# Rendering pipeline (WebGL) and performance notes

This page explains how Cellucid renders millions of points efficiently, and what changes are likely to cause performance regressions or WebGL failures.

It is written for contributors who might touch:
- `cellucid/assets/js/rendering/*`
- state→viewer synchronization code (`cellucid/assets/js/app/state/*`)
- overlay features (smoke, connectivity edges, vector fields)

## At a glance

**Audience**
- Computational users: read “Common failure modes” + “Troubleshooting”.
- Developers: read the whole page before changing renderer APIs or buffer formats.

**Time**
- 30–60 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`
- {doc}`06_state_datastate_and_events`

---

## Renderer entry point

The viewer is created in `cellucid/assets/js/rendering/viewer.js` via:
- `createViewer({ canvas, labelLayer, viewTitleLayer, onViewFocus })`

Hard constraints:
- **WebGL2 only**. If `canvas.getContext('webgl2')` fails, the app throws early.
- The context is created **single-sampled** (`antialias: false`). Antialiasing is
  not a context attribute any more; the app owns it, so do not add the attribute
  back. See the multisample resolve target below.

---

## Major render subsystems (what gets drawn)

`viewer.js` orchestrates multiple render passes/subsystems:

- **Scatter points (main cloud)**
  - Backend: `HighPerfRenderer` (`cellucid/assets/js/rendering/high-perf-renderer.js`)
  - Receives: positions (float32, one component per embedding axis), colors
    (RGBA uint8 packed), and alpha/transparency.
  - Notes:
    - Spatial indexing is dimension-aware (binary tree / quadtree / octree), used for picking/LOD.
    - Positions and colors are **two separate vertex streams** — see
      “Point vertex streams” below.

### Point vertex streams

Full detail is drawn from two `ARRAY_BUFFER`s bound into one VAO:

| Attribute | Stream | Layout | Bytes per point | Usage |
| --- | --- | --- | --- | --- |
| 0 `a_position` | `buffers.positions` | Float32, one component per embedding axis (x, xy or xyz), tightly packed | 4, 8 or 12 | `STATIC_DRAW` |
| 1 `a_color` | `buffers.colors` | RGBA8, normalized | 4 | `DYNAMIC_DRAW` |

Positions are never padded to XYZ. Every view's positions — `DataState`, the
viewer's per-view cache, the renderer, the spatial index, the LOD compact
store, snapshot geometry and the edge position texture — hold exactly as many
components per cell as the embedding has, and that count is the view's
dimension level (`rendering/position-layout.js`). The shaders keep a `vec3`
`a_position`; the attribute is bound with 1, 2 or 3 components and WebGL
supplies 0 for the axes a 1D or 2D stream does not carry, which is exactly the
value of an axis the embedding does not have. A 2D embedding therefore costs
8 bytes a point of positions instead of 12 on both the CPU and the GPU, and is
normalized in place on load instead of being expanded into a padded copy. CPU
code that needs a 3D point reads it with `readPositionXYZ`, which writes the
missing axes as zero without allocating. Because a view's dimension *is* its
position layout, a view changes dimension only by publishing that dimension's
positions (`updatePositions` / `setViewPositions`).

Alpha is not taken from the color stream while a dataset is loaded: it lives in
a separate R8 texture (one texel per point), so a visibility change touches
neither stream.

Each full-detail store is allocated straight from the published typed array
(`bufferData(positions)`, `bufferData(colors)`), so nothing is repacked and the
renderer holds no client copy of either. A recolor is one in-place
`bufferSubData` of the RGBA array; the positions are not re-sent.

Level of detail keeps one **compact prefix** for the live geometry (the live
positions have one dimension, so there is exactly one): the largest reduced
level's points in hierarchical order, held as the same two streams (a compact
position store with the live dimension's components and a compact color
store, one VAO), plus one R32UI index texture, shared by every reduced level,
that maps a compact rank back to its cell for the alpha lookup. That is the
whole LOD generation; there is no element buffer. The
compact stores are gathered through one fixed 65,536-point staging block
(`GATHER_CHUNK_POINTS`), positions and colors in separate passes, colors one
32-bit word at a time — which is why color arrays handed to the renderer must
start on a 4-byte boundary. A recolor rewrites only that compact color store,
in place, the next time the live dimension renders.

A kept (snapshot) view of the live geometry draws its reduced levels from the
same compact positions, with its own compact color store: 4 bytes per
maximum-prefix point, gathered from its published colors when they or the
live owner change, and retired with the view. A kept view with its own
geometry draws its own tree through a per-view element buffer.

The compact copy is deliberate. Drawing reduced levels through an original-ID
element buffer straight from the full-detail streams would remove the copy and
its recolor cost, but it turns sequential vertex fetch into a gather, and that
measured 2.0–2.4× slower per reduced-level frame at 5M points (Chromium, Apple
M1 Pro, ANGLE Metal). Kept views drew that way until they got their own compact
color stores: 12.0 ms against the live pane's 3.8 ms for the same
1,333,334-point level of a 10,000,000-point embedding. Splitting the streams did not cost vertex throughput:
the split full-detail draw measured within noise of the old 16-byte
interleaved record at 5M points.

- **Smoke / volumetric density**
  - Backend: `SmokeRenderer` (`cellucid/assets/js/rendering/smoke-cloud/`)
  - `viewer.buildSmokeVolumeGPU(...)` builds the visible-point volume through
    one WebGL2 transaction: bounded position + colour batches → R32F density
    splat atlas (eight trilinear corners per cell) and RGBA32F colour splat
    atlas (one texel per cell) → logarithmic GPU maximum reduction of the
    density → normalized RGBA8 atlas → RGBA8 3D texture with mipmaps. Red is
    the normalized density; green, blue, and alpha are the cells' mean display
    colour premultiplied by that density, so filtering and mipmapping average
    colour by the cells that produced it (`density-shaders.js`).
  - A flat cloud — a 1D or 2D embedding, or a planar 3D one, thinner than two
    voxels along an axis — is spread along that axis into a soft slab of fixed
    half-thickness (`SMOKE_SLAB_HALF_THICKNESS`, about four voxels of the 128³
    grid) by both splat passes, so a ray never steps over a voxel-thin sheet.
  - From each published volume, `buildDerivedVolumesGPU(...)` derives two more
    for the fixed key-light direction: a lighting volume (RG16F, half the grid:
    r = optical depth towards the key light, g = sky occlusion) that replaces a
    per-sample shadow ray, and an occupancy volume (R8, one texel per 4³ block
    of voxels: the largest density a sample in the block can read) that lets
    the ray march leap over empty blocks. Both are rebuilt only when the
    density volume or the light direction changes; the live density and
    absorption scale them in the shader.
  - The product path performs no `readPixels` synchronization. It publishes the
    new texture and bounds before completing its notification, preserves the
    prior volume on failure, and returns only a frozen non-owning summary.
  - The volume is a function of visibility, positions, and colours. The viewer
    reports every live position or colour republication (a field, colormap,
    category colour or log-scale change, a dimension switch, a new dataset)
    through `setSmokeSourceChangeHandler`, and the render controls mark the
    volume dirty exactly as they do for a visibility change.
  - The exact UI grid inventory is 32³, 48³, 64³, 96³, and 128³. The builder
    accepts exact integer sizes from 8 through 128 and rejects values outside
    that range before allocation. At 128³, its 2D atlas is 1536×1408—below
    WebGL2's required 2048 minimum `MAX_TEXTURE_SIZE`—and the published RGBA8
    volume is 8 MiB. The transient float atlases peak below 50 MiB, which is
    why the grid is not the memory lever it was.
  - Per frame, the ray march takes `Ray quality` samples per voxel (0.5–2,
    capped at `MAX_SMOKE_MARCH_STEPS` = 512 per ray) into an offscreen target
    of 0.25–1.00× the canvas, never more than 2560×1440 pixels
    (`MAX_SMOKE_TARGET_PIXELS`), and upsamples it with an alpha-gated
    bilateral composite. The noise volumes are 32³–128³. The ray-start jitter
    is a fixed blue-noise tile, so a still camera with `Animation speed` 0
    produces identical frames, and `SmokeRenderer.needsContinuousFrames()`
    then lets the viewer stop drawing.
  - A build requires a clean GL error state, `EXT_color_buffer_float`,
    `EXT_float_blend`, `MAX_3D_TEXTURE_SIZE` ≥ 128 and `MAX_TEXTURE_SIZE` ≥
    1536. These are checked once per context, up front and without allocating,
    by `detectSmokeSupport()` (`rendering/smoke-cloud/smoke-capabilities.js`,
    exposed as `viewer.getSmokeSupport()`); every build entry point refuses with
    a `SmokeUnsupportedError` whose message is user-facing and whose `detail`
    names the missing capability. Choosing smoke on such a device, or restoring
    a session saved in smoke, keeps Points and shows that message. There is no
    alternate density implementation. `OES_texture_float_linear` is not needed:
    every float32 texture smoke uses is NEAREST and read with `texelFetch`.
  - Per-build framebuffers, vertex arrays, staging buffers, atlases, and
    reduction textures are deleted on success and failure, and every touched
    GL binding/capability is restored. Programs plus the corner/quad buffers
    are cached per context, disposed with `SmokeRenderer`, and invalidated on
    context loss. The viewer requires a page reload after context loss.
  - The user notification duration covers the synchronous renderer call,
    including validation and publication. The `[GPU Splat] … submission`
    console duration starts after preflight and measures command submission.
    Neither duration is a GPU-completion measurement.

- **Connectivity edges**
  - Backend: instanced lines + edge textures (in `viewer.js`)
  - Visibility is derived from filter alpha and LOD; edges can be capped for UI performance.

- **Highlight tools**
  - Backend: `HighlightTools` (`cellucid/assets/js/rendering/highlight-renderer.js`)
  - Handles interactive selection and highlight overlay rendering.

- **Centroids**
  - Small, separate centroid shader program (centroid count is typically small).

- **Overlays**
  - Overlay framework: `OverlayManager` + overlay context
  - Vector field / velocity overlay: `VelocityOverlay` (`cellucid/assets/js/rendering/overlays/velocity/velocity-overlay.js`)

- **Multisample resolve target (antialiasing)**
  - Backend: `createSceneMsaaTarget` (`cellucid/assets/js/rendering/scene-msaa-target.js`)
  - The scene is drawn into an app-owned multisampled renderbuffer and blitted
    onto the default framebuffer at the end of each frame — the same sample
    count, coverage arithmetic and resolve the browser performed for a
    multisampled default drawing buffer, and the application owns the switch.
  - `SCENE_MSAA_TARGET_SAMPLES = 4`, or the device maximum when that is lower.
    A device reporting fewer than 2 samples has no target at all and the
    setting is reported as unavailable.
  - With antialiasing off, `beginFrame` binds framebuffer zero and
    `resolveFrame` does nothing: the frame costs exactly what it cost before
    this module existed, and nothing is allocated until it is switched on.

:::{warning}
Two caller contracts follow from this and are easy to break:

- Anything drawing into the scene must target `getSceneFramebuffer()` rather
  than assuming framebuffer zero. The smoke renderer and the velocity overlay
  both bind a framebuffer of their own mid-frame and must come back to the
  scene, not to the default.
- Anything reading pixels back must read **after** `resolveFrame`, from the
  default framebuffer, or read the resolved copy.
:::

---

## Viewer public API (how state becomes pixels)

The viewer exposes a small API surface used by `DataState` and `main.js`.

### Core dataset load

- `viewer.setData({ positions, colors, transparency, dimensionLevel })`
  - Initializes point buffers and builds spatial indices.
  - `positions` hold exactly `dimensionLevel` Float32 values per point.
  - `colors` are RGBA uint8 packed; `transparency` is a separate `Uint8Array` holding exactly 0 or 1 per cell, used for filtering.

### Incremental updates (critical for performance)

- `viewer.updateColors(colors)`
  - Rewrites the 4-byte-per-point color stream in place on the next frame, plus
    4 bytes per compact LOD point when the live dimension, which owns the LOD
    resources, next renders. Positions, the alpha texture and the LOD topology
    are not touched.

- `viewer.updateTransparency(alphaArray)` (a `Uint8Array` of 0/1; a Float32 or out-of-range array is rejected)
  - Updates alpha/visibility texture/buffer.
  - Also triggers highlight/overlay rebuild hooks (e.g., selection needs the new visibility mask).

- `viewer.updatePositions(positions, dimensionLevel)`
  - Used for dimension switching (same point count, different embedding).
  - Must preserve the same indices; `positions` must hold exactly
    `pointCount * dimensionLevel` values, and other lengths are rejected.
  - It is a new data publication: both streams are re-allocated straight from
    their arrays (4 bytes per embedding axis + 4 bytes a point, nothing
    repacked), the alpha texture is
    republished, and spatial indices are rebuilt when LOD or frustum culling
    needs them. With LOD on, the spatial-index build dominates its cost.

Design rule:
- Prefer incremental updates (`updateColors`, `updateTransparency`, `updatePositions`) over full `setData` reloads.
Full reloads are the fastest way to introduce jank on large datasets.

---

## Multiview rendering (live + snapshots)

Cellucid can render:
- one live view
- **at most three** snapshot views (`MAX_SNAPSHOTS = 3`; a fourth raises a
  `RangeError`)

Key idea:
- A snapshot view is not “a second dataset”; it is a second **view context** with its own camera/dimension/filter/field choices, rendered from the same underlying point identity index space.

The viewer:
- tracks per-view dimension levels (`viewer.setViewDimension(...)` is called from state when available)
- shares positions with snapshot views: a snapshot of the live geometry borrows
  the live positions array and GPU position store, with no copy, and a snapshot
  with its own geometry adopts the array it was published with. A published
  positions array is never written afterwards.
- uses optimized paths for alpha/transparency sharing so it does not re-upload N full buffers per view

### The presented layout: one source of truth for where views are

The canvas fills the window and the sidebar overlays its left edge. The
**presented area** is the part of the canvas the sidebar leaves uncovered, and
every view is laid out relative to it by one pure function,
`resolvePresentedLayout()` in `assets/js/rendering/viewer/presented-view-layout.js`.

- `ui/modules/sidebar-controls.js` measures the sidebar and publishes the
  fraction of the canvas it covers from the left
  (`computeLeftOcclusionRatio()` in `rendering/viewport-center.js`) to
  `viewer.setViewportLeftOcclusionRatio()`. It republishes on a toggle, during
  a resize drag (from cached geometry, without a layout read), on window and
  canvas resizes, and on responsive breakpoints. The open/close slide is a CSS
  transform no `ResizeObserver` sees, so while it runs the presented area stays
  the one beside the *open* sidebar (it shrinks at once on opening and grows
  only when the closing slide has ended): panes are never under the moving
  sidebar and keep equal areas in every frame. On the phone layout the open
  sidebar is a full-screen sheet the user closes to reach the plot, so it
  covers nothing. The hidden sidebar is `inert`.
- The viewer resolves the layout from the canvas backing-store and CSS sizes
  (so device-pixel-ratio changes are included), the occlusion ratio, the
  layout mode, the focused view, and the view inventory. It is memoised: it is
  recomputed only when one of those inputs changes, together with one
  projection matrix per pane, and is otherwise one frozen object.
  `viewer.getPresentedLayout()` returns it and
  `viewer.onPresentedLayoutChange(listener)` delivers each new one.
- **Grid** (`Grid compare` with two or more views): the presented area is split
  into `cols × rows` cells of identical integer framebuffer size, starting at
  the first framebuffer column right of the sidebar. No cell is under the
  sidebar, every cell has the same area, and each renders with a symmetric
  projection and its own scissor rectangle.
- **Single view**: the pane keeps the whole canvas (the scene continues under
  the sidebar), and its projection is shifted horizontally so the principal
  point - the camera target in Orbit and Planar, the look direction in
  Free-fly - sits at the centre of the presented area. This is what a single
  view always showed in Orbit and Planar; Free-fly now uses the same centre,
  so switching modes no longer shifts the picture.

Each pane carries its WebGL viewport (`viewport*`, framebuffer pixels,
bottom-left origin), its render rectangle (`css*`) and presented rectangle
(`visibleCss*`) in canvas-local CSS pixels, and `projectionCenterNdcX`. Every
mapping between the screen and a view reads these panes: the render loop's
viewports, scissors and per-pane projections (so frustum culling and LOD see
the real pane shape), `getViewportInfoAtScreen()` and therefore click/tap
focus, wheel zoom (which focuses the pane under the cursor), picking (nothing
outside a pane's presented rectangle can be picked), hover delivery to
Jupyter, every highlight tool (lasso polygons are clipped to the pane's
presented rectangle; proximity and KNN indicators are clipped like the WebGL
scissor), projectile aiming and the charge ring, centroid labels, view title
chips, overlays, `getPresentedViewState()` (figure export and its preview),
and the selection bar, which is centred in `presentedArea`.

### Layout arithmetic, and why it changes the cost

```
cols = viewCount <= 3 ? viewCount : ceil(sqrt(viewCount))
rows = ceil(viewCount / cols)
paneWidth = floor((canvasWidth - presentedLeft) / cols)
paneHeight = floor(canvasHeight / rows)
```

So: 2 views → 1×2, 3 → 1×3, 4 → 2×2. Figure export arranges multi-panel
figures with the same column rule (`getPresentedGridColumns()`).

`paneHeight` is what the point vertex shaders receive as the viewport-height
uniform, and perspective point size is proportional to it. Point **area** is
therefore proportional to `1/rows²`, and total fill across all panes is
proportional to `viewCount / rows²`, i.e. to `cols / rows`.

| Views | cols × rows | Relative total fill |
|---|---|---|
| 1 | 1 × 1 | 1 |
| 2 | 2 × 1 | 2 |
| 3 | 3 × 1 | 3 |
| 4 | 2 × 2 | 1 |

This is the **fill** part of the cost only. Every view also pays for every
point it draws, about 5–7 ns per on-screen point sprite on an Apple M1 Pro
through ANGLE Metal, whatever its size. At the default point size on a large
dataset that per-point part dominates and views add up: a 10,000,000-point 2D
embedding took 30 ms for one view and 147 ms for four. The fill ordering above
matters only when points are large enough to overlap heavily, and only while
size attenuation is non-zero; at zero attenuation point size ignores the
viewport. Level of detail is the lever for multiview: with it, the same four
views took 13.5 ms.

:::{warning}
Do not “optimise” the viewport-height uniform to the canvas height for
consistency across panes. Feeding the full canvas height into a 2×2 grid would
quadruple fill cost and make four views genuinely four times a single view.
:::

---

## Overlays (vector fields) and why they’re tricky

Vector field overlays are conceptually simple (“animate particles along vectors”), but interact with:
- dimension switching (vectors are dimension-specific),
- filtering (particles should respect visibility),
- multiview (each view has its own dimension level and visibility),
- performance (particle systems can be expensive).

Implementation notes:
- The overlay is “opt-in”: it is initialized only when enabled/needed.
- When visibility changes, `viewer.updateTransparency` marks the overlay “visibility dirty” so it can lazily rebuild spawn sources.

Developer advice:
- Avoid coupling overlay state to UI controls directly; route it through state events and viewer methods.

---

## What a frame actually costs

Before changing anything in the draw path, know the current baseline: **a settled
scene is not drawn at all, and a drawn frame that changes nothing transfers no data
to the GPU.**

```
render()  (self-rescheduling rAF; draws only when something can have changed)
  │
  ├─ frame gate ── invalidated, animating, or resized? ── no, settled ─▶ (no draw)
  │                                    │ yes
  ├─ dirty-flag gate ── buffer dirty? LOD dimension dirty? ─── no ──┐
  │                                    │ yes                        │
  │                                    └─▶ flushBufferUpdates()     │
  │                                                                 │
  ├─ frustum cache ── MVP + dimension + bounds unchanged? ─ yes ─────┤ (no work)
  │                                    │ no                         │
  │                                    └─▶ re-extract planes        │
  │                                                                 │
  ├─ index publication ── admitted leaf set unchanged? ──── yes ─────┤ (reuse EBO)
  │                                    │ no                         │
  │                                    └─▶ upload element buffer    │
  │                                                                 │
  └─ useProgram · uniforms · bind alpha texture · drawArrays ◀───────┘
```

The frame gate decides whether a frame is drawn at all. Any mutating viewer or
renderer call, any viewer input that changes something (a drag, free-look, a
selection gesture, `Alt`, navigation keys outside text fields, releasing a held
key), and any canvas resize marks the next frame for drawing — hovering or typing
into a sidebar field draws nothing; held navigation keys, inertia, smooth zoom, a grid fade, smoke mode
while it animates or is still waiting for its noise volumes, enabled overlays, visible highlights (their glow animates), and live projectiles
keep frames drawing, and so does a camera motion still waiting for its
`camera-settled` event. After the last such cause the viewer keeps drawing until
both one second has passed and 12 frames have been drawn
(`FRAME_SETTLE_MS`, `FRAME_SETTLE_MIN_FRAMES`), so multi-frame convergences
such as LOD hysteresis finish even where frames arrive more than a second apart,
and then stops until something changes. A frame that throws marks the next
frame for drawing, so a partial frame is redrawn. The drawing buffer is preserved,
so the canvas keeps showing — and reading back — the last drawn frame.

Inside a drawn frame, three independent gates keep the work proportional to what
changed, each cheap to check:

1. **Buffer dirty flags** — set by colour updates, cleared on flush.
2. **Frustum cache** — an exact element-wise comparison of the model-view-projection
   matrix, the dimension and the bounds. An idle camera re-extracts nothing.
3. **Index publication watermark** — a camera move that changes the matrix but
   not the *ordered set of admitted spatial nodes* reuses the published element
   buffer. Even the post-upload error check is gated on the allocation
   watermark, because running it every frame was measured in milliseconds.
   Highlight buffers follow the same rule: `gl.getError()` runs on the first
   allocation and when the store grows past the size already verified, not on
   every republication.

The culled full-detail element buffer lists the visible points in **ascending
index order**, not leaf by leaf. Every spatial leaf spans the whole point-ID
range, so leaf order made the GPU fetch both vertex streams scattered across
memory: at 5,000,000 points with everything in view, a culled static frame took
30.0 ms against 9.4 ms unculled (Chromium, Apple M1 Pro, ANGLE Metal). In
ascending order it takes 9.7 ms. Storage order is also the unculled draw order,
so a culled frame blends translucent points exactly like an unculled one, which
`tests/browser/frustum-culling-parity.spec.mjs` checks pixel for pixel. A small
visible set is sorted directly; a large one is marked in a one-byte-per-point
array and read back with one linear scan (`_writeAscendingVisibleIndices`).
Reduced LOD levels were already emitted in compact-rank order. Do not go back to
concatenating leaf index arrays.

Binding the alpha texture is a bind and two uniform writes; it never uploads.

The practical rule: **if you add work to the draw path, it must sit behind a gate
of its own.** “It is only a few microseconds” is how the idle frame stops being
free.

### Where the GPU time goes, and which optimisations are already closed

Measured with GPU timer queries on **one Apple M1 Pro through ANGLE Metal**,
1440×1000 at DPR 1, 10,000,000 synthetic points, LOD and frustum culling off,
counterbalanced arms with a byte-identical twin in every set. The absolute
milliseconds belong to that machine; the ratios are the durable part.

| Stage | Measured | Verdict |
| --- | --- | --- |
| Vertex | 2.62 ms, flat at every point size — 4.2% of the frame | fetch-bound; hoisting the redundant matrix product and the per-vertex `tan()` measured **no effect** |
| Fragment | 190× the fragments costs 1.82× the time | only ~5% of rasterised fragments survive the depth test to be shaded |
| Per-sprite rasterisation | ~5.7 ns/sprite, ~175M sprites/s | **this is the floor** |

Fill only reaches parity with the per-sprite cost at `gl_PointSize` ≈ 15 px.
Below that the frame is primitive-bound; above it, fill-bound.

**Do not re-open these.** Each was measured and rejected, and the numbers are in
the maintenance ledger rather than being folklore:

- **Fragment-shader simplification of any kind, in 3D views.** Fog, lighting,
  the alpha `texelFetch`, the round `discard`, the mere presence of `discard`,
  and the whole `full` → `light` → `ultralight` ladder (43 → 9 → 2 statements in
  the translated Metal) each stayed inside its own twin band at every point
  size. These measurements are 3D: depth testing rejects hidden fragments
  before shading. **2D views are different** — they draw without depth, so every
  overlapping, fog-translucent fragment is shaded and blended in order, per
  sample under MSAA. There, `full` and `light` now use a single `discard`, and
  `light` has a no-depth variant (`HP_FS_LIGHT_NO_DEPTH`, picked by
  `_pointProgramName()` for dimension ≤ 2) that writes a zero-alpha fragment
  instead of discarding, which leaves the pixels byte-identical. At 1,000,000
  cells in 2D with antialiasing, `light` went from 30.7 to 17.0 ms and `full`
  from 33.8 to 30.6 ms.
- **A different primitive.** `GL_POINTS` beat instanced quads and both triangle
  arrangements by 1.41–2.12× at one sample and by up to 1.90× with MSAA enabled
  at 4 samples, with the quad arms' rasterised and shaded fragment counts matched to
  five decimal places. Halving the primitive count with one triangle per point
  was the *worst* arm.
- **Order-independent transparency, a depth prepass, storage reordering, and the
  vertex-ALU hoist.** All measured; none won.

What did move: **antialiasing**, worth about a fifth of the frame at point size
0.75 and a third at large sizes in 3D — and about 10× in a dense 2D view, where
nothing is depth-rejected (1,000,000 cells: 33.8 ms on against 3.4 ms off). It is a live user setting (`#hp-antialias`) and
no longer a context-creation attribute: `scene-msaa-target.js` owns a
multisampled renderbuffer resolved by a blit each frame, so a change applies to
the next frame with no reload. The stored preference is three-state — `auto`
(the default, and what an absent key means), `on`, `off` — and `auto` resolves
from the dataset's cell count against `AUTOMATIC_ANTIALIAS_CELL_LIMIT = 5_000_000`,
re-decided on every dataset publication.

**Shader quality selects shading, never geometry.** All three point vertex
shaders — and the copy `HP_VS_HIGHLIGHT` keeps of the formula — clamp
`gl_PointSize` to the same `[0.5, 128.0]`. A higher ceiling in the lighter
shaders once made them draw *larger* sprites than `full`, so the "faster"
quality was up to 1.93× slower at close range. `tests/point-size-clamp-contract.test.mjs`
pins the three bounds together; keep them identical in any fork.

**The point-size curve lives in one module**,
`cellucid/assets/js/rendering/point-size-scale.js`: the slider curve (anchored at
0.25 at position 0 and 200 at position 100, with the domain extended down to
position -24 so stored sessions and published presets keep their meaning), the
viewer's bound, and `pointSizeSliderPositionForCellCount`, the size a dataset
opens at. It used to be a `0.25` in the sidebar module and a second `0.25` in
`viewer.setPointSize`; do not reintroduce a second copy — lowering the reachable
minimum in one place and not the other makes the slider throw at its own low end.

The reachable minimum is 0.050, well below the `0.5` clamp floor, and that is not
a contradiction: sizes under about 0.5 device pixels are only distinguishable
where perspective size scaling lifts the rendered size back above the clamp.

---

## Performance footguns (common mistakes)

### 0) Rejecting invisible points

Filtered-out points are not skipped by shrinking them. Setting `gl_PointSize` to
zero does **not** remove a point: the driver’s minimum aliased point size is one
pixel, so the point is still assembled, rasterised and shaded, and only then
discarded in the fragment stage — and the ultralight fragment shader has no
`discard` at all, so hidden points were also writing depth and occluding visible
ones.

The point vertex shaders now push hidden vertices outside the clip volume as
well as zeroing the size, which removes them before rasterisation. This is what
makes a heavily filtered view cheaper than an unfiltered one.

**Every** vertex shader that draws per-cell geometry does this now, including
the highlight-dot and centroid shaders, which previously only zeroed the point
size. The centroid case was not theoretical: a category with no visible members
has its alpha set to zero by the field summary, and centroids are drawn with
depth writes enabled, so those invisible sprites were occluding real points.

If you add or fork a point vertex shader, it must do the same. This is enforced
rather than documented: `tests/high-perf-hidden-point-rejection-contract.test.mjs`
sweeps every module in the tree that carries shader source, so a new shader
cannot be added without the rejection.

### 1) Accidental hot-path allocations

Avoid:
- allocating new arrays/objects inside per-frame loops
- creating new typed arrays on every slider change when a reusable scratch buffer would work

Prefer:
- reusing scratch buffers (float32/uint8) sized to `pointCount`
- debouncing expensive work in UI modules

Smoke density never makes a full duplicate position or colour array. It
retains the dataset-owned positions, colours, alpha, and optional
outlier-quantile arrays, validates them in place, and streams visible cells
through two transient buffers fixed at 262,144 cells: positions in the source
layout (at most 3 MiB, for a 3D embedding) and RGBA8 colours (1 MiB), 4 MiB per
batch. Large inputs are submitted in bounded batches.
Do not replace this with a buffer sized to the complete dataset or introduce
GPU→CPU volume readback. The grid cap and early release of completed reduction
textures bound the other transient memory.

### 2) Full buffer re-uploads when only alpha changed

Filtering often only needs alpha changes.
Do not rebuild/re-upload positions or full color buffers for filtering.
For the same reason, do not merge the position and color streams back into one
interleaved record: a recolor would then have to repack and re-upload the
positions beside every color, four times the bytes it needs.

Alpha itself is not re-uploaded wholesale either. Visibility lives in a texture
whose width is the device’s maximum texture size, and `updateAlphas` compares
the incoming array against the uploaded bytes row by row, coalesces the dirty
rows into maximal runs, and issues one sub-image upload per run. Above a bounded
number of disjoint runs it collapses to a single bounding range, so the cost is
never worse than the whole-texture upload it replaced. A republication where
nothing moved uploads zero bytes and returns a signal saying so.

Two invariants to preserve:

- The comparison must stay row-major over the same layout the texture uses;
  changing the packing without changing the comparison silently disables it.
- The “did anything move?” return value is load-bearing. It is what stops a
  no-op filter edit from repacking highlight geometry, which costs 16 bytes per
  selected cell per view — megabytes on a large selection.

### 2b) Cell-index inputs

Visibility APIs that accept a list of cell indices must handle typed arrays
directly and must not build a `Set` of them. A `Set` caps out at 2²⁴ entries, so
a selection above ~16.7 million cells was not slow — it *threw*. Presence is now
tracked in a `Uint8Array` bitmask sized to `pointCount`, which is O(1) per index,
allocation-free per call, and has no ceiling below the dataset size.

(There is a second, unrelated 16.7 million threshold in the shaders: float32
cannot address texel coordinates exactly above 2²⁴, which is why the index
textures are integer-typed. Do not conflate the two.)

### 3) Doing DOM work in render-critical flows

The renderer should not depend on DOM reads/writes at render frequency.
UI changes should be event-driven and decoupled from the render loop.

---

## Common failure modes

### WebGL2 unavailable

Symptom:
- immediate error: “WebGL2 is required but not supported…”

Causes:
- old browser
- corporate GPU policy / disabled WebGL
- remote desktop environments

Mitigation:
- test with a modern browser first
- provide a clear error message (already does)

### WebGL context lost

Symptoms:
- canvas freezes
- DevTools shows “WebGL context lost”

Causes:
- GPU memory pressure (large datasets + high settings)
- browser tab throttling
- driver instability

Mitigation:
- reduce smoke settings / disable heavy overlays
- lower dataset size for reproduction
- reload the page after context loss; restored contexts are not reused by the
  current viewer
- use the performance troubleshooting docs: {doc}`../n_benchmarking_performance/index`

---

## Troubleshooting (renderer-level)

### Symptom: “Points disappear / everything is invisible”

Likely causes:
- `categoryTransparency` got set to ~0 for all points (filters/outlier threshold).
- Alpha array length mismatch caused viewer to reject update.

How to confirm:
- In console: check `window._cellucidState.filteredCount` and `window._cellucidState.categoryTransparency`.
- Put a breakpoint in `viewer.updateTransparency`.

Fix:
- Ensure filter logic updates transparency correctly and recomputes counts.

### Symptom: “Dimension switch breaks rendering”

Likely causes:
- `updatePositions` received the wrong length (`pointCount * dimensionLevel`
  required).
- Embedding fetch failed and returned empty/invalid buffer.

How to confirm:
- DevTools → Network: check embedding requests.
- Console: look for `Viewer live positions must contain exactly … values`.

Fix:
- Validate embedding buffers before calling viewer updates.
- Ensure the dimension manager returns consistent point ordering.

---

Next: {doc}`08_ui_modules_map`.
