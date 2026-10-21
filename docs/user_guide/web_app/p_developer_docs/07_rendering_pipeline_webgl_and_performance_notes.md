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
viewer's per-view cache, the renderer, the spatial index and snapshot geometry — hold exactly as many
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

**Per-point textures are linear textures** (`rendering/linear-texture.js`): a
`TEXTURE_2D_ARRAY` whose rows are a power of two wide, stacked into as many
layers of `MAX_TEXTURE_SIZE` rows as the count needs; a count that fits one
row gets one row exactly as wide as the count, so its first upload covers the
whole image (Firefox clears an image a first upload covers only in part). The
alpha texture, the permutation and the permute pass's source all use it, and a
shader finds item `i` with a shift, a mask and one division (`linearTexel`). A single 2D texture
used to cap them at `MAX_TEXTURE_SIZE²` items: 67,108,864 on a GPU that
reports 8,192 (Firefox on macOS), so a 3D position gather (three words a
point) failed there above 22.4M points.

**Unchanged panes keep their pixels.** The scene target (the multisampled
target, or the default framebuffer with antialiasing off) is not cleared
between frames. Every change that can alter what a pane draws advances a scene
generation (`invalidateFrame`); a camera gesture alone does not
(`invalidateCameraFrame`), and time-driven effects other than the camera's own
motion (fades, smoke, the velocity overlay, the highlight pulse, projectiles,
LOD still converging) advance it every frame they run. A grid pane whose
generation, camera, projection, rectangle and navigation mode are unchanged
is skipped and keeps what it drew; a redrawn pane clears its own rectangle.
Orbiting one pane of four redraws one pane, and a frame that would repeat what
a single view already shows (the settle frame after a change) draws nothing.
A new target, size or layout draws every pane afresh.

**Input lands in the frame that draws it.** Animation-frame callbacks run in
the order they were registered, and the viewer re-registers its own at the
start of every frame, so a control that waited for its own animation frame ran
after the viewer and showed every slider step a frame late. Work that changes
what is drawn (coalesced slider and legend input, camera-path playback) is
queued with `runBeforeNextFrame` (`rendering/frame-phases.js`) and runs at the
top of the viewer's frame, before it decides what to draw; `runAfterNextFrame`
runs once the frame is drawn. While no viewer drives frames, each callback
runs from an animation frame of its own.

**No draw submits more than 2^24 vertices.** Firefox refuses a draw over
30,000,000 vertex ids with `INVALID_OPERATION` and draws nothing, so every draw
whose count grows with the data -- the point passes, the permute pass, the
highlight pass and figure export -- goes through `drawPointsInSegments`, which
issues consecutive ranges of at most 16,777,216 vertices. Draws execute in
submission order, so the image (and transform-feedback output) is that of one
draw; a 100M-point frame costs at most 6 extra calls.

**Draw order.** A geometry's points are stored on the GPU in the draw order of
the spatial tree over them (`SpatialIndex#drawOrder`, computed by
`high-perf/draw-order.js`): grouped by LOD band, then by draw unit (a tree
leaf in 3D, one of 64 coarse nodes in 1D and 2D), then by ascending cell
index. Band `b` holds the cells the hierarchical order first
admits at LOD level `b`; the last band holds the full-detail remainder. So:

- a reduced LOD level is one `drawArrays` **prefix** of the stored points;
- full detail is `drawArrays` over all of them;
- frustum culling collects the runs of consecutive tree leaves inside the
  frustum (a subtree wholly inside is one run, read from its node) and draws
  each run as one contiguous range per band. Where the ranges exceed the
  budget -- 1,024 with `WEBGL_multi_draw` (one `multiDrawArraysWEBGL` call),
  256 without it (Firefox, one `drawArrays` per range) -- the gaps between runs
  that hide the fewest points, counted band by band, are bridged first
  (`SpatialIndex#writeDrawRanges`), so the drawn total stays close to the
  exact cull. No per-point index is ever written or uploaded after the order
  is built, and a camera that does not move re-culls nothing.

Neighbouring points are neighbours in memory and on screen, which is also what
makes the stored order faster at full detail: 10,000,000 points drew in 18 ms
against 68 ms in source order on Chromium (15 against 57 on WebKit, 14 against
57 on Firefox), a reduced level of 8,000,000 points in 15 ms against 14 ms for
the compact copy it replaced, and a zoomed-in culled frame in 7–8 ms against
21–31 ms for the per-point index buffer it replaced (Apple M1 Pro; standalone
WebGL2 harness on the real draw order, translucent 2-px points, all cases in
one run on a loaded machine).

Positions are per geometry, so each geometry carries its order; colours are
per view, and each view's colour store is gathered into the order of the
geometry it shows. Alpha stays in source order in its R8 texture, read through
the order's R32UI permutation texture (`u_useLodIndexTex`: texel `k` is the
cell stored in slot `k`), so a visibility change still touches neither stream.
Picking, highlights, overlays, analysis and export all read source-order CPU
arrays; nothing but the point draw reads the stored order.

**Trees are built off the main thread.** The tree over a geometry and its
draw order (banded by LOD when LOD is on) are built together in a worker
(`high-perf/spatial-tree-worker.js`, the same `buildSpatialTreeJob` either
way), so loading a dataset, turning LOD or culling on, keeping a view or
showing another representation never runs a tree build on the main thread
(about 2.2 s at 10M points before). Until a view's tree lands, the view draws
full detail from its store as it is; when it lands, the renderer installs it
and asks the viewer for a frame (`onAsyncTreeReady`). A build whose positions
were replaced before it returned is dropped. Positions of any size go to the
worker, in transferred slices of 64 MiB, one task each, so no main-thread task
copies a whole geometry and the main thread never holds a second copy. A query that needs a tree at once (a click-time pick on a pane whose
index has not arrived, projectiles) still builds a plain tree directly; the
viewer asks the worker for every large pane's index when cell info is on and
when a selection tool is armed, so this is rare.

**A store adopts its tree's order in the frame that first needs it,**
whenever a tree exists, in every dimension (a hover tree counts too): in 3D
the order draws a 10M-point full-detail frame about 3x faster, in 2D about
3x faster too (21 ms against 67 ms, no vsync, Apple M1 Pro). Before a tree
exists a geometry is stored in source order, straight from the published
arrays with no client copy. The adoption uploads the permutation (4 bytes a
point), then rewrites the store's positions and the colours of every view
that draws it in place on the GPU (`GpuPermuter`,
`high-perf/gpu-permute.js`): the source-order words go up into a transient
R32UI texture, one texel a word, and one transform-feedback pass writes each
stored slot's words, fetched through the order's R32UI permutation texture,
bit for bit. At 10M 3D points the frame that adopts the order spends about
25 ms of CPU on it (positions 13, permutation 5, colours 5), once; turning
LOD on used to block the page for 2.3 s. A kept view's colour store is RGBA8,
4 bytes a point like the live one, and is written the same way, one word a
point (its alpha byte is never read: a kept view's alpha is its R8 texture).
It used to be packed RGB, 3 bytes a point, but Metal needs a 4-byte vertex
stride, so ANGLE kept a converted 4-byte copy of that store and converted it
again after every write. Nothing proportional to the point count
is gathered on the CPU, and no frame is kept alive for it. Colours are copied
as 32-bit words, which is why colour arrays handed to the renderer must start
on a 4-byte boundary.

A recolour rewrites the view's colour store in place: one upload of 4 bytes a
point and one GPU pass, whether or not the store is in a draw order. A kept
view's recolour does the same when the view alone owns its store and the
store holds the order it draws with; only a shared store, or one in another
order, is replaced.
A view that switches to another representation keeps the colour store it
leaves (keyed by that order); switching back while the colours are unchanged
draws from it again with no upload. Kept colour stores count against the
retention budget with retained geometry and are evicted first.

The 2D overlap order of translucent points therefore follows the spatial tree
(within a leaf, the later cell is drawn on top), and it is the same in the live
view, kept views, every LOD level, culled frames and figure export (which draws
in `getViewDrawOrder`). Nothing in Cellucid sets a draw order on purpose (no
"high values on top" or "bring to front" option exists), so no feature depends
on source order.

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
    bilateral composite. The noise volumes are 32³–128³, generated on the GPU
    over a few scheduling turns (`gpu-noise-generator.js`). A turn runs every
    cheap setup step that fits its time budget and stops at a program still
    compiling; each slice batch is a turn of its own, ending with a fence the
    next waits for. The batch doubles back (up to its caps) when the GPU
    finished it by the next turn, keeps its size when it was one turn late
    (the renderer's latency, whatever the batch), and halves only when it
    outlasted two turns, so a slow or software renderer never has seconds of
    noise work queued in front of the frames and input that follow, and is
    not left drawing one slice per two frames. The ray-start jitter
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
  - Backend: `rendering/connectivity-edge-layer.js` (which edges a pane shows:
    `rendering/connectivity-edge-selection.js`).
  - The graph stays on the CPU as loaded: `uint32` endpoints in storage order and
    one `float32` relative strength an edge (the `float64` weights are read once
    and not kept). Nothing is shuffled or copied.
  - A pane shows the first K of its visible edges — both endpoints pass the
    pane's filters and, with level of detail on, its current level — in one fixed
    random order: a keyed Feistel bijection of the edge indices, so the order
    costs no memory and every K-prefix is a uniform sample, the same in every
    session for the same graph. Selecting costs O(K / visible fraction).
    K is the edge slider's value: 250,000 (or every edge, if fewer) for a newly
    loaded graph (`CONNECTIVITY_DEFAULT_VISIBLE_TARGET`), the same in every
    pane. A pane added later selects on its first frame.
  - Only those K edges reach the GPU, one 28-byte instance record each (both
    endpoints' positions and the strength), so the vertex shader fetches no
    texture and tests nothing. GPU memory is K × 28 bytes a pane at any graph
    size; a camera move uploads nothing, a filter or LOD change re-selects and
    re-gathers K records, and a representation switch only re-gathers.
  - The visible-edge count for the slider's range is one sequential pass over the
    edges, run in time-boxed slices after the selection, so a filter change never
    waits for it.

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
  - Rewrites the 4-byte-per-point color stream in place on the next frame,
    gathered on the GPU into the live geometry's draw order when it has one, and drops the
    colour stores kept for other orders. Positions and the alpha texture are not
    touched.

- `viewer.updateTransparency(alphaArray)` (a `Uint8Array` of 0/1; a Float32 or out-of-range array is rejected)
  - Updates alpha/visibility texture/buffer.
  - Also triggers highlight/overlay rebuild hooks (e.g., selection needs the new visibility mask).

- `viewer.updatePositions(positions, dimensionLevel)`
  - Shows another representation in the live view (same point count, another
    embedding, of any dimension); it is the only way the live view changes
    dimension.
  - Must preserve the same indices; `positions` must hold exactly
    `pointCount * dimensionLevel` values, and other lengths are rejected.
  - It replaces only the positions: colours, transparency and the alpha
    texture are per cell, not per representation, and are kept. The position
    stream and spatial tree come from that representation's pooled geometry
    when a kept view shows it or it was shown lately (no upload, no tree
    build); otherwise the positions are uploaded once (4 bytes per embedding
    axis a point, nothing repacked) and spatial indices are rebuilt when LOD
    or frustum culling needs them. With LOD on, the spatial-index build
    dominates its cost. A failure restores the previous publication exactly.

Design rule:
- Prefer incremental updates (`updateColors`, `updateTransparency`, `updatePositions`) over full `setData` reloads.
Full reloads are the fastest way to introduce jank on large datasets.

---

## Multiview rendering (live + snapshots)

Cellucid can render:
- one live view
- **any number** of snapshot views. Each owns an RGBA8 colour store and an R8
  alpha texture (5 bytes a cell; positions are shared by geometry). Memory,
  not a count, limits them: when an allocation is refused, the stores of kept
  views not presented this frame are released first
  (`HighPerfRenderer._releaseDormantSnapshotStores`; the viewer reports the
  presented views every frame with `setPresentedViews`), and a released view
  is rebuilt from its own colours and alpha bytes when it is next drawn or
  edited (`_restoreSnapshotStores`), drawing exactly what it drew. Only a set
  of views shown at once that does not fit is refused, by the device.

Key idea:
- A snapshot view is not “a second dataset”; it is a second **view context** with its own camera/representation/filter/field choices, rendered from the same underlying point identity index space.

The viewer:
- tracks each view's layout key (`viewer.setViewEmbedding(viewId, { dimension, layoutKey })`, called from state); a layout change invalidates that view's LOD, frustum, picking and vector state
- shares geometry by representation: geometry is looked up by its positions
  array (`_findPooledGeometry(positions, dimension)`; there is one array per
  representation), so two views of one embedding share one GPU position
  stream and one spatial tree, and two embeddings of one dimension never do. A
  snapshot of the live geometry borrows the live positions array and GPU
  position store, with no copy. A published positions array is never written
  afterwards.
- retains geometry no view shows any more, so showing it again costs no
  upload and no tree build: within a budget of twice the largest geometry a
  view shows, between 64 MiB and 768 MiB, least recently used first
  (`evictRetainedGeometry`), and never at the cost of a new allocation:
  nothing kept is evicted while an allocation fits, and an allocation the
  driver refuses releases everything kept and runs once more
  (`_allocateGpu`, used by every GPU store and texture of the point renderer).
  A spatial tree whose build failed is built again on the next request once
  GPU memory has been released
- links cameras by representation: with **Link cameras** on, the views that
  show the same embedding share one camera and navigation mode
  (`viewer.setCamerasLocked()`, `viewer.getCameraLinkGroups()`); views of
  different embeddings always keep their own
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
- dimension and embedding switching (vectors are bound to one embedding),
- filtering (particles should respect visibility),
- multiview (each view has its own representation and visibility; a field is
  drawn in every view showing an embedding it is bound to),
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
render()  (rAF armed only while something can change; draws only then)
  │
  ├─ frame gate ── invalidated, animating, or resized? ── no, settled ─▶ (no draw; loop stops)
  │                                    │ yes
  ├─ dirty-flag gate ── colours dirty or in another order? ── no ──┐
  │                                    │ yes                        │
  │                                    └─▶ upload + GPU gather      │
  │                                                                 │
  ├─ frustum cache ── MVP + dimension + bounds unchanged? ─ yes ─────┤ (no work)
  │                                    │ no                         │
  │                                    └─▶ re-extract planes,       │
  │                                        collect leaf runs        │
  │                                                                 │
  ├─ draw ranges ── same runs and LOD band? ─────────────── yes ─────┤ (reuse)
  │                                    │ no                         │
  │                                    └─▶ write ≤ 1024 ranges (CPU)│
  │                                                                 │
  └─ useProgram · uniforms · bind alpha + permutation · draw ◀───────┘
```

The frame gate decides whether a frame is drawn at all. Once the scene has
settled the loop stops: no animation frame is requested until something asks
for one (`ensureFrameScheduled` from an invalidation, a held navigation key, a
canvas resize or queued frame-phase work), so an idle page does not wake
60–120 times a second. A frame that is drawn re-arms the loop, which runs its
settle frames and any time-driven effect, then stops again. Any mutating viewer or
renderer call (except a republication of what is already drawn: an unchanged
view label, or a highlight with the same cells, which is what a selection
preview sends on most pointer moves), any viewer input that changes something (a drag, free-look, a
selection gesture, `Alt`, navigation keys outside text fields, releasing a held
key), and any canvas resize marks the next frame for drawing — hovering or typing
into a sidebar field draws nothing; held navigation keys, inertia, smooth zoom, a grid fade, smoke mode
while it animates or is still waiting for its noise volumes, an overlay while it draws something that moves (velocity particles flowing in a presented view, or a spawn table being built for one; with no field for the layout, nothing visible, particles off at this level of detail, or in smoke mode it keeps no frames going), visible highlights in a quality whose highlight style animates (only `full`: its ring pulses and shimmers; every other quality draws a still ring, and a selection there costs nothing while the view is idle), and live projectiles
keep frames drawing, and so does a camera motion still waiting for its
`camera-settled` event and an adaptive LOD still stepping toward the level its
camera distance asks for (hysteresis moves one level per drawn frame; the
renderer reports it per view with `isViewLodConverging`). After the last such
cause the viewer draws two more frames (`FRAME_SETTLE_MIN_FRAMES`) and then
stops until something changes. Every multi-frame convergence is a cause of its
own, so the floor only covers the frame that shows a change. There is no
wall-clock window on top: one kept a fast display drawing 60–120 identical
frames after every change, and a renderer that takes seconds per frame now
spends two frames on a change, not a dozen. A frame that throws marks the next
frame for drawing, so a partial frame is redrawn. The drawing buffer is preserved,
so the canvas keeps showing — and reading back — the last drawn frame.

Inside a drawn frame, three independent gates keep the work proportional to what
changed, each cheap to check:

1. **Buffer dirty flags** — set by colour updates, cleared on flush. A colour
   store also follows its geometry's draw order, from a kept store when one
   holds the current colours.
2. **Frustum cache** — an exact element-wise comparison of the model-view-projection
   matrix, the dimension and the bounds. An idle camera re-extracts nothing.
3. **Draw ranges** — a camera move collects the leaf runs inside the frustum
   (a walk of the tree's nodes, not of its points) and writes at most 1,024
   ranges (256 without `WEBGL_multi_draw`); a level change with the same runs
   only rewrites the ranges. Nothing
   is uploaded: the ranges are offsets into the stored points.
   `gl.getError()` is a synchronous round trip that waits for every queued
   command, the previous frame's point draws included, so it runs only where
   a failure can appear: after an allocation (`_allocateGpu`, which retries
   only `OUT_OF_MEMORY`, once, after releasing memory kept for speed).
   Writes into storage WebGL has already accepted (a filter tick's alpha
   rows, a recolour, a kept view's alpha, a draw order's rows) fail
   synchronously or not at all and only check `isContextLost()`
   (`requireLiveContext`); a filter tick at 10M points in four views spent
   40 of its 75 ms waiting there. Highlight buffers check on the first
   allocation and when the store grows past the size already verified.

A culled frame draws a subsequence of the stored points in storage order, so
it blends translucent points exactly like an unculled frame, which
`tests/browser/frustum-culling-parity.spec.mjs` checks pixel for pixel.

Binding the alpha and permutation textures is two binds and a few uniform
writes; it never uploads.

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
  sample under MSAA. There, `full` and `light` each have a no-depth variant
  (`HP_FS_FULL_NO_DEPTH`, `HP_FS_LIGHT_NO_DEPTH`, picked by
  `_pointProgramName()` for dimension ≤ `NO_DEPTH_MAX_DIMENSION` = 2, which the
  figure export reads too) that writes a zero-alpha fragment where the depth
  variant discards, which leaves the pixels byte-identical
  (`tests/browser/point-shader-no-depth-parity.spec.mjs`). At 1,000,000 cells
  in 2D with antialiasing, `light` went from 30.7 to 17.0 ms and `full` from
  33.8 to 30.6 ms with the first of these.
- **What 0.9.2 moved out of the per-vertex and per-fragment work anyway.** None
  of it is measurable on the machine above, and all of it is work a smaller GPU
  pays for: the model-view product, the point-size projection factor
  (`u_pointScale`, which replaced `u_viewportHeight` and `u_fov`) and the Blinn
  half vector (`u_halfDir`) are formed once per draw by `writePointShaderFrame`;
  fog is evaluated once per point (below) instead of an `exp` per fragment; the
  specular `pow` is five squarings; the point fragment shaders are `mediump`.
  Depth is formed by the same two products in every quality, so the same point
  wins a near-tie whichever quality is selected.
- **The disc edge in 1D and 2D views is antialiased.** MSAA cannot smooth a
  point sprite's edge: the fragment shader runs once per pixel, and a discard
  or zero alpha removes every sample of it. The no-depth variants therefore
  compute the disc's coverage of the pixel — a one-pixel ramp across the rim,
  narrowing below a radius of four pixels so the smallest points stay crisp —
  from the sprite radius the vertex shader passes (`v_radiusPx`; screen-space
  derivatives are undefined along a sprite's edge). Depth passes keep the hard
  edge: a partially covered rim fragment would still write depth, and in a
  dense cloud every disc grew a ring of background colour where later points
  behind it were rejected.
- **Fog is one model, evaluated per vertex.** Points, highlight rings,
  connectivity edges, the orbit anchor and the velocity overlay splice the same
  GLSL (`FOG_GLSL`, `cellucid/assets/js/rendering/shaders/fog-shaders.js`) and
  read the same two uniforms, `u_fogNear` and `u_fogScale`
  (`log2(e) · density² · 0.6 / span`, derived per pane by `writeFogUniforms`;
  zero when fog is off). Before, each layer restated the law inline with its
  own terms (the anchor by eye depth, the velocity overlay squared, edges with
  their own alpha floor), so the layers of one scene faded at different rates.
  The transmittance is an `exp2` evaluated once per point or edge endpoint and
  passed on as a varying. Fog **thins** a translucent layer —
  `fogOpacity(T) = T·(0.1 + 0.9·T)` — instead of mixing its colour toward a
  constant `u_fogColor`: over a backdrop of the fog colour that is exactly the
  picture the colour-mixing law drew, and over anything else (a grid, a
  vignette, the points further back) a fogged point fades into what is really
  behind it instead of being tinted with a colour that is not. The figure
  export's vector strategies apply the same opacity, so a circle fades into the
  paper. With fog off a point is drawn at full opacity, and that measured
  cheaper than before the change for dense views (10,000,000 points, device
  pixel ratio 2: 2D 99.7 → 92.4 ms, 3D 87.8 → 80.4 ms a frame); with fog on the
  frame is unchanged. A compile-time variant for fog off would remove one
  `exp2` per vertex and was not added. The
  fog parameters are two `uniform1f` calls per program per pane; a shared
  uniform block would save microseconds a frame and add a buffer binding to
  every program, so it was not adopted.
- **Every vertex shader pins its inputs, one to location 0**
  (`tests/vertex-attribute-location-contract.test.mjs`). A draw whose attribute
  0 is not an enabled array makes Firefox on desktop OpenGL emulate it per
  vertex; with locations left to the linker, the instanced edges drew that way
  whenever no earlier pass had left location 0 enabled.
- **Connectivity edges are a 4-vertex triangle strip per edge**, not two
  6-vertex triangles: every vertex projects both endpoints, so four vertices
  are a third less vertex work than six. The endpoints and strength are
  per-instance attributes gathered for the shown edges only, so no vertex
  fetches a texture.
- **The smoke ray march does less arithmetic per sample**: the density-texture
  coordinate map and the volume normalisation are formed once per pixel, the
  three multiple-scattering octaves share one `exp` (the other two are its 4th
  and 16th powers), and the blue-noise jitter is fetched once. Same pixels;
  17.2 against 14.5 ms a frame on the 562,000-cell Suo sample (median of three
  interleaved runs, Apple M1 Pro, device pixel ratio 2).
- **A different primitive.** `GL_POINTS` beat instanced quads and both triangle
  arrangements by 1.41–2.12× at one sample and by up to 1.90× with MSAA enabled
  at 4 samples, with the quad arms' rasterised and shaded fragment counts matched to
  five decimal places. Halving the primitive count with one triangle per point
  was the *worst* arm.
- **Order-independent transparency, a depth prepass, storage reordering, and the
  vertex-ALU hoist.** All measured; none won.

What did move: **antialiasing**, worth about a fifth of the frame at point size
0.75 and a third at large sizes in 3D — and about twice the GPU time in a
dense 2D view, where nothing is depth-rejected (orbiting 1,000,000 cells: GPU
timer totals 25.5 ms on against 11.6 ms off with `full`, 17.5 against 9.9 ms
with `ultralight`; `scripts/render-scenarios.mjs --antialias on|off`). It is a live user setting (`#hp-antialias`) and
no longer a context-creation attribute: `scene-msaa-target.js` owns a
multisampled renderbuffer resolved by a blit each frame, so a change applies to
the next frame with no reload. The stored preference is three-state — `auto`
(the default, and what an absent key means), `on`, `off` — and `auto` resolves
from the dataset's cell count and the active view's dimension — against
`AUTOMATIC_ANTIALIAS_CELL_LIMIT = 5_000_000` in 3D and
`AUTOMATIC_ANTIALIAS_CELL_LIMIT_FLAT = 200_000` in 1D/2D, where nothing is
depth-rejected — re-decided on every dataset publication and on every
`view:representation-changed` of the active view.

**Shader quality selects shading, never geometry.** All three point vertex
shaders — and the copy `HP_VS_HIGHLIGHT` keeps of the formula — clamp
`gl_PointSize` to the same `[0.5, 128.0]`. A higher ceiling in the lighter
shaders once made them draw *larger* sprites than `full`, so the "faster"
quality was up to 1.93× slower at close range. `tests/point-size-clamp-contract.test.mjs`
pins the three bounds together; keep them identical in any fork.

**`sculpted`, the fourth quality, is opaque by design.** Its fragment
shaders (`HP_FS_SCULPTED`, `HP_FS_SCULPTED_NO_DEPTH`) light each sprite as a
sphere in linear light (the colour squared, encoded back with a square root,
so the round trip is exact) and multiply by `u_invCenterShade`, the
reciprocal of the lighting at the sphere's centre for the scene light, so the
centre shows exactly the legend colour (`sculptedCenterShade` mirrors the
shader's constants). The gloss fades in between sprite radii of two and six
pixels from the vertex shader's `v_radiusPx`, never from screen-space
derivatives, which are undefined along a sprite's edge. In a 3D pass the
points are opaque: `enterOpaquePointBlending` replaces the colour
(`blendFuncSeparate(ONE, ZERO, ONE, ONE_MINUS_SRC_ALPHA)`) and, when the pass
draws into the multisampled scene target (`renderParams.multisampled`, from
`sceneMsaaTarget.isMultisampledThisFrame()`), enables
`SAMPLE_ALPHA_TO_COVERAGE`, so the antialiased rim becomes a per-sample mask
that keeps its own depth — the only way to soften an edge with depth writes
on without the background showing through a nearer point.
`restorePointDrawBaseline(..., true)` undoes both after the draw
(`tests/browser/sculpted-point-pass-state.spec.mjs`). Single-sampled, the
coverage floor of 0.5 gives the same hard edge as every other quality. Because
the points are opaque, fog cannot thin them; they alone still mix toward
`u_fogColor`, which the figure export sets to its paper. In 1D and 2D views
they blend like the other qualities. One list names the qualities
(`cellucid/assets/js/rendering/shader-quality.js`: `SHADER_QUALITIES`,
`SHADER_QUALITY_TRAITS`), and one table their highlight styles
(`HIGHLIGHT_STYLES_BY_QUALITY` in `high-perf-shaders.js`); `sculpted` draws
`full`'s torus ring with `animated: false`.

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

### Symptom: “Switching a view's embedding breaks rendering”

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
