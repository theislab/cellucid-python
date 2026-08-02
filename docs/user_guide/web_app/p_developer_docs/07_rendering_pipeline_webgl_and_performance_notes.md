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

Hard constraint:
- **WebGL2 only**. If `canvas.getContext('webgl2')` fails, the app throws early.

---

## Major render subsystems (what gets drawn)

`viewer.js` orchestrates multiple render passes/subsystems:

- **Scatter points (main cloud)**
  - Backend: `HighPerfRenderer` (`cellucid/assets/js/rendering/high-perf-renderer.js`)
  - Receives: positions (float32), colors (RGBA uint8 packed), and alpha/transparency.
  - Notes:
    - Spatial indexing is dimension-aware (binary tree / quadtree / octree), used for picking/LOD.

- **Smoke / volumetric density**
  - Backend: `SmokeRenderer` (`cellucid/assets/js/rendering/smoke-cloud/`)
  - `viewer.buildSmokeVolumeGPU(...)` builds the visible-point volume through
    one WebGL2 transaction: bounded position batches → R32F splat atlas →
    logarithmic GPU maximum reduction → normalized R8 atlas → R8 3D texture.
  - The product path performs no `readPixels` synchronization. It publishes the
    new texture and bounds before completing its notification, preserves the
    prior volume on failure, and returns only a frozen non-owning summary.
  - The exact UI grid inventory is 32³, 48³, 64³, 96³, and 128³. The builder
    accepts exact integer sizes from 8 through 128 and rejects values outside
    that range before allocation. At 128³, its 2D atlas is 1536×1408—below
    WebGL2's required 2048 minimum `MAX_TEXTURE_SIZE`.
  - A build requires a clean GL error state, `EXT_color_buffer_float`,
    `EXT_float_blend`, and sufficient 2D/3D texture limits. These are explicit
    preconditions; there is no alternate density implementation.
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

---

## Viewer public API (how state becomes pixels)

The viewer exposes a small API surface used by `DataState` and `main.js`.

### Core dataset load

- `viewer.setData({ positions, colors, outlierQuantiles, transparency, dimensionLevel })`
  - Initializes point buffers and builds spatial indices.
  - `colors` are RGBA uint8 packed; `transparency` is a separate float array used for filtering.

### Incremental updates (critical for performance)

- `viewer.updateColors(colors)`
  - Updates color buffer without reloading positions.

- `viewer.updateTransparency(alphaArray)`
  - Updates alpha/visibility texture/buffer.
  - Also triggers highlight/overlay rebuild hooks (e.g., selection needs the new visibility mask).

- `viewer.updatePositions(positions)`
  - Used for dimension switching (same point count, different embedding).
  - Must preserve the same indices; mismatched lengths are rejected.

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
- caches per-view positions for snapshot views
- uses optimized paths for alpha/transparency sharing so it does not re-upload N full buffers per view

### Layout arithmetic, and why it changes the cost

```
cols = viewCount <= 3 ? viewCount : ceil(sqrt(viewCount))
rows = ceil(viewCount / cols)
paneHeight = floor(canvasHeight / rows)
```

So: 2 views → 1×2, 3 → 1×3, 4 → 2×2.

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

This is why a 2×2 grid costs about what a single view costs while submitting
four times the vertices, and why the row layouts are the expensive ones. It
holds only while size attenuation is non-zero; at zero attenuation point size
ignores the viewport and every view is straightforwardly additive.

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
frame transfers no data to the GPU.**

```
render()  (self-rescheduling rAF; runs every frame, unconditionally)
  │
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

Three independent gates, each cheap to check:

1. **Buffer dirty flags** — set by colour updates, cleared on flush.
2. **Frustum cache** — an exact element-wise comparison of the model-view-projection
   matrix, the dimension and the bounds. An idle camera re-extracts nothing.
3. **Index publication watermark** — a camera move that changes the matrix but
   not the *ordered set of admitted spatial nodes* reuses the published element
   buffer. Even the post-upload error check is gated on the allocation
   watermark, because running it every frame was measured in milliseconds.

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

- **Fragment-shader simplification of any kind.** Fog, lighting, the alpha
  `texelFetch`, the round `discard`, the mere presence of `discard`, and the
  whole `full` → `light` → `ultralight` ladder (43 → 9 → 2 statements in the
  translated Metal) each stayed inside its own twin band at every point size.
- **A different primitive.** `GL_POINTS` beat instanced quads and both triangle
  arrangements by 1.41–2.12× at one sample and by up to 1.90× under the shipped
  4× MSAA, with the quad arms' rasterised and shaded fragment counts matched to
  five decimal places. Halving the primitive count with one triangle per point
  was the *worst* arm.
- **Order-independent transparency, a depth prepass, storage reordering, and the
  vertex-ALU hoist.** All measured; none won.

What did move: **antialiasing**, now a user setting (`#hp-antialias`), worth
about a fifth of the frame at the default point size and a third at large sizes.

**Shader quality selects shading, never geometry.** All three point vertex
shaders — and the copy `HP_VS_HIGHLIGHT` keeps of the formula — clamp
`gl_PointSize` to the same `[0.5, 128.0]`. A higher ceiling in the lighter
shaders once made them draw *larger* sprites than `full`, so the "faster"
quality was up to 1.93× slower at close range. `tests/point-size-clamp-contract.test.mjs`
pins the three bounds together; keep them identical in any fork.

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

Smoke density never makes a full duplicate position array. It retains the
dataset-owned positions, alpha, and optional outlier-quantile arrays, validates
them in place, and streams visible positions through a transient buffer fixed
at 262,144 XYZ points (3 MiB). Large inputs are submitted in bounded batches.
Do not replace this with a buffer sized to the complete dataset or introduce
GPU→CPU volume readback. The grid cap and early release of completed reduction
textures bound the other transient memory.

### 2) Full buffer re-uploads when only alpha changed

Filtering often only needs alpha changes.
Do not rebuild/re-upload positions or full color buffers for filtering.

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
- `updatePositions` received wrong length (`pointCount * 3` required).
- Embedding fetch failed and returned empty/invalid buffer.

How to confirm:
- DevTools → Network: check embedding requests.
- Console: look for `[Viewer] updatePositions: position count mismatch`.

Fix:
- Validate embedding buffers before calling viewer updates.
- Ensure the dimension manager returns consistent point ordering.

---

Next: {doc}`08_ui_modules_map`.
