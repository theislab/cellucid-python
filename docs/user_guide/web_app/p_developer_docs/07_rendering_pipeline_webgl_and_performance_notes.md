# Rendering pipeline (WebGL) and performance notes

This page explains how Cellucid renders millions of points efficiently, and what changes are likely to cause performance regressions or WebGL failures.

The decisions behind the renderer, the numbers that justify them and what was tried
and dropped are in the decision record {doc}`28_decisions_rendering_layers`; this page
says what the code does now.

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
    - Spatial indexing is dimension-aware (binary tree / quadtree / octree), used for
      picking/LOD. Trees are built in a worker (see "Trees are built off the main
      thread" below).
    - Positions and colors are **two separate vertex streams** — see
      “Point vertex streams” below.
- **Smoke / volumetric density**
  - Backend: `SmokeRenderer` (`cellucid/assets/js/rendering/smoke-cloud/`)
  - `viewer.buildSmokeVolumeGPU(...)` builds the visible-point volume through
    one WebGL2 transaction: bounded position + color batches → R32F density
    splat atlas (eight trilinear corners per cell) and RGBA32F color splat
    atlas (one texel per cell) → logarithmic GPU maximum reduction of the
    density → normalized RGBA8 atlas → RGBA8 3D texture with mipmaps. Red is
    the normalized density; green, blue, and alpha are the cells' mean display
    color premultiplied by that density, so filtering and mipmapping average
    color by the cells that produced it (`density-shaders.js`).
  - A flat cloud — a 1D or 2D embedding, or a planar 3D one, thinner than two
    voxels along an axis — is spread along that axis into a soft slab of fixed
    half-thickness (`SMOKE_SLAB_HALF_THICKNESS`, about four voxels of the 128³
    grid) by both splat passes, so a ray never steps over a voxel-thin sheet.
  - From each published volume, `buildDerivedVolumesGPU(...)` derives two more
    for the fixed key-light direction: a lighting volume (RG16F, half the grid:
    r = optical depth toward the key light, g = sky occlusion) that replaces a
    per-sample shadow ray, and an occupancy volume (R8, one texel per 4³ block
    of voxels: the largest density a sample in the block can read) that lets
    the ray march leap over empty blocks. Both are rebuilt only when the
    density volume or the light direction changes; the live density and
    absorption scale them in the shader.
  - The product path performs no `readPixels` synchronization. It publishes the
    new texture and bounds before completing its notification, preserves the
    prior volume on failure, and returns only a frozen non-owning summary.
  - The volume is a function of visibility, positions, and colors. The viewer
    reports every live position or color republication (a field, colormap,
    category color or log-scale change, a dimension switch, a new dataset)
    through `setSmokeSourceChangeHandler`, and the render controls mark the
    volume dirty exactly as they do for a visibility change.
  - The exact UI grid inventory is 32³, 48³, 64³, 96³, and 128³. The builder
    accepts exact integer sizes from 8 through 128 and rejects values outside
    that range before allocation. At 128³, its 2D atlas is 1536×1408—below
    WebGL2's required 2048 minimum `MAX_TEXTURE_SIZE`—and the published RGBA8
    volume is 8 MiB. The transient float atlases peak below 50 MiB.
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
  - A build reads the GL error state once, after its allocations (each read
    is a synchronous round trip that waits for every queued command). A
    rebuild asked for inputs identical to the shown volume's (the same
    positions, colors, visibility and outlier-quantile arrays, held by weak
    reference, the colors and the visibility at the same write counts
    `getSmokeDensitySource` reports — every write to them goes through
    `per-cell-array-ownership.js`, which counts it — and the same threshold,
    grid size and gamma) does no work. Nothing per cell is copied for it.
    It requires `EXT_color_buffer_float`,
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
    context loss; a restored context gets a new `SmokeRenderer` with the same
    settings, and the volume is built again from its source.
  - The user notification duration covers the synchronous renderer call,
    including validation and publication. The `[GPU Splat] … submission`
    console duration starts after the input checks and measures command
    submission.
    Neither duration is a GPU-completion measurement.

- **Connectivity edges**
  - Backend: `rendering/connectivity-edge-layer.js` (which edges a pane shows:
    `rendering/connectivity-edge-selection.js`).
  - The graph stays on the CPU as the reader worker built it
    (`ConnectivityEdges`, `data/connectivity-edges.js`): its forward
    adjacency — each edge once, at its smaller cell, with its partner and its
    `float32` relative strength, and the forward offsets — 8 bytes an edge
    and 4 a cell. Nothing is stored in
    draw order and nothing is copied on the page.
  - A pane shows the first K of its visible edges — both endpoints pass the
    pane's filters and, with level of detail on, its current level — in one fixed
    random order: a keyed Feistel bijection of the edge indices (each round a
    table of `sqrt(E)` entries), so every K-prefix is a uniform sample, the same
    in every session for the same graph.
    K is the edge slider's value: 250,000 (or every edge, if fewer) for a newly
    loaded graph (`CONNECTIVITY_DEFAULT_VISIBLE_TARGET`), the same in every
    pane. A pane added later selects on its first frame.
  - The edge of rank `r` is `order.at(r)`; its destination and strength are
    read by edge, and its source is found by binary search within one block of
    128 edges (`blockSources`, the source of every 128th edge). The records of
    the first 4K ranks are kept in rank order (the shown-order prefix, 48 bytes
    a shown edge, made when K is set and kept across filter changes); records
    are written in edge order (each chunk's edges sorted, their sources found
    in one forward walk), so K records cost K order steps and one pass.
  - A selection takes the cheaper of two exact walks: reading ranks from 0
    until K are visible (the prefix's records sequentially, past it each rank's
    edge, a hidden destination costing no source search; O(K / visible
    fraction)), or ranking only the visible cells' forward edges and keeping
    the K smallest ranks in at most 2K slots (O(cells + visible cells'
    edges), O(K) memory). Both give identical records; the cell walk is the one that costs time when
    a filter hides most cells ({doc}`28_decisions_rendering_layers`, "Measurements that justify it").
  - Only those K edges reach the GPU: their 12-byte records, from which a
    transform-feedback pass writes one 28-byte instance record each (both
    endpoints' positions, missing axes 0, and the strength) by copying bits
    from the pane geometry's position texture (one texel a cell). The draw
    reads the positions as attributes and tests nothing: on Metal the same
    vertex arithmetic on texture-fetched positions rasterizes some lines
    differently, so the draw keeps attribute inputs and every frame matches
    a CPU gather pixel for pixel. GPU memory is K × 40 bytes a pane plus one
    position texture a geometry; the CPU never gathers positions; a camera
    move uploads nothing, a filter or LOD change re-selects, uploads K records
    and gathers, and a representation switch only gathers.
  - The visible-edge count for the slider's range is exact at once when the
    selection ranked the visible cells' edges; otherwise it is counted over the
    visible cells' forward edges in time-boxed slices after the selection, so a
    filter change never waits for it.
  - The k-NN tool does not need the page's graph: the graph is checked once
    (its edges built, the check recorded in `app/ui/modules/connectivity-controls.js`, the edges let go
    unless shown), and the tool reads the rows of each degree's frontier from
    the file in the reader worker
    (`getConnectivityRows`, 65,536 cells a read, the next read issued while one
    is walked; `KnnExpansion` in `selection-gesture-queries.js`). While a view
    shows the graph, so the page holds its edges, a frontier is answered from
    them instead (`forwardEdgeRows`, two passes over every edge) whenever that
    is measured faster than reading its rows from the file (`knnRowReader` in
    `app/ui/modules/connectivity-controls.js`). A released
    drag commits its step once its rows have arrived, after every step released
    before it.

- **Highlight tools**
  - Backend: `HighlightTools` (`cellucid/assets/js/rendering/highlight-renderer.js`)
  - Handles interactive selection and highlight overlay rendering: each
    highlighted cell a dot of `HIGHLIGHT_COLOR` inside a ring of
    `highlightRingColorFor` it (black), the same in figure export
    ({doc}`28_decisions_rendering_layers`).

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
    `resolveFrame` does nothing: the frame costs nothing extra, and nothing is
    allocated until it is switched on.

:::{warning}
Two caller contracts follow from this and are easy to break:

- Anything drawing into the scene must target `getSceneFramebuffer()` rather
  than assuming framebuffer zero. The smoke renderer and the velocity overlay
  both bind a framebuffer of their own mid-frame and must come back to the
  scene, not to the default.
- Anything reading pixels back must read **after** `resolveFrame`, from the
  default framebuffer, or read the resolved copy.
:::


### Point vertex streams

The rules of the point stores, one heading each. The mechanism behind each, what it
costs, what was measured and what was tried and dropped are in the decision record
{doc}`28_decisions_rendering_layers`, sections "The point stores in detail",
"Levels of detail, the draw order and the order's adoption" and
"The frame, shaders and sprites".

Full detail is drawn from two `ARRAY_BUFFER`s bound into one VAO:

| Attribute | Stream | Layout | Bytes per point | Usage |
| --- | --- | --- | --- | --- |
| 0 `a_position` | `buffers.positions` | Float32, one component per embedding axis (x, xy or xyz), tightly packed | 4, 8 or 12 | `STATIC_DRAW` |
| 1 `a_color` | `buffers.colors` | RGBA8, normalized | 4 | `DYNAMIC_DRAW` |

- Positions are never padded to XYZ: every view's positions hold exactly as many
  components per cell as the embedding has, and that count is the view's dimension
  level (`rendering/position-layout.js`). The shaders keep a `vec3` `a_position`,
  bound with 1, 2 or 3 components; CPU code that needs a 3D point reads it with
  `readPositionXYZ`. A view changes dimension only by publishing that dimension's
  positions (`updatePositions` / `setViewPositions`).
- Alpha is not taken from the color stream: it lives in a separate R8 texture (one
  texel per point), so a visibility change touches neither stream.

#### Visibility is the alpha texture; a color's alpha byte layers it


A cell is shown or hidden by its view's transparency alone (`u_alphaTex`); a filter
writes the transparency and never the colors. The color's alpha byte is the cell's
layer, bottom to top (`alpha-visibility.js`): `UNREAD_ALPHA_BYTE` (253) for a cell a
gene sample did not read, `BENEATH_ALPHA_BYTE` (254) for a field's no-value gray,
255 for the rest; `colorLayersOf` reads a colors array once. A layered view is
drawn one pass a layer over the same slots and ranges (`u_drawPass` 1 then 2, the
gray first; in 3D the gray pass yields in depth by one drawn point diameter,
`u_yieldDepth`), at every level of detail, for culled ranges and kept views, and a
figure draws the gray cells first too (`drawnSequence`). The cells a gene sample
did not read take a pass of their own first (`u_drawPass` 3,
`rendering/read-cell-layer.js`), and the cells read are drawn larger from their own
slot list (`high-perf/read-cell-index.js`). Picking follows the same layers
(`getCurrentPointLayers`), and figure export draws the same three layers
(`FigurePointLayers`). The look and its decisions: items 1 to 10 of "What it does
now" in {doc}`28_decisions_rendering_layers`.

#### Per-point textures are linear textures


`rendering/linear-texture.js`: a `TEXTURE_2D_ARRAY` of power-of-two-wide rows,
stacked into as many layers of `MAX_TEXTURE_SIZE` rows as the count needs. The alpha
texture, the permutation, the permute pass's source and the velocity overlay's
per-cell textures (`rendering/overlays/shared/packed-texture.js`) all use it; a
shader finds item `i` with `linearTexel`.

#### Unchanged panes keep their pixels

The scene target is not cleared between frames. Every change that can alter what a
pane draws advances a scene generation (`invalidateFrame`); a camera gesture alone
does not (`invalidateCameraFrame`). A grid pane whose generation, camera, projection,
rectangle and navigation mode are unchanged is skipped and keeps what it drew; a
redrawn pane clears its own rectangle. A new target, size or layout draws every pane
afresh.

#### Input lands in the frame that draws it


Work that changes what is drawn (coalesced slider and legend input, camera-path
playback) is queued with `runBeforeNextFrame` (`rendering/frame-phases.js`) and runs
at the top of the viewer's frame; `runAfterNextFrame` runs once the frame is drawn.
A control must not wait for an animation frame of its own: it would run after the
viewer's and show the change a frame late.

#### Draws and buffers are bounded by segments


No draw submits more than 2^24 vertices: every draw whose count grows with the data
(the point passes, the permute pass, the highlight pass, figure export) goes through
`drawPointsInSegments`, ranges of at most 16,777,216 vertices, in submission order.
No buffer holds more than one segment of points: every per-point vertex store is a
`PointStore` (`rendering/high-perf/point-store-segments.js`) of consecutive segments
of at most 16,777,216 slots, one buffer each, drawn through one vertex array a
segment (`SegmentedVertexArray`); ranges are split at segment boundaries
(`splitRangesBySegment`) and the point shader adds the segment's first slot
(`u_slotBase`) to `gl_VertexID` wherever it reads a per-slot texture. The picture is
the one-buffer picture byte for byte (`tests/browser/segmented-point-stores.spec.mjs`
forces 4,096-point segments). Growing a store appends: full segments never move
(`unchangedSegments`).

#### No CPU array is limited by one `ArrayBuffer`


A per-cell array that an engine refuses as an `ArrayBuffer` is made in a
WebAssembly memory (`utils/large-arrays.js`, `allocateLargeArray`; every reader
allocation goes through it via `allocateTypedArray`) and read as any other
`Float32Array`. Such an array crosses threads in pieces of at most 192 MiB
(`PieceSender` → `PieceReceiver`); embeddings and vector fields are read by row
range and made in pieces of 2^24 cells in the reader worker (`embedding-payload.js`).
Past 4 GiB no page can hold one array, and a read that needs one fails naming the
bytes it needed.

#### Levels of detail are fixed counts in priority order


A level of detail is the first `m` cells of the open cells' **priority order**
(`data/cell-priority.js`): rank = the 32-bit hash
`fmix32(position ^ fmix32(seed ^ 0x9e3779b9))` of the cell's original position,
under the file's own seed for a priority-ordered file and the subsample's seed
otherwise (the default seed for a whole file Cellucid did not write). The counts are
one fixed ladder (`high-perf/draw-order.js`, `lodLadder`): 1,000 cells, then 5/4 of
the level before (rounded down), and last every open cell, so a level's cells do not
depend on how many cells are open (100M cells have 53 levels in 100 bands, 10M have
43 in 47, 1,000 or fewer one). Every reduced level is a uniform random sample of the
open cells (`tests/lod-priority-order.test.mjs`,
`tests/browser/lod-priority-sample.spec.mjs`), exactly the uniform subsample of size
`m` with the same seed; the same prefix is the first pass of a progressive analysis
(`app/analysis/shared/progressive-cells.js`). The forced-level slider spans the
ladder of the cells shown. The priority comes from how the dataset was opened
(`cellPriorityOf`, set by `state.setCellPriority` → `viewer.setCellPriority` →
`HighPerfRenderer#setCellPriority` right after the scene is published):

- a priority-ordered file (`prepare()` / `cellucid_prepare()` write cells in
  priority order) opened whole or by a prefix, and a subsample of a cell-level
  method (`uniform`, `stratified`, `spatial`), are open **in** priority order
  (`OPEN_ORDER_PRIORITY`): open cell `i` has rank `i`;
- a `blocks` subsample opens whole runs of neighbors in run order, so its cells are
  ranked one by one under its seed; so is every other share not open in its own order
  (one opened in the file's cell order, a priority-ordered file's share in file
  order, and a step of a share grown from few chunks, its cells in chunk order);
- every other dataset (a file Cellucid did not write, opened whole; a synthetic
  benchmark) is ranked under the default seed.

A ranked cell's rank is always that of its **identity**, its original position
(`openOrder[i]`, else `kept[i]`, else `i`: what `data/cell-identity.js` names it by),
never of where it is open, so the same cell has the same rank in every order, step,
view, revisit and session of a share and its levels, picking and the levels a figure
exports draw the same cells first.

#### Cells open in priority order draw every level at once


Their stores are published in rank order, so a level is a prefix of the store as it
is (`high-perf/lod-levels.js`, `RankLevels`): no tree, no permutation, no gather. A
level is described by its count alone. Cells ranked one by one (`blocks`, a whole
foreign file) have no such prefix: their levels arrive with the tree.

#### One record of a level, one validator


Every consumer of the current level (picking, highlights, connectivity edges, the
velocity overlay's spawn sources, the highlight count, figure export) reads
`LodMembership` (`high-perf/renderer-contracts.js`), published by
`HighPerfRenderer#getCurrentLodMembership` and checked by the one
`requireLodMembership`: `count`, `lodLevel`, `pointCount`, `dimensionLevel`, a
`generationToken`, `admissionLevels` (null, or one byte a cell) and `indices` (or
null). A cell is drawn when
`admissionLevels === null ? cell < count : admissionLevels[cell] <= lodLevel`
(`lodMembershipAdmits`). Consumers check the fields they read, never the set of keys.

#### Draw order


A geometry's points are stored on the GPU in the draw order of the spatial tree over
them (`SpatialIndex#drawOrder`, computed by `high-perf/draw-order.js`): grouped by
**band**, then by draw unit (a tree leaf in 3D, one of 64 coarse nodes in 1D and
2D), then by ascending cell index. A band is a fixed run of ranks (`drawBands`): the
ladder's steps, each split again at every multiple of `DRAW_BAND_POINTS`
(2,097,152). So:

- a reduced LOD level is one `drawArrays` **prefix** of the stored points;
- full detail is `drawArrays` over all of them;
- frustum culling collects the runs of consecutive tree leaves inside the frustum
  and draws each run as one contiguous range per band. Where the ranges exceed the
  budget (1,024 with `WEBGL_multi_draw`, one `multiDrawArraysWEBGL` call; 256
  without it, one `drawArrays` per range) the gaps that hide the fewest points are
  bridged first (`SpatialIndex#writeDrawRanges`). No per-point index is written or
  uploaded after the order is built, and a camera that does not move re-culls
  nothing.

Positions are per geometry, so each geometry carries its order; colors are per view,
and each view's color store is gathered into the order of the geometry it shows.
Alpha stays in source order, read through the store order's R32UI slot map
(`u_orderedSlots`). Picking, highlights, overlays, analysis and export read
source-order CPU arrays; nothing but the point draw reads the stored order. The 2D
overlap order of translucent points follows the spatial tree, and is the same in the
live view, kept views, every LOD level, culled frames and figure export (which draws
in `getViewDrawOrder`); nothing sets a draw order on purpose, so no feature depends
on source order.

#### A store takes its tree's order band by band, without blocking a frame

A store in rank order takes its tree's order one band at a time while it is drawn
(`StoreOrder`, `high-perf/store-order.js`), so every level still holds exactly its own
cells; a frame moves at most `STORE_ORDER_POINTS_PER_FRAME` (2,097,152) points into
order, and culling starts once every band is in order. A store of cells not in
priority order takes its tree's order whole in one frame. A tree built under a
priority the cells no longer have is dropped when it lands
(`builtUnderStalePriority`).

#### Cells added or removed in place

An Apply of a larger or smaller share of the same method and seed
(`app/open-cells-resize.js`) keeps the stores of the cells kept
(`viewer.resizeLiveCells` → `HighPerfRenderer#resizeLiveCells`): added cells are
appended in rank order, removed cells are cut from the end, and nothing else is
uploaded or reordered. A progressive opening names the cells it grows to (`capacity`),
so each step up to it only writes: no allocation, no GPU copy and no `gl.getError`.
After the new tree is adopted the stores and every frame are byte-identical to a
direct open of the same cells (`tests/browser/lod-rank-order-resize.spec.mjs`).

#### Trees are built off the main thread

The tree over a geometry and its draw order are built together in a worker
(`high-perf/spatial-tree-worker.js`, `buildSpatialTreeJob`), one tree at a time. Until
a view's tree lands, the view draws levels by rank when its cells are open in priority
order, else full detail; when it lands, the renderer installs it and asks the viewer
for a frame (`onAsyncTreeReady`). A click on a large pane whose tree has not landed
scans its points (the same exact pick); a hover is told the lookup is being prepared;
projectiles collide nothing until it lands. On Chromium a build with levels is
admitted against the tab's JavaScript heap (`jsHeapSizeLimit`).

#### Cells an embedding does not place are left out of it, not refused

A modality that did not measure every cell (MuData), or any embedding with non-finite
rows, gives those cells NaN positions in every axis: the one mark
(`rendering/unplaced-positions.js`). The viewer accepts a row that is finite in every
axis or NaN in every axis (`assertCellPositions`), and nothing else. Bounds, projections,
picking and lasso leave the cells out by themselves; a spatial tree holds the placed
cells only, and each band of its draw order ends with its unplaced cells. NaN never
goes to the GPU (`positionsForGpu`), and every vertex shader that reads a cell's
position rejects it (`isUnplacedPosition`). They are still cells: filters, analyses and
every other embedding keep them.

#### Level of detail and culling turn on when they are needed


Both start off. Once per opening (a dataset, cells grown or cut, another
representation) the viewer measures three frames of every cell of every presented pane
(`rendering/full-detail-frame-cost.js`), once the publication it measures has
committed. Each frame is two shares, each as bounds the true value lies within: the
page's (from the clock) and the GPU's, read with the GPU timer
(`EXT_disjoint_timer_query_webgl2`, exact) where the browser has it, else with a fence
polled in tasks. A frame lasts the longer of the two shares. When the GPU's share
alone is proven over `INTERACTIVE_FRAME_BUDGET_MS` (33 ms: a UX target; a page slow
on its own thread is not made faster by drawing fewer cells), the render controls turn
both on (level of detail only where the ladder has a reduced level) and measure the
frames drawn with them, once their trees are built and adopted
(`verifyFullDetailFrameCost`). No automatic choice may make a frame slower, and the
verdict is what the bounds prove (`verdictOf`): `faster` (they stay and the
Performance tip says what was measured, `onFullDetailFrameCost`, `render-controls.js`),
`not-faster`, or `unresolved` (the bounds overlap); the last two turn them off again
for this framing and the tip says which. The next framing is measured afresh. The rule
does not branch on the dimension.

Only frames drawn with level of detail off are measured; an opening that starts with
it on, or where a person, a session or a caller of the viewer turns it on, is not
measured, and a late answer is dropped. Every measurement the viewer starts (cells
published, grown or cut, another representation, a restored graphics context, a new
framing) is said first (`onFullDetailFrameCost`, `kind` 'measuring'): the render
controls turn their own automatic choice off, the frames of every cell are measured,
and the choice is turned on and verified again as at the first opening. A person's own
change to either box stands for the page. A caller that sets either itself
(`setAdaptiveLOD`, `setFrustumCulling`) makes that opening's choice: every report and
every such call names its opening (`onFullDetailFrameCost`,
`onPerformanceSettingsChanged`), and the render controls leave the rest of that
opening's reports alone. The benchmark's pass profiler times passes with the same
extension; timer queries cannot nest, so the measurement times no frame while a pass
query is open and a pass that begins while the measurement times the frame is not
timed (`untimedPasses`).

#### A store adopts its tree's order in the frame that first needs it

Whenever a tree exists, in every dimension (a hover tree counts too). The order is an
input of the picture: translucent points overlap in it, so the frame after the
adoption is not the frame before it. `isPresentedGeometrySettled` is therefore false
while a tree is being built for a presented view's positions, or lies published but
not yet adopted, with culling and level of detail off too. The adoption uploads the
permutation (4 bytes a point), then rewrites the store's positions and the colors of
every view that draws it in place on the GPU (`GpuPermuter`,
`high-perf/gpu-permute.js`), bit for bit; above 16,777,216 points
(`GATHER_CHUNK_POINTS`) it takes the source one chunk at a time. Every vertex
attribute is released onto the gather's own anchor buffer, never onto no buffer
(`releaseVertexAttribute`, `rendering/gl-utils.js`). An adoption the GPU refuses says
so once, the view draws at full detail, and the store tries again when LOD or culling
is toggled. A kept view's color store is RGBA8, 4 bytes a point, written one word a
point; colors are copied as 32-bit words, so color arrays handed to the renderer must
start on a 4-byte boundary. A recolor rewrites the color store in place (one upload of
4 bytes a point and one GPU pass); only a shared store, or one in another order, is
replaced.

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
    color stores kept for other orders. Positions and the alpha texture are not
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
  - It replaces only the positions: colors, transparency and the alpha
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
- **any number** of snapshot views. Each owns an RGBA8 color store and an R8
  alpha texture (5 bytes a cell; positions are shared by geometry). Memory,
  not a count, limits them: when an allocation is refused, the stores of kept
  views not presented this frame are released first
  (`HighPerfRenderer._releaseDormantSnapshotStores`; the viewer reports the
  presented views every frame with `setPresentedViews`), and a released view
  is rebuilt from its published colors and transparency (both kept by
  reference) when it is next drawn or edited (`_restoreSnapshotStores`),
  drawing exactly what it drew. Only a set of views shown at once that does
  not fit is refused, by the device. A kept view holds no CPU copy of its R8
  image: as for the live view, a transparency change is compared with
  the published values and only the rows whose bytes moved are sent, through
  one bounded 1 MiB scratch; `null` transparency shows every cell.

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
  layout mode, the focused view, and the view inventory. It is memoized: it is
  recomputed only when one of those inputs changes, together with one
  projection matrix per pane, and is otherwise one frozen object.
  `viewer.getPresentedLayout()` returns it and
  `viewer.onPresentedLayoutChange(listener)` delivers each new one.
- **Field of view**: set with each presented layout from the presented area's
  aspect (`rendering/viewer/field-of-view.js`). It is 45° vertically on a
  landscape or square area; on a portrait one the vertical angle widens so
  that 45° spans the width. A single view takes the presented area's value and
  each pane of a grid the value of its own shape (`paneFieldOfView`); the
  pane's point-size projection factor, centroids, picking and figure export read
  that pane's value, so a camera's distance frames the same data on every screen
  shape (`tests/browser/camera-framing-aspect.spec.mjs`).
- **Grid** (`Grid compare` with two or more views): the presented area is split
  into `cols × rows` cells of identical integer framebuffer size, starting at
  the first framebuffer column right of the sidebar. No cell is under the
  sidebar, every cell has the same area, and each renders with a symmetric
  projection and its own scissor rectangle.
- **Single view**: the pane keeps the whole canvas (the scene continues under
  the sidebar), and its projection is shifted horizontally so the principal
  point - the camera target in Orbit and Planar, the look direction in
  Free-fly - sits at the center of the presented area, in every mode, so
  switching modes does not shift the picture.

Each pane carries its WebGL viewport (`viewport*`, framebuffer pixels,
bottom-left origin), its render rectangle (`css*`) and presented rectangle
(`visibleCss*`) in canvas-local CSS pixels, and `projectionCenterNdcX`. Every
mapping between the screen and a view reads these panes: the render loop's
viewports, scissors and per-pane projections (so frustum culling and LOD see
the real pane shape), `getViewportInfoAtScreen()` and therefore click/tap
focus, wheel zoom (which focuses the pane under the cursor), picking (nothing
outside a pane's presented rectangle can be picked), hover delivery to a
Python or R viewer, every highlight tool (lasso polygons are clipped to the pane's
presented rectangle; proximity and KNN indicators are clipped like the WebGL
scissor), projectile aiming and the charge ring, centroid labels, view title
chips, overlays, `getPresentedViewState()` (figure export and its preview),
and the selection bar, which is centered in `presentedArea`.

### Layout arithmetic, and why it changes the cost

```
cols = getPresentedGridColumns(viewCount, presentedWidth, presentedHeight)
rows = ceil(viewCount / cols)
paneWidth = floor((canvasWidth - presentedLeft) / cols)
paneHeight = floor(canvasHeight / rows)
```

`getPresentedGridColumns()` picks the arrangement whose panes have the largest
shorter side (ties go to the fewest empty cells, then to the wider arrangement),
so three panes stand side by side on a wide screen and one above the other on a
phone held upright. On a landscape area: 2 views → 2×1, 3 → 3×1, 4 → 2×2.
Figure export arranges multi-panel figures with the same column rule.

`paneHeight`, with the pane's field of view, forms the point-size projection
factor (`u_pointScale`), so perspective point size is proportional to it. Point
**area** is therefore proportional to `1/rows²`, and total fill across all panes
is proportional to `viewCount / rows²`:

| Views (landscape area) | cols × rows | Relative total fill |
|---|---|---|
| 1 | 1 × 1 | 1 |
| 2 | 2 × 1 | 2 |
| 3 | 3 × 1 | 3 |
| 4 | 2 × 2 | 1 |

This is the **fill** part of the cost only. Every view also pays for every point
it draws, whatever its size, and at the default point size on a large dataset
that per-point part dominates, so views add up (measurements in
{doc}`28_decisions_rendering_layers`, "Measurements that justify it"). The fill
ordering above matters only when points are large enough to overlap heavily, and
only while size attenuation is non-zero; at zero attenuation point size ignores
the viewport. Level of detail is the lever for multiview.

:::{warning}
Do not “optimize” the pane height of the point-size scale to the canvas height for
consistency across panes. Feeding the full canvas height into a 2×2 grid would
quadruple fill cost and make four views genuinely four times a single view.
:::

---|---|---|
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
Do not “optimize” the viewport-height uniform to the canvas height for
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
- Particles spawn from a table of at most `spawnTableSize` (65,536) visible
  cells. It is filled by walking the candidates (the level's cells, or every
  cell) in a fixed order: start, then a step coprime with their count near
  its golden ratio, so any run of them spreads over the cells. Each respawn
  epoch (one particle lifetime) refills the table from where the last one
  stopped. Over successive epochs every visible cell is a spawn source, each
  once a cycle, at 10M cells as at 100k, not a fixed 65,536 of them. A refill
  reads only as many cells as it takes and keeps the visibility texture; a
  table that already holds every visible cell is never refilled. A change of
  inputs (visibility, level, positions) starts a full build, which checks
  every transparency value first.

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
  ├─ dirty-flag gate ── colors dirty or in another order? ── no ──┐
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

The frame gate decides whether a frame is drawn at all. Once the scene has settled
the loop stops: no animation frame is requested until something asks for one
(`ensureFrameScheduled` from an invalidation, a held navigation key, a canvas
resize or queued frame-phase work). A frame that is drawn re-arms the loop, which
runs its settle frames and any time-driven effect, then stops again. Any mutating
viewer or renderer call (except a republication of what is already drawn), any
viewer input that changes something, and any canvas resize marks the next frame for
drawing; hovering or typing into a sidebar field draws nothing. Held navigation keys,
inertia, smooth zoom, a grid fade, smoke while it animates, an overlay while it draws
something that moves, visible highlights in a quality whose highlight style animates
(only `full`), live projectiles, a camera motion still waiting for its
`camera-settled` event and an adaptive LOD still stepping toward its level
(`isViewLodConverging`) keep frames drawing. A frame that only an ambient animation
asks for is paced by `rendering/ambient-frame-pacing.js`: it waits until the main
thread has been free as long as the last frames held it and, once the page's frames
run a display frame late or more while the GPU is still running the last frame, waits
as long again; a GPU that keeps up is never held back. While the highlight pulses,
each pane's base layer is kept in a framebuffer (`rendering/pulse-layer-cache.js`)
and a pulse-only frame draws only the layer over it; a multisampled frame resets and
draws again only the samples the rings cover (`rendering/pulse-sample-mask.js`).
After the last cause the viewer draws two more frames (`FRAME_SETTLE_MIN_FRAMES`) and
stops; every multi-frame convergence is a cause of its own, and there is no
wall-clock window on top. A frame that throws marks the next frame for drawing. The
drawing buffer is preserved, so the canvas keeps showing, and reading back, the last
drawn frame.

Inside a drawn frame, three independent gates keep the work proportional to what
changed, each cheap to check:

1. **Buffer dirty flags**: set by color updates, cleared on flush. A color store also
   follows its geometry's draw order, from a kept store when one holds the current
   colors.
2. **Frustum cache**: an exact element-wise comparison of the model-view-projection
   matrix, the dimension and the bounds. An idle camera re-extracts nothing.
3. **Draw ranges**: a camera move collects the leaf runs inside the frustum (a walk
   of the tree's nodes, not of its points) and writes at most 1,024 ranges (256
   without `WEBGL_multi_draw`); a level change with the same runs only rewrites the
   ranges. Nothing is uploaded: the ranges are offsets into the stored points.

`gl.getError()` and `gl.isTexture` are synchronous round trips that wait for every
queued command, so they run only where a failure can appear: after an allocation
(`_allocateGpu`, which retries only `OUT_OF_MEMORY`, once, after releasing memory kept
for speed). Writes into storage WebGL has already accepted (alpha rows, a recolor, a
kept view's alpha, a draw order's rows) only check `isContextLost()`
(`requireLiveContext`); highlight buffers, centroid buffers and the velocity
overlay's tables are rewritten in place when the new data fits and read the error
queue only for new storage. No owner may leave an error for the next owner to read: a
deleted object is never bound again (`deleteProgramReleasingUse`,
`deleteTextureTracked` and `restorableTexture` in `gl-utils.js`).

A culled frame draws a subsequence of the stored points in storage order, so it
blends translucent points exactly like an unculled frame, which
`tests/browser/frustum-culling-parity.spec.mjs` checks pixel for pixel. Binding the
alpha and permutation textures is two binds and a few uniform writes; it never
uploads.

The practical rule: **if you add work to the draw path, it must sit behind a gate
of its own.** “It is only a few microseconds” is how the idle frame stops being
free.

### Shader and draw rules

The mechanism of each, the numbers and what was rejected are in
{doc}`28_decisions_rendering_layers`, section "The frame, shaders and sprites".

- **Shader quality selects shading, never geometry.** Every vertex shader that sizes a
  cell's sprite calls the one GLSL function `pointSpriteSize(eyeDepth)`
  (`POINT_SPRITE_GLSL`, `rendering/point-sprite.js`): the perspective size, clamped to
  `[0.5, u_pointSizeMax]`. `u_pointSizeMax` is the viewer's cap, `POINT_SPRITE_MAX_PX`
  (160 px), or the GPU's `ALIASED_POINT_SIZE_RANGE[1]` where smaller (`pointSpriteCap`);
  picking (`drawnPointRadiusCss`), the read-cell layer's density estimate and the
  figure export's vector circles use the JavaScript `pointSpriteSize`, which
  `tests/point-sprite-contract.test.mjs` proves equal to the GLSL.
- **Nothing fades or shrinks by distance.** A cell near the camera is drawn whole,
  opaque, at its perspective size up to the cap, and is picked like any other
  (`tests/browser/point-near-camera-solid.spec.mjs`).
- **The point-size curve lives in one module**,
  `cellucid/assets/js/rendering/point-size-scale.js` (the slider curve, the viewer's
  bound, `pointSizeSliderPositionForCellCount`). Do not reintroduce a second copy of its
  `0.25` anchor: the slider would throw at its own low end.
- **2D views use no-depth point shaders** (dimension ≤ `NO_DEPTH_MAX_DIMENSION` = 2,
  chosen by `_pointProgramName()`, which the figure export reads too;
  `HP_FS_FULL_NO_DEPTH`, `HP_FS_LIGHT_NO_DEPTH`), pixel-identical to the depth
  variants (`tests/browser/point-shader-no-depth-parity.spec.mjs`) and with an
  antialiased disc edge (`v_radiusPx`).
- **Fog is one model, evaluated per vertex** (`FOG_GLSL`,
  `cellucid/assets/js/rendering/shaders/fog-shaders.js`; uniforms `u_fogNear`,
  `u_fogScale` from `writeFogUniforms`). Fog thins a translucent layer,
  `fogOpacity(T) = T·(0.1 + 0.9·T)`, instead of mixing toward `u_fogColor`.
- **Every vertex shader pins its inputs, one to location 0**
  (`tests/vertex-attribute-location-contract.test.mjs`).
- **Connectivity edges are a 4-vertex triangle strip per edge**, with per-instance
  attributes, so no vertex fetches a texture.
- **`sculpted`, the fourth quality, is opaque by design**: sphere lighting
  (`HP_FS_SCULPTED`, `HP_FS_SCULPTED_NO_DEPTH`), `enterOpaquePointBlending` and
  `SAMPLE_ALPHA_TO_COVERAGE` in a multisampled 3D pass, undone by
  `restorePointDrawBaseline(..., true)`
  (`tests/browser/sculpted-point-pass-state.spec.mjs`). One list names the qualities
  (`cellucid/assets/js/rendering/shader-quality.js`: `SHADER_QUALITIES`,
  `SHADER_QUALITY_TRAITS`), and one table their highlight styles
  (`HIGHLIGHT_STYLES_BY_QUALITY` in `high-perf-shaders.js`).
- **Antialiasing is a live user setting** (`#hp-antialias`), not a context attribute:
  `scene-msaa-target.js` resolves a multisampled renderbuffer by a blit each frame. The
  stored preference is `auto` (the default, and what an absent key means), `on` or
  `off`; `auto` resolves from the cell count and the active view's dimension, against
  `AUTOMATIC_ANTIALIAS_CELL_LIMIT = 5_000_000` in 3D and
  `AUTOMATIC_ANTIALIAS_CELL_LIMIT_FLAT = 200_000` in 1D/2D, on every dataset publication
  and every `view:representation-changed` of the active view.
- **Closed by measurement**, so not to be re-opened without a new one: fragment-shader
  simplification in 3D views, per-vertex arithmetic hoists, a different primitive than
  `GL_POINTS`, order-independent transparency, a depth prepass and storage reordering
  for fill. Per-sprite rasterization is the floor; below `gl_PointSize` ≈ 15 px a frame
  is primitive-bound, above it fill-bound.

---

## Performance footguns (common mistakes)

### 0) Rejecting invisible points

Filtered-out points are not skipped by shrinking them. Setting `gl_PointSize` to
zero does **not** remove a point: the driver’s minimum aliased point size is one
pixel, so the point is still assembled, rasterized and shaded, and only then
discarded in the fragment stage — and the ultralight fragment shader has no
`discard` at all, so a hidden point would also write depth and occlude visible
ones.

The point vertex shaders push hidden vertices outside the clip volume as well as
zeroing the size, which removes them before rasterization. This is what makes a
heavily filtered view cheaper than an unfiltered one.

**Every** vertex shader that draws per-cell geometry does this, including the
highlight-dot and centroid shaders: a category with no visible members has its alpha
set to zero by the field summary, and centroids are drawn with depth writes enabled,
so an invisible sprite would occlude real points.

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

Smoke density never makes a full duplicate position or color array. It
retains the dataset-owned positions, colors, alpha, and optional
outlier-quantile arrays, validates them in place, and streams visible cells
through two transient buffers fixed at 262,144 cells: positions in the source
layout (at most 3 MiB, for a 3D embedding) and RGBA8 colors (1 MiB), 4 MiB per
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
a selection above ~16.7 million cells was not slow — it *threw*. Presence is
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
- the startup failure page “This browser cannot run Cellucid” (the probe in
  `app/startup-failure.js` found no WebGL2), or “The browser’s graphics stopped
  working” when the browser has WebGL2 but refused the viewer a context
  (`contextRefusedError`, `rendering/gl-errors.js`). `createViewer` asks for one
  single-sampled `webgl2` context and throws when it gets none.

Causes:
- old browser
- corporate GPU policy / disabled WebGL
- remote desktop environments
- a graphics process reset, or too many pages drawing at once

Mitigation:
- test with a modern browser first
- the failure page already says what to do and links the requirements

### WebGL context lost

Symptoms:
- a dialog says the graphics were reset by the browser and are being restored
- DevTools shows “WebGL context lost”

Causes:
- a GPU process crash or reset, a driver timeout
- GPU memory pressure (large datasets + high settings)

What happens (`rendering/viewer.js`, "Lost and restored graphics"): the
viewer asks for the context back (`preventDefault` on `webglcontextlost`),
fences every GPU owner without a GL call, and captures the session at once
(`app/graphics-reset-session.js`: the same bundle Save session writes, read
from CPU state), which the dialog offers. When the browser restores the
context, every GPU object is made again from what the CPU holds
(`restoreLostGraphics`): the renderer with its settings (quality, level of
detail, culling, forced level, cell priority), the live cells, each kept view
with its own geometry, centroids and colors, the highlights
(`HighlightTools#resumeRendering`), the connectivity edges (the layer's
`suspend`/`resume`), the grid, studio and centroid programs, the antialiasing
target and smoke settings, the backdrop look with its ambient motion; the
velocity overlay and the smoke volume are drawn again by their owners
(`onGraphicsReset`, the smoke source change). A restored context has every
value of its state back at the initial one, so the viewer sets again what it
set on the context itself (`applyViewerGlState`: depth test, blending and the
background's clear color, without which a black, white or grid background
cleared to transparent black), and the state shadow forgets what it knew and
reads the context afresh (`primeGlState`, `gl-state-shadow.js`).
Filters, selections, kept views and analyses never leave the page. Only if
restoring fails does the dialog offer Reload, and the captured session is
written to disk first. `tests/browser/graphics-context-restore.spec.mjs`
proves it on every engine, on the studio and the black background: a steady
pulse frame before and after submits the same draws, programs and uniforms
(`helpers/gl-submissions.mjs`), and draws the same pixels, compared exactly in every
engine (`helpers/pixel-parity.mjs` allows no tolerance); a backdrop look comes back with
its ambient motion.

Mitigation:
- reduce smoke settings / disable heavy overlays
- lower dataset size for reproduction
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


Next: {doc}`08_ui_modules_map`.
