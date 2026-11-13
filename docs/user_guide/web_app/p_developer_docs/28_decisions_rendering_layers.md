# Decision record: rendering layers, the gene-sample look and levels of detail

This page answers: **"Why does a gene read for a sample of the open cells look
the way it does, why is the gray for 0 an option taken from the theme, how do the alpha
layers and levels of detail fit together, and how are rendering timings
taken?"** It records decisions, their evidence and what was replaced. The
rules of the renderer (buffers, passes, culling, frame cost) are in
{doc}`07_rendering_pipeline_webgl_and_performance_notes`; this page holds what the
code does beyond a rule, the reasons and the numbers. Update this page in the same
change as the code it describes.

## At a glance

**Audience**
- Developers changing point colors, layers, sizes,
  levels of detail, picking, figure export or smoke, or measuring the
  renderer.

**Time**
- 10 minutes.

**Prerequisites**
- {doc}`07_rendering_pipeline_webgl_and_performance_notes`, the sections
  "Point vertex streams" and "Levels of detail are fixed counts in priority
  order".

---

## What it does now

All paths below are under `cellucid/assets/js/` unless they name another
repository.

1. **A color's alpha byte is its layer, never its visibility.** Visibility
   is the view's transparency (an R8 texture). A field's RGBA8 color
   carries the layer in its alpha byte (`rendering/alpha-visibility.js`):
   - 255: colored, on top;
   - `BENEATH_ALPHA_BYTE` = 254: the no-value gray (a missing value, a
     non-positive value on a log scale, and a gene's 0 when gray for zero is
     on), drawn beneath;
   - `UNREAD_ALPHA_BYTE` = 253: a cell the gene sample did not read, drawn
     at the bottom.

   Every consumer states what an alpha byte means:
   `requireColorAlphaMeaning` returns `'layer'` for cells and `'opacity'`
   for markers of their own (centroids). `colorLayersOf` reads a colors
   array once and reports whether it is layered and how many cells are not
   read.
2. **The color kernel writes the layer.** `app/state/managers/continuous-colors.js`
   writes whole RGBA words. No value becomes `MISSING_WORD` (alpha 254), and
   a cell past the read rows becomes `UNREAD_WORD` (alpha 253). Both carry
   the RGB of the gray 211 for any reader that draws no layers. A sparse gene
   costs one word fill plus its stored entries (`data/column-forms.js`).
3. **Passes, bottom first** (`rendering/high-perf-renderer.js`): one pass for
   a single layer; with layers, the gray (pass 1) and then the colored cells
   (pass 2). When some cells are not read, those come first (pass 3):
   - flat (no lighting), at the level's size, without depth writes, in
     `u_unreadColor`;
   - then the gray in `u_zeroColor`;
   - then the colored cells at the read size. In 3D the gray pass yields
     in depth by one point diameter (`u_yieldDepth`).

   Each pass owns its blend state.
4. **Gray for zero is an option, off by default.** A gene's 0 is a value: it
   takes the color scale's lowest color, opaque. The legend's **Gray for
   zero** switch (`DataState.setContinuousZeroGray`, the field's `_zeroGray`,
   defaulted to false by `ensureContinuousMetadata`) gives the kernel a
   floor at 0 (`neutralAtOrBelow: 0` in `updateColorsContinuous`), so a
   gene's zeros take `MISSING_WORD` and are drawn beneath. It applies to a
   gene with no negative value only (`_zeroGrayApplies`; an obs field or
   scaled data has no floor, and the switch is hidden); on a log scale a 0
   has no logarithm and is the gray either way (the switch is disabled). The
   legend model says it as `zeroGray: {enabled, applies, drawn}` (`drawn`:
   the zeros are gray now). It is per field like the log scale: kept views
   clone it (`_cloneFieldState`), the var axis carries it by gene
   (`GENE_INTENT_KEYS`), a duplicated field keeps it, and the switch is the
   action `color.continuous.zero-gray`, one history step. It is a color of
   the field, held by the `colors` state domain (`state/domains/colors.js`:
   `zeroGray: true` in the field's entry when on, absent when off, since a
   value records only what differs from the field's opening state), so a step,
   an undo and a session file hold it. Python and R read a session's `colors`
   value only as data (an apply writes the `cells`, `field-overlays`,
   `highlights` and `gene-sets` values to an object and reports `colors` as
   presentation, not applied), and the shared golden history
   (`history-golden.cellucid-session`) holds no `zeroGray`. Everything that
   draws a colors array follows from the alpha byte alone (views, smoke,
   figures, the cell card). A figure whose zeros are gray says `Gray: 0`
   under its color bar (`buildColorScaleNote`).
5. **The gray for 0 is made from the background** (`rendering/read-cell-layer.js`,
   `zeroCellColor`). It is a neutral gray `ZERO_LIGHTNESS_CONTRAST`
   (15.44 CIE L\*, the step of gray 211 from white) away from the
   background's lightness, toward the middle. On white that is 211, as
   before. On the dark grid it is 56, and on black 38.
6. **The fade** (`unreadCellColor`). A cell not read is drawn in the
   background moved `UNREAD_LIGHTNESS_STEP` (8 L\*) toward the gray for 0,
   but never more than half the way. The map keeps its shape and recedes,
   and a cell not read never looks like a 0.
7. **The read size** (`readCellSizeMultiplier`). Each frame projects up to
   8,192 map cells and as many read cells, and measures the map's on-screen
   area with an occupancy grid. It then picks the diameter at which the
   shown read cells cover `READ_CELLS_COVERED_FRACTION` (40%) of that area.
   The multiplier stays between 1 and √(drawn ÷ read): a sample of every
   cell draws at the normal size.
8. **Read cells at every level of detail** (`rendering/high-perf/read-cell-index.js`).
   When the level draws at least as many cells as were read, the read cells
   are drawn from their own element buffer of slots (4 B per read cell),
   using `drawElements` over the store's own vertex arrays. Nothing is
   copied. When more cells were read than the level draws, the level draws
   its own share of them by rank, and the size multiplier is computed with
   that share (`readShare`).
9. **One set of layer colors for every surface.** `cellLayerColors` and
   `drawnCellRgb` are used by:
   - the views (`setLayerBackground`, and the CSS variables
     `--cell-zero-color` and `--cell-unread-color` that the legend dots use);
   - smoke, whose density splat substitutes the layer colors for bytes 254
     and 253 and rebuilds when the background changes
     (`rendering/smoke-cloud/smoke-density.js`, `rendering/shaders/density-shaders.js`);
   - figure export, which takes its colors from the figure's own paper
     (`FigurePointLayers` in
     `app/ui/modules/figure-export/utils/lod-membership.js`; `drawnSequence`
     orders the layers);
   - the cell info card's swatch (`app/ui/modules/cell-info.js`).

   The legend and the figure share one note, `geneSampleLayersNote`
   (`data/gene-read-plan.js`): "Faded: not read", and "Faded: not read ·
   gray: 0" when the zeros are gray (`zeroGray.drawn`). The legend writes it
   in one line beside a dot of each color, under the color bar, only for a
   gene read for a sample; under it, one line naming the sample links to
   Genes read for (`ui/modules/gene-read-control.js`, `createGeneReadLink`).
10. **Picking at the drawn size** (`rendering/picking.js`, fed by
   `getCurrentPointLayers`). An upper layer wins (in 3D, any read cell over
   a cell not read). Read cells are found at their drawn size and, when they
   are drawn at every level, outside the level too.
11. **A level of detail is the first `m` cells of the priority order**
    (`data/cell-priority.js`, `rendering/high-perf/draw-order.js`,
    `rendering/high-perf/lod-levels.js`). The priority is the subsampling
    rank of a cell's identity (its original position) under the opening's
    seed, so every level is a uniform random sample. Level counts are a fixed
    ladder: 1,000, then ×5/4, then every open cell. That fixed ladder is what
    lets growth append ranks in place. Adaptive selection stops at 2,000,000
    drawn points or an 8× reduction, whichever is smaller. The forced-level
    slider is not bounded. Points grow as √(reduction) × 0.2 + 0.8, up to a
    44× reduction.
12. **Levels of detail and culling turn on by measurement**
    (`rendering/full-detail-frame-cost.js`). They turn on only when the
    GPU's share of a frame of every cell is proven over
    `INTERACTIVE_FRAME_BUDGET_MS` (33 ms), timed exactly by the GPU timer
    where the browser offers one, else bounded by fence polls. They are kept
    only when the bounds prove frames with them faster by more than the
    frames' proven spread; otherwise they go off, and the tip says whether
    they were measured no faster or closer than the browser can time. See
    page 07, "Level of detail and culling turn on when they are needed".

### The point stores in detail

The rules are in {doc}`07_rendering_pipeline_webgl_and_performance_notes`, "Point
vertex streams". What the code does beyond the rule:

- **Positions.** A 2D embedding costs 8 bytes a point of positions instead of 12 on
  both the CPU and the GPU, and is normalized in place on load instead of being
  expanded into a padded copy. WebGL supplies 0 for the axes a 1D or 2D stream does
  not carry, which is exactly the value of an axis the embedding does not have;
  `readPositionXYZ` writes the missing axes as zero without allocating. Because a
  view's dimension is its position layout, it changes dimension only by publishing
  that dimension's positions.
- **Unread cells and the read size.** The unread pass is flat (no lighting), at the
  level's size, without depth writes, so every cell read shows over them in 3D too.
  `unreadCellColor` moves the background toward the gray for 0 by 8 CIE L*, or half the
  way where the two are closer, so the map keeps its shape, recedes, and never looks
  like a 0. `readCellSizeMultiplier` projects up to 8,192 cells of the map and as many
  read cells each frame (under a millisecond at any size), measures the map's area on
  screen as the cells of an occupancy grid that hold a probe, and takes the diameter at
  which the read cells shown would cover 40% of that area; it never goes below 1 (the
  normal size) nor above the square root of the drawn cells over the read ones, so a
  sample of every cell draws at the normal size. The read cells are drawn at every level
  of detail from their own slot list (`high-perf/read-cell-index.js`): an element buffer,
  4 bytes a read cell, drawn with `drawElements` over the store's own vertex arrays, so
  nothing is copied; it is used when no more cells are read than the level draws, and is
  made again when the colors, the store's order or its segments change. A level with
  fewer cells than that draws its own share of the read cells, by range. Picking (an
  upper layer wins; in 3D, any read cell over a cell not read) finds read cells at their
  drawn size and, when drawn at every level, outside the level too. Figure export draws
  the faded color made from the figure's own background in PNG, Hybrid and both vector
  strategies.
- **Linear textures.** A count that fits one row gets one row exactly as wide as the
  count, so its first upload covers the whole image (Firefox clears an image a first
  upload covers only in part). A shader finds item `i` with a shift, a mask and one
  division; the velocity overlay's per-cell textures are the vector field, positions,
  visibility and spawn table.
- **Panes that keep their pixels.** Time-driven effects other than the camera's own
  motion (fades, smoke, the velocity overlay, the highlight pulse, projectiles, LOD still
  converging) advance the scene generation every frame they run. Orbiting one pane of
  four redraws one pane, and a frame that would repeat what a single view already shows
  (the settle frame after a change) draws nothing.
- **Input and the frame order.** Animation-frame callbacks run in the order they were
  registered, and the viewer re-registers its own at the start of every frame; while no
  viewer drives frames, each `runBeforeNextFrame` / `runAfterNextFrame` callback runs
  from an animation frame of its own.
- **Draws.** `drawPointsInSegments` covers the point passes, the permute pass, the
  highlight pass and figure export; a 100M-point frame costs at most 6 extra calls. Draws
  execute in submission order, so the image (and transform-feedback output) is that of
  one draw.
- **Segments.** Full detail, a reduced level's prefix and frustum-culled ranges (one
  `multiDrawArraysWEBGL` a segment) are split at segment boundaries. Highlights drawn from
  the store bind one VAO a segment too, and figure export's own context uploads each of
  its passes (points, centroid overlay, highlight rings) as segmented stores; a store the
  export's GPU refuses fails the export by name instead of drawing nothing. The same
  points are submitted in the same order, so the picture is the one-buffer picture byte
  for byte (`tests/browser/segmented-point-stores.spec.mjs` compares every pixel and every
  stored word); a store of up to 16,777,216 points is one buffer. A refused allocation is
  one segment, and `PointStore`'s layout makes growing a store an append. Per-point
  textures are addressed by signed 32-bit shader integers, so they hold up to 2^31 − 1
  items; the highlight renderer's slot-word buffer (4 bytes a slot) stays one buffer.
- **Large arrays.** Only arrays the engine refuses take a memory, and every consumer reads
  them as any other `Float32Array`: picking, lasso, projectiles, figure export, kept
  views, sessions (the cell-order digest), the velocity overlay, edges. Memories made are
  positions, vector fields, layouts grown or reordered in place, copies for figure export
  and the read-only facades. A memory's buffer cannot be transferred, and cloning a typed
  array copies its whole buffer, so such an array crosses threads in pieces (the AnnData
  reader worker posts any array of a reply it cannot transfer that way, `{event: 'piece'}`
  before the reply; the generator worker its positions and colors; the figure-export
  client its uploads; the tree worker assembles the slices it is sent). Embeddings and
  vector fields are made in pieces of 2^24 cells, each posted as it is made, so the worker
  never holds a whole one and a float64 `obsm` of any size is read. Past 4 GiB (357.9M 3D
  cells) no page can hold one array.
- **Kept views.** A kept view's color store is RGBA8, 4 bytes a point like the live one,
  and is written the same way, one word a point (its alpha byte is never read: a kept
  view's alpha is its R8 texture). Nothing proportional to the point count is gathered on
  the CPU, and no frame is kept alive for it. A recolor rewrites the view's color store in
  place when the view alone owns its store and the store holds the order it draws with. A
  view that switches to another representation keeps the color store it leaves (keyed by
  that order).

### Levels of detail, the draw order and the order's adoption

The rules are in {doc}`07_rendering_pipeline_webgl_and_performance_notes`, "Levels of
detail are fixed counts in priority order" and the sections after it. What the code does
beyond the rule:

- **Levels.** A level's count, and therefore its cells, does not depend on how many cells
  are open. Every reduced level draws the cloud's density and its mix of categories as all
  cells do, up to sampling noise, and where a cell is or what it is never decides whether
  it is drawn: `tests/lod-priority-order.test.mjs` holds every level of a category-sorted
  atlas, culled and not, within five standard errors of the population's shares, and
  `tests/browser/lod-priority-sample.spec.mjs` does the same for the drawn indices in the
  browser. A level of `m` cells is exactly the uniform subsample of size `m` with the same
  seed. A ranked cell's rank is that of its identity, so its levels of detail, picking and
  the levels a figure exports draw the same cells first.
- **Levels of cells in priority order.** A level is its count alone: its cells are the open
  cells `[0, count)`, in whatever order the store holds them, so the same record serves
  every band of the tree's adoption (a consumer that iterates the cells in draw order,
  figure export, reads the store's order, `getViewDrawOrder`). The tree is built in the
  worker meanwhile, and its order is adopted afterwards. `LodMembership`:
  `admissionLevels` null means the level is the first `count` cells (`RankLevels`, and a
  tree over cells in priority order), else one byte a cell (a tree over cells ranked by
  hash); `indices` are the level's cells in a tree's draw order. Picking sizes a dot with
  the level's own size multiplier, as it is drawn (`getCurrentLODSizeMultiplier`).
- **Draw order.** Band boundaries are the same ranks whatever the cell count. A subtree
  wholly inside the frustum is one run, read from its node. When the ranges exceed the
  budget, the gaps between runs that hide the fewest points, counted band by band, are
  bridged first, so the drawn total stays close to the exact cull.
- **Ranked order.** A ranked order is two passes of a 16-bit radix sort over the ranks in
  the tree worker; ranks are a bijection of positions, so a dataset opened whole sorts the
  ranks alone and recovers each cell from its rank. Within a band the units and the
  ascending cell order are as for any tree, so frustum culling draws each unit's own share
  of the same sample, and full detail draws every cell. A tree over cells open in priority
  order is banded whether or not levels were asked for (it costs nothing there); levels
  built under another priority are dropped and built again. Colors are per view, and each
  view's color store is gathered into the order of the geometry it shows; alpha stays in
  source order in its R8 texture, read through the store order's R32UI slot map
  (`u_orderedSlots`: texel `k` is the cell stored in slot `k`, below that count), so a
  visibility change still touches neither stream. Nothing in Cellucid sets a draw order on
  purpose (no "high values on top" or "bring to front" option exists), so no feature
  depends on source order; within a leaf, the later cell is drawn on top.
- **Band by band.** A store holds slots `[0, adopted)` in the tree's order, slots of an
  earlier order below `orderedSlots`, and every later slot in rank order; the GPU reads the
  same mapping from one R32UI slot map (texel `s` is the cell of slot `s`) below
  `u_orderedSlots`, and takes slot `s` itself from there on. Each band is gathered in place
  from its own cells (`GpuPermuter#writeWordsRange`): its slot map texels are written, then
  its positions and the colors of every view presented on the store, copied from the store
  itself on the GPU (a band is at most 2,097,152 cells, one gather pass). That holds too for
  a band still in an earlier tree's order (cells kept by an Apply below a cut, or below the
  cells added): it holds the same ranks in another order, so for the gather its map texels
  first name the slot of that order holding each cell (`earlierSlotsOf`), then the cell; 4
  more bytes a cell instead of the band's positions and colors from page memory. The earlier
  order's slot map already holds that order below its `orderedSlots`, so it is taken over
  when it is large enough (a cut, or a new tree over the same cells) instead of a new map
  being written. Only a band wider than one pass (a store segment forced smaller in tests)
  is read from page memory.
- **Nothing in an adoption waits for the GPU.** The slot map and the gather's source texture
  are allocated with their check read in a later frame when nothing kept by choice could be
  evicted to make room (a refusal could then only be reported); the bands' checks are read
  once, when the store holds the order whole or has given it up (`deferWebGLErrorCheck`'s
  `readAfter`), because in Chromium every `getError` is a round trip that waits for all
  commands queued before it, a whole frame of a large view, even after the fence passed; the
  permute pass reads the generic value of its attribute from the state shadow
  (`gl-state-shadow.js`, read where the context is made and again where it is restored)
  instead of `getVertexAttrib`, which waited the same way.
- **Cells added or removed in place.** `PointStore#resize` copies the one partial segment on
  the GPU and writes only the new slots; where the store holds a tree's order past a cut, the
  band the cut falls in (at most 2,097,152 cells) is rewritten in rank order, so the first
  `n'` slots hold exactly the first `n'` cells. Colors are compared with the published ones
  over the kept cells and sent whole only if they differ (a color scale over every open cell
  can change). The alpha texture, one byte a cell, keeps its storage whenever that holds the
  new count, and only the rows that change are sent: the kept cells' rows whose bytes moved,
  and every row from the one where the old and new counts meet (`_resizeAlphaTexture`,
  `resizedLinearPlan`: the same shape, so the same shader addressing). A progressive
  opening's `capacity` comes from its plan's target; the highlight byte texture and the read
  cells' slot list follow the same rule: storage that holds the new count is written in
  place. The spec compares full detail, culled, a reduced level and culled frames in 2D and
  3D, a raise and a lower across 4,096-point segments.
- **Trees in the worker.** A build whose positions were replaced before it returned is
  dropped. The viewer asks the worker for every large pane's index when cell info is on and
  when a selection tool is armed. A build that runs on the page's thread is admitted with its
  whole peak, since a heap that cannot hold it collects forever instead of failing; a build
  in the worker is admitted with what the page receives of it, the decoded tree's nodes (2 B
  a point, measured 1.1–1.5); the draw order and leaf ids are typed arrays, outside the heap
  the limit bounds. In 3D the tree keeps no leaf ids (`dropDerivableLeafIds`, 4 bytes a cell
  that neither cross from the worker nor stay on the page; the draw order's runs are read
  with `forEachLeafRun`: one run of the order a band); a 1D or 2D unit is a coarse node
  whose cells keep source order, so its leaves keep their ids. Positions of any size go to
  the worker in transferred slices of 64 MiB, one task each, so no main-thread task copies a
  whole geometry and the page never holds a second copy; a geometry's positions are sent only
  once the worker has answered the build before (`_startNextSpatialTreeJob`). A build of the
  tree alone still waiting its turn makes the levels too when they are asked for meanwhile,
  instead of a second build of the same positions.
- **Unplaced cells.** A tree's draw order ends each band with its unplaced cells, after the
  last leaf's range (`buildDrawOrder`); bounds are read by `readPositionBounds`. Every upload
  of positions takes the stand-in `positionsForGpu` (a position store, a kept view's geometry
  store, the edge position texture, the velocity textures and the figure-export upload) and
  only for a geometry that has unplaced cells; the vertex shaders that reject it are points,
  highlight rings and edges. The smoke volume leaves them out; the edge slider's first K
  counts only edges whose both ends are placed in the pane (`connectivityVisibleCells`), so
  the count and the lines drawn agree; the velocity spawn table skips them.
- **The adoption in the first frame.** A cell lookup's tree (`prepareCellPicking`) counts. The
  source-order words go into a transient R32UI texture, one texel a word, and one
  transform-feedback pass writes each stored slot's words, fetched through the order's R32UI
  permutation texture. A store still in source order is its own source: its words are drawn
  into the texture on the GPU (one point a word, at its texel, into the texture as a color
  attachment), so nothing is uploaded from the page; only a store in another order, or whose
  first upload is still pending, uploads its words. Above `GATHER_CHUNK_POINTS` the gather
  uploads the source one chunk of points at a time from page memory, and every target range of
  at most a chunk of slots takes, in each pass, the slots whose source point is in that chunk,
  keeping the rest (read back as a vertex attribute, written to a chunk-sized scratch store and
  copied over the range on the GPU). Its transient memory is one chunk-sized source texture and
  one chunk-sized scratch store at any size; the scratch store is made once a write, not once a
  chunk. The vertex attribute that read a range lets go of it after the passes onto the
  gather's own anchor buffer, which also frees a dormant kept view's color store when it is
  released for memory; `tests/vertex-attribute-release-contract.test.mjs` rejects a pointer
  made with no buffer bound. The error queue is read once a chunk, after its passes and before
  the next chunk allocates; a write whose passes raised an error is run again with every call
  checked, so the error names the call, the pass and the segment. Writes made often (a view's
  highlight words, a draw order adopted a band a frame) are checked later
  (`deferWebGLErrorCheck`): a fence goes after the write, each frame asks whether the GPU
  passed it (which never waits), and the queue is read once it has. A failure then goes to the
  earliest outstanding write whose run, call by call, names it, and says which writes the read
  covered; a write whose source is its own store is never run again. A synchronous check that
  reads such an error first names the writes still outstanding, and the GPU gather hands one
  its own calls did not raise to them. The tree crosses from the worker as typed arrays
  (`encodeSpatialTree`). A recolor is one upload of 4 bytes a point and one GPU pass, whether
  or not the store is in a draw order; kept color stores count against the retention budget
  with retained geometry and are evicted first.
- **Automatic level of detail and culling.** The frames of every cell are measured three at a
  time, which the viewer asks for as soon as the measurement is armed, once the publication it
  measures has committed. A fence polled in tasks says the GPU finished after the last poll
  that found it unsignaled and by the poll that found it signaled, so its share is at most the
  time to that poll and at least the time to the poll before, less the page's share. The GPU
  timer is a `TIME_ELAPSED` query around the frame's commands, exact whenever it is read. A
  measurement's answer bounds its cheapest frame (the least lower bound and the least upper
  bound of its samples, so a first frame or a stall does not count), and its spread is how far
  the next-cheapest frame is proven to lie above it. Below the ladder's first count, 1,000
  cells, every level is every cell, so level of detail is not turned on there. `faster` is
  when the chosen frames' upper bound lies below the frames of every cell's lower bound by more
  than either spread; `not-faster` when the frames of every cell's upper bound is not above the
  chosen frames' lower bound by more than it; `unresolved` when the bounds overlap, which no
  exact instrument reaches. Not proven faster, they go off again for this framing and the tip
  says which ("no faster here", or "closer than this browser can time", with the bounds
  measured as ranges). Frames are measured with culling possibly on (a frame of the whole view
  still draws every cell). An opening that starts with level of detail on, or where it is
  turned on, is not measured even if it is turned off again, and an answer that lands after it
  was turned on is dropped: the choice is made. On a restore the measurement starts once the
  viewer accepts changes again. A query that finds no tree (a click-time pick on a pane whose
  index has not arrived) scans, as above, and projectiles wait; the viewer asks the worker for
  every large pane's index when cell info is on and when a selection tool is armed, so this is
  rare.

### The frame, shaders and sprites

The rules are in {doc}`07_rendering_pipeline_webgl_and_performance_notes`, "What a frame
actually costs". What the code does beyond the rule:

- **What keeps frames going.** Any mutating viewer or renderer call except a republication of
  what is already drawn (an unchanged view label, or a highlight with the same cells, which is
  what a selection preview sends on most pointer moves), any viewer input that changes
  something (a drag, free-look, a selection gesture, `Alt`, navigation keys outside text
  fields, releasing a held key), and any canvas resize marks the next frame for drawing.
  Hovering or typing into a sidebar field draws nothing. The causes that keep frames drawing:
  held navigation keys, inertia, smooth zoom, a grid fade, smoke mode while it animates or is
  still waiting for its noise volumes, an overlay while it draws something that moves
  (velocity particles flowing in a presented view, or a spawn table being built for one; with
  no field for the layout, nothing visible, particles off at this level of detail, or in smoke
  mode it keeps no frames going), visible highlights in a quality whose highlight style
  animates (only `full`: its ring pulses and shimmers; every other quality draws a still ring,
  and a selection there costs nothing while the view is idle), live projectiles, a camera
  motion still waiting for its `camera-settled` event and an adaptive LOD still stepping toward
  the level its camera distance asks for (hysteresis moves one level per drawn frame).
- **Ambient pacing.** A frame that only an ambient animation asks for (a backdrop's ambient
  motion, smoke that animates, flowing velocity particles, the highlight pulse) waits until the
  main thread has been free as long as the last frames held it, and, once the page's frames
  have run a display frame late or more while the GPU was still running the viewer's last frame
  (a fence read when a frame is due, then every half display interval until the GPU finishes,
  so a page that stopped meanwhile, a hidden tab or a long task, is not counted as GPU time),
  waits as long again, the viewer timing the frames in between without drawing; the display
  interval is measured on an idle GPU before the first ambient frame; an overlay's frames that
  converge (a spawn table being built or not yet taken up) are drawn at once. While it pulses,
  each pane's base layer, everything under the highlights, is kept in a framebuffer
  (`rendering/pulse-layer-cache.js`) and a pulse-only frame draws only the layer over it, pixel
  for pixel the frame drawn whole; a multisampled frame keeps the scene target's own samples
  instead, and a pulse-only frame resets and draws again only the samples the rings cover,
  which the ring pass of the frame before marked in the stencil.
- **Settling.** After the last cause the viewer draws two more frames, so the floor only covers
  the frame that shows a change; there is no wall-clock window on top.
- **Error reads in the frame.** After an allocation `_allocateGpu` retries only
  `OUT_OF_MEMORY`, once, after releasing memory kept for speed. Highlight buffers check once,
  after the store is sized, on the first allocation and when it grows past the size already
  verified. A highlight of many cells (when its color words, 4 bytes a slot of the view's point
  store, are fewer bytes than its own records, 16 bytes a cell) is drawn from the point store
  itself: positions from the store, the words written on the GPU from the highlight and
  visibility bytes, the slots not admitted rejected by the vertex shader; this form keeps its
  word store when it fits. The live and kept views' centroid buffers and the velocity overlay's
  per-view tables (visibility and spawn table, packed textures) are rewritten in place when the
  new data fits the storage they hold, with no read; only new storage reads the error queue,
  once (the velocity overlay reads it for both tables a slice after allocating them, before
  publishing). A filter change, a category toggle or a kept view's recolor with the velocity
  overlay on reads nothing; restoring a session of seven views reads once per new table, not
  four times per view. An error nobody reads where it happens is read by the next owner that
  does, so no owner may raise one: WebGL keeps a deleted program current
  (`deleteProgramReleasingUse` in `gl-utils.js` unbinds it first), and Chrome keeps reporting a
  texture deleted while bound to `TEXTURE_3D` or `TEXTURE_2D_ARRAY` on some units as still bound
  there; every texture is deleted through `deleteTextureTracked`, and every owner that saves and
  restores bindings (smoke density and noise builds, the draw-order permute, the edge gather,
  packed textures) restores through `restorableTexture`, which drops a deleted one.
- **The no-depth point shaders.** They write a zero-alpha fragment where the depth variant
  discards, which leaves the pixels byte-identical; the disc's coverage of the pixel is a
  one-pixel ramp across the rim, narrowing below a radius of four pixels so the smallest
  points stay crisp, computed from the sprite radius the vertex shader passes (`v_radiusPx`;
  screen-space derivatives are undefined along a sprite's edge). Depth passes keep the hard
  edge: a partially covered rim fragment would still write depth, and in a dense cloud every
  disc grew a ring of background color where later points behind it were rejected.
- **Per-draw frame values.** The model-view product, the point-size projection factor
  (`u_pointScale`) and the Blinn half vector (`u_halfDir`) are formed once per draw by
  `writePointShaderFrame`; the specular `pow` is five squarings; the point fragment shaders are
  `mediump`. Depth is formed by the same two products in every quality, so the same point wins
  a near-tie whichever quality is selected.
- **Fog.** `u_fogScale` is `log2(e) · density² · 0.6 / span`, derived per pane by
  `writeFogUniforms`, and zero when fog is off. The transmittance is an `exp2` evaluated once
  per point or edge endpoint and passed on as a varying. Over a backdrop of the fog color the
  thinning is exactly the picture the color-mixing law drew, and over anything else (a grid, a
  vignette, the points further back) a fogged point fades into what is really behind it; the
  figure export's vector strategies apply the same opacity, so a circle fades into the paper.
  With fog off a point is drawn at full opacity. The fog parameters are two `uniform1f` calls
  per program per pane.
- **The sprite cap.** `POINT_SPRITE_MAX_PX` is 160 px; a figure export sets `u_pointSizeMax`
  times its raster scale. The point-size slider is anchored at 0.25 at position 0 and the
  largest sprite drawn, `POINT_SPRITE_MAX_PX` (160), at position 100, so no slider travel lies
  past what can be drawn, with the domain extended down to position -24 below the anchor.
  `tests/browser/point-near-camera-solid.spec.mjs` checks, on the GPU, that a cell at 0.13 and
  at 0.52 of the orbit radius is drawn whole at exactly the cap, picked, and in the figure.
- **`sculpted`.** It lights each sprite as a sphere in linear light (the color squared, encoded
  back with a square root, so the round trip is exact) and multiplies by `u_invCenterShade`,
  the reciprocal of the lighting at the sphere's center for the scene light, so the center
  shows exactly the legend color (`sculptedCenterShade` mirrors the shader's constants). The
  gloss fades in between sprite radii of two and six pixels from the vertex shader's
  `v_radiusPx`. In a 3D pass `enterOpaquePointBlending` replaces the color
  (`blendFuncSeparate(ONE, ZERO, ONE, ONE_MINUS_SRC_ALPHA)`) and, when the pass draws into the
  multisampled scene target (`renderParams.multisampled`, from
  `sceneMsaaTarget.isMultisampledThisFrame()`), enables `SAMPLE_ALPHA_TO_COVERAGE`, so the
  antialiased rim becomes a per-sample mask that keeps its own depth: the only way to soften an
  edge with depth writes on without the background showing through a nearer point.
  Single-sampled, the coverage floor of 0.5 gives the same hard edge as every other quality.
  The points are opaque, so fog cannot thin them; they alone still mix toward `u_fogColor`,
  which the figure export sets to its paper. In 1D and 2D views they blend like the other
  qualities. `sculpted` draws `full`'s torus ring with `animated: false`.

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-21 | GPU-heavy work adapts to the measured GPU speed; macOS Firefox CI GL is about 2,500× slower than Metal | No browser preference changes it; giving a spec more time hides the cause | Fixed timeouts; one GPU budget for every machine |
| 2026-10-21 | Pixel-parity specs copy cameras only after the scene settles (`getFrameScheduleState().settled`), assert camera equality, and never use a pixel tolerance | A camera copied mid-ease differs by a few edge pixels on slow GL | Tolerances that hide which input differed |
| 2026-11-12 | A leap of the smoke ray march over an empty block lands on a sample the march without leaps takes, computed the same way (`rendering/shaders/smoke-shaders.js`) | The leap restarted the march at a new sample position, so the march with leaps and the march without sampled different places: 9,953 pixels differed, by at most 8, in Chromium (`smoke-volume-reconstruction`, once its hosted-runner bound was removed, record 38). Now the two give the same pixels, the skipped samples all reading zero | A leap that restarts the sample lattice; a bound on the difference |
| 2026-10-28 | Non-expressing cells are drawn beneath expressing ones, in a neutral color | Sparse expression must not be covered by the gray around it | One pass with arbitrary overlap |
| 2026-10-28 | A color's alpha byte is its layer; filters write only the transparency | A filter never copies or re-uploads colors a kept view shares | Alpha as visibility |
| 2026-10-28 | Levels of detail are fixed counts in priority order; ranks by count, with no admission table | A level is a uniform sample, and growth is in place (per band 10–20 ms and 48 MB → 0 at 10M; 66–79 ms and 229 MB → 0 at 50M) | Spatially chosen levels; per-dataset ladders |
| 2026-10-28 | Levels of detail and culling turn on automatically only when measured too slow, and are kept only when measured faster | Small data never needs them; a 2D full framing got slower with culling | Turning them on by cell count or by dimension |
| 2026-10-29 | Gene on a sample: unread cells fade almost to the background, distinct from the gray for 0; read cells drawn on top at a size set by the sample's on-screen density; no extra downloads | The expression pattern must read at once on a map of all cells | Drawing unread cells as no-value gray; reading more to fill the map |
| 2026-10-30 | The gray for 0 and missing values comes from the theme (fixed L\* contrast) | On a dark background gray 211 outshone the darker colormap cells (seen on a real atlas) | One fixed gray on every theme |
| 2026-11-07 | Ambient animations (smoke, flowing velocity particles, backdrop motion, the highlight pulse) are paced by how late the page's own frames run behind the GPU, confirmed by a fence read only when a frame is due, with the display interval measured | On a slow GPU (CI software GL; hosted macOS Firefox) an endless animation saturated the GPU and every click waited behind it (overlay on: 17.7 page fps under SwiftShader, clicks of 10–15 s on hosted macOS Firefox); a GPU that keeps up is never held back | Pacing every frame by fence latency, which throttled fast GPUs |
| 2026-10-31 | Gray for zero is a per-field switch, off by default; off, a gene's 0 takes the scale's lowest color (as before 2026-10-28); on, the theme's gray beneath. The verbose legend note under the color bar is removed; a gene sample keeps one line naming its states | Owner, on HLCA and Suo: a 0 is a value, and the gray should be a choice, explained by an (i), not a fixed rule and a paragraph | Gray for 0 always on; a note under every gene's color bar |
| 2026-10-30 | A level ranks a cell by its identity under the share's seed, whatever order it is open in | The same cells come first in every order, step, view, revisit and session | Ranking by open position |
| 2026-10-30 | Rendering and app timings only on a real GPU (Chromium `--use-angle=metal` on macOS); software renderers are refused unless asked for (`--software-gl`, never a timing) | SwiftShader draws on the CPU and holds the main thread on synchronous GL calls | Any timing from SwiftShader or llvmpipe |
| 2026-11-01 | A store's draw order is an input of the picture. A view is not settled (`isPresentedGeometrySettled`) while a tree for its positions is being built or is published but not adopted, a cell lookup's tree included; a pixel-parity spec draws the frames it compares in one order (or chooses level of detail and culling by hand) and records it | `unplaced-cells.spec.mjs` drew its reference in a lookup tree's order (the cell info card's idle request ran after the spec changed the positions) and the unplaced frame in source order: 82,845 pixels apart on Linux Chromium and 83,440 on Linux Firefox, reproduced exactly on SwiftShader by building the tree for the reference. `graphics-context-restore.spec.mjs` captured before the automatic choice turned culling on and after it (3 pixels, macOS Firefox) | Comparing frames drawn in different orders; treating a lookup tree as changing nothing on screen |
| 2026-11-01 | The page's own automatic level of detail and culling are verified again whenever their measurement basis changes (a restored context, cells grown or cut, another representation, a new framing), exactly as at an opening: off, frames of every cell measured, on, verified, kept or off. One path: the viewer says `measuring` before each measurement and the render controls turn their choice off then (`applyDatasetPerformance`, which did it for a new dataset only, is removed). A person's or a caller's choice is left alone | Kept across a re-measure, the choice was never checked on the new basis: after a restore its measurement measured frames with culling on as frames of every cell and turned nothing new on, so nothing verified it (found in `graphics-context-restore.spec.mjs` on macOS Firefox); a re-measure with level of detail on was not measured at all | Keeping an unverified choice; a restore-only special case |
| 2026-11-01 | A frame-cost sample is bounds, not a number: the page's share (exact, from the clock) and the GPU's (exact from the GPU timer where the browser offers one, else between two fence polls). Too slow is the GPU's share proven over the budget; the verdict on the automatic choice is `faster`, `not-faster` or `unresolved` from the bounds and the proven spread, and the tip says which | The rule that added the larger of the spread and the poll interval to the cheapest sample counted the resolution twice and called a frame measured 150 ms faster "no faster" (Linux Chromium); a frame slow on the page's thread is not made faster by drawing fewer cells | A noise constant; adding an instrument's resolution to a bound that already holds it; a tip that says "no faster" where the numbers show faster |
| 2026-11-01 | A measurement is armed only once the state it measures is committed (a publication of cells or positions, a restored context), asks for its frame at once, and keeps the scene unsettled until it has answered (`frameCostOwesFrames`) | A publication that failed had armed it first and its `measuring` report turned the page's own culling off (`camera-defaults.spec.mjs`); a measurement armed for the next framing of a still scene waited for whatever drew next and could turn level of detail on under it | Arming before a commit; a measurement that answers at a time something else decides |
| 2026-11-01 | Every browser spec runs on a fixed device (`tests/browser/helpers/device.mjs`, the `device` option of the `test` fixture), so the page's own measurements of its device (level of detail, culling, a progressive opening) never change under a spec; what the device holds and the CI rounds that showed it are in {doc}`38_decisions_testing_and_hosted_runners` | Five successive CI failures were each patched in one spec; the cause was one | Per-spec holds of the automatic choices; specs that pass or fail by the runner's GPU |
| 2026-11-11 | One sprite rule for every cell (`rendering/point-sprite.js`): the cap is 160 px within the GPU's `ALIASED_POINT_SIZE_RANGE`, for every point shader, the highlight ring, the centroids, picking and the figure export. Nothing else depends on the camera's distance: a cell near the camera is drawn whole, opaque and picked. Owner decisions | Close-ups (camera paths, `camera.framing`) drew the cells between the camera and the subject at the 128 px clamp, and four shader copies of the clamp could drift apart; the owner chose the larger cap and rejected any fade of near cells | A near fade (in 3-D orbit panes opacity falling by a smoothstep from 0.5 to 0.2 of the orbit radius, as an ordered dither; near cells under half seen not picked): shipped and removed the same day, the owner: near cells must not fade out as the camera approaches; shrinking near cells (a shrinking ball reads as a far one); a cap of the screen's share (flattens depth among near cells); a fixed limit above the GPU's |
| 2026-11-11 | A highlighted cell is a dot of `HIGHLIGHT_COLOR` (cyan) inside a ring of `highlightRingColorFor` it: black or white, whichever has the higher WCAG contrast with the dot (black around the cyan). One rule in `rendering/shaders/high-perf-shaders.js` for every shader quality (torus and flat rings), 2-D and 3-D, every view, the figure's raster (PNG, Hybrid) and its vector points (`highlightRingRadii`, the shader's ring geometry). Coordinator decision, reversible | Garcia colored by `LVL1` with its 10,993 germ cells highlighted looked exactly as before: dot and ring were both cyan (the flat ring 0.75 of it) and the germ cells sky blue (#4CAFEE), 0.12 OKLab from the dot and 0.07 from the ring (two colors under 0.08 read as one); the black ring is 0.73 from it. The ring cannot know the cells' colors (a highlight is drawn over every coloring), so the pair must read against all of them: any color is at least 0.42 OKLab from the cyan or the black, five times the distance two colors need to be told apart (`tests/highlight-ring-contrast.test.mjs`) | Dimming the other cells while a highlight shows (every highlight page would gray the map for every coloring; the figure keeps it as its opt-in Emphasize selection); a ring in a contrasting color per cell (the highlight pass does not carry the cells' colors); a dot in the page's color (the page color is the category's: the same failure) |
| 2026-11-01 | Every fence poll carries `SYNC_FLUSH_COMMANDS_BIT` (`fenceHasPassed`, `rendering/gl-utils.js`) | Firefox warns at the hundredth poll of an unsignalled fence made without it, whatever `flush` followed the fence; a GPU that takes seconds a frame reaches that count (`renderer-smoke.spec.mjs` velocity slider, macOS Firefox); reproduced on an Apple M1 Pro with slow fragment work, where the flagged poll stays silent | `getSyncParameter(SYNC_STATUS)` or an unflagged `clientWaitSync` as a poll |

### Stores, orders, trees and the frame

A date written "by 2026-11-12" is the day the entry was first written into this record,
where the day of the decision itself was not recorded.

| When | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | Positions are never padded to XYZ: a stream holds the embedding's own components | A 2D embedding costs 8 bytes a point instead of 12 on both the CPU and the GPU and is normalized in place instead of copied | A padded copy of every 2D and 1D stream |
| by 2026-11-12 | A color store is RGBA8, 4 bytes a point, in the live and the kept views | A kept view's store was once packed RGB, 3 bytes a point, but Metal needs a 4-byte vertex stride, so ANGLE kept a converted 4-byte copy of that store and converted it again after every write | Packed RGB for kept views |
| by 2026-11-12 | Position and color streams are not merged into one interleaved record | A recolor would have to repack and re-upload the positions beside every color, four times the bytes it needs | An interleaved vertex record |
| by 2026-11-12 | Per-point textures are linear textures (`TEXTURE_2D_ARRAY`) | A single 2D texture capped them at `MAX_TEXTURE_SIZE²` items: 67,108,864 on a GPU that reports 8,192 (Firefox on macOS), so a 3D position gather (three words a point) failed there above 22.4M points, and the velocity overlay refused datasets above 67.1M cells there (268M where the engine reports 16,384) | One 2D texture for per-cell data |
| by 2026-11-12 | No draw submits more than 2^24 vertices (`drawPointsInSegments`) | Firefox refuses a draw over 30,000,000 vertex ids with `INVALID_OPERATION` and draws nothing; a 100M-point frame costs at most 6 extra calls | One draw a pass |
| by 2026-11-12 | No buffer holds more than one segment of points (16,777,216 slots: 192 MiB of 3D positions, 64 MiB of colors) | A WebGL buffer is capped by bytes (Chromium 2^31 − 2^21, Firefox 2^31 − 1, WebKit 2^32 − 1), so one buffer of 3D positions (12 bytes a point) would stop at 178.8M points on Chromium, and Direct3D 11 caps one resource lower still. The highlight renderer's slot-word buffer (4 bytes a slot) stays one buffer, which binds at 536M points on Chromium | One buffer a store |
| by 2026-11-12 | A per-cell array an engine refuses as an `ArrayBuffer` is made in a WebAssembly memory | Engines cap one `ArrayBuffer` (measured: Chromium 2,145,386,496 bytes, WebKit 2^32, Firefox 2^34), so 3D positions on the page would stop at 178.8M cells on Chromium; a WebAssembly memory holds up to 2^32 bytes on all three engines, and a typed array views it whole. Chromium keeps about 120 memories alive in one page, so only arrays the engine refuses take one. A memory's buffer cannot be transferred, and cloning a typed array copies its whole buffer, so such an array crosses threads in pieces of at most 192 MiB | A page limit at the engine's `ArrayBuffer` cap |
| by 2026-11-12 | A store takes its tree's order band by band; a tree built under a stale priority is dropped when it lands (`builtUnderStalePriority`) | A culling tree requested before a priority-ordered dataset's priority was set landed after it, was taken whole in one frame, and made the next tree's order be taken whole too (a 3.1 s frame at 100M 3D) | Taking an order whole in one frame for cells in rank order |
| by 2026-11-12 | Adoption checks are read late (`deferWebGLErrorCheck`) and the permute pass reads its attribute's generic value from the state shadow | In Chromium every `getError` is a round trip that waits for all commands queued before it; on Firefox and WebKit the read waits for the GPU (10-13 ms after a 10M-point write; Chromium under 1 ms) | A synchronous `getError` after every write |
| by 2026-11-12 | The vertex attribute that read a gather range lets go of it onto the gather's own anchor buffer, never onto no buffer (`releaseVertexAttribute`) | Chromium's command-buffer client never sends a vertex attribute pointer made with no `ARRAY_BUFFER` bound to its GPU process, so the vertex array kept the store's last segment, and the next gather writing that segment with it bound drew `INVALID_OPERATION` ("a transform feedback buffer that would be written to is also bound to a non-transform-feedback target"): a kept view of more than one segment (16.7M cells) could not be recolored. Releasing also frees a dormant kept view's color store when it is released for memory: on Chromium its vertex arrays kept every segment allocated (128 MiB measured for 32M cells) | A pointer made with no buffer bound |
| by 2026-11-12 | The gather takes its source one chunk at a time above 16,777,216 points (`GATHER_CHUNK_POINTS`) | A whole-source texture was 1.2 GB beside the 1.2 GB store at 100M 3D points; Chromium refuses it (`INVALID_OPERATION` on its storage; the 800 MB texture of 100M 2D points is accepted), the refusal threw in every frame, and the view drew all 100M points at every forced level with no message. The transient memory is now 2 × 192 MiB at any size. An adoption the GPU refuses says so once, draws at full detail and tries again when LOD or culling is toggled | A whole-source gather texture; a silent refusal |
| by 2026-11-12 | The gather uploads from page memory with a texture, not from a `PIXEL_UNPACK_BUFFER` | An upload from the store bound as `PIXEL_UNPACK_BUFFER` is the obvious copy, but Firefox on macOS warns on every one | `PIXEL_UNPACK_BUFFER` uploads |
| by 2026-11-12 | The tree crosses from the worker as typed arrays (`encodeSpatialTree`) | A structured clone of a 10M-point octree's 38,568 node objects took 35 ms on the main thread; decoding and installing it takes 10-11 ms | Structured clone of tree nodes |
| by 2026-11-12 | Trees are built in the worker, one at a time, from slices | Building on the main thread took about 2.2 s at 10M points. Levels of detail come up for 100,000,000 points, where the whole-build projection refused them (6.2 GB projected over a 4.1 GiB cap). Four panes of four 3D representations used to post all four position copies into the worker's queue at once; now it holds one | A tree build on the main thread; every position copy queued at once |
| by 2026-11-12 | An Apply of a larger or smaller share writes into storage sized for its target (`capacity`) | A step that allocated or read the error queue waited for every command before it: at 10M cells the read-backs of each step cost 144 ms, with 100 ms more for the new texture | A new allocation, a GPU copy or a `gl.getError` a step |
| by 2026-11-12 | One record of a level (`LodMembership`) with one validator; consumers check the fields they read, never the set of keys | Six hand-made copies that pinned the key set made picking and the KNN edges throw at every rank-order level | Hand-made copies of the record |
| by 2026-11-12 | Input that changes what is drawn is queued with `runBeforeNextFrame` | Animation-frame callbacks run in registration order and the viewer re-registers its own at the start of every frame, so a control that waited for its own animation frame ran after the viewer and showed every slider step a frame late: 50 to 90 ms of latency at 10M cells in four views | An animation frame of a control's own |
| by 2026-11-12 | The frame loop draws two more frames after the last cause, with no wall-clock window | A wall-clock window kept a fast display drawing 60-120 identical frames after every change; a renderer that takes seconds per frame now spends two frames on a change, not a dozen | A time window after a change |
| by 2026-11-12 | `getError` and `isTexture` run only where a failure can appear | `gl.getError()` is a synchronous round trip that waits for every queued command, the previous frame's point draws included: a filter tick at 10M points in four views spent 40 of its 75 ms waiting there. A highlight word store sized anew is cleared by the GPU process before the check can answer (160 MB: 42 ms), so it keeps its store when it fits. Re-binding a deleted texture raised `INVALID_OPERATION` (0x502) that a later dataset switch reported as its own | A read of the error queue at every write |
| by 2026-11-12 | Hidden vertices are pushed outside the clip volume as well as sized to zero, in every vertex shader that draws per-cell geometry, with no exception list | The driver's minimum aliased point size is one pixel (`ALIASED_POINT_SIZE_RANGE` starts at 1.0 on ANGLE Metal, Apple M1 Pro), so a zero-size point is still assembled, rasterized and shaded: a shipped-shader point with alpha 0 lit 1 pixel with multisampling off and 4 with it on, where the same point moved outside the clip volume lit none; at ten million points with everything hidden the GPU frame cost 12 times what it needed and hiding 90% of the dataset rendered no faster than showing all. The ultralight fragment shader has no `discard`, so hidden points were writing depth and occluding visible ones, and a centroid of a category with no visible members occluded real points (1024x768, single-sample, 4,096 visible points with 4,096 hidden in front: without the rejection hidden points took 4,096 pixels from visible ones, with it none, and the framebuffer digest is the same with and without it) | Zero size as the only rejection |
| by 2026-11-12 | One sprite-size function for every shader quality (`pointSpriteSize`) | A higher ceiling in the lighter shaders once made them draw larger sprites than `full`, so the "faster" quality was up to 1.93× slower at close range, and a ring sized by its own copy of the clamp sat inside its dot | A copy of the clamp in a shader |
| 2026-11-11 | The sprite cap is 160 px (128 before) and the slider ends at the cap (before any session was saved with it; a position stored under the old anchor means a size at most 3.4% smaller) | The slider's top was 200: a size above the 160-pixel cap draws the cap for every cell nearer the camera than the perspective scale, so the travel from position 96.5 to 100 drew 160 px for every cell while the readout climbed to 200. The curve was once a `0.25` in the sidebar module and a second `0.25` in `viewer.setPointSize`; lowering the reachable minimum in one place and not the other made the slider throw at its own low end. The reachable minimum, 0.053, is well below the `0.5` clamp floor because sizes under about 0.5 device pixels are distinguishable only where perspective size scaling lifts the rendered size back above the clamp | Two copies of the curve; a slider range past what can be drawn |
| by 2026-11-12 | Fog is one model: one GLSL, two uniforms, thinning instead of color mixing | Before, each layer restated the law inline with its own terms (the anchor by eye depth, the velocity overlay squared, edges with their own alpha floor), so the layers of one scene faded at different rates. Thinning is exactly the picture the color-mixing law drew over a backdrop of the fog color, and fades a fogged point into what is really behind it elsewhere. With fog off a point is drawn at full opacity | A law restated per layer |
| by 2026-11-12 | The 2D no-depth point shaders antialias the disc edge by coverage | MSAA cannot smooth a point sprite's edge: the fragment shader runs once per pixel, and a discard or zero alpha removes every sample of it | A hard disc edge in 1D and 2D |
| by 2026-11-12 | Antialiasing is the app's multisampled target, not a context attribute | It is worth about a fifth of the frame at point size 0.75 and a third at large sizes in 3D, and about twice the GPU time in a dense 2D view where nothing is depth-rejected, the one measured lever; `antialias` is a context-creation attribute and the app has no path to a second context, so the setting is the app's own multisampled renderbuffer (`rendering/scene-msaa-target.js`): live, with three states, applying to the next frame with no reload (the control used to say it applied on the next load), and `auto` follows the cell count and the dimension (`AUTOMATIC_ANTIALIAS_CELL_LIMIT_FLAT`, 200,000 cells for a 1D or 2D view, where the antialiased frame fits a 120 Hz budget on far smaller GPUs). A click is the choice even where the browser keeps no storage (it was rolled back there), and a session restore must not dispatch `change`: that made a session the last writer of a device preference, so antialiasing could never stay off after a reload | A context-creation attribute |
| by 2026-11-12 | `sculpted` is opaque by design | The only way to soften an edge with depth writes on, without the background showing through a nearer point, is alpha to coverage on a multisampled target | A translucent `sculpted` |
| by 2026-11-12 | Smoke volume memory does not depend on the grid | The transient float atlases peak below 50 MiB at 128³ | A grid cap as a memory lever |
| by 2026-11-12 | Fixed-function GL state (`isEnabled`, `getParameter` of `ACTIVE_TEXTURE`, `VIEWPORT`, `SCISSOR_BOX`, `COLOR_WRITEMASK`, the blend functions and equations, `COLOR_CLEAR_VALUE`, the enable capabilities, a generic attribute's value) is read from a shadow of the context (`rendering/gl-state-shadow.js`), never asked of it per frame | In Firefox and WebKit these reads are round trips to the GPU process that first run every queued command: behind one 64 MB texture upload 6-8 ms in Firefox and 16-24 ms in WebKit, 0 ms for a binding or a pixel-store parameter in every engine. Chromium answers `getVertexAttrib(CURRENT_VERTEX_ATTRIB)` the same way. The backdrop look asked three of the waiting kind per pane per frame, the GPU gather one per draw-order change | Reading state from the context in a frame |
| by 2026-11-12 | The shadow is primed again after `WEBGL_lose_context` and a restore (`primeGlState`), and the viewer sets its clear color again in `restoreLostGraphics` | In Firefox 159 the shadow kept answering with the values from before the loss (clear color the background, the context's own transparent black), and the highlight pulse's sample reset drew the former where whole frames cleared to the latter; a restored page cleared every frame to transparent black under a black, white or grid background | Keeping the shadow across a context restore |
| by 2026-11-12 | The viewer draws only after something can have changed: pointer and focus input, cell lookups and viewer calls that change nothing drawn (preview borrowing, listener registration, input preferences) are not mutations | The viewer redrew every animation frame, a full GPU pass sixty times a second for an untouched dataset; each such call invalidated the frame and so kept drawing for a settle window | Drawing every animation frame; a mutation count that includes calls that change nothing |
| by 2026-11-12 | A backdrop fade moves toward its target in either direction | A grid whose target dropped (Grid light 0.85 to Grid dark 0.75) stayed above it, and the scheduler redrew every frame for good | A fade that only rises |
| by 2026-11-12 | Adaptive level of detail names its convergence as a frame cause (`'lod'`), and the floor after a change covers only the frame that shows it | The viewer covered a multi-level jump of level of detail with a blanket floor of twelve frames after every change | A blanket floor of frames after a change |
| by 2026-11-12 | The highlight pulse draws only the layer over the points; the pane's base layer is kept and reused (`rendering/pulse-layer-cache.js`) | Every pulse frame drew the whole pane, every point included | Drawing the whole pane for each pulse frame |
| by 2026-11-12 | A scene is not settled while the opening's frame-cost measurement awaits a sampled frame's fence, or is armed and the loop can draw (`getFrameScheduleState().settled`, `frameCostOwesFrames`); a disposed viewer cancels its deferred WebGL error checks; a context loss inside a frame is told from a defect | On a GPU that takes seconds a frame the measurement's next frame comes seconds after the loop stopped, and `idle-frame-scheduling.spec.mjs` counted those frames as idle draws (hosted macOS Firefox, a nightly run). The polls of the deferred error checks otherwise went on after the last frame's fence (`viewer-terminal-disposal.spec.mjs`, same runner). Rethrowing a loss reached `window.onerror` and reported a handled loss as a crash, which a page-error assertion in the terminal-disposal spec caught | Counting a measurement's frames as idle draws; rethrowing a context loss |
| by 2026-11-12 | The late-frame count of ambient pacing covers the display frames up to the first frame after the GPU's last busy read (one more for rounding) | A page hidden for 10 s right after a slow ambient frame held its ambient animation for 10 s once shown | Counting the whole gap as GPU lateness |
| by 2026-11-12 | An ambient pacer decides a frame from the frame just drawn, not the one before the last | Deciding from the frame before the last drew pulse frames back to back when the next animation frame runs before queued tasks: on Firefox 450 ms frames came in pairs and the page was held 900 ms at a time | Deciding a frame from older frames |
| by 2026-11-12 | Memory, not a count, decides what stays on the GPU for kept views; the layout change to one view presents a frame before views are dormant | There used to be room for four views; on WebKit, 1 run in 7, nothing was dormant, so the refused allocation had nothing to release and was thrown | A fixed number of views on the GPU |
| by 2026-11-12 | The camera's held keys are credited with the share of each frame they were held (`rendering/navigation-keys.js`), not "is it down now" at frame time | Sampling at frame time lost a tap shorter than a frame and credited a held key with the whole frame it was seen in: on a slow frame (1 to 8 s a frame for 562K cells) equal 400 ms holds moved the camera -0.586 against +0.034 rad | Sampling key state once per frame |
| by 2026-11-12 | A view's geometry is its layout, not its dimension: a switch to another layout of the same dimension publishes that layout's positions and drops every per-view cache derived from the previous geometry | Keyed by dimension the switch was a no-op, and a view kept level-of-detail, frustum and spawn state describing the embedding it no longer showed | Caches keyed by dimension |
| by 2026-11-12 | The positions of an array are scanned once for their per-axis bounds (`POSITION_BOUNDS`), and a representation switch back to a layout the viewer has shown costs no pass over the cells | The viewer, the projectile system, the renderer and every snapshot of the same geometry each scanned the same array: three to four full passes on each dataset publication, 100 ms at 10,000,000 2-D cells (machine not stated); a full scan is 30 to 60 ms at 10M 3-D cells (`assertCellPositions`), so a switch back to a shown layout skips it | A scan of the positions by each consumer |
| by 2026-11-12 | One embedding, one set of displayed coordinates: the payload route hands its array over to the dimension manager and forgets it | The adapter handed the manager the very array it kept in its own cache, so the in-place display normalization rewrote the adapter's cache, a reload after an eviction normalized that array again, and the eviction released nothing while the adapter still held it | Two owners of one array |
| by 2026-11-12 | A row is left out when any value of it is NaN or infinite (`data/embedding-payload.js`), for positions and vectors alike, and the notice (`app/unplaced-cells-notice.js`) names both: "NaN or infinite in the file" | A vector field whose row held an infinity said "(NaN in the file)" | A notice that names NaN alone |
| by 2026-11-12 | The centroid pass returns exactly the categories with cells | It had returned one entry per category of the inventory, so a category without cells got a centroid at the origin with `n_points` 0 and every consumer built a point and a label for it: a field with 70,000 categories of which three hold cells drew 69,997 labels stacked at (0, 0, 0) | One centroid per category of the inventory |
| by 2026-11-12 | Adaptive level of detail is bounded by a point budget and a ratio together (2,000,000 drawn points or an 8× reduction, whichever is smaller) | A fixed ratio made a large dataset unusable on Auto: the coarsest step of the ladder is 44× whatever the cell count, so pulling the camera back on 18.1M cells discarded 97.7% of them. An absolute budget alone broke the other end: no level of a 200,000-cell cloud holds two million points, so Auto answered full detail at every distance and the forced-level slider was the only control that did anything | A ratio alone; a budget alone |
| by 2026-11-12 | Every vertex shader pins its attributes, one of them to location 0 | With the locations left to the linker, attribute 0 was whichever input it chose, and a program drew without it whenever no earlier pass had left location 0 enabled (the instanced connectivity edges did once the default background stopped drawing the grid first); Firefox on desktop OpenGL then emulates the attribute for every vertex and warns | Locations chosen by the linker |
| by 2026-11-12 | The field of view applies 45 degrees to the presented area's shorter side, and a pane of a kept view takes its own field of view the same way | It was 45 degrees vertically on every screen, so on a phone held upright a camera framed to the data's height cut the data off at both sides | A vertical field of view fixed on every screen |
| 2026-11-11 | The rows of the grid for a pane count are chosen by the largest shorter pane side; each pane's own field of view follows its own shape | A 390 px phone laid three panes of about 120 px in one row. A kept view cropped the Pancreas 2-D UMAP when two views side by side made portrait panes under the presented area's landscape field of view | Always side by side; one field of view for every pane |

### Velocity overlay, shader quality and the analyzers

| When | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | One tone map for the velocity trail keeps the hue: the luminance is mapped and the color scaled with it. Every pass reading the trail clamps it to a finite ceiling (a half-float sum can overflow), and every blur takes its weights from one normalized Gaussian | A flow blown out at the dials' far ends drew magenta and cyan rims on a light backdrop. Mapped channel by channel, a saturated color's largest channels reached white first and its smallest lagged: a purple flow (the slow end of viridis, green near 0) turned magenta before white and a teal one cyan; its contrast curve and the bloom's screen blend assumed display values on HDR ones. The hand-typed bloom kernel summed to 1.0008, a 0.08% gain every pass, and its zero last tap read two texels for nothing. The trail held no NaN or infinity | Per-channel tone mapping; hand-typed kernels |
| by 2026-11-12 | The newest operation owns the velocity overlay's loading block; of several fields only the one bound to the shown embedding is chosen, none is guessed | A field chosen while the overlay was turning on left it "loading" for good, from every way in (`setting.set`, a session or history restore): the field change took the operation over and returned because a load was under way, and the load, no longer the owner, returned without clearing its loading state; the switch and select grayed and no particle drew | An operation that returns while another owns the block |
| by 2026-11-12 | The velocity spawn table is refilled each respawn epoch from where the last one stopped, so every visible cell becomes a spawn source (owner decision, the 50M to 100M scale work) | A fixed 65,536 sources | A fixed spawn set |
| by 2026-11-12 | A field whose every vector is zero is said, and a field that could not be read says its failure once | A zero field turned on silently and drew nothing, though the read knows its longest vector; a failed read was said twice (the loading notice failed with it and the overlay's caller said it again) | Two voices for one failure |
| by 2026-11-12 | Shader quality is a look choice, not a speed control: no tooltip says it reduces GPU work, and the bottleneck analyzer prints no `Shader overhead:` line | Full, Light and Ultra-light are within noise at every point size, term by term, because about 95% of rasterized fragments are killed by the depth test before the fragment shader runs (depth writes are on and submission order is uncorrelated with depth). A 1D or 2D view draws without depth, so every overlapping fragment is shaded and blended; any frame-time difference quoted for it is scoped to 2-D views and names the machine. The tooltip had said lowering it "reduces GPU work for smoother interaction", and the analyzer printed a delta of full minus ultralight; both were removed. Point size is the real control (per-sprite rasterization) | Presenting quality as a speed lever |
| by 2026-11-12 | The bottleneck analyzer draws through the renderer's public `render()`, and its contract is checked against the renderer's; its level-of-detail sweep orders levels by points drawn | `pointSizeMax` was added to the contract and the analysis hung for the person's Analyze performance. The sweep ordered levels by index (level 0 is the smallest sample), so a frame that did not follow the point count scored about 10,000 and was named "Vertex Processing" on every run, smooth or not | A private draw path for the analyzer; levels ordered by index |

## Measurements that justify it

- **Ambient pacing** (2026-11-07, Apple M1 Pro; page fps is the page's own rAF rate): SwiftShader with the velocity overlay, page 17.3–17.9 → 34–38 fps, overlay 17.9 → 17.5–17.8 fps, clicks 111–148 → 30–121 ms; SwiftShader smoke, page 25 → 41 fps; default headless Chromium overlay, page 17.6 → 45 fps, clicks 650 → 22–94 ms; native Metal (Chromium 60, Firefox 120, WebKit 60 fps) unchanged.
- **Real atlas walk** (2026-10-30, Chromium on Metal, Apple M1 Pro, HLCA
  full, 2,282,447 cells): opening 12.9–17.0 s, a gene 4.6–5.4 s on a sample
  of 1,759–2,273 cells, read-size multiplier 8.6–10.35 at a base point of
  0.32, and a PNG figure with the three layers in 0.77–0.8 s. On the light
  theme the pattern showed at once. On the dark theme gray 211 outshone the
  darker viridis cells, which led to the 2026-10-30 gray decision. The
  rerun had 0 page errors.
- **The zero-gray rule** is checked in `cellucid/tests/read-cell-layer.test.mjs`.
  It gives a fixed contrast on every background, 211 on white and a
  mid-dark gray on dark, and on every theme, backdrop looks included, at
  least half of every colormap stands out more than the gray does.
- **Growth without re-allocation** (2026-10-30, real-GPU page profile at
  10M). Each step of a progressive opening used to rebuild the alpha texture
  and read GPU errors back synchronously, about 250 ms of page time per
  step. Stores and textures are now sized for the growth target and kept.
- **Page color cost on Metal** (10M): a gene shown in about 90 ms before
  sparse forms crossed to the page, 60–78 ms after on every source, with no
  long task. Local 10M: 78–84 → 65–66 ms. See {doc}`27_decisions_memory_model`
  and {doc}`23_decisions_partial_reading` for the reading side.
- **macOS Firefox CI GL** (2026-10-21, GitHub-hosted runners, a 1024²
  fragment pass): Chromium headed on ANGLE Metal 22 ms; WebKit 31–154 ms;
  Chromium headless on SwiftShader about 500 ms; Linux Firefox on llvmpipe
  about 330 ms; Linux Chromium and every Windows engine 930–1,770 ms;
  macOS Firefox (OpenGL on the paravirtual GPU) 55,000–69,000 ms.
- **Specs:** `gene-expression-layering` passes in 2D and 3D at every level
  of detail on Chromium, Firefox and WebKit, with the passes and alpha bytes
  verified by a draw probe. The growth specs (open cells in place,
  progressive subsample, subsample control, unplaced cells) passed 15/15 on
  all three engines on 2026-10-30.
- **Software-GL caveat.** The orbit timings that motivated the automatic
  levels of detail were taken on headless Chromium with software GL: 3D 1M
  433 → 213 ms, 3D 10M 1,200 → 600 ms, 2D 1M 438 → 721 ms. Under the
  2026-10-30 rule they show direction only. Re-measure on a real GPU before
  citing them as absolute numbers.
- **Where the GPU time goes** (one Apple M1 Pro through ANGLE Metal, 1440×1000 at device
  pixel ratio 1, 10,000,000 synthetic points, LOD and frustum culling off, GPU timer queries,
  counterbalanced arms with a byte-identical twin in every set; the absolute milliseconds belong
  to that machine, the ratios are the durable part):

  | Stage | Measured | Verdict |
  | --- | --- | --- |
  | Vertex | 2.62 ms, flat at every point size, 4.2% of the frame | fetch-bound; hoisting the redundant matrix product and the per-vertex `tan()` had no effect |
  | Fragment | 190× the fragments costs 1.82× the time | only about 5% of rasterized fragments survive the depth test to be shaded |
  | Per-sprite rasterization | about 5.7 ns a sprite, about 175M sprites a second | the floor |

  Fill reaches parity with the per-sprite cost at `gl_PointSize` ≈ 15 px: below that the frame
  is primitive-bound, above it fill-bound.
- **Fragment-shader simplification, 3D** (same machine): fog, lighting, the alpha `texelFetch`,
  the round `discard`, the mere presence of `discard`, and the whole `full` → `light` →
  `ultralight` ladder (43 → 9 → 2 statements in the translated Metal) each stayed inside its own
  twin band at every point size, because depth testing rejects hidden fragments before shading.
  In 2D, which draws without depth, every overlapping fog-translucent fragment is shaded and
  blended in order, per sample under MSAA: at 1,000,000 cells with antialiasing the no-depth
  variants took `light` from 30.7 to 17.0 ms and `full` from 33.8 to 30.6 ms.
- **What 0.9.2 moved out of the per-vertex and per-fragment work**: none of it measurable on the
  machine above, and all of it work a smaller GPU pays for. `u_pointScale` replaced
  `u_viewportHeight` and `u_fov`; the model-view product, `u_pointScale` and the Blinn half
  vector `u_halfDir` are formed once per draw by `writePointShaderFrame`; fog is evaluated once
  per point instead of an `exp` per fragment.
- **Fog off** (10,000,000 points, device pixel ratio 2): 2D 99.7 → 92.4 ms, 3D 87.8 → 80.4 ms a
  frame; with fog on the frame is unchanged. A compile-time variant for fog off would remove one
  `exp2` per vertex and was not added; a shared uniform block would save microseconds a frame and
  add a buffer binding to every program, so it was not adopted.
- **Edges as a strip**: every vertex projects both endpoints, so four vertices are a third less
  vertex work than six.
- **The smoke ray march** does less arithmetic per sample (the density-texture coordinate map and
  the volume normalization once per pixel, the three multiple-scattering octaves share one `exp`,
  the blue-noise jitter fetched once): 17.2 against 14.5 ms a frame on the 562,000-cell Suo
  sample, same pixels (median of three interleaved runs, Apple M1 Pro, device pixel ratio 2).
- **A different primitive**: `GL_POINTS` beat instanced quads and both triangle arrangements by
  1.41–2.12× at one sample and by up to 1.90× with MSAA enabled at 4 samples, with the quad arms'
  rasterized and shaded fragment counts matched to five decimal places. Halving the primitive
  count with one triangle per point was the worst arm. Order-independent transparency, a depth
  prepass, storage reordering and the vertex-ALU hoist were all measured; none won.
- **Antialiasing** is worth about a fifth of the frame at point size 0.75 and a third at large
  sizes in 3D, and about twice the GPU time in a dense 2D view, where nothing is depth-rejected:
  orbiting 1,000,000 cells, GPU timer totals 25.5 ms on against 11.6 ms off with `full`, 17.5
  against 9.9 ms with `ultralight` (`scripts/render-scenarios.mjs --antialias on|off`).
- **Multiview** (Apple M1 Pro through ANGLE Metal): every view pays about 5–7 ns for every
  on-screen point sprite, whatever its size; a 10,000,000-point 2D embedding took 30 ms for one
  view and 147 ms for four, and 13.5 ms for the same four views with level of detail.
- **The stored order** (Apple M1 Pro; standalone WebGL2 harness on the real draw order,
  translucent 2-px points, all cases in one run on a loaded machine): 10,000,000 points drew in 18
  ms against 68 ms in source order on Chromium (15 against 57 on WebKit, 14 against 57 on Firefox),
  a reduced level of 8,000,000 points in 15 ms against 14 ms for the compact copy it replaced, and
  a zoomed-in culled frame in 7–8 ms against 21–31 ms for the per-point index buffer it replaced.
  Adopting the order of a 10M-point frame draws it about 3× faster in 3D and in 2D (21 ms against
  67 ms, no vsync, Apple M1 Pro).
- **A ranked order** is two radix passes: 82 ms for 10M cells in Node on an Apple M1 Pro.
- **Adoption at 10M 3D points**: 14–19 ms of CPU, once (it uploaded 120 MB of positions and 40 MB
  of colors: 32–35 ms); turning LOD on used to block the page for 2.3 s. At 100M 3D points before
  the chunked gather every forced level drew 100M points at 383 ms a frame; now the levels draw
  2.3M–51M points.
- **Cells and arrays past one `ArrayBuffer`**, measured in Chromium: 178,956,971 3D cells (2^29 +
  1 floats, refused as an `ArrayBuffer`) generated in the worker, assembled in a memory on the
  page, their bounds (1.0 s) and cell-order digest (3.0 s) read, and their tree built in the
  worker from slices (50 s).
- **Connectivity edges**: a k = 15 graph of 50M cells holds 3.0 GB of forward adjacency. The
  selection gives the records of the old ordered-records layout bit for bit, proved against it on
  synthetic graphs (648 selections, 12.3M records). At 10M cells (70M edges, K = 250,000): no
  hidden cell 0 ms after the prefix (about 0.1 s once, when K is set); a filter leaving half the
  cells 4–11 ms; one leaving a fifth about 0.2 s (the cell walk).
- **The read size**: on a whole-screen map, 2,000 of 2.3M cells draw about ten times larger; 1M
  of 10M are dense already and draw normally.
- **Frame-cost instruments** (2026-11-01): the GPU timer (`EXT_disjoint_timer_query_webgl2`) is
  in Chromium on a hardware GPU; Chromium's software GL, Firefox and WebKit do not have it. A
  fence polled in tasks answers within 6 to 10 ms in WebKit and 4 to 6 ms in Chromium on an idle
  Apple M1 Pro, and within up to 150 ms on a loaded hosted runner. A 2D view at full framing is
  measured no faster with culling (721 against 438 ms a frame at 1M cells on software GL) and
  keeps every cell.
- **Waits and costs that set a path** (Apple M1 Pro through ANGLE Metal unless stated; no date is
  kept for most of them):

  | What | Measured | What it set |
  |---|---|---|
  | The first read of a shadowed GL value during a frame of a large view | Waits for that whole frame: 0.2 s for `colorMask`, first read by the GPU permute's source draw during a 100M-point adoption; Chromium's `getVertexAttrib` was asked once a pass and waited for the frame before it, 2.3 s over the adoption | `primeGlState` where the context is made |
  | The error queue read after an allocation at the start of a 100M 3-D adoption | 0.7 s | `_allocateGpu` checks in a later frame when there is nothing to evict |
  | `getError()` read after in-place writes by an overlay | An overlay rebuilding the tables of seven restored views spent 112 ms there; two error reads in the allocating slice cost up to 90 ms a view on the Suo atlas | `requireLiveContext` for sub-uploads; one read, a slice later, for both textures |
  | A highlight of every cell at 10M points | Gathering millions of positions on the CPU 72 ms, the store cleared by the GPU process before it is written 54 ms more; three `getError()` checks in a publication cost 104 ms of a 199 ms highlight | A highlight draws from the view's point store when that writes no more bytes; one check, after the store is sized |
  | A shader program compiled on first use | Holds that frame for about 100 ms (the pulse's reset program is compiled ahead in `prepare()`); the permute and highlight programs took 90 to 190 ms the first time a machine compiled them (a 10M tree's adoption, or a highlight of millions of cells) | Programs are compiled ahead of first use (`prewarm` in `high-perf/gpu-permute.js`) |
  | The gather's source texture | Metal has no three-channel 32-bit format, so an RGB32UI texture would be converted on the CPU: 50 ms against 13 ms for 10M 3-D points | One-channel R32UI, one texel a word |
  | Tree transfer to the main thread | Tens of thousands of per-leaf transfers made the receiving task cost about 50 ms at 10M points; a structured clone of the node objects costs about a microsecond a node on the receiving thread (35 ms in one task for the 38,568 nodes of a 10M-point octree) | One shared leaf-id buffer; `encodeSpatialTree` |
  | Adopting a draw order band by band | With 8M 3-D cells, 8M in one frame held the main thread 29 ms (Chromium), 54 ms (Firefox) and 104 ms (WebKit), mostly writing the slot map's texels; at 2^21 points a frame the adoption takes a median 12.5 / 3 / 23 ms (the first, which allocates the map and the source, up to 22 / 40 / 65 ms) | `STORE_ORDER_POINTS_PER_FRAME`, 2^21 points a frame |
  | Frustum-culled draws without `WEBGL_multi_draw` | At 10M points, 252 coarsened ranges drew a zoomed view in 7 ms on all three engines, against 21 to 31 ms for the per-point index path they replaced | `MAX_CULLED_DRAW_RANGES` (256) |
  | Separate position and color streams | With 5M points on Chromium/Metal the split full-detail draw was within noise of the interleaved one (positions and colors used to share one 16-byte record, so every color change repacked and re-uploaded the positions) | Two streams |
  | The draw order's units in 1D and 2D | In 3D the leaf order makes a 10M-point full-detail frame about 4 times cheaper than source order. With no depth, every point blended and multisampling, leaf order is the slowest order: 10M 2-D points drew in 300 to 500 ms against 145 to 195 ms in source order (Chromium, benchmark harness); in the app a 10M-point 2-D full-detail frame drew in 21 ms with 64 coarse units against 67 ms in source order (no vsync, GPU-bound) | 64 coarse units for 1D and 2D (`DRAW_UNIT_DEPTH_BY_DIMENSION`) |
  | Firefox's limit on one draw | It refuses any draw over 30,000,000 vertex ids (`webgl.max-vert-ids-per-draw`) with `INVALID_OPERATION`: 30,000,000 is accepted and 30,000,001 is not (probed on Chromium, Firefox and WebKit); Chromium and WebKit drew 200M points in one call (the 50M to 100M scale probe) | The segment limit of `drawPointsInSegments` |
  | Ambient animation on the main thread | Under a software renderer or a headless browser that reads every WebGL frame back, a marker run of 400 genes (32 ms of statistics) took over 40 s, because its worker's replies only ran between frames | Ambient pacing |
  | Pacing under SwiftShader | A pipeline hides the cost of the first two or three frames of a burst: the lateness of three frames came two frames after them (why the loop times the frames it runs without drawing). Pacing every frame by a fence took smoke on Metal from 60 to 45 frames a second. The velocity overlay of `data-mapping-multi.h5ad` on 2026-11-07: page 34 to 38 frames a second instead of 17, a change shown in 60 ms instead of up to 112 | Pacing on the page's own lateness |
  | The display interval | 2026-11-07: two animation-frame periods measured 14.5 ms in WebKit (60 Hz) and 7.5 ms in Firefox (120 Hz); WebKit's timestamps are the callbacks' own times, up to a third of an interval apart from the display's | The interval is the mean of eight periods (`INTERVAL_PERIODS`); only the first and last timestamp count in it, so the two ends' third each over eight periods leaves it within a twelfth of an interval, far inside the half a late frame is told by (`LATE_AT_INTERVALS`, 1.5) |
  | The frame-cost measurement on a GPU that takes seconds a frame | The frames it asked for were most of an opening's time (macOS Firefox on hosted runners) | It stops where no presented pane would draw fewer cells (`nothing-fewer`) |
  | The smoke noise volumes | 256³ cost 146 MiB of GPU memory as a pair (18 MiB at 128³); GPU generation (`startCloudNoiseGenerationGPU`) is about 100 to 200 times faster than CPU generation (machine not stated). A generation turn that took one phase each made a generation some twenty-five animation frames long, well over a minute for a 32-cubed volume where a frame takes seconds (a software renderer) | `MAX_SMOKE_NOISE_RESOLUTION` 128; phases batched by `batchBudgetMs` |
  | A discard in the no-depth fragment shader | Two separate discards cost 12% of a dense antialiased 2-D frame (1M points, 4x multisampling): 35.6 ms against 31.4 ms for the same pixels | One discard in `HP_FS_FULL_NO_DEPTH` |
  | A point-size ceiling in the lighter shaders | A ceiling of 192 measured 1.23 times the GPU frame time at point size 40 and 1.90 times at 60 | One sprite-size function for every quality |
  | Edge draw inputs | The same vertex arithmetic on texture-fetched positions rasterized some 3-D connectivity lines differently from attribute-fetched ones on Metal (measured 2026-10-23) | The edge draw keeps attribute inputs, and the GPU gather copies bits into instance records |
  | The backdrop figure lattice (`backdrop-looks.spec.mjs`) | The share of a figure PNG's lattice darker than 192, 0.16 in WebKit and 0.12 in Chromium, and the PNG that held no point (measured in {doc}`36_decisions_figure_export`, Measurements) | The spec measures the lattice, not the count of colors |
  | The no-depth point shaders against a discarding twin | They were compared within one rounding step per channel, after hosted macOS Firefox drew one pixel of Light differently once it gained a rim; since 2026-11-12 they are compared exactly (record 38) | The tolerance of `point-shader-no-depth-parity.spec.mjs` |
  | Antialiasing on a flat view | Orbiting a 1,000,000-cell 2-D view, the GPU timer total went from 11.6 to 25.5 ms a frame with `full` (2.2 times) and from 9.9 to 17.5 ms with `ultralight` (`scripts/render-scenarios.mjs --antialias`), against 20 to 35% in 3-D, where depth testing rejects most fragments first | `AUTOMATIC_ANTIALIAS_CELL_LIMIT_FLAT` |

## Tried and reverted

- **Pacing every ambient frame by a GPU fence's latency** (2026-11-07). A fence measures latency, not occupancy, so it held back GPUs that keep up: smoke on Metal went from 60 to 45 fps. Replaced by pacing on the page's own lateness, confirmed by the fence only when a frame is due.
- **Unread cells as the no-value gray** (design of 2026-10-28). Cells not
  read were NaN and drawn in the gray beneath, with the difference stated
  only in text. On a large map a 0 and "not read" looked the same. It was
  replaced on 2026-10-29 by the fade and the larger read cells.
- **One gray (211) on every theme.** It was replaced on 2026-10-30 by
  `zeroCellColor`. On light themes the gray is unchanged.
- **`hasBeneathCells`** was replaced by `colorLayersOf` when the third layer
  arrived.
- **Ranking a share's cells by open position.** A share not open in its priority
  order, chunk steps included, was ranked by open index under the default seed, so the
  cells drawn first changed with the order the cells were open in. Since 2026-10-30
  ranks follow identity.
- **An earlier fixed ladder of 18 reductions** was replaced by fixed counts
  in priority order. The coarsest reduction, 44, survives as the farthest
  camera's size cap.
- **Automatic culling on a 2D full framing** made frames slower (above). It
  led to the rule that an automatic choice is verified and turned off again
  when it is no faster.
- **A growth fix counted capacity as held bytes.** `PointStore#byteLength`
  briefly reported capacity, so recoloring was refused during growth. It
  now reports the cells held, and `allocatedByteLength` reports the storage.
- **A layering spec set up wrongly.** Its "2D" view still had depth and its
  3D sheet faced away at an angle. The spec now uses `setViewEmbedding` and a
  sheet facing the camera. The lesson: prove the inputs of a pixel
  comparison before trusting it.
- **Fragment-shader simplification in 3D views, per-vertex arithmetic hoists**: measured no
  effect ("Measurements that justify it"); not to be re-opened without a new measurement.
- **Instanced quads, two triangles or one triangle a point** instead of `GL_POINTS`:
  measured slower (above).
- **Order-independent transparency, a depth prepass, storage reordering for fill, the
  vertex-ALU hoist**: measured, none won.
- **A compile-time variant for fog off, a shared fog uniform block**: not added (above).
- **A packed RGB kept color store** (3 bytes a point): replaced by RGBA8 because Metal needs a
  4-byte vertex stride and ANGLE kept a converted copy that it converted again after every write.
- **A whole-source gather texture**: Chromium refuses it at 100M 3D points (see "Stores, orders,
  trees and the frame").
- **A near fade or shrink of cells by distance, a higher point-size ceiling in the lighter
  shaders, a slider range past the cap**: see the 2026-11-11 rows and "Stores, orders, trees and
  the frame". The near fade (an ordered-dither opacity from 0.5 to 0.2 of the focus
  distance in 3-D orbit panes; cells nearer than 1.5 of a camera 3 from the target were left
  unpicked while it lasted) shipped and was removed the same day, 2026-11-11, on the owner's
  decision that near cells must not fade; the sprite cap went from 128 to 160 px the same day at
  the owner's request. `tests/browser/point-near-camera-solid.spec.mjs` fails if any
  distance-based fade or shrink returns (it also caught the orbit anchor's compass, whose red
  mark showed as one red pixel beside the near cell once the target came close to the eye:
  "310 px against the cap 160" on WebKit). It would need the owner's word.
- **A frame-time verdict on shader quality from the wall-clock interval between
  `requestAnimationFrame` callbacks.** That interval is quantised to the display refresh and
  mixes CPU, GPU, compositor and scheduling; an implementation that measured one quality twice
  and used the twin difference as a noise floor was observed live to pass a 22 ms "separation"
  over a 7.5 ms floor on a machine whose frames were 100% janky and whose true effect is zero.
  The instrument now emits no verdict (`tests/shader-quality-look-not-speed-contract.test.mjs`),
  and the feature-toggle analysis restores no shader quality and has no shader-quality sweep:
  the three qualities measure the same on instrumentation that can resolve a difference.
- **A verdict that overrode a caller's own setting**: a verdict landing after a caller's
  `setAdaptiveLOD` or `setFrustumCulling` turned the caller's culling off or level of detail on
  again, mid-benchmark and mid-test (a hosted runner's software GL measures most openings too
  slow). Every report and every such call now names the opening it belongs to.
- **A rule that added the spread and the poll interval** to the cheapest sample: see the
  2026-11-01 rows.

## Known limits and open items

- **Highlight rings** are drawn at the level's point size, not the read
  size, so a ring around an enlarged read cell is too small
  (`rendering/highlight-renderer.js` sizes rings by
  `getCurrentLODSizeMultiplier` only).
- **Smoke** substitutes layer colors, but a volume has no layers or sizes:
  read cells are not enlarged there.
- **The zero-gray walk on a real atlas** (light, dark, zoom, coarse level,
  kept view, figure) has not been done. The first attempt stalled while the
  atlas was still opening, before the rendering code ran.
- **Planned renderer work:**
  - draw PNG and Hybrid figures from the GPU stores in the page's context,
    so they get layers, LOD and read cells for free;
  - let kept views share the live color store and R8 alpha (copy-on-write);
  - give edges and velocity one shared position texture;
  - reduce typed-array allocation sites;
  - split `rendering/viewer.js` and `rendering/high-perf-renderer.js`.
- **`full-detail-frame-cost` on Firefox** read a 65 ms frame once while two
  heavy jobs ran. Measure on a quiet machine.
- **Automatic level of detail and culling: a verdict is what the bounds
  prove** (2026-11-01). A fence polled in tasks says only that the GPU
  finished between two polls: 6 to 10 ms apart in WebKit and 4 to 6 ms in
  Chromium on an idle Apple M1 Pro, up to 150 ms on a loaded hosted runner.
  Three samples' spread alone (0 to 2 ms) once called choices measured 2 to
  8 ms faster inside one poll interval faster (culling for 120 cells, 67
  against 61 ms, macOS WebKit). The fix of
  that added the larger of the spread and the poll interval to the cheapest
  sample, which counted the resolution twice: a 150 ms poll called frames
  measured 174 ms against 324 ms "no faster" (Linux
  Chromium), an untrue tip. Each sample is now bounds (the GPU's share
  between the poll before and the poll that saw it done, less the page's
  share, which the clock measures exactly), the answer bounds the cheapest
  frame, the spread is only what the samples prove, and the verdict compares
  bounds: faster, no faster, or unresolved, said as such. Where the browser
  offers the GPU timer (Chromium on a hardware GPU), a frame's GPU time is
  exact and no verdict is unresolved. Too slow is the GPU's share alone over
  the budget: a frame slow on the page's thread is not made faster by
  drawing fewer cells. Still open, a product question: frames of every cell
  are measured right after an opening and the choice's frames later, so
  work of the opening that delays the polls widens the first bounds.
  Measuring both in the same window would remove that, at the price of a
  few frames drawn at full detail once the choice is on.
- **The browser suite on hosted runners** (2026-11-01): the page's own measurements made
  specs depend on the runner's GPU, which is why every spec runs on a fixed device, and a few
  pixels differ on hosted macOS and Windows Firefox that no local GPU shows. Both stories are in
  {doc}`38_decisions_testing_and_hosted_runners`.
- **Unchanged since the 2026-10-28 commit, still open:**
  - after a session restore, a kept view with labels on shows no centroid
    labels until labels are switched off and on;
  - on a phone, opening the panel while an analysis sheet covers 75–100%
    of the screen puts the panel underneath it; which one comes to the front
    is still to be decided.
- **Changes in flight:** colors and the legend were being changed on
  2026-10-30, and `data/cell-priority.js` changed that day (ranks by
  identity). Re-check the legend note text, the layer colors and the rank
  rule against the code before relying on this page for them.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| Alpha-byte layers | `cellucid/assets/js/rendering/alpha-visibility.js` | `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/alpha-word-kernels-contract.test.mjs` |
| Gray for zero (the switch, its floor and legend model) | `cellucid/assets/js/app/state/managers/color-manager.js` (`setContinuousZeroGray`, `updateColorsContinuous`), `cellucid/assets/js/app/state/domains/colors.js` (`zeroGray`), `cellucid/assets/js/app/actions/definitions/coloring.js` (`color.continuous.zero-gray`), `cellucid/assets/js/app/ui/modules/legend-renderer.js` | `cellucid/tests/legend-continuous-filter-independence.test.mjs`, `cellucid/tests/session-filter-restore-full-state-contract.test.mjs` (restored; a value without it is off), `cellucid/tests/state-domains-contract.test.mjs` (the colors value), `cellucid/tests/gene-sample-contract.test.mjs` (figure note) |
| Gray for 0, fade, read size | `cellucid/assets/js/rendering/read-cell-layer.js`, `cellucid/assets/js/app/state/core/constants.js` | `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/browser/gene-expression-layering.spec.mjs`, `cellucid/tests/browser/gene-sample.spec.mjs` |
| Color kernel (layer words) | `cellucid/assets/js/app/state/managers/continuous-colors.js`, `cellucid/assets/js/data/column-forms.js` | `cellucid/tests/continuous-color-table-contract.test.mjs`, `cellucid/tests/page-pipeline-contract.test.mjs` |
| Layer passes, read slot list | `cellucid/assets/js/rendering/high-perf-renderer.js`, `cellucid/assets/js/rendering/high-perf/read-cell-index.js`, `cellucid/assets/js/rendering/shaders/high-perf-shaders.js` | `cellucid/tests/high-perf-draw-order-contract.test.mjs`, `cellucid/tests/browser/render-blend-contract.spec.mjs` |
| Priority order and levels of detail | `cellucid/assets/js/data/cell-priority.js`, `cellucid/assets/js/rendering/high-perf/draw-order.js`, `cellucid/assets/js/rendering/high-perf/lod-levels.js`, `cellucid/assets/js/rendering/high-perf/store-order.js` | `cellucid/tests/lod-priority-order.test.mjs`, `cellucid/tests/browser/lod-priority-sample.spec.mjs`, `cellucid/tests/browser/lod-rank-order-resize.spec.mjs` |
| Picking at drawn size | `cellucid/assets/js/rendering/picking.js`, `cellucid/assets/js/rendering/viewer.js` | `cellucid/tests/picking-spatial-query.test.mjs` |
| Sprite size and cap | `cellucid/assets/js/rendering/point-sprite.js` (spliced into `shaders/high-perf-shaders.js` and `shaders/centroid-shaders.js`; read by `picking.js`, `read-cell-layer.js` and the figure export's `app/ui/modules/figure-export/utils/point-size.js` and `app/ui/modules/figure-export/utils/webgl-point-rasterizer.js`) | `cellucid/tests/point-sprite-contract.test.mjs`, `cellucid/tests/browser/point-near-camera-solid.spec.mjs` |
| Figure export layers | `cellucid/assets/js/app/ui/modules/figure-export/utils/lod-membership.js`, `cellucid/assets/js/app/ui/modules/figure-export/utils/webgl-point-rasterizer.js`, `cellucid/assets/js/app/ui/modules/figure-export/utils/point-projector.js`, `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-point-markup.js` | `cellucid/tests/figure-export-lod-membership-contract.test.mjs`, `cellucid/tests/figure-export-viewer-layer-fidelity.test.mjs`, `cellucid/tests/browser/figure-export.spec.mjs` |
| Smoke layer colors | `cellucid/assets/js/rendering/smoke-cloud/smoke-density.js`, `cellucid/assets/js/rendering/shaders/density-shaders.js` | `cellucid/tests/read-cell-layer.test.mjs` |
| Legend note, cell card | `cellucid/assets/js/data/gene-read-plan.js`, `cellucid/assets/js/app/ui/modules/legend-renderer.js`, `cellucid/assets/js/app/ui/modules/gene-read-control.js`, `cellucid/assets/js/app/ui/modules/cell-info.js` | `cellucid/tests/legend-renderer-terminal-ownership.test.mjs`, `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/gene-read-control-contract.test.mjs` |
| Automatic LOD and culling | `cellucid/assets/js/rendering/full-detail-frame-cost.js`, `cellucid/assets/js/rendering/viewer.js` (`armFullDetailFrameCost`, `frameCostOwesFrames`), `cellucid/assets/js/app/ui/modules/render-controls.js`, `cellucid/assets/js/app/ui/modules/benchmark/pass-profiler.js` (shares the GPU timer) | `cellucid/tests/full-detail-frame-cost.test.mjs`, `cellucid/tests/render-controls-exact-contract.test.mjs`, `cellucid/tests/benchmark-render-probes-contract.test.mjs`, `cellucid/tests/browser/full-detail-frame-cost.spec.mjs`, `cellucid/tests/browser/graphics-context-restore.spec.mjs` (verified again after a restore), `cellucid/tests/browser/helpers/device.mjs` (the simulated GPU timer every spec runs on) |
| Real-GPU launch for timings | `cellucid/scripts/lib/browser-launch.mjs`, `cellucid/assets/js/app/ui/modules/benchmark/renderer-identity.js` | `cellucid/tests/browser-launch-contract.test.mjs` |
| Pixel parity after settle | `cellucid/tests/browser/snapshot-compact-lod-parity.spec.mjs`, `cellucid/tests/browser/figure-export-multiview-parity.spec.mjs`, `cellucid/tests/browser/frustum-culling-parity.spec.mjs`, `cellucid/tests/browser/unplaced-cells.spec.mjs`, `cellucid/tests/browser/graphics-context-restore.spec.mjs` | (the specs themselves) |
| Settled includes every pending draw order | `cellucid/assets/js/rendering/high-perf-renderer.js` (`isPresentedGeometrySettled`) | `cellucid/tests/high-perf-draw-order-contract.test.mjs` (a lookup tree without culling or LOD) |
| Fence polls | `cellucid/assets/js/rendering/gl-utils.js` (`fenceHasPassed`): deferred error checks, the frame-cost measurement, `whenGpuIdle`, smoke noise batches, backdrop bakes | `cellucid/tests/webgl-deferred-error-check-contract.test.mjs`, `cellucid/tests/full-detail-frame-cost.test.mjs`, `cellucid/tests/smoke-noise-generation-transaction.test.mjs` (their fakes refuse an unflagged poll) |
| Point stores, segments, large arrays | `cellucid/assets/js/rendering/high-perf/point-store-segments.js`, `cellucid/assets/js/rendering/gl-utils.js` (`drawPointsInSegments`), `cellucid/assets/js/utils/large-arrays.js`, `cellucid/assets/js/rendering/linear-texture.js` | `cellucid/tests/browser/segmented-point-stores.spec.mjs`, `cellucid/tests/large-arrays-contract.test.mjs` |
| Adopting a tree's order | `cellucid/assets/js/rendering/high-perf/gpu-permute.js`, `cellucid/assets/js/rendering/high-perf/store-order.js`, `cellucid/assets/js/rendering/gl-errors.js` (`deferWebGLErrorCheck`) | `cellucid/tests/vertex-attribute-release-contract.test.mjs`, `cellucid/tests/webgl-deferred-error-check-contract.test.mjs`, `cellucid/tests/browser/lod-rank-order-resize.spec.mjs` |
| The frame loop and its pacing | `cellucid/assets/js/rendering/viewer.js`, `cellucid/assets/js/rendering/frame-phases.js`, `cellucid/assets/js/rendering/ambient-frame-pacing.js`, `cellucid/assets/js/rendering/pulse-layer-cache.js` | `cellucid/tests/frame-phases-contract.test.mjs`, `cellucid/tests/browser/frustum-culling-parity.spec.mjs` |
| Shader quality, fog, point-size curve | `cellucid/assets/js/rendering/shader-quality.js`, `cellucid/assets/js/rendering/shaders/fog-shaders.js`, `cellucid/assets/js/rendering/point-size-scale.js`, `cellucid/assets/js/rendering/shaders/high-perf-shaders.js` | `cellucid/tests/vertex-attribute-location-contract.test.mjs`, `cellucid/tests/high-perf-hidden-point-rejection-contract.test.mjs`, `cellucid/tests/browser/point-shader-no-depth-parity.spec.mjs`, `cellucid/tests/browser/sculpted-point-pass-state.spec.mjs` |

## Related pages

- {doc}`07_rendering_pipeline_webgl_and_performance_notes` (the renderer:
  passes, culling, frame cost, footguns)
- {doc}`12_figure_export_architecture`
- {doc}`19_subsampling_specification` (the uniform rank behind levels of
  detail) and {doc}`20_read_path_architecture`
- {doc}`23_decisions_partial_reading` (the gene sample and its size)
- {doc}`27_decisions_memory_model`
- {doc}`26_decisions_served_store_and_widgets`
- {doc}`01_codebase_map_and_entry_points`
- {doc}`../d_fields_coloring_legends/03_color_by_behavior` (what a user sees:
  the gray for 0, the fade, the legend note)
