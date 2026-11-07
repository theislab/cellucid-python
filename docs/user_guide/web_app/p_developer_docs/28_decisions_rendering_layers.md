# Decision record: rendering layers, the gene-sample look and levels of detail

This page answers: **"Why does a gene read for a sample of the open cells look
the way it does, why is the grey for 0 an option taken from the theme, how do the alpha
layers and levels of detail fit together, and how are rendering timings
taken?"** It records decisions, their evidence and what was replaced. The
renderer itself (buffers, passes, culling, frame cost) is described in
{doc}`07_rendering_pipeline_webgl_and_performance_notes`, and this page does
not repeat it. Update this page in the same change as the code it describes.

## At a glance

**Audience**
- Developers changing point colours, layers, sizes,
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

1. **A colour's alpha byte is its layer, never its visibility.** Visibility
   is the view's transparency (an R8 texture). A field's RGBA8 colour
   carries the layer in its alpha byte (`rendering/alpha-visibility.js`):
   - 255: coloured, on top;
   - `BENEATH_ALPHA_BYTE` = 254: the no-value grey (a missing value, a
     non-positive value on a log scale, and a gene's 0 when gray for zero is
     on), drawn beneath;
   - `UNREAD_ALPHA_BYTE` = 253: a cell the gene sample did not read, drawn
     at the bottom.

   Every consumer states what an alpha byte means:
   `requireColorAlphaMeaning` returns `'layer'` for cells and `'opacity'`
   for markers of their own (centroids). `colorLayersOf` reads a colours
   array once and reports whether it is layered and how many cells are not
   read.
2. **The colour kernel writes the layer.** `app/state/managers/continuous-colors.js`
   writes whole RGBA words. No value becomes `MISSING_WORD` (alpha 254), and
   a cell past the read rows becomes `UNREAD_WORD` (alpha 253). Both carry
   the RGB of the grey 211 for any reader that draws no layers. A sparse gene
   costs one word fill plus its stored entries (`data/column-forms.js`).
3. **Passes, bottom first** (`rendering/high-perf-renderer.js`): one pass for
   a single layer; with layers, the grey (pass 1) and then the coloured cells
   (pass 2). When some cells are not read, those come first (pass 3):
   - flat (no lighting), at the level's size, without depth writes, in
     `u_unreadColor`;
   - then the grey in `u_zeroColor`;
   - then the coloured cells at the read size. In 3D the grey pass yields
     in depth by one point diameter (`u_yieldDepth`).

   Each pass owns its blend state.
4. **Gray for zero is an option, off by default.** A gene's 0 is a value: it
   takes the colour scale's lowest colour, opaque. The legend's **Gray for
   zero** switch (`DataState.setContinuousZeroGray`, the field's `_zeroGray`,
   defaulted to false by `ensureContinuousMetadata`) gives the kernel a
   floor at 0 (`neutralAtOrBelow: 0` in `updateColorsContinuous`), so a
   gene's zeros take `MISSING_WORD` and are drawn beneath. It applies to a
   gene with no negative value only (`_zeroGrayApplies`; an obs field or
   scaled data has no floor, and the switch is hidden); on a log scale a 0
   has no logarithm and is the grey either way (the switch is disabled). The
   legend model says it as `zeroGray: {enabled, applies, drawn}` (`drawn`:
   the zeros are grey now). It is per field like the log scale: kept views
   clone it (`_cloneFieldState`), the var axis carries it by gene
   (`GENE_INTENT_KEYS`), a duplicated field keeps it, and the switch is the
   action `color.continuous.zero-gray`, one history step. It is a colour of
   the field, held by the `colors` state domain (`state/domains/colors.js`:
   `zeroGray: true` in the field's entry when on, absent when off, since a
   value records only what differs from the field's opening state), so a step,
   an undo and a session file hold it. Python and R read a session's `colors`
   value only as data (an apply writes the `cells`, `field-overlays`,
   `highlights` and `gene-sets` values to an object and reports `colors` as
   presentation, not applied), and the shared golden history
   (`history-golden.cellucid-session`) holds no `zeroGray`. Everything that
   draws a colours array follows from the alpha byte alone (views, smoke,
   figures, the cell card). A figure whose zeros are grey says `Gray: 0`
   under its colour bar (`buildColorScaleNote`).
5. **The grey for 0 is made from the background** (`rendering/read-cell-layer.js`,
   `zeroCellColor`). It is a neutral grey `ZERO_LIGHTNESS_CONTRAST`
   (15.44 CIE L\*, the step of grey 211 from white) away from the
   background's lightness, toward the middle. On white that is 211, as
   before. On the dark grid it is 56, and on black 38.
6. **The fade** (`unreadCellColor`). A cell not read is drawn in the
   background moved `UNREAD_LIGHTNESS_STEP` (8 L\*) toward the grey for 0,
   but never more than half the way. The map keeps its shape and recedes,
   and a cell not read never looks like a 0.
7. **The read size** (`readCellSizeMultiplier`). Each frame projects up to
   8,192 map cells and as many read cells, and measures the map's on-screen
   area with an occupancy grid. It then picks the diameter at which the
   shown read cells cover `READ_CELLS_COVERED_FRACTION` (40 %) of that area.
   The multiplier stays between 1 and √(drawn ÷ read): a sample of every
   cell draws at the normal size.
8. **Read cells at every level of detail** (`rendering/high-perf/read-cell-index.js`).
   When the level draws at least as many cells as were read, the read cells
   are drawn from their own element buffer of slots (4 B per read cell),
   using `drawElements` over the store's own vertex arrays. Nothing is
   copied. When more cells were read than the level draws, the level draws
   its own share of them by rank, and the size multiplier is computed with
   that share (`readShare`).
9. **One set of layer colours for every surface.** `cellLayerColors` and
   `drawnCellRgb` are used by:
   - the views (`setLayerBackground`, and the CSS variables
     `--cell-zero-color` and `--cell-unread-color` that the legend dots use);
   - smoke, whose density splat substitutes the layer colours for bytes 254
     and 253 and rebuilds when the background changes
     (`rendering/smoke-cloud/smoke-density.js`, `rendering/shaders/density-shaders.js`);
   - figure export, which takes its colours from the figure's own paper
     (`FigurePointLayers` in
     `app/ui/modules/figure-export/utils/lod-membership.js`; `drawnSequence`
     orders the layers);
   - the cell info card's swatch (`app/ui/modules/cell-info.js`).

   The legend and the figure share one note, `geneSampleLayersNote`
   (`data/gene-read-plan.js`): "Faded: not read", and "Faded: not read ·
   gray: 0" when the zeros are grey (`zeroGray.drawn`). The legend writes it
   in one line beside a dot of each colour, under the colour bar, only for a
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

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-21 | GPU-heavy work adapts to the measured GPU speed; macOS Firefox CI GL is about 2,500× slower than Metal | No browser preference changes it; giving a spec more time hides the cause | Fixed timeouts; one GPU budget for every machine |
| 2026-10-21 | Pixel-parity specs copy cameras only after the scene settles (`getFrameScheduleState().settled`), assert camera equality, and never use a pixel tolerance | A camera copied mid-ease differs by a few edge pixels on slow GL | Tolerances that hide which input differed |
| 2026-10-28 | Non-expressing cells are drawn beneath expressing ones, in a neutral colour | Sparse expression must not be covered by the grey around it | One pass with arbitrary overlap |
| 2026-10-28 | A colour's alpha byte is its layer; filters write only the transparency | A filter never copies or re-uploads colours a kept view shares | Alpha as visibility |
| 2026-10-28 | Levels of detail are fixed counts in priority order; ranks by count, with no admission table | A level is a uniform sample, and growth is in place (per band 10–20 ms and 48 MB → 0 at 10M; 66–79 ms and 229 MB → 0 at 50M) | Spatially chosen levels; per-dataset ladders |
| 2026-10-28 | Levels of detail and culling turn on automatically only when measured too slow, and are kept only when measured faster | Small data never needs them; a 2D full framing got slower with culling | Turning them on by cell count or by dimension |
| 2026-10-29 | Gene on a sample: unread cells fade almost to the background, distinct from the grey for 0; read cells drawn on top at a size set by the sample's on-screen density; no extra downloads | The expression pattern must read at once on a map of all cells | Drawing unread cells as no-value grey; reading more to fill the map |
| 2026-10-30 | The grey for 0 and missing values comes from the theme (fixed L\* contrast) | On a dark background grey 211 outshone the darker colormap cells (seen on a real atlas) | One fixed grey on every theme |
| 2026-11-07 | Ambient animations (smoke, flowing velocity particles, backdrop motion, the highlight pulse) are paced by how late the page's own frames run behind the GPU, confirmed by a fence read only when a frame is due, with the display interval measured | On a slow GPU (CI software GL; hosted macOS Firefox) an endless animation saturated the GPU and every click waited behind it (overlay on: 17.7 page fps under SwiftShader, clicks of 10–15 s on hosted macOS Firefox); a GPU that keeps up is never held back | Pacing every frame by fence latency, which throttled fast GPUs |
| 2026-10-31 | Gray for zero is a per-field switch, off by default; off, a gene's 0 takes the scale's lowest colour (as before 2026-10-28); on, the theme's grey beneath. The verbose legend note under the colour bar is removed; a gene sample keeps one line naming its states | Owner, on HLCA and Suo: a 0 is a value, and the grey should be a choice, explained by an (i), not a fixed rule and a paragraph | Gray for 0 always on; a note under every gene's colour bar |
| 2026-10-30 | A level ranks a cell by its identity under the share's seed, whatever order it is open in | The same cells come first in every order, step, view, revisit and session | Ranking by open position |
| 2026-10-30 | Rendering and app timings only on a real GPU (Chromium `--use-angle=metal` on macOS); software renderers are refused unless asked for (`--software-gl`, never a timing) | SwiftShader draws on the CPU and holds the main thread on synchronous GL calls | Any timing from SwiftShader or llvmpipe |
| 2026-11-01 | A store's draw order is an input of the picture. A view is not settled (`isPresentedGeometrySettled`) while a tree for its positions is being built or is published but not adopted, a cell lookup's tree included; a pixel-parity spec draws the frames it compares in one order (or chooses level of detail and culling by hand) and records it | `unplaced-cells.spec.mjs` drew its reference in a lookup tree's order (the cell info card's idle request ran after the spec changed the positions) and the unplaced frame in source order: 82,845 pixels apart on Linux Chromium and 83,440 on Linux Firefox (push run 41292953076), reproduced exactly on SwiftShader by building the tree for the reference. `graphics-context-restore.spec.mjs` captured before the automatic choice turned culling on and after it (3 pixels, macOS Firefox, nightly run 41269083996) | Comparing frames drawn in different orders; treating a lookup tree as changing nothing on screen |
| 2026-11-01 | The page's own automatic level of detail and culling are verified again whenever their measurement basis changes (a restored context, cells grown or cut, another representation, a new framing), exactly as at an opening: off, frames of every cell measured, on, verified, kept or off. One path: the viewer says `measuring` before each measurement and the render controls turn their choice off then (`applyDatasetPerformance`, which did it for a new dataset only, is removed). A person's or a caller's choice is left alone | Kept across a re-measure, the choice was never checked on the new basis: after a restore its measurement measured frames with culling on as frames of every cell and turned nothing new on, so nothing verified it (found in `graphics-context-restore.spec.mjs` on macOS Firefox, nightly run 41269083996); a re-measure with level of detail on was not measured at all | Keeping an unverified choice; a restore-only special case |
| 2026-11-01 | A frame-cost sample is bounds, not a number: the page's share (exact, from the clock) and the GPU's (exact from the GPU timer where the browser offers one, else between two fence polls). Too slow is the GPU's share proven over the budget; the verdict on the automatic choice is `faster`, `not-faster` or `unresolved` from the bounds and the proven spread, and the tip says which | The rule that added the larger of the spread and the poll interval to the cheapest sample counted the resolution twice and called a frame measured 150 ms faster "no faster" (run 41307458156, Linux Chromium); a frame slow on the page's thread is not made faster by drawing fewer cells | A noise constant; adding an instrument's resolution to a bound that already holds it; a tip that says "no faster" where the numbers show faster |
| 2026-11-01 | A measurement is armed only once the state it measures is committed (a publication of cells or positions, a restored context), asks for its frame at once, and keeps the scene unsettled until it has answered (`frameCostOwesFrames`) | A publication that failed had armed it first and its `measuring` report turned the page's own culling off (`camera-defaults.spec.mjs`, run 41307458156); a measurement armed for the next framing of a still scene waited for whatever drew next and could turn level of detail on under it | Arming before a commit; a measurement that answers at a time something else decides |
| 2026-11-01 | Every browser spec runs on a fixed device (`tests/browser/helpers/device.mjs`, the `device` option of the `test` fixture): a simulated GPU timer (1 ms a frame by default), a held opening cost (the machine's `shownCellMs`, `utils/machine.js`; nothing measured by default), 4 cores, 8 GiB where reported. A spec about a measured choice sets its own; a spec about the real GPU's time says `gpu: 'real'` | Five CI rounds each patched one spec that the runner's own measurements changed under it (`drawEveryCellByHand`, `drawEveryCellFromTheStart`, zooming in first, waiting for the start sample's record); the cause was one | Per-spec holds of the automatic choices; specs that pass or fail by the runner's GPU |
| 2026-11-01 | Every fence poll carries `SYNC_FLUSH_COMMANDS_BIT` (`fenceHasPassed`, `rendering/gl-utils.js`) | Firefox warns at the hundredth poll of an unsignalled fence made without it, whatever `flush` followed the fence; a GPU that takes seconds a frame reaches that count (`renderer-smoke.spec.mjs` velocity slider, macOS Firefox, nightly run 41269083996); reproduced on an Apple M1 Pro with slow fragment work, where the flagged poll stays silent | `getSyncParameter(SYNC_STATUS)` or an unflagged `clientWaitSync` as a poll |

## Measurements that justify it

- **Ambient pacing** (2026-11-07, Apple M1 Pro; page fps is the page's own rAF rate): SwiftShader with the velocity overlay, page 17.3–17.9 → 34–38 fps, overlay 17.9 → 17.5–17.8 fps, clicks 111–148 → 30–121 ms; SwiftShader smoke, page 25 → 41 fps; default headless Chromium overlay, page 17.6 → 45 fps, clicks 650 → 22–94 ms; native Metal (Chromium 60, Firefox 120, WebKit 60 fps) unchanged.
- **Real atlas walk** (2026-10-30, Chromium on Metal, Apple M1 Pro, HLCA
  full, 2,282,447 cells): opening 12.9–17.0 s, a gene 4.6–5.4 s on a sample
  of 1,759–2,273 cells, read-size multiplier 8.6–10.35 at a base point of
  0.32, and a PNG figure with the three layers in 0.77–0.8 s. On the light
  theme the pattern showed at once. On the dark theme grey 211 outshone the
  darker viridis cells, which led to the 2026-10-30 grey decision. The
  rerun had 0 page errors.
- **The zero-grey rule** is checked in `cellucid/tests/read-cell-layer.test.mjs`.
  It gives a fixed contrast on every background, 211 on white and a
  mid-dark grey on dark, and on every theme, backdrop looks included, at
  least half of every colormap stands out more than the grey does.
- **Growth without re-allocation** (2026-10-30, real-GPU page profile at
  10M). Each step of a progressive opening used to rebuild the alpha texture
  and read GPU errors back synchronously, about 250 ms of page time per
  step. Stores and textures are now sized for the growth target and kept.
- **Page colour cost on Metal** (10M): a gene shown in about 90 ms before
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

## Tried and reverted

- **Pacing every ambient frame by a GPU fence's latency** (2026-11-07). A fence measures latency, not occupancy, so it held back GPUs that keep up: smoke on Metal went from 60 to 45 fps. Replaced by pacing on the page's own lateness, confirmed by the fence only when a frame is due.
- **Unread cells as the no-value grey** (design of 2026-10-28). Cells not
  read were NaN and drawn in the grey beneath, with the difference stated
  only in text. On a large map a 0 and "not read" looked the same. It was
  replaced on 2026-10-29 by the fade and the larger read cells.
- **One grey (211) on every theme.** It was replaced on 2026-10-30 by
  `zeroCellColor`. On light themes the grey is unchanged.
- **`hasBeneathCells`** was replaced by `colorLayersOf` when the third layer
  arrived.
- **Ranking a share's cells by open position.** Under the default seed, the
  cells drawn first changed with the order the cells were open in. Since
  2026-10-30 ranks follow identity.
- **An earlier fixed ladder of 18 reductions** was replaced by fixed counts
  in priority order. The coarsest reduction, 44, survives as the farthest
  camera's size cap.
- **Automatic culling on a 2D full framing** made frames slower (above). It
  led to the rule that an automatic choice is verified and turned off again
  when it is no faster.
- **A growth fix counted capacity as held bytes.** `PointStore#byteLength`
  briefly reported capacity, so recolouring was refused during growth. It
  now reports the cells held, and `allocatedByteLength` reports the storage.
- **A layering spec set up wrongly.** Its "2D" view still had depth and its
  3D sheet faced away at an angle. The spec now uses `setViewEmbedding` and a
  sheet facing the camera. The lesson: prove the inputs of a pixel
  comparison before trusting it.

## Known limits and open items

- **Highlight rings** are drawn at the level's point size, not the read
  size, so a ring around an enlarged read cell is too small
  (`rendering/highlight-renderer.js` sizes rings by
  `getCurrentLODSizeMultiplier` only).
- **Smoke** substitutes layer colours, but a volume has no layers or sizes:
  read cells are not enlarged there.
- **The zero-grey walk on a real atlas** (light, dark, zoom, coarse level,
  kept view, figure) has not been done. The first attempt stalled while the
  atlas was still opening, before the rendering code ran.
- **Planned renderer work:**
  - draw PNG and Hybrid figures from the GPU stores in the page's context,
    so they get layers, LOD and read cells for free;
  - let kept views share the live colour store and R8 alpha (copy-on-write);
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
  against 61 ms, macOS WebKit, Web contracts run 41255961977). The fix of
  that added the larger of the spread and the poll interval to the cheapest
  sample, which counted the resolution twice: a 150 ms poll called frames
  measured 174 ms against 324 ms "no faster" (run 41307458156, Linux
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
- **The page's own measurements made the browser suite depend on the
  runner** (2026-11-01). Five CI rounds (Web contracts runs 41213981681 to
  41307458156) each patched one spec that saw the page change state under it
  because a runner's software GL measured otherwise than a developer's Mac:
  level of detail and culling turned on, verified and off at times the
  runner decided; a failed publication turned the page's own culling off
  (it armed the measurement before the publication committed, now after);
  the start sample's opening-cost record made an "at once" opening
  progressive. Every spec now runs on a fixed device
  (`tests/browser/helpers/device.mjs`: a simulated GPU timer, a held opening
  cost, fixed cores and memory); specs about a measured choice set their
  own. See page 14, "What the hosted CI runners differ in".
- **Pixel differences on CI that no local GPU shows** (2026-11-01, not
  reproduced, causes unproven):
  - `figure-export.spec.mjs` "the point pass drawn on the page equals the
    one drawn in the worker": on macOS Firefox (nightly run 41269083996) the
    two rasters differed at 7 of 15,346 drawn pixels, by at most 3 levels,
    all on antialiased rims, alpha included (so not a premultiplication
    round trip: both bitmaps are read through the same 2D canvas). Both are
    drawn by the same code with the same requested attributes in two WebGL2
    contexts (an `OffscreenCanvas` on the page and one in the worker), so
    the two contexts are the only input left. The spec now compares both
    contexts' descriptions (`rasterContext`: renderer, vendor, granted
    attributes, samples) before the pixels and attaches them with the first
    32 differing pixels: the next failure names the input. If the contexts
    differ there, the contract is one context's pixels, not the same pixels
    from any two.
  - `segmented-point-stores.spec.mjs` "flatumap", reduced level: on Windows
    Firefox (same run) 2 pixels, both on row 359, differed between one store
    and segments of 4,096 points, with store bytes, cameras, level and cells
    drawn equal. The spec now attaches the context (renderer, granted
    attributes, samples) and every draw of both frames (program, mode,
    first, count, vertex array, viewport) of a configuration that differs.
    On macOS Firefox both patterns timed out (600 s, the slow GL above).
- **Unchanged since the 2026-10-28 commit, still open:**
  - after a session restore, a kept view with labels on shows no centroid
    labels until labels are switched off and on;
  - on a phone, opening the panel while an analysis sheet covers 75–100 %
    of the screen puts the panel underneath it; which one comes to the front
    is still to be decided.
- **Changes in flight:** colours and the legend were being changed on
  2026-10-30, and `data/cell-priority.js` changed that day (ranks by
  identity). Re-check the legend note text, the layer colours and the rank
  rule against the code before relying on this page for them.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| Alpha-byte layers | `cellucid/assets/js/rendering/alpha-visibility.js` | `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/alpha-word-kernels-contract.test.mjs` |
| Gray for zero (the switch, its floor and legend model) | `cellucid/assets/js/app/state/managers/color-manager.js` (`setContinuousZeroGray`, `updateColorsContinuous`), `cellucid/assets/js/app/state/domains/colors.js` (`zeroGray`), `cellucid/assets/js/app/actions/definitions/coloring.js` (`color.continuous.zero-gray`), `cellucid/assets/js/app/ui/modules/legend-renderer.js` | `cellucid/tests/legend-continuous-filter-independence.test.mjs`, `cellucid/tests/session-filter-restore-full-state-contract.test.mjs` (restored; a value without it is off), `cellucid/tests/state-domains-contract.test.mjs` (the colours value), `cellucid/tests/gene-sample-contract.test.mjs` (figure note) |
| Grey for 0, fade, read size | `cellucid/assets/js/rendering/read-cell-layer.js`, `cellucid/assets/js/app/state/core/constants.js` | `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/browser/gene-expression-layering.spec.mjs`, `cellucid/tests/browser/gene-sample.spec.mjs` |
| Colour kernel (layer words) | `cellucid/assets/js/app/state/managers/continuous-colors.js`, `cellucid/assets/js/data/column-forms.js` | `cellucid/tests/continuous-color-table-contract.test.mjs`, `cellucid/tests/page-pipeline-contract.test.mjs` |
| Layer passes, read slot list | `cellucid/assets/js/rendering/high-perf-renderer.js`, `cellucid/assets/js/rendering/high-perf/read-cell-index.js`, `cellucid/assets/js/rendering/shaders/high-perf-shaders.js` | `cellucid/tests/high-perf-draw-order-contract.test.mjs`, `cellucid/tests/browser/render-blend-contract.spec.mjs` |
| Priority order and levels of detail | `cellucid/assets/js/data/cell-priority.js`, `cellucid/assets/js/rendering/high-perf/draw-order.js`, `cellucid/assets/js/rendering/high-perf/lod-levels.js`, `cellucid/assets/js/rendering/high-perf/store-order.js` | `cellucid/tests/lod-priority-order.test.mjs`, `cellucid/tests/browser/lod-priority-sample.spec.mjs`, `cellucid/tests/browser/lod-rank-order-resize.spec.mjs` |
| Picking at drawn size | `cellucid/assets/js/rendering/picking.js`, `cellucid/assets/js/rendering/viewer.js` | `cellucid/tests/picking-spatial-query.test.mjs` |
| Figure export layers | `cellucid/assets/js/app/ui/modules/figure-export/utils/lod-membership.js`, `cellucid/assets/js/app/ui/modules/figure-export/utils/webgl-point-rasterizer.js`, `cellucid/assets/js/app/ui/modules/figure-export/utils/point-projector.js`, `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-point-markup.js` | `cellucid/tests/figure-export-lod-membership-contract.test.mjs`, `cellucid/tests/figure-export-viewer-layer-fidelity.test.mjs`, `cellucid/tests/browser/figure-export.spec.mjs` |
| Smoke layer colours | `cellucid/assets/js/rendering/smoke-cloud/smoke-density.js`, `cellucid/assets/js/rendering/shaders/density-shaders.js` | `cellucid/tests/read-cell-layer.test.mjs` |
| Legend note, cell card | `cellucid/assets/js/data/gene-read-plan.js`, `cellucid/assets/js/app/ui/modules/legend-renderer.js`, `cellucid/assets/js/app/ui/modules/gene-read-control.js`, `cellucid/assets/js/app/ui/modules/cell-info.js` | `cellucid/tests/legend-renderer-terminal-ownership.test.mjs`, `cellucid/tests/read-cell-layer.test.mjs`, `cellucid/tests/gene-read-control-contract.test.mjs` |
| Automatic LOD and culling | `cellucid/assets/js/rendering/full-detail-frame-cost.js`, `cellucid/assets/js/rendering/viewer.js` (`armFullDetailFrameCost`, `frameCostOwesFrames`), `cellucid/assets/js/app/ui/modules/render-controls.js`, `cellucid/assets/js/app/ui/modules/benchmark/pass-profiler.js` (shares the GPU timer) | `cellucid/tests/full-detail-frame-cost.test.mjs`, `cellucid/tests/render-controls-exact-contract.test.mjs`, `cellucid/tests/benchmark-render-probes-contract.test.mjs`, `cellucid/tests/browser/full-detail-frame-cost.spec.mjs`, `cellucid/tests/browser/graphics-context-restore.spec.mjs` (verified again after a restore), `cellucid/tests/browser/helpers/device.mjs` (the simulated GPU timer every spec runs on) |
| Real-GPU launch for timings | `cellucid/scripts/lib/browser-launch.mjs`, `cellucid/assets/js/app/ui/modules/benchmark/renderer-identity.js` | `cellucid/tests/browser-launch-contract.test.mjs` |
| Pixel parity after settle | `cellucid/tests/browser/snapshot-compact-lod-parity.spec.mjs`, `cellucid/tests/browser/figure-export-multiview-parity.spec.mjs`, `cellucid/tests/browser/frustum-culling-parity.spec.mjs`, `cellucid/tests/browser/unplaced-cells.spec.mjs`, `cellucid/tests/browser/graphics-context-restore.spec.mjs` | (the specs themselves) |
| Settled includes every pending draw order | `cellucid/assets/js/rendering/high-perf-renderer.js` (`isPresentedGeometrySettled`) | `cellucid/tests/high-perf-draw-order-contract.test.mjs` (a lookup tree without culling or LOD) |
| Fence polls | `cellucid/assets/js/rendering/gl-utils.js` (`fenceHasPassed`): deferred error checks, the frame-cost measurement, `whenGpuIdle`, smoke noise batches, backdrop bakes | `cellucid/tests/webgl-deferred-error-check-contract.test.mjs`, `cellucid/tests/full-detail-frame-cost.test.mjs`, `cellucid/tests/smoke-noise-generation-transaction.test.mjs` (their fakes refuse an unflagged poll) |

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
  the grey for 0, the fade, the legend note)
