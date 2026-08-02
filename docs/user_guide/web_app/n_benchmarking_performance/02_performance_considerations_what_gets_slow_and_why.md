# Performance considerations (what gets slow and why)

**Audience:** computational users + power users (still readable for everyone)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- Which actions are “hot paths” (and why they get slow)
- The biggest multipliers (`n_cells`, `n_views`, pixels, categories)
- The highest-leverage knobs in the UI (and which bottleneck they target)
- Safe performance workflows that preserve scientific intent

**Prerequisites:**
- A dataset loaded (small is fine; large helps you *feel* the differences)

---

## Start here: “what is the expensive thing?”

Cellucid slowdowns usually come from one of these:

- **Compute something for every cell** (CPU, often `O(n_cells)` or worse),
- **Draw too much every frame** (GPU, often scales with pixels and views),
- **Load too much data** (I/O, often scales with bytes requested and latency).

This page maps *common user actions* to their likely bottleneck and the most effective knob.

If you haven’t read the mental model yet, start with {doc}`01_performance_mental_model`.

---

## The “cost model” by interaction (what scales with what)

### Loading data (dataset + fields + genes)

**What tends to scale with bytes:**
- initial dataset load (embeddings + metadata),
- loading gene expression chunks on demand,
- remote server latency (each request has overhead even if data is small).

**Big practical reality:**
- loading a large `.h5ad` directly in the browser is not truly lazy and can freeze/crash.

Best practices and alternatives:
- {doc}`../b_data_loading/01_loading_options_overview`
- {doc}`../b_data_loading/03_browser_file_picker_tutorial` (browser-only limits)
- {doc}`../b_data_loading/04_server_tutorial` (recommended for large files)

---

### Rendering points (navigation smoothness)

**What tends to scale with GPU work per frame:**
- number of points *actually drawn* (and whether LOD/frustum culling reduces it),
- point size, because it sets how many pixels each sprite covers,
- window size and device pixel ratio (retina screens are expensive),
- antialiasing, which multiplies the cost of the drawing buffer itself.

Three things that are commonly assumed to cost per frame, and do not:

- **Hidden points.** A cell that filters have hidden is rejected before
  rasterisation — it is pushed outside the clip volume in the vertex shader, so
  no fragment is ever produced for it. A heavily filtered view is therefore
  *cheaper* to navigate than the same dataset unfiltered, not more expensive.
- **Standing still.** When nothing changes, the frame transfers no data to the
  GPU at all. Holding a view steady does not accumulate cost.
- **Shader quality.** `Full`, `Light` and `Ultra-light` measured the same on the
  machine below, at every point size tested. The setting changes what a point
  looks like, not how fast the frame draws — see the next section.

High-leverage knobs:
- **Make the window smaller** while exploring (fewer pixels) — *usually the biggest immediate win*.
- **Filter down to what you are actually looking at** — hidden points stop costing anything to draw.
- **Enable Level-of-Detail (LOD)** (points mode) — reduces draw cost when zoomed out.
- **Reduce `Point size (log):`** once points are big enough to overlap — the one
  rendering setting that measurably moves frame time, and it is measured below.
- **Turn off `Antialiasing (smooth point edges)`** under **Visualization →
  Image quality** — measured below, and the largest single lever left on the
  plain path once filtering and window size have been used.
- **Change the layout shape, not just the view count** — a two- or three-view
  row is the expensive arrangement; one view and the 2×2 four-view grid shade
  about the same number of pixels. See {doc}`01_performance_mental_model`.

#### Where the frame time actually goes, and why shader quality is not a lever

All the figures in this section and the next were measured on **one machine** —
an Apple M1 Pro driving WebGL2 through ANGLE Metal, 1440×1000 at device pixel
ratio 1, 10,000,000 synthetic points, LOD and frustum culling off, camera held
fixed. Treat the absolute milliseconds as belonging to that machine. What
carries across hardware is the *shape*: which term scales with what.

Split by GPU timer query, at the shipped default point size:

| Stage | Cost | Share of the frame |
|---|---|---|
| Vertex processing | 2.62 ms, and **flat** — identical at every point size | 4.2% |
| Fragment shading | rises 1.82× for **190×** the fragments | small, and mostly hidden |
| Per-sprite rasterisation | the remainder | the bulk of the frame |

Two consequences follow, and both contradict advice that used to appear on this
page:

- **Simplifying the fragment shader buys nothing.** Depth writes are on and
  submission order is uncorrelated with depth, so at the default point size only
  about **5%** of rasterised fragments survive to be shaded at all. Removing fog,
  removing lighting, removing the alpha texture fetch, removing the round
  discard, and stepping the whole quality ladder down — 43 shader statements to 9
  to 2 in the compiled Metal — each stayed inside its own measurement noise band,
  at every point size tested. The fragment shader runs on a rounding error's
  worth of the fill, so making it cheaper changes nothing you can see.
- **The floor is the per-sprite rasterisation footprint**, not shading and not
  vertex work: binning, rasteriser setup and raster operations, measured at about
  **5.7 ns per sprite (~175 million sprites/s)** and reproduced three times by
  three independent instruments. `GL_POINTS` is the cheapest way to reach it —
  instanced quads and every triangle arrangement tried came out **1.41× to 2.12×
  slower**, and the quad arms were shading a fragment count matched to the point
  arm to five decimal places, so the gap is the primitive and not extra work.
  Halving the primitive count with one triangle per point made it *worse*, so
  there is no arrangement left to try. At 10 million points that floor alone
  accounts for most of a ~50 ms frame, which is the ~19 FPS in the antialiasing
  table below.

:::{note}
The **Analyze Performance** tool used to report a "Shader overhead" figure and
suggest switching to `Ultra-light`. Both are gone, and the sidebar's
`Shader quality:` tooltip no longer claims a speed benefit either — it describes
what the setting changes about how points *look*. The tool reports point-size
response in that row now. Reach for point size, window size, filtering or
antialiasing.
:::

#### Point size: the control that does respond

Point size is the one rendering setting with a measurable effect on frame time,
and the reason is geometric rather than arithmetic: a sprite's area grows with
the square of its size, so it changes how much the rasteriser has to cover
rather than how much arithmetic runs per pixel.

It does not respond everywhere on its range. Below roughly **`gl_PointSize` 15**
— about the middle of the `Point size (log):` slider, or the default size after
zooming in about sixfold — the frame is bound by the per-sprite cost, and
shrinking points further buys very little. Above it the frame becomes fill-bound
and point size becomes the dominant term. In practice:

- If points already overlap heavily on screen, **reduce point size first**.
- If points are tiny specks and the frame is still slow, point size is not your
  problem — you are paying per sprite, so **draw fewer sprites** (filter, or
  enable LOD).

Both thresholds are hardware rates, so the crossover position moves with the
GPU; the ordering does not.

#### Antialiasing: what it costs and what it buys

Antialiasing is on by default and can be turned off under **Visualization →
Image quality**. It is not a cosmetic setting, so both halves of the trade are
stated here.

Measured on the same machine as above — 10,000,000 synthetic points, one Apple
M1 Pro through ANGLE Metal, 1440×1000 at device pixel ratio 1, camera held
fixed, four independent browser windows per configuration — and on the frames
actually presented to the screen, not in an offscreen buffer:

| Point size | Frame time, antialiasing on | off | Change |
|---|---|---|---|
| 0.75 (the default) | 52.7 ms (19.0 FPS) | 42.5 ms (23.6 FPS) | **19% faster** |
| 2.0 | 69.7 ms | 48.9 ms | **30% faster** |
| 4.0 | 93.1 ms | 61.6 ms | **34% faster** |

A GPU timer query over the same runs agreed within three points at every size
(22%, 32%, 36%), so read the saving as **about a fifth at the default point size
and about a third at large ones**. Note that the saving *grows* with point size:
antialiasing costs per pixel covered, so it is worth least in exactly the
millions-of-tiny-dots case that most needs the frames.

What you give up, on the same frames: **18% of pixels change at the default
point size** and 32% with `Ultra-light (square points)`. Small dots are almost
entirely edge, so this is the setting that changes what a dense cloud looks
like, not merely how fast it is drawn. Larger points are affected far less,
which is the opposite of the speed table — at point size 4 the picture barely
changes while the saving is largest.

The browser fixes antialiasing when it creates the drawing buffer, so the
setting cannot take effect while the page is open. The panel says so, and
saves your choice for the next load.

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`
- {doc}`../a_orientation/02_system_requirements` (WebGL context lost = GPU memory pressure)

---

### Volumetric “smoke” mode (cinematic density rendering)

Smoke mode is intentionally heavy and can dominate GPU time.

Separate build cost from frame cost:

- **Visible-cell count** scales density-build validation and GPU submission
  approximately linearly.
- **Grid density** scales density-build work and GPU memory steeply—roughly
  with the cube of the selected grid width.
- **Render resolution** and **Ray quality** are the main steady-state
  per-frame costs because they control rendered pixels and ray-march work.

Practical consequences:
- The exact grid choices are 32³, 48³, 64³, 96³, and 128³; 128³ is the
  initial and maximum setting on every supported WebGL2 implementation.
- Smoke can look “blank” if density is low or if few points are visible (filters/outliers).
- Smoke can cause “context lost” if grid/quality are too high for your GPU.

---

### Filtering (visibility recomputation)

Filtering is a classic CPU hot path — and the CPU half is the *whole* cost.

Separate the two halves, because they scale completely differently:

- **Deciding what is visible (CPU).** Every filter change re-derives visibility
  for every cell, testing each enabled filter. This is roughly
  **`O(n_cells × n_enabled_filters)`** and it does **not** get cheaper when the
  filter hides almost everything — narrowing a range to eight cells costs the
  same pass as narrowing it to eight hundred thousand.
- **Telling the GPU (upload).** Only the parts of the visibility texture whose
  bytes actually moved are re-uploaded, in at most a few dozen contiguous
  regions. A change that touches one small region of the dataset moves
  kilobytes; a change scattered across the whole dataset touches every region
  and costs a full re-upload. A filter edit that changes nothing at all uploads
  nothing and does not disturb the highlight buffers either.

So the upload is rarely what you feel. The thing you feel is the per-cell
recomputation, repeated once per slider frame.

The most common performance mistake is **scrubbing sliders** while Live filtering is on.

Best practices:
- for continuous filters: turn **Live filtering** off and click **FILTER** once,
- avoid stacking many enabled filters if you can disable/remove no-op filters,
- tune filters in a single view first, then create snapshots after you stabilize the filter state.

Related docs:
- {doc}`../e_filtering/05_performance_considerations`
- {doc}`../d_fields_coloring_legends/03_color_by_behavior` (Live filtering vs FILTER)

---

### Highlighting and selection (large sets)

Highlighting is often cheap for moderate sizes, but there are “cliffs”:

- huge selections/highlight groups (hundreds of thousands → millions of cells),
- many highlight groups/pages that must be merged/resolved per draw,
- frequent selection updates (dragging lasso over millions of points).

What it feels like:
- selection tools lag while dragging,
- highlight updates stutter,
- sessions get very large (saving/restoring takes longer).

Related docs:
- {doc}`../f_highlighting_selection/05_edge_cases_highlighting` (large highlights)
- {doc}`../l_sessions_sharing/index` (sessions size and restore behavior)

---

### Analysis (group computations)

Analysis cost depends on the mode, but the big drivers are usually:

- group sizes (how many cells are in A and B),
- the number of features involved (e.g., genes),
- repeated reruns caused by changing inputs frequently.

Practical workflow:
- define groups carefully (avoid “everything vs everything” when you don’t need it),
- keep analysis scope small while iterating,
- export results for offline analysis when needed.

Related docs:
- {doc}`../h_analysis/01_analysis_mental_model`
- {doc}`../h_analysis/10_troubleshooting_analysis` (includes performance symptoms)

---

### Figure export (large artifacts)

Export can be expensive because it may:
- render at high resolution (many pixels),
- capture many points/legends/text elements,
- produce very large SVGs for huge datasets.

Practical guidance:
- treat “publication export” as a separate step from exploration,
- expect export time to scale with output size and scene complexity,
- if an export mode warns about huge outputs, follow the recommendation (often “optimized/hybrid”).

Related docs:
- {doc}`../k_figure_export/index`

---

## Symptom → likely cause → best knob (cheat sheet)

Use this when you just need the “one thing that helps”.

| Symptom | Most likely bottleneck | Best first knob (safe) | Why it works |
|---|---|---|---|
| Camera movement is choppy | GPU | Reduce window size | fewer pixels shaded per frame |
| Choppy only in a two- or three-view row | GPU | Go back to one view, or add a fourth | the 2×2 grid halves point diameter; a row does not |
| Everything is slower on retina | GPU (pixels) | Reduce window size | fewer pixels shaded per frame |
| Filter sliders stutter | CPU | Live filtering off → `FILTER` once | avoids repeated per-cell recompute |
| Lag increases with more filters | CPU | Disable/remove no-op filters | fewer tests per cell |
| Switching genes is slow | I/O | Use server mode / keep data local | reduces bytes/latency |
| Overlay is slow | GPU | Lower overlay density / disable bloom | reduces per-frame GPU load |
| App crashes / “context lost” | GPU memory | Reduce smoke/overlays/snapshots | reduces VRAM pressure |

---

## “Safe workflows” that keep you fast (and honest)

These habits reduce lag without compromising scientific meaning.

### 1) Work in one view, then snapshot

When exploring:
- keep one live view while tuning filters/fields,
- create snapshots only when you want a stable comparison.

Why this works:
- one view is unambiguously the cheapest state to iterate in, and it keeps the
  comparison you eventually build deliberate rather than accidental. Cellucid
  caps you at three snapshots beside the live view, so a snapshot is a scarce
  slot, not a scratchpad.

### 2) Stage expensive actions

If an action is known to be expensive (filter recompute, big analysis, export):
- make the *decision* first (what range? what groups? what output size?),
- then apply once.

Examples:
- continuous filtering: Live filtering off → move sliders → `FILTER`.
- analysis: define groups → run once → export results if needed.

### 3) Reduce pixels intentionally

If you are GPU-bound:
- shrink the browser window while exploring,
- then return to full-size when you need a “final view”.

This sounds silly, but it is one of the fastest no-risk fixes for GPU-bound workloads.

### 4) Separate “exploration settings” from “presentation settings”

For performance and reproducibility:
- keep a conservative preset (fast, stable),
- keep a separate “presentation preset” (higher quality, slower),
- don’t mix them while debugging performance.

This is especially important for smoke mode and vector overlays.

---

## When the fix is “change the data” (not the UI)

Some problems are best solved upstream:

- **Category explosion** (tens of thousands of categories) makes legends and palettes unusable.
- **Pathological filters** (dozens of stacked filters) means you’re doing data cleaning in the UI.
- **Huge feature sets** (very large var/gene tables) can make gene search and loading heavy.

If you hit these repeatedly, treat it as a data preparation issue:
- pre-filter obvious QC failures in Python,
- choose a smaller set of fields to export,
- collapse/rename categories to a human-usable set.

---

## Interface reference

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel-controls.png
:alt: The Performance Benchmark panel before any run, showing the six point-count presets from 100K to 20M, the custom point count field, the data pattern dropdown, and the Load Synthetic Data, Copy Situation Report and Analyze Performance buttons.
:width: 516px

Every knob in this section that produces a number lives in one place: the
**Performance Benchmark** accordion in the left sidebar. Use it to find the
shape of *your* machine’s limits — the answer is different on every GPU, so
there is no number here to copy. {doc}`05_benchmark_tools` explains each
control and what its readout means.
```

---

## Next steps

- {doc}`03_large_dataset_best_practices` (a “do this in order” playbook for very large datasets)
- {doc}`06_edge_cases_performance` (performance cliffs and surprising slowdowns)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
