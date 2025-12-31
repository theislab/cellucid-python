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
- number of points drawn (and whether LOD/frustum culling reduces it),
- shader quality (lighting/fog/post-processing),
- point size and screen coverage (many large points = more pixel shading),
- number of views (each view renders its own frame),
- window size and device pixel ratio (retina screens are expensive).

High-leverage knobs:
- **Reduce views** (clear snapshots) — *often the biggest immediate win*.
- **Enable Level-of-Detail (LOD)** (points mode) — reduces draw cost when zoomed out.
- **Lower shader quality / simplify rendering** (points mode).
- **Make the window smaller** while exploring (fewer pixels).

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`
- {doc}`../a_orientation/02_system_requirements` (WebGL context lost = GPU memory pressure)

---

### Volumetric “smoke” mode (cinematic density rendering)

Smoke mode is intentionally heavy and can dominate GPU time.

**The two biggest cost drivers** are documented in {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`:
- **Grid density** (3D volume resolution)
- **Ray quality** (raymarch steps)

Practical consequences:
- Smoke can look “blank” if density is low or if few points are visible (filters/outliers).
- Smoke can cause “context lost” if grid/quality are too high for your GPU.

---

### Filtering (visibility recomputation)

Filtering is a classic CPU hot path.

High-level scaling:
- visibility recomputation is roughly **`O(n_cells × n_enabled_filters)`**.

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
| Camera movement is choppy | GPU | Clear snapshots / reduce views | fewer renders per frame |
| Everything is slower on retina | GPU (pixels) | Reduce window size | fewer pixels shaded per frame |
| Filter sliders stutter | CPU | Live filtering off → `FILTER` once | avoids repeated recompute |
| Lag increases with more filters | CPU | Disable/remove no-op filters | less work per cell |
| Switching genes is slow | I/O | Use server mode / keep data local | reduces bytes/latency |
| Overlay is slow | GPU | Lower overlay density / disable bloom | reduces per-frame GPU load |
| App crashes / “context lost” | GPU memory | Reduce smoke/overlays/views | reduces VRAM pressure |

---

## “Safe workflows” that keep you fast (and honest)

These habits reduce lag without compromising scientific meaning.

### 1) Work in one view, then snapshot

When exploring:
- keep one live view while tuning filters/fields,
- create snapshots only when you want a stable comparison.

Why this works:
- you avoid multiplying GPU work and state updates by `n_views`.

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

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: performance-knobs-overview
Suggested filename: benchmarking_performance/02_performance-knobs-overview.png
Where it appears: User Guide → Web App → Benchmarking and Performance → 02_performance_considerations_what_gets_slow_and_why.md
Capture:
  - UI location: show the most relevant “performance knobs” in a single screenshot (or two)
  - Include: Visualization → renderer settings (LOD), render mode, and an example of Live filtering toggle + FILTER button
  - Include: at least one visible snapshot indicator (so view count is tangible)
Crop:
  - Include: enough left sidebar to see the knobs + enough canvas to see multiple views
Redact:
  - Remove: private dataset identifiers
Alt text:
  - Left sidebar showing performance-related knobs such as render mode, LOD, and Live filtering controls.
Caption:
  - Most performance issues can be fixed by reducing views, reducing pixels, lowering GPU-heavy modes, or avoiding repeated filter recomputation.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for performance-related knobs overview.
:width: 100%

Most performance issues can be fixed by reducing views, reducing pixels, lowering GPU-heavy modes, or avoiding repeated filter recomputation.
```

---

## Next steps

- {doc}`03_large_dataset_best_practices` (a “do this in order” playbook for very large datasets)
- {doc}`06_edge_cases_performance` (performance cliffs and surprising slowdowns)
- {doc}`07_troubleshooting_performance` (symptom → diagnosis → fix)
