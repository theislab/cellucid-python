# Benchmark tools (if exposed)

**Audience:** developers + performance-minded users  
**Time:** 10–25 minutes  
**What you’ll learn:**
- How to run the built-in **synthetic rendering benchmark** (when available)
- How to interpret the metrics (FPS, frame time, GPU memory, LOD)
- How to use “Copy Situation Report” for bug reports and regressions
- What this benchmark does *not* measure (and how to benchmark those separately)

**Prerequisites:**
- Cellucid open in a desktop browser (Chrome/Edge recommended)

---

## Where the benchmark lives (UI location)

If your Cellucid build exposes it, you will find:

- Left sidebar → **Performance Benchmark** accordion

It includes:
- quick presets (100K → 20M points),
- a custom point count input,
- synthetic “data patterns” (clusters, batch effects, atlas-like, etc.),
- a **Load Synthetic Data** button,
- live stats (FPS, frame time, GPU memory, LOD),
- a **Copy Situation Report** button,
- and an **Analyze Performance** tool that suggests what to change.

:::{important}
The synthetic benchmark **loads synthetic data into the viewer**.

Use it when:
- you want to test rendering limits without needing a real dataset, or
- you want to compare hardware/browser builds.

Do **not** use it to benchmark:
- real dataset loading,
- gene expression fetching,
- analysis runtime,
- or figure export.

For those, use {doc}`04_benchmarking_methodology_and_metrics`.
:::

---

## What this benchmark measures (and why it’s useful)

This benchmark is primarily a **rendering stress test**:

- how many points your GPU can render smoothly,
- how well LOD/frustum culling helps,
- and how GPU memory scales with point count.

It is especially useful for answering questions like:
- “Can this laptop handle a 1M-cell dataset interactively?”
- “Did a renderer change reduce FPS on the same machine?”
- “Is the bottleneck likely GPU vs CPU?”

---

## Quick start (3 minutes)

1) Open **Performance Benchmark** in the left sidebar.
2) Click a preset like **500K** or **1M**.
3) Click **Load Synthetic Data**.
4) Orbit/pan/zoom for ~10 seconds and watch:
   - **FPS**
   - **Frame Time**
   - **GPU Memory**
5) If it feels slow, click **Analyze Performance** and read the “What to do” list.

If you want a copy/paste summary for a bug report, click **Copy Situation Report**.

---

## Recommended workflow (reproducible)

If you’re comparing two builds or changes, make the benchmark comparable:

1) Fix window size (pixels matter).
2) Fix view count (single view only; clear snapshots).
3) Start in points mode (smoke mode is a different workload).
4) Use the same preset count and the same pattern.
5) Let FPS stabilize for ~5–10 seconds, then record the metric.

If you want statistics, repeat the run 3–5 times and report the median.

---

## Data patterns (how to choose)

The pattern dropdown controls the spatial distribution of synthetic points.

Practical guidance:
- **Uniform Random**: worst-case for LOD/frustum culling (points everywhere).
- **Gaussian Clusters**: closer to “clustered embedding” intuition.
- **Batch Effects / Atlas-like**: useful for “many groups across space” intuition.
- **Flat UMAP**: approximates a dense 2D manifold style.

If you are testing *rendering limits*, Uniform Random is often the harshest.
If you are testing *realistic single-cell scenes*, Clusters/Atlas-like are usually more representative.

---

## How to interpret the stats

### FPS and frame time

- **60 FPS** is the “smooth” target on most displays.
- **Frame Time** is usually shown in milliseconds:
  - 16.7ms/frame ≈ 60 FPS
  - 33ms/frame ≈ 30 FPS

Interpretation:
- If FPS drops as you increase point count: you’re hitting a GPU throughput limit.
- If FPS is unstable with long spikes (high p95/max frame time): you are seeing “jank” (stutters).

### GPU memory

GPU memory is the main predictor of “WebGL context lost” on large datasets.

Interpretation:
- If GPU memory climbs close to your GPU’s limit, you can see context loss.
- If memory is fine but FPS is low, you’re compute/pixel-bound rather than memory-bound.

### LOD level and visible points

LOD exists to reduce the number of rendered points when zoomed out.

Interpretation:
- If LOD reduces visible points and FPS improves, the renderer is successfully adapting.
- If visible points stays high and FPS stays low, you may need:
  - LOD enabled in Visualization, or
  - fewer views/pixels, or
  - lower shader quality.

Related docs:
- {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` (LOD, culling, shader quality)

---

## “Copy Situation Report” (use this for bug reports)

The **Copy Situation Report** button is meant to produce a compact block of text you can paste into:
- a GitHub issue,
- a performance regression note,
- or a lab internal thread.

Best practice:
- include the report *and* your scenario description (“what I was doing when it was slow”).

If you’re filing an issue, follow {doc}`09_reporting_performance_bugs`.

---

## “Analyze Performance” (why is it slow?)

The **Analyze Performance** button runs a short measurement and produces:
- a verdict (e.g., “GPU-limited”, “pixel-bound”, etc.),
- a list of detected problems,
- and a prioritized “What to do” checklist.

Use it as a triage assistant:
- If it says “pixel-bound”: reduce window size and/or reduce view count.
- If it says “shader overhead”: lower shader quality.
- If it says “LOD overhead”: enable LOD and avoid forcing high LOD levels.

---

## Safety notes and footguns

- Synthetic benchmark data is not your real dataset. Reload the page to restore your original data.
- Don’t start with 20M points on a laptop. Ramp up: 100K → 500K → 1M → 5M.
- If you hit “WebGL context lost”, reduce point count and/or reduce GPU load (fewer views, smaller window).

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: benchmark-panel-overview
Suggested filename: benchmarking_performance/05_benchmark-panel-overview.png
Where it appears: User Guide → Web App → Benchmarking and Performance → 05_benchmark_tools_if_exposed.md
Capture:
  - UI location: Performance Benchmark accordion expanded
  - Show: presets, pattern dropdown, Load Synthetic Data, stats panel, and Analyze Performance button
Crop:
  - Exclude: browser chrome, personal bookmarks, account avatars
Redact:
  - Not needed (synthetic data), but avoid personal info in the browser UI
Alt text:
  - Performance Benchmark panel showing synthetic point presets and performance statistics.
Caption:
  - The built-in benchmark generates synthetic data to stress-test rendering and helps diagnose whether you are GPU-, pixel-, or shader-limited.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Performance Benchmark panel.
:width: 100%

The built-in benchmark generates synthetic data to stress-test rendering and helps diagnose whether you are GPU-, pixel-, or shader-limited.
```

---

## Next steps

- {doc}`04_benchmarking_methodology_and_metrics` (benchmark scenarios beyond synthetic rendering)
- {doc}`06_edge_cases_performance` (performance cliffs and limits)
- {doc}`09_reporting_performance_bugs` (how to file a performance issue with the right context)
