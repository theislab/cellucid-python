# Screenshots (benchmarking and performance)

This page is a **screenshot capture checklist** for the Benchmarking and Performance section.

It exists so you (or a collaborator) can capture screenshots once, systematically, without hunting through every page.

---

## How placeholders work

All pages in this section can use:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Recommended screenshot storage:

- `cellucid-python/docs/_static/screenshots/benchmarking_performance/`

Each screenshot placeholder block should include an HTML comment describing:
- what state to capture,
- what to crop/redact,
- suggested filename conventions,
- and suggested alt text + caption content.

General guidance lives in:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## Recommended screenshot set (benchmarking + performance)

### Overview / mental model (1 screenshot)

1) **Three bottlenecks diagram** (GPU vs CPU vs I/O)

Page:
- `n_benchmarking_performance/01_performance_mental_model`

Suggested filename:
- `benchmarking_performance/01_three-bottlenecks-diagram.png`

---

### Performance knobs (1 screenshot)

2) **Performance knobs overview** (views + render mode + LOD + Live filtering)

Page:
- `n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`

Suggested filename:
- `benchmarking_performance/02_performance-knobs-overview.png`

---

### Views multiplier (1–2 screenshots)

3) **Single view vs grid view comparison**

Page:
- `n_benchmarking_performance/03_large_dataset_best_practices`

Suggested filename:
- `benchmarking_performance/03_grid-vs-single-view.png`

Nice-to-have:
- A second screenshot showing the same view with overlays enabled (to emphasize “pixels × views”).

---

### Benchmarking tools (1 screenshot)

4) **Browser performance monitor or task manager next to Cellucid**

Page:
- `n_benchmarking_performance/04_benchmarking_methodology_and_metrics`

Suggested filename:
- `benchmarking_performance/04_devtools-performance-monitor.png`

---

### In-app benchmark (if exposed) (1 screenshot)

5) **Performance Benchmark accordion overview**

Page:
- `n_benchmarking_performance/05_benchmark_tools_if_exposed`

Suggested filename:
- `benchmarking_performance/05_benchmark-panel-overview.png`

Nice-to-have:
- A screenshot after clicking “Analyze Performance” showing the verdict + suggested fixes.

---

### Failure modes (high value support screenshots)

6) **WebGL context lost overlay**

Page:
- `n_benchmarking_performance/07_troubleshooting_performance` (also referenced in `a_orientation/02_system_requirements`)

Suggested filename:
- `benchmarking_performance/70_context-lost.png`

7) **Live filtering OFF with FILTER enabled**

Page:
- `n_benchmarking_performance/07_troubleshooting_performance` (and also relevant to fields/filtering pages)

Suggested filename:
- `benchmarking_performance/71_live-filtering-off.png`

8) **Category explosion example** (legend unusable / too many categories)

Page:
- `n_benchmarking_performance/06_edge_cases_performance`

Suggested filename:
- `benchmarking_performance/72_category-explosion.png`

---

## Nice-to-have screenshots (optional but very helpful)

9) **Retina vs external monitor / window-size effect**

Goal:
- demonstrate that “smaller window → smoother” is real.

Suggested filename:
- `benchmarking_performance/80_pixels-matter.png`

10) **Smoke mode “too slow” example**

Goal:
- show smoke mode controls with a conservative preset and a warning about grid density.

Suggested filename:
- `benchmarking_performance/81_smoke-performance-knobs.png`

11) **Overlay “bloom strength 0” as a fast fix**

Goal:
- show overlay advanced settings with bloom set to 0.00.

Suggested filename:
- `benchmarking_performance/82_overlay-bloom-off.png`

