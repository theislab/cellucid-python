# Benchmarking and Performance

These pages help you keep Cellucid **responsive on large datasets**, run **repeatable benchmarks**, and report performance issues in a way that can actually be fixed.

They are written for **mixed audiences**:
- **Wet lab / non-technical**: a clear “is this normal?” checklist and safe defaults.
- **Computational users**: scaling rules, workflow patterns, and the data characteristics that dominate performance.
- **Power users / developers**: how to measure, compare, and debug performance issues without guesswork.

:::{important}
Cellucid performance is usually limited by **one** of these three resources:

- **GPU-bound** (rendering): too many points/views, high visual quality, volumetric smoke, vector field overlays.
- **CPU-bound** (computation): repeated filtering, heavy analysis, large category accounting.
- **Network/storage-bound** (I/O): remote server latency, loading large expression chunks, slow disk/remote mounts.

The fastest fix is almost always: **identify which one you hit first**, then change the *one* knob that actually affects it.
:::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```

---

## Fast path (if it feels slow right now)

Use this table when you’re in the middle of work and just need relief.

| What you see | Most likely bottleneck | Do this first (safe) | Next page |
|---|---|---|---|
| Navigation (orbit/pan/zoom) feels choppy | GPU | Switch to **Points** mode, clear snapshots, close other GPU-heavy tabs | `01_performance_mental_model` |
| Slider scrubbing (filters) stutters | CPU | Turn **Live filtering** off and use **FILTER** to apply once | `02_performance_considerations_what_gets_slow_and_why` |
| Loading takes forever or tab freezes | I/O or memory | Prefer **Server Mode** for large `.h5ad`/`.zarr` and large exports | `03_large_dataset_best_practices` |
| Fans spin up / laptop throttles | GPU + thermals | Lower quality knobs (smoke/grid, overlays), reduce window size, fewer views | `06_edge_cases_performance` |
| “WebGL context lost” / blank canvas | GPU memory | Reload, then reduce GPU load (no smoke, fewer views, lower overlay density) | `07_troubleshooting_performance` |
| “It got slower after a change” | Regression | Run a small benchmark loop and compare numbers | `04_benchmarking_methodology_and_metrics` |

If you suspect an environment/browser issue (not your data), also read:
- {doc}`../a_orientation/02_system_requirements`

---

## Recommended reading order

1) `01_performance_mental_model` (learn the 3-bottleneck model once)
2) `02_performance_considerations_what_gets_slow_and_why` (what actually scales with `n_cells`, `n_views`, etc.)
3) `03_large_dataset_best_practices` (how to stay fast on huge datasets)
4) `07_troubleshooting_performance` (symptom → diagnosis → fix)
5) If you are measuring/optimizing: `04_benchmarking_methodology_and_metrics` → `05_benchmark_tools_if_exposed`
6) If you’re reporting an issue: `09_reporting_performance_bugs`

---

## “Performance” is often a cross-feature issue (related pages)

The biggest performance knobs often live outside this section:

- Filtering hot paths: {doc}`../e_filtering/05_performance_considerations`
- Smoke mode (volumetric) quality/performance knobs: {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke`
- Vector field overlay tuning: {doc}`../i_vector_field_velocity/05_performance_and_quality`
- Large-file loading constraints (browser vs server mode): {doc}`../b_data_loading/01_loading_options_overview`
- Figure export can be expensive for huge SVGs: {doc}`../k_figure_export/04_quality_knobs_and_best_practices`

---

## Pages in this section

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`cpu;1.5em;sd-mr-1` Performance mental model
:link: 01_performance_mental_model
:link-type: doc

How to think about GPU vs CPU vs I/O bottlenecks, and how to triage quickly.
:::

:::{grid-item-card} {octicon}`graph;1.5em;sd-mr-1` What gets slow (and why)
:link: 02_performance_considerations_what_gets_slow_and_why
:link-type: doc

The concrete “cost model”: which interactions scale with `n_cells`, `n_views`, categories, and overlays.
:::

:::{grid-item-card} {octicon}`rocket;1.5em;sd-mr-1` Large dataset best practices
:link: 03_large_dataset_best_practices
:link-type: doc

A practical, step-by-step workflow for staying fast on very large exports.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Benchmarking methodology
:link: 04_benchmarking_methodology_and_metrics
:link-type: doc

How to measure load time, FPS/latency, and memory reliably—and compare changes fairly.
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Benchmark tools (if exposed)
:link: 05_benchmark_tools_if_exposed
:link-type: doc

How to use in-app benchmark tools when available, plus interpretation and common pitfalls.
:::

:::{grid-item-card} {octicon}`alert;1.5em;sd-mr-1` Edge cases
:link: 06_edge_cases_performance
:link-type: doc

Surprising “performance cliffs” like category explosion, too many views, GPU memory limits, and huge exports.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 07_troubleshooting_performance
:link-type: doc

Symptom-based debugging for lag, freezes, context loss, and slow loads.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Screenshots
:link: 08_screenshots
:link-type: doc

Screenshot capture checklist for this section (orientation, success, common failures).
:::

:::{grid-item-card} {octicon}`issue-opened;1.5em;sd-mr-1` Reporting bugs
:link: 09_reporting_performance_bugs
:link-type: doc

A copy/paste template for reporting performance issues with the right context and metrics.
:::

::::
