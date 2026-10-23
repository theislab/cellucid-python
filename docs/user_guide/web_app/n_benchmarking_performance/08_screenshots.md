# Verified performance captures

This page holds the pictures and the measured evidence behind the rest of the
section. Everything here is a record of **one machine on one day**. None of it
is a target, a guarantee, or something to compare your own machine against
directly.

---

## The benchmark panel

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel-controls.png
:alt: The Performance Benchmark panel before any run, showing the six point-count presets from 500K to 40M, the custom point count field, the data pattern dropdown, and the Load Synthetic Data, Copy Situation Report and Analyze Performance buttons; the live performance readout is omitted.
:width: 516px

The **Benchmark** accordion as it appears before a run: six presets
from 500K to 40M, a custom point count, the data pattern selector, and the three
actions. The live metric grid that opening the panel shows below **Load
Synthetic Data** is omitted from the capture; its fields are documented in
{doc}`05_benchmark_tools`.
```

:::{note}
There is deliberately no screenshot of the metric grid *after* a run.

Those tiles show a live frame rate, so any capture of them records how busy the
capture machine happened to be, not what Cellucid does. A trial capture on a
loaded build machine read “FPS 2” at 100,000 points — a true statement about
that machine and a badly misleading one about the software. The metric names,
units and meanings are written out instead in {doc}`05_benchmark_tools`, where
they can be checked against the interface and do not decay.
:::

---

## Acceptance samples (historical, one machine)

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel.png
:alt: Cellucid Performance Benchmark rendering a 10-million-point GLB model surface while showing the live synthetic-data controls and renderer metrics.
:width: 1440px

Build 2026-07-26.4 rendering the 10-million-point Model surface pattern in
Chromium on an Apple M1 Pro. Read this as evidence that a ten-million-point
scene renders and reports coherent statistics, not as a frame rate to expect.
```

The scripted 1440 × 1000 acceptance series recorded the following samples on
that machine and that build:

| Browser engine | Points | Generation | FPS sample | Frame-time sample |
|---|---:|---:|---:|---:|
| Chromium | 10,000,000 | 1,592 ms | 18 | 56.33 ms |
| Firefox | 10,000,000 | 1,838 ms | 18 | 55.76 ms |
| WebKit | 10,000,000 | 1,735 ms | 19 | 53.74 ms |
| Chromium | 20,000,000 | 3,082 ms | 8 | 123.00 ms |

Every row loaded its pattern exactly once, published the exact requested point
count, returned WebGL error `0`, and emitted no browser error or warning.

:::{warning}
**These figures predate the renderer changes of 2026-08-01 and have not been
re-measured since.**

Two of those changes move the numbers in opposite directions depending on the
scene. Hidden points are now rejected before rasterisation, and synthetic
generation moved off the main thread. A synthetic acceptance run draws every
point it generates, so the first change should barely touch these particular
rows; the generation column is the one most likely to have moved.

Treat the table as a historical acceptance record. If you need a current number,
run the benchmark yourself — {doc}`04_benchmarking_methodology_and_metrics`
explains how to make the result comparable.
:::

Compare measurements only on the same machine, browser build, viewport, and
dataset. Numbers from different engines are different benchmarks.

---

## Layout shapes

Cellucid renders any number of views: the live view plus any number of
snapshots. Two and three views are laid out as a single row; four or more as a
grid of ⌈√n⌉ columns (four as 2×2). That difference matters more than the count, because point size follows
the height of the pane a point is drawn in.

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels fill the area to the right of the sidebar and show the Pancreas 3D embedding from one linked camera — the left panel, the live view, coloured by cell_type, and the right panel, the kept view, recoloured by S_score on a Viridis scale — each under a title chip naming its field. The sidebar, scrolled to Compare Views in an earlier interface, lists badge 1 cell_type and the selected badge 2 S_score above the Layout select set to Grid compare, and Coloring & Filtering shows S_score in the Continuous obs select.
:width: 1440px

Two views share the full canvas height and split the width to the right of the
sidebar equally, so each point keeps its full size and the pair draws every
cell twice — roughly twice the point shading of a single view. This is the most
expensive layout per unit of insight — tune in one view first, then split.
(Earlier interface: the view rows now also name the embedding each view shows.)
```

```{figure} ../../../_static/screenshots/web_app/dark-theme-multiview.png
:alt: The whole window in the dark theme on the Grid (dark) background, in Grid compare with the sidebar open: two equal panels to the right of the sidebar show the Pancreas 3D embedding coloured by cell_type from one linked camera, each under a cell_type title chip. The sidebar, scrolled to Visualization in an earlier interface, reads Theme Dark and Background Grid (dark), and Compare Views below it lists the two view badges.
:width: 1440px

Theme and background apply to every panel of the comparison layout alike, and
neither choice changes the cost of the frame (earlier interface).
```
