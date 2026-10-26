# Verified performance captures

This page holds the pictures and the measured evidence behind the rest of the
section. Everything here is a record of **one machine on one day**. None of it
is a target, a guarantee, or something to compare your own machine against
directly.

---

## The benchmark panel

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel-controls.png
:alt: The Benchmark panel before any run: Synthetic data with six presets from 500K to 40M, Point count, Pattern and Embeddings, the Load synthetic data button, and Diagnostics with Analyze performance and Copy situation report; the live performance readout is omitted.
:width: 532px

The **Benchmark** accordion as it appears before a run: six presets
from 500K to 40M, a custom point count, the data pattern selector, and the three
actions. The live metric grid that opening the panel shows below **Load
synthetic data** is omitted from the capture; its fields are documented in
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
:alt: The whole window after the Benchmark panel generated one million points of the Model surface pattern: the synthetic point cloud on the canvas and the panel's presets, point count, pattern and embeddings controls in the sidebar; the live readout is omitted.
:width: 1440px

One million points of the **Model surface** pattern, made by **Load
synthetic data**, beside the panel that made them. The live readout is
omitted: it would report the machine that took the picture, not Cellucid.
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
:alt: The whole window in Grid compare with the sidebar open: two equal panels to the right of the sidebar show the Pancreas 3D embedding from the same linked camera — view 1 coloured by cell_type, view 2 recoloured by S_score on a Viridis scale — each under a title chip naming its embedding and field. The sidebar shows Compare views with the two view rows, each marked linked and Orbit, View 2 shows UMAP · 3D, Keep view, Link cameras, Layout Grid compare and Clear, and Coloring & filtering with S_score in the Continuous obs select.
:width: 1440px

Two views share the full canvas height and split the width to the right of the
sidebar equally, so each point keeps its full size and the pair draws every
cell twice — roughly twice the point shading of a single view. This is the most
expensive layout per unit of insight — tune in one view first, then split.
```

```{figure} ../../../_static/screenshots/web_app/dark-theme-multiview.png
:alt: The whole window in the dark theme on the Grid (dark) background, in Grid compare with the sidebar open: two equal panels show the Pancreas 3D embedding coloured by cell_type from the same linked camera, each under a title chip. The sidebar shows Visualization with Theme Dark and Background Grid (dark), and Compare views below it with the two view rows.
:width: 1440px

Theme and background apply to every panel of the comparison layout alike, and
neither choice changes the cost of the frame.
```
