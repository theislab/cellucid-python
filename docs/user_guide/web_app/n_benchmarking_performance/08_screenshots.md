# Verified performance captures

## Benchmark panel

```{figure} ../../../_static/screenshots/benchmarking_performance/benchmark-panel.png
:alt: Performance Benchmark panel showing synthetic-data controls and live renderer metrics.
:width: 100%

The benchmark panel reports points, FPS, frame time, GPU memory estimate, LOD
level, visible points, and timing distribution for the current view.
```

The values in this capture describe one browser run and are not universal
performance claims. Compare measurements on the same machine, browser, viewport,
and dataset.

## View-count cost

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: Cellucid showing two side-by-side views of the same dataset.
:width: 100%

Each additional visible view adds rendering work, so tune one view before
creating a small comparison layout.
```

```{figure} ../../../_static/screenshots/web_app/dark-theme-multiview.png
:alt: Two Cellucid views displayed side by side with the dark theme.
:width: 100%

Theme and background changes apply consistently across the active comparison
layout.
```
