# Verified filtering captures

Every image below was captured from the running application on the
**Pancreatic endocrinogenesis (scVelo)** sample (3,696 cells), at a 1440×1000
viewport and twice the device pixel ratio, and is regenerated from a named
scenario in `docs/_tooling/screenshots/scenarios.mjs`. The scenario id is given
with each figure so you can reproduce the state yourself in about a minute.

---

## The panel, with nothing filtered

```{figure} ../../../_static/screenshots/filtering/panel-clusters-unfiltered.png
:alt: The Coloring and Filtering panel showing a Categorical obs select holding clusters, empty Continuous obs and Gene Expression selects, a Display options box with the outlier slider at 100 percent, Show centroids and Show labels checkboxes, Show All and Hide All buttons, eight category rows with colour swatches and cell counts, and an Active filters block reading Showing all 3,696 points above No filters active.
:width: 516px

The whole control surface filtering is built from, in its resting state.
Everything above the dashed **Display options** box chooses *what colours the
points*; everything inside it and below it controls *which points are drawn*.
Scenario `filtering-panel-unfiltered`.
```

---

## The four filter surfaces

```{figure} ../../../_static/screenshots/filtering/type-category-visibility.png
:alt: The categorical legend for clusters with Show All and Hide All buttons above eight rows; Ductal, Ngn3 high EP and Ngn3 low EP are unchecked with grey swatches and read 0 of their total cells, while the five checked rows show a plain cell count.
:width: 428px

**Category visibility.** Unchecking a row hides its cells and greys the row.
Scenario `filtering-type-category-visibility`.
```

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: The Filtering block of the continuous legend, headed Filtering and Visible range, with the help line Adjust limits, then click FILTER or enable Live filtering, a Live filtering toggle reading Off, the hint Drag sliders then click Filter, Min and Max sliders with numeric readouts, and black FILTER and RESET buttons with the mouse pointer on FILTER.
:width: 428px

**Numeric range.** With `Live filtering` off, `FILTER` becomes clickable and
applies the range once. Scenario `filtering-type-numeric-range`.
```

```{figure} ../../../_static/screenshots/filtering/filter-button-disabled-live-on.png
:alt: The same Filtering block with the Live filtering toggle reading On, the hint Changes apply as you drag, and a pale grey FILTER button under the mouse pointer beside a black RESET button.
:width: 428px

**The same block with `Live filtering` on.** `FILTER` is greyed out and cannot
be clicked, because the range is already being applied on every slider step.
Scenario `filtering-button-disabled-live-on`.
```

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: A control labelled Outlier filter (latent space) with its slider handle a short distance in from the right end and the number 95 percent printed beside it, the mouse pointer resting on the handle.
:width: 436px

**Latent-space outliers.** Scenario `filtering-outlier-slider-95`.
```

---

## What the filter list looks like as you build a stack

```{figure} ../../../_static/screenshots/filtering/active-filters-empty.png
:alt: A block headed ACTIVE FILTERS (SELECTED VIEW ONLY) with a small circled i button, the line Showing all 3,696 points, and a bordered area containing the words No filters active.
:width: 480px

Nothing filtered. Scenario `filtering-active-filters-empty`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-one-row.png
:alt: The same block reading Showing 1,876 of 3,696 points above one row with a ticked blue checkbox, the truncated text clusters, hiding Ductal, Ngn3, and a circled cross at the right edge.
:width: 480px

One category filter. Scenario `filtering-active-filters-one-row`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The same block reading Showing 709 of 3,696 points above three ticked rows, the first reading S_score minus 0.15 to 1.14 followed by a cell count, the second reading clusters colon hiding Ductal, Epsilon truncated, the third reading clusters colon outlier less than or equal to 95 percent with no count after it.
:width: 480px

Three filters of three different kinds, stacked. Only the first two carry a
`visible / available` cell count. Scenario
`filtering-active-filters-three-rows`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-row-disabled.png
:alt: The same three rows with the middle checkbox cleared, that row greyed and struck through while the other two stay black, the mouse pointer on the cleared checkbox, and the line above now reading Showing 1,472 of 3,696 points.
:width: 480px

The middle filter disabled rather than removed. Scenario
`filtering-active-filters-row-disabled`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-remove-button.png
:alt: A one-row filter list with the mouse pointer resting on the circled cross at the right edge of the row.
:width: 480px

The remove control. Scenario `filtering-active-filters-remove`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-scope-tooltip.png
:alt: The Active filters heading with its circled i button pressed and highlighted, and a pop-up panel opened to the right reading that filters stay with the selected panel when its coloring changes and other panels remain unchanged.
:width: 908px

Why the heading says *selected view only*. Scenario
`filtering-scope-tooltip`.
```

---

## What the filters do to the picture

```{figure} ../../../_static/screenshots/filtering/window-clusters-all-visible.png
:alt: The whole application window with the sidebar on the left and the Pancreas 3D embedding on a light grid, all eight legend rows ticked with their cell counts, and Showing all 3,696 points below them.
:width: 1440px

Before. Scenario `filtering-window-all-visible`.
```

```{figure} ../../../_static/screenshots/filtering/window-clusters-three-hidden.png
:alt: The same window after unchecking three rows; the orange, purple and blue point groups have gone from the canvas, those three legend rows read 0 of 916, 0 of 642 and 0 of 262 cells, and the sidebar reads Showing 1,876 of 3,696 points above one filter row.
:width: 1440px

After. The hidden cells are **not drawn faintly** — they are not drawn at all.
Scenario `filtering-window-categories-hidden`.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-100-off.png
:alt: The whole window coloured by clusters with the outlier slider at its right end reading 100 percent, loose points scattered around and between the coloured groups, and Showing all 3,696 points.
:width: 1440px

Outlier filter off. Scenario `filtering-window-outlier-100`.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-95-applied.png
:alt: The same window with the outlier slider reading 95 percent and the mouse pointer on it, noticeably fewer loose points around the groups, every legend count reduced, and Showing 3,518 of 3,696 points above a row reading clusters colon outlier less than or equal to 95 percent.
:width: 1440px

Outlier filter at 95%. Scenario `filtering-window-outlier-95`.
```

```{figure} ../../../_static/screenshots/filtering/window-continuous-range-applied.png
:alt: The whole window coloured by S_score on a Viridis scale with the Min slider moved to the right, only a scattered subset of cells drawn, and the sidebar reporting a reduced count above one S_score range row.
:width: 1440px

A numeric range on a continuous field. Scenario
`filtering-window-continuous-range`.
```

```{figure} ../../../_static/screenshots/filtering/window-zero-visible-cells.png
:alt: The whole window with a completely blank grid where the embedding was, every legend row unchecked and reading 0 of its total, and the sidebar reading Showing 0 of 3,696 points above one clusters hiding row.
:width: 1440px

Everything filtered out. This is a filter state, not a crash. Scenario
`filtering-window-zero-visible`.
```

---

## States that surprise people

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under another filter, where Beta, Delta, Epsilon and Pre-endocrine are greyed with unticked disabled checkboxes and read 0 cells, while Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay black with small counts.
:width: 428px

Rows with nothing left to show are disabled, not merely empty. Scenario
`filtering-legend-category-unavailable`.
```

```{figure} ../../../_static/screenshots/filtering/window-gene-range-active.png
:alt: The whole window coloured by the gene Ins1 with the Min slider raised, a small cluster of surviving points on an otherwise empty grid, and the sidebar listing one Ins1 range row.
:width: 1440px

A gene range filter while the gene is the active field. Scenario
`filtering-window-gene-range-active`.
```

```{figure} ../../../_static/screenshots/filtering/window-gene-range-scope-lost.png
:alt: The same window after switching the Categorical obs select to clusters; the full coloured embedding is back, the filter list reads No filters active and the line above reads Showing all 3,696 points.
:width: 1440px

The same session one field-change later. The gene filter stopped applying.
Scenario `filtering-window-gene-range-scope-lost`.
```
