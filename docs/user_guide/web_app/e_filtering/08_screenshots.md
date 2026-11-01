# Verified filtering captures

Every image below was captured from the running application on the
**Pancreatic endocrinogenesis (scVelo)** sample (3,696 cells), at a 1440×1000
viewport and twice the device pixel ratio, and is regenerated from a named
scenario in `docs/_tooling/screenshots/scenarios.mjs`. The scenario id is given
with each figure so you can reproduce the state yourself in about a minute.

---

## The panel, with nothing filtered

```{figure} ../../../_static/screenshots/filtering/panel-clusters-unfiltered.png
:alt: The Coloring & filtering panel with the clusters field selected: Color by with its three field selectors, Show cell info, the Legend group with Show centroids, Show labels and the Outlier filter at 100 percent, the Show all and Hide all buttons above the eight-category legend, and Filters reading Showing all 3,696 cells above No filters active.
:width: 532px

The whole control surface filtering is built from, in its resting state.
Everything above the **Legend** heading chooses *what colours the points*;
the **Legend** group, the category rows and **Filters** below it control
*which points are drawn*.
Scenario `filtering-panel-unfiltered`.
```

---

## The four filter surfaces

```{figure} ../../../_static/screenshots/filtering/type-category-visibility.png
:alt: The categorical legend for clusters: Show all and Hide all above eight rows, three of them (Ductal, Ngn3 high EP, Ngn3 low EP) unchecked and reading 0 of their cell counts, each row carrying a colour swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

**Category visibility.** Unchecking a row hides its cells and greys the row.
Scenario `filtering-type-category-visibility`.
```

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: The Visible range block of the continuous S_score legend with Live filtering switched off: Minimum and Maximum sliders with their values, and the Filter button beside Reset enabled, with the pointer on Filter.
:width: 488px

**Numeric range.** With `Live filtering` off, `Filter` appears beside `Reset`
and applies the range once; with it on, only `Reset` is there. Scenario
`filtering-type-numeric-range`.
```


```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: The Outlier filter with its (i), its slider moved left of the end stop and its readout showing 95 percent, with the pointer on the slider.
:width: 496px

**Latent-space outliers.** Scenario `filtering-outlier-slider-95`.
```

---

## What the filter list looks like as you build a stack

```{figure} ../../../_static/screenshots/filtering/active-filters-empty.png
:alt: The Filters block with its (i), reading Showing all 3,696 cells above the words No filters active.
:width: 496px

Nothing filtered. Scenario `filtering-active-filters-empty`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-one-row.png
:alt: The Filters block reading Showing 1,876 of 3,696 cells above one row: a ticked checkbox, the text clusters: hiding Ductal, cut off with an ellipsis, and a remove cross at its right.
:width: 496px

One category filter. Scenario `filtering-active-filters-one-row`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The Filters block holding three rows of different kinds at once: an S_score range row ending in its cell count, cut off with an ellipsis, a clusters hiding row, and a clusters: outlier ≤ 95% row with no count.
:width: 496px

Three filters of three different kinds, stacked. Only the first two carry a
`visible / available` cell count. Scenario
`filtering-active-filters-three-rows`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-row-disabled.png
:alt: The same three-row Filters list with the clusters hiding row unchecked under the pointer: that row is greyed and struck through while the other two stay black, and the count line above has risen.
:width: 496px

The middle filter disabled rather than removed. Scenario
`filtering-active-filters-row-disabled`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-remove-button.png
:alt: A one-row Filters list with the pointer resting on the remove cross at the right edge of the row.
:width: 496px

The remove control. Scenario `filtering-active-filters-remove`.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-scope-tooltip.png
:alt: The (i) of the Filters block pressed, and the information box it opens: every filter hiding cells in the selected view, from the legend, the range sliders and the outlier filter, that unticking one pauses it and removing it drops it, that filters belong to the view they were set in, and Learn more links.
:width: 712px

The (i) of **Filters**: what the list holds, and that filters belong to the
view they were set in. Scenario `filtering-scope-tooltip`.
```

---

## What the filters do to the picture

```{figure} ../../../_static/screenshots/filtering/window-clusters-all-visible.png
:alt: The whole application window with the Pancreas embedding coloured by clusters, every legend row checked, and the Filters block reading Showing all 3,696 cells above No filters active.
:width: 1440px

Before. Scenario `filtering-window-all-visible`.
```

```{figure} ../../../_static/screenshots/filtering/window-clusters-three-hidden.png
:alt: The same window after unchecking Ductal, Ngn3 high EP and Ngn3 low EP: those points are absent from the canvas, their legend rows are unchecked and read 0 of their cell counts, and the Filters block reads Showing 1,876 of 3,696 cells above one clusters row.
:width: 1440px

After. The hidden cells are **not drawn faintly** — they are not drawn at all.
Scenario `filtering-window-categories-hidden`.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-100-off.png
:alt: The whole window coloured by clusters with the Outlier filter at its default 100 percent: scattered cells sit between and around the clusters, and the Filters block reads Showing all 3,696 cells.
:width: 1440px

Outlier filter off. Scenario `filtering-window-outlier-100`.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-95-applied.png
:alt: The same window with the Outlier filter at 95 percent: fewer scattered cells between the clusters, every legend count reduced, and a clusters: outlier ≤ 95% row in the Filters block above Showing 3,506 of 3,696 cells.
:width: 1440px

Outlier filter at 95%. Scenario `filtering-window-outlier-95`.
```

```{figure} ../../../_static/screenshots/filtering/window-continuous-range-applied.png
:alt: The whole window coloured by S_score with the Minimum slider raised: only the higher-scoring cells are drawn, and the Filters block reads a reduced count above one S_score range row.
:width: 1440px

A numeric range on a continuous field. Scenario
`filtering-window-continuous-range`.
```

```{figure} ../../../_static/screenshots/filtering/window-zero-visible-cells.png
:alt: The whole window with a completely empty canvas after pressing Hide all in the clusters legend: every row unchecked and reading 0 of its count, and the Filters block reading Showing 0 of 3,696 cells above one clusters hiding row.
:width: 1440px

Everything filtered out. This is a filter state, not a crash. Scenario
`filtering-window-zero-visible`.
```

---

## States that surprise people

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under a narrow S_score range filter: Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay ticked with small counts, while Beta, Delta, Epsilon and Pre-endocrine are greyed out with disabled checkboxes and read 0 cells.
:width: 488px

Rows with nothing left to show are disabled, not merely empty. Scenario
`filtering-legend-category-unavailable`.
```

```{figure} ../../../_static/screenshots/filtering/window-gene-range-active.png
:alt: The whole window coloured by the gene Ins1 with the Minimum slider raised, so only the cells expressing it are drawn and an Ins1 range row is listed in the Filters block.
:width: 1440px

A gene range filter, listed under Filters; it keeps applying when the view is
coloured by another field. Scenario `filtering-window-gene-range-active`.
```
