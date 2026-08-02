# Verified coloring and legend captures

All captures on this page come from the **Pancreatic endocrinogenesis (scVelo)**
sample (3,696 cells, 3,753 genes) at a 1440×1000 viewport and twice the device
pixel ratio. Each names the scenario in
`docs/_tooling/screenshots/scenarios.mjs` that regenerates it.

---

## Choosing a field

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: Three stacked controls labelled CATEGORICAL OBS holding clusters, CONTINUOUS OBS holding None and GENE EXPRESSION showing the placeholder Search genes, each followed by four small icon buttons that are greyed out on the two empty rows, with the mouse pointer pressing the first dropdown.
:width: 488px

The three selectors and their four per-field actions. Scenario
`fields-selector-rows`.
```

```{figure} ../../../_static/screenshots/web_app/gene-search-dropdown-open.png
:alt: The GENE EXPRESSION row with the letters Rp typed into the search box and an open dropdown beneath listing matching mouse gene names one per row.
:width: 480px

Gene search matches a case-insensitive substring anywhere in the name.
Scenario `fields-gene-search-open`.
```

```{figure} ../../../_static/screenshots/web_app/gene-search-no-match.png
:alt: The same row with CD19 typed in and a dropdown reading No gene matches, a sentence giving the number of gene names this dataset publishes and warning it may be a subset of the source data, and a link reading Why a gene may be missing.
:width: 480px

The empty state, which names this dataset's own published gene count. Scenario
`fields-gene-search-no-match`.
```

---

## The two legends

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: A legend with the hint Click a swatch to pick a color, Show All and Hide All buttons, and eight alphabetically ordered rows each with a checkbox, a coloured swatch, a category name, a cell count, a pencil icon and a bin icon.
:width: 428px

Categorical. Scenario `fields-legend-categorical`.
```

```{figure} ../../../_static/screenshots/web_app/legend-continuous-sscore.png
:alt: A legend headed COLOR SCALE (VIRIDIS) with a purple to yellow gradient bar labelled minus 0.37 and 1.14, a Log color scale toggle reading Off, a Rescale colorbar to slider range toggle reading On, and a FILTERING block with a Live filtering toggle reading On, Min and Max sliders and FILTER and RESET buttons.
:width: 428px

Continuous. Scenario `fields-legend-continuous`.
```

```{figure} ../../../_static/screenshots/web_app/legend-colormap-menu-open.png
:alt: The same continuous legend with the gradient bar expanded into a grid of small named gradient swatches, Viridis marked active, followed by Plasma, Inferno, Magma, Cividis, Turbo and many more.
:width: 428px

The colormap menu, open. Scenario `fields-colormap-menu`.
```

---

## What that does to the canvas

```{figure} ../../../_static/screenshots/web_app/window-orientation-map.png
:alt: The whole application window with the sidebar scrolled to Coloring and Filtering, the Pancreas 3D embedding drawn on a light grid, and the eight-row clusters legend whose swatch colours match the coloured groups of points on the canvas.
:width: 1440px

Colouring by a categorical field. The legend swatches and the point colours are
the same palette — this is the correspondence the rest of this section relies
on. Scenario `orientation-window-map`.
```

```{figure} ../../../_static/screenshots/web_app/window-color-by-gene-ins1.png
:alt: The whole application window coloured by expression of the gene Ins1 on a dark-purple to yellow Viridis scale, with most of the embedding dark and one region bright, and the continuous legend open in the sidebar.
:width: 1440px

Colouring by a gene. The same continuous legend as a continuous obs field, and
the same controls. Scenario `fields-window-gene-coloring`.
```

```{figure} ../../../_static/screenshots/web_app/window-centroids-and-labels.png
:alt: The whole application window coloured by clusters with Show centroids and Show labels both ticked, and the category names Alpha, Beta, Delta, Ductal, Epsilon, Ngn3 high EP, Ngn3 low EP and Pre-endocrine printed on the canvas at the centre of each coloured group.
:width: 1440px

Centroids and labels on. Scenario `fields-window-centroids`.
```

---

## The wide overview capture shared with other sections

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The Cellucid web app with the sidebar scrolled to the Session panel showing the Suo dataset summary, and a large multicoloured single-cell embedding filling the canvas on a light grid.
:width: 1440px

The wide **Suo** overview used across several sections of this documentation.
It is a 1× capture of a different sample, and its sidebar is scrolled to
**Session** rather than to the legend, so use the Pancreas captures above when
you need to read a control. It is kept here because six other pages reference
it.
```
