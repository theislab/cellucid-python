# Verified coloring and legend captures

All captures on this page come from the **Pancreatic endocrinogenesis (scVelo)**
sample (3,696 cells, 3,753 genes) at a 1440×1000 viewport and twice the device
pixel ratio. Each names the scenario in
`docs/_tooling/screenshots/scenarios.mjs` that regenerates it.

---

## Choosing a field

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: The Color by heading with its (i), above the three field selectors stacked in order — Categorical obs holding clusters, Continuous obs holding None, and an empty Gene expression search box — the Categorical obs row followed by six small icon buttons, the Continuous obs row by five and the gene row by four, with the pointer pressing the first selector.
:width: 504px

The three selectors and their per-field actions. Scenario
`fields-selector-rows`.
```

```{figure} ../../../_static/screenshots/web_app/gene-search-dropdown-open.png
:alt: The Gene expression search box holding the letters Rp with its list open below, naming the genes that contain them.
:width: 496px

Gene search matches a case-insensitive substring anywhere in the name.
Scenario `fields-gene-search-open`.
```

```{figure} ../../../_static/screenshots/web_app/gene-search-no-match.png
:alt: The Gene expression list after typing a symbol this sample does not publish: No gene matches "CD19", a sentence saying it searched the 3,753 genes of this dataset named by its var index and that another var column (symbols) can name them instead, a Choose gene names… button, and a link Why a gene may be missing.
:width: 496px

The empty state, which names this dataset's own published gene count. Scenario
`fields-gene-search-no-match`.
```

---

## The two legends

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: The categorical legend for clusters: the Show all and Hide all buttons above eight rows, each holding a checkbox, a color swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

Categorical. Scenario `fields-legend-categorical`.
```

```{figure} ../../../_static/screenshots/web_app/legend-continuous-sscore.png
:alt: The continuous legend for S_score: the Color scale (Viridis) bar with its two end values, Log color scale off and Rescale to range on, each with its (i), and below a rule the Visible range block with Live filtering on, Minimum and Maximum sliders, and Reset (Filter appears only with Live filtering off).
:width: 488px

Continuous. Scenario `fields-legend-continuous`.
```

```{figure} ../../../_static/screenshots/web_app/legend-colormap-menu-open.png
:alt: The colormap list open under the S_score color bar: seven named gradients, Viridis, Plasma, Inferno, Magma, Cividis, Turbo and Twilight, with Viridis marked as the active one.
:width: 488px

The colormap menu, open. Scenario `fields-colormap-menu`.
```

---

## What that does to the canvas

```{figure} ../../../_static/screenshots/web_app/window-orientation-map.png
:alt: The whole application window with the Pancreas sample loaded, the sidebar showing Coloring & filtering with the embedding colored by clusters, and the eight-row legend beside the matching colors on the canvas.
:width: 1440px

Coloring by a categorical field. The legend swatches and the point colors are
the same palette — this is the correspondence the rest of this section relies
on. Scenario `orientation-window-map`.
```

```{figure} ../../../_static/screenshots/web_app/window-color-by-gene-ins1.png
:alt: The whole window colored by expression of the gene Ins1 on a Viridis scale from 0 to 1260, almost every cell dark and a small patch at the foot of the lower-right arm bright, with the continuous legend in the sidebar.
:width: 1440px

Coloring by a gene. The same continuous legend as a continuous obs field, and
the same controls plus one switch of its own, **Gray for zero**. Scenario `fields-window-gene-coloring`.
```

```{figure} ../../../_static/screenshots/web_app/window-centroids-and-labels.png
:alt: The Pancreas embedding colored by clusters with Show centroids and Show labels both ticked, so each cluster carries a marker at its center with its category name printed on it.
:width: 1440px

Centroids and labels on. Scenario `fields-window-centroids`.
```

---

## The wide overview capture shared with other sections

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample colored by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colors.
:width: 1440px

The overview used across several sections of this documentation: the
Pancreas sample colored by `cell_type`, with **Coloring & filtering** at the
top of the sidebar. Scenario `orientation-app-overview-cell-type`.
```
