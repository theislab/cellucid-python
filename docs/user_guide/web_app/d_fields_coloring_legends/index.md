# Fields, Coloring, and Legends (Obs/Var)

These pages explain how Cellucid turns **per-cell metadata (obs)** and **gene expression (var)** into:
- **Colors** in the viewer (“Color by”)
- **Filters** (hide categories; numeric range filters; latent outlier filtering when available)
- **Legends** (the control surface for interpreting and editing what you see)

They are written for mixed audiences:
- **Wet lab / non-technical**: click-by-click, what to look for, and safe defaults.
- **Computational users**: exact semantics, edge cases (NaN/Inf, log scale), and performance notes.
- **Power users**: how editing categories creates derived fields, how restore/purge works, and how to debug “why did points disappear?”

:::{important}
Cellucid treats “coloring” and “filtering” as coupled:

- Choosing a field determines the **legend UI** you see.
- Legend interactions often create **real filters** (not just visual tweaks).
:::

---

## Fast path (what most users want)

| You want to… | Do this in the UI | Why |
|---|---|---|
| Color by clusters / sample / batch | **Coloring & Filtering → Categorical obs** | Best for discrete labels; legend lets you hide categories and edit labels/colors |
| Color by a score (QC, pseudotime, marker score) | **Coloring & Filtering → Continuous obs** | Gives a continuous colorbar; you can filter by min/max and rescale the palette |
| Color by a gene | **Coloring & Filtering → Gene Expression** | Search genes, load on demand, and use the same continuous legend controls |

---

## The three routes, on screen

All three rows of the table above live in the same place, one under the other:

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: Three stacked controls labelled CATEGORICAL OBS holding clusters, CONTINUOUS OBS holding None, and GENE EXPRESSION holding the placeholder Search genes, each followed by four small icon buttons, with the mouse pointer pressing the first dropdown.
:width: 488px

Only one of the three can be active at a time. Choosing a field in one selector
clears the other two.
```

Which one you choose decides which legend you get, and the legend is where most
of the work happens:

```{figure} ../../../_static/screenshots/web_app/legend-categorical-clusters.png
:alt: A legend with the hint Click a swatch to pick a color, Show All and Hide All buttons, and eight rows each with a checkbox, a coloured swatch, a category name, a cell count, a pencil icon and a bin icon.
:width: 428px

**Categorical** → a list of rows you can tick, recolour, rename and merge.
```

```{figure} ../../../_static/screenshots/web_app/legend-continuous-sscore.png
:alt: A legend headed COLOR SCALE (VIRIDIS) with a purple to yellow gradient bar under it, the numbers minus 0.37 and 1.14 at its ends, a Log color scale toggle reading Off with the note Zero negative and NaN values use the None color, a Rescale colorbar to slider range toggle reading On with the note On colors track sliders Off full data range, then a FILTERING section with a VISIBLE RANGE subtitle, a Live filtering toggle reading On, Min and Max sliders with numeric readouts, and FILTER and RESET buttons.
:width: 428px

**Continuous or gene** → a colour bar, two colour-mapping toggles, and a
filtering block.
```

---

## Recommended reading order

1) {doc}`01_field_types_and_sources` (what “obs” and “var” mean)
2) {doc}`02_field_selector_ux` (how selection/rename/delete/restore actually works)
3) {doc}`03_color_by_behavior` (how colors and filters are computed)
4) {doc}`04_legend_behavior` (what you can do inside the legend)
5) {doc}`05_troubleshooting_fields_legends` (when colors/legends feel “wrong”)
6) {doc}`06_screenshots` (verified interface captures)
7) {doc}`07_genes_in_the_built_in_samples` (why a real gene can be missing from a sample)

---

## Pages in this section

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Field Types & Sources
:link: 01_field_types_and_sources
:link-type: doc

Obs vs var, categorical vs continuous, and what “gene expression” means in Cellucid.
:::

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Field Selector UX
:link: 02_field_selector_ux
:link-type: doc

How to pick a field, search genes, duplicate/rename/delete fields, and restore from Deleted Fields.
:::

:::{grid-item-card} {octicon}`sliders;1.5em;sd-mr-1` Color-by Behavior
:link: 03_color_by_behavior
:link-type: doc

How categorical vs continuous coloring works, how missing values look, and what changes are “real filters”.
:::

:::{grid-item-card} {octicon}`eye;1.5em;sd-mr-1` Legend Behavior
:link: 04_legend_behavior
:link-type: doc

All legend interactions: colormaps, log scale, range sliders, category colors, rename/merge/delete categories.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 05_troubleshooting_fields_legends
:link-type: doc

Symptom → diagnosis → fix for missing fields, slow genes, “everything is gray”, and more.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Verified captures
:link: 06_screenshots
:link-type: doc

Current coloring and categorical-legend states captured from the running app.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Genes in the built-in samples
:link: 07_genes_in_the_built_in_samples
:link-type: doc

Which genes each sample publishes, why the rest are not there, and what to do
when you need one that is missing.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
