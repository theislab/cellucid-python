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

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: fields-coloring-legends-panel-overview
Suggested filename: web_app/fields_legends/01_coloring-filtering-panel-overview.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → index.md
Capture:
  - UI location: left sidebar → “Coloring & Filtering”
  - State prerequisites: dataset loaded; pick a categorical obs field (so the categorical legend is visible)
  - Include: the three field selectors (Categorical obs / Continuous obs / Gene Expression) and the Display options box with the legend
  - Optional: also show the “Active filters” block underneath to hint that legend interactions create real filters
Crop:
  - Include: left sidebar from “Coloring & Filtering” down through “Active filters”
  - Include: enough canvas to show the colored point cloud
  - Exclude: any private dataset identifiers, file paths, browser profile UI
Annotations:
  - Call out: #1 Categorical obs selector, #2 Continuous obs selector, #3 Gene Expression search, #4 Display options (legend), #5 Active filters summary
Alt text:
  - The Coloring & Filtering panel with field selectors and a legend controlling how points are colored and filtered.
Caption:
  - Explain that the field selector chooses what “color-by” means, and the legend is where you interpret and adjust it.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Coloring & Filtering panel.
:width: 100%

The Coloring & Filtering panel is where you pick the active field (obs or gene expression) and use the legend to interpret and adjust coloring and filters.
```

---

## Recommended reading order

1) `d_fields_coloring_legends/01_field_types_and_sources` (what “obs” and “var” mean)
2) `d_fields_coloring_legends/02_field_selector_ux` (how selection/rename/delete/restore actually works)
3) `d_fields_coloring_legends/03_color_by_behavior` (how colors and filters are computed)
4) `d_fields_coloring_legends/04_legend_behavior` (what you can do inside the legend)
5) `d_fields_coloring_legends/05_troubleshooting_fields_legends` (when colors/legends feel “wrong”)
6) `d_fields_coloring_legends/06_screenshots` (checklist for screenshots you’ll add later)

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

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Screenshots
:link: 06_screenshots
:link-type: doc

Screenshot checklist with exact capture instructions and placeholder figures to replace.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
