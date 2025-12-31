# Analysis (Modes, Outputs, Interpretation)

Cellucid’s Analysis panel helps you **turn selections/highlights into interpretable results**:
summary plots, comparisons between groups, and exportable tables/figures.

These pages are written for mixed audiences:
- **Wet lab / non-technical**: click-by-click recipes and “what success looks like”.
- **Computational users**: definitions, assumptions, and interpretation pitfalls.
- **Power users**: caching/state semantics, edge cases, and debugging playbooks.

:::{important}
Analysis results depend on **app state**.

At minimum, always confirm:
- which **view/snapshot** is active,
- what is currently **visible** (filters),
- what groups you’re comparing (**highlight pages** / derived pages),
- and whether the required data (e.g., **gene expression**) is available in your dataset.
:::

:::{tip}
If you’re here because “analysis is empty / missing results”, start with {doc}`10_troubleshooting_analysis`.
:::

---

## Fast path (pick what you’re trying to answer)

| You want to… | Start with | When it’s the right choice | Next page |
|---|---|---|---|
| Get a quick read on a selection/group | Quick | You want a fast “is this different?” summary | `03_analysis_mode_quick_insights` |
| Inspect distributions and details | Detailed Analysis | You need shape/variance/outliers, not just a single number | `04_analysis_mode_detailed_analysis` |
| See relationships between genes/fields | Correlation Analysis | You’re asking “do these move together?” | `05_analysis_mode_correlation_analysis` |
| Compare genes between two groups | Differential Expression (DE) | You have gene expression and two pages to compare (A vs B) | `06_analysis_mode_differential_expression_de` |
| Score a curated gene program | Gene Signature | You have a predefined gene set (pathway/module) | `07_analysis_mode_gene_signature` |
| Find markers for many groups at once | Marker Genes (Genes Panel) | You want one-vs-rest markers for every group in a categorical field | `08_analysis_mode_genes_panel` |

---

## Related concepts (worth skimming first)

- Defining groups: {doc}`../f_highlighting_selection/index`
- Restricting what’s visible: {doc}`../e_filtering/index`
- Loading gene expression (and why it can be missing): {doc}`../b_data_loading/index`
- Exporting plots and publication figures: {doc}`../k_figure_export/index`
- Saving/restoring analysis state (sessions): {doc}`../l_sessions_sharing/index`

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-panel-overview
Suggested filename: analysis/01_analysis-panel-overview.png
Where it appears: User Guide → Web App → Analysis → index.md
Capture:
  - UI location: Analysis panel open (show the mode toggle and at least one result area)
  - State prerequisites: dataset loaded; at least one highlight group or selection available
  - Action to reach state: select a group of cells, open Analysis, run any mode once
Crop:
  - Include: the Analysis panel + enough canvas to show the selected/colored points
  - Exclude: private dataset names, file paths, browser bookmarks, personal avatars
Redact:
  - Remove: any identifying dataset IDs if private
Alt text:
  - Analysis panel open in the left sidebar with mode selector visible.
Caption:
  - Analysis starts from the current view and your defined groups; confirm you are in the expected state before interpreting results.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Analysis panel overview.
:width: 100%

The Analysis panel summarizes and compares groups defined by selection/highlights, and depends on the active view and current filtering state.
```

---

## Recommended reading order

1) `01_analysis_mental_model` (what analysis operates on + what is cached)
2) `02_analysis_ui_overview` (where things live; copy/restore behaviors)
3) Mode pages (`03`–`08`) for the specific analysis you need
4) `09_exporting_analysis_results` (tables/plots, what is reproducible)
5) `10_troubleshooting_analysis` (symptom → diagnosis → fix)
6) `11_screenshots` (capture checklist)

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
