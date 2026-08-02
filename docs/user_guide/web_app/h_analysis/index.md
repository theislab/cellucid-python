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
**New here? Read {doc}`11_screenshots` first.** It walks one analysis from a
mouse click to an exported table, in eight numbered steps with a screenshot at
each one, and assumes nothing.

If you’re here because “analysis is empty / missing results”, start with {doc}`10_troubleshooting_analysis`.
:::

---

## Fast path (pick what you’re trying to answer)

| You want to… | Start with | When it’s the right choice | Next page |
|---|---|---|---|
| Get a quick read on a selection/group | Quick | You want a fast “is this different?” summary | {doc}`03_analysis_mode_quick_insights` |
| Inspect distributions and details | Detailed | You need shape/variance/outliers, not just a single number | {doc}`04_analysis_mode_detailed_analysis` |
| See relationships between genes/fields | Correlation | You’re asking “do these move together?” | {doc}`05_analysis_mode_correlation_analysis` |
| Compare genes between two groups | Differential Expression | You have gene expression and two pages to compare (A vs B) | {doc}`06_analysis_mode_differential_expression_de` |
| Score a curated gene program | Gene Signature | You have a predefined gene set (pathway/module) | {doc}`07_analysis_mode_gene_signature` |
| Find markers for many groups at once | Marker Genes | You want one-vs-rest markers for every group in a categorical field | {doc}`08_analysis_mode_genes_panel` |

The names in the second column are the labels the app itself uses, exactly as
they appear in the panel.

---

## Related concepts (worth skimming first)

- Defining groups: {doc}`../f_highlighting_selection/index`
- Restricting what’s visible: {doc}`../e_filtering/index`
- Loading gene expression (and why it can be missing): {doc}`../b_data_loading/index`
- Exporting plots and publication figures: {doc}`../k_figure_export/index`
- Saving/restoring analysis state (sessions): {doc}`../l_sessions_sharing/index`

---

## Interface reference

```{figure} ../../../_static/screenshots/analysis/analysis-panel-modes.png
:alt: Six stacked rows reading Quick (Automatic insights), Detailed (Full control over options), Correlation (Explore variable relationships), Differential Expression (Find DE genes between groups), Gene Signature (Compute signature scores) and Marker Genes (Discover markers across groups). Each row carries a small copy-window icon and a chevron.
:width: 808px

The Analysis panel exposes six explicit modes, each with its own validated inputs
and result area. The icon on each row undocks that mode into a floating window;
the chevron opens it in place.
```

---

## Recommended reading order

1) {doc}`11_screenshots` — **one analysis from start to finish**, in eight
   screenshotted steps. Start here if you have not used the panel before.
2) {doc}`01_analysis_mental_model` (what analysis operates on + what is cached)
3) {doc}`02_analysis_ui_overview` (where things live; sidebar vs expanded; copy/restore behaviors)
4) Mode pages (`03`–`08`) for the specific analysis you need
5) {doc}`09_exporting_analysis_results` (tables/plots, what is reproducible)
6) {doc}`10_troubleshooting_analysis` (symptom → diagnosis → fix)

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
