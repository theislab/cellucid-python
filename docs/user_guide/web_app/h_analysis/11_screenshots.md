# The analysis lifecycle, end to end

This page is one continuous walkthrough: a cell group is made, a mode is chosen,
it is configured, it runs, the result is read, one setting is changed, and the
answer leaves the browser. Every screenshot below is a frame of that single
session, in order, on the same dataset.

If you have never used the app before, start here and do the steps as you read
them. Nothing on this page assumes you have read the mode pages first.

:::{admonition} What was captured, exactly
:class: note

- **Dataset:** the built-in **Pancreatic endocrinogenesis (scVelo)** sample —
  3,696 cells, 3,753 genes, 8 cell types. Pick it from **Sample dataset** in the
  Session panel.
- **View:** 2D embedding, coloured by `cell_type`, light theme.
- **Window:** 1440 × 1000 CSS pixels, captured at 2× pixel density.
- **Group used throughout:** one annotation group — the 262 `Ngn3 low EP` cells —
  as **Page 1**, with its complement **Rest of Page 1** (3,434 cells).

Your own numbers will differ if you select a different group. The shapes of the
screens will not.
:::

---

## Step 1 — Turn a click into a group

Analysis never operates on "what I am looking at". It operates on **pages**:
named, explicit sets of cells. Until at least one page has cells in it, every
analysis mode is empty, and that is the single most common reason a first-time
user thinks analysis is broken.

Making a page takes two actions, and **both are required**:

1. In **Coloring & filtering**, choose a **Categorical obs** field —
   `cell_type` here. Annotation selection reads the clicked cell's value on the
   active field, so with no active field the tool never starts.
2. In **Highlighting**, leave the mode on **Annotation**, then
   **Alt+click** any cell. Every cell sharing that cell's annotation value is
   picked up at once. The panel now reads `Step 1: 262 cells`.

```{figure} ../../../_static/screenshots/analysis/life-01-confirm-a-selection.png
:alt: The Highlighting panel with the Annotation tool chosen and a pending step reading "Step 1: 262 cells" with Confirm, undo, redo and Cancel, above Page 1, which still holds no cells. The pointer is on Confirm.
:width: 532px

A selection is a **pending step** until you confirm it. Note that `Page 1` still
reads `(0)` and the group list below still says "No cells highlighted" — the 262
cells are staged, not saved.
```

Press **Confirm**. Only now does Page 1 acquire cells.

```{figure} ../../../_static/screenshots/analysis/life-02-two-pages-ready.png
:alt: The same Highlighting panel after Confirm: Page 1 holds 262 cells, and the list below it reads 262 cells highlighted with Clear, above the confirmed group with its checkbox and remove button.
:width: 532px

After **Confirm**, `Page 1` holds 262 cells and the group appears in the list
below with its own enable and remove controls. Analysis can now see it.
```

:::{tip}
You get a second page for free. Every base page has a derived complement called
`Rest of <page>` — here `Rest of Page 1`, the other 3,434 cells. That is enough
to run a two-group comparison without selecting anything else.
:::

---

## Step 2 — Choose a mode

Open the **Analysis** section. It lists its modes in three groups (Genes, Cells,
Explore), all collapsed. Opening one closes the others, so exactly one mode is
active at a time.

```{figure} ../../../_static/screenshots/analysis/analysis-panel-modes.png
:alt: The Analysis panel with every analysis collapsed, in three groups: Genes (Marker genes, Differential expression, Dot plot, Gene signature), Cells (Composition) and Explore (Quick, Detailed, Correlation); each row carries a chevron, an (i) and a copy-window button.
:width: 532px

The modes. Each row carries a **copy** button, which undocks that mode into a
floating window, and a chevron, which opens it in place.
```

| Mode | Runs | Needs | Answers |
| --- | --- | --- | --- |
| **Quick** | automatically | one page with cells | "What is in this group?" |
| **Detailed** | automatically | one variable, one or more pages | "How does this variable differ between my groups?" |
| **Correlation** | on **Run** | two continuous or gene variables | "Do these two move together?" |
| **Differential Expression** | on **Run** | two groups of cells (pages or categories), gene data | "Which genes separate group A from group B?" |
| **Gene Signature** | on **Run** | one or more gene sets | "How strongly does this program score per cell?" |
| **Marker Genes** | on **Run** | a categorical field, gene data | "What marks every group at once?" |
| **Dot plot**, **Composition** | on **Run** | genes and a field / samples and conditions | see {doc}`13_analysis_mode_dot_plot`, {doc}`12_analysis_mode_composition` |

Quick and Detailed recompute as you change a control. Every other mode is a form
that ends in **Run**: essentials first, everything else under a collapsed **More
options** (see {doc}`02_analysis_ui_overview`).

We will follow **Differential Expression**, because it is the mode with the most
moving parts. The others are the same shape.

---

## Step 3 — Configure it

Open **Differential Expression**. The essentials ask for three things:
**Compare** (Pages, here), **Group A** and **Group B**.

```{figure} ../../../_static/screenshots/analysis/de-page-selection.png
:alt: The top of the Differential expression form: Compare set to Pages, Group A reading Page 1 and Group B reading Rest of Page 1.
:width: 464px

**Group A vs Group B.** Two different pages that share no cell. The base page
and its `Rest of` complement are the simplest valid pair. Until the pair is
valid, Run is disabled and the line above it says why.
```

**Which page goes in A?** The sign convention follows A: a gene with positive
log₂ fold change is higher in **Group A**. Put the group you are characterising
in A and the background in B, and "upregulated" then means what you expect;
**Swap A and B** exchanges them.

**More options** can stay collapsed. It holds the cells to compare within, the
test (Wilcoxon rank-sum by default; two t-tests), the p-value convention, the
multiple-testing correction, the expression matrix, the fold change,
**Replicates** (pseudobulk), the filters and the labels — every
default is in the table on {doc}`06_analysis_mode_differential_expression_de`.
There are no performance settings: the run sizes itself.

---

## Step 4 — Let it run

Press **Run**. The button reads `Running…`, and a box under it says what the
run covers, what it is doing, and offers **Cancel**.

```{figure} ../../../_static/screenshots/analysis/de-running.png
:alt: The Differential expression panel while it runs: the Run button disabled, and under it a box naming the comparison (Page 1 against Rest of Page 1, with their cell counts), the step under way with how many genes of how many are done, a progress bar and a Cancel button.
:width: 488px

The first step, `Testing genes`: every gene is read and tested, counted as it goes.
Step 2, `Applying multiple testing correction`, adjusts the false-discovery rate
across all of them at once and is usually over in an instant.
```

**Cancel** stops the run, its reads and its computations, and keeps nothing; the
box then reads `Differential expression cancelled` with **Run again**.
Changing the form while a run is in flight does not stop it: the result arrives
dimmed, marked **Stale · Run**, so a result that no longer matches the form is
never mistaken for one that does.

---

## Step 5 — Read the preview

When the run finishes, a small volcano plot appears below the form, with an
**⤢ Expand** button under it.

```{figure} ../../../_static/screenshots/analysis/de-sidebar-result.png
:alt: The sidebar result of a differential expression run: a small volcano plot with the genes up in Page 1 in red on the right, those down in blue on the left and the rest grey, one of its gene labels shown, and the Expand button below it under the pointer.
:width: 440px

The sidebar preview is a **glance**, not a workspace. At 224 pixels wide the plot
is confirmation that the run produced something; the numbers live in the full
view. There is no legend here because the band under the axis fits the tick
labels and the axis title *or* the legend, and an unlabelled axis is worse than
an unkeyed volcano.
```

:::{important}
The sidebar shows the plot and one line naming what was compared. Summary
statistics, the gene table and every export live **exclusively** in the expanded
view. If you are looking for a number and cannot find it, you have not expanded
yet.
:::

The **⤢ Expand** button is the only way in. The preview itself is not clickable,
and the button is a real button: you can reach it with `Tab` and open the full
view with `Enter` or `Space`.

---

## Step 6 — Open the full view

```{figure} ../../../_static/screenshots/analysis/de-volcano-expanded.png
:alt: A large modal headed Differential expression: Page 1 vs Rest of Page 1: a volcano plot with seven labelled genes top-left over a legend reading Up (395), Down (698) and Not significant (2660), an Export row of PNG, SVG and CSV buttons above the Plot options on the right, a Summary statistics table bottom-left, and a Statistical analysis table of the top genes bottom-right.
:width: 1440px

The expanded view has four regions: the **plot**, the **options** column (which
also carries the export row), **summary statistics**, and **statistical
analysis**. Every mode uses this same four-region layout. The dotted seams between
them are drag handles.
```

Read it in this order.

**The volcano.** Horizontal position is log₂ fold change: right means higher in
Group A, left means higher in Group B. Vertical position is `−log₁₀` of the
adjusted p-value, so higher up means more statistically confident. Genes that
pass the filters are coloured — red on the right, blue on the left — and
everything else is grey. The dashed lines are the filters themselves.

**Summary statistics.** Worth reading carefully.

```{figure} ../../../_static/screenshots/analysis/de-fdr-denominator.png
:alt: The statistics table of an expanded differential expression result: Genes tested, Not tested, Significant at the adjusted-p and |log2FC| thresholds, Up in A, Down in A, the p-value convention, Multiple testing (Benjamini–Hochberg), the expression matrix and the fold change.
:width: 850px

The multiple-testing correction is applied over **genes tested**, not over the
size of the gene panel. Both numbers are shown so the false-discovery rate can be
reconstructed from this table alone.
```

- **Genes tested** — genes that produced a p-value. This is the number the
  correction divided by.
- **Not tested (fewer than 2 cells with a value)** — genes one of the two groups
  had fewer than 2 measured cells for. They have no p-value and are not part of
  the correction.
- **Significant** — tested genes passing **both** **Max adjusted p** and **Min
  |log2FC|**. The row label restates the two values in force.
- **Up in A / Down in A** — significant genes split by the sign of the fold
  change.
- **Multiple testing**, **Expression matrix**, **Fold change** — what the run
  used (and **p-values** for Wilcoxon, **Replicates** for a pseudobulk run).

**The gene table.** Every significant gene with its score, log₂ fold change,
adjusted p and the share of each group’s cells expressing it, highest score first:
the genes up in A lead, those up in B close the table.

In this run, the largest movers include `Iapp`, `Gcg` and `Ghrl`, all strongly
*lower* in Group A. That is the expected answer: Page 1 is a progenitor
population that has not yet switched on the hormone genes that define the mature
endocrine cells filling the rest.

---

## Step 7 — Change one setting and watch the answer move

Nothing is recomputed when you change a filter. The test statistics were fixed
the moment the run finished; the filters only decide which genes get called
significant, get coloured, and get listed. Changing them is instant and lossless,
in the form or at the top of the expanded view’s options.

Change **Min |log2FC|** and compare.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-a.png
:alt: A volcano plot at a log2 fold-change threshold of 0.5, with broad red and blue clouds of coloured points either side of a narrow grey band, and a legend below the axis title counting each colour.
:width: 1192px

**|log₂FC| ≥ 0.5** — a permissive cutoff. 1,588 genes are called: 750 up and
838 down.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-b.png
:alt: The same volcano plot at a log2 fold-change threshold of 3, with the threshold lines pushed far apart and almost every point now grey, and a legend below the axis title counting each colour.
:width: 1192px

**|log₂FC| ≥ 3** — a strict cutoff. The same run now reports 190 significant
genes instead of 1,588, all of them downregulated.
```
:::
::::

This is worth internalising: **"how many genes are significant" is a statement
about your thresholds as much as about your data.** If a volcano looks empty,
check the cutoffs before concluding there is nothing there — that is the single
most reported analysis surprise.

The other controls behave the same way:

- **Max adjusted p** — any value from 0 to 1 (0.05).
- **Min % expressing in either group** — hides genes expressed in few cells of
  both groups.
- **Labelled genes** — how many gene names the plot is allowed to draw; labels
  that would collide are dropped, and the plot says how many fit.
- **Point size**, **Show threshold lines**, **Color scheme** — appearance only.

---

## Step 8 — Take it with you

The **Export:** row sits at the top of the options column, in the expanded view
only.

```{figure} ../../../_static/screenshots/analysis/export-toolbar.png
:alt: The export row of an expanded analysis result: an Export label followed by PNG, SVG and CSV buttons.
:width: 438px

Three formats, the same three in every mode: **PNG**, **SVG**, **CSV**. There is
no PDF export.
```

- **PNG** — a raster image of the plot as drawn.
- **SVG** — the same plot as vector art, for figure assembly.
- **CSV** — the numbers. For Differential Expression the file is
  `differential_expression.csv`, with scanpy’s column names: `names, scores,
  logfoldchanges, pvals, pvals_adj, pct_nz_group, pct_nz_reference, auroc,
  mean_a, mean_b, n_a, n_b`.

:::{warning}
The CSV holds **every gene**, not just the ones passing your current filters,
and it does **not** record which groups A and B were, which test you chose, or
what the filters were. Write that down alongside the file
or you will not be able to reconstruct the run later. See
{doc}`09_exporting_analysis_results` for every mode's schema.
:::

To keep the setup and the result, save a session — see
{doc}`../l_sessions_sharing/index`. Settings are restored, and so are the
results of differential expression, marker genes, signatures, correlation, dot
plot and composition; Quick and Detailed views are recomputed from their settings.

---

## The same arc in the other modes

Steps 2 and 5–8 are identical everywhere. Only the middle changes.

| Mode | What you configure | What the preview shows | What the full view adds |
| --- | --- | --- | --- |
| {doc}`Quick <03_analysis_mode_quick_insights>` | nothing | composition bars and a statistics table | — (Quick has no expanded view) |
| {doc}`Detailed <04_analysis_mode_detailed_analysis>` | one variable, plot type, pages | the plot | summary table + statistical tests + export |
| {doc}`Correlation <05_analysis_mode_correlation_analysis>` | X, Y, pages; method and colour-by under More options | the scatter | per-page `r`, `r²`, `p`, `n` + strength reading + export |
| {doc}`Gene Signature <07_analysis_mode_gene_signature>` | gene sets, pages; method and display under More options | the distribution | per-page summary + between-page tests + gene chips + export |
| {doc}`Marker Genes <08_analysis_mode_genes_panel>` | Group by, compare against; groups, test, filters, ranking under More options | a ranked list, heatmap or clustered heatmap | the sortable per-group table, heatmap options + CSV, all genes CSV, GMT |

---

## Next steps

- The vocabulary underneath all of this: {doc}`01_analysis_mental_model`
- Where every control lives: {doc}`02_analysis_ui_overview`
- When a result looks wrong: {doc}`10_troubleshooting_analysis`
