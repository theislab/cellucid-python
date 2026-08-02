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
  3,696 cells, 3,753 genes, 8 cell types. Pick it from **Sample datasets** in the
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

1. In **Coloring & Filtering**, choose a **Categorical obs** field —
   `cell_type` here. Annotation selection reads the clicked cell's value on the
   active field, so with no active field the tool never starts.
2. In **Highlighting**, leave the mode on **Annotation based**, then
   **Alt+click** any cell. Every cell sharing that cell's annotation value is
   picked up at once. The panel now reads `Step 1: 262 cells`.

```{figure} ../../../_static/screenshots/analysis/life-01-confirm-a-selection.png
:alt: The Highlighting panel with Annotation based selected, a line reading "Step 1: 262 cells", and Confirm, undo, redo and Cancel buttons above an empty page list labelled "Page 1 (0)". A pointer rests on Confirm.
:width: 540px

A selection is a **pending step** until you confirm it. Note that `Page 1` still
reads `(0)` and the group list below still says "No cells highlighted" — the 262
cells are staged, not saved.
```

Press **Confirm**. Only now does Page 1 acquire cells.

```{figure} ../../../_static/screenshots/analysis/life-02-two-pages-ready.png
:alt: The same Highlighting panel after confirming. The page tab reads "Page 1 (262)" and the group list below shows one entry, "Annotation (262 cells)", with an enable toggle and a remove button.
:width: 540px

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

Open the **Analysis** section. It lists six modes, all collapsed. Opening one
closes the others, so exactly one mode is active at a time.

```{figure} ../../../_static/screenshots/analysis/analysis-panel-modes.png
:alt: Six stacked rows reading Quick (Automatic insights), Detailed (Full control over options), Correlation (Explore variable relationships), Differential Expression (Find DE genes between groups), Gene Signature (Compute signature scores) and Marker Genes (Discover markers across groups). Each row carries a small copy-window icon and a chevron.
:width: 808px

The six modes. Each row carries a **copy** button, which undocks that mode into a
floating window, and a chevron, which opens it in place.
```

| Mode | Runs | Needs | Answers |
| --- | --- | --- | --- |
| **Quick** | automatically | one page with cells | "What is in this group?" |
| **Detailed** | automatically | one variable, one or more pages | "How does this variable differ between my groups?" |
| **Correlation** | automatically | two continuous or gene variables | "Do these two move together?" |
| **Differential Expression** | on **Run** | two different, non-empty pages, gene data | "Which genes separate group A from group B?" |
| **Gene Signature** | automatically | a list of gene names | "How strongly does this program score per cell?" |
| **Marker Genes** | on **Discover Markers** | a categorical field, gene data | "What marks every group at once?" |

Four modes recompute as you type or change a control. Two do not, because they
are expensive: **Differential Expression** waits for `Run Differential
Expression`, and **Marker Genes** waits for `Discover Markers`.

We will follow **Differential Expression**, because it is the mode with the most
moving parts. The other five are the same shape.

---

## Step 3 — Configure it

Open **Differential Expression**. The form asks for three things.

```{figure} ../../../_static/screenshots/analysis/de-page-selection.png
:alt: The Differential Expression form with a heading "Select pages to compare:", a Page A select set to "Page 1", the word "vs", a Page B select set to "Rest of Page 1", a Statistical method select reading Wilcoxon with the description "Rank-based test, robust to outliers and non-normal distributions.", a collapsed Performance Settings row, and a Run Differential Expression button.
:width: 488px

**Page A vs Page B.** The two must be different and both non-empty. The base
page and its `Rest of` complement are the simplest valid pair.
```

**Which page goes in A?** The sign convention follows A: a gene with positive
log₂ fold change is higher in **Page A**. Put the group you are characterising in
A and the background in B, and "upregulated" then means what you expect.

**Which statistical method?** Two, and only two:

- **Wilcoxon** (default) — a Mann-Whitney rank-sum test. It asks whether one
  group's values tend to rank above the other's. It does not assume a bell curve
  and is not thrown off by a few extreme cells. Use this unless you have a
  reason not to.
- **t-test** — Welch's t-test on the means. Faster, but it assumes the values
  are roughly normally distributed, which single-cell expression usually is not.

**Performance Settings** is collapsed by default and can stay that way. It exists
for very large datasets, where the limit is how fast gene vectors can be fetched
and held in memory rather than the arithmetic.

```{figure} ../../../_static/screenshots/analysis/de-performance-settings.png
:alt: The Performance Settings disclosure expanded, showing four selects: Batch size set to "500 genes (recommended)", Memory budget set to 512 MB, Network parallelism set to 12 parallel, and Compute parallelism set to Auto.
:width: 488px

The four performance controls, expanded. Every one of them has a working default;
they change how fast the run goes, never what it computes.
```

---

## Step 4 — Let it run

Press **Run Differential Expression**. The button reads `Running...` and a
progress bar appears with a named phase.

```{figure} ../../../_static/screenshots/analysis/de-running.png
:alt: The Differential Expression panel mid-run. The button is disabled and reads "Running...", and a progress bar above it shows the phase "Loading & Computing" with a gene count and an estimated time.
:width: 488px

Phase 1 of 2, `Loading & Computing`: every gene is fetched and tested. Phase 2,
`Multiple Testing Correction`, applies the false-discovery-rate adjustment across
all of them at once and is usually over in an instant.
```

Changing anything in the form while a run is in flight cancels it. That is
deliberate — a result that no longer matches the form on screen would be worse
than no result.

---

## Step 5 — Read the preview

When the run finishes, a small volcano plot appears below the form, with an
**⤢ Expand** button under it.

```{figure} ../../../_static/screenshots/analysis/de-sidebar-result.png
:alt: A small volcano plot in the sidebar with a red cloud on the right, a blue cloud on the left and a grey cloud in the middle, its axes labelled and no legend, and an "⤢ Expand" button below it with a pointer on it.
:width: 448px

The sidebar preview is a **glance**, not a workspace. At 224 pixels wide the plot
is confirmation that the run produced something; the numbers live in the full
view. There is no legend here because the band under the axis fits the tick
labels and the axis title *or* the legend, and an unlabelled axis is worse than
an unkeyed volcano.
```

:::{important}
The sidebar shows the plot only. Summary statistics, the ranked gene table, the
threshold controls and every export live **exclusively** in the expanded view.
If you are looking for a number and cannot find it, you have not expanded yet.
:::

The **⤢ Expand** button is the only way in. The preview itself is not clickable,
and the button is a real button: you can reach it with `Tab` and open the full
view with `Enter` or `Space`.

---

## Step 6 — Open the full view

```{figure} ../../../_static/screenshots/analysis/de-volcano-expanded.png
:alt: A large modal headed "DIFFERENTIAL EXPRESSION: PAGE 1 VS REST OF PAGE 1", divided into four regions: a volcano plot with labelled genes top-left over a legend reading Up (396), Down (701) and Not significant (2656) drawn below its axis title, an Export row (PNG, SVG, CSV) and a PLOT OPTIONS column on the right, a SUMMARY STATISTICS table bottom-left, and a STATISTICAL ANALYSIS panel bottom-right holding the Top Differentially Expressed Genes table.
:width: 1440px

The expanded view has four regions: the **plot**, the **options** column (which
also carries the export row), **summary statistics**, and **statistical
analysis**. Every mode uses this same four-region layout. The dotted seams between
them are drag handles.
```

Read it in this order.

**The volcano.** Horizontal position is log₂ fold change: right means higher in
Page A, left means higher in Page B. Vertical position is `−log₁₀` of the
p-value, so higher up means more statistically confident. Genes that clear both
thresholds are coloured — red on the right, blue on the left — and everything
else is grey. The dashed lines are the two thresholds themselves.

**Summary statistics.** Five rows, and they are worth reading carefully.

```{figure} ../../../_static/screenshots/analysis/de-fdr-denominator.png
:alt: A five-row table. Genes tested (FDR denominator) 3,753. Not tested (< 10 cells with a value) 0. Significant (FDR ≤ 0.05, |log₂FC| ≥ 1) 1,097. Upregulated 396. Downregulated 701.
:width: 906px

The multiple-testing correction is applied over **genes tested**, not over the
size of the gene panel. Both numbers are shown so the false-discovery rate can be
reconstructed from this table alone.
```

- **Genes tested (FDR denominator)** — genes that produced a p-value. This is the
  number the correction divided by.
- **Not tested (< 10 cells with a value)** — genes skipped before testing,
  because one of the two groups had fewer than 10 cells with a measured value for
  that gene. They have no p-value and no adjusted p-value, and they are not part
  of the correction. It reads `0` here because this export stores a value for
  every cell; on a sparser dataset it will not.
- **Significant** — tested genes passing **both** the current p-value threshold
  and the current fold-change threshold. The row label restates the two cutoffs
  in force, so it changes when you change them.
- **Upregulated / Downregulated** — significant genes split by the sign of the
  fold change, relative to Page A.

**The ranked table.** `Top Differentially Expressed Genes` lists significant genes
sorted by **absolute** fold change, largest first, so the strongest movers in
either direction come first. The `Top 5 / 10 / 20 / 100 / All` select changes how
many rows are drawn. The p-value column is headed `Adj. P-value` while the
FDR checkbox is ticked and `P-value` when it is not, and values below 0.001 are
printed as `<0.001` rather than as a long exponent.

In this run, the three largest movers are `Iapp` (−11.8), `Gcg` (−11.0) and
`Ghrl` (−10.7) — all strongly *lower* in Page A. That is the expected answer: Page
A is a progenitor population that has not yet switched on the hormone genes that
define the mature endocrine cells filling Page B.

---

## Step 7 — Change one setting and watch the answer move

Nothing in the expanded view is recomputed when you move a threshold. The test
statistics were fixed the moment the run finished; the thresholds only decide
which genes get called significant, get coloured, and get listed. Moving them is
instant and lossless.

Drag **log2 FC threshold** and compare.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-a.png
:alt: A volcano plot at a log2 fold-change threshold of 0.5, with large red and blue clouds of coloured points on both sides of a narrow grey band, and a legend below the axis title reading Up (750), Down (841) and Not significant (2162).
:width: 1160px

**|log₂FC| ≥ 0.5** — a permissive cutoff. 1,591 genes are called: 750 up and
841 down.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-b.png
:alt: The same volcano plot at a log2 fold-change threshold of 3, with the threshold lines pushed far apart and almost every point now grey, and a legend below the axis title reading Up (0), Down (192) and Not significant (3561).
:width: 1160px

**|log₂FC| ≥ 3** — a strict cutoff. The same run now reports 192 significant
genes instead of 1,591, all of them downregulated.
```
:::
::::

This is worth internalising: **"how many genes are significant" is a statement
about your thresholds as much as about your data.** If a volcano looks empty,
check the cutoffs before concluding there is nothing there — that is the single
most reported analysis surprise.

The other controls in the column behave the same way:

- **P-value threshold** — one of `0.001`, `0.01`, `0.05`, `0.1`.
- **Use FDR-adjusted p-values** — on by default. Untick it to threshold on raw
  p-values instead. Do not untick it to get more hits; with thousands of genes
  tested, raw p-values will hand you hundreds of false positives.
- **Maximum gene labels** — how many gene names the plot is allowed to draw. The
  plot prints how many it actually managed to fit — `7 of 15 gene labels shown`
  in the capture above — because labels that would collide are dropped.
- **Point size**, **Show threshold lines**, **Show legend**, **Color scheme** —
  appearance only.

---

## Step 8 — Take it with you

The **Export:** row sits at the top of the options column, in the expanded view
only.

```{figure} ../../../_static/screenshots/analysis/export-toolbar.png
:alt: A row labelled EXPORT: followed by three pill buttons reading PNG, SVG and CSV.
:width: 312px

Three formats, the same three in every mode: **PNG**, **SVG**, **CSV**. There is
no PDF export.
```

- **PNG** — a raster image of the plot as drawn.
- **SVG** — the same plot as vector art, for figure assembly.
- **CSV** — the numbers. For Differential Expression the file is
  `differential_expression.csv` and its columns are, in order:
  `gene, meanA, meanB, log2FoldChange, pValue, adjustedPValue`.

:::{warning}
The CSV holds **every gene that was tested**, not just the ones passing your
current thresholds, and it does **not** record which pages A and B were, which
method you chose, or what the thresholds were. Write that down alongside the file
or you will not be able to reconstruct the run later. See
{doc}`09_exporting_analysis_results` for every mode's schema.
:::

To keep the *setup* rather than the numbers, save a session — see
{doc}`../l_sessions_sharing/index`. Settings are restored; results are recomputed.

---

## The same arc in the other modes

Steps 2 and 5–8 are identical everywhere. Only the middle changes.

| Mode | What you configure | What the preview shows | What the full view adds |
| --- | --- | --- | --- |
| {doc}`Quick <03_analysis_mode_quick_insights>` | nothing | composition bars and a statistics table | — (Quick has no expanded view) |
| {doc}`Detailed <04_analysis_mode_detailed_analysis>` | one variable, plot type, pages | the plot | summary table + statistical tests + export |
| {doc}`Correlation <05_analysis_mode_correlation_analysis>` | X, Y, colour-by, method, pages | the scatter | per-page `r`, `r²`, `p`, `n` + strength reading + export |
| {doc}`Gene Signature <07_analysis_mode_gene_signature>` | gene list, scoring, normalisation, plot | the distribution | per-page summary + between-page tests + gene chips + export |
| {doc}`Marker Genes <08_analysis_mode_genes_panel>` | grouping field, mode, clustering | a heatmap, in all three modes | the ranked per-group table, a wider heatmap with dendrograms + export |

---

## Next steps

- The vocabulary underneath all of this: {doc}`01_analysis_mental_model`
- Where every control lives: {doc}`02_analysis_ui_overview`
- When a result looks wrong: {doc}`10_troubleshooting_analysis`
