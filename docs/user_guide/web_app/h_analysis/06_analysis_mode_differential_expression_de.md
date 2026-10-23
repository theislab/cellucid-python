# Analysis mode: Differential Expression (DE) (Page A vs Page B)

**Audience:** computational users + wet lab scientists doing marker discovery  
**Time:** 25–60 minutes  
**What you’ll learn:**
- How DE in Cellucid is defined (it compares **two highlight pages**)
- The exact statistics used (Wilcoxon/Mann–Whitney U or Welch’s t-test)
- How log2 fold change is computed (including pseudocount)
- How to read the volcano plot + top-genes table without common mistakes
- How to tune performance settings for large datasets

**Prerequisites:**
- A dataset loaded
- **Gene expression available** (var / genes)
- At least one highlight page (you need two “page options” to compare)

---

## What DE is for (in Cellucid)

Differential Expression (DE) answers:

> “Which genes differ between Page A and Page B?”

Where:
- **Page A** and **Page B** are highlight pages (or derived “Rest of \<page\>” wildcards).

Common uses:
- validate that a selected group has the expected marker genes,
- generate a first-pass marker list for follow-up analysis,
- sanity-check that an annotation is consistent with expression.

DE in Cellucid is **exploratory**:
- it does not fit covariate models (batch, donor, condition),
- it does not replace a full DE pipeline when you need complex designs,
- it assumes your page definitions are meaningful groups.

---

## Inputs (how you pick A and B)

In **Analysis → Differential Expression**, you select:

- **Page A** (dropdown)
- **Page B** (dropdown)

If you have many pages, the selector shows:
- base pages under **Pages**
- derived options under **Wildcards**, including **Rest of \<page\>**

```{figure} ../../../_static/screenshots/analysis/de-page-selection.png
:alt: The Differential Expression form with the heading "Select pages to compare:", a Page A select showing "Page 1", the word "vs", a Page B select showing "Rest of Page 1", a Statistical method select reading Wilcoxon with its description line, a collapsed Performance Settings row, and a Run Differential Expression button.
:width: 488px

The configured form. Pages with **zero** cells are listed but disabled, so an
invalid pair cannot be assembled by accident.
```

Until two non-empty options exist, the selector replaces itself with a notice:

```{figure} ../../../_static/screenshots/analysis/de-needs-two-groups.png
:alt: The Differential Expression panel showing the notice "Differential expression requires two non-empty cell groups." above two greyed rows reading "Page 1 (0 cells)" and "Rest of Page 1 (3,696 cells)", with the Statistical method select and the Run button still below.
:width: 488px

`Differential expression requires two non-empty cell groups.` — shown whenever
fewer than two of the available pages have any cells in them. The counts beside
each row tell you which one is empty.
```

### What “Rest of \<page\>” means here

`Rest of Page A` means “all cells not in Page A”.
This enables one-vs-rest DE without manually building a background page.

### Filtering note

DE uses page membership (stored cell indices), not canvas visibility.
If you filtered cells out after creating a page, they may still be included in DE.
See {doc}`01_analysis_mental_model` for the visibility vs membership distinction.

### Minimum group size

DE requires enough cells in both groups:
- the implementation enforces a minimum of **10 cells** per group (default).
- the run is refused before any gene is loaded when either chosen group has fewer
  cells than that: a page of five cells stops with
  `Differential expression requires at least 10 cells per group; received 5 and …`.
  A page with zero cells cannot be selected at all.
- the same minimum is applied again **per gene**, on the cells that actually
  carry a measured value for that gene. A gene that falls below the minimum on
  either side is skipped: it gets no p-value, no adjusted p-value, and is
  excluded from the FDR correction. The summary reports how many genes this
  happened to — see
  [Which genes the correction is applied to](#which-genes-the-correction-is-applied-to).

---

## Statistics (exact implementation)

Cellucid supports two statistical methods:

### 1) Wilcoxon (default) = Mann–Whitney U

- This is a rank-based, non-parametric two-sample test.
- It compares the distribution of expression values between the two groups.

Implementation notes:
- Cellucid computes the U statistic from a full, tie-aware rank pass over the two selected groups.
- The rank pass walks each gene's non-zero values only: every zero is tied with every other zero, so the zeros are one rank block whose size is the zero count. The ranks, the U statistic, the tie correction and the p-value are exactly the numbers a pass over every cell gives.
- The U statistic is not histogram-binned or sampled.
- The p-value is exact when both comparison groups contain fewer than 50 finite values and there are no ties.
- Otherwise, the p-value uses the exact U statistic with a tie- and continuity-corrected asymptotic reference distribution.

Practical implication:
- large groups retain the exact U statistic, but use asymptotic inference. Both tests cost about the same: the work grows with a gene's non-zero values, not with the number of cells.

### 2) t-test = Welch’s t-test

- Parametric test comparing means.
- Uses Welch’s formulation (does not assume equal variances).

It is less robust than Wilcoxon to heavy non-normality and outliers. Its means and variances are the exact two-pass sample moments of each group's values in cell order.

For analysis-wide scope and assumptions, see {doc}`01_analysis_mental_model`.

---

## Effect size: log2 fold change (log2FC)

For each gene, Cellucid computes:

`log2FC = log2((meanA + 0.01) / (meanB + 0.01))`

Where:
- `meanA` is the mean expression in **Page A**
- `meanB` is the mean expression in **Page B**
- `0.01` is a fixed pseudocount to avoid division by zero

That is the whole formula whenever the compared cells are non-negative, which
is every count, normalized or log1p matrix. Scaled or z-scored expression can
be negative, and a ratio of negative means has no meaning (it is undefined when
the means straddle `-0.01` and its sign inverts when both are negative). For
such a gene both means are first shifted by `m`, the most negative value over
the cells of both pages:

`log2FC = log2((meanA − m + 0.01) / (meanB − m + 0.01))`

so the fold change is always finite and its sign is always the sign of
`meanA − meanB`. Its magnitude depends on that shift, so read it as a ranking
of direction and strength rather than a biological fold for scaled data.

Interpretation:
- `log2FC > 0` → gene is higher in **Page A**
- `log2FC < 0` → gene is higher in **Page B**

:::{important}
The sign depends on the A/B ordering.

If the biology you expect is “markers of Page B”, consider swapping A and B so that “upregulated” corresponds to Page B.
:::

---

## Multiple testing correction (FDR)

Cellucid computes a Benjamini–Hochberg adjusted p-value per gene:
- this is exposed as **adjustedPValue** and visualized when “Use FDR-adjusted p-values” is enabled in the volcano plot options.

Default behavior:
- the volcano plot uses **adjusted p-values** by default,
- and many summary counts are framed as “FDR < 0.05”.

### Which genes the correction is applied to

Benjamini–Hochberg needs a denominator: the number of tests in the family being
corrected. Cellucid uses **only the genes that produced a p-value**.

A gene produces no p-value when one of the two groups has fewer than the
minimum number of cells with a measured value for that gene (10 by default, see
[Minimum group size](#minimum-group-size)). Such a gene was never tested, so it
is excluded from the denominator, its **adjustedPValue** stays empty, and it can
never be counted as significant.

This matters whenever your pages are small or your panel is sparse: with 20,000
panel genes of which 8,000 fail the minimum cell check, the correction runs at
**m = 12,000**, not 20,000. Correcting at 20,000 would have produced adjusted
p-values roughly 1.7× larger.

:::{important}
When you quote an FDR in a figure caption or a methods section, quote the number
of **genes tested** — the summary row labelled *Genes tested (FDR denominator)* —
not the size of your gene panel. They are the same number only when every gene
passed the minimum cell check.
:::

---

## Outputs (what you see)

### Summary stats (quick sanity check)

You’ll see:

| Row | What it counts |
| --- | --- |
| **Genes tested (FDR denominator)** | Genes that produced a p-value. This is the `m` Benjamini–Hochberg used. |
| **Not tested (< N cells with a value)** | Genes skipped before testing because a group had fewer than `N` cells with a measured value (`N` = the minimum group size, 10 by default). These have no p-value and no adjusted p-value. |
| **Significant (FDR ≤ 0.05, \|log₂FC\| ≥ 1)** | Tested genes passing the current p-value **and** fold-change thresholds. The row label is rebuilt from the two cutoffs in force, so it changes when you change them. Always a subset of *Genes tested*. |
| **Upregulated** | Significant genes with log2FC > 0. |
| **Downregulated** | Significant genes with log2FC < 0. |

*Genes tested* + *Not tested* is the size of the gene panel that was scanned.

Upregulated/downregulated are defined by the sign of log2FC (Page A relative to Page B).

```{figure} ../../../_static/screenshots/analysis/de-fdr-denominator.png
:alt: A five-row table reading Genes tested (FDR denominator) 3,753; Not tested (< 10 cells with a value) 0; Significant (FDR ≤ 0.05, |log₂FC| ≥ 1) 1,097; Upregulated 396; Downregulated 701.
:width: 906px

Both halves of the panel are reported, which is what makes the false-discovery
rate reconstructable from this table alone. *Not tested* reads `0` here because
this dataset stores a value for every cell; on a sparser matrix it will not.
```

### While it runs

```{figure} ../../../_static/screenshots/analysis/de-running.png
:alt: The Differential Expression panel mid-run. The button reads "Running…" and is disabled; under it the cost line "Tests 2,000 genes in 1,000 vs 3,000 cells.", then a box titled "Differential expression" with "2,000 genes in Beta cells (1,000 cells) vs Rest of Beta cells (3,000 cells)", "Testing genes: 0 of 2,000 genes", a progress bar, and a Cancel button.
:width: 464px

The box under the button says what the run covers, which step runs
(`Testing genes`, then `Applying multiple testing correction`), how many genes
are done of how many, how much has been read and how long is left. **Cancel**
stops the run; so does changing any control.
```

```{figure} ../../../_static/screenshots/analysis/de-cancelled.png
:alt: The Differential Expression panel after Cancel. The run button reads "Run differential expression" again; below it a box titled "Differential expression cancelled" reads "Cancelled after 2.3 s. Nothing was kept." with a Run again button.
:width: 464px

A cancelled run keeps nothing and is not an error. **Run again** runs the form
as it stands. See {doc}`02_analysis_ui_overview` for what Cancel stops.
```

### Volcano plot

```{figure} ../../../_static/screenshots/analysis/de-sidebar-result.png
:alt: A small volcano plot in the sidebar with a red cloud on the right, a blue cloud on the left and a grey cloud in the middle, its axes labelled and no legend, and an "⤢ Expand" button below it with a pointer on it.
:width: 448px

The sidebar shows the volcano **and nothing else**. Summary statistics, the ranked
gene table and every threshold control live behind `⤢ Expand`. The legend is
absent on purpose: at this width the band under the axis has room for the tick
labels and the axis title or for the legend, and the labels win. Set **Show
legend** to any explicit position in the expanded view and it is drawn there
regardless.
```

```{figure} ../../../_static/screenshots/analysis/de-volcano-expanded.png
:alt: The expanded Differential Expression modal headed "DIFFERENTIAL EXPRESSION: PAGE 1 VS REST OF PAGE 1" with a labelled volcano plot over a legend reading Up (396), Down (701) and Not significant (2656) drawn below its axis title, an Export row of PNG SVG CSV, a PLOT OPTIONS column holding the p-value threshold, log2 FC threshold, FDR checkbox and label controls, a SUMMARY STATISTICS table, and a Top Differentially Expressed Genes table with a Top 5 select.
:width: 1440px

The expanded view, where DE is actually read. Four regions: plot, options and
export, summary statistics, and the ranked gene table. There is room for the
legend here, so it is drawn — below the axis title, never across it.
```

Axes:
- **x-axis:** log2 fold change (A vs B)
- **y-axis:** `-log10(p)` (raw or adjusted, depending on the toggle)

Coloring (conceptual):
- significant up (positive log2FC + passes thresholds)
- significant down (negative log2FC + passes thresholds)
- not significant (fails thresholds)

Threshold controls (in the plot options):
- p-value threshold (0.001 / 0.01 / 0.05 / 0.1)
- log2FC threshold (slider, 0 to 3 in steps of 0.25; default 1.0 = 2-fold)
- use adjusted vs raw p-values (**on** by default)
- maximum gene labels (0 to 50 in steps of 5; default 15)
- point size, threshold lines, legend position, colour scheme

None of these re-runs the test. The statistics were fixed when the run finished;
the thresholds only decide which genes are called significant, coloured and
listed. That is why moving them is instant — and why the same run can report very
different counts.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-a.png
:alt: A volcano plot at a log2 fold-change threshold of 0.5, with broad red and blue clouds of coloured points either side of a narrow grey band, and a legend below the axis title reading Up (750), Down (841) and Not significant (2162).
:width: 1156px

**|log₂FC| ≥ 0.5** — 1,591 genes significant, 750 up and 841 down.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-b.png
:alt: The same volcano plot at a log2 fold-change threshold of 3, with the dashed threshold lines pushed far apart and nearly every point grey, and a legend below the axis title reading Up (0), Down (192) and Not significant (3561).
:width: 1156px

**|log₂FC| ≥ 3** — the same run, now 192 significant genes, all downregulated.
```
:::
::::

:::{warning}
Do not untick **Use FDR-adjusted p-values** to get more hits. With thousands of
genes tested at once, raw p-values below 0.05 will include hundreds of genes that
are there by chance alone — that is exactly what the adjustment exists to remove.
:::

### Top genes table (expanded view only)

“Top Differentially Expressed Genes” is drawn in the **expanded view**, never in
the sidebar — the sidebar holds the volcano and the `⤢ Expand` button and nothing
else. The table is:
- filtered by the current p-value threshold (adjusted while **Use FDR-adjusted
  p-values** is ticked, raw when it is not),
- sorted by **absolute** log2FC, so the strongest movers in either direction come
  first,
- limited to the row count chosen in the `Top 5 / 10 / 20 / 100 / All` select
  beside it.

`All` shows every significant gene. To get the genes that did *not* pass the
thresholds as well, export the CSV.

---

## Performance settings (when DE is slow)

DE is gene-heavy: it can touch thousands of genes and hundreds of thousands of cells.
```{figure} ../../../_static/screenshots/analysis/de-performance-settings.png
:alt: The Performance Settings disclosure expanded under the Differential Expression form, showing four selects: Batch size set to "500 genes (recommended)", Memory budget set to 512 MB, Network parallelism set to 12 parallel, and Compute parallelism set to Auto.
:width: 488px

The four performance controls. They change how fast a run goes and never what it
computes, so a slow run and a fast run of the same configuration return the same
numbers.
```

Use **Performance Settings** (collapsible) to tune:

- **Batch size**: how many genes to preload ahead of compute (higher can be faster but uses more memory).
- **Memory budget**: limits how many gene vectors can be in flight without crashing the browser. The default is a sixteenth of the device memory the browser lets the page estimate (512 MB on an 8 GB machine, 1 GB on a 16 GB one); the choices run from a thirty-second to half of it.
- **Network parallelism**: how many gene loads a run keeps asked for at once; the reader reads the genes asked for together as one batch.
- **Compute parallelism**: how many genes are computed at once. *Auto* uses every analysis worker: one per CPU core, less one core for the page and one for the reader that feeds the workers (8 workers on a 10-core machine). The choices go up to that number.

How a run is computed, and what bounds it:

- Two named groups of a large dataset are read **by cell**, every gene in one read, when that reads less than every gene's column. The genes are then tested in runs, a run of genes per worker.
- Otherwise each gene's column is read once and moved to a worker without a copy. The two groups were sent to every worker once at the start of the run, so the page never gathers a group's values, including for a *rest of* page.
- The tests walk each gene's non-zero values only (see the method notes above), so testing is a small share of a run: what remains is reading the genes.

Measured on a 561,947-cell, 8,192-gene dataset on a 10-core laptop (Wilcoxon, every gene): two cell types of 4,737 and 3,995 cells in 0.34 s; the two largest types (55,096 and 48,640 cells) in 3.7 s; one 266-cell type against the rest of the dataset in 2.1 s. In each case almost all of the time is the read of the expression matrix.

Recommended workflow for very large datasets:
- start with smaller pages (more targeted groups),
- export results and do full DE offline for publication workflows.

---

## Edge cases and “gotchas”

- **Pages overlap** (same cell in A and B): DE is not meaningful; fix page definitions.
- **Tiny pages**: unstable means and unreliable p-values; many genes will fail minimum cell checks.
- **All-zero/constant genes**: log2FC ~ 0 and p-values ~ 1 (expected).
- **Expression scale ambiguity**: DE uses whatever scale you exported (counts vs log1p vs normalized).
- **Missing gene expression**: mode will fail or be empty if var/gene data is unavailable.

---

## Troubleshooting (DE)

### Symptom: “Differential expression requires two non-empty cell groups.”

This notice replaces the page selector when fewer than two of the available pages
have any cells. Pressing **Run** with an incomplete pair instead reports
`Select two non-empty cell groups to run differential expression`.

Cause:
- you don’t have enough non-empty pages/wildcards to populate A and B.

Fix:
- create highlight pages in Highlighted Cells (see {doc}`../f_highlighting_selection/index`),
- or use `Rest of \<page\>` once you have at least one base page.

### Symptom: “DE is extremely slow”

Likely causes:
- huge dataset (n_cells) and/or many genes,
- exact Wilcoxon ranking costs,
- insufficient memory budget leading to throttling.

How to confirm:
- try a smaller page (fewer cells),
- switch method to t-test and compare runtime,
- watch the progress phases; if it stalls on loading, it’s an I/O issue.

Fix:
- reduce scope (smaller pages),
- adjust Performance Settings (lower batch size and/or memory budget to avoid crashes; increase cautiously for speed),
- prefer offline DE for publication-scale runs.

### Symptom: “Volcano plot looks wrong / missing points”

Common causes:
- thresholds hide most points (p-value threshold or log2FC threshold),
- adjusted p-values are much larger than raw (expected after BH),
- many genes have NaN due to insufficient valid cells.

How to confirm the last one:
- read the **Not tested (< N cells with a value)** row in the summary stats. If
  it is a large fraction of the panel, most of your genes never entered the
  analysis and the volcano plot has nothing to draw for them.

Fix:
- lower log2FC threshold,
- switch between adjusted and raw p-values to understand correction effects,
- increase group sizes / fix page overlap.

---

## Next steps

- {doc}`07_analysis_mode_gene_signature` (score a curated gene set)
- {doc}`09_exporting_analysis_results` (export DE tables/plots)
- {doc}`10_troubleshooting_analysis` (DE slow, volcano missing points)
