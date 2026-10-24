# Analysis mode: Differential Expression (DE) (Group A vs Group B)

**Audience:** computational users + wet lab scientists doing marker discovery  
**Time:** 25–60 minutes  
**What you’ll learn:**
- How DE in Cellucid is defined (it compares **two groups of cells**: pages, or categories of a field)
- The exact statistics used (Wilcoxon rank-sum or Welch’s t-test, in scanpy’s and SciPy’s conventions)
- How log2 fold change is computed (as stored, or with `log1p` undone)
- How to read the volcano plot and the gene table without common mistakes
- How Cellucid sizes a run by itself

**Prerequisites:**
- A dataset loaded
- **Gene expression available** (var / genes)
- Two highlight pages, or a categorical field to take the groups from

---

## What DE is for (in Cellucid)

Differential Expression (DE) answers:

> “Which genes differ between Group A and Group B?”

Common uses:
- validate that a selected group has the expected marker genes,
- compare a condition against another inside one cell type (disease against healthy),
- generate a first-pass marker list for follow-up analysis.

DE in Cellucid is **exploratory**:
- it does not fit covariate models (batch, donor, condition),
- it does not replace a full DE pipeline when you need complex designs,
- it treats cells as independent observations (see {doc}`01_analysis_mental_model`).

---

## The form

Open **Analysis → Differential expression**. The essentials are always shown;
everything else is under **More options**, in collapsed groups whose headers
count what you changed from the defaults (see {doc}`02_analysis_ui_overview`).

| Where | Options |
| --- | --- |
| Essentials | **Compare** (Pages / Categories), **Group A**, **Group B**, **Swap A and B** |
| More options → Cells | **Within**, **Min cells per group** |
| More options → Test | **Replicates**, **Method**, **Wilcoxon p-values**, **Multiple testing**, **Expression matrix**, **Fold change from** |
| More options → Filters | **Max adjusted p**, **Min \|log2FC\|**, **Min % expressing in either group**, the gene filters |
| More options → Labels | **Labelled genes** |

Nothing runs until you press **Run**. When the groups cannot be compared, Run is
disabled and the line above it says why (see
[When Run is disabled](#when-run-is-disabled)). Changing a filter or the labels
after a run redraws the volcano at once; changing anything else marks the result
**Stale · Run** until you run again.

### Defaults

<!-- options:differential -->
| Label | Default |
|---|---|
| Compare | Pages |
| Group A | The first page |
| Group B | The second page |
| Field | The first categorical field with two or more categories |
| Group A | None |
| Group B | None |
| Within | All cells |
| Within field | The first categorical field |
| Within categories | All |
| Min cells per group | 10 |
| Replicates | None (cells) |
| Method | Wilcoxon rank-sum |
| Wilcoxon p-values | Exact or continuity-corrected (SciPy, R) |
| Multiple testing | Benjamini–Hochberg |
| Expression matrix | The mapped matrix |
| Fold change from | Automatic |
| Max adjusted p | 0.05 |
| Min \|log2FC\| | 1 |
| Min % expressing in either group | 0 |
| Exclude genes | None |
| Exclude pattern | None |
| Only genes in | All genes |
| Gene list | None |
| Labelled genes | 15 |
<!-- /options:differential -->

---

## Groups: pages or categories

**Compare → Pages.** Group A and Group B are highlight pages, or a derived
**Rest of \<page\>** (every cell not in that page), which gives one-vs-rest DE
without building a background page. Pages with no cells are not offered.

**Compare → Categories.** Choose a **Field**, then the categories that make up
Group A and those that make up Group B, each in a searchable checklist with its
cell count. Several categories in one group are pooled: Group A can be
`Disease` and Group B `Healthy`, or Group A three clusters and Group B two
others. A category cannot be in both groups.

**Swap A and B** exchanges the two groups. The sign of every fold change flips
with it.

### Within: compare inside some cells

**Within** (More options → Cells) keeps both groups inside a set of cells:

- **All cells** (default), **Visible cells** (what the current filters show), or
  any page;
- **Categories of a field**: choose a **Within field** and **Within
  categories**, for example `cell_type` = `T cell`, to compare disease against
  healthy inside T cells only.

Groups are always made of page membership or category labels; filters matter
only when **Within** is set to **Visible cells**. See
{doc}`01_analysis_mental_model` for visibility versus membership.

### Minimum group size: flagged, not refused

**Min cells per group** (10 by default, at least 2) does not stop a run. A group
smaller than it is tested all the same, and the result carries a **low n** chip
naming the group (the notification says the same). Only two things refuse a
run: overlapping groups, and a group with fewer than 2 cells.

Each gene is then tested when **both groups hold at least 2 cells with a
measured value** for it; a gene below that on either side is reported as not
tested (see
[Which genes the correction is applied to](#which-genes-the-correction-is-applied-to)).

### When Run is disabled

| Line above Run | Meaning |
| --- | --- |
| `Choose two pages` | A page group is empty in the form. |
| `Group A and Group B are the same page` | Choose another page for one side. |
| `Choose categories for Group A and Group B` | One side has no category. |
| `A category is in both groups: …` | Take it out of one side. |
| `The groups share N cells` | The pages overlap; a cell can be on one side only. |
| `Each group needs at least 2 cells` | No test says anything about one cell. |
| `Reading the categories…` | The field is still being read. |

---

## Statistics (exact implementation)

The test options (More options → Test) follow scanpy’s `rank_genes_groups`, so a
Cellucid result can be compared with one computed in Python column for column.

### Method

- **Wilcoxon rank-sum** (default). The score is the z of the rank sum,
  `(R − n(N+1)/2) / sqrt(n·m·(N+1)/12 · f)` — scanpy’s `scores` — where `n` and
  `m` are the group sizes, `N = n + m` and `f = 1 − Σ(t³ − t)/(N³ − N)` the tie
  factor. **AUROC** is `U / (n·m)`: the chance that a random cell of Group A
  ranks above a random cell of Group B.
- **t-test (Welch)**. The score is Welch’s t; means and variances are the exact
  sample moments of each group.
- **t-test, overestimated variance**. Welch’s t with Group B’s variance divided
  by Group A’s size instead of its own (scanpy’s `t-test_overestim_var`).

### Replicates (pseudobulk)

**Replicates → None (cells)** (default) treats every cell as an observation, as
scanpy’s `rank_genes_groups` does. Cells of one donor or sample are not
independent, so with many cells a cell-level test calls even small differences
significant. Choose a **sample field** instead and the run is a **pseudobulk**
test:

- each sample’s cells of Group A form one observation, and its cells of Group B
  another (a sample with cells on both sides gives one to each); cells with no
  sample are left out;
- a gene’s value for an observation is its mean over those cells, read in one pass
  over the gene’s non-zero entries;
- the test is **Welch’s t** between the two groups’ per-sample means, as
  `scipy.stats.ttest_ind(..., equal_var=False)` on those means; the score is t,
  `log2FC` follows **Fold change from** on the mean of the per-sample means, and
  `% A` / `% B` are the share of the group’s cells above 0.

Each group needs cells in at least 2 samples, else Run says so. **Method** and
**Wilcoxon p-values** do not apply and are hidden. No variance moderation (limma,
edgeR, DESeq2) and no pairing is applied; a gene with fewer than 2 samples with a
value on a side is not tested. The result line reads `Pseudobulk t-test`, and the
expanded view states the field and the samples on each side.

### Wilcoxon p-values

| Choice | p-value | Matches |
| --- | --- | --- |
| **Exact or continuity-corrected (SciPy, R)** (default) | exact when both groups have fewer than 50 values and no ties; otherwise normal, with tie and continuity correction | SciPy `mannwhitneyu(method='auto')`, R `wilcox.test` |
| **Normal approximation (scanpy)** | `2Φ(−\|z\|)` with `f = 1` | scanpy `rank_genes_groups` (default `tie_correct=False`) |
| **Normal, tie-corrected** | `2Φ(−\|z\|)` with the tie factor | scanpy with `tie_correct=True` |

The score uses the tie factor except under **Normal approximation**, as scanpy
does.

Implementation notes:
- The rank pass walks each gene’s non-zero values only: every zero is tied with
  every other zero, so the zeros are one rank block whose size is the zero
  count. Ranks, U, the tie factor and the p-value are exactly the numbers a pass
  over every cell gives. Nothing is binned or sampled.
- Both methods cost about the same: the work grows with a gene’s non-zero
  values, not with the number of cells.

### Expression matrix

The genes are read from the matrix the data mapping shows (**Mapped**, the
default), or from another matrix of the file — `X`, a layer, or `raw` with its own
genes — for this run alone, without changing what is shown. Every gene of the
chosen matrix is tested (after the gene filters). See
[Expression matrix](08_analysis_mode_genes_panel.md#expression-matrix) in Marker
genes for how each matrix names its genes and when only **Mapped** is offered.

### Fold change from

| Choice | `log2FC` |
| --- | --- |
| **Automatic** (default) | **Undo log1p** when the data mapping says the expression matrix read holds `log1p` values, otherwise **As stored** |
| **Undo log1p (scanpy)** | `log2((expm1(meanA) + 1e-9) / (expm1(meanB) + 1e-9))`, `expm1` in the base the matrix was logged with — scanpy’s `logfoldchanges` |
| **As stored** | `log2((meanA + 0.01) / (meanB + 0.01))` on the values as stored |

`meanA` and `meanB` are the means over each group’s cells. Whether a matrix holds
`log1p` values is set under **Session › Data mapping › Genes** (guessed from
`uns['log1p']`); see {doc}`../b_data_loading/13_data_mapping`.

Under **As stored**, scaled or z-scored expression can be negative, and a ratio
of negative means has no meaning. For such a gene both means are first shifted
by `m`, the most negative value over the cells of both groups:
`log2((meanA − m + 0.01) / (meanB − m + 0.01))`, so the fold change is always
finite and its sign is the sign of `meanA − meanB`. Read it as a direction and a
strength, not as a biological fold.

Interpretation:
- `log2FC > 0` → higher in **Group A**
- `log2FC < 0` → higher in **Group B**

### % expressing

`% A` and `% B` are the share of each group’s cells whose value is not zero
(scanpy’s `pct_nz_group` and `pct_nz_reference`).

### Multiple testing

**Benjamini–Hochberg** (default) or **Bonferroni** (`min(1, p · m)`), applied over
the genes tested.

### Which genes the correction is applied to

Either correction needs a denominator: the number of tests in the family. Cellucid
uses **only the genes that produced a p-value**.

A gene produces no p-value when one of the two groups has fewer than 2 cells with a
measured value for it. Such a gene was never tested, so it is excluded from the
denominator, its adjusted p-value stays empty, and it can never be counted as
significant. With 20,000 genes of which 8,000 could not be tested, the correction
runs at **m = 12,000**.

:::{important}
When you quote an FDR in a figure caption or a methods section, quote the number
of **genes tested** shown in the summary, not the size of your gene panel.
:::

For analysis-wide scope and assumptions, see {doc}`01_analysis_mental_model`.

---

## Filters and labels

The filters decide which tested genes are called significant, coloured and
listed; **Labelled genes** decides how many of them carry a name on the volcano.

| Option | A gene passes when |
| --- | --- |
| **Max adjusted p** | its adjusted p-value is ≤ the value (0.05) |
| **Min \|log2FC\|** | its absolute fold change is ≥ the value (1, a 2-fold change) |
| **Min % expressing in either group** | `% A` or `% B` is ≥ the value (0) |

None of these tests again: the statistics were fixed when the run finished, so
changing a filter — in the form or in the expanded view — redraws at once. That is
also why the same run can report very different counts.

The gene filters in the same group decide which genes are **tested**, so changing
them marks the result stale: **Exclude genes** (mitochondrial, ribosomal,
hemoglobin, each with how many genes it matches), **Exclude pattern** (a regular
expression on gene names), and **Only genes in** (all genes, the highly variable
genes when the data mapping names them, or a gene list).

---

## Outputs (what you see)

### The result line

Above the volcano a line says what was compared:
`Group A (n) vs Group B (n) · Wilcoxon · 3,753 genes`, ending with the matrix's
name when the run read one other than the mapped matrix, and a **low n** chip
when a group is under **Min cells per group**.

### While it runs

Run reads `Running…`. A box under it says what the run covers
(`3,753 genes in Beta (1,062 cells) vs Rest of Beta (2,634 cells)`), the step
(`Testing genes`, then `Applying multiple testing correction`), how many genes are
done, how much has been read and how long is left, and offers **Cancel**. See
{doc}`02_analysis_ui_overview`.

```{figure} ../../../_static/screenshots/analysis/de-cancelled.png
:alt: The Differential Expression panel after Cancel. The run button reads "Run differential expression" again; below it a box titled "Differential expression cancelled" reads "Cancelled after 2.3 s. Nothing was kept." with a Run again button.
:width: 464px

A cancelled run keeps nothing and is not an error. **Run again** runs the form
as it stands.
```

### Volcano plot

```{figure} ../../../_static/screenshots/analysis/de-sidebar-result.png
:alt: A small volcano plot in the sidebar with a red cloud on the right, a blue cloud on the left and a grey cloud in the middle, its axes labelled and no legend, and an "⤢ Expand" button below it with a pointer on it.
:width: 448px

The sidebar shows the volcano and its result line. The summary, the gene table
and the export live behind `⤢ Expand`.
```

Axes:
- **x-axis:** log2 fold change (A vs B)
- **y-axis:** `−log10(adjusted p)`

Coloring: up (positive log2FC, passes the filters), down (negative log2FC, passes
the filters), not significant.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-a.png
:alt: A volcano plot at a log2 fold-change threshold of 0.5, with broad red and blue clouds of coloured points either side of a narrow grey band, and a legend below the axis title reading Up (750), Down (841) and Not significant (2162).
:width: 1156px

**Min |log₂FC| 0.5** — 1,591 genes significant, 750 up and 841 down.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-b.png
:alt: The same volcano plot at a log2 fold-change threshold of 3, with the dashed threshold lines pushed far apart and nearly every point grey, and a legend below the axis title reading Up (0), Down (192) and Not significant (3561).
:width: 1156px

**Min |log₂FC| 3** — the same run, now 192 significant genes, all down.
```
:::
::::

### The expanded view

`⤢ Expand` opens the four-region view (see {doc}`02_analysis_ui_overview`):

- **Options:** the filters and **Labelled genes** at the top (the same values as
  the form), then the volcano’s look: point size, threshold lines, colour scheme.
- **Summary statistics:**

  | Row | What it counts |
  | --- | --- |
  | **Genes tested** | Genes that produced a p-value: the correction’s `m`. |
  | **Not tested (fewer than 2 cells with a value)** | Genes a group had fewer than 2 measured cells for. |
  | **Significant (adj p ≤ …, \|log2FC\| ≥ …)** | Tested genes passing the two filters; the label states the values in force. |
  | **Up in A** / **Down in A** | Significant genes with log2FC > 0 / < 0. |
  | **Multiple testing**, **Expression matrix**, **Fold change** | What the run used. |

- **Statistical analysis:** every significant gene with its score, log2FC,
  adjusted p, `% A` and `% B`, highest score first (so the genes up in A lead and
  those up in B close the table). Click a gene to colour the cells by it, as
  choosing it under **Coloring & filtering → Gene expression** does; the sidebar's
  gene search then names it.

To get every gene, including those that did not pass, export the CSV
({doc}`09_exporting_analysis_results`).

---

## How Cellucid sizes a run

There are no performance settings: every run sizes itself from what it measures,
and keeps adjusting while it goes. Nothing it decides changes a number in the
result — a fast run and a slow run of the same form return the same statistics.

- **Read-ahead.** A run reads genes ahead of the workers: as many as the workers
  consume while one gene is fetched, from the link’s measured round trip and
  throughput and each gene’s measured size and test time. A slow, distant link
  reads further ahead; a local file barely at all.
- **Workers.** Every analysis worker is used (one per CPU core, less one for the
  page and one for the file reader), each with one gene running and one waiting.
- **Memory.** Read-ahead is held to 30% of the analysis memory budget (a sixteenth
  of the device memory the browser reports, 512 MB on an 8 GB machine), and what
  the workers hold to half of it. Memory pressure, or an allocation the browser
  refuses, halves the read-ahead at once.
- **Adapting.** Every quarter second (and at least every 32 genes) the run checks
  itself: workers waiting on reads double the read-ahead; a doubling that did
  not make the run faster is undone and not tried again (past that point more
  reads only wait in the reader); a drop in throughput halves it.
- **Requests.** How many reads go out at once is the file reader’s: its scheduler
  paces them to what the link carries.

How a run reads the data:

- Two small groups of a large dataset are read **by cell**, every gene in one
  read, when that reads less than every gene’s column.
- Otherwise each gene’s column is read once and moved to a worker without a copy;
  the groups were sent to every worker once at the start, so the page never
  gathers a group’s values.

Measured with `npm run bench:analysis` on 2026-10-24 (Suo sample, 561,947
cells, 5,103 genes; Wilcoxon, every gene; a 10-core, 16 GB laptop, Chromium,
median of 5 cold runs): the two largest `LVL1` categories as pages in 1.0 s, the
largest against the rest in 1.1 s; over a simulated 4G link (170 ms round
trip, 9 Mbit/s) the two largest in 53 s, almost all of it reading the matrix.

For publication-scale designs (donors, batches, covariates), export and run a
replicate-aware DE pipeline offline.

---

## Edge cases and “gotchas”

- **Tiny groups**: tested and flagged **low n**; their p-values are unstable and
  many genes may be untestable.
- **All-zero or constant genes**: log2FC near 0 and p-values near 1 (expected).
- **Expression scale**: DE uses the matrix chosen under **Session › Data mapping**;
  check **Fold change from** when comparing with scanpy.
- **Missing gene expression**: the mode needs var / gene data.

---

## Troubleshooting (DE)

### Symptom: Run is disabled

Read the line above it; [When Run is disabled](#when-run-is-disabled) lists every
reason. With **Compare → Pages**, create two pages in Highlighted Cells (see
{doc}`../f_highlighting_selection/index`), or use **Rest of \<page\>** once you
have one.

### Symptom: “DE is slow”

The read of the expression matrix is almost all of a run. Narrow it:
- compare fewer cells (**Within**, or smaller pages),
- test fewer genes (**Only genes in** → highly variable genes, or a list),
- watch the running box: a run stuck on reading is limited by the link.

### Symptom: “Volcano plot looks wrong / missing points”

Common causes:
- the filters hide most genes (**Max adjusted p**, **Min |log2FC|**, **Min %**),
- many genes were not tested: read **Not tested** in the summary,
- A and B are the other way round: press **Swap A and B**.

---

## Next steps

- {doc}`07_analysis_mode_gene_signature` (score a curated gene set)
- {doc}`08_analysis_mode_genes_panel` (markers of every category at once)
- {doc}`09_exporting_analysis_results` (export DE tables/plots)
- {doc}`10_troubleshooting_analysis` (DE slow, volcano missing points)
