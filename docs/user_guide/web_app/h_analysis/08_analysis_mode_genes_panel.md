# Analysis mode: Marker Genes (Genes Panel)

**Audience:** everyone (especially computational users; wet lab-friendly for interpretation)  
**Time:** 25–60 minutes  
**What you’ll learn:**
- How Marker Genes finds the markers of every category of a field at once
- How the reference, the groups and the cells are chosen
- What happens to groups that are too small
- How to read score, log2FC, adjusted p, % expressing and AUROC
- How filters and ranking re-apply without testing again

**Prerequisites:**
- A dataset loaded
- Gene expression available
- At least one **categorical obs** field with two or more categories (e.g., `cell_type` or `cluster`)

---

## What Marker Genes is for

Marker Genes is the “many groups” marker discovery tool. Use it when you want:
- the marker genes of *each* cell type or cluster of a field,
- a list or heatmap of those markers across groups,
- marker gene sets to score ({doc}`07_analysis_mode_gene_signature`) or to plot
  ({doc}`13_analysis_mode_dot_plot`).

How it differs from DE:
- {doc}`06_analysis_mode_differential_expression_de` compares **two groups** (A vs B).
- Marker Genes compares **every group** of a field against the rest, or against
  one reference group, in one run.

The tests are scanpy’s `rank_genes_groups`, with the same column names in the
CSV; see [Test](#test).

---

## The form

| Where | Options |
| --- | --- |
| Essentials | **Group by**, **Compare against** (Rest / Group), **Reference group** |
| More options → Groups and cells | **Groups to test**, **Cells**, **Rest includes**, **Min cells per group**, **Small groups** |
| More options → Test | **Method**, **Wilcoxon p-values**, **Multiple testing**, **Expression matrix**, **Fold change from** |
| More options → Filters | **Direction**, **Max adjusted p**, **Min \|log2FC\|**, **Min % expressing in group**, **Max % in reference**, the gene filters |
| More options → Ranking | **Rank by**, **Top genes per group** |

Nothing runs until you press **Run**. When nothing can be tested, Run is disabled
and the line above it says why (see
[When Run is disabled](#when-run-is-disabled)). Every option of **Direction**,
the filters above the gene filters, and **Ranking** re-applies to the shown
result at once; any other change marks it **Stale · Run**.

### Defaults

<!-- options:genesPanel -->
| Label | Default |
|---|---|
| Group by | The first categorical field with two or more categories |
| Compare against | Rest |
| Reference group | The largest group |
| Groups to test | All |
| Cells | All cells |
| Rest includes | All other cells |
| Min cells per group | 10 |
| Small groups | Skip and list |
| Method | Wilcoxon rank-sum |
| Wilcoxon p-values | Exact or continuity-corrected (SciPy, R) |
| Multiple testing | Benjamini–Hochberg |
| Expression matrix | The mapped matrix |
| Fold change from | Automatic |
| Direction | Up |
| Max adjusted p | 0.05 |
| Min \|log2FC\| | 0 |
| Min % expressing in group | 0 |
| Max % in reference | 100 |
| Exclude genes | None |
| Exclude pattern | None |
| Only genes in | All genes |
| Gene list | None |
| Rank by | Score (z or t) |
| Top genes per group | 10 |
<!-- /options:genesPanel -->

---

## Group by and compare against

**Group by** is a categorical field; each of its categories with cells is a group,
largest first. Cells whose category is missing belong to no group. An
out-of-range category code stops the run with a data-contract error.

**Compare against**:
- **Rest** (default): each group against every other cell (scanpy’s
  `reference='rest'`).
- **Group**: each group against one **Reference group** (scanpy’s
  `reference='<group>'`), by default the largest. Only the reference and the
  groups to test take part.

---

## Groups and cells

- **Groups to test** — a searchable checklist of the field’s groups with their
  cell counts (All / None / Invert). By default every group is tested.
- **Cells** — **All cells**, **Visible cells** (what the current filters show) or
  a page. The groups and their counts are taken inside these cells only.
- **Rest includes** — shown when comparing against the rest with a subset of
  groups: **All other cells** (default) or **Selected groups only** (the rest is
  the other chosen groups).

---

## Small groups

A group with fewer cells than **Min cells per group** (10 by default, at least 2)
is never an error. **Small groups** decides what happens to it:

| Choice | The small group |
| --- | --- |
| **Skip and list** (default) | is not tested and stays in the rest; the result shows **n skipped** and the notification names them |
| **Pool into one group** | is tested together with the other small groups as `Other small groups (k)` (**n pooled**); a pool still under the minimum is skipped |
| **Test anyway, flagged** | is tested when it has at least 2 cells, and marked **low n** |

A group is then compared against its reference for one gene when **both sides
hold at least `max(2, Min cells per group)` cells with a measured value** — 2
under **Test anyway**. When either side falls short, the gene is **not tested**
for that group: it has no p-value, is left out of that group’s correction, and
never passes a filter. Each group is its own family, so two groups of one run can
have different numbers of genes tested; the expanded view states each group’s
count.

### When Run is disabled

| Line above Run | Meaning |
| --- | --- |
| `This dataset has no categorical field to group by` | Marker Genes needs a categorical field. |
| `Reading the groups…` | The field is still being read. |
| `At least 2 groups with cells are needed to compare. Found: N.` | Too few groups in the chosen cells. |
| `Choose the reference group` / `The reference group has fewer than N cells` | Pick a larger reference. |
| `Choose groups to test` | Groups to test is empty. |
| `No group to test has N cells or more` | Lower Min cells per group, or choose Test anyway. |

---

## Test

The same options as Differential expression; formulas and conventions are in
{doc}`06_analysis_mode_differential_expression_de` (Statistics).

- **Method**: Wilcoxon rank-sum (default), t-test (Welch), or t-test with
  overestimated variance.
- **Wilcoxon p-values**: exact or continuity-corrected (SciPy, R; default), the
  normal approximation (scanpy’s default), or normal, tie-corrected.
- **Multiple testing**: Benjamini–Hochberg (default) or Bonferroni, over the genes
  tested in each group.
- **Expression matrix**: the matrix the genes are read from (see
  [Expression matrix](#expression-matrix)).
- **Fold change from**: Automatic (undo `log1p` when the data mapping says the
  matrix read holds `log1p` values), Undo log1p (scanpy), or As stored.

For each group and gene the result holds:

| Value | Meaning |
| --- | --- |
| **Score** | Wilcoxon: the z of the group’s rank sum (scanpy’s `scores`); t-test: Welch’s t |
| **log2FC** | Group against its reference, as **Fold change from** says |
| **adj p** | Adjusted p-value |
| **% in** / **% ref** | Share of the group’s / the reference’s cells whose value is not zero |
| **AUROC** | Wilcoxon: the chance a random group cell ranks above a random reference cell |

Each gene is read once and moved to an analysis worker without a copy; the worker
tests every group in one pass over the gene’s non-zero values, and also keeps
every gene’s mean in every group, so the heatmap needs no second read. A run sizes
itself; see {doc}`06_analysis_mode_differential_expression_de`
(How Cellucid sizes a run).

### Expression matrix

**Expression matrix** chooses, for this analysis alone, which matrix of the file
the genes are read from. Nothing shown changes: the embedding, the gene search and
every other analysis keep the matrix the data mapping shows.

| Choice | Reads |
| --- | --- |
| **Mapped (…)** (default) | the matrix the data mapping shows, named in brackets; follows the mapping when it changes |
| **X**, **layers['…']**, **raw.X · N genes** | that matrix, when it is not the mapped one |

Each matrix names its genes on its own axis: `X` and every layer share `var`;
`raw` has `raw.var`, which often holds genes `X` lacks (`X` is frequently the
highly variable subset). The names come from the mapping's gene-name column where
that axis has the column, else from its index. The genes tested are that matrix's
genes, after the gene filters; **Only genes in → Highly variable** keeps the genes
the mapped matrix marks. **Fold change from → Automatic** reads the chosen
matrix's own `log1p` record. The same option is in Differential expression and
Dot plot, where a gene the chosen matrix lacks is named as not found.

The list holds only **Mapped** when the file has no other usable matrix, or when
the open dataset cannot change its expression matrix; its (i) says so. Every way of
opening a dataset — a file opened here, a link, `cellucid serve`, a notebook — reads
the file in the browser as AnnData, so each offers the file's other matrices.

The result line names a matrix other than the mapped one
(`12 groups · Wilcoxon · BH · 33,538 genes · raw`), the expanded view's statistics name
the matrix read, and clicking a gene that only that matrix has says so instead of
colouring by it. Results and sessions record the choice.

---

## Filters

Filters decide which tested genes each group shows. They apply to the full
result, so changing one never tests again.

| Option | A gene passes when |
| --- | --- |
| **Direction** | Up: score above 0 (higher in the group, default) · Down: below 0 · Both |
| **Max adjusted p** | adj p ≤ the value (0.05) |
| **Min \|log2FC\|** | \|log2FC\| ≥ the value (0) |
| **Min % expressing in group** | % in ≥ the value (0) |
| **Max % in reference** | % ref ≤ the value (100) |

When no group has a gene left, the list says how many genes were tested and offers
**Relax filters**, which sets Direction to Both and the four thresholds to their
loosest.

### Excluded and selected genes

These decide which genes are **tested**, so changing them marks the result stale:

- **Exclude genes** — Mitochondrial (`MT-`), Ribosomal (`RPS`, `RPL`) and
  Hemoglobin (`HBA1`, `HBB`, …) genes, each with how many genes of the dataset it
  matches.
- **Exclude pattern** — a regular expression on gene names; an invalid one is
  marked in place.
- **Only genes in** — All genes, **Highly variable** (when the data mapping names
  a highly variable column), or **A gene list** (typed or pasted; names matched as
  in Gene signature).

---

## Ranking

- **Rank by**: Score (scanpy’s order, default), Adjusted p (smallest first),
  log2 fold change, AUROC, or % in group − % in reference. With Direction Down
  the most negative come first; with Both the largest magnitude.
- **Top genes per group**: how many markers each group shows (10).

Ties are broken by gene order in the dataset.

---

## Outputs (what you see)

Above the result a line says what was compared
(`12 groups · Wilcoxon · BH · 3,753 genes`, or `… vs <reference>`), with chips for
**skipped**, **pooled** and **low n** groups (hover for their names and sizes) and
**Cached** with **Recompute** when the answer came from memory. A run on an
expression matrix other than the mapped one ends the line with that matrix's name.

**Show** switches the view without testing again:
- **List** — each tested group and its markers, in rank order;
- **Heatmap** (default) — the shown markers (rows) against the tested groups
  (columns), from the group means the run kept;
- **Clustered** — the heatmap with genes and groups ordered by hierarchical
  clustering.

**Click a gene** — in the list or in the expanded view's table — to colour the
cells by it: the same as choosing it under **Coloring & filtering → Gene
expression**, whose search box then names it, so the sidebar and the embedding
always agree. A gene the shown expression matrix lacks (a run on another
[Expression matrix](#expression-matrix)) colours nothing; a notification says so.

The result also hands its markers on: **Score as signature** opens Gene signature
with one set per group, and **Dot plot** opens a dot plot of them by the same
field.

### Expanded view

`⤢ Expand` opens the four-region view:

- **Options:** the filters, ranking and Top genes at the top, then the heatmap’s
  look: Transform (z-score, log(1 + x), none), Cluster genes, Cluster groups,
  Distance (correlation, Euclidean, cosine), Linkage (average, complete, single),
  the dendrograms and the colours. Clustering choices re-cluster at once.
- **Summary statistics:** Group by, groups tested, method, p-values, multiple
  testing, expression matrix, fold change, genes tested, duration.
- **Statistical analysis:** one group at a time (choose it above the table), with
  `n of m tested genes pass the filters`, the markers sortable by any column,
  **Copy genes**, **All genes CSV** and **GMT**.

Clustering is exact agglomerative clustering, in double precision, in the
analysis workers: about 0.15 s for 3,000 genes × 69 groups, rows and columns, on
an M1 Pro. Clustering n genes holds one 8-byte distance per pair (about 4 MB for
1,000 genes, 256 MB for 8,192), which must fit the analysis memory budget (a
sixteenth of the device memory); a larger request is refused with the size it
needs. Show fewer genes, or turn that axis’s clustering off.

### Cache

A finished run is kept in memory (up to a quarter of the analysis memory budget,
least recently used first out). The same run — same field and grouping, cells,
groups, reference, test and expression matrix — is answered from it and marked **Cached**; press
**Recompute** to run it again. Changing the field’s categories, the data mapping
or the dataset makes it a different run. A cancelled run is never cached.

---

## Export

From the expanded view:

- **CSV** — the shown markers, one line per group and marker, with scanpy’s
  column names.
- **All genes CSV** — every gene of every group.
- **GMT** — each group’s shown markers as a gene set.

Columns and formats are in {doc}`09_exporting_analysis_results`.

---

## Edge cases and pitfalls

- **No markers found** — read the line above the table: when no gene could be
  tested for a group, the filters are not the cause; otherwise press **Relax
  filters**.
- **Markers of a mixture** — a group that mixes cell types finds the markers of
  its majority.
- **Scale** — compare with scanpy only under the same matrix and **Fold change
  from**; see {doc}`../b_data_loading/13_data_mapping`.

---

## Troubleshooting (Marker Genes)

### Symptom: Run is disabled

Read the line above it; [When Run is disabled](#when-run-is-disabled) lists every
reason.

### Symptom: a group is missing from the result

It was smaller than **Min cells per group** and skipped (see the **skipped**
chip). Choose **Pool into one group** or **Test anyway, flagged**, or lower the
minimum.

### Symptom: “Analysis is slow”

Almost all of a run is reading the expression matrix: test fewer genes (**Only
genes in**), or fewer cells (**Cells**). A large matrix stored column-major (CSC)
is read without a gene index pass first.

---

## Interface reference

```{figure} ../../../_static/screenshots/analysis/marker-genes-ranked.png
:alt: The Marker genes panel grouped by cell_type, with its Run button, and beneath it the result header, the Show choice (List, Heatmap, Clustered) on Heatmap, and a heatmap of eight cell-type rows against their top marker genes with a Z-score colour bar.
:width: 488px

The marker heatmap in the sidebar. Every group is named, but only a handful of
the gene names fit across the bottom, so expand when you need to read which gene
a column is.
```

```{figure} ../../../_static/screenshots/analysis/marker-genes-clustered.png
:alt: A heatmap with cell-type groups down one axis and gene names across the other, dendrograms attached to both axes, and a colour bar for the Z-scored expression.
:width: 1192px

**Clustered**, in the expanded view. Genes and groups are each reordered by
hierarchical clustering, and the two dendrograms are drawn along the axes.
```

---

## Next steps

- {doc}`06_analysis_mode_differential_expression_de` (two-group DE with volcano plot)
- {doc}`07_analysis_mode_gene_signature` (score the markers as gene sets)
- {doc}`09_exporting_analysis_results` (what each mode exports)
