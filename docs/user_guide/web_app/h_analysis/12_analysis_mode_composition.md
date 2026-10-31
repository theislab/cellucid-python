# Analysis mode: Composition

**Audience:** wet lab + computational users (differential abundance of cell types)
**What you’ll learn:**
- How Composition counts cells per cell type, sample and condition
- Which test it runs, and why samples (not cells) are the replicates
- What the plots, the table and the CSV files hold

**Prerequisites:**
- A dataset loaded
- A categorical field of cell types; for a test, a field of samples and a field of conditions

---

## What Composition answers

*Do cell-type proportions differ between conditions, accounting for samples?*
For example: is the share of regulatory T cells higher in disease than in
healthy donors?

Composition counts the cells of every cell type in every sample (or in every
condition, without samples), shows the proportions as stacked bars, and — with
a condition — tests every cell type for a difference between the conditions.

It reads no gene expression: it needs only the categorical fields, and a run
is one pass over the cells in a worker.

---

## Cell types, samples and condition

The three fields above **Run**:

- **Cell types**: the field whose categories are counted. The default is the
  field the cells are colored by, when it is categorical and not named like a
  sample; otherwise the first
  categorical field named like a cell type or cluster (`cell_type`,
  `annotation`, `leiden`, `louvain`, `cluster`, `label`, `LVL1`…); otherwise
  the first not named like a sample.
- **Samples**: the field that names each cell’s sample, donor or library — the
  biological replicate. The default is the first other field whose name
  contains `sample`, `donor`, `patient`, `individual`, `library` or `batch`;
  otherwise None.
- **Condition**: the field whose groups are compared (disease, treatment,
  timepoint). None by default, which shows the composition without a test.

Each role takes a different field: a field chosen for one is not offered for
the others.

A sample whose cells fall into two conditions is counted once per condition
(as two units); the result says so with a `split samples` chip.

```{figure} ../../../_static/screenshots/analysis/composition-result.png
:alt: The Composition panel on the Kanemaru heart sample with Cell types LVL1, Samples sample_ID and Condition sex, and the stacked bars beneath: one bar per sample, grouped under its condition, split by cell type.
:width: 488px

The Kanemaru heart sample: **Cell types** `LVL1`, **Samples** `sample_ID` and
**Condition** `sex`. The line under **Run** says what was counted and tested;
the stacked bars beneath show one bar per sample, grouped under its condition.
```

---

## Cells

Under **More options → Cells**:

- **Cells**: all cells, the cells the view shows (its filters applied), or the
  cells of one highlight page.
- **Only these cell types**: count only the chosen cell types. The proportions
  are then of the chosen types (they add up to 100% without the others), which
  is what you want when comparing, say, T-cell subsets within T cells.

Cells with a missing value in any of the three fields are not counted.

---

## Tests

With a samples field and a condition, each cell type is tested across the
samples of each condition. **More options → Test** chooses how:

| Test | Values per sample | Two conditions | More conditions |
|---|---|---|---|
| Logit proportions (propeller) — default | `logit((count + 0.5) / (total + 0.5·K))`, K the types with cells | Welch’s t-test | one-way ANOVA |
| Ranks of proportions | `count / total` | Wilcoxon rank-sum | Kruskal-Wallis |
| Centred log-ratio | `log(count + 0.5) − mean over types of log(count + 0.5)` | Welch’s t-test | one-way ANOVA |

- The logit test follows propeller (Phipson et al., *Bioinformatics* 2022), on
  its logit-transformed pseudo-proportions exactly (speckle's
  `getTransformedProps`: `(count + 0.5) / rowSums(count + 0.5)`, the total
  plus 0.5 a type). propeller then moderates the variances with limma’s
  empirical Bayes; Cellucid does not, so with few samples its p-values are
  more conservative than propeller’s.
- The Wilcoxon rank-sum test is exact when both conditions have fewer than 50
  samples and no ties, otherwise asymptotic with tie and continuity
  correction (as R’s `wilcox.test`).
- **Multiple testing**: Benjamini-Hochberg (default) or Bonferroni over the
  cell types.
- A test across samples needs at least two samples in every condition; with
  fewer, the type gets no p-value and the result shows `Too few samples`.
- A cell type with **no cells in scope** (a category the chosen cells, a page
  or a subsample left out) is not tested: its transformed proportions would
  depend only on how many cells each sample has, so it would look changed
  whenever conditions differ in size. It stays out of the centred log-ratio's
  mean, the tests and the multiple-testing correction; the table reads
  `No cells in scope` for it, and the result carries a
  `N types without cells` chip.

These match SciPy (`ttest_ind(equal_var=False)`, `f_oneway`, `mannwhitneyu`,
`kruskal`, `false_discovery_control`) on the same values; the parity is held by
the web app’s test suite.

---

## Without samples

With a condition but no samples field, the cells are the only replicates.
Composition then runs a χ² test of cell type × condition and, for each type, a
χ² test of that type against all others across the conditions (Fisher’s exact
test for a 2 × 2 table with an expected count below 5, two-sided as R's
`fisher.test`; Yates’s correction on a 2 × 2 χ², as
`scipy.stats.chi2_contingency` and R's `chisq.test`). These are the Detailed
view's tests too, from one kernel, so one table gets one answer in both.

Treating cells as replicates ignores how much samples differ from one another,
so these p-values overstate the evidence. The result carries a
`Cells as replicates` chip; choose a samples field whenever the data has one.

---

## Display

Under **More options → Display** (these redraw the result without counting
again):

- **Bars per**: one bar per sample (grouped under its condition) or one bar per
  condition.
- **Cell type order**: as in the field, most abundant first, or largest change
  (smallest adjusted p) first.
- **Show**: proportions or cell counts.

---

## Results and export

- **Sidebar**: the stacked bars, and chips for what limits the test.
- **Expanded view** (⤢): the stacked bars or, with samples, **Per sample** — one
  point per sample and type over a box per condition, with `*` (adj. p ≤ 0.05),
  `**` (≤ 0.01) and `***` (≤ 0.001) above the types. The summary lists the
  cells, types, samples, conditions and test; the table lists per type the mean
  proportion in each condition, the log₂ ratio of the second condition’s mean
  over the first’s (two conditions), the statistic, p and adjusted p.
- **CSV** (in the expanded view) downloads `composition_results.csv` when
  there is a condition: `cell_type`, `mean_proportion_<condition>` per
  condition, `log2_ratio`, `test`, `statistic`, `pvalue`, `pvalue_adjusted`,
  `n_samples_<condition>` per condition. Without a condition it downloads the
  counts.
- **Counts CSV** (under the summary) downloads `composition_counts.csv`:
  `cell_type` and one column of cell counts per sample (or condition), ready
  for propeller, scCODA or any other tool.

```{figure} ../../../_static/screenshots/analysis/composition-expanded.png
:alt: The expanded Composition view on the Kanemaru heart sample: the stacked bars per sample by sex, the options beside them, and the table of each cell type's proportions and test result.
:width: 1440px

The same result expanded: the cell-type legend above the stacked bars, the
summary (cells, cell types, samples, conditions, test, multiple testing) with
**Counts CSV**, and per cell type the mean proportion in each condition, the
log₂ ratio and the statistic.
```

Sessions save the options and the shown result.

---

## Options and defaults

<!-- options:composition -->
| Label | Default |
|---|---|
| Cell types | The coloring field when it is not named like a sample, else the first categorical field named like a cell type or cluster, else the first not named like a sample |
| Samples | The first other field named like sample, donor, patient, individual, library or batch; else None |
| Condition | None |
| Cells | All cells |
| Only these cell types | All |
| Test | Logit proportions (propeller) |
| Multiple testing | BH |
| Bars per | Sample |
| Cell type order | As in the field |
| Show | Proportions |
<!-- /options:composition -->

---

## Limits and edge cases

- On a subsample, the result says how the open cells were chosen
  ({doc}`../b_data_loading/21_subsampled_loading`). Composition does not run
  when the open cells are stratified by the cell-type field at unequal shares
  (balanced, or with a minimum per category): the proportions would be the
  subsample's design, not the data's. Open every cell, or a uniform or
  proportional subsample.

- Composition compares proportions, which are linked: when one type grows the
  others shrink. A significant change in one type can be the reflection of a
  change in another. The centred log-ratio test reads every type relative to
  the sample’s geometric mean, which is one way to see this.
- Neighbourhood-level differential abundance (Milo) and Bayesian models
  (scCODA) are not in the web app; export the counts to run them.
