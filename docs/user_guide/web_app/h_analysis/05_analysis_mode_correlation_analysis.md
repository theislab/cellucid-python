# Analysis mode: Correlation (X vs Y across pages)

**Audience:** computational users (still usable for wet lab with guidance)  
**Time:** 15–35 minutes  
**What you’ll learn:**
- What Cellucid’s correlation computes (Pearson/Spearman + p-value)
- How scope works (correlation is computed **per selected page**)
- How missing values are handled (paired values only)
- How plotting differs from computation on very large pages (plot downsampling)

**Prerequisites:**
- A dataset loaded
- At least one highlight page (recommended)
- For gene-based correlations: gene expression available

---

## What Correlation mode is for

Correlation answers the question:

> “Across cells in a page, do these two variables move together?”

Examples:
- `CXCL8` vs `S100A8` (gene–gene co-expression)
- `pct_mito` vs `n_counts` (QC relationships)
- `CXCL8` vs `neutrophil_score` (gene vs score)

What it is *not*:
- causality,
- a model that controls for batch/covariates,
- a replacement for regression/GLMs in downstream analysis.

---

## Inputs (what you choose)

### X and Y variables

In the UI you select:
- **X Axis Variable:** (continuous obs or gene expression)
- **Y Axis Variable:** (continuous obs or gene expression)

Cellucid prevents X and Y from being identical (it won’t run if they are the same variable).

### Pages

You select pages under **Compare pages:**.
Correlation is computed **separately for each selected page**.

### Optional: Color by

If you choose a categorical obs field under **Color by:**
- points are colored by category (useful to reveal substructure),
otherwise
- points are colored by page (using the page colors you set).

---

## What is computed (exact semantics)

### Paired values only

Correlation is computed on the set of cells that have:
- a finite X value **and**
- a finite Y value

Cells missing either value are excluded from that page’s calculation.

Minimum size:
- a page with fewer than 3 paired finite values **fails the run**. The whole
  correlation stops with an error naming the page and the count it found —
  `Correlation page "<pageId>" requires at least 3 paired values; received N` —
  rather than printing a coefficient for the other pages and a blank for that one.

### Pearson vs Spearman

- **Pearson** measures linear correlation on the raw values.
- **Spearman** computes Pearson correlation on the ranked values (rank correlation).

### Reported statistics

Cellucid computes (per page):
- `r` (correlation coefficient)
- `r²`
- `p-value` (two-sided Student t distribution using `n-2` degrees of freedom)
- `n` (number of paired values used)
- `slope` and `intercept` for the trend line (linear regression)

:::{note}
**Spearman's p-value is an approximation, and its ties are not corrected for.**

The coefficient itself is exact: values are converted to ranks with tied values
sharing their average rank, and Pearson's formula is applied to those ranks.
The p-value, however, comes from the *same* Student-t transform used for
Pearson — `t = r √((n−2)/(1−r²))` on `n−2` degrees of freedom — with no
tie correction and no exact permutation branch. `slope` and `intercept` are
always the least-squares fit to the **raw** values, so the trend line drawn on
the scatter plot does not change when you switch to Spearman.

That transform is a good approximation for a moderate `n` with few ties. It gets
optimistic in exactly the case single-cell data produces most often: a sparse
gene where thousands of cells share the value 0, so the ranks are dominated by
one enormous tie group. Read `r` in that situation and treat the p-value as
ordering information, not as a calibrated error rate.
:::

:::{important}
The p-value is per page and is not multiple-testing corrected.
If you run many correlations (many variable pairs, many pages), do not treat “p < 0.05” as a global discovery threshold.

For analysis-wide scope and assumptions, see {doc}`01_analysis_mental_model`.
:::

---

## Plotting vs computation (large datasets)

To keep the scatter plot responsive, Cellucid may downsample points for plotting:
- the plot shows at most ~50,000 points per page,
- downsampling is deterministic (stride-based) so repeated renders are stable.

The reported statistics (`r`, `p-value`, `n`) are computed on the **full paired dataset**, not just the plotted subset.

---

## Fast path (practical)

1) Create/select a page you want to analyze
2) Open **Analysis → Correlation**
3) Choose X and Y variables
   - Start with continuous obs fields (QC) to confirm the workflow.
4) Select one or more pages under **Compare pages**
5) Optional: set **Color by** to a categorical obs field
6) Choose Pearson vs Spearman
7) Interpret
   - inspect scatter shape,
   - confirm `n` is large enough to trust,
   - consider whether the relationship is driven by a hidden grouping (Color by helps).

---

## Interpretation pitfalls (don’t skip)

### Zero inflation (gene expression)

Many single-cell expression matrices are sparse:
- many cells have exactly 0 for a gene,
- correlation can be dominated by the “zero cloud”.

Practical tips:
- try Spearman if you suspect non-linearity,
- stratify by cell type (Color by or separate pages).

### Confounding and mixture effects

If your page contains a mixture of cell types, you can see strong correlation simply because:
- cell type A expresses both genes high,
- cell type B expresses both genes low.

That can be biologically real, but it’s not “within cell type” correlation.

### Outliers

Pearson is sensitive to outliers.
If a small number of extreme values drive the relationship, consider:
- Spearman,
- or plotting with a log axis (if valid) and inspecting density contours.

### Log scales

Correlation does not automatically log-transform.
If you enable log axes in the plot options:
- values ≤ 0 cannot be shown on a log axis,
- and the visual can differ from the underlying correlation (which is still computed on raw values).

---

## Troubleshooting (Correlation mode)

### Symptom: the run fails with `requires at least 3 paired values; received N`

The message names the page that failed and the `N` it actually found.

Likely causes:
- one or both variables are missing for most cells in that page,
- you picked a gene that is not present in this dataset,
- the page is tiny.

How to confirm:
- read `N` off the error — it is the paired finite count, not the page size,
- try a different variable pair (e.g., two QC metrics),
- check page cell counts in Compare pages.

Fix:
- deselect the offending page, or use a larger one,
- pick variables that exist and are populated for those cells,
- ensure gene expression is available for gene variables.

### Symptom: `r = 0` with `p = 1`, unexpectedly

This is not a failure or a missing value — it is the defined answer when one
variable has no variance inside that page. With a constant X or Y the
correlation denominator is zero, so `r` is set to exactly `0`, and `t` and the
p-value follow as `0` and `1`.

Likely causes:
- one variable is constant within the page (an all-zero gene is the usual one),
- after rank conversion, Spearman sees a single tie group and the same thing happens.

Fix:
- choose a gene or field with variance in *this* page, not in the dataset overall,
- split pages by a categorical field and re-run.

---

## Interface reference

```{figure} ../../../_static/screenshots/analysis/correlation-results.png
:alt: The Correlation panel with X Axis Variable set to Continuous obs then S_score, Y Axis Variable set to G2M_score, both pages selected, Color by set to cell_type, a Correlation method select reading Pearson (linear) under a pointer, and a coloured scatter plot with trend lines beneath.
:width: 488px

Correlation on the two cell-cycle scores, both pages selected and points coloured
by `cell_type`. **Color by** changes only the colouring; the coefficient is always
computed per page, never per colour.
```

```{figure} ../../../_static/screenshots/analysis/correlation-expanded.png
:alt: The expanded Correlation modal headed "CORRELATION: S_SCORE VS G2M_SCORE" with a coloured scatter and per-page statistic boxes, a PLOT OPTIONS column, a SUMMARY STATISTICS table with columns Page, R, R², P and N, and a STATISTICAL ANALYSIS table with columns Page, Direction, Strength and R.
:width: 1440px

Expanded Correlation. **SUMMARY STATISTICS** gives one row per page: `r`, `r²`,
the p-value in exponential notation, and `n` — the number of **paired finite**
values that actually entered the calculation. **STATISTICAL ANALYSIS** restates
the same `r` as a direction and a strength band, with the band boundaries printed
underneath.
```

:::{important}
Read `n` and `r` together. In the capture above, `Rest of Page 1` reports
`r = 0.4671` over `n = 3,434` with `p = 1.063e-185`. That vanishing p-value says
only that the correlation is not zero; the `r² = 0.2182` beside it says the score
explains about 22% of the variance. With tens of thousands of cells, a
scientifically uninteresting `r` of 0.02 will also carry an astronomically small
p-value. **The p-value answers “is it non-zero?”, never “is it big?”**
:::

---

## Next steps

- {doc}`06_analysis_mode_differential_expression_de` (gene-level A vs B differences)
- {doc}`04_analysis_mode_detailed_analysis` (single-variable comparisons with tests)
- {doc}`10_troubleshooting_analysis` (analysis-wide issues: missing expression, slow compute)
