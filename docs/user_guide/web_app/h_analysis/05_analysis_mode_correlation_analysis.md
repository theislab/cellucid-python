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
- if there are fewer than 3 paired values, Cellucid reports “Insufficient paired values”.

### Pearson vs Spearman

- **Pearson** measures linear correlation on the raw values.
- **Spearman** computes Pearson correlation on the ranked values (rank correlation).

### Reported statistics

Cellucid computes (per page):
- `r` (correlation coefficient)
- `r²`
- `p-value` (t-distribution approximation using `n-2` degrees of freedom)
- `n` (number of paired values used)
- `slope` and `intercept` for the trend line (linear regression)

:::{important}
The p-value is per page and is not multiple-testing corrected.
If you run many correlations (many variable pairs, many pages), do not treat “p < 0.05” as a global discovery threshold.
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

### Symptom: “Insufficient paired values”

Likely causes:
- one or both variables are missing for most cells in the selected page(s),
- you picked a gene that is all zeros or not present,
- pages are tiny.

How to confirm:
- try a different variable pair (e.g., two QC metrics),
- check page cell counts in Compare pages.

Fix:
- use larger pages,
- pick variables that exist and have variance,
- ensure gene expression is available for gene variables.

### Symptom: “Correlation is NaN or r = 0 unexpectedly”

Likely causes:
- one variable is effectively constant within the page,
- heavy missingness leaves a tiny effective `n`,
- Spearman ties (many identical values) reduce informativeness.

Fix:
- choose a more variable gene/field,
- split pages by a categorical field and re-run.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-correlation-success
Suggested filename: analysis/12_correlation-success.png
Where it appears: User Guide → Web App → Analysis → 05_analysis_mode_correlation_analysis.md
Capture:
  - UI location: Analysis → Correlation, expanded/modal view preferred
  - State prerequisites: pick X and Y variables with non-trivial variance; pick at least one page
  - Action to reach state: set X/Y → pick pages → (optional) set Color by → open expanded plot view
Crop:
  - Include: X/Y selectors + correlation plot + displayed r/p-value (if visible) + legend
Alt text:
  - Correlation mode showing a scatter plot of two variables with trend line and correlation statistics.
Caption:
  - Correlation compares X vs Y within each selected page using paired values only; plot points may be downsampled, but statistics use all paired data.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Correlation mode success state.
:width: 100%

Correlation compares X vs Y within each selected page using paired values only; plot points may be downsampled, but statistics use all paired data.
```

---

## Next steps

- {doc}`06_analysis_mode_differential_expression_de` (gene-level A vs B differences)
- {doc}`04_analysis_mode_detailed_analysis` (single-variable comparisons with tests)
- {doc}`10_troubleshooting_analysis` (analysis-wide issues: missing expression, slow compute)
