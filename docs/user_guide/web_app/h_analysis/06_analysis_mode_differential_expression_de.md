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

- **Page A:** (dropdown)
- **Page B:** (dropdown)

If you have many pages, the selector shows:
- base pages under **Pages**
- derived options under **Wildcards**, including **Rest of \<page\>**

### What “Rest of \<page\>” means here

`Rest of Page A` means “all cells not in Page A”.
This enables one-vs-rest DE without manually building a background page.

### Filtering note

DE uses page membership (stored cell indices), not canvas visibility.
If you filtered cells out after creating a page, they may still be included in DE.
See {doc}`01_analysis_mental_model` for the visibility vs membership distinction.

### Minimum group size

DE requires enough cells in both groups:
- the implementation enforces a minimum of **10 valid cells** per group (default).
- if either side is below the minimum, genes may show as invalid/NaN.

---

## Statistics (exact implementation)

Cellucid supports two statistical methods:

### 1) Wilcoxon (default) = Mann–Whitney U

- This is a rank-based, non-parametric two-sample test.
- It compares the distribution of expression values between the two groups.

Implementation notes:
- For datasets with **≤ 5,000 cells total**, Cellucid computes U using an exact rank pass (sorting).
- For datasets with **> 5,000 cells**, it uses a **histogram-based approximation** of U for performance.
  - values are binned on a capped log1p scale,
  - the number of bins is controlled by **Wilcoxon bins** in Performance Settings (default 128).

Practical implication:
- on large datasets, Wilcoxon p-values are **approximate** (but designed to be stable and fast).

### 2) t-test = Welch’s t-test

- Parametric test comparing means.
- Uses Welch’s formulation (does not assume equal variances).

This can be faster than Wilcoxon on some datasets, but is less robust to heavy non-normality/outliers.

---

## Effect size: log2 fold change (log2FC)

For each gene, Cellucid computes:

`log2FC = log2((meanA + 0.01) / (meanB + 0.01))`

Where:
- `meanA` is the mean expression in **Page A**
- `meanB` is the mean expression in **Page B**
- `0.01` is a fixed pseudocount to avoid division by zero

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

---

## Outputs (what you see)

### Summary stats (quick sanity check)

You’ll see:
- Genes tested
- Significant (FDR < 0.05)
- Upregulated
- Downregulated

Upregulated/downregulated are defined by the sign of log2FC (Page A relative to Page B).

### Volcano plot

Axes:
- **x-axis:** log2 fold change (A vs B)
- **y-axis:** `-log10(p)` (raw or adjusted, depending on the toggle)

Coloring (conceptual):
- significant up (positive log2FC + passes thresholds)
- significant down (negative log2FC + passes thresholds)
- not significant (fails thresholds)

Threshold controls (in the plot options):
- p-value threshold (0.001 / 0.01 / 0.05 / 0.1)
- log2FC threshold (slider, default 1.0 = 2-fold)
- use adjusted vs raw p-values
- label top N genes
- point size and other rendering toggles

### Top genes table (sidebar convenience)

The sidebar shows a small “Top Differentially Expressed Genes” table:
- filtered by adjusted p-value threshold (default 0.05),
- sorted by absolute log2FC,
- limited to a small number of rows (for readability).

For the full ranked list, export CSV.

---

## Performance settings (when DE is slow)

DE is gene-heavy: it can touch thousands of genes and hundreds of thousands of cells.
Use **Performance Settings** (collapsible) to tune:

- **Batch size**: how many genes to preload ahead of compute (higher can be faster but uses more memory).
- **Memory budget**: limits how many gene vectors can be in flight without crashing the browser.
- **Network parallelism**: concurrent gene loads (relevant for server/remote loading).
- **Compute parallelism**: number of concurrent genes computed (bounded by worker pool size and memory).
- **Wilcoxon bins**: accuracy/speed tradeoff for approximate Wilcoxon on large datasets.

Recommended workflow for very large datasets:
- start with smaller pages (more targeted groups),
- consider t-test if Wilcoxon is too slow,
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

### Symptom: “Need at least 1 page for comparison”

Cause:
- you don’t have enough pages/wildcards to populate A and B.

Fix:
- create highlight pages in Highlighted Cells (see {doc}`../f_highlighting_selection/index`),
- or use `Rest of \<page\>` once you have at least one base page.

### Symptom: “DE is extremely slow”

Likely causes:
- huge dataset (n_cells) and/or many genes,
- Wilcoxon exact/approx costs,
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

Fix:
- lower log2FC threshold,
- switch between adjusted and raw p-values to understand correction effects,
- increase group sizes / fix page overlap.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-de-success
Suggested filename: analysis/13_de-success.png
Where it appears: User Guide → Web App → Analysis → 06_analysis_mode_differential_expression_de.md
Capture:
  - UI location: DE mode showing a volcano plot and a ranked gene table
  - State prerequisites: two reasonably sized groups with clear marker differences
  - Action to reach state: define two groups → Analysis → DE → run
Crop:
  - Include: group selectors + volcano plot + top of the gene table
Alt text:
  - Differential Expression mode showing a volcano plot and ranked genes table.
Caption:
  - DE compares two groups under the current scope; interpret results carefully when group sizes are small or confounded.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Differential Expression success state.
:width: 100%

Differential Expression shows a volcano plot and ranked genes for Group A vs Group B.
```

---

## Next steps

- {doc}`07_analysis_mode_gene_signature` (score a curated gene set)
- {doc}`09_exporting_analysis_results` (export DE tables/plots)
- {doc}`10_troubleshooting_analysis` (DE slow, volcano missing points)
