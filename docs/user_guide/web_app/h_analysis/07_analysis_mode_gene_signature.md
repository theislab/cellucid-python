# Analysis mode: Gene Signature (Gene Signature Score)

**Audience:** wet lab + computational users (gene program scoring)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- What a gene signature score represents (in Cellucid terms)
- How Cellucid computes the score (mean/sum over genes; per-cell)
- How normalization options change interpretation (z-score/min-max)
- How to avoid the most common signature failures (missing genes, formatting, scale confusion)

**Prerequisites:**
- A dataset loaded
- Gene expression available
- A gene list (“signature”) to score

---

## What a gene signature is (user-facing)

A gene signature is a curated gene list meant to represent a biological state, pathway, or cell program.

Cellucid turns a signature into a **per-cell score**:
- one number per cell,
- which can be compared across highlight pages (groups).

Intuition for wet lab users:
- a high score means “many genes in the program are high (on the dataset’s expression scale)”,
- a low score means “the program is not expressed strongly in these cells”.

---

## Inputs (genes + pages)

### 1) Signature genes input

In **Analysis → Gene Signature**, you paste genes into **Signature Genes**.

Important formatting rule:
- Gene lists are parsed as **comma-separated** values (e.g., `CD3E, CD4, IL7R`).

:::{important}
Newline-only lists are not reliably parsed.

If you paste one gene per line, also include commas or convert to comma-separated format.
:::

Gene matching rules (practical):
- matching is **exact** to the dataset’s gene keys,
- no alias mapping is applied (e.g., symbol ↔ Ensembl),
- case sensitivity depends on your dataset’s gene keys (treat it as case-sensitive).

### 2) Pages (“Compare pages”)

You select pages under **Compare pages:**.
By default, if you have pages and haven’t selected anything, the UI will often start by selecting all pages.

Gene Signature supports derived pages:
- **Rest of \<page\>** for one-vs-rest comparisons.

---

## Scoring algorithm (exact)

For each selected page, Cellucid computes a score for each cell in that page.

Let `G` be the set of genes you entered.

For each cell `i`:
- collect expression values `x_{i,g}` for all genes `g ∈ G` that exist and have finite values
- compute:
  - **Mean expression**: `score_i = mean_g(x_{i,g})` over available genes
  - **Sum expression**: `score_i = sum_g(x_{i,g})` over available genes

Notes:
- Genes missing from the dataset are skipped.
- If a cell has zero valid gene values (e.g., all genes missing), its score becomes NaN and is excluded from summary statistics/plots.

:::{note}
Signature scores are computed on the expression values present in your dataset.

If your dataset stores log-normalized expression, the score is on that scale.
If your dataset stores counts, the score is on the count scale.
:::

### About the “Median expression” option

The UI exposes a “Median expression” option, but the current backend aggregation is mean/sum-based.
Until this is updated, treat “Median” as experimental and prefer **Mean expression** for reproducibility.

---

## Normalization options (what they do)

After scoring, you can normalize the scores:

- **None**: keep scores as computed
- **Z-score**: transform scores to `(x - μ) / σ` using μ and σ computed across **all selected pages combined**
- **Min-Max (0–1)**: transform scores to `(x - min) / (max - min)` using global min/max across **all selected pages combined**

Practical implication:
- normalization makes scores more comparable across pages *within the current analysis run*,
- but it also changes the meaning of “high score” (especially z-score).

---

## Outputs and interpretation

You’ll typically see:

- A distribution plot per page (violin/box/histogram depending on your selection)
- A per-page summary table (mean/median/std + number of cells with valid scores)
- A “Genes in Signature” list (chips) so you can verify inputs

Interpretation guidance:
- Compare **effect size** (how separated are distributions) before staring at p-values.
- If two pages differ strongly, check whether the signature is really a program or just a proxy for cell type composition.

Statistical tests:
- When you have ≥2 pages selected, the modal can show the same style of distribution-comparison tests used for continuous variables in Detailed mode.
Treat these as exploratory.

---

## Export (CSV)

Gene Signature exports a CSV named like `gene_signature_scores.csv` containing:
- `page`
- `score`

Important limitations:
- exported rows do not include cell indices/IDs, so the file is not directly joinable back to AnnData without additional context.

If you need per-cell mapping:
- compute the signature in Python and attach it to `adata.obs`, or
- export a per-cell table through your own workflow.

---

## Edge cases and pitfalls

- **Most genes not found**: score becomes meaningless (based on a tiny subset).
- **First gene missing**: some pages may show “No data available” behavior; put a known-present gene early in the list.
- **Duplicate genes in the list**: duplicates effectively up-weight that gene (it is added multiple times).
- **Huge signatures (hundreds of genes)**: can be slow and may stress browser memory.
- **Housekeeping-dominated signatures**: scores track library size/QC rather than biology.

---

## Troubleshooting (Gene Signature)

### Symptom: “Everything is empty / no valid values”

Likely causes:
- gene expression is not available in this dataset/loading method,
- gene keys don’t match (symbols vs Ensembl, case mismatch),
- signature input formatting is wrong (newline-only or extra punctuation).

How to confirm:
- try a single known gene (e.g., `MS4A1`) as a “signature” and see if it produces non-empty output.

Fix:
- correct gene identifiers to match your dataset,
- ensure comma-separated input,
- load the dataset with gene expression (see {doc}`../b_data_loading/index`).

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-gene-signature-success
Suggested filename: analysis/14_gene-signature-success.png
Where it appears: User Guide → Web App → Analysis → 07_analysis_mode_gene_signature.md
Capture:
  - UI location: Gene Signature mode with a non-empty signature score output
  - State prerequisites: a signature with mostly-present genes; output visible in plot/legend
  - Action to reach state: open Gene Signature → enter signature → run
Crop:
  - Include: signature input area + score output
Alt text:
  - Gene Signature mode showing a scored signature output.
Caption:
  - Gene signature scoring produces a per-cell continuous score that you can compare across groups, but interpretation depends on normalization and gene matching.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Gene Signature success state.
:width: 100%

Gene Signature mode scores a gene list and shows a per-cell signature value you can compare across groups.
```

---

## Next steps

- {doc}`08_analysis_mode_genes_panel` (marker discovery across categorical groups)
- {doc}`09_exporting_analysis_results` (export signature results)
- {doc}`10_troubleshooting_analysis` (missing genes, missing expression)
