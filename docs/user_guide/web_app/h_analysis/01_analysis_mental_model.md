# Analysis mental model (Pages, Variables, Scope)

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 20–35 minutes  
**What you’ll learn:**
- The 3 core concepts: **pages**, **variables**, and **scope**
- What analysis includes/excludes (especially filters and missing values)
- What “Rest of \<page\>” means and when to use it
- What is cached vs recomputed (and why results can feel “sticky”)
- What is reproducible, and what is “UI-only”

**Prerequisites:**
- A dataset loaded in the web app
- For gene-based modes: gene expression available (var / genes)
- For most analysis: at least one **highlight page** (see below)

---

## Mental model (one sentence)

Cellucid analysis compares **highlight pages** (named groups of cells) by summarizing **variables** (obs fields or gene expression) under a clearly defined scope.

---

## The 3 nouns you must understand

### 1) Dataset = “all cells that exist”

Your dataset defines:
- the complete set of cells (indices `0..n_cells-1`),
- obs fields (categorical + continuous),
- and optionally gene expression (var / genes).

### 2) Page = “a named group of cells”

In Cellucid, the Analysis section is built around **highlight pages**:

- A **highlight page** is a named container that holds one or more **highlight groups**.
- A page’s effective membership is the **union of cell indices** across its **enabled** highlight groups.
- Pages are the “groups” you compare in most analysis modes.

:::{important}
If you do not have highlight pages, most analysis modes have nothing to compare.

Create pages in {doc}`../f_highlighting_selection/index` (the **Highlighting** panel just above **Analysis** in the sidebar).
:::

### 3) Variable = “one value per cell”

Analysis variables come from the same conceptual sources as the rest of the app:
- **Categorical obs**: labels like `cell_type`, `sample`, `cluster`.
- **Continuous obs**: numbers like `n_counts`, `pct_mito`, `score`.
- **Gene expression**: per-cell expression values for a chosen gene.

In the analysis UI you’ll often pick:
- a single variable (Detailed), or
- two variables (Correlation), or
- many genes (Gene Signature / Marker Genes).

---

## What analysis includes (and what it does not)

### Analysis operates on pages, not on “what I can currently see”

It’s common to confuse these:

- **Filtering** controls *visibility* on the canvas.
- **Pages** store *membership* (a concrete set of cell indices).

In practice:
- If you filtered cells out, they may still be in your page (because the page stores indices, not visibility).
- If you want analysis to exclude cells, define pages accordingly (e.g., create the page while a filter is active), or, in Marker genes, Differential expression, Dot plot and Composition, set **Cells** / **Within** to **Visible cells** — the one choice that follows the filters.

### Groups that are not pages

Two modes take their groups from a categorical field instead of pages:
- **Marker genes** compares every category of a **Group by** field against the
  rest or a reference category;
- **Differential expression** can compare categories (**Compare → Categories**),
  and can keep both groups inside some cells (**Within**: visible cells, a page,
  or categories of another field).

:::{tip}
Debug rule: if analysis outputs surprise you, first inspect your pages:
open the **Highlighting** panel → verify which page is active → verify group counts.
:::

### Missing values are dropped per-variable

When analysis extracts values:
- cells with missing/non-finite values for that variable are excluded from that variable’s computation,
- and different variables can end up using different effective `n` (especially in correlation).

This is why you may see:
- a correlation run stop with `requires at least 3 paired values` for a page that
  looked large enough ({doc}`05_analysis_mode_correlation_analysis`),
- or apparent “shrinking” group sizes for certain fields.

---

## “Rest of \<page\>” (derived pages)

Some analysis modes let you compare a page to its complement:

- **Rest of \<page\>** means: *all cells in the dataset that are not in \<page\>*.
- It is computed as a **set complement** relative to the dataset’s total cell count.

Typical uses:
- **One-vs-rest DE**: `Page A` vs `Rest of Page A`
- **One-vs-rest correlation/summaries**: compare a group against the background distribution

Important caveats:
- If your page is extremely large, “rest of” can be small (and vice versa).
- If pages overlap heavily, “rest of Page A” can still contain many cells from another page (because “rest of” only excludes A).

---

## Caching and recomputation (why results can be fast—or confusing)

Cellucid analysis uses caching at multiple levels:

- **Field loading cache**: obs/gene fields are loaded on demand and kept in memory for responsiveness.
- **Gene read-ahead**: a gene-heavy run reads genes ahead of its workers, as far as it measures it needs, within the analysis memory budget (see {doc}`06_analysis_mode_differential_expression_de`, How Cellucid sizes a run).
- **Marker genes cache**: the same marker run is answered from memory and marked **Cached**, with **Recompute** beside it.
- **Gene moments**: each gene's mean and spread over all cells, read once per dataset for Gene signature's Mean z-score and Module score.
- **Page data cache**: the per-cell values of each page and variable, kept in their own numeric type and bounded at 100 entries and 256 MiB, least recently used first out.
- **Density curves**: a density plot's curve is reused across cosmetic option changes and CSV export, and recomputed only when the data or the bandwidth changes.

Quick and Detailed recompute when their inputs change. Every other mode computes
only when you press **Run**; after that, a change that would change the numbers —
options, selected pages, page membership, and the filters when **Cells** /
**Within** is **Visible cells** — marks the result **Stale · Run**, even in a
mode you are not looking at, and a change that only changes what is shown
(display filters, ranking, the plot) redraws it at once. A copied (floating) analysis window starts with the settings, not the
result.

Changing the gene names or the expression matrix under **Session › Data
mapping** replaces the gene list every gene result was computed from, so every
analysis window closes and every cache above is cleared; reopen the analysis to
compute it on the new values (see {doc}`../b_data_loading/13_data_mapping`).

:::{note}
Caching is a performance feature, not a scientific guarantee.

If you need a stable artifact, export it (CSV/SVG/PNG) and/or save a session bundle.
:::

---

## Reproducibility: what your numbers “mean”

Cellucid does not assume your expression values are raw counts or log-normalized.

All analysis operates on the values of the expression matrix in use:
- If you exported log-normalized expression, DE and signatures use that scale.
- If you exported counts, they use counts.
- That is the matrix chosen under **Session › Data mapping › Genes** (`X`, a
  layer such as `log1p`, or `raw.X`). Marker genes, Differential expression and
  Dot plot can read another matrix of the file for one run with their
  **Expression matrix** option (see
  [Expression matrix](08_analysis_mode_genes_panel.md#expression-matrix)).

This matters for interpretation:
- fold changes follow **Fold change from**: by default, `log1p` is undone when the
  data mapping says the matrix holds `log1p` values, as scanpy does; otherwise
  they are on the stored scale (see {doc}`06_analysis_mode_differential_expression_de`),
- and correlation can differ dramatically between raw vs normalized values.

---

## Statistical interpretation and assumptions

Cellucid’s statistical results are exploratory comparisons of the cells represented by the selected pages.

Keep these scope limits in mind:

- Missing or non-finite values are removed according to each analysis mode, so the effective sample size can differ by variable.
- Overlapping pages are not independent because the same cell can appear in more than one group.
- Cells from the same donor or sample are not independent biological replicates. The interactive tests do not model donor, batch, condition, or other covariates; Differential expression's **Replicates** and Composition's **Samples** make samples, not cells, the observations (see {doc}`06_analysis_mode_differential_expression_de` and {doc}`12_analysis_mode_composition`).
- P-values are not globally adjusted unless a mode explicitly reports an adjustment, such as Benjamini–Hochberg values in Differential Expression or Marker Genes.

For categorical comparisons, Cellucid uses Pearson’s chi-squared test only
when every expected cell count is at least 5. It automatically switches a
sparse 2×2 table to a two-sided Fisher exact test. A larger sparse table is
reported as **N/A** because Cellucid does not implement the corresponding
general exact test; combine scientifically compatible sparse categories or
use a suitable exact workflow outside Cellucid.

For two-page continuous-distribution comparisons, Cellucid reports Welch’s t-test and Mann–Whitney U. Mann–Whitney uses an exact two-sided p-value when both groups contain fewer than 50 finite values and there are no ties; otherwise it reports a tie- and continuity-corrected asymptotic p-value. The result card shows which method was used. Multi-page Kruskal–Wallis inference corrects for tied ranks.

Interpret plots and effect sizes alongside p-values. For confirmatory inference, use a replicate-aware statistical workflow that models the experimental design.

---

## Common pitfalls (read this once; it saves hours)

- **Pages overlap**: the same cell can exist in multiple pages. Comparisons are then not independent.
- **Groups are tiny**: many statistics are unstable below ~50–100 cells. Differential expression tests a small group and flags it **low n**; Marker genes skips, pools or flags it, as **Small groups** says.
- **Gene expression missing**: if the dataset was loaded without gene expression, gene-based modes will be unavailable/empty.
- **Gene names mismatch**: a name is matched exactly, or else to the one gene that matches it ignoring case (“Tp53” → “TP53”); a name the dataset does not have is left out and reported as **not found**, never guessed.
- **Interpreting p-values as biology**: DE and correlation tests are exploratory; batch/technical confounding can dominate.

---

## Next steps

- {doc}`02_analysis_ui_overview` (where Analysis lives; how modes and windows behave)
- {doc}`03_analysis_mode_quick_insights` (Quick mode: composition + stats at a glance)
- {doc}`10_troubleshooting_analysis` (symptom → diagnosis → fix)
