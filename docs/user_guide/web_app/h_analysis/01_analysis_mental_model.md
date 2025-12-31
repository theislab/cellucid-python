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

Create pages in {doc}`../f_highlighting_selection/index` (the “Highlighted Cells” UI just above Analysis in the sidebar).
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
- If you want analysis to exclude cells, you must define pages accordingly (e.g., create the page while a filter is active, or maintain a dedicated “filtered group” page).

:::{tip}
Debug rule: if analysis outputs surprise you, first inspect your pages:
open the Highlighted Cells section → verify which page is active → verify group counts.
:::

### Missing values are dropped per-variable

When analysis extracts values:
- cells with missing/non-finite values for that variable are excluded from that variable’s computation,
- and different variables can end up using different effective `n` (especially in correlation).

This is why you may see:
- “Insufficient paired values” in correlation,
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
- **Bulk gene cache**: gene-heavy modes can preload and reuse gene vectors to avoid repeated I/O.
- **Marker genes cache** (Marker Genes mode): optional caching of marker results so reruns are faster.

What triggers recomputation depends on the mode, but common triggers include:
- changing selected pages or variables,
- changing page membership (adding/removing highlighted cells),
- changing plot options (thresholds, transforms, plot type),
- opening a copied (floating) analysis window (it recomputes; results are not copied).

:::{note}
Caching is a performance feature, not a scientific guarantee.

If you need a stable artifact, export it (CSV/SVG/PNG) and/or save a session bundle.
:::

---

## Reproducibility: what your numbers “mean”

Cellucid does not assume your expression values are raw counts or log-normalized.

All analysis operates on the values present in your dataset:
- If you exported log-normalized expression, DE and signatures use that scale.
- If you exported counts, they use counts.

This matters for interpretation:
- fold changes are computed on the stored scale (see DE page for exact formula),
- and correlation can differ dramatically between raw vs normalized values.

---

## Common pitfalls (read this once; it saves hours)

- **Pages overlap**: the same cell can exist in multiple pages. Comparisons are then not independent.
- **Pages are tiny**: many statistics are unstable below ~50–100 cells; DE below the app’s minimum size may be skipped.
- **Gene expression missing**: if the dataset was loaded without gene expression, gene-based modes will be unavailable/empty.
- **Gene IDs mismatch**: “TP53” vs “Tp53” vs Ensembl IDs; gene signature lists must match the dataset’s gene keys.
- **Interpreting p-values as biology**: DE and correlation tests are exploratory; batch/technical confounding can dominate.

---

## Next steps

- {doc}`02_analysis_ui_overview` (where Analysis lives; how modes and windows behave)
- {doc}`03_analysis_mode_quick_insights` (Quick mode: composition + stats at a glance)
- {doc}`10_troubleshooting_analysis` (symptom → diagnosis → fix)
