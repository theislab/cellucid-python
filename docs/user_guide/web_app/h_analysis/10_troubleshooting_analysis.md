# Troubleshooting (analysis)

**Audience:** everyone using the Analysis panel  
**Time:** find your symptom in 1–2 minutes; fixes take 5–15  
**What you’ll get:** a result where there was none, or a clear reason there cannot be one.

This page is organized by **symptom → diagnosis → fix**.

If you need help understanding what analysis should operate on in the first place, read {doc}`01_analysis_mental_model` first.

---

## Before you debug: 60‑second triage checklist

Most “analysis is broken” reports reduce to one of these.

1) **Do you have a non-empty group to analyze?**  
   - Open Highlighted Cells and confirm at least one page has a non-zero cell count.
2) **Are you on the page you think you are?**  
   - Quick mode can follow the *active* page (Dynamic).
   - Compare-style modes can have their own internal “pages to compare” selector.
3) **Are you confusing visibility with membership?**  
   - Filters change visibility; pages store membership (indices).  
   - Many analysis modes operate on page membership, not “what is currently visible”.
4) **Is gene expression actually available?**  
   - If DE/signature/marker genes look empty, confirm gene search returns genes and expression loads.
5) **Is this a mode that ignores pages?**  
   - Marker Genes (Genes Panel) groups by a categorical obs field, not by highlight pages.
6) **Is Run disabled, or the result marked Stale?**  
   - The line above Run says why it cannot run; **Stale · Run** means the form changed since the result was computed.

If any of these fail, fix them first; most downstream symptoms disappear.

---

## Symptom: “Run is disabled” or a chip on the result

Every mode with **Run** says what it needs in a muted line above the button, and
labels a result with chips for what the run decided. None of them is an error.

| You see | Meaning | What to do |
| --- | --- | --- |
| A line above Run (`Choose two pages`, `Add genes to a set`, …) | The form cannot run yet | Do what it says; each mode’s page lists every reason |
| **Stale · Run** | The options, pages or page contents changed after this result | Press it to run again, or change the form back |
| **n skipped** (Marker genes) | Groups under **Min cells per group** were not tested and stayed in the rest | Hover for their names; choose **Pool into one group** or **Test anyway, flagged** under Small groups |
| **n pooled** (Marker genes) | Small groups were tested together as one | — |
| **low n** | A group was tested with fewer cells than the minimum | Read its p-values with care |
| **n not found** (Gene signature) | Gene names the dataset does not have were left out | Hover for them; fix the names or **Remove not found** |
| **n constant** (Gene signature, Mean z-score) | Genes the same in every cell were left out | — |
| **Cached** (Marker genes) | The same run was answered from memory | **Recompute** runs it again |

---

## Symptom: “Analysis is empty / shows no results”

```{figure} ../../../_static/screenshots/analysis/analysis-empty-no-pages.png
:alt: The Quick panel reading "Page 1: 0 cells" with a Composition section stating "No composition data available." and a Statistics section stating "No statistics data available."
:width: 488px

What this looks like. The heading `Page 1: 0 cells` is the whole diagnosis: a page
exists, and nothing has been confirmed into it. Selecting cells is not enough —
press **Confirm** in the Highlighting panel. See {doc}`11_screenshots` step 1.
```

### Likely causes (ordered)

1) **No pages are defined / selected**
   - You have not created any highlight pages yet, or the analysis mode requires you to select pages and you haven’t.
2) **Your pages are empty**
   - A page exists but contains 0 cells (common after clearing selection or overwriting a page).
3) **All cells you expected are hidden (filters)**
   - Your selection might be fine, but visibility is zero (or near zero), making the plot/interaction feel “dead”.
4) **Required data is missing**
   - Expression-based modes require gene expression; correlation requires X/Y variables with finite values.
5) **You are on the wrong view/snapshot**
   - You created pages in one snapshot/view but are now looking at another.

### How to confirm

1) Open Highlighted Cells
   - Do you see at least one page with a non-zero cell count?
2) Open Filters
   - Do you have “0 visible cells” (or a very small number) due to stacked filters?
3) For DE / Gene Signature / Marker Genes
   - Try searching for a well-known gene in the gene selector.
   - If no genes appear, gene expression is likely unavailable for this dataset/loading method.
4) If using multiple snapshots/views
   - Switch to the view where you created the pages and re-check.

### Fix

Safe fixes first:

1) Create a page (group) to analyze
   - Make a selection in the viewer.
   - Save it into a highlight page (see {doc}`../f_highlighting_selection/index`).
2) If your page is empty
   - Re-fill it by saving a selection into it again.
3) If all cells are filtered out
   - Disable filters (temporarily) to confirm the diagnosis.
   - Then decide whether you want analysis on:
     - the *full page membership*, or
     - only currently visible cells (in which case: build a new page from the visible subset).
4) If gene expression is missing
   - Open a dataset that includes expression (`X`, a layer or `raw.X`; see {doc}`../b_data_loading/index`).
   - Re-load the dataset.

### Prevention

- Treat “page creation” as the first step of any analysis workflow.
- Name pages with intent: `T_cells__day7`, `Rest_of_T_cells`, `QC_fail`, etc.
- Before interpreting results, always confirm page size and active view.

---

## Symptom: “DE is extremely slow”

### Likely causes (ordered)

1) **Reading the expression matrix**
   - Almost all of a run is reading genes; over a slow or distant link that is the limit.
2) **Large gene count**
   - DE tests every gene by default; 20k–40k genes is common.
3) **Huge groups**
   - Testing grows with each gene’s non-zero values in the two groups.
4) **Browser memory pressure**
   - Other tabs or large plots leave the run less memory; it then reads less far ahead.

There are no performance settings to tune: a run measures the link, the genes and
the workers and sizes itself (see
{doc}`06_analysis_mode_differential_expression_de`, How Cellucid sizes a run).

### How to confirm

- Read the running box under Run: whether it is reading (how much is read and how fast) or `Applying multiple testing correction`. A run you do not want to wait for can be stopped with **Cancel**, and run again later.
- If the whole app slows down, suspect memory pressure.

### Fix

1) Test fewer genes
   - **More options → Filters → Only genes in**: highly variable genes, or a gene list.
2) Test fewer cells
   - **Within** a cell type or a page, or smaller pages.
3) Refresh if the browser is in a degraded memory state
   - Close extra windows and tabs, then reload and run again.

### Prevention

- For exploratory DE, start with targeted groups (hundreds–thousands of cells).
- Export CSV artifacts and move heavy analysis to Python/R for large studies.

---

## Symptom: “Volcano plot looks wrong / missing points”

Before anything else, check the thresholds. The same completed run produces both
of these, and neither is a bug:

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-a.png
:alt: A volcano plot at a log2 fold-change threshold of 0.5, with broad red and blue clouds of coloured points either side of a narrow grey band, and a legend below the axis title reading Up (750), Down (841) and Not significant (2162).
:width: 1156px

`|log₂FC| ≥ 0.5` — 1,591 genes called.
```
:::

:::{grid-item}
```{figure} ../../../_static/screenshots/analysis/volcano-threshold-b.png
:alt: The same volcano plot at a log2 fold-change threshold of 3, with the dashed threshold lines pushed far apart and nearly every point grey, and a legend below the axis title reading Up (0), Down (192) and Not significant (3561).
:width: 1156px

`|log₂FC| ≥ 3` — 192 genes called, from the identical run.
```
:::
::::

### Likely causes (ordered)

1) **Filters hide most genes**
   - The volcano colours and lists the genes passing **Max adjusted p**,
     **Min |log2FC|** and **Min % expressing**; **Labelled genes** names the top ones.
2) **Weak signal after correction**
   - The volcano always uses adjusted p-values; with weak signal few genes pass.
3) **Non-finite values**
   - Genes with NaN/Inf p-values or fold-changes may not be plotted or may collapse.
4) **Constant / all-zero genes**
   - These can produce undefined statistics depending on method and filtering.
5) **You compared pages in the wrong direction**
   - Page A vs Page B is directional; log2FC sign flips if you swap them.

### How to confirm

1) Confirm Group A and Group B
   - Re-check them at the top of the form; **Swap A and B** flips the direction.
2) Relax the filters
   - Temporarily raise **Max adjusted p** and lower **Min |log2FC|** (they redraw without testing again).
   - If points “appear”, it was a thresholding/view issue.
3) Check the gene table
   - If the gene table has rows but the plot looks sparse, it’s a plotting/threshold issue.
   - If the table is also sparse/empty, it’s likely data/signal/page definition.

### Fix

- Relax the filters to inspect overall behavior, then tighten them for reporting.
- Use the adjusted p-values for “call significant” claims; the CSV also has the raw ones for ranking.
- If you see NaN-heavy behavior:
  - try the other method (Wilcoxon vs t-test),
  - ensure gene expression is present and not all-zero for the compared pages.

### Prevention

- Always sanity-check DE direction by verifying known markers have the expected sign.
- Save a session for each “final” DE configuration so you can reproduce the plot + thresholds later.

---

## Symptom: “Correlation reports r = 0, or refuses to run”

Correlation has exactly two disappointing outcomes, and they mean different
things. It does not produce `NaN` from the UI.

### Likely causes (ordered)

1) **The run failed with `requires at least 3 paired values`**
   - Correlation counts only cells where both X and Y are finite. A page that
     drops below 3 such cells stops the whole run, and the error names the page
     and the count it found.
2) **A clean `r = 0` with `p = 1`**
   - One variable is constant inside that page (an all-zero gene is the usual
     case). The correlation denominator is zero, so `r` is defined as exactly
     `0` and the p-value follows as `1`. This is an answer, not a failure.
3) **Severe missingness / sparse expression**
   - For genes, many values can be 0 or NaN depending on how expression is encoded/loaded.
4) **You expected filtering to change the analyzed set**
   - If you “filtered out” cells but your page membership still includes them, you may be correlating across the larger set.

### How to confirm

- Read the error text, or read `n` in the correlation results table. `n` is the
  paired finite count, not the page size, and the gap between them is the
  diagnosis.
- Swap to two continuous QC obs fields with known variance (e.g., `n_counts` vs `pct_mito`).
  If that works, the issue is your chosen variables (gene missing, constant field, etc.).

### Fix

1) Pick variables that have variance **inside the page you selected**, not just across the dataset
2) Use larger pages (more cells), or deselect the page named in the error
3) If expression is sparse, consider:
   - Spearman (rank-based),
   - splitting by cell type (separate pages),
   - or correlating gene vs QC covariate instead of two sparse genes.

### Prevention

- Validate that the gene exists and has non-zero expression in the relevant subset before interpreting correlation.
- See {doc}`05_analysis_mode_correlation_analysis` for the exact pairing rule and
  the caveat on Spearman's p-value.

---

## Symptom: “Analysis windows don’t restore”

### Likely causes (ordered)

1) **You expected results to persist, but only settings persist**
   - Analysis windows restore layout and settings; results typically recompute.
2) **Dataset identity mismatch**
   - Sessions apply to a specific dataset identity; if you loaded a different dataset (or a differently-written version), restoration may be partial.
3) **Local browser storage cleared / blocked**
   - Some persistence relies on browser storage; private browsing / strict policies can interfere.
4) **You’re in a different view/snapshot**
   - Windows may restore, but the “active view” context differs from what you expect.

### How to confirm

- Open the session restore page and confirm the session load completed without errors.
- After restore:
  - confirm windows exist (even if plots are blank),
  - check whether page selections/variables match your saved configuration.

### Fix

1) Check what a session restores
   - Settings come back for every analysis; differential expression, marker
     genes, signatures and correlation also bring back their results. Quick
     and Detailed views compute again from their settings.
2) Ensure you loaded the same dataset
   - Use the same dataset file (version) and dataset id.
3) If windows fail to appear entirely
   - try restoring in a fresh tab,
   - verify in a clean profile with extensions disabled (extensions can block storage).

### Prevention

- Save sessions only after you are happy with page naming and variable selection.
- Keep dataset files versioned so “the same dataset” stays the same across collaborators.

---

## Symptom: “Marker Genes skipped a group”

### Likely causes (ordered)

1) **The group is smaller than Min cells per group** (10 by default)
   - Under the default **Small groups → Skip and list**, it is not tested, stays in the rest, and the result shows **n skipped**.
2) **Cells narrows the groups**
   - With **Cells** set to visible cells or a page, groups are counted inside those cells only.

### How to confirm

- Hover the **n skipped** chip: it names each group and its cell count.
- Open **Groups to test**: every group shows its cell count, and small ones are flagged.

### Fix

1) **More options → Groups and cells → Small groups**: **Pool into one group**, or **Test anyway, flagged** (any group with at least 2 cells).
2) Lower **Min cells per group** (at least 2).
3) Choose a coarser **Group by** field.

### Prevention

- Fields meant for marker discovery work best with ≥ 50 cells per category; read **low n** results with care.

---

## Symptom: “Gene Signature returns nothing / looks wrong”

### Likely causes (ordered)

1) **Genes not found in your dataset**
   - The names you pasted may be a different vocabulary from the dataset’s (symbols vs Ensembl ids, human vs mouse names). The result says **n not found**; Run is disabled when no gene of any set is found.
2) **Gene expression missing**
   - Signature scoring requires expression.
3) **Rescaling misunderstood**
   - z-score and 0–1 rescale every drawn score together and can change “intuitive” expectations; the CSV keeps the scores as computed.
4) **Method mismatch**
   - Mean scores track expression level; Module score and UCell are relative (see {doc}`07_analysis_mode_gene_signature`).

### How to confirm

- Hover **n not found** for the names left out.
- Try a 1–2 gene set with a known marker gene you can find in the field selector.

### Fix

1) Use the names the dataset publishes (read a few in the field selector), then **Remove not found**.
2) Set **Rescale** to None to sanity-check raw behavior.
3) Reduce to a minimal set to debug, then expand it.

### Prevention

- Keep a “known-good” gene list for each dataset, in that dataset’s own names, and import it.

---

## Symptom: “I changed filters but analysis didn’t change”

### Likely causes (ordered)

1) **You changed visibility, but analysis uses page membership**
2) **The result is from before the change**
   - Modes with **Run** keep their result until you run again; a change to pages or page contents marks it **Stale · Run**.

### How to confirm

- Check {doc}`01_analysis_mental_model` and confirm which modes require clicking Run.
- Create a new page from the currently visible cells and re-run analysis on that page.

### Fix

- If you want analysis to respect filtering, build pages after filtering (or rebuild them), or set **Cells** / **Within** to **Visible cells** in Marker genes and Differential expression.
- Re-run the analysis after changing the group definition.

### Prevention

- Adopt a convention: “filters are for looking; pages are for computing”.

---

## Symptom: “Exported CSV/PNG is missing / downloads don’t start”

Likely causes:
- your browser blocks downloads,
- you didn’t open the expanded view (where export controls live),
- you don’t have a computed result yet.

Fix:
- follow {doc}`09_exporting_analysis_results`,
- verify downloads in a clean profile with extensions disabled.

---

## Related pages

- {doc}`../e_filtering/index` (0 visible cells is often a filtering issue)
- {doc}`../f_highlighting_selection/index` (defining groups reliably)
- {doc}`../b_data_loading/index` (loading methods that include gene expression)
- {doc}`09_exporting_analysis_results` (export artifacts to debug outside the UI)
