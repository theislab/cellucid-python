# Troubleshooting (analysis)

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
   - Marker Genes (Genes Panel) groups by a categorical obs field across the dataset, not by highlight pages.

If any of these fail, fix them first; most downstream symptoms disappear.

---

## Troubleshooting template (use this structure)

When adding a new troubleshooting entry, use:

### Symptom
What the user sees (include on-screen text if possible).

### Likely causes (ordered)
3–7 plausible causes, each testable.

### How to confirm
Concrete checks in the UI (what panel to open, what count to look at).

### Fix
Step-by-step actions (safe fixes first).

### Prevention
What to do earlier to avoid the issue next time.

---

## Symptom: “Analysis is empty / shows no results”

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
2) Open Filtering / Active filters
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
   - Use a loading method/export that includes expression (see {doc}`../b_data_loading/index`).
   - Re-load the dataset.

### Prevention

- Treat “page creation” as the first step of any analysis workflow.
- Name pages with intent: `T_cells__day7`, `Rest_of_T_cells`, `QC_fail`, etc.
- Before interpreting results, always confirm page size and active view.

---

## Symptom: “DE is extremely slow”

### Likely causes (ordered)

1) **Huge groups (page sizes)**
   - Runtime scales with the number of cells in Page A and Page B.
2) **Large gene count**
   - DE typically tests all genes; 20k–40k genes is common.
3) **Gene expression loading overhead**
   - If expression is loaded lazily/over the network, the data transfer dominates.
4) **Browser memory pressure**
   - Large expression batches + large plots can trigger GC churn or WebGL memory pressure.
5) **Too aggressive performance settings**
   - High parallelism + high network concurrency can saturate memory and thrash.

### How to confirm

- Check Page A / Page B sizes in the DE UI.
- If the UI shows a progress tracker, see whether it stalls on “Loading & Computing” vs “Multiple Testing Correction”.
- If you see browser slowdown across the whole app, suspect memory pressure.

### Fix

Do the least invasive thing that answers your question:

1) Reduce the scope
   - Compare smaller pages (e.g., subset to one sample/batch first).
   - If you need global DE, consider exporting and running DE offline in Python/R.
2) Use performance settings deliberately
   - Open the DE “Performance Settings” section.
   - If you see stalls/crashes, try:
     - lower network concurrency,
     - lower memory budget (to avoid spikes),
     - lower compute parallelism.
3) Avoid “repeat reruns” during setup
   - Choose pages and method first, then run once.
4) Refresh if the browser is in a degraded memory state
   - Close extra modals/windows, then reload the page and run again.

### Prevention

- For exploratory DE, start with smaller pages (hundreds–thousands of cells) to validate your hypothesis.
- Export CSV artifacts and move heavy analysis to Python/R for large studies.

---

## Symptom: “Volcano plot looks wrong / missing points”

### Likely causes (ordered)

1) **Thresholds hide most genes**
   - Volcano plots are often drawn with filters:
     - p/FDR threshold
     - |log2FC| threshold
     - “label top N” choices
2) **Adjusted p-values vs raw p-values confusion**
   - If you switch to FDR (adjusted) and your dataset has weak signal, many genes will become non-significant.
3) **Non-finite values**
   - Genes with NaN/Inf p-values or fold-changes may not be plotted or may collapse.
4) **Constant / all-zero genes**
   - These can produce undefined statistics depending on method and filtering.
5) **You compared pages in the wrong direction**
   - Page A vs Page B is directional; log2FC sign flips if you swap them.

### How to confirm

1) Confirm Page A and Page B
   - Re-check the page selector at the top of DE.
2) Open the expanded view and adjust thresholds
   - Temporarily relax thresholds (e.g., higher p/FDR threshold; lower |log2FC| threshold).
   - If points “appear”, it was a thresholding/view issue.
3) Check the gene table
   - If the gene table has rows but the plot looks sparse, it’s a plotting/threshold issue.
   - If the table is also sparse/empty, it’s likely data/signal/page definition.

### Fix

- Relax thresholds to inspect overall behavior, then tighten for reporting.
- Decide explicitly: use raw p-values for exploratory ranking; use FDR for “call significant” claims.
- If you see NaN-heavy behavior:
  - try the other method (Wilcoxon vs t-test),
  - ensure gene expression is present and not all-zero for the compared pages.

### Prevention

- Always sanity-check DE direction by verifying known markers have the expected sign.
- Save a session for each “final” DE configuration so you can reproduce the plot + thresholds later.

---

## Symptom: “Correlation is NaN”

### Likely causes (ordered)

1) **One variable is constant (no variance)**
   - Pearson/Spearman are undefined or uninformative when X or Y has no variance.
2) **Too few paired finite values**
   - Correlation is computed only on cells where both X and Y are finite.
3) **Severe missingness / sparse expression**
   - For genes, many values can be 0 or NaN depending on how expression is encoded/loaded.
4) **You expected filtering to change the analyzed set**
   - If you “filtered out” cells but your page membership still includes them, you may be correlating across the larger set.

### How to confirm

- Swap to two continuous QC obs fields with known variance (e.g., `n_counts` vs `pct_mito`).
  - If that works, the issue is your chosen variables (gene missing, constant field, etc.).
- Check the reported `n` (paired sample size) in the correlation results table.
  - If `n` is tiny (or missing), you don’t have enough paired data.

### Fix

1) Pick variables with variance and enough finite values
2) Use larger pages (more cells)
3) If expression is sparse, consider:
   - Spearman (rank-based),
   - splitting by cell type (separate pages),
   - or correlating gene vs QC covariate instead of two sparse genes.

### Prevention

- Validate that the gene exists and has non-zero expression in the relevant subset before interpreting correlation.

---

## Symptom: “Analysis windows don’t restore”

### Likely causes (ordered)

1) **You expected results to persist, but only settings persist**
   - Analysis windows restore layout and settings; results typically recompute.
2) **Dataset identity mismatch**
   - Sessions apply to a specific dataset identity; if you loaded a different dataset (or a differently-exported version), restoration may be partial.
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

1) Re-run analyses after restore
   - Treat a restored session as “settings are back; compute again”.
2) Ensure you loaded the same dataset export
   - Use the same export folder/version and dataset id.
3) If windows fail to appear entirely
   - try restoring in a fresh tab,
   - try a different browser profile (extensions can block storage).

### Prevention

- Save sessions only after you are happy with page naming and variable selection.
- Keep dataset exports versioned so “the same dataset” stays the same across collaborators.

---

## Symptom: “Marker Genes fails with ‘Group \"X\" has only N cells. Minimum required: 10.’”

### Likely causes (ordered)

1) **Your Group By field has tiny categories**
   - Common when grouping by overly granular labels, or after heavy filtering in preprocessing.
2) **You have many rare categories**
   - Even one category below the minimum can block the whole run in the current implementation.

### How to confirm

- In the Marker Genes UI, inspect the chosen Group By field:
  - do you see many groups with very small counts?
- If you have access to the data table in Python, compute per-category counts.

### Fix

1) Choose a different Group By field
   - Prefer stable, reasonably sized groups (e.g., `cell_type` over a very granular label).
2) Collapse/merge rare categories upstream
   - Merge rare labels into `Other` during preprocessing/export.
3) If you only need a subset of groups
   - consider re-exporting a dataset that excludes those rare categories.

### Prevention

- When creating categorical fields intended for marker discovery, keep a minimum group size in mind (e.g., ≥ 50 cells per category for stable marker behavior).

---

## Symptom: “Gene Signature returns nothing / looks wrong”

### Likely causes (ordered)

1) **Gene list formatting**
   - In the current UI, genes are parsed by commas (`,`). Newline-only lists can be interpreted as one “gene”.
2) **Genes not found in your dataset**
   - Gene IDs may be symbols vs Ensembl, or different casing.
3) **Gene expression missing**
   - Signature scoring requires expression.
4) **Normalization misunderstanding**
   - Z-score/min-max normalization changes the scale and can invert “intuitive” expectations.

### How to confirm

- Try a 1–2 gene signature with a known marker gene you can visualize in the field selector.
- Confirm gene expression exists by searching for that gene in other places (field selector).
- Re-enter the gene list as comma-separated: `CD3D, CD3E, LST1`

### Fix

1) Fix gene list formatting (comma-separated)
2) Use the correct gene namespace
   - If your dataset uses Ensembl IDs, provide Ensembl IDs.
3) Turn off normalization temporarily to sanity-check raw behavior
4) Reduce to a minimal signature to debug
   - then expand the list once you see non-empty results

### Prevention

- Keep a “known-good” gene list snippet for each dataset (symbols vs Ensembl) and reuse it.

---

## Symptom: “I changed filters but analysis didn’t change”

### Likely causes (ordered)

1) **You changed visibility, but analysis uses page membership**
2) **You are looking at a cached analysis result**
   - Some modes cache results until you re-run.

### How to confirm

- Check {doc}`01_analysis_mental_model` and confirm which modes require clicking Run.
- Create a new page from the currently visible cells and re-run analysis on that page.

### Fix

- If you want analysis to respect filtering, build pages after filtering (or rebuild them).
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
- try a different browser profile (extensions can block downloads).

---

## Screenshot placeholders (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-missing-expression
Suggested filename: analysis/90_missing-expression.png
Where it appears: User Guide → Web App → Analysis → 10_troubleshooting_analysis.md
Capture:
  - Show an analysis mode that requires expression with a clear empty/disabled state
  - Include any warning text (“Gene expression not available”, etc.)
Redact:
  - Remove: private dataset identifiers
Alt text:
  - Analysis mode showing an empty or disabled state due to missing gene expression.
Caption:
  - If gene expression is not available in the loaded dataset, expression-based analysis modes may be empty or disabled.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for missing gene expression in analysis.
:width: 100%

Expression-based analysis modes can be empty/disabled when gene expression was not exported or is not available in the current loading method.
```

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-no-groups-defined
Suggested filename: analysis/91_no-groups-defined.png
Where it appears: User Guide → Web App → Analysis → 10_troubleshooting_analysis.md
Capture:
  - Show the analysis panel when no selection/highlight group is defined
  - Include any “select cells”/“create group” guidance shown by the UI
Alt text:
  - Analysis panel showing an empty state because no groups are defined.
Caption:
  - Analysis requires at least one defined group (selection or highlight group); create one before running analysis.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for analysis with no groups defined.
:width: 100%

If nothing is selected and no highlight groups exist, analysis will often show an empty state until you define a group.
```

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-all-cells-filtered-out
Suggested filename: analysis/92_all-cells-filtered-out.png
Where it appears: User Guide → Web App → Analysis → 10_troubleshooting_analysis.md
Capture:
  - Show “0 visible cells” (or equivalent) and the analysis output being empty
  - Include the Active filters panel if possible
Alt text:
  - Active filters showing 0 visible cells leading to empty analysis outputs.
Caption:
  - If all cells are filtered out, analysis outputs can be empty; disable filters or reset visibility first.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for all cells filtered out affecting analysis.
:width: 100%

Analysis can appear empty when all cells are filtered out; check Active filters and restore visibility before re-running.
```

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-de-slow-warning
Suggested filename: analysis/93_de-slow-warning.png
Where it appears: User Guide → Web App → Analysis → 10_troubleshooting_analysis.md
Capture:
  - If the UI shows a progress indicator, spinner, or warning for slow DE, capture it
  - Include the group sizes if visible
Alt text:
  - Differential expression analysis showing a progress indicator or slow-compute warning.
Caption:
  - DE runtime can scale with dataset and group sizes; prefer smaller scopes or export tables for offline computation when needed.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for DE being slow.
:width: 100%

DE runtime can scale with dataset and group size; if DE is slow, reduce scope or switch to offline analysis workflows.
```

---

## Related pages

- {doc}`../e_filtering/index` (0 visible cells is often a filtering issue)
- {doc}`../f_highlighting_selection/index` (defining groups reliably)
- {doc}`../b_data_loading/index` (loading methods that include gene expression)
- {doc}`09_exporting_analysis_results` (export artifacts to debug outside the UI)
