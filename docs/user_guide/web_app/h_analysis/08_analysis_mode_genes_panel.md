# Analysis mode: Marker Genes (Genes Panel)

**Audience:** everyone (especially computational users; wet lab-friendly for interpretation)  
**Time:** 25–60 minutes  
**What you’ll learn:**
- How Marker Genes discovers one-vs-rest markers for many groups at once
- How grouping works (categorical obs field → groups)
- How to interpret log2FC, p-values/FDR, and percent-expressing
- How Ranked vs Clustered vs Custom modes differ
- How caching and performance settings affect runtime

**Prerequisites:**
- A dataset loaded
- Gene expression available
- At least one **categorical obs** field that can define groups (e.g., `cell_type` or `cluster`)

---

## What Marker Genes is for

Marker Genes is the “many groups” marker discovery tool.

Use it when you want:
- marker genes for *each* cell type/cluster in a categorical field,
- a heatmap view of those markers across groups,
- and an exportable marker table.

How it differs from DE:
- {doc}`06_analysis_mode_differential_expression_de` compares **two pages** (A vs B).
- Marker Genes computes **one-vs-rest** markers for **every group** in a categorical obs field.

---

## Inputs (what you choose)

### 1) Group By (categorical obs field)

You pick a categorical obs field under **Group By:** (e.g., `cell_type`).

Cellucid then builds groups as:
- one group per category label (e.g., `B cell`, `T cell`, …),
- across the full dataset (not highlight pages).

Important behaviors:
- cells with missing/invalid category codes are excluded from grouping,
- groups are sorted by size (largest first) for stable UI.

### 2) Mode

Marker Genes supports three modes:

- **Ranked Genes**: show ranked marker lists per group
- **Clustered**: build a marker heatmap and cluster genes/groups
- **Custom Genes**: skip marker discovery; visualize a user-supplied gene list across groups

### 3) Statistical method

- **Wilcoxon** (default): rank-based test, robust to outliers
- **t-test**: Welch’s t-test (mean comparison)

### 4) Use cached results

If enabled, Cellucid can reuse cached marker results for the same dataset + group-by field + settings.
This can make repeated runs much faster.

### 5) Performance Settings (collapsible)

Marker discovery can be heavy. The same performance controls used by DE apply here:
- batch size, memory budget, network parallelism, compute parallelism, Wilcoxon bins.

---

## Statistics (what is computed)

Marker discovery is **one-vs-rest** per group.

For each group `g` and each gene:
- “in-group” = cells in group `g`
- “out-group” = all other cells (with valid group labels)

For each gene and group, Cellucid computes (conceptually):
- `meanInGroup`, `meanOutGroup`
- `log2FoldChange = log2((meanInGroup + 0.01) / (meanOutGroup + 0.01))`
- `pValue` (Wilcoxon U or Welch t-test)
- `adjustedPValue` via Benjamini–Hochberg (computed per group across genes)
- `percentInGroup`, `percentOutGroup` = percent of cells with expression `> 0`

Markers are filtered by thresholds (typically controlled in the expanded view):
- p-value/FDR threshold (default ~0.05)
- |log2FC| threshold (default ~1.0)
- and whether to use adjusted p-values by default

---

## Outputs (what you see)

### Ranked Genes mode

- Select a group from a dropdown.
- View the top markers for that group.
- Use Expand (modal) for the full table and exports.

### Clustered mode

- A heatmap of genes (rows) vs groups (columns).
- Optional clustering of rows/columns (distance + linkage choices).
- Plot options usually include:
  - p-value threshold / log2FC threshold
  - use adjusted vs raw p-values
  - transform (e.g., z-score/log1p) and colorscale choices

### Custom Genes mode

- You provide genes.
- Cellucid builds an expression matrix for those genes across groups.
- No p-values are computed (because this mode is visualization, not discovery).

---

## Export (CSV)

Marker Genes exports depend on what you’re viewing:

- **Heatmap CSV**: `gene` column + one column per group (matrix values)
- **Ranked markers CSV**: `group,gene,rank,log2FoldChange,pValue,adjustedPValue,meanInGroup,meanOutGroup,percentInGroup,percentOutGroup`

Use exports when you need:
- reproducible reports,
- downstream filtering in R/Python,
- figure preparation outside the app.

---

## Edge cases and pitfalls

### Small groups can block the whole run

Marker discovery enforces a minimum group size (default ~10 cells).
If any group is below the minimum, the run can fail.

Workarounds:
- choose a different categorical field (coarser grouping),
- merge rare categories into “Other” in preprocessing,
- create a derived categorical field that excludes rare categories.

### No markers found

Common causes:
- thresholds too strict,
- groups are extremely similar,
- gene expression scale is inappropriate for the test (e.g., already heavily transformed),
- missing gene expression.

### Caching confusion

If “Use cached results” is enabled:
- reruns may return instantly with the same markers.
Disable caching if you suspect you changed the underlying dataset or want a fresh recompute.

---

## Troubleshooting (Marker Genes)

### Symptom: “No categorical fields available”

Cause:
- dataset has no categorical obs annotations.

Fix:
- export/load a dataset with cluster/cell type labels (obs categorical).

### Symptom: “Group ‘X’ has only N cells. Minimum required: 10.”

Cause:
- at least one category is too small for robust marker discovery.

Fix:
- merge rare categories or choose a different group-by field (see Edge cases above).

### Symptom: “Analysis is very slow / browser becomes unresponsive”

Fix:
- reduce Performance Settings (lower batch size, lower parallelism) to avoid memory pressure,
- keep fewer tabs/windows open,
- consider server mode for large datasets (data loading + memory stability).

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: analysis-genes-panel-success
Suggested filename: analysis/15_genes-panel-success.png
Where it appears: User Guide → Web App → Analysis → 08_analysis_mode_genes_panel.md
Capture:
  - UI location: Genes Panel open showing several genes and an active selection
  - State prerequisites: add 5–10 genes to the panel; click one so its expression is visible
  - Action to reach state: add genes → open panel → pick a gene
Crop:
  - Include: the Genes Panel list + gene expression legend (if it changes)
Alt text:
  - Genes Panel listing multiple genes with one active gene selected.
Caption:
  - Genes Panel speeds up iteration by keeping a curated gene list available for one-click switching.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Genes Panel success state.
:width: 100%

Marker Genes discovers one-vs-rest markers per group and visualizes them as ranked lists or a clustered heatmap.
```

---

## Next steps

- {doc}`06_analysis_mode_differential_expression_de` (two-page DE with volcano plot)
- {doc}`09_exporting_analysis_results` (what each mode exports)
