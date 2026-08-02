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
- cells with the declared missing-code sentinel are excluded from grouping,
- an out-of-range or otherwise undeclared category code stops the run with an explicit data-contract error,
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

If enabled, Cellucid can reuse cached marker results for the same dataset + group-by field + settings,
and only when the group-by field still assigns exactly the same cells to the same groups.
Editing that field — merging two categories, or moving one to unassigned — changes the grouping, so the
next run is recomputed rather than read from the cache.
This can make repeated runs much faster.

### 5) Performance Settings (collapsible)

Marker discovery can be heavy. The same performance controls used by DE apply here:
- batch size, memory budget, network parallelism, and compute parallelism.

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

Wilcoxon always uses exact, tie-aware ranks. For low-cardinality values (including quantized prepared data), Cellucid aggregates exact counts per distinct value; high-cardinality values use an exact index sort. Its p-value is exact only when both comparison groups contain fewer than 50 finite values and there are no ties; otherwise inference uses the exact U statistic with a tie- and continuity-corrected asymptotic reference distribution. See {doc}`01_analysis_mental_model` for analysis-wide scope and assumptions.

Markers are filtered by thresholds (typically controlled in the expanded view):
- p-value/FDR threshold (default ~0.05)
- |log2FC| threshold (default ~1.0)
- and whether to use adjusted p-values by default

### Genes that could not be tested

The minimum group size is checked twice. The first check is on the group as a
whole and refuses the run (see
[Small groups can block the whole run](#small-groups-can-block-the-whole-run)).
The second is **per gene, per group**, on the cells that actually carry a
measured value for that gene.

A group is compared against the rest for a given gene only when **both sides
hold at least `max(2, minimum group size)` cells with a measured value** — 10 by
default, and never fewer than two, because neither the rank test nor Welch says
anything about a single observation. When either side falls short, no comparison
is run and the gene is reported as **not tested** for that group. It is not a
failure: it is a property of the data, and a run that hits it completes
normally.

An untested gene carries no p-value, so:

- **it is excluded from the Benjamini–Hochberg denominator `m`.** Each group is
  its own family — one one-vs-rest comparison across every gene in the panel —
  and `m` is the number of genes that produced a p-value *in that group*, not
  the size of the panel. Counting untested genes would deflate every adjusted
  p-value in the group: an untested gene is not a gene with p = 1;
- **it never enters the Top-N ranking**, because it has no rank to compete with;
- **it has no adjusted p-value**, and can never be counted as significant.

This is the same rule the differential-expression path applies — see
{doc}`06_analysis_mode_differential_expression_de` — so an FDR quoted from
Marker Genes and one quoted from DE mean the same thing.

:::{important}
`m` is per group, and two groups in the same run can have different denominators.
When you quote an FDR in a figure caption or a methods section, quote the
**genes tested** count shown for that group, not the size of your gene panel.
They are the same number only when every gene passed the per-gene cell check for
that group.
:::

---

## Outputs (what you see)

### Ranked Genes mode

- The **sidebar preview is a heatmap**, the same plot Clustered mode draws. The
  difference between the two modes is that Ranked Genes applies no clustering, so
  no dendrograms are computed and the Clustering settings are hidden.
- The **ranked table lives in the expanded view**: press `⤢ Expand`, then choose a
  group from the dropdown above the table. The dropdown lists each group by its
  category name — the same label the heatmap draws on its own axis, and the same
  string the ranked CSV's `group` column carries.
- Exports are in the expanded view too.

In the expanded view, a line above the marker table reports the denominator that
group's adjusted p-values were divided by:

```text
1,204 genes tested (FDR denominator) · 96 not tested (fewer than 10 cells with a measured value in this group or in the rest)
```

- The `· … not tested (…)` half appears only when that group has at least one
  untested gene.
- The number in *fewer than N cells* is the per-side requirement actually
  applied, `max(2, minimum group size)`.
- The line is hidden where there is no denominator to report: Custom Genes mode
  runs no comparison, and partial results shown while a run is still streaming
  have not corrected anything yet.

When a group has no markers to show, the message distinguishes the two reasons:

- **“No gene could be tested for this group, so it has no markers.”** — the
  denominator is zero. Nothing was compared, so no threshold could have excluded
  anything. Lowering the thresholds will not help; see
  [Genes that could not be tested](#genes-that-could-not-be-tested).
- **“No markers for selected group at the current thresholds.”** — genes *were*
  tested, and none passed the current p-value and |log2FC| cutoffs. Relaxing the
  thresholds may help.

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

**The CSV follows the mode the result was computed in**, not the plot on screen —
which matters because all three modes draw the same heatmap:

- **Ranked Genes** exports the per-group ranking behind `⤢ Expand`:
  `marker_genes_ranked.csv`, columns
  `group,gene,rank,log2FoldChange,pValue,adjustedPValue,meanInGroup,meanOutGroup,percentInGroup,percentOutGroup`
- **Clustered** and **Custom Genes** export the wide matrix the heatmap is drawn
  from: `marker_genes_heatmap.csv`, a `gene` column followed by one column per
  group.

The mode is read from the result, so the file always matches the run you are
looking at even if the form has moved on since. Group labels and gene names are
quoted where they need to be, so a comma inside a category name cannot split a
row.

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

This is the only version of the check that stops a run. A *gene* that falls
below the minimum on either side is reported as not tested and the run
continues; see [Genes that could not be tested](#genes-that-could-not-be-tested).

### No markers found

Common causes:
- thresholds too strict,
- groups are extremely similar,
- gene expression scale is inappropriate for the test (e.g., already heavily transformed),
- missing gene expression.

Read the *genes tested* line above the table before changing anything. If it
reports zero tested genes for the group, the thresholds are not the cause: no
gene in the panel had enough cells with a measured value on both sides, and the
group's own message says so rather than blaming the thresholds. A sparse matrix
or a heavily subset gene panel is the usual reason.

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

## Interface reference

```{figure} ../../../_static/screenshots/analysis/marker-genes-ranked.png
:alt: The Marker Genes panel with Mode set to Ranked Genes and the help line "Show top markers sorted by significance", a Wilcoxon method select, a Discover Markers button, and beneath it a heatmap of eight cell-type rows against many gene columns with a Z-score colour bar, its gene axis carrying evenly spaced, separated names.
:width: 488px

**Ranked Genes** after **Discover Markers**. Note that the sidebar draws a
heatmap here, not a list — the ranked table is behind `⤢ Expand`. Every group is
named down the left, but only a handful of the gene names fit across the bottom,
so expand when you need to read which gene a column is.
```

```{figure} ../../../_static/screenshots/analysis/marker-genes-clustered.png
:alt: A heatmap with cell-type groups down one axis and gene names across the other, dendrograms attached to both axes, and a colour bar for the Z-scored expression.
:width: 1160px

**Clustered** mode, shown in the expanded view. Genes and groups are each
reordered by hierarchical clustering, and the two dendrograms are drawn along
the axes. The wider box prints many more of the gene names than the sidebar can,
which is the practical reason to expand this one.
```

:::{note}
Both Ranked Genes and Clustered produce the same kind of plot. Choose **Ranked
Genes** when you only want the per-group tables and do not want to pay for
clustering; choose **Clustered** when the point is to see which genes and which
groups travel together.
:::

---

## Next steps

- {doc}`06_analysis_mode_differential_expression_de` (two-page DE with volcano plot)
- {doc}`09_exporting_analysis_results` (what each mode exports)
