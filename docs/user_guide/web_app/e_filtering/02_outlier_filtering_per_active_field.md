# Outlier filtering (per active field)

**Audience:** everyone (wet lab + computational users)  
**Time:** 10–15 minutes  
**What you’ll learn:**
- What the Outlier filter (latent space) slider removes (and what it does *not* remove)
- How the percent threshold maps to “top tail” outliers
- What to do when the outlier slider is missing or appears to “do nothing”

**Prerequisites:**
- A dataset loaded
- A categorical obs field selected that provides latent outlier quantiles (often clusters)

---

## What this filter is (and is not)

The outlier filter is a **per-cell outlierness score** that comes from your dataset export.

It is typically computed for **categorical fields** (e.g., clustering labels) as:

- compute a centroid per category in a latent space (a vector embedding used during export),
- compute each cell’s distance to its category centroid,
- convert that distance into a within-category **quantile** in `[0, 1]`.

Interpretation:
- `0.50` means “around the median distance in this category”
- `0.95` means “in the farthest ~5% of this category”

This filter is **not**:
- a min/max range filter on a continuous field, and
- not a “global outlier” score across the entire dataset (it is category-relative when generated that way).

---

## Fast path (wet lab / non-technical)

Use this when your clusters/regions look “fuzzy” and you want to hide fringe points.

1) Pick a categorical field (clusters/cell types) in **Coloring & Filtering → Categorical obs**.
2) Look in **Display options** for **Outlier filter (latent space)**:
   - if it’s visible, your dataset provides outlier stats for this field.
3) Set the slider to `95%`.
4) Confirm in **Active filters** that you see a line like:
   - `<field>: outlier ≤ 95%`
5) If you removed too much:
   - increase the slider (e.g., `98%`), or reset to `100%`.

What success looks like:
- the “halo” / far-away cells diminish,
- the main structure becomes clearer,
- “Showing X of Y points” drops slightly (often a few percent).

---

## Practical path (computational users): exact semantics

### What the percent means

The slider is a **threshold** in `[0, 1]` displayed as a percent.

- `100%` means **no outlier filtering**.
- `p%` means: **hide cells whose outlier quantile is above `p/100`.**

Example:
- `95%` hides approximately the top 5% most outlier-like cells *within each category*, for fields where the outlier quantiles were computed that way.

### What happens with missing outlier quantiles

Outlier quantiles can be missing (e.g., `NaN`) for practical reasons:
- categories with too few cells (below a minimum threshold during export),
- or categories/cells that could not be scored.

Cells with missing outlier quantiles are **not removed by outlier filtering**.

### When the outlier slider appears

The outlier slider is shown only when the **active field** provides outlier quantiles.

Important consequences:
- If you switch to a field without outlier stats, the control disappears and outlier filtering stops applying.
- If you create a derived field by merging/deleting categories, the derived field typically does **not** carry outlier stats (so the slider may disappear).

### Data requirements (how outliers are produced)

Outlier quantiles are generated during export in `cellucid-python`:

- `cellucid.prepare(..., latent_space=<n_cells × n_dims>, obs=<DataFrame>)`
- For each categorical obs field, Cellucid computes per-cell quantiles based on latent-space distances to category centroids.

If your dataset is loaded from a source that does not provide outlier quantiles for a field, the UI simply hides the outlier control for that field.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: filtering-outlier-slider
Suggested filename: filtering/01_outlier-filter-slider.png
Where it appears: User Guide → Web App → Filtering → 02_outlier_filtering_per_active_field.md
Capture:
  - UI location: left sidebar → Coloring & Filtering (or equivalent) → outlier slider
  - State prerequisites: dataset loaded; categorical field selected that exposes the outlier slider; slider moved below 100%
  - Action to reach state: select a clustering-like categorical field and adjust the outlier slider to ~95%
Crop:
  - Include: the slider + its label + the Active filters area showing the outlier filter entry
  - Include: enough canvas to see that outliers are gone
Redact:
  - Remove: private dataset names and local paths if visible
Alt text:
  - Outlier filtering slider set below 100% with the outlier filter listed in Active filters.
Caption:
  - Show where the slider is and how to confirm it is active.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the outlier filter slider.
:width: 100%

Use the outlier filter slider to hide extreme values for the active field, and confirm it is enabled in the Active filters list.
```

---

## Troubleshooting (quick)

### “I don’t see the outlier slider”

Most common causes:
- you selected a field without outlier stats (try a categorical clusters field),
- you created/edited a derived categorical field (outlier stats often don’t carry over),
- your dataset source doesn’t include outlier quantiles for this field.

### “The slider does nothing”

Common causes:
- it’s still at `100%`,
- your categories are small and many cells have missing outlier quantiles (they won’t be filtered),
- the dataset has very tight clusters so few cells exceed `95%`.

### “It removes almost everything”

Common causes:
- threshold is too low (try `98–99%`),
- the latent space used for outlier scoring does not match the structure you’re viewing (export choice issue).

---

## Next steps

- `03_filter_stack_ui_active_filters` (how to disable/remove the outlier filter)
- `06_edge_cases_filtering` (what happens when all cells are filtered out)
