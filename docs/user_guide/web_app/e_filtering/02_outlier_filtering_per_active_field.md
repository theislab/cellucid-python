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
3) Drag the slider from `100%` down to `95%`.

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: A control labelled Outlier filter (latent space) whose slider handle sits a short distance in from the right end, with the number 95 percent printed to its right and the mouse pointer resting on the handle.
:width: 436px

The control lives inside the dashed **Display options** box, directly above
**Show centroids**. It starts at `100%`, which means no outlier filtering at
all.
```

4) Confirm in **Active filters** that you see a line like:
   - `<field>: outlier ≤ 95%`
5) If you removed too much:
   - increase the slider (e.g., `98%`), or reset to `100%`.

### What success looks like

```{figure} ../../../_static/screenshots/filtering/window-outlier-100-off.png
:alt: The whole application window coloured by clusters with the outlier slider at its right end reading 100 percent, loose single points scattered in the space around and between the coloured groups, and the sidebar reading Showing all 3,696 points above No filters active.
:width: 1440px

**At `100%`** — no filtering. Note the loose points drifting between the
coloured groups, and the legend counts: Alpha 481, Ductal 916, Pre-endocrine
592.
```

```{figure} ../../../_static/screenshots/filtering/window-outlier-95-applied.png
:alt: The same window with the slider moved left and reading 95 percent, the mouse pointer on it, visibly fewer loose points between the coloured groups, every legend count reduced to Alpha 458, Ductal 872 and Pre-endocrine 564, and the sidebar reading Showing 3,518 of 3,696 points above a row reading clusters colon outlier less than or equal to 95 percent.
:width: 1440px

**At `95%`** — the same camera, the same colouring, 178 fewer cells drawn.
Every category has lost roughly its own top 5%: Alpha 481 → 458, Ductal
916 → 872, Pre-endocrine 592 → 564. That per-category proportionality is the
signature of this filter, and it is what distinguishes it from a global
threshold.
```

So: the “halo” of far-away cells diminishes, the main structure becomes
clearer, and “Showing X of Y points” drops by a few percent — not by half.

---

## Practical path (computational users): exact semantics

### What the percent means

The slider is a **threshold** in `[0, 1]` displayed as a percent.

- `100%` means **no outlier filtering**.
- `p%` means: **hide cells whose outlier quantile is above `p/100`.**

Example:
- `95%` hides approximately the top 5% most outlier-like cells *within each category*, for fields where the outlier quantiles were computed that way.

### When a row appears in Active filters

The outlier row is written only when the threshold is **below `99.99%`**. A
slider parked at `100%` produces no row, which is correct — it is filtering
nothing — but it is also why “the slider does nothing” is such a common
report.

The row also carries **no `visible / available` cell count**, unlike category
and range rows. Its only numeric effect is on the `Showing X of Y points` line.

### What happens with missing outlier quantiles

A quantile can be missing, and a cell with a missing quantile is **never
removed by outlier filtering**.

The export writes a missing quantile as `NaN` (or, in the quantized encoding,
as the terminal integer marker: `255` for 8-bit, `65535` for 16-bit). The web
app converts that to exactly `-1` when it loads the field, so its visibility
loop can tell “unavailable” apart from the valid `[0, 1]` interval without
carrying a non-finite value into a hot loop.

Why a quantile goes missing:

- the cell's category holds fewer than `centroid_min_points` cells (**default
  `10`**), so no reliable centroid could be computed for it;
- or the category could not be scored at all.

:::{note}
An export in which *every* quantile for a field is missing is rejected outright
at preparation time, with an error naming `centroid_min_points`. So a shipped
dataset either has usable outlier statistics for a field or has none at all —
never an all-missing set that would silently do nothing.
:::

### When the outlier slider appears

The outlier slider is shown only when the **active field** provides outlier quantiles.

Important consequences:
- If you switch to a field without outlier stats, the control disappears and outlier filtering stops applying.
- The control never appears for a **continuous** field or a **gene**, whatever the dataset — outlier quantiles are a per-category statistic.
- If you create a derived field by merging/deleting categories, the derived field typically does **not** carry outlier stats (so the slider may disappear).

### Data requirements (how outliers are produced)

Outlier quantiles are generated during export in `cellucid-python`:

- `cellucid.prepare(..., latent_space=<n_cells × n_dims>, obs=<DataFrame>)`
- For each categorical obs field, Cellucid computes per-cell quantiles based on latent-space distances to category centroids.

If your dataset is loaded from a source that does not provide outlier quantiles for a field, the UI simply hides the outlier control for that field.

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

- {doc}`03_filter_stack_ui_active_filters` (how to disable/remove the outlier filter)
- {doc}`06_edge_cases_filtering` (what happens when all cells are filtered out)
