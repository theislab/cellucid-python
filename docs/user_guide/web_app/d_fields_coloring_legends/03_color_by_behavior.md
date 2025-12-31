# Color-by behavior

**Audience:** everyone (users interpreting colors) + computational users (exact semantics)  
**Time:** 20–30 minutes  
**What you’ll learn:**
- What “color by” actually does in Cellucid (and what it does *not* do)
- How **categorical** vs **continuous** fields map to colors
- What the neutral gray “None” color means (missing / invalid values)
- How legend interactions create **real filters** (visibility changes)
- How log scale, rescaling, and outlier filtering affect what you see

---

## First principle: coloring is data → colors, but the legend can also filter

Cellucid uses the active field to compute an RGBA color for every cell (point).

However, many legend interactions change **visibility**, not just color:
- unchecking a category hides those cells,
- narrowing a numeric range hides cells outside the range,
- enabling outlier filtering hides outliers (when supported).

This is why the UI groups it as **Coloring & Filtering**.

---

## Categorical color-by (discrete labels)

### What is colored, exactly?

For a categorical obs field, each cell has a category code:
- If the code is valid, the cell is colored with that category’s RGB color.
- If the code is missing/invalid, the cell uses a neutral gray “None” color.

Category colors come from a fixed categorical palette by default, and can be overridden per category.

### What does “hiding a category” do?

Unchecking a category in the legend is not a “dim”:
- the category becomes **fully transparent** (invisible),
- and the cells are treated as **filtered out** for downstream counts and summaries.

This means you can:
- filter by one categorical field (e.g., hide `doublets`)
- while coloring by another field (e.g., clusters)

### Category order (important for interpretation)

The legend list is sorted alphabetically (human-friendly), but:
- colors are attached to category **identity** (the category index), not the legend row order.

So “the first visible row” is not a special color; it’s just the first in alphabetical sort order.

<!-- SCREENSHOT PLACEHOLDER
ID: categorical-color-by-hiding-a-category
Suggested filename: web_app/fields_legends/20_categorical-hide-category.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 03_color_by_behavior.md
Capture:
  - Select a categorical obs field with multiple categories
  - Take the screenshot after unchecking one category so those cells disappear
  - Ensure the Active filters box shows a line like “<field>: hiding …”
Crop:
  - Include: legend with one checkbox unchecked + visible point cloud
  - Include: Active filters summary if possible
Annotations:
  - Call out: the unchecked category row and the fact that the points disappear (not just recolor)
Alt text:
  - Categorical legend with one category unchecked and the corresponding points hidden.
Caption:
  - Explain that unchecking a category is a real filter that removes those cells from view.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for hiding a categorical category.
:width: 100%

Unchecking a category in a categorical legend hides those cells (it is a real visibility filter, not just a recolor).
```

---

## Continuous color-by (numeric values)

Continuous fields (continuous obs and gene expression) map numbers → colors using:
- a **colormap** (palette), and
- a **color domain** (min/max used for color mapping).

### Neutral “None” gray (missing values)

Continuous values use the neutral gray “None” color when:
- the value is missing (`NaN`/null), or
- log scale is enabled and the value is ≤ 0 (log is undefined there).

This makes missing/invalid values obvious without forcing you to filter them out.

### Filtering vs coloring (two separate concepts)

Continuous legends expose a **visible range** (min/max sliders).

- The *filter* decides which cells are visible at all.
- The *color domain* decides how the remaining values map onto the colormap.

If you filter but keep the color domain wide, visible points can look “washed out”.
That’s why Cellucid also provides **Rescale colorbar to slider range**.

### Log scale (coloring only)

The **Use log scale** toggle affects the *color mapping*, not the filter math:

- With log scale **off**: values map linearly across the color domain.
- With log scale **on**: values map by `log10(value)` across the log-domain.

Important consequences:
- Values `<= 0` cannot be placed on a log axis, so they render as the neutral “None” gray.
- Your **Min/Max filter sliders still operate in the original units** (not log units).

If you enable log scale on a field with no positive values, the colorbar will fall back to a safe dummy range and most/all points will appear gray (because they are non-positive). That’s a data issue, not a rendering bug.

---

## Rescale colorbar to slider range (contrast control)

The toggle **Rescale colorbar to slider range** controls the color domain:

- **Off**: colors are mapped using the field’s full data range  
  (useful for consistent comparisons across different filters / views).
- **On**: colors are mapped using your current slider Min/Max range  
  (useful for bringing out contrast inside a narrow subset).

Practical example:
1) Filter a QC score to a tight range (e.g., “good cells only”).
2) Turn **Rescale colorbar to slider range** on so the remaining variation is visible.

---

## Live filtering vs FILTER (performance + predictability)

Continuous legends support two ways to apply the numeric filter:

- **Live filtering: On** (default)  
  Moving sliders immediately updates visibility.
- **Live filtering: Off**  
  Move sliders freely, then click **FILTER** to apply once.

When you might turn Live filtering off:
- very large datasets where recomputing visibility on every slider tick feels slow,
- remote desktop / low-power machines where UI responsiveness matters.

The **RESET** button restores the full numeric range for that field.

<!-- SCREENSHOT PLACEHOLDER
ID: continuous-filter-live-vs-filter-button
Suggested filename: web_app/fields_legends/21_continuous-filter-live-toggle.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 03_color_by_behavior.md
Capture:
  - Select a continuous field
  - Turn Live filtering OFF so the FILTER button becomes enabled
  - Set Min/Max sliders to a narrow range (but don’t click FILTER yet)
  - Optional: include Active filters summary to show whether the filter is already applied
Crop:
  - Include: Filtering section with Live filtering toggle + FILTER/RESET buttons + sliders
Annotations:
  - Call out: Live filtering toggle state and the enabled FILTER button
Alt text:
  - Continuous filtering controls showing Live filtering turned off and the FILTER button enabled.
Caption:
  - Explain that turning Live filtering off lets you stage slider changes and apply them in one step.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Live filtering vs FILTER behavior.
:width: 100%

Turn off Live filtering to adjust sliders without immediate recomputation, then click FILTER to apply the range once.
```

---

## How range sliders work (what “0–100” means)

The Min/Max sliders are normalized to **0–100%** of the field’s numeric range:

- `0` means the field’s minimum value
- `100` means the field’s maximum value

The UI shows the corresponding numeric value next to each slider so you don’t have to reason in percent.

If you drag Min above Max (or Max below Min), Cellucid clamps them so the range stays valid.

---

## Gene expression specifics (var fields)

Gene expression fields behave like continuous fields with a few practical differences:

- Genes are selected via search and **loaded on demand** (first selection can be slower).
- Many genes are sparse: a large fraction of cells can be exactly zero.
  - With **log scale on**, zeros become “None” gray.
  - With log scale off, zeros map to the low end of the colormap.

Current behavior to know:
- **Gene filters apply only while that gene is the active field.**  
  If you switch to an obs field, the gene range filter is not applied to visibility.

If you want “filter by gene, color by clusters” workflows, treat this as a current limitation and use:
- categorical/continuous obs filters, or
- analysis outputs exported as obs fields (then filter as obs).

---

## Outlier filter (latent space)

Some fields provide a per-cell **latent outlier quantile** (a dataset-provided statistic).
When available, the UI shows:

- **Outlier filter (latent space)** slider in Display options

Semantics:
- `100%` means **no outlier filtering**.
- Lower values remove more outliers:
  - `95%` removes cells above the 95th percentile of the outlier score (top ~5% most outlier-like).

Important details:
- Outlier filtering only applies when the current active field supports it (that’s when the control is visible).
- Derived/user-defined categorical fields created by merges/deletions typically do **not** carry latent outlier stats, so the outlier control may disappear after you “edit” a field.

<!-- SCREENSHOT PLACEHOLDER
ID: outlier-filter-slider-visible
Suggested filename: web_app/fields_legends/22_outlier-filter-slider.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 03_color_by_behavior.md
Capture:
  - Select a categorical field that shows the Outlier filter slider
  - Set it to a non-100 value (e.g., 95%) so the effect is visible
  - Include the on-screen percent readout next to the slider
Crop:
  - Include: Outlier filter slider + legend + enough of canvas to see outliers disappear
Annotations:
  - Call out: the percent setting and what it means (“remove top tail of outlier-like cells”)
Alt text:
  - Outlier filter slider set below 100% in the Display options box.
Caption:
  - Explain that the slider removes outliers when the dataset provides latent outlier quantiles for the active field.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the outlier filter slider.
:width: 100%

When available, Outlier filter (latent space) removes the most outlier-like cells; 100% keeps all cells.
```

---

## Filter stacking (why points can disappear “unexpectedly”)

Cellucid supports a filter stack across **obs fields**:
- hidden categories (from categorical legends),
- numeric ranges (from continuous legends),
- and (when active and supported) outlier filtering.

Key workflow pattern:
1) Filter using one field (e.g., continuous QC metric).
2) Switch the active field to something else (e.g., clusters) to color by it.

The filters remain active (and are summarized in **Active filters**) until you reset them.

:::{important}
Clearing the active field (choosing **None**) makes points gray, but it does **not** clear existing obs filters.
:::

---

## Quick checks when “colors look wrong”

Before assuming something is broken:

1) Look at **Coloring & Filtering**:
   - Is the active field actually the one you intended (categorical vs continuous vs gene)?
2) Look at the **legend**:
   - Is log scale on?
   - Is rescale on (which changes contrast)?
3) Look at **Active filters**:
   - Are you hiding categories on another field?
   - Did you narrow a continuous range earlier?
4) Look for the neutral “None” gray:
   - It usually means missing values, invalid codes, or log-scale-incompatible values.

For a full symptom catalog, see `d_fields_coloring_legends/05_troubleshooting_fields_legends`.

---

## Next steps

- For the exact UI controls (what each toggle/button does), go to `d_fields_coloring_legends/04_legend_behavior`.
- For selection, rename, delete, and restore semantics, go to `d_fields_coloring_legends/02_field_selector_ux`.
