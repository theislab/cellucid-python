# Common filter types (document every filter the UI exposes)

**Audience:** computational users + power users (still readable for wet lab)  
**Time:** 20–30 minutes  
**What you’ll learn:**
- Every filter type Cellucid currently exposes in the UI
- Where each filter comes from (legend vs dedicated slider)
- Exact semantics: inclusivity, missing values, and scope

**Prerequisites:**
- A dataset loaded

---

## Quick reference (what exists today)

Cellucid’s current “visibility stack” is built from these filter types:

| Filter type | Created in the UI by… | Applies to… | Typical use |
|---|---|---|---|
| Category visibility filter | Categorical legend checkboxes | obs categorical fields | Hide clusters/samples/batches you don’t want to see |
| Numeric range filter | Continuous legend Min/Max sliders (FILTER / Live filtering) | obs continuous fields | QC gating (e.g., `n_counts`, `pct_mito`, scores) |
| Gene expression range filter | Gene legend Min/Max sliders | the active gene (var) | “Show only cells expressing gene X above Y” |
| Outlier filter (latent space) | Display options → Outlier slider (when available) | fields with outlier stats | Hide fringe/outlier-like cells within categories |

:::{important}
Not all filters apply at all times:

- Gene filters apply only while that gene is the active field.
- Outlier filtering applies only while the active field supports outlier stats (that’s when the slider is visible).
:::

---

## 1) Category visibility filter (categorical legend)

### UI name / where it lives

- Categorical legend (checkboxes next to category labels)
- Buttons: `Show All`, `Hide All`

### What it does

For a categorical obs field (clusters, sample, batch, etc.):

- a category checkbox unchecked → cells in that category become **not visible** (when category filtering is enabled)
- checked → cells in that category can be visible (subject to other filters)

### Semantics (exact)

- Category filtering is **inclusive**: “visible categories are kept”.
- The filter stack is an **AND** across all enabled filters.
- Missing/invalid category codes are typically **not removed** by hiding categories (they behave like “None/unknown”).

### “Show All / Hide All” behavior

In the categorical legend, `Show All` / `Hide All` affect **only categories that are currently available** after other filters.

Practical implication:
- if a category has `0 available` cells (because of other filters), its checkbox may be disabled and Show/Hide All won’t toggle it.

### Active filters entry

Shows up in **Active filters** when **at least one category is hidden**, even if you later disable the filter:

- `Field: hiding A, B, C +N more`

### Enable/disable vs remove

- **Disable** the category filter in Active filters:
  - cells from hidden categories can reappear,
  - but they may appear visually “de-emphasized” (e.g., greyed) because the legend still considers them hidden.
- **Remove** the category filter:
  - resets the field to “show all categories”
  - and restores default visibility state.

---

## 2) Numeric range filter (continuous legend)

### UI name / where it lives

In the continuous legend:

- “Filtering” section with Min/Max sliders
- `Live filtering` toggle
- `FILTER` and `RESET` buttons

### What it does

For a continuous obs field (QC metric/score):

- sets a visible range `[min, max]`
- cells are visible only if their value is inside the range

### Semantics (exact)

- **Inclusive range**: values `min ≤ v ≤ max` are kept.
- **Missing values**: `NaN` (and `null`) are treated as “outside range” and are filtered out whenever the filter is active.
- The filter stack is an **AND** with other filters (categories, other continuous fields, etc.).

### What “0–100” means in the sliders

The UI sliders are normalized to `0–100%` of the field’s numeric range:

- `0` corresponds to the field’s minimum value
- `100` corresponds to the field’s maximum value

The UI shows the corresponding numeric values so you don’t have to reason in percent.

### Live filtering vs FILTER button

- `Live filtering = On`: changing sliders recomputes visibility immediately.
- `Live filtering = Off`: changing sliders updates the displayed values, and you must click `FILTER` to apply.

### RESET

`RESET` returns the filter to the full range (equivalent to “no filter”).

### Active filters entry

Shows up in **Active filters** when the range is not the full range:

- `Field: min – max`

(with `(disabled)` if you turned the filter off via Active filters).

---

## 3) Gene expression range filter (var field)

Gene expression filters behave like numeric range filters, but with an important scope limitation:

- **Gene filters apply only while that gene is the active field.**

If you switch back to an obs field (clusters/QC/etc.), the gene filter no longer affects visibility.

Practical workflows:
- If you need “filter by gene, color by clusters”, consider exporting the gene score as an obs field and filtering it as obs.

---

## 4) Outlier filter (latent space)

### UI name / where it lives

When the active field provides outlier stats, Display options shows:

- `Outlier filter (latent space)` slider with a percent readout

### What it does

At `p%`, the outlier filter hides cells whose latent outlier quantile is above `p/100`.

Interpretation:
- `100%` = no outlier filtering
- `95%` ≈ removes the top 5% most outlier-like cells (often per category)

### Data requirements

Outlier quantiles are dataset-provided statistics (typically computed during export for categorical fields).

If the dataset/field does not provide outlier quantiles, the outlier control is hidden.

### Active filters entry

When active and set below 100%, Active filters shows a line like:

- `Field: outlier ≤ 95%`

### Enable/disable vs remove

- **Disable** keeps the slider value but stops it from filtering.
- **Remove** resets the outlier slider to `100%`.

See `02_outlier_filtering_per_active_field` for a deeper explanation and troubleshooting.

---

## Next steps

- `05_performance_considerations` (what gets slow and why)
- `06_edge_cases_filtering` (failure modes and surprising states)
