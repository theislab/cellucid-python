# Common filter types

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
| Numeric range filter | Continuous legend Min/Max sliders (Filter / Live filtering) | obs continuous fields | QC gating (e.g., `n_counts`, `pct_mito`, scores) |
| Gene expression range filter | Gene legend Min/Max sliders | the active gene (var) | “Show only cells expressing gene X above Y” |
| Outlier filter | Legend → **Outlier filter** slider | categorical obs fields of the file, once a latent space is chosen | Hide fringe/outlier-like cells within categories |

:::{important}
Not all filters apply at all times:

- Gene filters apply only while that gene is the active field.
- Outlier filtering applies only while the categorical field it was set on is the active field.
:::

---

## 1) Category visibility filter (categorical legend)

### UI name / where it lives

- Categorical legend (checkboxes next to category labels)
- Buttons: `Show all`, `Hide all`

### What it does

For a categorical obs field (clusters, sample, batch, etc.):

- a category checkbox unchecked → cells in that category become **not visible** (when category filtering is enabled)
- checked → cells in that category can be visible (subject to other filters)

### Semantics (exact)

- Category filtering is **inclusive**: “visible categories are kept”.
- The filter stack is an **AND** across all enabled filters.
- Missing/invalid category codes are typically **not removed** by hiding categories (they behave like “None/unknown”).

```{figure} ../../../_static/screenshots/filtering/type-category-visibility.png
:alt: A legend headed by the hint Click a swatch to pick a color, then Show All and Hide All buttons, then eight rows; Ductal, Ngn3 high EP and Ngn3 low EP have empty checkboxes, grey swatches and grey labels and read 0 of 916, 0 of 642 and 0 of 262 cells, while Alpha, Beta, Delta, Epsilon and Pre-endocrine stay ticked and coloured with a plain cell count.
:width: 428px

Three categories hidden. Note what a hidden row keeps: its swatch position, its
label, and its total (`0 / 916 cells`). Nothing was deleted.
```

### “Show all / Hide all” behavior

In the categorical legend, `Show all` / `Hide all` affect **only categories that are currently available** after other filters.

Practical implication:
- if a category has `0 available` cells (because of other filters), its checkbox is disabled and Show/Hide All skip it entirely.

### Active filters entry

Shows up in **Filters** when **at least one category is hidden**, even if you later disable the filter:

- `Field: hiding A, B, C +N more`

### Enable/disable vs remove

- **Disable** the category filter in Filters:
  - cells from hidden categories can reappear,
  - but they may appear visually “de-emphasized” (e.g., greyed) because the legend still considers them hidden.
- **Remove** the category filter:
  - resets the field to “show all categories”
  - and restores default visibility state.

---

## 2) Numeric range filter (continuous legend)

### UI name / where it lives

In the continuous legend:

- the **Visible range** section with `Minimum` / `Maximum` sliders
- `Live filtering` toggle
- `Filter` and `Reset` buttons

### What it does

For a continuous obs field (QC metric/score):

- sets a visible range `[min, max]`
- cells are visible only if their value is inside the range

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: A block headed Filtering with the subtitle Visible range, a help line reading Adjust limits then click FILTER or enable Live filtering, a Live filtering toggle reading Off with the hint Drag sliders then click Filter, a Min slider at the far left reading minus 0.37 and a Max slider at the far right reading 1.14, and black FILTER and RESET buttons with the mouse pointer on FILTER.
:width: 428px

The complete numeric-range control, with `Live filtering` switched **off** so
`Filter` is clickable.
```

### Semantics (exact)

- **Inclusive range**: values `min ≤ v ≤ max` are kept.
- **Missing values**: `NaN` is treated as “outside range” and is filtered out whenever the filter is active. A `NaN` in a continuous obs column is a missing value (both writers store one as `NaN`); gene expression written by `prepare()` or `cellucid_prepare()` never holds one.
- The filter stack is an **AND** with other filters (categories, other continuous fields, etc.).

### What “0–100” means in the sliders

The UI sliders are normalized to `0–100%` of the field’s numeric range in whole
steps of 1:

- `0` corresponds to the field’s minimum value
- `100` corresponds to the field’s maximum value

The UI shows the corresponding numeric values so you don’t have to reason in percent.

The two sliders cannot cross: dragging `Minimum` past `Maximum` pushes `Maximum` along with
it, and vice versa.

### Live filtering vs Filter button

`Live filtering` starts **On**.

- `Live filtering = On`: changing sliders recomputes visibility immediately —
  and `Filter` is **greyed out and unclickable**, with the tooltip
  `Turn off Live filtering to use Filter`.
- `Live filtering = Off`: changing sliders updates the displayed values only,
  and you must click `Filter` to apply.

```{figure} ../../../_static/screenshots/filtering/filter-button-disabled-live-on.png
:alt: The same Filtering block with the Live filtering toggle reading On, the hint below it reading Changes apply as you drag, and a pale grey FILTER button under the mouse pointer beside a black RESET button.
:width: 428px

With `Live filtering` on, `Filter` is disabled. This is the intended state, not
a fault: there is nothing for the button to do.
```

:::{warning}
**The toggle does not stay off.** It is rebuilt with the legend, so it returns
to `On` whenever the legend re-renders: a different active field, a different
colormap, a different kept view, or a loaded session. It is also not stored in
a `.cellucid-session`. Check it again after any of those before you drag on a
large dataset — see {doc}`05_performance_considerations`.
:::

### Reset

`Reset` returns the filter to the full range (equivalent to “no filter”). It
also snaps both sliders back to `0` and `100`.

### Active filters entry

Shows up in **Filters** when the range is not the full range:

- `Field: min – max`, both numbers rounded to two decimal places

(with `(disabled)` appended if you turned the filter off via Filters —
though on a narrow sidebar that suffix is usually past the ellipsis, so read
the greyed, struck-through row instead).

```{figure} ../../../_static/screenshots/filtering/window-continuous-range-applied.png
:alt: The whole application window coloured by S_score on a Viridis scale with the Min slider moved to the right, only a scattered subset of the cells drawn on the grid, and the sidebar reporting a reduced visible count above a single S_score range row.
:width: 1440px

An `S_score` range filter in effect. Coloring and filtering are driven by the
same field here, which is the common case — but they are separate: the colour
scale still spans the full data range unless you also turn off **Rescale
colorbar to slider range**.
```

---

## 3) Gene expression range filter (var field)

Gene expression filters behave like numeric range filters — same block, same
`Live filtering` / `Filter` / `Reset` controls — but with an important scope
limitation:

- **Gene filters apply only while that gene is the active field.**

If you switch back to an obs field (clusters/QC/etc.), the gene filter no longer affects visibility.

```{figure} ../../../_static/screenshots/filtering/window-gene-range-active.png
:alt: The whole application window coloured by the gene Ins1, the Gene Expression box holding Ins1, the Min slider raised, only a small group of surviving points on an otherwise empty grid, and the sidebar listing a single Ins1 range row.
:width: 1440px

**While `Ins1` is the active field**, its range filter is listed and enforced.
```

```{figure} ../../../_static/screenshots/filtering/window-gene-range-scope-lost.png
:alt: The same window one action later with the Categorical obs select holding clusters and Gene Expression empty; the full coloured embedding is back, the filter list reads No filters active and the line above reads Showing all 3,696 points.
:width: 1440px

**After switching the active field to `clusters`**, the gene row has left the
list and every cell is drawn again. Nothing was removed by you — the filter is
simply out of scope.
```

Practical workflows:
- If you need “filter by gene, color by clusters”, consider exporting the gene score as an obs field and filtering it as obs.

---

## 4) Outlier filter

### UI name / where it lives

When the active field is a categorical obs field of the file, the **Legend**
box shows, below **Show centroids**:

- an `Outlier filter` slider with a percent readout,
- a line under it, when there is something to say: that a latent space must
  be chosen first (with a **Choose latent space** button), that the field's
  categories are too small to measure, what the first use reads, or the
  progress of that read (with **Cancel**).

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: A control labelled Outlier filter (latent space) whose slider handle sits a short distance in from the right end with the number 95 percent printed beside it and the mouse pointer on the handle.
:width: 436px

The slider runs `0`–`100` in whole steps and starts at `100`.
```

### What it does

At `p%`, the outlier filter hides cells whose latent outlier quantile is above `p/100`.

Interpretation:
- `100%` = no outlier filtering
- `95%` ≈ removes the top 5% most outlier-like cells (often per category)

### Data requirements

The filter measures in a latent space (an `obsm` array such as `X_scVI` or
`X_pca`), named by the file (`prepare(latent_key=...)`,
`cellucid_prepare(latent_key = ...)`), by the launch (`--latent-key`,
`latent_key=`), or by you under **Data mapping › Latent space**. Without one
the control is shown disabled
([Turning it on](02_outlier_filtering_per_active_field.md#turning-it-on-choose-a-latent-space)).

Its quantiles are either stored in the file for that latent space (read at
once, one field at a time) or computed in the browser on the slider's first
move below `100%`, with progress and **Cancel**
([Stored or computed](02_outlier_filtering_per_active_field.md#stored-or-computed)).

A field you create by merging or deleting categories has no outlier
quantiles, so the control is hidden for it.

### Active filters entry

When active and set below 100%, Filters shows a line like:

- `Field: outlier ≤ 95%`

### Enable/disable vs remove

- **Disable** keeps the slider value but stops it from filtering.
- **Remove** resets the outlier slider to `100%`.

See {doc}`02_outlier_filtering_per_active_field` for a deeper explanation and troubleshooting.

---

## Next steps

- {doc}`05_performance_considerations` (what gets slow and why)
- {doc}`06_edge_cases_filtering` (failure modes and surprising states)
