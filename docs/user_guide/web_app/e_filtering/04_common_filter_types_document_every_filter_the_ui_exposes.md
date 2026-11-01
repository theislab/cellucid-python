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
Every filter applies until it is removed, whatever the view is coloured by:
a gene's range filter and a categorical field's outlier filter keep hiding
cells after you colour by another field, and **Filters** lists each of them.
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
:alt: The categorical legend for clusters: Show all and Hide all above eight rows, three of them (Ductal, Ngn3 high EP, Ngn3 low EP) unchecked and reading 0 of their cell counts, each row carrying a colour swatch, a label, a cell count, a rename pencil and a delete bin.
:width: 488px

Three categories hidden. Note what a hidden row keeps: its swatch position, its
label, and its total (`0 / 916 cells`). Nothing was deleted.
```

### “Show all / Hide all” behavior

In the categorical legend, `Show all` / `Hide all` affect **only categories that are currently available** after other filters.

Practical implication:
- if a category has `0 available` cells (because of other filters), its checkbox is disabled and **Show all** / **Hide all** skip it entirely.

### Filters list entry

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

- the **Visible range** section with `Minimum` / `Maximum`, each a slider
  and a number box
- `Live filtering` toggle
- `Filter` and `Reset` buttons

### What it does

For a continuous obs field (QC metric/score):

- sets a visible range `[min, max]`
- cells are visible only if their value is inside the range

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: The Visible range block of the continuous S_score legend with Live filtering switched off: Minimum and Maximum sliders with their values, and the Filter button beside Reset enabled, with the pointer on Filter.
:width: 488px

The complete numeric-range control, with `Live filtering` switched **off** so
`Filter` is clickable.
```

### Semantics (exact)

- **Inclusive range**: values `min ≤ v ≤ max` are kept.
- **Missing values**: `NaN` is treated as “outside range” and is filtered out whenever the filter is active. A `NaN` in a continuous obs column is a missing value (both writers store one as `NaN`); gene expression written by `prepare()` or `cellucid_prepare()` never holds one.
- The filter stack is an **AND** with other filters (categories, other continuous fields, etc.).

### Sliders and number boxes

Each bound has a slider and a number box beside it, and the box always holds
the exact bound applied:

- **Type a bound** in the box (a decimal point or a decimal comma; `1e-3`
  works too) and press Enter or leave the box: it is applied exactly as
  typed, so `≥ 1500` on `total_counts` from 500 to 100,000 is one entry.
  A value past the data's minimum or maximum is that end; text that is not
  a number is not applied, and the box shows the bound again.
- **Drag a slider** for the feel: it moves in thousandths of the track, and
  the bound it sets is the shortest decimal that tells one step from the
  next (`1500`, `50250`, `0.0485`), which the box then shows. The track's
  ends are the field's minimum and maximum exactly. With **Log color scale**
  on, the track is spaced in log units from the smallest positive value to
  the maximum, with the zeros (and anything below the smallest positive
  value) at its first step, so most of the track is no longer the top decade.
- Only the bound you move or type changes: the other keeps its exact value,
  also one restored from a session. The two cannot cross: setting `Minimum`
  past `Maximum` pushes `Maximum` to it, and vice versa.

### Live filtering vs Filter button

`Live filtering` starts **On**.

- `Live filtering = On`: changing sliders recomputes visibility immediately,
  so there is no `Filter` button; only `Reset` sits under the sliders.
- `Live filtering = Off`: changing sliders updates the displayed values only;
  `Filter` appears beside `Reset`, and you click it to apply.

:::{warning}
**The toggle does not stay off.** It is rebuilt with the legend, so it returns
to `On` whenever the legend re-renders: a different active field, a different
colormap, a different kept view, or a loaded session. It is also not stored in
a `.cellucid-session`. Check it again after any of those before you drag on a
large dataset — see {doc}`05_performance_considerations`.
:::

### Reset

`Reset` returns the filter to the full range (equivalent to “no filter”). It
also returns both sliders to the track's ends and the boxes to the field's
minimum and maximum.

### Filters list entry

Shows up in **Filters** when the range is not the full range:

- `Field: min – max`, each bound to six significant figures (scientific
  notation below 0.01 and from 1,000,000), so a typed or slider bound reads
  exactly as applied; figure provenance writes it the same way

(with `(disabled)` appended if you turned the filter off via Filters —
though on a narrow sidebar that suffix is usually past the ellipsis, so read
the greyed, struck-through row instead).

```{figure} ../../../_static/screenshots/filtering/window-continuous-range-applied.png
:alt: The whole window coloured by S_score with the Minimum slider raised: only the higher-scoring cells are drawn, and the Filters block reads a reduced count above one S_score range row.
:width: 1440px

An `S_score` range filter in effect. Coloring and filtering are driven by the
same field here, which is the common case — but they are separate: with
**Rescale to range** on (the default) the colour scale follows the slider
range; turn it off to keep the colours spanning the full data range.
```

---

## 3) Gene expression range filter (var field)

Gene expression filters behave like numeric range filters — same block, same
`Live filtering` / `Filter` / `Reset` controls — and like them they keep
applying when you colour by another field: "filter by gene, colour by
cluster" works. The filter is listed under **Filters** as
`<gene>: <min> – <max>`; remove or disable it there, or colour by the gene
again to change its range.

```{figure} ../../../_static/screenshots/filtering/window-gene-range-active.png
:alt: The whole window coloured by the gene Ins1 with the Minimum slider raised, so only the cells expressing it are drawn and an Ins1 range row is listed in the Filters block.
:width: 1440px

Coloured by `Ins1`, its range filter is listed and hides every cell outside
the range; coloured by another field afterwards, it is still listed and still
hides them until you remove or disable it.
```

---

## 4) Outlier filter

### UI name / where it lives

When the active field is a categorical obs field of the file, the **Legend**
box shows, below **Show centroids**:

- the `Outlier filter` label with its (i), and then one of:
  - a single **Choose latent space** button, with no slider, until a latent
    space is chosen;
  - a line saying the field's categories are too small to measure;
  - the slider with a percent readout and at most one short line: what the
    first use reads (for example "First use reads X_scVI (12 MB)."), or the
    progress of that read with **Cancel**.

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: The Outlier filter with its (i), its slider moved left of the end stop and its readout showing 95 percent, with the pointer on the slider.
:width: 496px

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
the control is its label, its (i) and one **Choose latent space** button
([Turning it on](02_outlier_filtering_per_active_field.md#turning-it-on-choose-a-latent-space)).

Its quantiles are either stored in the file for that latent space (read at
once, one field at a time) or computed in the browser on the slider's first
move below `100%`, with progress and **Cancel**
([Stored or computed](02_outlier_filtering_per_active_field.md#stored-or-computed)).

A field you create by merging or deleting categories has no outlier
quantiles, so the control is hidden for it.

### Filters list entry

When active and set below 100%, Filters shows a line like:

- `Field: outlier ≤ 95%`

It keeps applying when you colour by another field, as every filter does, and
stays listed under **Filters** until you remove or disable it.

### Enable/disable vs remove

- **Disable** keeps the slider value but stops it from filtering.
- **Remove** resets the outlier slider to `100%`.

See {doc}`02_outlier_filtering_per_active_field` for a deeper explanation and troubleshooting.

---

## Next steps

- {doc}`05_performance_considerations` (what gets slow and why)
- {doc}`06_edge_cases_filtering` (failure modes and surprising states)
