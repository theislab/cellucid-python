# Edge cases (filtering)

**Audience:** computational users + power users  
**Time:** 10–20 minutes  
**What you’ll learn:**
- The most common “surprising” filtering states and how to recognize them
- How missing values affect visibility (NaN and small categories)
- When some filters apply only conditionally (active gene / outlier-supporting fields)
- How multiview changes what “my filters” means
- What degenerate fields look like: one category, one value, or none of a kind

**Prerequisites:**
- A dataset loaded

---

## 1) All cells filtered out (0 visible points)

### What you see

- The canvas looks empty.
- Filters shows: `Showing 0 of N cells`.

```{figure} ../../../_static/screenshots/filtering/window-zero-visible-cells.png
:alt: The whole window with a completely empty canvas after pressing Hide all in the clusters legend: every row unchecked and reading 0 of its count, and the Filters block reading Showing 0 of 3,696 cells above one clusters hiding row.
:width: 1440px

An empty canvas produced by pressing `Hide all`. This is a filter state, not a
crash and not a failed load — and the two lines in the sidebar say so exactly:
the count says `0 of 3,696`, and there is a filter row to remove.
```

### Why it happens

One or more filters reduced visibility to zero:
- a continuous Min/Max range excludes everything,
- categories were hidden until none remain,
- outlier threshold is very low (and outlier filtering is active),
- you are looking at a snapshot with a different filter stack than you expect.

### What to do

Use a deterministic “unwind” sequence:

1) In **Filters**, disable filters one-by-one.
2) When you find the culprit, either:
   - relax it (widen range / show categories / raise outlier percent), or
   - remove it (×), or
   - reset it in the legend (`Reset`, `Show all`, outlier `100%`).

---

## 2) NaN / missing values (continuous filters)

### What you see

- Applying a continuous filter removes more cells than expected.
- Even with Min≈min and Max≈max, some cells never come back while the filter is active.

### Why it happens

Continuous range filters treat `NaN` as “outside range”, so those cells become
invisible whenever the filter is active.

:::{important}
A file written by `prepare()` or `cellucid_prepare()` stores a missing value
in a continuous obs column as `NaN`, so this is ordinary data: those cells are
the ones the column does not measure. Gene expression cannot hold `NaN` in
such a file — both writers refuse it, counting the offending values and
naming the first affected cells.
:::

### What to do

- Fill or impute the missing values upstream if those cells should stay visible, or leave the filter off for that column.
- If you want to exclude them, current behavior already does this.

---

## 3) “Some categories are disabled / I can’t click them”

### What you see

- In the categorical legend, some category rows are greyed out, their checkbox and colour well disabled, and their count reads `0 cells`.
- Hovering the row shows the tooltip `No cells available in this category after other filters`.

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under a narrow S_score range filter: Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay ticked with small counts, while Beta, Delta, Epsilon and Pre-endocrine are greyed out with disabled checkboxes and read 0 cells.
:width: 488px

Four categories with nothing left in them, under a narrow `S_score` range
filter. The rows are **disabled**, not merely empty, and `Show all` / `Hide
all` will skip them.
```

### Why it happens

Cellucid computes:
- **available** counts = after other filters (ignoring this field’s category visibility),
- **visible** counts = after all filters (including this field).

If upstream filters remove all cells in a category, that category becomes unavailable in the current view.

Note the count format changes with the state, which is a useful tell:

| Row reads | Meaning |
|---|---|
| `481 cells` | Nothing else is filtering this category; visible equals available |
| `122 / 916 cells` | 916 cells survive the other filters; 122 of those are drawn |
| `0 cells` | Nothing survives the other filters; the row is disabled |

### What to do

1) Relax upstream filters (continuous ranges, other category hides, gene/outlier filters if active).
2) Revisit the categorical legend.
3) Use `Show all` to make categories visible again.

---

## 4) “I disabled a category filter, but those categories look gray”

### What you see

- You disable a category filter in Filters.
- Cells return, but some categories are colored gray/de-emphasized.

### Why it happens

Cellucid tracks:
- which categories are “visible” in the legend, and
- whether that visibility is currently enforced as a filter.

Disabling the category filter stops it from hiding cells, but the legend may still treat those categories as “hidden” for display purposes (so they appear gray).

### What to do

- Use `Show all` in the categorical legend (or remove the filter) to fully restore.

---

## 5) Outlier slider hidden, replaced by a button, or reading data

### What you see

- After merging/deleting categories (creating a derived field), the outlier
  slider is gone; or
- there is no slider, only a **Choose latent space** button or a line; or
- moving it below `100%` shows a progress bar instead of filtering at once.

### Why it happens

- A field you create by merging or deleting categories has **no** outlier
  quantiles: they exist only for the file's own categorical `obs` columns, so
  the control is hidden for it.
- The filter measures in a latent space, and none is chosen: the control is
  one **Choose latent space** button
  ([Turning it on](02_outlier_filtering_per_active_field.md#turning-it-on-choose-a-latent-space)).
- No category of the field has `centroid_min_points` cells (10 by default):
  there is nothing to measure, and the line says so.
- The file stores no outlier quantiles for this latent space, so the first
  move reads the latent space once and computes them; the line under the
  slider said beforehand how much it reads
  ([Stored or computed](02_outlier_filtering_per_active_field.md#stored-or-computed)).

### What to do

- Switch back to the original (non-derived) field, or add the categorical
  field to `obs` upstream and write the file again.
- Choose a latent space under **Data mapping › Latent space**.
- Lower `centroid_min_points` when writing or serving the file, to measure
  smaller categories.
- Wait for the read, or press **Cancel** (the slider returns to `100%`); write
  the file with `prepare(latent_key=...)` or `cellucid_prepare(latent_key = ...)`
  to store the quantiles so it works at once.

---

## 6) “My gene filter (or outlier filter) stops applying when I switch fields”

### What you see

- You filter by gene expression or outliers.
- You switch to another field.
- The filtered-out cells reappear.

```{figure} ../../../_static/screenshots/filtering/window-gene-range-scope-lost.png
:alt: The same window after switching the active field to clusters: every cell is drawn again, the Ins1 row has left the Filters block, and the count line reads Showing all 3,696 cells.
:width: 1440px

One field change after an `Ins1` range filter was hiding almost everything.
The row is gone from the list, so the app is not hiding the filter from you —
the filter genuinely stopped applying.
```

### Why it happens

Current behavior:
- gene filters apply only while that gene is the active field,
- outlier filtering applies only while the categorical field it was set on is the active field (its threshold is kept for when you come back).

### What to do

- Treat gene/outlier filtering as “active-field scoped” tools.
- If you need persistent gating while coloring by something else, write the gate into the file as an obs field (e.g., a score/boolean) and filter it as obs.

---

## 7) Multiview: “filters look different across snapshots”

### What you see

- In Live + Snapshots mode, different views show different visible counts.
- Filters changes when you click different snapshots.

### Why it happens

Each view/snapshot can have its own filter state.

### What to do

- Click the view you care about first, then inspect Filters.
- If you want two views to match, apply the same filter stack to both (or recreate the snapshot from the filtered live view).

---

## 8) Degenerate fields: one category, or one value

### A categorical field with a single category

Nothing special happens, and that is the point. The legend renders exactly one
row with the whole dataset's count beside it, `Show all` / `Hide all` act on
that one row, and unchecking it takes you straight to `Showing 0 of N cells`.
There is no separate “all cells are in one group” state to recognise — it is
the ordinary zero-visible case from section 1 above, reached in one click.

### A continuous field where every cell has the same value

A constant continuous field is legitimate: the writers publish it with
equal bounds rather than rejecting it, so you will meet one eventually (a gene
detected in no written cell is the common case).

Its **Visible range** block is degenerate. Both readouts print the same number,
because the field's minimum *is* its maximum, and there is no interval for the
sliders to divide.

:::{note}
The sliders and `Reset` still work, but there is nothing to select: every
slider position maps onto the one value, so the filter always keeps every cell.
Use a categorical field, or the outlier filter, to gate cells instead — or
filter on this field upstream, before you write the file, where a constant column is easy
to spot.
:::

### A dataset with no continuous obs fields at all

The **Continuous obs** dropdown shows `(no continuous obs fields)` and holds
nothing to select. There is no failure and no warning: the numeric-range filter
type simply does not exist for that dataset, and the only filters available are
category visibility, outliers, and gene ranges.

The categorical dropdown behaves the same way in the mirror case, showing
`(no categorical obs fields)`.

---

## 9) An older screenshot shows fewer cells than the app does now

### What you see

- A figure you exported earlier, or a screenshot in an old report, shows fewer
  visible cells than the same filter state produces today.

### Why it happens

Hidden points used to be drawn at the driver's minimum point size instead of
being rejected outright. In the **Ultra-light (square points)** shader quality
the fragment shader does not discard, so those hidden points wrote depth and
could occlude visible points sitting behind them. The filtered view was
silently missing real cells.

Hidden points are now rejected in clip space, so nothing they used to cover is
hidden any more. The fix reveals **more** cells, never fewer, and only in
Ultra-light quality — the other modes rendered identically before and
after.

### What to do

- Re-export any figure captured in Ultra-light quality with a filter active.
- If you are comparing an old figure with a new one, expect the new one to show
  the same or more cells at the same settings.

---

## 10) Slider resolution

### What you see

- Sliders feel “chunky” (small changes don’t affect visibility).
- Very fine thresholds are hard to set.

### Why it happens

UI sliders are `0–100` percent steps (1% increments) of the field's full
range.

### What to do

- Use UI filtering for coarse/interactive gating.
- Do fine-grained gating upstream in Python or R when precision matters.

---

## Next steps

- {doc}`07_troubleshooting_filtering` (symptom-based debugging)
