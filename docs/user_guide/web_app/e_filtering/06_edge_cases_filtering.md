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
- Active filters shows: `Showing 0 of N points`.

```{figure} ../../../_static/screenshots/filtering/window-zero-visible-cells.png
:alt: The whole application window with a completely blank light grid where the embedding was, every legend row unchecked and greyed and reading 0 of its own total, and the sidebar reading Showing 0 of 3,696 points above a single clusters hiding row.
:width: 1440px

An empty canvas produced by pressing `Hide All`. This is a filter state, not a
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

1) In **Active filters**, disable filters one-by-one.
2) When you find the culprit, either:
   - relax it (widen range / show categories / raise outlier percent), or
   - remove it (×), or
   - reset it in the legend (`RESET`, `Show All`, outlier `100%`).

---

## 2) NaN / missing values (continuous filters)

### What you see

- Applying a continuous filter removes more cells than expected.
- Even with Min≈min and Max≈max, some cells never come back while the filter is active.

### Why it happens

Continuous range filters treat `NaN` as “outside range”, so those cells become
invisible whenever the filter is active.

:::{important}
**This cannot happen with a prepared export.** Both writers refuse to prepare a
continuous obs field or a gene that contains a non-finite value, and the refusal
counts them — how many are `NaN` or infinite, and the first affected cells — so
the export never carries one. You can only reach this state by opening a raw `.h5ad` or
Zarr archive directly in the browser, which reads the values as they are.
:::

### What to do

- If you are opening a raw archive: fill or impute missing values upstream, or prepare an export instead, which will tell you exactly which field is at fault.
- If you want to exclude them, current behavior already does this.

---

## 3) “Some categories are disabled / I can’t click them”

### What you see

- In the categorical legend, some category rows are greyed out, their checkbox and colour well disabled, and their count reads `0 cells`.
- Hovering the row shows the tooltip `No cells available in this category after other filters`.

```{figure} ../../../_static/screenshots/filtering/legend-category-unavailable.png
:alt: The clusters legend under another filter, where Beta, Delta, Epsilon and Pre-endocrine are drawn grey with unticked disabled checkboxes and read 0 cells, while Alpha, Ductal, Ngn3 high EP and Ngn3 low EP stay black and ticked with small counts of 1, 122, 10 and 16 cells.
:width: 428px

Four categories with nothing left in them, under a narrow `S_score` range
filter. The rows are **disabled**, not merely empty, and `Show All` / `Hide
All` will skip them.
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
3) Use `Show All` to make categories visible again.

---

## 4) “I disabled a category filter, but those categories look gray”

### What you see

- You disable a category filter in Active filters.
- Cells return, but some categories are colored gray/de-emphasized.

### Why it happens

Cellucid tracks:
- which categories are “visible” in the legend, and
- whether that visibility is currently enforced as a filter.

Disabling the category filter stops it from hiding cells, but the legend may still treat those categories as “hidden” for display purposes (so they appear gray).

### What to do

- Use `Show All` in the categorical legend (or remove the filter) to fully restore.

---

## 5) Outlier slider disappears after editing a categorical field

### What you see

- Outlier filtering was available on a categorical field.
- After merging/deleting categories (creating a derived field), the outlier slider is gone.

### Why it happens

Derived categorical fields typically do **not** carry latent outlier statistics, because those stats are tied to the original category definitions at export time.

### What to do

- Switch back to the original (non-derived) field if you need outlier filtering, or
- re-export data with the desired categorical field defined upstream.

---

## 6) “My gene filter (or outlier filter) stops applying when I switch fields”

### What you see

- You filter by gene expression or outliers.
- You switch to another field.
- The filtered-out cells reappear.

```{figure} ../../../_static/screenshots/filtering/window-gene-range-scope-lost.png
:alt: The whole application window with the Categorical obs select holding clusters and the Gene Expression box empty, the full coloured embedding drawn on the grid, the filter list reading No filters active and the line above it reading Showing all 3,696 points.
:width: 1440px

One field change after an `Ins1` range filter was hiding almost everything.
The row is gone from the list, so the app is not hiding the filter from you —
the filter genuinely stopped applying.
```

### Why it happens

Current behavior:
- gene filters apply only while that gene is the active field,
- outlier filtering applies only while the active field supports outlier stats (that’s when the control is visible).

### What to do

- Treat gene/outlier filtering as “active-field scoped” tools.
- If you need persistent gating while coloring by something else, export the gate as an obs field (e.g., a score/boolean) and filter it as obs.

---

## 7) Multiview: “filters look different across snapshots”

### What you see

- In Live + Snapshots mode, different views show different visible counts.
- Active filters changes when you click different snapshots.

### Why it happens

Each view/snapshot can have its own filter state.

### What to do

- Click the view you care about first, then inspect Active filters.
- If you want two views to match, apply the same filter stack to both (or recreate the snapshot from the filtered live view).

---

## 8) Degenerate fields: one category, or one value

### A categorical field with a single category

Nothing special happens, and that is the point. The legend renders exactly one
row with the whole dataset's count beside it, `Show All` / `Hide All` act on
that one row, and unchecking it takes you straight to `Showing 0 of N points`.
There is no separate “all cells are in one group” state to recognise — it is
the ordinary zero-visible case from section 1 above, reached in one click.

### A continuous field where every cell has the same value

A constant continuous field is a legitimate export: the writers publish it with
equal bounds rather than rejecting it, so you will meet one eventually (a gene
detected in no exported cell is the common case).

Its **Filtering** block is degenerate. Both readouts print the same number,
because the field's minimum *is* its maximum, and there is no interval for the
sliders to divide.

:::{warning}
Do not use the range sliders on a constant field. There is no range to select,
and the slider positions no longer correspond to values in the data. Use a
categorical field, or the outlier filter, to gate cells instead — or filter on
this field upstream, before export, where a constant column is easy to spot.
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
Ultra-light quality — the other two modes rendered identically before and
after.

### What to do

- Re-export any figure captured in Ultra-light quality with a filter active.
- If you are comparing an old figure with a new one, expect the new one to show
  the same or more cells at the same settings.

---

## 10) Quantization and slider resolution (export-time edge case)

### What you see

- Sliders feel “chunky” (small changes don’t affect visibility).
- Very fine thresholds are hard to set.

### Why it happens

Two effects can compound:

- UI sliders are `0–100` percent steps (1% increments).
- Exports can quantize continuous fields and/or outlier quantiles (8-bit/16-bit), reducing numeric precision.

### What to do

- Use UI filtering for coarse/interactive gating.
- Do fine-grained gating upstream in Python before export when precision matters.

---

## Next steps

- {doc}`07_troubleshooting_filtering` (symptom-based debugging)
