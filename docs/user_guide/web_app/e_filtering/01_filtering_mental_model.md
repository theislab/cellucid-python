# Filtering mental model

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- What “filtering” means in Cellucid (visibility, not deletion)
- Where filters come from (legend controls, outlier slider, gene range filters)
- How the filter stack behaves (AND semantics, enable/disable vs remove)
- How filtering interacts with views/snapshots, highlights, and exports

**Prerequisites:**
- A dataset loaded in the web app (any embedding)

---

## Mental model (one sentence)

In Cellucid, a filter is a **rule that decides whether each cell is visible right now**.

:::{important}
Filtering is **not deletion**.

- Filtered-out cells still exist in the dataset.
- They can reappear instantly when you disable/remove/reset filters.
:::

---

## What “hidden” looks like

This is the single most useful thing to know before you touch a control.

```{figure} ../../../_static/screenshots/filtering/window-clusters-all-visible.png
:alt: The whole application window, sidebar on the left and the Pancreas 3D embedding on a light grid, with eight ticked legend rows carrying cell counts and the line Showing all 3,696 points beneath them.
:width: 1440px

**Before.** All 3,696 cells drawn, every category ticked.
```

```{figure} ../../../_static/screenshots/filtering/window-clusters-three-hidden.png
:alt: The same window after unchecking three legend rows; the orange, purple and blue point groups are absent from the canvas, those rows are greyed and read 0 of 916, 0 of 642 and 0 of 262 cells, and the sidebar reads Showing 1,876 of 3,696 points above one filter row.
:width: 1440px

**After.** Three categories unchecked. The cells are *gone from the picture*,
not faded — and the legend row for each still remembers how many cells the
category holds (`0 / 916 cells`), which is how you know nothing was deleted.
```

---

## Where filters come from (the three surfaces)

Cellucid creates filters through three UI surfaces, all of them inside
**Coloring & Filtering**.

### 1) Categorical legend — checkboxes

```{figure} ../../../_static/screenshots/filtering/type-category-visibility.png
:alt: A legend with Show All and Hide All buttons above eight rows; three rows are unchecked with grey swatches and read 0 of their total cells while five ticked rows read a plain cell count, and each row carries a pencil and a bin icon.
:width: 428px

Unchecking a row creates a **category visibility** filter. `Show All` and
`Hide All` act on every row you are currently allowed to toggle.
```

### 2) Continuous legend — Min/Max sliders

```{figure} ../../../_static/screenshots/filtering/type-numeric-range.png
:alt: A block headed Filtering with the subtitle Visible range, the help line Adjust limits then click FILTER or enable Live filtering, a Live filtering toggle reading Off, the hint Drag sliders then click Filter, Min and Max sliders with numeric readouts, and black FILTER and RESET buttons with the pointer on FILTER.
:width: 428px

Narrowing the range creates a **numeric range** filter. The same block appears
for a continuous obs field and for a gene.
```

### 3) Display options — Outlier filter (latent space)

```{figure} ../../../_static/screenshots/filtering/outlier-slider-95.png
:alt: A control labelled Outlier filter (latent space) whose slider handle sits a short distance in from the right end, with 95 percent printed beside it and the mouse pointer on the handle.
:width: 436px

Lowering the percentage creates an **outlier** filter that hides the cells
furthest from their own category's centre. The control is only present when the
active field carries outlier statistics.
```

There is also a scope rule that catches everyone once:

- **Gene expression range filters apply only while that gene is the active
  field.**
- **Outlier filters apply only while the field they belong to is the active
  field.**

If you ever feel “the plot is lying to me”, the debugging rule is:
**trust the Active filters panel** (it lists what currently affects visibility in the active view).

---

## Fast path (wet lab / non-technical): “hide what I don’t want to see”

1) Pick what you want to look at in **Coloring & Filtering**:
   - clusters / cell types (categorical), or
   - QC metrics / scores (continuous), or
   - a gene (gene expression).
2) Apply a simple filter:
   - **categorical:** uncheck categories you want to hide,
   - **continuous:** narrow the Min/Max slider range,
   - **outliers (if available):** set outlier slider to `95%` to hide the most outlier-like cells.
3) Look at **Active filters**:
   - confirm your filter is listed,
   - confirm “Showing X of Y points” matches your expectation.
4) If you “lost” cells:
   - disable filters one-by-one (checkbox next to each filter),
   - or remove them (×) to reset fully.

---

## Practical path (computational users): exact semantics

### A filter stack is an AND

Within a view, Cellucid’s filter stack is effectively:

> **A cell is visible only if it passes every enabled filter.**

This means:
- adding a filter can only keep the same number of cells or reduce it;
- disabling/removing a filter can only keep the same number of cells or increase it.

### What each filter type considers “missing”

These are easy-to-miss gotchas:

- **Continuous range filters**: a cell whose value is `NaN` is treated as
  “outside range” and becomes invisible whenever that filter is active.
  :::{note}
  A **prepared export cannot contain `NaN` in a continuous obs field or a
  gene** — both writers refuse it, counting the offending values and naming the
  first affected cells. You can only meet this case by opening a raw `.h5ad`
  or Zarr archive in the browser, where the values are read as they are.
  :::
- **Categorical visibility filters**: a cell whose category code is the
  missing-value marker is **skipped** by category filters, so it stays visible
  no matter which categories you untick. It is drawn in the neutral “None” grey
  and has no legend row of its own.
- **Outlier filtering**: a cell without an outlier quantile is never removed by
  outlier filtering. This happens for every cell in a category smaller than
  `centroid_min_points` (default `10`) at export time.

### Enable/disable vs remove (two different workflows)

Cellucid supports both:

- **Disable a filter** (checkbox in Active filters):
  - keeps the filter configured (useful for “what if?” comparisons),
  - but stops it from affecting visibility.
- **Remove a filter** (× in Active filters):
  - resets it to its default “no filtering” state (e.g., show all categories, full numeric range, outlier = `100%`),
  - and re-enables it, so a removed filter usually disappears from the list entirely.

```{figure} ../../../_static/screenshots/filtering/active-filters-row-disabled.png
:alt: Three filter rows where the middle one has an empty checkbox and is greyed and struck through while the other two stay black with ticked checkboxes, the mouse pointer resting on the cleared checkbox, and the line above reading Showing 1,472 of 3,696 points.
:width: 480px

**Disabled**, not removed. The row keeps its configuration and its `×`, and it
is drawn greyed and struck through so it is obvious at a glance which rows are
still doing work. Ticking the box puts the filter back exactly as it was.
```

```{figure} ../../../_static/screenshots/filtering/active-filters-remove-button.png
:alt: A single filter row with the mouse pointer resting on the circled cross button at its right edge.
:width: 480px

**Removed.** The `×` resets that filter to its no-op state — show all
categories, full numeric range, outlier back to `100%` — and the row then
vanishes because there is nothing left to report.
```

:::{tip}
Use **disable** first when debugging. Use **remove** once you’re sure you don’t need that filter state anymore.
:::

### Per-view vs global (live view vs snapshots)

Filtering is **view-scoped**:

- In the normal “single view” workflow, you have one active view (live) so it feels global.
- In **Live + Snapshots** (small multiples), each view/snapshot keeps its own filter state.

Practical implication:
- If a snapshot looks “wrong”, click/select that snapshot first and then check **Active filters** (the panel reflects the currently active view).

---

## Deep path (power users / developers): what filtering costs

Filtering recomputes visibility by scanning all cells and checking each enabled filter.

High-level performance model:
- Visibility recomputation is roughly **O(n_cells × n_enabled_filters)**.
- Cells that are already hidden short-circuit: the mask is tested first, and no
  filter predicate is evaluated for a cell that is hidden anyway.
- Some UI interactions are intentionally optimized:
  - continuous filtering can be applied only on demand (turn off Live filtering → click FILTER),
  - outlier slider updates are throttled to at most one recomputation every 50 ms while dragging.

:::{note}
**Filtering now makes drawing cheaper, not more expensive.** Hidden points are
rejected before rasterisation, so a heavily filtered view draws faster than an
unfiltered one. This has not always been true; see
{doc}`05_performance_considerations` for what changed and why it also fixed a
correctness bug.
:::

For practical guidance, see {doc}`05_performance_considerations`.

---

## Common misconceptions (FAQ-style)

### “Filtering deleted my cells.”

No—filters are reversible visibility rules. Remove or disable them to bring cells back.

### “I cleared the active field (set it to None), but my cells are still missing.”

Clearing the active field affects **coloring** (points become neutral gray), but it does **not** clear filters.

### “My gene filter worked, then stopped working when I switched fields.”

Current behavior: gene-range filters apply only while that gene is the active field (see {doc}`04_common_filter_types_document_every_filter_the_ui_exposes`).

---

## Next steps

- {doc}`02_outlier_filtering_per_active_field` (how outlier filtering is computed and when to use it)
- {doc}`03_filter_stack_ui_active_filters` (how to inspect/disable/remove filters safely)
- {doc}`07_troubleshooting_filtering` (symptom → diagnosis → fix)
