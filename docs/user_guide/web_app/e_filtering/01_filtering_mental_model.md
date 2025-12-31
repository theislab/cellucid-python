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

## Where filters come from (the “three surfaces”)

Cellucid creates filters through three main UI surfaces:

1) **Categorical legend** (checkboxes)
   - Hiding categories creates a “category visibility” filter.
2) **Continuous legend** (Min/Max sliders)
   - Adjusting the visible range creates a numeric range filter.
3) **Display options → Outlier filter (latent space)** (percent slider, when available)
   - Reduces visual clutter by hiding the most outlier-like cells for fields that provide latent outlier stats.

There is also a special case:
- **Gene expression range filters apply only while that gene is the active field.**

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

- **Continuous range filters**: cells with `NaN` (and `null`) are treated as “outside range” and become invisible whenever that filter is active.
- **Categorical visibility filters**: cells with missing/invalid category codes are typically not removed by “hide category” rules (they behave like “None/unknown” and can remain visible).
- **Outlier filtering**: cells without an outlier quantile value (e.g., rare categories below a minimum size during export) are not removed by outlier filtering.

### Enable/disable vs remove (two different workflows)

Cellucid supports both:

- **Disable a filter** (checkbox in Active filters):
  - keeps the filter configured (useful for “what if?” comparisons),
  - but stops it from affecting visibility.
- **Remove a filter** (× in Active filters):
  - resets it to its default “no filtering” state (e.g., show all categories, full numeric range, outlier = `100%`).

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

## Deep path (power users / developers): why filtering can be “expensive”

Filtering recomputes visibility by scanning all cells and checking each enabled filter.

High-level performance model:
- Visibility recomputation is roughly **O(n_cells × n_enabled_filters)**.
- Some UI interactions are intentionally optimized:
  - continuous filtering can be applied only on demand (turn off Live filtering → click FILTER),
  - outlier slider updates are throttled to avoid too many recomputations while dragging.

For practical guidance, see `05_performance_considerations`.

---

## Common misconceptions (FAQ-style)

### “Filtering deleted my cells.”

No—filters are reversible visibility rules. Remove or disable them to bring cells back.

### “I cleared the active field (set it to None), but my cells are still missing.”

Clearing the active field affects **coloring** (points become neutral gray), but it does **not** clear filters.

### “My gene filter worked, then stopped working when I switched fields.”

Current behavior: gene-range filters apply only while that gene is the active field (see `04_common_filter_types_document_every_filter_the_ui_exposes`).

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: filtering-mental-model-active-filters
Suggested filename: filtering/00_active-filters-and-count.png
Where it appears: User Guide → Web App → Filtering → 01_filtering_mental_model.md
Capture:
  - UI location: left sidebar → Active filters
  - State prerequisites: dataset loaded; at least one categorical hide and one continuous range filter active
  - Action to reach state: hide a category and apply a continuous min/max range, then open Active filters
Crop:
  - Include: “Showing X of Y points” line + at least two filter entries
  - Include: enough canvas to show that some points are missing
Alt text:
  - Active filters panel showing multiple filters and the current visible point count.
Caption:
  - Active filters is the source of truth for what currently affects visibility in the active view.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Active filters panel.
:width: 100%

Active filters lists the current filter stack and shows “Showing X of Y points” for the active view.
```

---

## Next steps

- `02_outlier_filtering_per_active_field` (how outlier filtering is computed and when to use it)
- `03_filter_stack_ui_active_filters` (how to inspect/disable/remove filters safely)
- `07_troubleshooting_filtering` (symptom → diagnosis → fix)
