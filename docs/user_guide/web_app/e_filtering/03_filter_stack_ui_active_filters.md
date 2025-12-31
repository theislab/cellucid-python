# Filter stack UI (Active filters)

**Audience:** everyone (wet lab + computational users)  
**Time:** 5–10 minutes  
**What you’ll learn:**
- How to find and interpret the Active filters list
- How to disable vs remove a filter safely
- How to debug “why did points disappear?” using the filter stack

**Prerequisites:**
- A dataset loaded
- At least one filter applied (any type)

---

## Fast path: “I want my cells back”

1) Open the **Active filters** panel.
2) Temporarily disable filters one-by-one (keep notes on which one changes the count).
3) Remove the filter only after you’re sure you don’t need it.
4) If you still see fewer cells than expected, check for:
   - a second filter you forgot about,
   - a per-view filter (if applicable),
   - or a field-based legend toggle that created a filter.

---

## Practical path: what each control means (exact semantics)

The Active filters UI is the **source of truth** for current visibility.

It has two parts:

- A **count line**:
  - `Showing all N points` (no filtering), or
  - `Showing X of N points` (some points hidden)
- A **filter list** (“stack”):
  - each entry is one active filter in the current view

### Empty state

If no filters are active, the list shows:
- `No filters active`

### Filter rows (what the UI contains)

Each filter row includes:

- a **checkbox**:
  - checked = filter enabled (affects visibility)
  - unchecked = filter disabled (does not affect visibility)
- a **text summary** (human-readable):
  - continuous range: `Field: min – max`
  - categorical hides: `Field: hiding A, B, ...`
  - outlier filter: `Field: outlier ≤ 95%`
- a **remove button** `×`:
  - removes the filter and resets it to its default “no filter” state

Disabled filters remain listed and are labeled `(disabled)` so you can re-enable them.

### “Remove” vs “Reset” vs “Disable”

There are three common ways to stop a filter from affecting visibility:

1) **Disable** (checkbox in Active filters)
   - best for debugging and A/B comparisons
   - preserves the filter configuration (range/category choices)
2) **Reset** (in the legend UI)
   - continuous: click `RESET` in the continuous legend
   - outliers: set outlier slider back to `100%`
   - categorical: click `Show All` in the categorical legend
3) **Remove** (× in Active filters)
   - does the reset for you (full range / show all / outlier=100%)
   - also re-enables the filter (it becomes a “no-op” and typically disappears from the list)

:::{tip}
When debugging, prefer **disable** first (safe and reversible). Use **remove** once you understand the cause.
:::

### Per-filter counts: “visible / available cells”

Some filters show a count suffix like:

`• 12,345 / 56,789 cells`

Interpretation:

- **available** = how many cells would be visible *if this filter did not exist*, but all other current filters still apply.
- **visible** = how many of those available cells remain visible *with this filter enabled*.

This is extremely useful for diagnosing which filter is “doing the work” in a large filter stack.

:::{note}
Not every filter type currently shows per-filter counts (for example, outlier filters may only affect the overall “Showing X of Y points” line).
:::

### View/snapshot awareness (easy pitfall)

Active filters is **view-scoped**:

- In single-view mode, you only have the live view so it feels global.
- In Live + Snapshots mode, each snapshot can have a different filter stack.

If you’re in grid mode, click the view/snapshot you care about first, then check Active filters.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: filtering-active-filters-stack
Suggested filename: filtering/02_active-filters-stack.png
Where it appears: User Guide → Web App → Filtering → 03_filter_stack_ui_active_filters.md
Capture:
  - UI location: left sidebar → Active filters
  - State prerequisites: at least 2 filters applied, with one disabled
  - Action to reach state: apply an outlier filter + a categorical hide filter, then disable one
Crop:
  - Include: entire Active filters list + visible count summary
  - Include: enough canvas to show the filtered point cloud
Alt text:
  - Active filters panel showing multiple stacked filters, with one disabled.
Caption:
  - Demonstrate stacked filters and how disable differs from remove.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Active filters stack.
:width: 100%

The Active filters panel shows every filter currently affecting visibility; disable a filter to test its effect before removing it.
```

---

## Next steps

- `04_common_filter_types_document_every_filter_the_ui_exposes` (what filters exist)
- `07_troubleshooting_filtering` (common “filter stack confusion” failures)
