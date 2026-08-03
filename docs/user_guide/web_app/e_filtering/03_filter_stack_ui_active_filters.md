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

The Active filters UI is the **source of truth** for current visibility — for
every filter you can create in the UI. The one thing it cannot show is a
per-cell visibility mask set from a notebook with `viewer.set_visibility()`,
which is applied before all of these and has no row here; that is the only way
`Showing X of N points` and `No filters active` can disagree. See
{doc}`01_filtering_mental_model` and
{doc}`../b_data_loading/05_jupyter_tutorial`.

:::{important}
Its label in the app is **`Active filters (selected view only):`**, not just
“Active filters”. The parenthesis is the whole story: this block describes one
view, and in multiview a different panel can have a completely different list.
:::

It has two parts:

- A **count line**:
  - `Showing all N points` (no filtering), or
  - `Showing X of N points` (some points hidden)
- A **filter list** (“stack”):
  - each entry is one active filter in the current view

### Empty state

If no filters are active, the list shows `No filters active`:

```{figure} ../../../_static/screenshots/filtering/active-filters-empty.png
:alt: A block headed ACTIVE FILTERS (SELECTED VIEW ONLY) with a small circled i button beside the heading, the line Showing all 3,696 points, and a bordered area containing the words No filters active.
:width: 480px

The resting state. In the browser, `Showing all N points` and `No filters
active` always agree — if they ever disagree, you are either looking at a
different view or a notebook has set a visibility mask (see above).
```

Click the small `i` next to the heading to see why the label says
*selected view only*:

```{figure} ../../../_static/screenshots/filtering/active-filters-scope-tooltip.png
:alt: The Active filters heading with its circled i button pressed and ringed, and a pop-up panel opened to its right reading that filters stay with the selected panel when its coloring changes while other panels remain unchanged.
:width: 908px

The app's own explanation of filter scope, in its own words.
```

### Filter rows (what the UI contains)

```{figure} ../../../_static/screenshots/filtering/active-filters-one-row.png
:alt: The Active filters block reading Showing 1,876 of 3,696 points above a single row containing a ticked blue checkbox, the text clusters colon hiding Ductal, Ngn3 cut off with an ellipsis, and a circled cross button at the right edge.
:width: 480px

One row, left to right: the enable checkbox, the summary text, the `×` that
removes it.
```

Each filter row includes:

- a **checkbox**:
  - checked = filter enabled (affects visibility)
  - unchecked = filter disabled (does not affect visibility)
  - its tooltip reads `Disable this filter` / `Enable this filter`
- a **text summary** (human-readable):
  - continuous range: `Field: min – max`, both bounds to two decimal places
  - categorical hides: `Field: hiding A, B, C, D` — at most four names, then
    ` +N more`
  - outlier filter: `Field: outlier ≤ 95%`
- a **remove button** `×` (tooltip `Remove this filter`):
  - removes the filter and resets it to its default “no filter” state

:::{tip}
**The row text is usually too long for the sidebar and is cut off with an
ellipsis.** The full, untruncated text is the row's tooltip — hover the row to
read it. This is also why a `(disabled)` suffix is often invisible; see below.
:::

Disabled filters stay in the list. Internally their summary text gains a
` (disabled)` suffix, but the reliable visual cue is that **the whole row is
greyed and struck through**:

```{figure} ../../../_static/screenshots/filtering/active-filters-row-disabled.png
:alt: Three filter rows where the middle row has an empty checkbox and is drawn grey and struck through while the rows above and below stay black with ticked checkboxes, the mouse pointer on the cleared checkbox.
:width: 480px

The middle row is disabled. Its configuration is intact and one click restores
it.
```

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

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The Active filters block reading Showing 709 of 3,696 points above three ticked rows; the first reads S_score minus 0.15 to 1.14 followed by a cell count, the second reads clusters colon hiding Ductal, Epsilon cut off with an ellipsis, and the third reads clusters colon outlier less than or equal to 95 percent with nothing after it.
:width: 480px

Three filters of three different kinds on the same 3,696 cells. The first two
carry a `visible / available` suffix; the outlier row carries none.
```

:::{note}
**Category and numeric-range filters carry per-filter counts. Outlier filters
do not** — an outlier row only ever moves the overall `Showing X of Y points`
line. That is a property of the row, not a sign that the outlier filter is
inactive.
:::

### View/snapshot awareness (easy pitfall)

Active filters is **view-scoped**:

- In single-view mode, you only have the live view so it feels global.
- In Live + Snapshots mode, each snapshot can have a different filter stack.

If you’re in grid mode, click the view/snapshot you care about first, then check Active filters.

---

## Next steps

- {doc}`04_common_filter_types_document_every_filter_the_ui_exposes` (what filters exist)
- {doc}`07_troubleshooting_filtering` (common “filter stack confusion” failures)
