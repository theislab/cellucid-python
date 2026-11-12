# Filter stack UI (Filters panel)

**Audience:** everyone (wet lab + computational users)  
**Time:** 5–10 minutes  
**What you’ll learn:**
- How to find and interpret the Filters list
- How to disable vs remove a filter safely
- How to debug “why did points disappear?” using the filter stack

**Prerequisites:**
- A dataset loaded
- At least one filter applied (any type)

---

## Fast path: “I want my cells back”

1) Open the **Filters** panel.
2) Temporarily disable filters one-by-one (keep notes on which one changes the count).
3) Remove the filter only after you’re sure you don’t need it.
4) If you still see fewer cells than expected, check for:
   - a second filter you forgot about,
   - a per-view filter (if applicable),
   - or a field-based legend toggle that created a filter.

---

## Practical path: what each control means (exact semantics)

The Filters UI is the **source of truth** for current visibility — for
every filter you can create in the UI. The one thing it cannot show is a
per-cell visibility mask set from a notebook with `viewer.cells_visibility_set()`,
which is applied before all of these and has no row here; that is the only way
`Showing X of N cells` and `No filters active` can disagree. See
{doc}`01_filtering_mental_model` and
{doc}`../b_data_loading/05_jupyter_tutorial`.

:::{important}
Its label in the app is **`Filters`**, and its `i` help says what the label
does not: it lists the filters of the **selected view**. In multiview a
different panel can have a completely different list.
:::

It has two parts:

- A **count line**:
  - `Showing all N cells` (no filtering), or
  - `Showing X of N cells` (some cells hidden);
  - on a subsample, N is the open cells of every cell: `Showing all 240 of
    1,200 cells`, `Showing 120 of 240 of 1,200 cells`
    ({doc}`../b_data_loading/21_subsampled_loading`)
- A **filter list** (“stack”):
  - each entry is one active filter in the current view

### Empty state

If no filters are active, the list shows `No filters active.`:

```{figure} ../../../_static/screenshots/filtering/active-filters-empty.png
:alt: The Filters block with its (i), reading Showing all 3,696 cells above the words No filters active.
:width: 496px

The resting state. In the browser, `Showing all N cells` and `No filters
active` always agree — if they ever disagree, you are either looking at a
different view or a notebook has set a visibility mask (see above).
```

Click the small `i` next to the **Filters** heading to see what the list
holds and the scope rule:

```{figure} ../../../_static/screenshots/filtering/active-filters-scope-tooltip.png
:alt: The (i) of the Filters block pressed, and the information box it opens: every filter hiding cells in the selected view, from the legend, the range sliders and the outlier filter, that unticking one pauses it and removing it drops it, that filters belong to the view they were set in, and Learn more links.
:width: 712px

The app's own explanation: what the list holds, what unticking and removing
a row do, and that filters belong to the view they were set in.
```

### Filter rows (what the UI contains)

```{figure} ../../../_static/screenshots/filtering/active-filters-one-row.png
:alt: The Filters block reading Showing 1,876 of 3,696 cells above one row: a ticked checkbox, the text clusters: hiding D… cut off with an ellipsis, the cell count 1,876, and a remove cross at its right.
:width: 496px

One row, left to right: the enable checkbox, the summary text, the cell count,
the `×` that removes it.
```

Each filter row includes:

- a **checkbox**:
  - checked = filter enabled (affects visibility)
  - unchecked = filter disabled (does not affect visibility)
  - its tooltip reads `Turn this filter off` / `Turn this filter on`
- a **text summary** (human-readable):
  - continuous range: `Field: min – max`, each bound to six significant
    figures, so a typed or slider bound reads exactly as applied
  - categorical hides: the smaller side named, `Field: hiding A, B, C, D` or,
    when fewer categories are shown than hidden, `Field: showing only A, B`
    (`Field: hiding every category` when none is) — at most four names, then ` +N more`;
    **Missing** counts like a category when some open cell has no category
  - categories keeping only some of their cells (filter by cell count, the
    legend's count): `Field: keeping 200 of Beta` — at most two, then
    ` +N more`; on the same row as the field's hides, after them
    (`cell_type: hiding Alpha; keeping 200 of Beta`)
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
grayed and struck through**:

```{figure} ../../../_static/screenshots/filtering/active-filters-row-disabled.png
:alt: The same three-row Filters list with the clusters hiding row unchecked under the pointer: that row is grayed and struck through while the other two stay black, and the count line above has risen.
:width: 496px

The middle row is disabled. Its configuration is intact and one click restores
it.
```

### “Remove” vs “Reset” vs “Disable”

There are three common ways to stop a filter from affecting visibility:

1) **Disable** (checkbox in Filters)
   - best for debugging and A/B comparisons
   - preserves the filter configuration (range/category choices)
2) **Reset** (in the legend UI)
   - continuous: click `Reset` in the continuous legend
   - outliers: set outlier slider back to `100%`
   - categorical: click `Show all` in the categorical legend (it also shows
     every cell of a category that kept only some); for one category, click
     its count and type `all`
3) **Remove** (× in Filters)
   - does the reset for you (full range / show all / every cell of every category / outlier=100%)
   - also re-enables the filter (it becomes a “no-op” and typically disappears from the list)
   - the legend redraws to match: range sliders return to the ends, the category checkboxes are ticked again, the counts show every cell, and the outlier slider returns to `100%`

Each of these is a step of the history, so **Undo** brings the filter back
as it was ({doc}`../m_state_history/01_the_history`).

:::{tip}
When debugging, prefer **disable** first (safe and reversible). Use **remove** once you understand the cause.
:::

### Per-filter counts: “visible of available cells”

Some filters show, in a column of their own at the row's right end (never
cut off by a long description), the number of cells they keep, `12,345`;
pointing at it says of how many: `12,345 of 56,789 cells pass this filter`.

Interpretation:

- **available** = how many cells would be visible *if this filter did not exist*, but all other current filters still apply.
- **visible** = how many of those available cells remain visible *with this filter enabled*.

This is extremely useful for diagnosing which filter is “doing the work” in a large filter stack.

```{figure} ../../../_static/screenshots/filtering/active-filters-three-rows.png
:alt: The Filters block reading Showing 768 of 3,696 cells above three rows of different kinds: S_score: -0.146 – 1.… and clusters: hiding Du…, each cut off with an ellipsis and followed by the cell count 768, and clusters: outlier ≤ 95% with no count.
:width: 496px

Three filters of three different kinds on the same 3,696 cells. The first two
carry a cell count in a column of their own; the outlier row carries none.
```

:::{note}
**Category and numeric-range filters carry per-filter counts. Outlier filters
do not** — an outlier row only ever moves the overall `Showing X of Y cells`
line. That is a property of the row, not a sign that the outlier filter is
inactive.
:::

### View/snapshot awareness (easy pitfall)

Filters is **view-scoped**:

- In single-view mode, you only have the live view so it feels global.
- In Live + Snapshots mode, each snapshot can have a different filter stack.

If you’re in grid mode, click the view/snapshot you care about first, then check Filters.

---

## Next steps

- {doc}`04_common_filter_types_document_every_filter_the_ui_exposes` (what filters exist)
- {doc}`07_troubleshooting_filtering` (common “filter stack confusion” failures)
