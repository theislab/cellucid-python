# Edge cases (highlighting and selection)

This page is a “weird but expected” checklist for highlighting/selection. It is not a bug tracker; it’s the set of behaviors that commonly surprise users.

If you’re actively stuck, use {doc}`06_troubleshooting_highlighting` (it’s organized by symptoms and fixes).

## Edge case index (scan first)

- [Zero visible cells](#highlighting-with-zero-visible-cells)
- [Changing the active field](#highlighting-after-field-changes-active-field-changed)
- [No active field at all](#annotation-based-with-no-field-at-all)
- [Cells in multiple groups](#cells-in-multiple-highlight-groups)
- [Switching pages mid-selection](#switching-highlight-pages-while-a-selection-is-active)
- [Huge highlight groups (1M+)](#huge-highlight-groups-performance-memory-and-session-size)
- [Visual ambiguity (it “doesn’t look highlighted”)](#visual-ambiguity-highlight-style-vs-color-by)
- [3D lasso pitfalls](#3d-selection-projection-pitfalls-lasso)
- [Multi-view pitfalls](#multi-view-live--snapshots)
- [KNN graph not available](#knn-selection-when-the-neighbor-graph-is-missing)
- [Reload/dataset changes](#dataset-changes-and-persistence)

---

## Highlighting with zero visible cells

**What you see**
- You have highlight groups, but the canvas looks empty or unhighlighted.
- The highlight count can read something like “0 of N highlighted cells visible”.

**Why this happens**
- Filters can hide *all* cells in the current view ({doc}`../e_filtering/index`).
- In very large datasets, LOD/downsampling can reduce the number of drawn points when zoomed out (so “visible” can drop even if you did not change filters).

**What to do**
- First, disable filters (or reset filter stack) to confirm whether highlights exist but are hidden.
- If you’re in multi-view, confirm which panel you’re looking at: different snapshot views can have different visibility.

**What is expected**
- Highlight group membership still exists even when no cells are visible; it is a membership layer, not a visibility layer.

---

## Highlighting after field changes (active field changed)

**What changes**
- The *look* of the plot changes (because color-by changed).
- Annotation-based selection behavior can change (because it depends on the active field’s kind: categorical vs continuous).

**What does not change**
- Highlight pages and highlight groups (membership) do not change when you switch fields.

**Common gotcha**
- You are in **Annotation based** mode, and you switch to a field that is not selectable (or has missing values). `Alt+click` may appear to “stop working” because the tool has no meaningful selection rule for that field.

**Best practice**
- If you are mid-selection, confirm/cancel before switching the active field to avoid mixing steps driven by different semantics.

## Annotation based with no field at all

Set both **Categorical obs** and **Continuous obs** to `None` and the tool has
nothing to work from. It says so rather than failing silently:

```{figure} ../../../_static/screenshots/highlighting_selection/read-no-active-field-notice.png
:alt: The Highlighting panel in Annotation based mode with no coloured field, showing the line that annotation selection needs an active field chosen under Coloring and Filtering.
:width: 516px

With no coloured field the usual `Alt+click` invitation is replaced by an
instruction. In this state the viewer never starts the gesture, so there is
nothing to debug beyond choosing a field.
```

The other three tools are unaffected: lasso, proximity and KNN are geometric or
graph-based and do not read the coloured field at all.

---

## Cells in multiple highlight groups

**Current behavior**
- Cellucid’s highlight rendering is effectively boolean: if a cell is in *any enabled group* on the active page, it is highlighted.
- There is no “priority” or “which group wins” visual rule because groups do not have distinct colors/styles yet.

**What this implies**
- Groups are combined as a union for rendering.
- Disabling a group only removes its contribution; if a cell is also in another enabled group, it stays highlighted.

**How to reason about it**
- If you care about mutually exclusive labels, use pages (separate workspaces),
  then turn them into a real categorical column with `Create Categorical`, which
  makes you decide explicitly what happens to a cell that is in two of them —
  see {doc}`03_highlight_ui`.

---

## Switching highlight pages while a selection is active

This is a common source of “my selection vanished”.

**What can happen**
- The selection candidate set may still exist (step controls may still show a step count).
- The preview highlight can disappear or change because switching pages recomputes permanent highlighting.
- If you click **Confirm** after switching pages, the new group is created in the *currently active page* (which may be different from where you started).

**Best practice**
- Confirm or cancel before switching pages.

**If you already did it**
- Switch back to the original page and look for a new group (if you confirmed).
- If you didn’t confirm, try switching tools (or re-enter the same tool) to see whether the candidate set UI reappears.

---

## Huge highlight groups (performance, memory, and session size)

Huge highlights are possible, but you should expect performance constraints.

**What gets slower**
- Selection preview during lasso/proximity drags (may repeatedly scan many points).
- Recomputing highlight buffers when you add/remove/enable/disable groups.
- Saving/loading sessions (highlight memberships are stored in the session bundle).

**What you’ll notice**
- “Confirm” feels delayed (UI pauses before the group appears).
- Toggling group checkboxes becomes sluggish.
- Session files become large and take longer to load/restore.

**Recommended workflows for 1M+ cells**
- Prefer filters to reduce the visible set before selecting.
- Use fewer groups; keep groups coarse and page-level organization clean.
- Avoid long continuous drags that force many preview recalculations.
- Save sessions intentionally (expect large files).

**Developer note (why)**
Highlight groups store explicit arrays of indices, and highlight recomputation scales with the total number of stored indices across enabled groups.

---

## Visual ambiguity: highlight style vs color-by

**Symptoms**
- “Nothing is highlighted” (but the count says there are highlighted cells).
- “I can’t tell what’s highlighted” (low contrast).

**Why**
- Highlighting is rendered as an overlay (halo/ring) on top of whatever color-by is showing.
- Depending on point size, background, and color palette, the highlight overlay can be subtle.

**What to try**
- Zoom in: highlights are easiest to see when points are not subpixel-sized.
- Raise **Point size (log)** under **Visualization** to make the drawn points —
  and therefore their halos — larger.
- Switch to a simpler color-by (e.g., a single color) just to visually validate highlights.
- Change background mode to increase contrast.

:::{warning}
**Point size** changes only what is *drawn*. It has no effect on what can be
*clicked*: picking uses a fixed radius in data space and never reads the drawn
size. If an `Alt+click` is missing cells, zoom in — do not reach for the point
size slider. See {doc}`06_troubleshooting_highlighting`.
:::

---

## 3D selection projection pitfalls (lasso)

**The core issue**
Lasso is a 2D polygon drawn on the screen. In 3D, “inside the polygon” is ambiguous because many points at different depths can project into the same 2D region.

**What you’ll see**
- You lasso a region that looks like one cluster, but additional points are included.
- You lasso “the same region” after rotating and get a different result.

**The intended fix**
- Use multi-step intersection:
  1) lasso from angle A (`Alt`)
  2) rotate substantially
  3) lasso again with `Alt` to intersect

This approximates a 3D volume by combining multiple 2D “cuts”.

---

## Multi-view (live + snapshots)

**What is global**
- Highlight pages and groups are global (dataset-scoped).
- Confirmed highlights apply to all views.

**What is per-view**
- Visibility (filters/transparency) can differ per view.
- Some snapshots may “share” the live view’s visibility, while others may have their own.

**Common surprises**
- “It highlights in the live view but not in a snapshot”: snapshot is filtering out those cells.
- “The count says visible < total”: you’re seeing visibility differences (filters or LOD).

---

## KNN selection when the neighbor graph is missing

**What you see**
- KNN drag mode shows a message like “Neighbor graph not available”, or KNN selection does nothing.

**Why**
- KNN selection requires a connectivity graph.
- The graph may not be present in the dataset export, or it may not be loaded yet.

**What to do**
- Enable/show edges (if the UI has a “Show edges” toggle).
- Re-export your dataset including connectivities.
- Use lasso/proximity/annotation-based selection if connectivity is not available.

---

## Dataset changes and persistence

### Reloading the page
If you reload without restoring a session, highlights typically reset (highlighting is session-scoped, not automatically URL-scoped).

### Loading a different dataset
Highlights do not transfer between datasets because they are stored as cell indices for a specific dataset ordering.

### Restoring a session on the “wrong” dataset
When dataset fingerprints do not match, the current reader rejects and rolls
back the complete session. It never salvages layout or silently omits
highlights. Load the exact dataset identity, then restore the unchanged session.

---

## Related pages
- {doc}`01_highlight_mental_model`
- {doc}`04_selection_synchronization`
- {doc}`06_troubleshooting_highlighting`
