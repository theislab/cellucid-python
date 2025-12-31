# Highlight UI (modes, pages, groups)

This page is the “UI map” for highlighting:
- which buttons exist,
- what state they change,
- and what the app will do when you click **Confirm**, **Clear**, or switch pages.

## At A Glance

**Audience**
- Wet lab / beginner: learn where highlighting lives and how to “save a selection”.
- Computational users: learn how pages/groups interact with filters, multi-view, and analysis.
- Developers: learn which UI actions mutate page/group state and what persists in sessions.

**Time**
- ~15–25 minutes

**Prerequisites**
- A dataset loaded
- Familiarity with selection tools helps:
  - `02_selection_tools_document_each_tool`

**What you’ll learn**
- What “Highlight mode”, “Highlight pages”, and “Highlighted groups” each control
- How to create/rename/recolor/delete pages
- How to enable/disable/remove groups (and what “Clear” clears)
- How to interpret the highlighted cell count (“visible” vs “total”)

---

## Where the highlight UI lives

In the left sidebar, open the accordion section labeled **Highlighting**.

Inside this section you will typically see three sub-areas:
1) **Highlight mode** (toolbelt for selection)
2) **Highlight pages** (tabs)
3) **Highlighted groups** (the list + count + Clear)

---

## 1) Highlight mode (choose how you build a selection)

The **Highlight mode** buttons control *which selection tool is active*:
- **Annotation based**
- **KNN drag**
- **Proximity drag**
- **Lasso**

All of these tools create a temporary **candidate set** and then let you **Confirm** it into a new highlight group.

### Step controls appear after you start selecting
Once you perform a selection gesture (e.g., finish a lasso), the UI shows “step controls” under the mode description:
- **Confirm**: create a new group from the current candidate set
- **Undo / Redo**: roll back steps (limited history)
- **Cancel**: discard the candidate set

These controls are tool-specific but behave consistently across tools.

### Modifier keys (all tools)
- `Alt` starts the selection gesture.
- `Shift+Alt` usually means “add” (union).
- `Ctrl+Alt` or `Cmd+Alt` usually means “remove” (subtract).

For the annotation-based tool on categorical fields, `Alt` acts like “replace” (details in `02_selection_tools_document_each_tool`).

---

## 2) Highlight pages (tabs = alternative “workspaces”)

Highlight pages are the top-level organization unit.

Each page contains its own list of highlight groups.
Only one page is **active** for highlight rendering at a time.

### Create a new page
Click the **`+`** button next to the page tabs.

Expected behavior:
- a new page is created,
- the app switches to it immediately,
- the group list becomes empty (because the new page has no groups yet).

### Rename a page
Double-click the page name in the tab strip.

Typical editing behavior:
- press `Enter` to commit,
- press `Escape` to cancel,
- clicking away (blur) commits the current text.

### Recolor a page
Each page tab has a color indicator (a small color swatch).

Click it to open a color picker and choose a new color.

The page color is used for:
- quick visual identification of pages,
- some analysis/UI components that reference pages.

### Delete a page
If you have more than one page, an **`×`** delete control appears on tabs.

Rules:
- you cannot delete the last remaining page,
- if you delete the currently active page, the UI switches to another page automatically.

### Combine pages (union / intersection)
You can create a derived page by dragging one page tab onto another.

On drop, the UI offers:
- **Intersection (∩)**: cells that are in *both* pages (enabled groups only)
- **Union (∪)**: cells that are in *either* page (enabled groups only)

This creates a new page whose name reflects the operation (e.g. `Page A ∩ Page B`).

:::{note}
If the intersection is empty, the derived page can be created with zero groups (an “empty result”). This is expected.
:::

---

## 3) Highlighted groups (the actual persistent selections)

Highlight groups live inside the currently active page.

### Create a group (the normal way)
To create a new highlight group:
1) choose a highlight mode (lasso/proximity/KNN/annotation)
2) build a candidate set (possibly multi-step)
3) click **Confirm**

Result:
- a new group appears in the list
- highlighted cells are emphasized in the canvas

### Enable/disable a group
Each group has a checkbox:
- checked = group contributes to highlighting and analysis
- unchecked = group is ignored (but kept)

This is the safest way to temporarily remove a group from consideration without deleting it.

### Remove a group
Click the **`×`** next to a group to remove it from the page.

There is no undo for deleting groups; if you need safety, save a session first.

### Clear all groups (what “Clear” clears)
The **Clear** button clears all groups on the **active page**.

It does not delete other pages.

Keyboard shortcut:
- press `x` (when you are not typing in a text field) to clear highlights.

:::{important}
If you have multiple highlight pages, “Clear” can feel like it “didn’t work” because you may be looking at a different page. Always confirm the active page tab.
:::

---

## Reading the highlighted cell count (“visible” vs “total”)

Above the group list, Cellucid displays a summary like:
- “No cells highlighted”
- “12,345 cells highlighted”
- “1,234 of 12,345 highlighted cells visible”

Interpretation:
- **total** = membership in enabled groups on the active page
- **visible** = how many of those highlighted cells are currently visible in the current view (after filters and, for very large datasets, LOD/downsampling)

Common reasons “visible < total”:
- filters hide some highlighted cells (`../e_filtering/index`)
- multi-view snapshots may have stricter filters than the live view
- LOD/downsampling is active (large datasets)

---

## What you cannot (yet) do in the highlight UI

Depending on what you’ve used in other tools, you might look for these features:
- rename highlight groups
- recolor highlight groups
- “add selection to an existing group” after confirming
- export highlighted cell lists directly from the UI

In the current UI, group membership is edited by:
- rebuilding the candidate set (union/subtract/intersect), then confirming a new group, and optionally removing the old one.

For exporting, the most robust current path is:
- save a session bundle (`../l_sessions_sharing/index`) and treat it as the persisted artifact of your highlights.

---

## Advanced (may be dev-phase): create a categorical field from pages

Some builds include a “Category builder” UI under Highlighting that lets you:
- drag highlight pages into a builder,
- choose an overlap resolution strategy,
- create a new categorical field from page membership.

If you don’t see it, it may be disabled or hidden in your build.

---

## Edge cases and troubleshooting
- Edge cases: `05_edge_cases_highlighting`
- Troubleshooting: `06_troubleshooting_highlighting`

## Related pages
- `01_highlight_mental_model`
- `02_selection_tools_document_each_tool`
- `04_selection_synchronization`
