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
  - {doc}`02_selection_tools_document_each_tool`

**What you’ll learn**
- What “Highlight mode”, “Highlight pages”, and “Highlighted groups” each control
- How to create/rename/recolor/delete pages
- How to enable/disable/remove groups (and what “Clear” clears)
- How to interpret the highlighted cell count (“visible” vs “total”)
- How to turn your pages into a real categorical obs column with `Create Categorical`

---

## Where the highlight UI lives

In the left sidebar, open the accordion section labeled **Highlighting**.

Inside this section you will see four sub-areas, top to bottom:
1) **Highlight mode** (toolbelt for selection)
2) **Highlight pages** (tabs)
3) **Highlighted groups** (the list + count + Clear)
4) `Create Categorical` (a collapsed accordion that turns pages into an obs column)

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

For the annotation-based tool on categorical fields, `Alt` acts like “replace” (details in {doc}`02_selection_tools_document_each_tool`).

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
You can create a derived page by dragging one page tab onto another, or by
pressing the `⋈` control on a tab.

On drop, the UI offers:
- **Intersection (∩)**: cells that are in *both* pages (enabled groups only)
- **Union (∪)**: cells that are in *either* page (enabled groups only)

This creates a new page whose name reflects the operation (e.g. `Page A ∩ Page B`).

```{figure} ../../../_static/screenshots/highlighting_selection/compare-highlight-pages.png
:alt: The Highlighting panel with two page tabs, a blue Page 1 carrying 870 and an active red Page 2 carrying 250, the group list showing one Proximity group, and the pointer on the button that adds another page.
:width: 516px

Two pages. Each tab carries a colour swatch, the page name, and that page's
total membership. Only the active tab's groups are listed below and only they
are rendered — which is why “Clear did nothing” is almost always “Clear worked,
on the other page”.
```

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
- filters hide some highlighted cells ({doc}`../e_filtering/index`)
- multi-view snapshots may have stricter filters than the live view
- LOD/downsampling is active (large datasets)

---

## 4) Create Categorical (turn pages into a real obs column)

At the bottom of the **Highlighting** section there is a collapsed accordion
headed `Create Categorical`, described as *Build a new categorical obs column
from highlight pages*. It turns the pages you have built into an ordinary
categorical observation field — one you can then colour by, filter on, and use
as the field for annotation-based selection.

It is always present. It is greyed out only until a dataset is loaded.

### The short version

1) Expand `Create Categorical`.
2) Drag page tabs into the `Drop pages here` box (or pick one in
   `Add a highlight page:` and press **Add**). Each page becomes one category,
   named after the page.
3) Type a name in `Column name:`.
4) Press **Create**.

A notification confirms it: `Created "<your name>" (<N> categories)`, and the
new field appears alongside the dataset's own obs fields.

### The two questions it asks you

The builder only shows these when they actually apply, so a set of
non-overlapping pages that covers every cell will show neither.

- **Cells in more than one page.** A section appears headed
  `Assign overlapping cells to:` with four choices: **First page** (the default),
  **Last page**, **New label** — one shared category, `Overlap` unless you rename
  it — and **Each intersection**, which makes a separate category per overlap
  combination and asks you to name each one. Only enabled groups count towards a
  page's membership here, exactly as for page combine.
- **Cells in no page at all.** A section reports `<n> cells not in any page` and
  offers `Label for uncovered cells:`, `Unassigned` by default. Clear that box if
  you would rather those cells carry no label.

:::{note}
**Limits and refusals** (the **Create** button stays disabled and prints the reason):

- at least one page must be added — `Add at least one highlight page`
- the column name must be non-empty, untrimmed-whitespace-free, at most 256
  characters, and must not contain `:`
- the name must not collide with a visible obs field —
  `A visible observation field named "<name>" already exists`
- **Each intersection** supports at most 12 pages, because it enumerates every
  combination of them
- a dataset carries at most 20 user-created fields in total

Fields made this way are saved in a `.cellucid-session` and restored with it, so
a collaborator opening your session on the same dataset gets the same column.
:::

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
- save a session bundle ({doc}`../l_sessions_sharing/index`) and treat it as the persisted artifact of your highlights.

---

## Edge cases and troubleshooting
- Edge cases: {doc}`05_edge_cases_highlighting`
- Troubleshooting: {doc}`06_troubleshooting_highlighting`

## Related pages
- {doc}`01_highlight_mental_model`
- {doc}`02_selection_tools_document_each_tool`
- {doc}`04_selection_synchronization`
