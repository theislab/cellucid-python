# Screenshots

This page is a **screenshot production checklist** for highlighting and selection.

All blocks below use the standard “screenshot placeholder” spec from:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

Recommended destination folder for real images:
- `cellucid-python/docs/_static/screenshots/highlighting_selection/`

---

## Orientation (where is Highlighting in the UI?)

<!-- SCREENSHOT PLACEHOLDER
ID: highlighting-orientation
Suggested filename: highlighting_selection/01_highlighting-ui-entry-point.png
Where it appears: Highlighting and Selection → Orientation
Capture:
  - UI location: left sidebar → Highlighting accordion section
  - State prerequisites: dataset loaded; points visible
  - Action to reach state: open Highlighting (ensure it is expanded)
Crop:
  - Include: highlight mode buttons, highlight page tabs, highlight count, and 1–2 groups if present
  - Include: enough of the canvas to orient the reader
Redact:
  - Remove: private dataset name/ID if needed
Annotations:
  - Callouts: highlight mode toolbelt; page tabs; group list; Clear button
Alt text:
  - Left sidebar showing the Highlighting section with mode buttons, page tabs, and highlighted groups.
Caption:
  - The Highlighting section contains the selection toolbelt (modes), page tabs, and the highlighted group list.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Highlighting UI entry point.
:width: 100%

The Highlighting section contains the selection toolbelt (modes), page tabs, and the highlighted group list.
```

---

## Lasso selection (before/after)

<!-- SCREENSHOT PLACEHOLDER
ID: lasso-before
Suggested filename: highlighting_selection/02_lasso-before.png
Where it appears: Highlighting and Selection → Lasso → Before
Capture:
  - UI location: canvas + Highlighting section visible (optional but helpful)
  - State prerequisites: dataset loaded; lasso mode selected; no confirmed group yet
  - Action to reach state: select Lasso mode; hold Alt; start drawing a lasso but do not release
Crop:
  - Include: the lasso polygon/line and enough points to show the intended region
Alt text:
  - 2D embedding with a lasso outline being drawn around a cluster.
Caption:
  - Hold Alt and draw a lasso; the outline shows the tentative region during selection.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing a lasso outline before finishing the selection.
:width: 100%

Hold Alt and draw a lasso; the outline shows the tentative region during selection.
```

<!-- SCREENSHOT PLACEHOLDER
ID: lasso-after-step-controls
Suggested filename: highlighting_selection/03_lasso-after-step-controls.png
Where it appears: Highlighting and Selection → Lasso → After
Capture:
  - UI location: Highlighting section + canvas
  - State prerequisites: you released the lasso, candidate set is non-empty, step controls are visible
  - Action to reach state: finish lasso (release mouse)
Crop:
  - Include: step controls (Confirm/Undo/Redo/Cancel) and the canvas showing preview-highlighted points
Alt text:
  - Lasso selection completed with step controls visible and selected points preview-highlighted.
Caption:
  - After releasing the lasso, Cellucid shows step controls so you can confirm, undo/redo, or cancel the selection.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing lasso step controls after completing a lasso step.
:width: 100%

After releasing the lasso, Cellucid shows step controls so you can confirm, undo/redo, or cancel the selection.
```

<!-- SCREENSHOT PLACEHOLDER
ID: lasso-after-confirm
Suggested filename: highlighting_selection/04_lasso-after-confirm.png
Where it appears: Highlighting and Selection → Lasso → Confirmed group
Capture:
  - UI location: Highlighting group list + canvas
  - State prerequisites: you clicked Confirm and created a new highlight group
  - Action to reach state: finish lasso → Confirm
Crop:
  - Include: the new group entry in the list and the highlighted region in the canvas
Alt text:
  - A new highlight group appears in the group list after confirming a lasso selection.
Caption:
  - Confirm converts the current selection candidate set into a persistent highlight group on the active page.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing a confirmed highlight group after lasso selection.
:width: 100%

Confirm converts the current selection candidate set into a persistent highlight group on the active page.
```

---

## Highlight groups list (2–3 groups + enable/disable)

<!-- SCREENSHOT PLACEHOLDER
ID: highlight-groups-list
Suggested filename: highlighting_selection/05_highlight-groups-list.png
Where it appears: Highlighting and Selection → Highlighted groups
Capture:
  - UI location: Highlighted groups list
  - State prerequisites: at least 2–3 groups exist; at least one group disabled (checkbox unchecked)
  - Action to reach state: create multiple groups; disable one via checkbox
Crop:
  - Include: group labels, checkboxes, remove (×), and the Clear button
Alt text:
  - Highlighted groups list showing multiple groups with checkboxes and a Clear button.
Caption:
  - The group list lets you enable/disable groups (without deleting) and remove or clear groups on the active page.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the highlighted groups list with multiple groups.
:width: 100%

The group list lets you enable/disable groups (without deleting) and remove or clear groups on the active page.
```

---

## Highlight pages UI (create/rename/recolor)

<!-- SCREENSHOT PLACEHOLDER
ID: highlight-pages-tabs
Suggested filename: highlighting_selection/06_highlight-pages-tabs.png
Where it appears: Highlighting and Selection → Highlight pages
Capture:
  - UI location: highlight page tab strip
  - State prerequisites: at least 2 pages exist; page colors visible
  - Action to reach state: click + to add a second page; switch between them
Crop:
  - Include: page tabs, + add button, and the active tab styling
Alt text:
  - Highlight page tabs showing multiple pages and an add-page button.
Caption:
  - Highlight pages are alternative workspaces; switching pages changes which highlight groups are active.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the highlight page tabs.
:width: 100%

Highlight pages are alternative workspaces; switching pages changes which highlight groups are active.
```

<!-- SCREENSHOT PLACEHOLDER
ID: highlight-page-rename
Suggested filename: highlighting_selection/07_highlight-page-rename.png
Where it appears: Highlighting and Selection → Rename page
Capture:
  - UI location: highlight page tabs
  - State prerequisites: rename input visible (editing a tab name)
  - Action to reach state: double-click a page name to rename
Crop:
  - Include: the inline rename input and any relevant surrounding context
Alt text:
  - Inline text input for renaming a highlight page tab.
Caption:
  - Double-click a page tab name to rename the page (use clear names like “Control” or “Treated”).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing highlight page rename UI.
:width: 100%

Double-click a page tab name to rename the page (use clear names like “Control” or “Treated”).
```

---

## Page combine menu (union / intersection)

<!-- SCREENSHOT PLACEHOLDER
ID: highlight-page-combine-menu
Suggested filename: highlighting_selection/08_highlight-page-combine-menu.png
Where it appears: Highlighting and Selection → Combine pages
Capture:
  - UI location: page combine menu overlay
  - State prerequisites: at least 2 pages exist
  - Action to reach state: drag one page tab onto another and drop
Crop:
  - Include: the menu showing Intersection (∩) and Union (∪) options
Alt text:
  - Page combine menu offering intersection and union operations after dragging one page onto another.
Caption:
  - Drag-and-drop page combine creates a derived page by union or intersection of enabled highlights.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the page combine menu.
:width: 100%

Drag-and-drop page combine creates a derived page by union or intersection of enabled highlights.
```

---

## “Visible vs total” highlight count (filters hide highlighted cells)

<!-- SCREENSHOT PLACEHOLDER
ID: highlight-visible-vs-total
Suggested filename: highlighting_selection/09_visible-vs-total-count.png
Where it appears: Highlighting and Selection → Visibility interplay
Capture:
  - UI location: highlight count line (above group list)
  - State prerequisites: at least one highlight group exists AND a filter hides some highlighted cells
  - Action to reach state:
      1) create a highlight group
      2) apply a filter that hides part of that group
Crop:
  - Include: highlight count text showing “X of Y highlighted cells visible”
Alt text:
  - Highlight count showing that only a subset of highlighted cells are visible due to filtering.
Caption:
  - The highlight count distinguishes total highlighted membership from the subset currently visible in the view.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for visible vs total highlighted count.
:width: 100%

The highlight count distinguishes total highlighted membership from the subset currently visible in the view.
```

---

## Proximity and KNN (recommended)

<!-- SCREENSHOT PLACEHOLDER
ID: proximity-overlay
Suggested filename: highlighting_selection/10_proximity-overlay.png
Where it appears: Highlighting and Selection → Proximity drag
Capture:
  - UI location: canvas overlay during proximity drag
  - State prerequisites: proximity mode selected; Alt held; drag in progress
  - Action to reach state: Alt+click a cell, drag to expand radius
Crop:
  - Include: the proximity circle overlay and enough points to interpret selection
Alt text:
  - Proximity drag overlay showing a circular radius around a seed cell.
Caption:
  - Proximity drag selects cells within an embedding-space radius around a seed point.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for proximity drag overlay.
:width: 100%

Proximity drag selects cells within an embedding-space radius around a seed point.
```

<!-- SCREENSHOT PLACEHOLDER
ID: knn-degree-overlay
Suggested filename: highlighting_selection/11_knn-degree-overlay.png
Where it appears: Highlighting and Selection → KNN drag
Capture:
  - UI location: canvas overlay during KNN drag
  - State prerequisites: KNN mode selected; edges loaded; Alt held; drag in progress
  - Action to reach state: Alt+click a seed cell, drag to increase degree
Crop:
  - Include: degree indicator and enough points to interpret selection
Alt text:
  - KNN drag overlay showing a seed cell and a degree indicator during neighbor expansion.
Caption:
  - KNN drag grows selection along the neighbor graph; degree increases as you drag farther.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for KNN drag degree overlay.
:width: 100%

KNN drag grows selection along the neighbor graph; degree increases as you drag farther.
```

---

## Common failure states (high value)

<!-- SCREENSHOT PLACEHOLDER
ID: knn-graph-missing
Suggested filename: highlighting_selection/12_knn-graph-missing.png
Where it appears: Troubleshooting → KNN not available
Capture:
  - UI location: highlight mode description text
  - State prerequisites: KNN mode selected but neighbor graph not loaded/available
  - Action to reach state: select KNN mode on a dataset without connectivities (or before enabling edges)
Crop:
  - Include: the KNN mode description showing “neighbor graph not available” (or equivalent)
Alt text:
  - KNN selection mode showing a message that the neighbor graph is not available.
Caption:
  - KNN selection requires a connectivity graph; if it is missing, use other selection tools or re-export with connectivities.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for missing neighbor graph in KNN mode.
:width: 100%

KNN selection requires a connectivity graph; if it is missing, use other selection tools or re-export with connectivities.
```

<!-- SCREENSHOT PLACEHOLDER
ID: filtered-out-all-highlights
Suggested filename: highlighting_selection/13_filtered-out-all-highlights.png
Where it appears: Troubleshooting → “Nothing highlights”
Capture:
  - UI location: highlight count + filter summary
  - State prerequisites: at least one group exists, but filters hide all highlighted cells
  - Action to reach state: apply filter until visible highlighted count is 0
Crop:
  - Include: highlight count reading “0 of N visible” and enough filter UI to indicate filtering is active
Alt text:
  - Highlight count showing zero visible highlighted cells due to active filters.
Caption:
  - If highlight membership exists but no highlighted cells are visible, relax filters or switch to a view where those cells are visible.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for all highlighted cells being filtered out.
:width: 100%

If highlight membership exists but no highlighted cells are visible, relax filters or switch to a view where those cells are visible.
```
