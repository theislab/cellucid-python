# UX design

:::{warning}
Cross-highlighting is **planned** and **under development**.
This page documents intended UX behavior so it can be implemented consistently across plot types and views.
:::

## At a glance

**Audience**
- Wet lab / beginner: understand what you will be able to do once the feature is enabled.
- Computational users: understand how cross-highlighting interacts with filters, multi-view, and persistence.
- Developers: use this as the UX contract for implementation.

**Time**
- 10–25 minutes

**Prerequisites**
- Familiarity with:
  - {doc}`../h_analysis/index` (analysis plots)
  - {doc}`../f_highlighting_selection/index` (highlights/pages and preview highlight)
  - {doc}`../e_filtering/index` (visibility affects what you see)

**What you’ll learn**
- Hover vs click behavior (preview vs selection)
- Planned UI affordances (“Save as Page”, “Clear”, view targeting)
- How cross-highlighting should avoid interfering with manual selection tools

---

## Design goals (what “good” feels like)

Cross-highlighting should:
- make it effortless to answer “where are the cells behind this plot element?”
- feel fast and reversible (no “I clicked something and broke my state”)
- integrate with the existing highlight/page system (so users can persist the result)
- behave predictably in multi-view layouts (highlights show up where users expect)

---

## Primary interactions (planned)

### Hover → preview highlight (transient)

Intended behavior:
- Hovering a plot element previews the corresponding cells in the embedding.
- Moving the mouse away clears the preview.

Constraints:
- Preview should be **transient** and should not modify permanent highlight groups.
- Preview should not be triggered for plot types where mapping is unsafe (e.g., sampled scatterplots).

### Click → cross-highlight selection (still non-destructive)

Intended behavior:
- Clicking a plot element selects a cell set (cross-highlight selection).
- The viewer highlights those cells (typically via the preview highlight path).
- A toast/notification appears summarizing the selection.

The toast should include:
- a readable message (e.g., `12,345 cells in "Cluster A"`)
- **Save as Page** (persist into highlight system)
- **Clear** (remove the cross-highlight highlight)

---

## Planned UX flows (step-by-step)

### Flow A: “Find cells behind a bar, then save them”

1) Run an analysis that produces a barplot (e.g., category counts).
2) Click the bar for `Cluster A`.
3) Viewer highlights those cells.
4) Toast appears: `N cells in "Cluster A"`.
5) Click **Save as Page**.
6) A new highlight page/group appears in the Highlighting UI.
7) You can now:
   - switch fields, filters, and views, and
   - reuse that page/group for downstream analyses and sessions.

### Flow B: “Use hover to explore, without committing”

1) Hover different bins in a histogram.
2) Watch the highlighted region move in the embedding.
3) If you find something interesting, click to select and save.

### Flow C: “Multi-view: highlight in the correct panel”

Planned behavior:
- By default, cross-highlighting targets the **active view** (the panel the analysis is associated with).
- If you are in grid compare / snapshots, highlights should appear in the intended view, not “randomly everywhere”.

Developer note: this is why the plan introduces `viewId` propagation and per-view highlight buffers.

---

## Interaction with manual selection tools (avoid conflicts)

Manual selection tools (lasso/proximity/KNN/annotation) already use:
- a selection candidate set, and
- a preview highlight while selecting.

Cross-highlighting should not:
- wipe out an in-progress manual selection,
- “steal” the preview highlight channel in a way that confuses the user.

Planned approach (implementation-level intent):
- cross-highlighting routes through the same highlight pipeline, but should be treated as a distinct “mode” or “source” so the UI can communicate what is happening.

User-facing guideline:
:::{tip}
If you are mid-lasso and cross-highlighting appears, cancel/confirm the selection first. Cross-highlighting is meant to be used when you are not actively dragging a selection gesture.
:::

---

## Error/empty states (what users should see)

When cross-highlighting cannot map safely, the UX should be explicit:

- If a plot type is unsupported: show a small warning/toast like “Cross-highlighting unavailable for scatterplots (sampled)”.
- If the plot has no cell mapping: “This plot has no per-cell mapping; cannot highlight cells”.
- If filters hide all highlighted cells: toast still shows count, but optionally add “(0 visible in this view; check filters)”.

---

## Screenshot placeholders (capture later)

<!-- SCREENSHOT PLACEHOLDER
ID: cross-highlighting-toast-save-as-page
Suggested filename: cross-highlighting/01_toast-save-as-page.png
Where it appears: User Guide → Web App → Cross-highlighting → UX design
Capture:
  - UI location: analysis window with a barplot or histogram + the embedding canvas
  - State prerequisites: dataset loaded; analysis plot visible; cross-highlighting enabled
  - Action to reach state: click a bar/bin so a toast appears with “Save as Page” and “Clear”
Crop:
  - Include: plot element clicked, highlighted region in embedding, toast with actions
Redact:
  - Remove: private dataset identifiers, local paths
Alt text:
  - A selection toast appears after clicking an analysis plot element, with actions to save the selection as a page or clear it.
Caption:
  - Cross-highlighting is designed to be reversible: click a plot element to highlight cells, then choose whether to persist the result as a highlight page.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for cross-highlighting toast UX.
:width: 100%

Cross-highlighting toast after clicking a plot element, showing the planned “Save as Page” action.
```

---

## Next steps / related pages

- {doc}`04_performance_correctness_notes`
- {doc}`05_troubleshooting_cross_highlighting`
- {doc}`../f_highlighting_selection/03_highlight_ui`
