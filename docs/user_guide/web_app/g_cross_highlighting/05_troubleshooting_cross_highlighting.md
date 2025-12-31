# Troubleshooting (cross-highlighting)

:::{warning}
Cross-highlighting is **planned** and **under development**.

If you are reading this because “nothing happens when I click a plot”, it may simply be that your build does not yet include the required wiring.
This page still helps you distinguish “not implemented” from “implemented but broken”.
:::

## At a glance

**Audience**
- Everyone: symptom → diagnosis → fix

**Time**
- 10–30 minutes (depending on how deep you go)

**Prerequisites**
- Ability to open browser DevTools console (recommended)
- Familiarity with highlights and filters:
  - {doc}`../f_highlighting_selection/index`
  - {doc}`../e_filtering/index`

---

## Quick triage (30 seconds)

1) **Check whether the plot type is supported**
   - Barplots/histograms are the most likely to be supported first.
   - Scatterplots are commonly sampled → cross-highlighting is expected to be disabled initially.

2) **Check whether highlights are simply invisible**
   - Temporarily disable filters in the target view.
   - Switch to a different view/snapshot where those cells might be visible.

3) **Open the console**
   - Look for logs/warnings like `[CrossHighlight] ...`.

---

## Common symptoms → likely causes

| Symptom | Likely causes (ordered) |
|---|---|
| Clicking a bar/bin does nothing | Feature not wired (manager not connected to viewer), plot not registered, plot has no cell mapping, unsupported plot type |
| Hover preview does nothing | Hover preview not implemented for that plot type, preview API missing, hover events not bound |
| The “wrong” cells highlight | Dataset mismatch, scatterplot sampling, pageData mapping bug |
| Toast appears but “Save as Page” does nothing | onCreatePage callback not connected to highlight/page system |
| Toast says N cells selected but you see 0 highlighted | Filters hide them in the target view; you are looking at the wrong panel/view |
| App becomes slow after clicks | Very large selections, repeated large allocations, or highlight updates across many snapshot views |

---

## Symptom: “Clicking does nothing”

### Likely causes
- CrossHighlightManager is not connected to the viewer (embeddingView is `null`).
- The viewer lacks required API methods (`highlightCells`, `previewHighlight`, `clearPreview`, `clearHighlight`).
- The analysis plot never registered itself with the manager (no click handler attached).
- The plot type is intentionally unsupported (e.g., sampled scatterplots).

### How to confirm
- Open DevTools console and look for debug lines like:
  - `Registered plot: <plotId> (barplot)`
- Click the plot and check for errors like `viewer.highlightCells is not a function`.

### Fix / workaround (current state)
- If the feature is not yet wired in your build: there is no user workaround.
- Use manual selection tools instead:
  - {doc}`../f_highlighting_selection/02_selection_tools_document_each_tool`

---

## Symptom: “Toast says N cells, but I see nothing highlighted”

### Likely causes
- Filters in the target view hide the highlighted cells.
- You are looking at a different view than the one receiving the highlight (multi-view/snapshots).

### How to confirm
- Disable filters in the target view (or hit “Show all” for categorical filters).
- Switch to single view temporarily and try again.

### Fix
- Disable filters temporarily and repeat the click.
- Ensure you are interacting with plots associated with the view you care about (once view targeting is implemented).

---

## Symptom: “The wrong cells highlight” (most dangerous failure mode)

### Likely causes
- Scatterplot sampling: `pointIndex` does not map to the original dataset order.
- Dataset mismatch: the analysis result was computed on a different ordering/subset than the currently loaded dataset.
- A plot reused pageData from another view or stale cache.

### How to confirm
- If it’s a scatterplot: assume sampling unless you know it carries a sampling map.
- Check whether the highlighted cells make any sense across multiple fields (e.g., click “Cluster A” bar, but highlighted cells are not in Cluster A when you color by cluster).

### Fix (planned safety behavior)
- Disable cross-highlighting for that plot type until mapping is correct.
- File a bug with a minimal repro (see below).

---

## Symptom: “Save as Page does nothing”

### Likely causes
- The toast UI exists, but the callback that creates highlight pages is not wired.

### How to confirm
- Click “Save as Page” and look for a warning like “Page creation not available”.
- Check whether a new highlight page appears in the Highlighting UI.

### Workaround
- Use manual selection tools to recreate the selection and save it.

---

## Reporting a bug (minimum repro package)

Include:
- What build/version/commit you are using (or deployment date)
- Browser + OS
- Dataset size and loading path (export folder / server / Jupyter)
- Plot type + exact plot you interacted with (bar/hist/violin/scatter)
- Whether filters were active
- Whether you were in multi-view/snapshots
- Console logs (copy/paste) including any `[CrossHighlight]` warnings/errors
- Screenshot or short screen recording

---

## Next steps / related pages

- {doc}`06_reference_implementation_notes` (developer wiring checklist and file locations)
- `cellucid/markdown/CROSS-HIGHLIGHTING-FIX-PLAN.md` (root causes + phases)
