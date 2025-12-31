# Performance / correctness notes

:::{warning}
Cross-highlighting is **planned** and **under development**.
This page summarizes the key performance and correctness constraints from the implementation plan so the feature can be delivered without breaking rendering or highlighting semantics.
:::

## At a glance

**Audience**
- Computational users: understand what might get slow and how to work safely on large datasets.
- Developers: implement cross-highlighting without regressing FPS/memory or breaking multi-view highlighting.

**Time**
- 10–20 minutes

**Prerequisites**
- Familiarity with highlights + visibility:
  - {doc}`../f_highlighting_selection/01_highlight_mental_model`
  - {doc}`../e_filtering/index`

**What you’ll learn**
- The performance “hot paths” that must not be touched
- Correctness rules for multi-view/snapshots and dimensionality
- The main known footgun: scatterplot sampling

---

## Performance constraints (what must remain fast)

Cross-highlighting must be **event-driven**.
It should run on:
- plot hover/click handlers,
- selection confirmation actions,
- (optionally) debounced hover previews.

It must **not** run in the render loop.

Reference for the detailed performance rationale:
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md`

### Avoid large allocations

On big datasets, allocating a fresh `Uint8Array(pointCount)` per click can be expensive (tens of MB).

Planned approach:
- route cross-highlighting through the existing DataState preview highlight pipeline, which reuses buffers and is designed for frequent updates.

### Avoid O(n_cells) scans on interactions

Cross-highlighting should use precomputed arrays (like `pageData.cellIndices`), not scan the entire dataset in JS on each event.

If a plot element corresponds to “all cells in category X”, cross-highlighting should:
- identify those cells from the analysis result mapping (preferred), or
- use an indexed representation (if introduced later),
not do a full per-event scan.

---

## Correctness constraints (avoid “wrong highlights”)

### 1) Cell indices are the only stable identity

Do not try to map via plotted coordinates, hover positions, or ordering assumptions.
Cross-highlighting should operate on **cell indices**, and let the viewer handle coordinates for the current dimension.

### 2) Multi-view requires view targeting (viewId)

In multi-view (snapshots/grid compare), each view can have:
- different positions (different embedding),
- different visibility/transparency (different filters),
- different render buffers.

Planned requirement:
- cross-highlighting passes `viewId` through the highlight pipeline so the correct view buffer is used.

Reference:
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (render vs renderWithSnapshot consistency + per-view buffers)

### 3) “Nothing visible” can be correct (filters)

A cross-highlight selection can include many cells that are currently filtered out in the target view.

Planned UX mitigation:
- include “0 visible in this view; check filters” style hints in notifications.

### 4) Scatterplot sampling must not produce wrong mappings

Scatterplots are commonly sampled (e.g., 500 points max). Without a sampling map, `pointIndex` does not map to dataset index.

Planned safety rule:
- disable scatterplot cross-highlighting unless a sampling map is present.

Reference:
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md`
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (sampling map pattern)

---

## Practical performance tips (for users, once implemented)

- If hover preview feels laggy: move the cursor slower, or disable hover preview (if the UI adds a toggle).
- If clicking a bar selects “millions of cells”: expect a slower highlight update; consider filtering down first, then using cross-highlighting for smaller sets.
- If you are in grid compare with many snapshots: highlights must update multiple view buffers; temporarily reduce view count when exploring huge selections.

---

## What to capture when reporting a performance/correctness bug

- Dataset approximate size (`n_cells`, `n_genes`, #categories in the active field)
- What plot type you interacted with (bar/hist/violin/scatter)
- Whether you were in single view or grid compare / snapshots
- Whether filters were active in the target view
- Browser console logs (look for `[CrossHighlight]` warnings)
- A screen recording if the bug is “laggy hover/click”

---

## Next steps / related pages

- {doc}`05_troubleshooting_cross_highlighting`
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (developer deep dive)
