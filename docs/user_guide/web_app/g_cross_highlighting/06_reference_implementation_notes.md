# Reference (implementation notes)

:::{warning}
Cross-highlighting is **planned** and **under development**.
This page is developer-facing: it summarizes where the relevant code lives today and what the tracked plan says to implement next.
:::

## At a glance

**Audience**
- Developers and power users

**Time**
- 10–30 minutes (depending on whether you follow links into the plan docs)

**Prerequisites**
- Comfort reading JS modules and browser console logs

---

## Source-of-truth plan docs (internal)

These are the tracked design/implementation documents:
- `cellucid/markdown/CROSS-HIGHLIGHTING-FIX-PLAN.md` (root cause analysis + phased plan)
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md` (concrete checklist; includes scatterplot warning)
- `cellucid/markdown/CROSS-HIGHLIGHTING-NEW-MODULES.md` (proposed refactor into reusable modules)
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md` (render/snapshot consistency, viewId, dimensionality, memory)

---

## Where the code currently lives (today)

Current manager implementation:
- `cellucid/assets/js/app/analysis/ui/cross-highlighting.js`

Why it is currently considered “not implemented” end-to-end:
- the manager expects a viewer “embeddingView” with APIs like `highlightCells()` and `previewHighlight()`
- the viewer is not consistently wired to provide these methods (and/or the manager is not initialized at the right lifecycle point)
- plots must register themselves after rendering, but this integration is not uniformly present across analysis UIs

---

## High-level plan summary (what gets built in what order)

The plan in `cellucid/markdown/CROSS-HIGHLIGHTING-FIX-PLAN.md` breaks into phases:

### Phase 1 (critical): core wiring
- Add missing viewer APIs:
  - `highlightCells(cellIndices, options)`
  - `previewHighlight(cellIndices)`
  - `clearPreview()`
  - `clearHighlight()`
- Initialize/wire the CrossHighlightManager in the analysis module lifecycle (after the viewer exists).
- Register analysis plots after Plotly render so click/hover handlers are actually connected.

### Phase 2 (high): state integration
- Route highlights through DataState preview/highlight APIs so cross-highlighting behaves like existing selection/highlight features.
- Avoid per-click allocation of large arrays by reusing state-managed buffers.

### Phase 3 (medium): plot coverage + hover behavior
- Extend hover preview beyond barplot/histogram (with safe fallbacks).
- Keep scatterplots disabled unless they carry a sampling map (avoid wrong highlights).

### Phase 4 (medium): multi-view correctness
- Add `viewId` support so highlights target the intended view/snapshot and use per-view buffers/transparency.

---

## Proposed module refactor (planned)

The plan proposes moving cross-highlighting to a reusable UI module directory:
- `cellucid/assets/js/app/ui/modules/cross-highlight/`

With responsibilities split into:
- a manager (event + selection lifecycle),
- a viewer/state “bridge” (DataState-first integration),
- plot cell-extractor utilities (per plot type),
- small UI helpers (notifications/actions).

Reference:
- `cellucid/markdown/CROSS-HIGHLIGHTING-NEW-MODULES.md`

---

## Known critical limitation: scatterplot sampling

Scatterplots in the analysis system may sample to ~500 points for performance.
Without a sampling map, `pointIndex` does not map to the original dataset index.

Tracked mitigation:
- disable scatterplot cross-highlighting until a sampling map is threaded through plot registration.

References:
- `cellucid/markdown/CROSS-HIGHLIGHTING-IMPLEMENTATION-CHECKLIST.md`
- `cellucid/markdown/CROSS-HIGHLIGHTING-PERFORMANCE-GUIDE.md`

---

## Debugging checklist (developer quick hits)

- Confirm plots are registered (look for `[CrossHighlight] Registered plot:` logs).
- Confirm the manager has a non-null `embeddingView`.
- Confirm required viewer APIs exist (no `... is not a function` errors).
- Confirm hover preview clears on unhover (`plotly_unhover` bound).
- In multi-view, confirm `viewId` is propagated all the way to the highlight renderer path.

---

## Related docs

- {doc}`01_what_cross_highlighting_is_user_story` (UX intent)
- {doc}`05_troubleshooting_cross_highlighting` (symptoms mapped to likely wiring failures)
- {doc}`../p_developer_docs/index` (broader app architecture + debugging)
