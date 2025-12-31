# Extension points overview

This page describes the **supported extension points** for contributors: how to add features without breaking performance, persistence, or UI modularity.

If you want step-by-step “do exactly this” guides, jump to:
- {doc}`16_extension_point_add_ui_module`
- {doc}`17_extension_point_add_analysis_mode`
- {doc}`18_extension_point_add_export_renderer`

## At a glance

**Audience**
- Developers extending Cellucid: read fully.

---

## Golden rules (avoid common architectural mistakes)

1) **Keep boundaries crisp**
   - UI modules own DOM and call state/viewer APIs.
   - State owns typed arrays and emits events.
   - Viewer owns WebGL resources and render loops.

2) **Never add per-frame DOM work**
   - If it needs to update frequently, it belongs in the renderer or a debounced UI update.

3) **Avoid hot-path allocations**
   - Reuse typed arrays and scratch buffers.
   - Benchmark before/after for large datasets.

4) **Think about persistence on day one**
   - Should your feature persist in sessions?
   - If not, should it be explicitly excluded so users aren’t surprised?

5) **Document edge cases and failure modes**
   - Especially for: filters ↔ highlights ↔ multiview ↔ analysis ↔ export ↔ sessions.

---

## Extension point categories

### 1) Add a UI module (sidebar feature)

When you want:
- a new accordion section
- a new sidebar control group
- a new user-facing workflow control panel

Go to:
- {doc}`16_extension_point_add_ui_module`

Key files involved:
- `cellucid/index.html` (DOM)
- `cellucid/assets/js/app/ui/core/dom-cache.js` (DOM references)
- `cellucid/assets/js/app/ui/modules/*` (module implementation)
- `cellucid/assets/js/app/ui/core/ui-coordinator.js` (wiring)

### 2) Add an analysis mode / plot / compute operation

When you want:
- a new analysis UI mode (e.g. new “tab” in Page Analysis)
- a new plot type
- a new transform/statistical test

Go to:
- {doc}`17_extension_point_add_analysis_mode`

Key files involved:
- `cellucid/assets/js/app/analysis/comparison-module.js`
- `cellucid/assets/js/app/analysis/ui/analysis-ui-manager.js`
- `cellucid/assets/js/app/analysis/plots/types/*`
- `cellucid/assets/js/app/analysis/compute/*`

### 3) Add an export renderer

When you want:
- a new SVG renderer variant
- a new PNG exporter option
- a new file format (e.g. PDF-like vector output, if supported)

Go to:
- {doc}`18_extension_point_add_export_renderer`

Key files involved:
- `cellucid/assets/js/app/ui/modules/figure-export/figure-export-engine.js`
- `cellucid/assets/js/app/ui/modules/figure-export/renderers/*`

### 4) Add a new data source

When you want:
- a new way to load datasets (beyond local-demo/local-user/remote/github/jupyter)

Start with:
- `cellucid/assets/js/data/data-source-manager.js`
- `cellucid/assets/js/data/data-source.js` (contracts + validation)

Then wire it into:
- `cellucid/assets/js/app/main.js` and/or dataset connections UI

### 5) Add session persistence coverage

When you add a feature that should persist:
- implement or extend a contributor under `cellucid/assets/js/app/session/contributors/`
- ensure you have bounds checks and dataset mismatch semantics

See:
- {doc}`10_sessions_persistence_and_serialization`

---

## “Should this be state, UI, or viewer?”

Use this decision checklist:

- If it is a **user preference or UI selection** → UI module + (maybe) `DataState` field.
- If it affects **which points are visible or colored** → `DataState` (typed arrays) + viewer update calls.
- If it affects **how pixels are drawn** (shader, buffers, per-frame loops) → rendering layer.
- If it is a **derived artifact** (analysis result, export output) → keep it off the hot path; compute on-demand.

---

Next: {doc}`16_extension_point_add_ui_module`.
