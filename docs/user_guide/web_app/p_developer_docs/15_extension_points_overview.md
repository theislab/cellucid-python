# Extension points overview

This page describes the **supported extension points** for contributors: how to add features without breaking performance, persistence, or UI modularity.

If you want step-by-step “do exactly this” guides, jump to:
- {doc}`16_extension_point_add_ui_module`
- {doc}`17_extension_point_add_analysis_mode`
- {doc}`18_extension_point_add_export_renderer`
- {doc}`29_action_layer_domains_and_history` (“How to add an action”, “How to
  add a piece of state”)

## At a glance

**Audience**
- Developers extending Cellucid: read fully.

---

## Golden rules (avoid common architectural mistakes)

1) **Keep boundaries crisp**
   - UI modules own DOM and dispatch actions; an action's implementation
     calls the state/viewer APIs. A control calls no mutating state method.
   - State owns typed arrays and emits events.
   - Viewer owns WebGL resources and render loops.

2) **Never add per-frame DOM work**
   - If it needs to update frequently, it belongs in the renderer or a debounced UI update.

3) **Avoid hot-path allocations**
   - Reuse typed arrays and scratch buffers.
   - Benchmark before/after for large datasets.

4) **Think about the history and persistence on day one**
   - Every change a person can make is an action, so it is undoable and
     replays; what it changes belongs to a state domain, so a step, an undo
     and a session hold it.
   - If something should not be kept, exclude it explicitly so users aren’t
     surprised.

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
- a new analysis UI mode (e.g. new “tab” in the Analysis panel)
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

Pick the smallest of three extension points:

- **A new host for links** (a repository or cloud drive): one provider module
  under `cellucid/assets/js/data/sources/providers/` plus one line in
  `DEFAULT_PROVIDERS` (`data/sources/registry.js`, whose header states the
  provider contract). A provider only says where the bytes are; it never
  parses a format.
- **A new file format**: one entry in `cellucid/assets/js/data/input-formats.js`
  (regenerate `tests/fixtures/input-formats.json`, which Python and R read),
  plus a layout view under `data/anndata/views/` for a format read in place,
  or a transcoder under `data/anndata/transcode/` for one read once.
- **A new kind of source** (beyond local-demo/local-user/remote/link/jupyter):
  start with `cellucid/assets/js/data/data-source-manager.js` and
  `cellucid/assets/js/data/data-source.js` (contracts + validation), then wire
  it into `cellucid/assets/js/app/main.js` and the dataset connections UI.

### 5) Add an action

When you want:
- a new change a person, an agent, or a Python or R host can make (it becomes
  one undoable history step, an agent tool, a Python method and an R function)
- a new read an agent or a host can ask for

Go to:
- {doc}`29_action_layer_domains_and_history` (“How to add an action”)

Key files involved:
- `cellucid/assets/js/app/actions/definitions/<area>.js` (the definition)
- the module that owns the state (`getActions().provide(type, handler)`)
- `cellucid/scripts/sync-shared-files.mjs` (it runs
  `generate-action-catalog.mjs` and writes the packages' copies),
  `cellucid-python/scripts/generate_actions.py`,
  `cellucid-r/tools/generate-actions.R`,
  `cellucid/scripts/generate-relay-bundle.mjs` (regenerate every surface)

### 6) Add a piece of persisted state

When you add state that should be undone and saved:
- make it a state domain's value, or part of one
  (`cellucid/assets/js/app/state/domains/`): the history, the autosave and the
  session file hold it with no encoder of their own
- every action that can change it names the domain in `touches`

See:
- {doc}`29_action_layer_domains_and_history` (“How to add a piece of state”)
- {doc}`10_sessions_persistence_and_serialization`

---

## “Should this be state, UI, or viewer?”

Use this decision checklist:

- If it is a **user preference or UI selection** → UI module + (maybe) `DataState` field; a stable sidebar control is a setting (`settings` domain) with no code.
- If it affects **which points are visible or colored** → `DataState` (typed arrays) + viewer update calls.
- If it affects **how pixels are drawn** (shader, buffers, per-frame loops) → rendering layer.
- If it is a **derived artifact** (analysis result, export output) → keep it off the hot path; compute on-demand.

---

Next: {doc}`16_extension_point_add_ui_module`.
