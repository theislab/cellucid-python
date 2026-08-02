# UI modules map

This page is a **map of the UI layer**: which modules exist, what DOM they own, and how they communicate with state/rendering without creating tight coupling.

It is especially useful when:
- you want to add a new sidebar module,
- you want to find where a UI behavior lives,
- you are debugging “why did this update twice?” issues.

## At a glance

**Audience**
- Computational users: skim “UI coordinator” and “Which module owns what?”
- Developers: read fully; use the module table as a navigation map.

**Time**
- 20–40 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`

---

## UI coordinator (the orchestrator)

The UI layer is orchestrated by:
- `cellucid/assets/js/app/ui/core/ui-coordinator.js` (`initUI`)

`initUI` is intentionally thin:
- collects DOM references once (`collectDOMReferences`)
- initializes feature modules (`initX` functions)
- wires cross-module callbacks (e.g. “visibility changed → refresh legend + badges”)
- subscribes to a small set of `DataState` events

This separation exists to prevent:
- a single monolithic UI file that becomes unmaintainable
- accidental performance regressions from repeated DOM queries

---

## DOM cache (single source of truth for selectors)

Instead of every module calling `document.getElementById(...)` repeatedly, Cellucid uses:
- `cellucid/assets/js/app/ui/core/dom-cache.js`

`collectDOMReferences(document)` returns a nested object keyed by UI domains:
- `dom.dataset`
- `dom.fieldSelector`
- `dom.filter`
- `dom.highlight`
- `dom.render`
- `dom.camera`
- `dom.view`
- `dom.session`
- `dom.communityAnnotation`
- …

Design rule:
- When adding new UI elements, prefer adding them to the DOM cache rather than scattering selectors across modules.

---

## Module index (what exists today)

Most UI modules live under:
- `cellucid/assets/js/app/ui/modules/`

Most — not all — are imported and initialized by:
- `cellucid/assets/js/app/ui/core/ui-coordinator.js`

:::{important}
“Lives under `ui/modules/`” and “is initialized by the coordinator” are two
different facts, and three of the modules below are initialized by a *sibling
module* rather than by the coordinator:

- `dataset-connections.js` ← initialized by `dataset-controls.js`
- `field-selector-gene-expression.js` ← initialized by `field-selector.js`
- `field-selector-deleted-fields.js` ← initialized by `field-selector.js`

If you are tracing why a control is inert, check which of the two paths owns it
before assuming the coordinator does.

The `ui/` layer also has directories this page does not otherwise cover:
`ui/components/` (shared widgets — the coordinator initializes
`components/info-popovers.js` directly), `ui/category-builder/`,
`ui/onboarding/` (welcome modal and the global keyboard shortcuts), and
`ui/keyboard-move.js`.
:::

### Core “shell” modules

| Module | File | Owns |
|---|---|---|
| Sidebar toggle/resize | `ui/modules/sidebar-controls.js` | sidebar visibility + `--sidebar-width` |
| Stats display | `ui/modules/stats-display.js` | top stats line (dataset/field/cell counts) |

### Dataset + session modules

| Module | File | Owns |
|---|---|---|
| Dataset controls | `ui/modules/dataset-controls.js` | dataset dropdown + dataset info block |
| Dataset connections | `ui/modules/dataset-connections.js` | local-user picker, remote server connect, GitHub connect |
| Session controls | `ui/modules/session-controls.js` | save/load `.cellucid-session` buttons |

### Field/filter/highlight modules

| Module | File | Owns |
|---|---|---|
| Field selector | `ui/modules/field-selector.js` | categorical/continuous/gene selectors + copy/rename/delete |
| Gene expression helper | `ui/modules/field-selector-gene-expression.js` | gene dropdown/search behavior |
| Deleted fields panel | `ui/modules/field-selector-deleted-fields.js` | restore/purge UI |
| Filter controls | `ui/modules/filter-controls.js` | filter stack UI + count |
| Legend renderer | `ui/modules/legend-renderer.js` | legend DOM + outlier UI container |
| Highlight controls | `ui/modules/highlight-controls.js` | highlight pages + groups + selection mode tools |

### View/render/camera modules

| Module | File | Owns |
|---|---|---|
| View controls | `ui/modules/view-controls.js` | multiview (“Keep view”), layout mode, camera lock, badges |
| Dimension controls | `ui/modules/dimension-controls.js` | dimension dropdown + per-view dim cycling |
| Render controls | `ui/modules/render-controls.js` | points vs smoke, smoke sliders, connectivity toggles |
| Velocity overlay controls | `ui/modules/velocity-overlay-controls.js` | vector field overlay settings + field selection |
| Camera controls | `ui/modules/camera-controls.js` | navigation mode + camera settings panels |

### “Large feature” modules

These are bigger subsystems but still treated as modules:

| Module | File | Owns |
|---|---|---|
| Figure export | `ui/modules/figure-export/` | export UI + engine + SVG/PNG renderers |
| Community annotation | `ui/modules/community-annotation-controls.js` | auth + pull/publish + voting/moderation UI |
| Community annotation (split parts) | `ui/modules/community-annotation/` | the controls module's own submodules (browser integration, connection flow, consensus column, modal shell, auto-pull storage) |
| Community annotation dialogs | `ui/modules/community-annotation-voting-modal.js`, `ui/modules/community-annotation-modal-owner.js` | the voting dialog and its ownership/lifecycle |
| Cinematic camera | `ui/modules/cinematic-camera/` | camera-path authoring, transport bar, playback. Initialized by the coordinator; persisted as a mandatory session chunk (see {doc}`10_sessions_persistence_and_serialization`) |
| Performance benchmark | `ui/modules/benchmark/` | synthetic data generation contract + worker, the configuration-matrix harness, GL upload counters, frame recorder, camera path. The panel’s wiring lives in `main.js`; the measurement code is here and in `assets/js/dev/benchmark.js` |
| Highlight subsystem | `ui/modules/highlight/` | the four selection tools, their copy, and the highlight-page UI |
| Legend subsystem | `ui/modules/legend/` | categorical and continuous legend rendering |
| Visualization reset | `ui/modules/visualization-reset.js` | reset actions (clear filters/highlights/etc.), and the global `R` reset-camera key |
| Field interaction owner | `ui/modules/field-interaction-owner.js` | arbitrates which module owns a field interaction |
| Camera input contract | `ui/modules/camera-input-contract.js` | the validated shape of camera control input |

:::{note}
`ui/modules/cross-highlight/` exists but contains only a `.gitkeep`. There is no
cross-highlight UI module today — cross-view highlight synchronisation is
handled in the state layer (`app/state/managers/highlight-manager.js`) and in
`ui/modules/highlight/`. A reader who trusts a complete module map would go
looking for a module that is not there.
:::

:::{note}
Page Analysis (plots / DE / correlation / etc.) is initialized directly in `cellucid/assets/js/app/main.js` (not via `ui-coordinator.js`).
Its UI lives in `cellucid/assets/js/app/analysis/ui/`.
See {doc}`11_analysis_architecture`.
:::

---

## How modules communicate (state events + callbacks)

Most modules follow this pattern:

- Inputs (DOM events) call `DataState` methods (or viewer methods when appropriate).
- `DataState` emits events (`field:changed`, `visibility:changed`, etc.).
- The UI coordinator reacts to those events and triggers targeted UI updates (legend/stats/badges).

Why this matters:
- It prevents “UI module A reaches into UI module B’s internal state”.
- It makes it easier to reason about what needs to update when.

Example: visibility changes
- Filtering module updates filters → state recomputes visibility → state emits `visibility:changed`.
- UI coordinator debounces a visibility UI update and triggers:
  - filterControls render
  - highlight summary render
  - view badge render

---

## Session serialization interaction (what UI state is persisted)

The session system can capture/restore generic UI control values, but **some UI subtrees are intentionally excluded**.
Exclusion uses a DOM attribute:
- `data-state-serializer-skip="true"`

Examples:
- Figure export section is excluded (`cellucid/index.html` and the figure-export UI code).
- Benchmark section is excluded.
- Dataset connection controls are excluded (sessions assume the dataset is already loaded).
- Community annotation UI state is excluded (auth/network-driven).

Code pointers:
- `cellucid/assets/js/app/state-serializer/ui-controls.js`
- `cellucid/assets/js/app/state-serializer/README.md`

---

## Troubleshooting (UI module boundary issues)

### Symptom: “One click triggers multiple UI updates”

Likely causes:
- Two modules subscribed to the same `DataState` event and both re-render the same DOM region.
- A module is initialized twice (usually due to calling `initUI` twice on reload paths).

How to confirm:
- Set a breakpoint in the module init function and reload.
- Log subscriptions by temporarily instrumenting `EventEmitter.on`.

Fix:
- Keep DOM ownership crisp: one module “owns” a DOM subtree.
- Prefer coordinator wiring over cross-module subscriptions.

---

Next: {doc}`09_data_loading_pipeline_and_caching` (data sources) or {doc}`15_extension_points_overview` (how to add modules safely).
