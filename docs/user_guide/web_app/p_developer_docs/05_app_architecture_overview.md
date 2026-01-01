# App architecture overview

This page is the **big-picture mental model** for the Cellucid web app: how the app boots, how data moves through the system, and where performance-sensitive boundaries are.

If you are new to the codebase, start here before diving into any single module.

## At a glance

**Audience**
- Wet lab / non-technical: read “Big picture” + “What is state?” + “How to report a bug”.
- Computational users: read “Boot sequence” + “Data sources” + “State and persistence”.
- Developers: read the whole page, then continue to {doc}`06_state_datastate_and_events`.

**Time**
- 20–40 minutes

**Prerequisites**
- None (links point to code files in this workspace)

---

## Big picture (boot → state → UI → rendering)

```
cellucid/index.html
  ├─ defines: sidebar DOM + <canvas id="glcanvas">
  ├─ loads:   early theme/analytics bootstraps
  └─ runs:    assets/js/app/main.js  (ES module)

main.js (orchestrator)
  ├─ createViewer()     → rendering/viewer.js (WebGL2 renderer)
  ├─ createDataState()  → app/state/core/data-state.js (typed arrays + state managers)
  ├─ initUI()           → app/ui/core/ui-coordinator.js (DOM + modules)
  ├─ data sources       → data/data-source-manager.js (local/remote/GitHub/Jupyter)
  ├─ data loaders       → data/data-loaders.js (binary + manifests + h5ad/zarr adapters)
  ├─ sessions           → app/session/session-serializer.js (.cellucid-session bundles)
  └─ analysis           → app/analysis/* (compute backends + plots + UI)

DataState (state coordinator)
  ├─ emits events: visibility/field/highlight/page/dimension changes
  ├─ owns: typed arrays + per-view contexts
  └─ calls viewer: updateColors/updateTransparency/updatePositions/etc
```

Design goal: **keep “hot” work** (per-frame render and per-point math) inside the renderer and state managers, not in `main.js` or UI modules.

---

## Boot sequence (what happens on page load)

The high-level startup flow is:

1) `cellucid/index.html` loads `assets/js/app/main.js` as an ES module.
2) `main.js` creates a WebGL2 viewer: `createViewer({ canvas, labelLayer, viewTitleLayer, sidebar })`.
3) `main.js` creates the app state: `createDataState({ viewer, labelLayer })`.
4) `main.js` initializes notifications and analytics, then sets up the `DataSourceManager`.
5) The app decides what to load first:
   - URL params (`?remote=…`, `?github=…`, `?dataset=…`) may override the default demo.
   - Jupyter context (iframe) may register a Jupyter data source and auto-load.
6) The app loads dataset metadata and core buffers:
   - `obs_manifest.json` (field list + metadata)
   - optional `var_manifest.json` (gene expression availability)
   - points buffers (embedding positions, colors, outlier quantiles, etc.)
   - optional connectivity (KNN graph edges)
7) The app initializes the UI coordinator (`initUI`) which wires the sidebar modules to state/viewer.
8) The session serializer may auto-restore a “latest session” from the dataset exports directory (if configured).

---

## Data sources (how the frontend “gets data”)

Cellucid supports multiple data sources through a single coordinator:

- `cellucid/assets/js/data/data-source-manager.js` (`DataSourceManager`)
  - Registers sources and manages dataset switching.
  - Produces a **dataset base URL** (`baseUrl`) used by the low-level loaders.

Common sources:
- **local-demo**: datasets loaded from an exports base URL (configured via `<meta name="cellucid-exports-base-url" ...>` or `?exportsBaseUrl=...`; in production this typically points at a separate `cellucid-datasets` host)
- **local-user**: browser file picker source (user selects folder/h5ad/zarr)
- **remote**: connects to a `cellucid-python` server (lazy loading; best for large h5ad/zarr)
- **github-repo**: reads exports from a GitHub repo/path (sharing; no server)
- **jupyter**: a special source used when running inside an embedded notebook viewer

Deep dive: {doc}`09_data_loading_pipeline_and_caching`.

---

## What is “state” in Cellucid?

Cellucid has *many* interacting features. The app treats state as first-class:

### `DataState` is the core coordinator

`DataState` (in `cellucid/assets/js/app/state/core/data-state.js`) owns:
- core typed arrays (positions, colors, visibility/transparency)
- field registries (rename/delete/user-defined)
- highlight pages/groups and per-point highlight overlay arrays
- per-view contexts (live view + snapshot views)
- dimension switching state and vector-field availability
- caches for loaded obs/var field data (bounded LRU)

`DataState` emits events so UI/analysis modules can stay decoupled:
- `visibility:changed`
- `field:changed`
- `highlight:changed`
- `page:changed`
- `dimension:changed`
- `vectorFields:changed`

Deep dive: {doc}`06_state_datastate_and_events`.

### Viewer state is separate

The WebGL viewer (`cellucid/assets/js/rendering/viewer.js`) owns GPU resources and the render loop.
It is updated by the state layer using explicit methods like:
- `viewer.setData(...)`
- `viewer.updateColors(...)`
- `viewer.updateTransparency(...)`
- `viewer.updatePositions(...)`

Deep dive: {doc}`07_rendering_pipeline_webgl_and_performance_notes`.

---

## UI architecture (how the sidebar is wired)

The UI is intentionally modular:

- `cellucid/assets/js/app/ui/core/dom-cache.js` collects DOM element references once.
- `cellucid/assets/js/app/ui/core/ui-coordinator.js` initializes modules and wires callbacks.
- `cellucid/assets/js/app/ui/modules/*` contains feature modules (filters, highlights, sessions, export, etc.).

Deep dive: {doc}`08_ui_modules_map`.

---

## Persistence and reproducibility

Cellucid has multiple persistence-like mechanisms; they are different and must not be conflated:

1) **Session bundle** (`.cellucid-session`)
   - A downloadable, shareable file that can restore UI + state.
   - Code: `cellucid/assets/js/app/session/` and `cellucid/assets/js/app/state-serializer/`.

2) **Browser storage** (preferences/caches)
   - Theme, debug flags, some community-annotation caches, etc.
   - Documented in: {doc}`04_configuration_env_vars_and_feature_flags` and {doc}`../o_accessibility_privacy_security/02_privacy_model`.

3) **URL state** (deep links)
   - `?dataset=…`, `?remote=…`, `?github=…`, `?annotations=…`.
   - Code: `cellucid/assets/js/app/url-state.js`.

Deep dive: {doc}`10_sessions_persistence_and_serialization`.

---

## Performance invariants (rules contributors should treat as “hard”)

These are the most common sources of accidental regressions:

- **No per-frame DOM work**: the render loop should not query or update DOM every frame.
- **No hot-path allocations**: avoid allocating arrays/objects in per-point loops; reuse typed-array scratch buffers.
- **State → viewer updates should be explicit**: prefer `viewer.update*` methods rather than “implicit” coupling.
- **Batch expensive recomputations**: when restoring filters or applying many changes, use `DataState.beginBatch()` / `endBatch()` when available.

---

## Troubleshooting: how to report a bug in a way that can be fixed

If you are reporting an issue (even as a non-developer), include:

1) **Environment**
   - Browser + OS
   - Hosted (`https://www.cellucid.com`) vs local app (`http://localhost:8000`) vs Jupyter iframe

2) **Dataset identity**
   - dataset id (and where it came from: local-user folder, remote server, GitHub repo/path)

3) **Exact steps**
   - “click X, then select Y, then move slider Z”

4) **Expected vs actual**
   - What did you think would happen?
   - What happened instead?

5) **Console + network clues**
   - Turn on debug: `localStorage.setItem('CELLUCID_DEBUG','true'); location.reload();`
   - Copy the first error stack trace and any failed network requests

Deep playbook: {doc}`13_debugging_playbook`.

---

Next: {doc}`06_state_datastate_and_events`.
