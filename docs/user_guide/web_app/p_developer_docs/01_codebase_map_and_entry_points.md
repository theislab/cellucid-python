# Codebase map and entry points

This page answers: **“Where is the code?”** and **“What runs first?”** for the Cellucid web app.

If you are trying to run Cellucid locally first, you can skim this and jump to {doc}`02_local_development_setup`.

## At a glance

**Audience**
- Non-technical / wet lab: skim the “Big picture” diagram and “Key terms”.
- Computational users: focus on “Entry points” and “Where to look for X”.
- Developers: read the whole page; it is the fastest way to orient.

**Time**
- 10–20 minutes to get oriented

**Prerequisites**
- Ability to open files in this workspace (no special tooling required)

---

## Big picture: layers (data ↔ state ↔ UI ↔ rendering)

Cellucid is intentionally organized into layers so performance-critical code stays isolated:

- **Rendering layer** (`cellucid/assets/js/rendering/`)
  - Owns the WebGL renderer and per-frame loops.
  - Goal: “fast pixels” for millions of points.

- **Data layer** (`cellucid/assets/js/data/`)
  - Owns data sources (local demo, file picker, remote server, GitHub, Jupyter) and binary loaders.
  - Goal: “load only what you need” + predictable caching.

- **App layer** (`cellucid/assets/js/app/`)
  - **State** (`app/state/`): `DataState` and its managers (fields, filters, colors, highlights, views).
  - **UI** (`app/ui/`): DOM wiring and sidebar modules.
  - **Session + persistence** (`app/session/`, `app/state-serializer/`).
  - Goal: connect data + rendering + UI without hot-path regressions.

---

## Repo map (what lives where)

The frontend lives in `cellucid/` (this section is documented from the Python docs repo).

```
cellucid/
├── index.html                     # Single-page app shell (sidebar + canvas + scripts)
└── assets/
    ├── css/                       # Design system + themes
    ├── js/
    │   ├── app/                   # App layer (state + UI wiring + sessions + analysis)
    │   │   ├── main.js            # App entry point (bootstrap + dataset load orchestration)
    │   │   ├── state/             # DataState + managers (fields/filters/colors/highlights/views)
    │   │   ├── ui/                # UI coordinator + sidebar modules
    │   │   ├── session/           # Session bundle save/load (.cellucid-session)
    │   │   ├── state-serializer/  # Feature-scoped snapshot/restore helpers (used by session)
    │   │   ├── analysis/          # Analysis module (compute backends + plugins + UI)
    │   │   └── community-annotations/ # GitHub-backed annotation voting (auth + sync + cache)
    │   ├── data/                  # Data sources + binary loaders + h5ad/zarr adapters
    │   ├── rendering/             # WebGL viewer + shaders + overlays (smoke, highlight, etc.)
    │   └── utils/                 # Shared utilities (debug, theme, style manager, etc.)
    └── exports/                   # (Optional) demo datasets, session snapshots, etc.
```

:::{tip}
If you only read one file to understand “what runs when”, start with:
`cellucid/assets/js/app/main.js`.
It is the bootstrap orchestrator and intentionally contains *no* per-frame logic.
:::

---

## Entry points (what runs first)

### 1) `cellucid/index.html` (app shell)

`cellucid/index.html` contains:
- The `<canvas id="glcanvas">` used by the WebGL renderer.
- The sidebar markup (accordion sections, inputs, buttons).
- Early “must not throw” bootstraps (theme init, analytics init).
- The module entry point:
  - `<script type="module" src="assets/js/app/main.js"></script>`

### 2) `cellucid/assets/js/app/main.js` (bootstrap orchestrator)

`main.js` orchestrates:
- creation of the WebGL viewer (`createViewer` from `assets/js/rendering/viewer.js`)
- creation of app state (`createDataState` from `assets/js/app/state/index.js`)
- UI initialization (`initUI` from `assets/js/app/ui/core/ui-coordinator.js`)
- data source registration and initial dataset selection (`DataSourceManager`)
- dataset loading (obs/var/connectivity/embeddings) via `assets/js/data/data-loaders.js`
- session persistence (`createSessionSerializer` from `assets/js/app/session/index.js`)
- analysis initialization (comparison module)

It also exposes dev-friendly globals:
- `window._cellucidViewer`
- `window._cellucidState`
- `window._comparisonModule` (analysis)

### 3) `cellucid/assets/js/app/ui/core/ui-coordinator.js` (UI coordinator)

`initUI()` is a thin orchestrator that:
- collects DOM references once (`ui/core/dom-cache.js`)
- initializes UI modules (`ui/modules/*`)
- wires callbacks across modules without creating tight coupling

### 4) `cellucid/assets/js/app/state/core/data-state.js` (state coordinator)

`DataState` is the “source of truth” for:
- loaded field data (with bounded LRU caches),
- filters and visibility state,
- highlights and highlight pages,
- per-view contexts (live view + snapshot views),
- dimension switching and vector-field overlay availability.

The public API is assembled from manager mixins (see `app/state/managers/*`).

---

## “Where do I look for…?” (common developer questions)

| If you’re looking for… | Start with | Then |
|---|---|---|
| “How does the app decide what dataset to load?” | `cellucid/assets/js/data/data-source-manager.js` | `cellucid/assets/js/app/main.js` |
| “Where are the file format assumptions?” | `cellucid/assets/js/data/data-loaders.js` | {doc}`09_data_loading_pipeline_and_caching` |
| “Where are filters implemented?” | `cellucid/assets/js/app/state/managers/filter-manager.js` | `cellucid/assets/js/app/ui/modules/filter-controls.js` |
| “Where does color-by happen?” | `cellucid/assets/js/app/state/managers/color-manager.js` | `cellucid/assets/js/app/ui/modules/legend-renderer.js` |
| “Where do highlights live?” | `cellucid/assets/js/app/state/managers/highlight-manager.js` | `cellucid/assets/js/app/ui/modules/highlight-controls.js` |
| “Where is multi-view (‘Keep view’) handled?” | `cellucid/assets/js/app/state/managers/view-manager.js` | `cellucid/assets/js/app/ui/modules/view-controls.js` |
| “Where are sessions saved/loaded?” | `cellucid/assets/js/app/session/session-serializer.js` | {doc}`10_sessions_persistence_and_serialization` |
| “Where is figure export implemented?” | `cellucid/assets/js/app/ui/modules/figure-export/README.md` | {doc}`12_figure_export_architecture` |
| “Where is community annotation?” | `cellucid/assets/js/app/ui/modules/community-annotation-controls.js` | `cellucid/assets/js/app/community-annotations/REPO_SETUP.md` |

---

## Key terms (used consistently in code)

These show up everywhere (state, UI, analysis, sessions):

- **Dataset**: a logical unit with a stable `datasetId` (see {doc}`../b_data_loading/06_dataset_identity_why_it_matters`).
- **Field**: an `obs` categorical/continuous field, or a `var` “gene expression” field.
- **View**: the **live** view or a **snapshot** view (created by “Keep view”).
- **Visibility**: what points are currently shown after filters/outliers.
- **Highlight**: annotation overlay on top of visibility (can exist even for filtered-out points).
- **Highlight page**: independent collection of highlight groups; switching pages changes overlay state.
- **Dimension level**: 1D/2D/3D embedding selection for a view.

---

## Troubleshooting (code navigation problems)

### Symptom: “I changed a file but nothing happens”

Likely causes (ordered):
1) Browser cache is serving old modules.
2) You edited a file not loaded in your current code path (e.g. a different module).
3) Your local server is serving the wrong directory (not `cellucid/`).

How to confirm:
- In DevTools → Network, enable “Disable cache” and hard-reload.
- In DevTools → Sources, open the file and confirm your changes are visible.

Fix:
- Run a local server from `cellucid/` and reload (see {doc}`02_local_development_setup`).

### Symptom: “Module script failed to load”

Likely causes:
- The server is serving `.js` files with the wrong MIME type.
- You opened the app via `file://` instead of HTTP.

Fix:
- Use a local HTTP server (see {doc}`02_local_development_setup`).

---

Next: {doc}`02_local_development_setup`.
