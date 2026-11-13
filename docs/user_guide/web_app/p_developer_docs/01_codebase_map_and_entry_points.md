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
  - Owns the data sources (sample catalog, files and folders picked or dropped,
    a `cellucid serve` server, links through the providers in `data/sources/`,
    the notebook/R viewer), the input-format registry (`data/input-formats.js`)
    and the AnnData reader worker.
  - The container readers it uses live in `cellucid/assets/js/formats/`
    (HDF5, Zarr, ZIP, R serialization, MatrixMarket text, byte I/O).
  - Goal: “load only what you need” + predictable caching.

- **App layer** (`cellucid/assets/js/app/`)
  - **State** (`app/state/`): `DataState` and its managers (fields, filters,
    colors, highlights, views); the owners of the state domains
    (`app/state/domains/`).
  - **Actions** (`app/actions/`): the one way state changes. Every change and
    read is an action of one catalog; controls, agents, Python and R hosts
    and the history all dispatch them ({doc}`29_action_layer_domains_and_history`).
  - **History** (`app/history/`): trackers of step trees, undo, redo,
    branches, the autosave and saved states, the history in words.
  - **UI** (`app/ui/`): DOM wiring and sidebar modules.
  - **Sessions** (`app/session/`): the `.cellucid-session` file, session
    format 4, which holds the history ({doc}`10_sessions_persistence_and_serialization`).
  - Goal: connect data + rendering + UI without hot-path regressions.

---

## Repo map (what lives where)

The frontend lives in `cellucid/` (this section is documented from the Python docs repo).

```
cellucid/
├── index.html                     # Single-page app shell (sidebar + canvas + scripts)
├── relay/                         # The Cellucid relay: a Cloudflare Worker that makes
│                                  #   Cellucid a connector in Claude and ChatGPT
│                                  #   (remote MCP + OAuth; see 32_agents_in_the_tab)
└── assets/
    ├── css/                       # Design system + themes
    ├── external/                  # Vendored third-party libraries, served from
    │                              #   this origin (plotly, gl-matrix, and
    │                              #   fflate-0.8.2/, the ES module the readers'
    │                              #   inflater imports, unchanged, license and
    │                              #   digest pinned: tests/wasm-codecs-contract).
    │                              #   Nothing is loaded from a CDN; the page CSP
    │                              #   would block it. plotly.min.js carries one
    │                              #   local, test-pinned patch (violin density;
    │                              #   see its header) to re-apply on upgrade.
    ├── wasm/                      # codecs.wasm: zstd, lz4 and c-blosc decoders,
    │                              #   built from pinned release tarballs by
    │                              #   scripts/build-wasm-codecs.mjs (sources in
    │                              #   scripts/wasm-codecs/); LICENSES/ beside it.
    ├── fonts/  img/               # Static assets
    └── js/
        ├── analytics/             # Google Analytics helpers. Inert unless the
        │                          #   page is served from a production hostname
        │                          #   (see ui/core/ga-init.js). Reached only by
        │                          #   the dynamic import in app/optional-analytics.js,
        │                          #   so a blocker refusing it cannot stop the app.
        │                          #   schema.js declares every event and parameter
        │                          #   (the GA4 event dictionary and admin setup).
        ├── app/                   # App layer (state + actions + history + UI + sessions + analysis)
        │   ├── main.js            # App entry point: the start-up sequence, in order
        │   ├── dataset/           # The dataset shown: publication.js (stage, commit,
        │   │                      #   restore, settle the interface), reload.js (one
        │   │                      #   dataset change as a transaction), published-state.js
        │   │                      #   (the catalog's starting view)
        │   ├── actions/           # The action catalog (catalog.js, schema.js,
        │   │                      #   definitions/<area>.js, definitions/params.js),
        │   │                      #   the dispatcher, the state domain list (domains.js),
        │   │                      #   a control's dispatch (control-dispatch.js),
        │   │                      #   outcomes as JSON (outcome.js), and names to
        │   │                      #   current objects (targets/: fields.js, pages.js,
        │   │                      #   views.js, names.js)
        │   ├── history/           # Trackers and step trees (tree.js), the value pool
        │   │                      #   (value-pool.js, value-parts.js), the history
        │   │                      #   (history.js), a dataset's history (dataset-history.js),
        │   │                      #   autosave and saved states (autosave.js,
        │   │                      #   autosave-store.js), current.js, exports/
        │   ├── state/             # DataState (state/core) + managers
        │   │                      #   (fields/filters/colors/highlights/views);
        │   │                      #   state/domains/: one owner per state domain
        │   │                      #   (cells.js ... layout.js; cell-positions.js,
        │   │                      #   fields.js, view-list.js, fingerprint.js);
        │   │                      #   state/actions/: the state's own action handlers
        │   ├── ui/                # UI coordinator (ui/core) + sidebar modules (ui/modules;
        │   │                      #   connectivity-controls.js: the edges and their
        │   │                      #   controls; benchmark/benchmark-controls.js: the
        │   │                      #   Benchmark panel; smoke-density.js: the smoke
        │   │                      #   volume rebuilt;
        │   │                      #   ui/modules/state/: the State panel;
        │   │                      #   ui/modules/my-datasets/: My datasets)
        │   ├── session/           # Session files (.cellucid-session, format 4: the history)
        │   ├── state-serializer/  # The sidebar's control inventory (ui-controls.js,
        │   │                      #   controlDomainOf: which domain holds a control)
        │   ├── library/           # My datasets (my-datasets.js, IndexedDB cellucid-library)
        │   │                      #   and kept file handles (local-file-handles.js)
        │   ├── viewer-host/       # A Python or R host's side: host-view.js (the view's
        │   │                      #   end, hover and click), host-channel.js (session
        │   │                      #   bundles, commands, results, history events),
        │   │                      #   its actions (command-handler.js), pointer
        │   │                      #   delivery, results
        │   ├── agent/             # Agents in this tab: what they are told
        │   │                      #   (instructions.js), one tool call and the
        │   │                      #   relay's tab end (bridge/: tool-call.js,
        │   │                      #   agent-calls.js, data-guard.js, mcp-tools.js),
        │   │                      #   the Connect panel (connect/), permissions,
        │   │                      #   confirm cards and activity lines (presence/),
        │   │                      #   the chat bar (chat/) and its model providers
        │   │                      #   (providers/), chats and the history
        │   │                      #   (history-link.js), the per-dataset switch
        │   │                      #   "names and summaries only" (data-scope.js)
        │   ├── analysis/          # Analysis module (compute backends + plugins + UI)
        │   ├── community-annotations/ # GitHub-backed annotation voting (auth + sync + cache)
        │   ├── dockable-accordions/   # Sidebar panel docking/undocking
        │   ├── notification-center/   # The toast/notification surface
        │   └── registries/            # Shared registration points
        ├── data/                  # Data sources, input-format registry, link providers
        │                          #   (data/sources/), AnnData reader worker + views;
        │                          #   data/anndata/write/: the AnnData writers (the
        │                          #   element layer the transcoders share, and a copy
        │                          #   of the open dataset as a .zarr.zip);
        │                          #   data/anndata/reader.js: the AnnData reader's state
        │                          #   and lifecycle, its method groups in
        │                          #   data/anndata/reader/ and the functions they share
        │                          #   beside it (uns.js, probe.js, cell-maps.js, ...)
        ├── formats/               # Container readers: hdf5/, zarr/, zip/, rds/,
        │                          #   text/ (MatrixMarket), io/ (byte sources, decode pool);
        │                          #   writers: zarr/writer.js (Zarr v2), zip/writer.js
        │                          #   (stored ZIP, ZIP64), hdf5/writer.js
        ├── dev/                   # benchmark.js — the measurement harness behind
        │                          #   the Benchmark panel. Lazily
        │                          #   imported; never on the product render path.
        ├── rendering/             # WebGL viewer + shaders + overlays (smoke, highlight, etc.)
        ├── utils/                 # Shared utilities (debug, theme, style manager, exact records; cancellation abort.js, failures errors.js and error-wire.js)
        └── workers/               # worker-script.js: where each worker starts from
                                   #   (its source on this machine, its bundle elsewhere)
    └── workers/                   # Generated, never edited: one bundle per worker
                                   #   (<name>.<sha256>.js + source map) and
                                   #   scripts.js, the URL of each; written by
                                   #   scripts/generate-web-asset-inventory.mjs
```

:::{note}
Sample datasets are **not** in this repository — they live in
`cellucid-datasets` as a dataset catalog (AnnData files plus `datasets.json`)
and are reached through the sample catalog URL (the
`cellucid-sample-catalog` meta element in `index.html`). See {doc}`09_data_loading_pipeline_and_caching`.
:::

:::{tip}
If you only read one file to understand “what runs when”, start with:
`cellucid/assets/js/app/main.js`.
It is the bootstrap orchestrator and intentionally contains *no* per-frame logic.
Its parts live in their own modules (below).
:::

---

## Entry points (what runs first)

### 1) `cellucid/index.html` (app shell)

`cellucid/index.html` contains:
- A restrictive content security policy. `default-src 'self'`; the only external
  script origins allowed are Google’s tag manager and, for Google Drive files
  only, Google's sign-in and picker (`accounts.google.com`,
  `apis.google.com`). Adding a CDN dependency means editing this policy, which
  is the intended friction. It grants `'wasm-unsafe-eval'`, so the readers can
  compile the WebAssembly decoders (`assets/wasm/codecs.wasm`), and no
  JavaScript evaluation. `blob:` is allowed for images (`img-src`: a figure
  export's PNG and a picture drawn for the report or an agent), for
  `connect-src` and for workers (`worker-src`); `connect-src` also allows
  `https:`, `wss:`, `http:` and `ws:`, since a dataset, a provider's API or
  the relay can be at any address.
- The `<canvas id="glcanvas">` used by the WebGL renderer.
- The sidebar markup (accordion sections, inputs, buttons). Data-source,
  renderer and benchmark controls are authored **disabled** and are enabled
  later by `ui/core/deferred-control-readiness.js`, so a control can never be
  clicked before something is listening for it.
- Early “must not throw” bootstraps, classic scripts (not modules) that run
  before the module graph: `ui/core/preconnect-init.js`, `ui/core/ga-init.js`,
  `ui/core/theme-init.js`, `ui/core/viewport-init.js`,
  `ui/core/code-load-watch.js`, `ui/core/file-drop-guard.js` and
  `ui/core/press-guard.js` (a double click or a held Enter presses a button
  that commits work or starts and stops it once: `.btn-primary` and
  `[data-single-press]`, every button that turns into its own Stop, Cancel or
  Run again).
- The import map that gives every module the URL naming its bytes
  (`CELLUCID_MODULE_VERSIONS`, written by
  `scripts/generate-web-asset-inventory.mjs`; see
  {doc}`03_build_run_and_deployment`, checklist item 5).
- The module entry point:
  - `<script type="module" src="assets/js/app/main.js?v=…"></script>`

### 2) `cellucid/assets/js/app/main.js` (bootstrap orchestrator)

`main.js` runs the start-up sequence, in order. The parts that are not the
sequence itself are factories it calls (`app/dataset/`, `app/viewer-host/`,
`ui/modules/connectivity-controls.js`, `ui/modules/benchmark/benchmark-controls.js`,
`ui/modules/smoke-density.js`; {doc}`05_app_architecture_overview` says what each
owns). It orchestrates:
- creation of the WebGL viewer (`createViewer` from `assets/js/rendering/viewer.js`)
- creation of app state (`createDataState` from `assets/js/app/state/index.js`)
- UI initialization (`initUI` from `assets/js/app/ui/core/ui-coordinator.js`)
- data source registration and initial dataset selection (`DataSourceManager`)
- dataset loading (obs/var/connectivity/embeddings) via `assets/js/data/data-loaders.js`, as a dataset change of `app/dataset/reload.js` published by `app/dataset/publication.js`
- session files (`createSessionSerializer` from `assets/js/app/session/index.js`)
- each dataset's history (`createDatasetHistories` from
  `assets/js/app/history/dataset-history.js`) and the owners of the state
  domains (`installStateDomains` from `assets/js/app/state/domains/index.js`)
- My datasets (`recordDatasetOpened` from `assets/js/app/library/my-datasets.js`,
  once a dataset is published)
- analysis initialization (comparison module)

It also exposes dev-friendly globals:
- `window._cellucidViewer`
- `window._cellucidState`
- `window._cellucidHistories` (the dataset histories: `current()` is the
  history of the dataset shown)
- `window._comparisonModule` (analysis)

### 3) `cellucid/assets/js/app/ui/core/ui-coordinator.js` (UI coordinator)

`initUI()` (about 1,200 lines of wiring) does three things:
- collects DOM references once (`ui/core/dom-cache.js`)
- initializes UI modules (`ui/modules/*`)
- wires callbacks across modules without creating tight coupling

New features go into a module of their own under `ui/modules/`, not into the
coordinator.

### 4) `cellucid/assets/js/app/state/core/data-state.js` (state coordinator)

`DataState` is the “source of truth” for:
- loaded field data (with bounded LRU caches),
- filters and visibility state,
- highlights and highlight pages,
- per-view contexts (live view + snapshot views),
- the embedding each view shows, the data mapping, and vector-field overlay
  availability.

The public API is assembled from manager mixins (see `app/state/managers/*`).

---

## “Where do I look for…?” (common developer questions)

Each folder of `cellucid/assets/js/` that is edited often has a short
`README.md` (what it owns, entry points, invariants, tests):
`formats/`, `formats/io/`, `data/`, `data/anndata/`, `data/sources/`,
`rendering/`, `utils/`, `app/`, `app/actions/`, `app/history/`,
`app/history/exports/`, `app/library/`, `app/session/`, `app/state/domains/`,
`app/state-serializer/`, `app/agent/` (and its `bridge/`, `chat/`,
`presence/`), `app/ui/modules/state/`,
`app/ui/modules/my-datasets/`. `assets/js/ARCHITECTURE.md` holds the
layers, the glossary (one word per concept) and the full "where to change X"
index; {doc}`20_read_path_architecture` follows one read from a source to
the pixels, naming the module that owns each step; the decision records
({doc}`21_decisions_decoders_and_codecs` to {doc}`28_decisions_rendering_layers`
and {doc}`33_decisions_actions_domains_and_history` to
{doc}`38_decisions_testing_and_hosted_runners`, listed on {doc}`index`)
say why each subsystem is built the way it is.

| If you’re looking for… | Start with | Then |
|---|---|---|
| “How does the app decide what dataset to load?” | `cellucid/assets/js/data/data-source-manager.js` | `cellucid/assets/js/app/main.js` (the start-up steps), `cellucid/assets/js/app/dataset/reload.js` (the change itself) |
| “Which files and folders does the app open?” | `cellucid/assets/js/data/input-formats.js` (the format registry) | `cellucid/assets/js/data/anndata/views/` (HDF5/Zarr layouts read in place), `cellucid/assets/js/data/anndata/transcode/` (R objects, MatrixMarket), {doc}`09_data_loading_pipeline_and_caching` |
| “How are reads turned into requests, and how is the link measured?” | `cellucid/assets/js/formats/io/byte-reader.js` | `formats/io/link.js`, `formats/io/byte-source.js`, `formats/io/hosts.js`, {doc}`20_read_path_architecture` |
| “Which cells does a dataset open with (subsample, auto-fit, block size)?” | `cellucid/assets/js/data/anndata/subsample-plan.js` | `data/subsample.js`, `data/anndata/views/cell-subset.js`, {doc}`19_subsampling_specification` |
| “How is a gene read, priced and sampled?” | `cellucid/assets/js/data/anndata/matrix-access.js` | `data/anndata/subsample-plan.js` (`geneSampleOf`), `data/sparse-column-index.js`, `data/anndata/prefix-column-index.js`, `app/state/managers/field/gene-sample.js` |
| “What is kept in this browser across visits?” | `cellucid/assets/js/data/anndata/transcode/space.js` | `data/anndata/transcode/kept-format.js`, `data/anndata/transcode/kept-journal.js`, `formats/io/kept-parts.js`, `data/anndata/kept-column-index.js`, `app/ui/modules/kept-stores-dialog.js` |
| “How much memory may the app use?” | `cellucid/assets/js/data/memory-budget.js` | `data/memory-census.js`, `app/page-memory.js` |
| “Where are the codecs and decode workers?” | `cellucid/assets/js/formats/io/decode-pool.js` | `formats/io/decoder-choice.js`, `formats/io/wasm-codecs.js` |
| “Where is a link provider (Zenodo, GEO, CELLxGENE, …)?” | `cellucid/assets/js/data/sources/providers/` | `data/sources/registry.js` |
| “How does the Python or R viewer talk to the page?” | `cellucid/assets/js/data/viewer-protocol.js` | `data/viewer-channel.js`, `data/viewer-host-source.js`, `app/viewer-host/` (`command-handler.js`: a host's `action` command dispatched as an action) |
| “Where is an action defined, and what does it change?” | `cellucid/assets/js/app/actions/definitions/<area>.js` | `app/actions/catalog.js`, `app/actions/dispatcher.js`, the module that provides it (`app/actions/README.md`), {doc}`29_action_layer_domains_and_history` |
| “How does a step name a field, a page or a view?” | `cellucid/assets/js/app/actions/targets/` | `app/actions/definitions/params.js` |
| “What does a step, an undo or a session hold of a piece of state?” | `cellucid/assets/js/app/actions/domains.js` | `app/state/domains/<domain>.js`, `app/state/domains/README.md` |
| “Where are undo, redo, branches and the autosave?” | `cellucid/assets/js/app/history/history.js` | `app/history/tree.js`, `app/history/autosave.js`, `app/ui/modules/state/` (the State panel) |
| “Where is My datasets?” | `cellucid/assets/js/app/library/my-datasets.js` | `app/library/local-file-handles.js`, `app/ui/modules/my-datasets/` |
| “Where are filters implemented?” | `cellucid/assets/js/app/state/managers/filter-manager.js` | `cellucid/assets/js/app/ui/modules/filter-controls.js` |
| “Where does color-by happen?” | `cellucid/assets/js/app/state/managers/color-manager.js` | `cellucid/assets/js/app/ui/modules/legend-renderer.js` |
| “How are unread, zero and colored cells drawn?” | `cellucid/assets/js/rendering/alpha-visibility.js` | `rendering/read-cell-layer.js`, `rendering/high-perf-renderer.js` |
| “Where do highlights live?” | `cellucid/assets/js/app/state/managers/highlight-manager.js` | `cellucid/assets/js/app/ui/modules/highlight-controls.js` |
| “Where is multi-view (‘Keep view’) handled?” | `cellucid/assets/js/app/state/managers/view-manager.js` | `cellucid/assets/js/app/ui/modules/view-controls.js` |
| “Where are sessions saved/loaded?” | `cellucid/assets/js/app/session/session-serializer.js` | `app/session/history-file.js` (format 4), `app/session/README.md`, {doc}`10_sessions_persistence_and_serialization` |
| “Where is analysis computed?” | `cellucid/assets/js/app/analysis/` | {doc}`11_analysis_architecture` |
| “Where is figure export implemented?” | `cellucid/assets/js/app/ui/modules/figure-export/README.md` | {doc}`12_figure_export_architecture` |
| “Where is community annotation?” | `cellucid/assets/js/app/ui/modules/community-annotation-controls.js` | `cellucid/assets/js/app/community-annotations/REPO_SETUP.md` |

---

## Key terms (used consistently in code)

These show up everywhere (state, UI, analysis, sessions):

- **Dataset**: a logical unit with a stable `datasetId` (see {doc}`../b_data_loading/06_dataset_identity_why_it_matters`).
- **Field**: an `obs` categorical/continuous field, or a `var` “gene expression” field.
- **View**: the **live** view or a **kept** view (created by “Keep view”; the
  code calls it a snapshot view). An action names a view `"live"` or
  `"View N"`, its row in the Compare views list.
- **Action**: a named, parameterized change or read of the catalog
  (`color.by`, `filter.category.visibility`); every state change is one.
- **State domain**: one of the sixteen units state is captured, undone and
  saved in (`cells`, `mapping`, … `layout`).
- **Step**, **tracker**: one recorded action in the history; a tree of
  steps of one dataset.
- **Visibility**: what points are currently shown after filters/outliers.
- **Highlight**: annotation overlay on top of visibility (can exist even for filtered-out points).
- **Highlight page**: independent collection of highlight groups; switching pages changes overlay state.
- **Representation**: one embedding (an `obsm` key at 1 to 3 of its
  columns) a view can show, identified by its
  `layoutKey`. The data mapping turns any number on and names each
  dimension's default; a view's **dimension level** (1D/2D/3D) is that of the
  representation it shows.

---

## Troubleshooting (code navigation problems)

### Symptom: “I changed a file but nothing happens”

Likely causes (ordered):
1) Browser cache is serving old modules.
2) You edited a file not loaded in your current code path (e.g. a different module).
3) Your local server is serving the wrong directory (not `cellucid/`).
4) You edited a worker's code and opened the page under another name than
   `localhost`, `127.0.0.1` or `[::1]`: there it is a deployed site, and its
   workers start from their bundles in `assets/workers/` until
   `node scripts/generate-web-asset-inventory.mjs` rebuilds them
   (see {doc}`03_build_run_and_deployment`, checklist item 5).

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
