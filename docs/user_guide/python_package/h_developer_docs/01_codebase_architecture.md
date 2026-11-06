# Codebase architecture

This page explains the **architecture of the `cellucid` Python package** (repo folder: `cellucid-python/`): what the core subsystems are, how data flows through them, and which invariants you must preserve when changing code.

:::{admonition} Audience
:class: note

- If you are a wet-lab user trying to “make the viewer appear”: read the *Fast path* and *Troubleshooting* sections only.
- If you are a computational user: focus on the Cellucid-ready AnnData + server behavior.
- If you are a developer: read the whole page before editing `prepare_data/`, `anndata_server.py`, or `jupyter/`.
:::

---

## What the Python package is responsible for (and what it is not)

### Python package responsibilities

`cellucid-python` provides:

- **A Cellucid-ready AnnData writer**: `cellucid.prepare(...)` writes one standard AnnData (`.h5ad`, `.zarr`, `.zarr.zip`), or from a MuData a Cellucid-ready MuData (`.h5mu`, MuData Zarr), laid out for the browser's range reads.
- **One server**: `cellucid.anndata_server.AnnDataServer` serves the bytes of one dataset in any format the web app reads, classified by the input format registry (`input_formats.json`: AnnData, MuData, Loom, 10x `.h5`/`outs/`, MatrixMarket folders, R objects), which the web app reads itself — with range requests where the format allows it, once from the start for an R object or a MatrixMarket folder.
- **Notebook embedding + hooks**: `cellucid.jupyter.AnnDataViewer` / `show_anndata(...)` embed the web app in notebooks and enable Python ↔ frontend interaction.
- **Session tooling**: read `.cellucid-session` files (the history of every step: trackers, branches, steps, notes, turns, the journal, the state at each step; `read_session`, `cellucid.history`), and plan, then apply, the state at any step onto AnnData or MuData, each apply recorded (`apply_cellucid_session_to_anndata`, `ApplyPlan`, `ApplyReport`, `history_applied`).
- **Convenience utilities**: vector-field derivation helpers, hosted web UI caching helpers.

### Not Python package responsibilities

- Rendering, UI state machines, analysis panels, figure export internals live in the **web app repo** (`cellucid/`).
  Start with: {doc}`../../web_app/p_developer_docs/index`.
- Community annotation workflows are primarily web-app-driven (repo template lives in `cellucid-annotation/`).
- R is owned by `cellucid-r`; its `cellucid_prepare()` writes the same
  Cellucid-ready AnnData the web app reads.

---

## Source tree: the “one screen” map

Most of the public API lives in `cellucid-python/src/cellucid/`:

| File | What it does | Read when you… |
|---|---|---|
| `src/cellucid/__init__.py` | Lazy public API exports via `__getattr__` | want to change public imports without slowing CLI startup |
| `src/cellucid/prepare_data/` | `prepare(...)`: what `data=` (an AnnData or MuData) supplies (`_data.py`), input checks (`_inputs.py`), latent-space outlier quantiles (`_outliers.py`), the element writer with its chunks and codecs (`_store.py`), the Cellucid-ready MuData layout (`_mudata.py`), the public function (`_prepare.py`) | want to change what a Cellucid-ready AnnData holds or how it is laid out |
| `src/cellucid/_launch_defaults.py` | the launch-default check shared by `prepare()` (`uns['cellucid']['defaults']`) and the server (`/_cellucid/datasets`) | want to add or change a launch default |
| `src/cellucid/data_mapping.py` | the data-mapping rules (validation, unique names, the default UMAP rule, vector-field keys), mirrored by the web app's `data/data-mapping.js` and pinned by a shared conformance fixture | want to change what the viewer plots and reads |
| `src/cellucid/_contracts.py` | dataset-id and display-text rules, standard library only | want to change what counts as a valid id or name |
| `src/cellucid/_non_finite.py` | the counted NaN/Inf/float32-range diagnosis in refusals | want to change how a non-finite input is reported |
| `src/cellucid/_anndata_files.py` | classifies any file or folder by the input format registry (`input_formats.json`), reads AnnData and MuData metadata (never an array), finds arrays the browser cannot decode, builds missing consolidated Zarr metadata | want to change what input is accepted or how it is handed to the browser |
| `src/cellucid/prepare_data/_gene_major.py` | the one gene-major layout pass (`GeneMajorLayout`): any sparse source, in memory or on disk, in any row order, as `prepare()` stores it (canonical, float32 or float64), with the prefix counts; shared by `prepare()`, the backed path and the served store | want to change how an expression matrix is laid out |
| `src/cellucid/_source_arrays.py`, `src/cellucid/_scratch.py` | arrays of a source read by rows or ranges in any order (in memory, mapped in place, or decoded once into a scratch folder); the scratch folders, removed with their owner and by the next one when a process was killed | want to change how a source array is read |
| `src/cellucid/anndata_server.py` | the data server: checks launch choices against AnnData/MuData metadata, serves the bytes with exact `Range` answers (`sendfile`), lists a folder of files, publishes `/_cellucid/datasets`, hands `/_cellucid/ws` to the viewer channel | want to add/modify server endpoints |
| `src/cellucid/jupyter/` | the viewer API (`_viewer.py`: commands, hooks/event dispatch, results, session capture, notebook turns (`turn()`), `load_session()`; a browser-tab viewer), the live history (`_history.py`: `LiveHistory`, `viewer.history`, from `history_changed`; `LoadSessionRequest`), its actions (`_actions.py`: `act()` checks and converts an action's parameters by the catalogue and sends one `action` command; `_actions_generated.py`: one method per action, written by `scripts/generate_actions.py` from `src/cellucid/action_catalog.json`, the web app's catalogue byte for byte, SHA-256 pinned in `tests/test_action_catalog.py`), the notebook widget (`_widget.py` with `static/widget.js`: an anywidget relaying the app's reads and channel over the comm); `_wire.py` checks every message | want to change notebook UX or hooks behavior |
| `src/cellucid/_viewer_channel.py` | the viewer channel: WebSocket upgrade and `Origin` check, `hello`/`welcome`, the per-viewer registry and command outbox, events, session bundles streamed to disk | want to change how the viewer and Python talk |
| `src/cellucid/_websocket.py` | RFC 6455 framing over an upgraded HTTP connection (masked, fragmented, in-place unmasking) | want to change the WebSocket layer |
| `src/cellucid/_server_base.py` | the HTTP server base: CORS, `Host` checks, byte ranges, web-asset serving, `405` for `POST` | want to change security/transport behavior |
| `src/cellucid/session.py` | the `.cellucid-session` reader (session format 4, the web app's `app/session/history-file.js`): framing, manifest, the chunks in their order (eager before lazy; steps, parts, array, values, journal, chat), every chunk inflated to its declared size (or stored as it is), arrays checked by digest, parts joined once (`_PartJoiner`, `partJoiner` in `value-parts.js`: a part that names itself and a value past 2^29 − 24 characters refused before anything is expanded; the values the head shows from the eager chunks alone) and values by hash once joined as `app/history/value-pool.js` hashes them, every record against the history model; a file written by a newer Cellucid refused as newer (`newer_session_message`); `read_session()` / `CellucidSession` (the header, `trackers`, `tracker()`, `conversations()`, `save()`, `apply()`) | want to read more of a session, or follow a change of the format |
| `src/cellucid/history.py` | the reading model: `Tracker` (the step table, `step()`, `branch()`, `group()`, `turn()`, `path_to()`, `notes_table()`, `stale()`, `methods()`, `pseudocode()`), `Branch`, `Group`, `Step` (of each kind: root, action, outside, squashed; `state()`, `diff()`) and `Diff`, with the web app's words for diffs, the methods paragraph and pseudocode (`app/history/exports/`) | want to change how a history reads in Python |
| `src/cellucid/_history_model.py`, `src/cellucid/_history_rules.py` | the history model as frozen dataclasses, generated by `scripts/generate_history_model.py` from `history_model.schema.json` (the web app's `historyModelDocument` byte for byte, SHA-256 pinned in `tests/test_history_model.py`, the web and R suites); the rules the types cannot say (a tracker's whole tree, step kinds, squashed steps, groups, names, turns, the journal, what a scope says it holds), in the web app's words | want to follow a change of the history model (regenerate; never edit `_history_model.py` by hand) |
| `src/cellucid/session_codecs.py` | the codecs of a session's array chunks (uvarint, cell sets), mirroring `app/session/codecs/` | want to align a codec with the web app |
| `src/cellucid/anndata_session/` | the round trip: `pull`, results → `obs` columns, GMT, obs-column CSV; and applying one step of a session (`_apply.py`: the dataset checks, the plan, one row per change, the atomic write and its record in `uns['cellucid']['history_applied']`; the same rules as cellucid-r, both held to the golden history's expected plans) | want to change how results or a session's state reach AnnData |
| `src/cellucid/vector_fields.py` | CellRank drift → vector fields | want to add helpers for velocity/drift conventions |
| `src/cellucid/web_cache.py` | exact inventory-driven viewer generation establishment and verification | want to change source, staging, verification, or atomic publication behavior |
| `src/cellucid/cli.py` | the `cellucid` command (`cellucid serve`, `cellucid mcp`) | want to add or change a command-line option |
| `src/cellucid/mcp.py` | `cellucid mcp`, the local MCP server for agent apps (the official SDK: initialize, tools, calls, cancellation, the person asked in the app) | want to change what an agent app is answered |
| `src/cellucid/_mcp_core.py` | one agent connection without the transport: which viewer it drives, the tools (`mcp_tools.json`, generated by the web repository), the app's name, `cellucid_status`, the sentences an agent reads | want to change how a viewer is chosen or what is said when none is |
| `src/cellucid/_agent_link.py` | a viewer's end for agents (`ViewerAgentPort`), its publication on this computer (the agents folder and the loopback endpoint), and a published viewer driven from another process (`RemoteViewerLink`) | want to change how `cellucid mcp` reaches a viewer |
| `src/cellucid/subsample.py` | `subsample_cells()`: the subsampling specification, implemented as the web app's `data/subsample.js` and R's `subsample.R` (pinned by the shared conformance fixture) | want to change which cells a subsample keeps (change all three) |
| `src/cellucid/_file_versions.py` | strong validators of a file's version, read from the open file: the server's `ETag` and dataset list, and the web cache's re-check (R's `.file_version()` builds the same tag) | want to change when a file counts as changed |
| `src/cellucid/_console.py`, `src/cellucid/_plural.py` | encoding-safe console output; the noun a count takes | want to change user-facing console text |

---

## The core contract: the dataset as stored

Every dataset reaches the viewer as the file or folder it is on disk, read as
it is: the web app's own readers (HDF5, Zarr, R serialization, MatrixMarket
text) read it directly, from a server route, a link, or a file picked in the
browser. An in-memory AnnData, and an AnnData file not in priority order, are
served by Python as a Zarr store in priority order made as the viewer reads it
(`src/cellucid/_served_store.py`), and so is an in-memory MuData, as a MuData
store. A file written by
`prepare()` follows the Cellucid-ready AnnData contract
({doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`):
gene-major, chunked for range reads, and naming itself in `uns['cellucid']`.
A server only adds `/_cellucid/datasets`, which names the dataset and states
the launch choices — see
{doc}`09_server_mode_architecture_endpoints_and_security`.

---

## Data flow diagrams (high level)

### A) Preparing a dataset (fast + shareable)

1) Python: `prepare(...)` checks the inputs and writes one Cellucid-ready
   AnnData (CPU + disk work, once).
2) Anything that opens AnnData opens it: a server (below), the web app's file
   buttons, or a link to wherever it is hosted.

### B) AnnData server workflow (nothing converted)

1) Python: `AnnDataServer` classifies the path by the input format registry
   (`classify_anndata_path`), presents an in-memory AnnData, an AnnData file
   not in priority order, and a MuData, Loom or 10x file as a Zarr store in
   priority order made as the browser reads it (`_served_store.py`,
   `_served_formats.py`; an in-memory MuData as a MuData store), reads its
   metadata (an R object and a MatrixMarket folder are
   read by the browser alone), and checks the launch choices against it.
2) Browser: reads `/_cellucid/datasets`, then opens the store at
   `/_cellucid/anndata/<dataset_id><ending>` (the file's own ending, such as
   `.h5ad`, `.zarr.zip`, `.loom` or `.rds.gz`), `<dataset_id>.zarr/<key>` for a
   Zarr folder, or `<dataset_id>/<file>` for a folder of files, whose file list
   `/_cellucid/datasets` carries.
3) Browser: to show `FOXP1`, its HDF5/Zarr reader issues HTTP range requests
   for exactly that gene's bytes; Python answers each with `206` and the bytes
   on disk.

### C) Notebook embedding + hooks workflow

1) Python: `show_anndata()` prepares the dataset (an `AnnDataServer` that is
   never bound), creates a `viewerId` and `viewerToken`, and shows a `Widget`.
2) Notebook: the widget's output frames `https://www.cellucid.com/?viewer=embedded`
   and hands the page the id, the token and a `MessagePort`.
3) Frontend: the viewer reads the data through the port (the widget relays
   each request to the kernel), opens its channel and says hello with its
   `viewerId` and `viewerToken`.

A browser-tab viewer (`AnnDataViewer`) instead binds the server (port 8765, or
the next free port when another viewer holds it and no port was given), and
the page at `<server>/?viewer=true&viewerId=…&viewerToken=…&anndata=true`
opens its channel at `/_cellucid/ws`.
4) Python → frontend: commands on the channel; frontend → Python: events on
   the same channel.

---

## Architecture diagram (recommended)

If you want a single figure you can refer to in bug reports and onboarding, add an architecture diagram.

---

## Invariants (do not break these without a coordinated frontend change)

### 1) Cell indices are positional

Frontend selections/highlights are expressed as **row indices** (0-based positions).

Implications:
- If you reorder cells, the same indices refer to different biological cells.
- For “transfer state” between datasets, you need stable dataset identity + explicit mapping logic.

### 2) The Cellucid-ready AnnData is a shared contract

Users will:
- prepare a file once,
- host/share it,
- and expect it to keep opening.

If you change element names, layouts, `uns['cellucid']` keys, or launch-default
semantics, you must change cellucid-r's writer, the web app's reader, and the
specification page together.

### 3) Notebook environments are hostile to assumptions

Do not assume:
- the browser can reach `http://127.0.0.1:<port>` (remote kernels, Colab),
- mixed-content is allowed (HTTPS notebooks),
- network access is available (airgapped environments),
- or that the kernel hears the viewer while a cell runs (it does not: a
  widget's messages arrive between cells).

Notebook constraints drive the existence of:
- the notebook widget, which needs no port, address or proxy,
- for browser-tab viewers, same-origin local viewer generation delivery and
  the explicit `client_server_url` argument,
- and the debugging helpers in `viewer.debug_connection()`.

---

## Troubleshooting (quick routing)

### Symptom: “I changed Python and now the viewer fails to load data”

Likely causes:
- you broke an invariant of the Cellucid-ready AnnData (element names, encodings, codecs),
- Python and the configured web source publish different current contracts,
- you introduced a CORS/origin mismatch.

Start with:
- {doc}`07_prepare_writer_architecture` (writer invariants),
- {doc}`09_server_mode_architecture_endpoints_and_security` (HTTP + CORS checks),
- {doc}`12_debugging_playbook` (how to capture minimal repro).

### Symptom: “Hooks/events stopped working”

Start with:
- {doc}`11_hooks_events_protocol_and_schema`,
- {doc}`10_jupyter_embedding_architecture`,
- then run `viewer.debug_connection()` and follow {doc}`12_debugging_playbook`.
