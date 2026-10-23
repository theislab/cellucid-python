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

- **A Cellucid-ready AnnData writer**: `cellucid.prepare(...)` writes one standard AnnData (`.h5ad`, `.zarr`, `.zarr.zip`) laid out for the browser's range reads.
- **One server**: `cellucid.anndata_server.AnnDataServer` serves the bytes of one *AnnData* store (`.h5ad` or Zarr), which the web app reads itself with range requests.
- **Notebook embedding + hooks**: `cellucid.jupyter.AnnDataViewer` / `show_anndata(...)` embed the web app in notebooks and enable Python ↔ frontend interaction.
- **Session bundle tooling**: read `.cellucid-session` bundles and apply them back onto AnnData (`CellucidSessionBundle`, `apply_cellucid_session_to_anndata`).
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
| `src/cellucid/prepare_data/` | `prepare(...)`: input checks (`_inputs.py`), the element writer with its chunks and codecs (`_store.py`), the public function (`_prepare.py`) | want to change what a Cellucid-ready AnnData holds or how it is laid out |
| `src/cellucid/_launch_defaults.py` | the launch-default check shared by `prepare()` (`uns['cellucid']['defaults']`) and the server (`/_cellucid/datasets`) | want to add or change a launch default |
| `src/cellucid/_contracts.py` | dataset-id and display-text rules, standard library only | want to change what counts as a valid id or name |
| `src/cellucid/_non_finite.py` | the counted NaN/Inf/float32-range diagnosis in refusals | want to change how a non-finite input is reported |
| `src/cellucid/_anndata_files.py` | classifies an `.h5ad` / Zarr folder / zipped Zarr store, reads its metadata (never an array), finds arrays the browser cannot decode, builds missing consolidated Zarr metadata, writes the temporary `.h5ad` for an in-memory object | want to change what AnnData input is accepted or how it is handed to the browser |
| `src/cellucid/anndata_server.py` | AnnData server: checks launch choices against the store's metadata, serves the store's bytes with exact `Range` answers, publishes `/_cellucid/datasets` | want to add/modify server endpoints |
| `src/cellucid/jupyter/` | notebook embedding, postMessage commands, hooks/event routing, session capture | want to change notebook UX or hooks behavior |
| `src/cellucid/_server_base.py` | the HTTP server base: CORS, `Host` checks, byte ranges, web-asset serving, event + session upload endpoints | want to change security/transport behavior |
| `src/cellucid/session_bundle.py` | `.cellucid-session` streaming reader | want to parse sessions or add chunk-level tooling |
| `src/cellucid/session_codecs.py` | session chunk codecs (varint, cell sets, user-defined codes, typed JSON) | want to add/align codecs with the web app |
| `src/cellucid/anndata_session/` | apply session bundle → AnnData | want to change “sessions to AnnData bridge” semantics |
| `src/cellucid/vector_fields.py` | CellRank drift → vector fields | want to add helpers for velocity/drift conventions |
| `src/cellucid/web_cache.py` | exact inventory-driven viewer generation establishment and verification | want to change source, staging, verification, or atomic publication behavior |

---

## The core contract: an AnnData store

Every dataset reaches the viewer as an AnnData store, read as it is: the web
app's HDF5/Zarr reader reads the `.h5ad` file or Zarr store directly (from a
server route, a link, or a file picked in the browser). A file written by
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

1) Python: `AnnDataServer` classifies the store (`.h5ad`, Zarr folder, zipped
   Zarr; an in-memory object is written once to a temporary `.h5ad`), reads its
   metadata, and checks the launch choices against it.
2) Browser: reads `/_cellucid/datasets`, then opens the store at
   `/_cellucid/anndata/<dataset_id>.h5ad` (or `.zarr.zip`, or `.zarr/<key>`).
3) Browser: to show `FOXP1`, its HDF5/Zarr reader issues HTTP range requests
   for exactly that gene's bytes; Python answers each with `206` and the bytes
   on disk.

### C) Notebook embedding + hooks workflow

1) Python: starts an AnnData server and creates a `viewerId`.
2) Notebook: renders an iframe pointing at the server (possibly via proxy URL).
3) Python → frontend: commands via `postMessage` (includes `viewerToken`).
4) Frontend → Python: events via HTTP POST to `/_cellucid/events` (routed by `viewerId`).

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
- or that iframes preserve a stable origin.

Notebook constraints drive the existence of:
- same-origin local viewer generation delivery,
- the explicit `client_server_url` argument,
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
