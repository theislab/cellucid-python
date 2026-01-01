# Codebase architecture

This page explains the **architecture of the `cellucid` Python package** (repo folder: `cellucid-python/`): what the core subsystems are, how data flows through them, and which invariants you must preserve when changing code.

:::{admonition} Audience
:class: note

- If you are a wet-lab user trying to “make the viewer appear”: read the *Fast path* and *Troubleshooting* sections only.
- If you are a computational user: focus on the export format + server behavior.
- If you are a developer: read the whole page before editing `prepare_data.py`, `server.py`, `anndata_server.py`, or `jupyter.py`.
:::

---

## What the Python package is responsible for (and what it is not)

### Python package responsibilities

`cellucid-python` provides:

- **An export pipeline**: `cellucid.prepare(...)` writes a static directory that the web app can load fast and reproducibly.
- **Two server implementations**:
  - `cellucid.server.CellucidServer`: serves *pre-exported* directories (static files).
  - `cellucid.anndata_server.AnnDataServer`: serves *AnnData* as virtual export files (dynamic).
- **Notebook embedding + hooks**: `cellucid.jupyter.CellucidViewer` / `show(...)` / `show_anndata(...)` embed the web app in notebooks and enable Python ↔ frontend interaction.
- **Session bundle tooling**: read `.cellucid-session` bundles and apply them back onto AnnData (`CellucidSessionBundle`, `apply_cellucid_session_to_anndata`).
- **Convenience utilities**: vector-field derivation helpers, hosted web UI caching helpers.

### Not Python package responsibilities

- Rendering, UI state machines, analysis panels, figure export internals live in the **web app repo** (`cellucid/`).
  Start with: {doc}`../../web_app/p_developer_docs/index`.
- Community annotation workflows are primarily web-app-driven (repo template lives in `cellucid-annotation/`).
- R export (`cellucid-r`) is planned but not yet the authoritative path.

---

## Source tree: the “one screen” map

Most of the public API lives in `cellucid-python/src/cellucid/`:

| File | What it does | Read when you… |
|---|---|---|
| `src/cellucid/__init__.py` | Lazy public API exports via `__getattr__` | want to change public imports without slowing CLI startup |
| `src/cellucid/prepare_data.py` | `prepare(...)` export pipeline + manifest formats | want to change export layout, quantization, metadata |
| `src/cellucid/server.py` | static-file server for exported datasets | want to change serving behavior for export folders |
| `src/cellucid/anndata_adapter.py` | AnnData → virtual export format adapter | want to change lazy loading, caching, mapping of AnnData fields |
| `src/cellucid/anndata_server.py` | dynamic server that exposes adapter data via HTTP routes | want to add/modify endpoints for AnnData mode |
| `src/cellucid/jupyter.py` | notebook embedding, postMessage commands, hooks/event routing, session capture | want to change notebook UX or hooks behavior |
| `src/cellucid/_server_base.py` | shared CORS, web-asset proxy, event + session upload endpoints | want to change security/transport behavior |
| `src/cellucid/session_bundle.py` | `.cellucid-session` streaming reader | want to parse sessions or add chunk-level tooling |
| `src/cellucid/session_codecs.py` | session chunk codecs (varint, delta-varint, RLE codes) | want to add/align codecs with the web app |
| `src/cellucid/anndata_session.py` | apply session bundle → AnnData | want to change “sessions to AnnData bridge” semantics |
| `src/cellucid/vector_fields.py` | CellRank drift → vector fields | want to add helpers for velocity/drift conventions |
| `src/cellucid/web_cache.py` | hosted viewer UI asset caching (offline-safe notebooks) | want to change caching, invalidation, prefetch behavior |

---

## The core contract: “the viewer loads the export format”

Whether the data comes from:
- static files on disk (export folder), or
- an AnnData adapter generating bytes on demand,

the browser viewer uses the *same* conceptual file contract:

- points: `points_2d.bin(.gz)`, `points_3d.bin(.gz)`, …
- obs manifest + field binaries: `obs_manifest.json`, `obs/<field>...`
- var manifest + gene binaries: `var_manifest.json`, `var/<gene>...`
- (optional) connectivity + vectors: `connectivity_manifest.json`, `connectivity/...`, `vectors/...`
- metadata: `dataset_identity.json`

This contract is documented in detail here:
{doc}`08_export_format_spec_and_invariants`.

---

## Data flow diagrams (high level)

### A) Exported dataset workflow (fast + shareable)

1) Python: `prepare(...)` writes export folder (CPU + disk work).
2) Python: `CellucidServer` serves that folder over HTTP.
3) Browser: loads the UI (proxied/cached) + fetches dataset files over HTTP.

### B) AnnData server workflow (convenient + lazy)

1) Python: `AnnDataServer` opens AnnData (in-memory, backed `.h5ad`, or `.zarr`).
2) Browser: requests “files” like `/var/FOXP1.values.f32`.
3) Python: `AnnDataAdapter` reads that gene column lazily and returns bytes (optionally gzip).

### C) Notebook embedding + hooks workflow

1) Python: starts a server (exported or AnnData) and creates a `viewerId`.
2) Notebook: renders an iframe pointing at the server (possibly via proxy URL).
3) Python → frontend: commands via `postMessage` (includes `viewerToken`).
4) Frontend → Python: events via HTTP POST to `/_cellucid/events` (routed by `viewerId`).

---

## Architecture diagram (recommended)

If you want a single figure you can refer to in bug reports and onboarding, add an architecture diagram.

<!-- SCREENSHOT PLACEHOLDER
ID: python-architecture-overview
Suggested filename: developer/python-architecture-overview.svg
Where it appears: Python Package → Developer Docs → Codebase architecture
Capture:
  - This is a diagram (not a screenshot). Export as SVG if possible.
Content:
  - Boxes: prepare/export, exported server, AnnData adapter/server, Jupyter viewer, web app
  - Arrows: HTTP file fetches, postMessage commands, HTTP POST events
  - Note: show “hosted-asset proxy cache” as a separate box (because it explains offline behavior)
Alt text:
  - Block diagram showing how cellucid-python exports/serves data and embeds the web app in notebooks.
Caption:
  - High-level architecture: the browser viewer always speaks the export-format protocol; Python either serves files from disk or synthesizes them from AnnData.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder for a high-level cellucid-python architecture diagram.
:width: 100%

High-level architecture: the viewer loads the export-format contract; Python serves it from disk or synthesizes it from AnnData.
```

---

## Invariants (do not break these without a coordinated frontend change)

### 1) Cell indices are positional

Frontend selections/highlights are expressed as **row indices** (0-based positions).

Implications:
- If you reorder cells, the same indices refer to different biological cells.
- For “transfer state” between datasets, you need stable dataset identity + explicit mapping logic.

### 2) Export format compatibility is real user data

Users will:
- export once,
- host/share the folder,
- and expect it to keep working.

If you change:
- filenames,
- manifest schema,
- quantization encoding,
- or reserved sentinel values,

you must coordinate with the web app and document the compatibility story.

### 3) Notebook environments are hostile to assumptions

Do not assume:
- the browser can reach `http://127.0.0.1:<port>` (remote kernels, Colab),
- mixed-content is allowed (HTTPS notebooks),
- network access is available (airgapped environments),
- or that iframes preserve a stable origin.

Notebook constraints drive the existence of:
- the hosted-asset proxy,
- `CELLUCID_CLIENT_SERVER_URL`,
- and the debugging helpers in `viewer.debug_connection()`.

---

## Troubleshooting (quick routing)

### Symptom: “I changed Python and now the viewer fails to load data”

Likely causes:
- you broke an invariant in the export format (paths, gzip flags, types),
- the viewer is still using cached assets (stale UI build),
- you introduced a CORS/origin mismatch.

Start with:
- {doc}`08_export_format_spec_and_invariants` (format checks),
- {doc}`09_server_mode_architecture_endpoints_and_security` (HTTP + CORS checks),
- {doc}`12_debugging_playbook` (how to capture minimal repro).

### Symptom: “Hooks/events stopped working”

Start with:
- {doc}`11_hooks_events_protocol_and_schema`,
- {doc}`10_jupyter_embedding_architecture`,
- then run `viewer.debug_connection()` and follow {doc}`12_debugging_playbook`.
