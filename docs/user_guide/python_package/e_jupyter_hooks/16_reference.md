# Reference (hooks, commands, schemas, endpoints)

This page is the “source of truth” reference for the notebook hooks system.

Primary sources in this repo:
- Python: `cellucid-python/src/cellucid/jupyter/`
- Python: `cellucid-python/src/cellucid/_viewer_channel.py` (the viewer channel)
- Python: `cellucid-python/src/cellucid/_server_base.py`
- Web app: `cellucid/assets/js/data/viewer-host-source.js`
- Web app: `cellucid/assets/js/app/main.js`

If you are new, start with:
- {doc}`02_quickstart_minimal_roundtrip`

---

## Public Python API (what you can call)

### Creating a viewer

```python
from cellucid import show_anndata
```

- `show_anndata(...) -> cellucid.Widget` (a notebook widget; same API as
  `AnnDataViewer`), with this closed call surface:

  ```text
  show_anndata(
      data,
      height=600,
      *,
      latent_key=None,
      gene_id_column=None,
      normalize_embeddings=None,
      centroid_min_points=None,
      dataset_name=None,
      dataset_id=None,
      obs_keys=None,
      embedding_keys=None,
      knn_key=None,
      expression_source=None,
      vector_field_default=None,
      modality=None,
      web_source_url="https://www.cellucid.com",
  )
  ```

  `dataset_name` and `dataset_id` default to `uns['cellucid']['name']` and
  `uns['cellucid']['id']`, and are required when the data declares neither.
  Every launch choice (`latent_key` through `modality`) is optional:
  each one named is checked against the data's metadata before the viewer is
  shown and handed to it as a launch default; one left unset is taken
  from `uns['cellucid']['defaults']`, else the viewer's default rule.
  `obs_keys` names the exact `adata.obs` columns the viewer shows; every column
  is shown when it is omitted. `vector_field_default` names the vector field the
  overlay opens on. `web_source_url` is where the browser loads the app from.
  The widget opens no port; for a viewer in a browser tab of its own (a
  script, an SSH tunnel) construct `AnnDataViewer`, which also takes `port`,
  `client_server_url`, `web_cache_dir` and `allowed_hosts`.

### Basic properties

- `viewer.state`, `viewer.data_mapping`, `viewer.results`, `viewer.selection`, `viewer.obs_names`
- `AnnDataViewer` only (a widget has no address and raises `AttributeError`):
  - `viewer.server_url` → the underlying data server URL (usually `http://127.0.0.1:<port>`)
  - `viewer.viewer_url` → the browser-tab URL (includes `viewer=true`, `viewerId`, `viewerToken`, and `anndata=true`)

### Display and lifecycle

- `viewer.display()` → show the widget in the cell's output (an
  `AnnDataViewer` shows a link to its browser tab)
- `viewer.stop()` → freeze the displayed view, stop answering it (an
  `AnnDataViewer`: stop its server), unregister hooks, and report any freeze
  or shutdown failure
- `cellucid.jupyter.cleanup_all()` → stop all active viewers (also registered via `atexit`)

### Commands (Python → viewer)

- `viewer.send_message(message: dict)` (low-level entry point for the closed
  command schemas below)
- `viewer.highlight_cells(cell_indices: list[int], color: str = "#ff0000")`
- `viewer.clear_highlights()`
- `viewer.set_color_by(field: str)`
- `viewer.set_visibility(cell_indices: list[int] | None = None, visible: bool = True)`
- `viewer.reset_view()`
- `viewer.push(values, key=None)` (a column as a field) and
  `viewer.push_gene_set(name, genes, description="")`

### Results back (one call)

- `summary = viewer.pull(adata)` → every highlight group, user-defined field,
  score and gene set the viewer holds, written onto `adata.obs` /
  `adata.uns['cellucid']['gene_sets']`
- `viewer.subset(adata, which=None)`, `viewer.gene_sets()`
- Function-level, also for the web app's exports: `cellucid.pull(adata, source)`
- Details: {doc}`../g_api_reference_coverage/api/round_trip`

### Hooks (viewer → Python)

Decorators:
- `@viewer.on_ready`
- `@viewer.on_selection`
- `@viewer.on_hover`
- `@viewer.on_click`
- `@viewer.on_command_error` (a command this notebook sent was rejected by the
  viewer; registering one also silences the default `sys.stderr` notice)
- `@viewer.on_data_mapping` (the viewer resolved or changed its data mapping;
  `viewer.data_mapping` holds the latest one)
- `@viewer.on_message` (catches all nine current event types)

Programmatic registration:
- `viewer.register_hook(event: str, callback) -> callback`
- `viewer.unregister_hook(event: str, callback) -> bool`
- `viewer.clear_hooks(event: str | None = None)`

Synchronous API (a notebook widget hears the viewer only between cells, so
the waiting calls raise at once in the cell that waits; see
{doc}`03_supported_environments_matrix`):
- `viewer.state` (latest event snapshot)
- `viewer.wait_for_event(event: str, timeout: float | None = 30.0, predicate=None) -> dict`
- `viewer.wait_for_ready(timeout: float | None = 30.0) -> dict`
- `viewer.data_mapping` (the latest reported data mapping, or `None`) and
  `viewer.data_mapping_origin`

Data mapping (the viewer reports it because this package lists `data_mapping` in `GET /_cellucid/protocol`):
- `viewer.set_data_mapping(**fields) -> None` — replace top-level mapping
  fields, validate, and send `setDataMapping` without waiting (at once, or as
  soon as the viewer reports a mapping); the viewer's answer becomes
  `viewer.data_mapping`, a refusal a `command_error`

### Session bundles (durable state)

- `request = viewer.request_session_bundle()`, then
  `bundle = request.result(timeout: float | None = 60.0)` -- the one way, in
  a notebook widget (in a later cell) and a browser-tab viewer alike
- `bundle.save("path/to/file.cellucid-session")`
- `adata2 = viewer.apply_session_to_anndata(adata, request, inplace=False)`
- `adata2 = bundle.apply_to_anndata(adata,
  expected_dataset_id="my-study-v1", inplace=False)`
- Function-level: `cellucid.apply_cellucid_session_to_anndata(...)`

### Connectivity / cache utilities

- `viewer.debug_connection(timeout=5.0)` → structured connectivity report
  with `channel_connected` and `recent_events["command_errors"]`, the most
  recent commands the viewer rejected. A widget's report adds
  `connected_outputs` and `web_source_url`; an `AnnDataViewer`'s adds its
  server's probes, including `anndata_store_probes`, keyed by every exact
  listed dataset id (a `HEAD` and a one-byte `Range` read of the served
  AnnData store)
- `AnnDataViewer` only:
  - `viewer.ensure_web_ui_cached(force=True, show_progress=True)` → establish
    the complete source UI generation
  - `viewer.ensure_web_ui_cached(force=False, show_progress=True)` → verify
    the selected existing generation without network access or mutation
  - `viewer.clear_web_cache()` → clear cached viewer UI assets

---

## Event schemas (viewer → Python)

Browser events are delivered on the viewer channel (`/_cellucid/ws`), which
the viewer authenticated with its `viewerId`/`viewerToken` pair when it opened
it. Most are JSON text frames; `selection` is built by the server from a binary
frame of cell indices, and `session_bundle` from a binary frame holding the
bundle (see {doc}`05_architecture_viewer_channel`).

Python receives payloads **without** `type`; no event carries `viewerId` or
`viewerToken`. Every listed field is required and undeclared fields are
rejected.

### `ready`

```python
{"n_cells": int, "n_obs": int, "dimensions": int}  # cells open, cells in the dataset
```

### `selection`

```python
{"cells": list[int], "source": str}
```

`cells` holds unique row indices in ascending order.

### `hover`

```python
{"cell": int | None, "position": dict | None}
```

### `click`

```python
{"cell": int, "button": int, "shift": bool, "ctrl": bool}
```

### `pong`

```python
{"requestId": str, "t": int}
```

### `debug_snapshot`

```python
{
  "requestId": str,
  "ts": str,
  "locationHref": str,
  "origin": str,
  "serverUrl": str,
  "connected": bool,
  "userAgent": str | None,
}
```

### `session_bundle`

```python
{
  "requestId": str,
  "status": "ok",
  "bytes": int,
  "path": str,
}
```

### `command_error`

```python
{
  "command": str,   # one of the command types below
  "reason": str,    # one line: the viewer sends at most 240 characters, Python accepts 500
}
```

### `data_mapping`

```python
{
  "requestId": str | None,  # a str exactly when origin is "notebook"
  "origin": str,            # "default", "best_guess", "session", "user" or "notebook"
  "mapping": dict,          # the complete data mapping, validated and canonical
}
```

### `@viewer.on_message` envelope

For each current event type `X`, `@viewer.on_message` receives:

```python
{"event": "X", **payload}
```

---

## Command schemas (Python → viewer)

Commands are sent on the viewer channel as JSON text frames, exactly the
fields below; a command carries no `viewerId` or `viewerToken`. A command that
names cells (`highlight`, `setVisibility`) is one binary frame instead: a
little-endian `uint32` header length, the JSON header (exactly the header
fields below, padded with spaces to a multiple of 4 bytes), then the cells as
little-endian `uint32`, strictly ascending. `viewer.send_message()` takes them
with their cells under `"cells"` (a list of ints or a NumPy integer array, in
any order). A command sent before the viewer connects waits and is sent, in
order, when it does.

### `highlight` (binary)

```python
{"type": "highlight", "color": "#ff00ff", "label": None}  # header; then the cells, at least one
```

### `clearHighlights`

```python
{"type": "clearHighlights"}
```

### `setColorBy`

```python
{"type": "setColorBy", "field": "cell_type"}
```

### `setVisibility` (binary)

```python
{"type": "setVisibility", "visible": False, "all": False}  # header; then the cells, at least one
{"type": "setVisibility", "visible": True, "all": True}    # every cell; no cells follow
```

### `resetCamera`

```python
{"type": "resetCamera"}
```

### `setDataMapping` (sent by `viewer.set_data_mapping()`)

```python
{"type": "setDataMapping", "requestId": "...", "mapping": {...}}  # a complete data mapping
```

Answered by a `data_mapping` event of origin `notebook` with the same
`requestId`, or a `command_error` naming `setDataMapping`.

### Session/diagnostics (internal)

```python
{"type": "requestSessionBundle", "requestId": "..."}
{"type": "ping", "requestId": "..."}
{"type": "debug_snapshot", "requestId": "..."}
{"type": "freeze"}
```

---

## Server endpoints (debugging + integration)

All of these are served by the Python data server that backs a browser-tab
viewer (`AnnDataViewer`, `serve_anndata`, `cellucid serve`). A notebook widget
answers the same requests over the notebook's connection instead, with no
server:

### `GET /_cellucid/health`

Returns `status` (`"ok"`), `type` (`"anndata"`), `version`, `format` (the
input format id, such as `"h5ad"`), `n_cells` and `n_genes` (`null` for a
format only the browser reads). Used for connectivity probing (and by the
viewer once, to tell a stopped Python from a proxy that refuses WebSockets);
it is not used to discover or select a proxy.

### `GET /_cellucid/info`

Returns `version`, `type`, `format`, `host`, `port`, `n_cells`, `n_genes` and
`runtime`: the program serving the viewer, `"jupyter"` for a viewer started in
a notebook kernel and `"python"` otherwise, which the viewer names.

### `GET /_cellucid/datasets`

Returns the one dataset the server serves: its id, name, store route and launch choices.

### `GET /_cellucid/protocol`

Returns where the viewer channel is and what it carries: `channel` (its path,
`_cellucid/ws`), `events` (the inbound event types it validates) and
`commands` (the commands it may send, which is also the set of names a
`command_error` may report). The lists grow between releases; there is no
version number. A web build asks this route before it opens the channel and
sends only the events listed. See
{doc}`../h_developer_docs/09_server_mode_architecture_endpoints_and_security`.

### `GET /_cellucid/ws` (WebSocket)

The viewer channel ({doc}`05_architecture_viewer_channel`).

Behavior:
- the upgrade passes the `Host` check (421) and an `Origin` check (403): the
  server's own host, a loopback name, or a host in `allowed_hosts`
- the first frame must be `{"type": "hello", "viewerId", "viewerToken"}`; the
  answer is `{"type": "welcome"}`, or close code 4401 for an unknown viewer or
  a wrong token
- commands: JSON text frames, and binary frames for `highlight`, `setField`
  and `setVisibility`; events: JSON text frames (at most 1 MiB each) and binary
  frames for `selection`, `result` and `session_bundle` (uint32 little-endian header
  length, JSON header, body)
- close codes: 1001 the viewer or its server was stopped, 4409 a newer channel of the same
  viewer replaced this one, 1008 a frame broke the protocol
- the server pings every 20 s; a session bundle is streamed to a temporary
  file, with no size cap but the temporary directory's free space

The server answers no `POST` (405).

---

## Connectivity and generation arguments

| Argument | Meaning | When to use |
|---|---|---|
| `web_source_url` | Where the app comes from: loaded by the browser for a widget; the exact origin publishing `cellucid-web-assets.json` for a browser-tab server | A copy of the app you host; testing |
| `web_cache_dir` | `AnnDataViewer`: directory for the active verified web generation | Controlled writable location or explicit inspection |
| `client_server_url` | `AnnDataViewer`: exact browser-facing HTTP(S) server base | A reverse proxy or tunnel in front of the server |
| `allowed_hosts` | `AnnDataViewer`, `serve_anndata`: extra exact host names the data server answers to (`Host`) and accepts viewer pages from (`Origin`); one bare host name per entry, no port, scheme or wildcard | Required behind a reverse proxy, which forwards the browser's `Host` verbatim |

---

## Known limitations (important)

- **Indices are positional** (row index), not stable cell IDs.
- `event["cells"]` is a Python `list[int]`, about 36 bytes per cell in memory
  (4 bytes a cell on the wire).
- Hook callbacks run on the viewer channel's thread, one event at a time; keep them fast and avoid heavy work inline.
- Session bundle application requires the exact `expected_dataset_id` and
  matching cell/variable counts; any fingerprint mismatch raises before
  AnnData is mutated.

---

## Related pages

- Overview: {doc}`01_overview_bidirectional_communication`
- Architecture: {doc}`05_architecture_viewer_channel`
- Commands: {doc}`06_python_to_frontend_commands`
- Events: {doc}`07_frontend_to_python_events`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
