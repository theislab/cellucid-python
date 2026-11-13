# Reference (hooks, commands, schemas, endpoints)

This page is the “source of truth” reference for the notebook hooks system.

Primary sources in this repo:
- Python: `cellucid-python/src/cellucid/jupyter/` (the actions: `_actions.py`, `_actions_generated.py`, from `src/cellucid/action_catalog.json`)
- Python: `cellucid-python/src/cellucid/_viewer_channel.py` (the viewer channel)
- Python: `cellucid-python/src/cellucid/_server_base.py`
- Web app: `cellucid/assets/js/data/viewer-host-source.js`
- Web app: `cellucid/assets/js/app/viewer-host/host-channel.js` and `host-view.js` (the host's side of the channel and of the view's end), `cellucid/assets/js/app/main.js` (the start-up that makes them)

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
      var_gene_id_column=None,
      normalize_embeddings=None,
      centroid_min_points=None,
      dataset_name=None,
      dataset_id=None,
      obs_keys=None,
      embedding_defaults=None,
      knn_key=None,
      expression_source=None,
      vector_field_default=None,
      modality=None,
      web_source_url="https://www.cellucid.com",
      cache_bytes=None,
  )
  ```

  `dataset_name` and `dataset_id` default to `uns['cellucid']['name']` and
  `uns['cellucid']['id']`; data that declares neither is named after its file,
  or in memory as `AnnData (N cells)` (`MuData (N cells)`), with an id made
  from that name and printed (sessions are kept under it).
  Every launch choice (`latent_key` through `modality`) is optional:
  each one named is checked against the data's metadata before the viewer is
  shown and handed to it as a launch default; one left unset is taken
  from `uns['cellucid']['defaults']`, else the viewer's default rule.
  `obs_keys` names the exact `adata.obs` columns the viewer shows; every column
  is shown when it is omitted. `vector_field_default` names the vector field the
  overlay opens on. `web_source_url` is where the browser loads the app from.
  `cache_bytes` bounds what the kernel keeps of what it makes for a store
  served in priority order (a 64th of the machine's memory when None).
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

### Actions (Python → viewer)

Every action of the viewer's catalog is a method, generated from the
catalog this package ships (`cellucid.jupyter.action_catalog()`); each
returns an `ActionRequest` at once ({doc}`06_python_to_frontend_commands`):

- `viewer.act(type: str, /, **params) -> ActionRequest` -- any action, by its
  dotted type
- one method per action, its type with underscores:
  `viewer.color_by(field=...)`, `viewer.highlight_cells(cells=..., color=..., label=...)`,
  `viewer.highlight_clear_all()`, `viewer.cells_visibility_set(cells=..., visible=...)`,
  `viewer.camera_reset(view="live")`, `viewer.mapping_set(...)`,
  `viewer.gene_set_save(name=..., genes=...)`, `viewer.fields_list()`,
  `viewer.history_get()`, `viewer.history_undo()`, `viewer.history_redo()`,
  `viewer.history_jump(step=...)` (these three move the history and record no
  step), … (the full list: {doc}`../g_api_reference_coverage/api/viewers`)
- `request.result(timeout: float | None = 60.0)` → the outcome (a read's
  result); raises `ActionError` (`.reason`) when the viewer did not apply it;
  warns `ActionNotRecordedWarning` (`.unrecorded`) when a change made no
  history step; `request.done`, `request.step` (`{"id", "label"}` of the
  history step, or `None`), `request.unrecorded` (why a change made no step:
  `"unchanged"`, `"retracted"`, `"off"`, `"paused"`, `"elsewhere"`,
  `"no-history"`; `None` beside a step and for a read), `request.cancel()`
- `viewer.push(values, key=None) -> ActionRequest` (a pandas column or an
  array as a field: the `field.push` action)
- `viewer.send_message(message: dict)` (low-level entry point for the closed
  command schemas below)

### The history, live

- `with viewer.turn(label: str) as turn_id:` → every action sent inside the
  block is one turn of the viewer's history (`turn_id` is `python-` and 12
  hexadecimal digits); the State panel groups its steps and undoes them
  together. Turns do not nest (`RuntimeError`).
- `viewer.history` → a `LiveHistory`: `last_change` (the latest
  `history_changed` payload, or `None`), `tracker`, `head`, `branch`,
  `on_change(callback)`, `read()` (the same as `viewer.request_session()`).
  The history is changed through the generated `history_*` actions:
  `viewer.history_undo()`, `viewer.history_jump(step=...)`,
  `viewer.history_branch_switch(...)`, `viewer.history_note_add(...)`, …
- `loading = viewer.load_session(path)` → a `LoadSessionRequest`: the file's
  trackers are added beside the person's own; `loading.result(timeout: float
  | None = 60.0)` → `{"tracker": {"id", "name"}, "head"}` (raises
  `RuntimeError` with the viewer's reason when it refused), `loading.done`,
  `loading.cancel()`

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
- `@viewer.on_command_error` (a command this notebook sent was refused by the
  viewer; registering one also silences the default `sys.stderr` notice)
- `@viewer.on_data_mapping` (the viewer resolved or changed its data mapping;
  `viewer.data_mapping` holds the latest one)
- `@viewer.on_history_change` (the viewer's history changed: a step, a move,
  an edit, a tracker, a loaded session; `viewer.state.history` holds the
  latest one)
- `@viewer.on_message` (catches every current event type)

`session_loaded` has no decorator of its own: `viewer.load_session(...)`
collects it, and `viewer.register_hook("session_loaded", callback)` or
`@viewer.on_message` sees it.

Programmatic registration:
- `viewer.register_hook(event: str, callback) -> callback`
- `viewer.unregister_hook(event: str, callback) -> bool`
- `viewer.clear_hooks(event: str | None = None)`

`event` is one of the viewer's events (`cellucid.jupyter.HOOK_EVENTS`: every
event of {doc}`07_frontend_to_python_events` except an agent app's own) or
`message`; any other name raises `ValueError`, as a misspelt one would never
fire. R's `cellucid_on()` takes the same names.

Synchronous API (a notebook widget hears the viewer only between cells, so
the waiting calls raise at once in the cell that waits; see
{doc}`03_supported_environments_matrix`):
- `viewer.state` (latest event snapshot)
- `viewer.wait_for_event(event: str, timeout: float | None = 30.0, predicate=None) -> dict`
- `viewer.wait_for_ready(timeout: float | None = 30.0) -> dict`
- `viewer.data_mapping` (the latest reported data mapping, or `None`) and
  `viewer.data_mapping_origin`

Data mapping (the viewer reports it because this package lists `data_mapping` in `GET /_cellucid/protocol`):
- `viewer.mapping_set(**groups) -> ActionRequest` — the `mapping.set` action:
  change the named groups (`embeddings`, `knn`, `gene_names`, `expression`,
  `log1p`, `highly_variable`, `latent`, `vector_fields`, `modality`); the
  viewer reports the new mapping as a `data_mapping` event of origin
  `notebook`, which becomes `viewer.data_mapping`

### Sessions (the history, durable)

- `request = viewer.request_session()`, then
  `session = request.result(timeout: float | None = 60.0)` (a
  {class}`~cellucid.CellucidSession`) -- the one way, in a notebook widget
  (in a later cell) and a browser-tab viewer alike
- `session.save("path/to/file.cellucid-session")`, `session.trackers`,
  `session.tracker(name_or_id=None)` (a {class}`~cellucid.history.Tracker`:
  `steps` as a table, `step(id).state()`, `branch(name)`, `notes_table()`,
  `methods()`, `pseudocode()`), `session.conversations(tracker)`
- `plan = viewer.apply_session_to_anndata(adata, request, *, tracker=None,
  branch=None, step=None, dry_run=True, overwrite=False, keep_session=False,
  timeout=60.0, cleanup_file=True)` → an `ApplyPlan` (one row per change);
  `report = plan.apply()` writes it into `adata`
- `plan = session.apply(adata, *, tracker=None, branch=None, step=None,
  dry_run=True, overwrite=False, keep_session=False)`
- Function-level: `cellucid.apply_cellucid_session_to_anndata(adata, session, ...)`;
  `cellucid.history_applied(adata)` lists every apply recorded on `adata`
- Details: {doc}`11_session_bundles_apply_session_to_anndata`

### Connectivity / cache utilities

- `viewer.debug_connection(timeout=5.0)` → structured connectivity report
  with `channel_connected` and `recent_events["command_errors"]`, the most
  recent commands the viewer refused. A widget's report adds
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

### `action_result` (binary)

```python
{
  "requestId": str,            # the action's request
  "step": {"id": str, "label": str} | None,   # the history step it made; None for a read
  "unrecorded": str | None,    # why a change made no step: "unchanged", "retracted", "off",
                               # "paused", "elsewhere", "no-history"; None beside a step and for a read
  "outcome": object,           # the frame's body, UTF-8 JSON: a change's outcome or a read's result
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
  "command": str,          # one of the command types below
  "requestId": str | None, # the failed command's request; None for one that names none
  "reason": str,           # one line: the viewer sends at most 1000 characters, Python accepts 4000
}
```

### `data_mapping`

```python
{
  "origin": str,            # "default", "best_guess", "session", "user" or "notebook"
  "mapping": dict,          # the complete data mapping, validated and canonical
}
```

### `history_changed`

```python
{
  "kind": str,              # "step", "moved", "edited", "tracker" or "loaded"
  "tracker": str,           # the id of the tracker in use
  "head": str,              # the step shown
  "step": str | None,       # the step recorded (kind "step"), else None
  "branch": str | None,     # the id of the active branch
}
```

### `session_loaded`

```python
{
  "requestId": str,                       # the load_session command's request
  "tracker": {"id": str, "name": str},    # the tracker the viewer shows now
  "head": str,                            # its step shown
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
fields below; a command carries no `viewerId` or `viewerToken`. An `action`
whose parameters hold NumPy arrays (cells, a column's values) is one binary
frame instead: a little-endian `uint32` header length, the JSON header
(`type`, `requestId`, `action`, `params` without the arrays, `turn`,
`turnLabel`, and `arrays`: each array's `path` in `params`, `type` and
`length`), padded with spaces to a multiple of 4 bytes, then each array's
elements, little-endian, in order. `load_session` is always binary: its
header, padded the same way, then the session file. A command sent before the
viewer connects waits and is sent, in order, when it does.

### `action`

```python
{"type": "action", "requestId": "...", "action": "color.by", "params": {"field": {"source": "obs", "key": "cell_type"}},
 "turn": None, "turnLabel": None}
# binary header, sent inside `with viewer.turn("Marker genes"):`; then the cells, uint32 little-endian
{"type": "action", "requestId": "...", "action": "highlight.cells", "params": {"color": "#e41a1c"},
 "turn": "python-0f3a9c21b7e4", "turnLabel": "Marker genes",
 "arrays": [{"path": ["cells"], "type": "uint32", "length": 3}]}
```

`turn` and `turnLabel` are both text (the turn's id, and what it was asked to
do) or both null. Answered by one `action_result` with the same `requestId`,
or a `command_error` naming `action` and the `requestId`.

### `load_session` (binary)

```python
# binary header; then the .cellucid-session file
{"type": "load_session", "requestId": "...", "name": "my-session.cellucid-session"}
```

`name` is the file's name, without a folder. Answered by one
`session_loaded` with the same `requestId`, or a `command_error` naming
`load_session` and the `requestId`.

### Session/diagnostics

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
`runtime`: what opened the viewer, `"jupyter"` for a viewer started in a
notebook kernel, `"python"` in a script or session, `"cellucid-serve"` for
`cellucid serve` and `"cellucid-mcp"` for `cellucid mcp --open`, which every
sentence of the viewer names.

### `GET /_cellucid/datasets`

Returns the one dataset the server serves: its id, name, store route and launch choices.

### `GET /_cellucid/protocol`

Returns where the viewer channel is and what it carries: `channel` (its path,
`_cellucid/ws`), `events` (the inbound event types it validates),
`commands` (the commands it may send, which is also the set of names a
`command_error` may report) and `binaryCommands` (those it may send as binary
frames: `action` and `load_session`). The lists grow between releases; there is no
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
- commands: JSON text frames, and binary frames for an `action` whose
  parameters hold arrays and for `load_session`; events: JSON text frames (at most 1 MiB each) and
  binary frames for `selection`, `result`, `action_result` and
  `session_bundle` (uint32 little-endian header length, JSON header, body)
- close codes: 1001 the viewer or its server was stopped, 4409 a newer channel of the same
  viewer replaced this one, 1008 a frame broke the protocol
- the server pings every 20 s; a session file is streamed to a temporary
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
| `cache_bytes` | Every entry point: bytes a store served in priority order keeps of what it makes (renumbered genes), least recently used first out; None keeps a 64th of the machine's memory, 0 nothing (cellucid-r's `cache_bytes` alike) | A large object on a machine whose memory other work needs, or a large machine with room to spare |

---

## Known limitations (important)

- **Cells are positional** (row index of `obs`), not stable cell IDs; name
  cells by obs name where the order of your object may differ.
- `event["cells"]` is a Python `list[int]`, about 36 bytes per cell in memory
  (4 bytes a cell on the wire).
- Hook callbacks run on the viewer channel's thread, one event at a time; keep them fast and avoid heavy work inline.
- Applying a session checks only what the chosen step needs (the cell
  count, the dataset id when `uns['cellucid']['id']` names one, the data
  mapping and gene count, the cell order when `adata` has the session's
  embedding, the fields the step uses); any mismatch raises before AnnData is
  changed.

---

## Related pages

- Overview: {doc}`01_overview_bidirectional_communication`
- Architecture: {doc}`05_architecture_viewer_channel`
- Commands: {doc}`06_python_to_frontend_commands`
- Events: {doc}`07_frontend_to_python_events`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
