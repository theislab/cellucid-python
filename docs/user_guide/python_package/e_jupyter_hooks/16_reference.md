# Reference (hooks, commands, schemas, endpoints)

This page is the “source of truth” reference for the notebook hooks system.

Primary sources in this repo:
- Python: `cellucid-python/src/cellucid/jupyter.py`
- Python: `cellucid-python/src/cellucid/_server_base.py`
- Web app: `cellucid/assets/js/data/jupyter-source.js`
- Web app: `cellucid/assets/js/app/main.js`
- Dev notes: `cellucid/markdown/HOOKS_DEVELOPMENT.md`

If you are new, start with:
- {doc}`02_quickstart_minimal_roundtrip`

---

## Public Python API (what you can call)

### Creating a viewer

```python
from cellucid import show, show_anndata
```

- `show("./export_dir", height=600) -> CellucidViewer`
- `show_anndata(adata_or_path, height=600, **adapter_kwargs) -> AnnDataViewer`

### Basic properties

- `viewer.server_url` → the underlying data server URL (usually `http://127.0.0.1:<port>`)
- `viewer.viewer_url` → the iframe URL (includes `jupyter=true`, `viewerId`, `viewerToken`, and `anndata=true` for AnnDataViewer)

### Display and lifecycle

- `viewer.display()` → (re)display iframe in a notebook output cell
- `viewer.stop()` → freeze view (best effort) then stop server + unregister hooks
- `cellucid.jupyter.cleanup_all()` → stop all active viewers (also registered via `atexit`)

### Commands (Python → viewer)

- `viewer.send_message(message: dict)` (low-level)
- `viewer.highlight_cells(cell_indices: list[int], color: str = "#ff0000")`
- `viewer.clear_highlights()`
- `viewer.set_color_by(field: str)`
- `viewer.set_visibility(cell_indices: list[int] | None = None, visible: bool = True)`
- `viewer.reset_view()`

### Hooks (viewer → Python)

Decorators:
- `@viewer.on_ready`
- `@viewer.on_selection`
- `@viewer.on_hover`
- `@viewer.on_click`
- `@viewer.on_message` (catches all events)

Programmatic registration:
- `viewer.register_hook(event: str, callback) -> callback`
- `viewer.unregister_hook(event: str, callback) -> bool`
- `viewer.clear_hooks(event: str | None = None)`

Synchronous “pull” API:
- `viewer.state` (latest event snapshot)
- `viewer.wait_for_event(event: str, timeout: float | None = 30.0, predicate=None) -> dict`
- `viewer.wait_for_ready(timeout: float | None = 30.0) -> dict`

### Session bundles (durable state)

- `bundle = viewer.get_session_bundle(timeout: float | None = 60.0)`
- `bundle.save("path/to/file.cellucid-session")`
- `adata2 = viewer.apply_session_to_anndata(adata, inplace=False, **kwargs)`
- `adata2 = bundle.apply_to_anndata(adata, inplace=False, **kwargs)`
- Function-level: `cellucid.apply_cellucid_session_to_anndata(...)`

### Connectivity / cache utilities

- `viewer.debug_connection(...)` → structured connectivity report
- `viewer.ensure_web_ui_cached(force=False, show_progress=True)` → prefetch viewer UI assets
- `viewer.clear_web_cache()` → clear cached viewer UI assets

---

## Event schemas (viewer → Python)

All events are delivered via `POST /_cellucid/events` and routed by `viewerId`.

Python receives payloads **without** the `type` and `viewerId` keys.

### `ready`

```python
{"n_cells": int, "dimensions": int}
```

### `selection`

```python
{"cells": list[int], "source": str}
```

### `hover`

```python
{"cell": int | None, "position": dict | None}
```

### `click`

```python
{"cell": int, "button": int, "shift": bool, "ctrl": bool}
```

### `console` (best-effort debug forwarding)

```python
{"level": "warn" | "error", "message": str, "ts": str, ...optional_fields...}
```

### `session_bundle`

```python
{
  "requestId": str,
  "status": "ok" | "error",
  "bytes": int,          # on success
  "path": str,           # on success (temp file path on server machine)
  "error": str           # on error
}
```

### Diagnostic events

Used by `viewer.debug_connection()`:
- `pong`
- `debug_snapshot`

### `@viewer.on_message` envelope

For any event type `X`, `@viewer.on_message` receives:

```python
{"event": "X", **payload}
```

---

## Command schemas (Python → viewer)

All commands are sent via `postMessage` into the iframe.

You typically send only the message-specific fields; the embedding layer injects:
- `viewerId`
- `viewerToken`

### `highlight`

```python
{"type": "highlight", "cells": [0, 10, 42], "color": "#ff00ff"}
```

### `clearHighlights`

```python
{"type": "clearHighlights"}
```

### `setColorBy`

```python
{"type": "setColorBy", "field": "cell_type"}
```

### `setVisibility`

```python
{"type": "setVisibility", "cells": [0, 10, 42], "visible": False}
```

### `resetCamera`

```python
{"type": "resetCamera"}
```

### Session/diagnostics (internal)

```python
{"type": "requestSessionBundle", "requestId": "..."}
{"type": "ping", "requestId": "..."}
{"type": "debug_snapshot", "requestId": "..."}
{"type": "freeze"}
```

---

## Server endpoints (debugging + integration)

All of these are served by the Python data server that backs the viewer:

### `GET /_cellucid/health`

Used for:
- connectivity probing
- notebook proxy selection

### `GET /_cellucid/info`

Returns server metadata (version, mode, etc.).

### `GET /_cellucid/datasets`

Returns dataset listings (one dataset in AnnData mode; one or many in exported mode).

### `POST /_cellucid/events`

Hooks endpoint.

Behavior:
- expects JSON
- requires `viewerId`
- request size limit: **1MB** (to guard against accidental giant payloads)
- responds with JSON:
  - `{"status": "ok", "delivered": true}` if routed
  - `{"status": "not_found", "delivered": false}` if viewerId not registered

### `POST /_cellucid/session_bundle?viewerId=...&requestId=...`

Session bundle upload endpoint (Jupyter no-download capture).

Behavior:
- requires a pre-registered pending request
- streams upload to a temp file
- hard size cap: **512MB**
- validates a MAGIC header (`CELLUCID_SESSION\n`)

---

## Environment variables

| Variable | Meaning | When to use |
|---|---|---|
| `CELLUCID_WEB_PROXY_CACHE_DIR` | Where the hosted viewer UI assets are cached | Offline use, persistent cache, shared environments |
| `CELLUCID_CLIENT_SERVER_URL` | Override the browser-facing server URL used in iframe embeds | Remote kernels, custom proxies, unusual notebook frontends |

---

## Known limitations (important)

- **Indices are positional** (row index), not stable cell IDs.
- `/_cellucid/events` has a **1MB request limit**; huge selections can be rejected.
- Hook callbacks can run on a server thread; keep them fast and avoid heavy work inline.
- Session bundle application is only safe when applying to a dataset with matching row order (use mismatch policies to guard).

---

## Related pages

- Overview: {doc}`01_overview_bidirectional_communication`
- Architecture: {doc}`05_architecture_message_routing_http_vs_postmessage`
- Commands: {doc}`06_python_to_frontend_commands`
- Events: {doc}`07_frontend_to_python_events`
- Troubleshooting: {doc}`14_troubleshooting_hooks`

