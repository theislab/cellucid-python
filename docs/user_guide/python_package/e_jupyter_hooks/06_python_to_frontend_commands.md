# Python → frontend commands

This page documents how your notebook Python code can **control** the embedded Cellucid viewer.

Under the hood:
- Python injects a small JS snippet into the notebook output.
- That JS uses `postMessage(...)` to send a command into the iframe.
- Commands include a per-viewer secret (`viewerToken`) for authentication.

If you’re new, start with the quickstart: {doc}`02_quickstart_minimal_roundtrip`.

## At a glance

**Audience**
- Wet lab / beginner: use the “Fast path recipes”.
- Computational users: read “Edge cases” + “Performance”.
- Developers: read “Message schema” + {doc}`05_architecture_message_routing_http_vs_postmessage`.

## Fast path recipes (copy/paste)

### Highlight a few cells

```python
viewer.highlight_cells([0, 10, 42], color="#ff00ff")
```

### Clear all highlights

```python
viewer.clear_highlights()
```

### Reset camera

```python
viewer.reset_view()
```

### Color by an obs field

```python
viewer.set_color_by("cell_type")
```

### Hide a set of cells (advanced)

```python
viewer.set_visibility([0, 10, 42], visible=False)
```

```{important}
Most commands require the viewer to be **displayed** and (for best results) **ready**:

```python
viewer.wait_for_ready(timeout=60)
```

If you call commands before the iframe exists, the Python side will warn and do nothing.
```

## Command catalog (public Python API)

All of these are methods on the `viewer` object returned by `show(...)` or `show_anndata(...)`.

| Python method | Purpose | Frontend message type |
|---|---|---|
| `viewer.send_message(message)` | Low-level escape hatch | `message["type"]` |
| `viewer.highlight_cells(cell_indices, color="#ff0000")` | Highlight specific cells by index | `"highlight"` |
| `viewer.clear_highlights()` | Remove highlights | `"clearHighlights"` |
| `viewer.set_color_by(field)` | Change active color-by field | `"setColorBy"` |
| `viewer.set_visibility(cells=None, visible=True)` | Hide/show specific cells (or all) | `"setVisibility"` |
| `viewer.reset_view()` | Reset camera | `"resetCamera"` |

Diagnostics / utilities (also useful for debugging):

| Python method | Purpose | Frontend message type |
|---|---|---|
| `viewer.debug_connection()` | End-to-end connectivity report | `"ping"`, `"debug_snapshot"` (internal) |
| `viewer.get_session_bundle()` | Request a session bundle upload | `"requestSessionBundle"` |
| `viewer.stop()` | Freeze the view, stop the server, unregister hooks, and report any failure | `"freeze"` |

## Message schema (what actually crosses the boundary)

### Important: you usually do **not** include `viewerId`/`viewerToken` yourself

When you call `viewer.send_message({...})`, the embedding layer automatically injects:
- `viewerId`: routes the message to the correct iframe
- `viewerToken`: authenticates the command inside the iframe

You normally send only the message-specific fields below.

### `highlight`

```json
{
  "type": "highlight",
  "cells": [0, 10, 42],
  "color": "#ff00ff"
}
```

Notes:
- `cells` should be a list of 0-based integer indices.
- `color` must be a hex string (`#RRGGBB`); invalid colors are rejected.

### `clearHighlights`

```json
{ "type": "clearHighlights" }
```

### `setColorBy`

```json
{ "type": "setColorBy", "field": "cell_type" }
```

Notes:
- `field` is the exact current field ID exposed by the loaded dataset.
- Use an existing field ID; the command does not translate aliases.

### `setVisibility`

```json
{ "type": "setVisibility", "cells": [0, 10, 42], "visible": false }
```

Notes:
- `"cells": null` targets all cells and is what
  `viewer.set_visibility(None, visible=...)` emits.
- If you need reproducible filtering, prefer the UI’s filter system; cell-level visibility is best treated as an interactive convenience.

### `resetCamera`

```json
{ "type": "resetCamera" }
```

### `requestSessionBundle` (internal, used by `viewer.get_session_bundle()`)

```json
{ "type": "requestSessionBundle", "requestId": "..." }
```

### `ping` / `debug_snapshot` (internal, used by `viewer.debug_connection()`)

```json
{ "type": "ping", "requestId": "..." }
```

```json
{ "type": "debug_snapshot", "requestId": "..." }
```

## Performance notes

- **Huge cell lists** (hundreds of thousands) can be slow because they must cross a notebook boundary.
  - If you need to highlight “almost everything”, consider highlighting the *small* complement instead.
- **High-frequency commands** (sending a highlight every mousemove) are not recommended; prefer debouncing and batching.
- If you want to drive many interactions, consider using hooks events to *decide* what to do, but send commands only when needed.

## Edge cases and footguns

- **Indexing is positional**: indices refer to row positions, not stable cell IDs.
- **Out-of-range indices** violate the current command contract; validate
  indices against the dataset cell count before sending them.
- **Viewer not displayed**: `viewer.send_message(...)` raises `RuntimeError`.
- **Multiple viewers**: make sure you are sending commands to the same `viewer` instance that rendered the iframe you are looking at.

## Troubleshooting

If “nothing happens” when sending commands:

1. Confirm the viewer is displayed and ready:
   ```python
   viewer.wait_for_ready(timeout=60)
   ```
2. Run:
   ```python
   viewer.debug_connection()
   ```
   - If ping/pong fails, you have a connectivity/proxy problem.
3. If ping/pong succeeds but a specific command is ignored:
   - Re-establish the exact configured source generation:
     ```python
     viewer.clear_web_cache()
     ```
   - Re-run the cell to create a fresh viewer.
4. Use the full guide: {doc}`14_troubleshooting_hooks`

## Next steps

- Viewer → Python events: {doc}`07_frontend_to_python_events`
- Robust callback patterns: {doc}`08_writing_robust_callbacks`
- Full reference: {doc}`16_reference`
