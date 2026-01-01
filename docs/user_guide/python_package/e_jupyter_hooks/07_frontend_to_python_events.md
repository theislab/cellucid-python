# Frontend → Python events

This page documents what the embedded Cellucid viewer can send back to Python, and how to consume those events safely.

Under the hood:
- the iframe sends JSON to `POST /_cellucid/events`
- the Python server routes the event by `viewerId`
- the viewer triggers callbacks registered via `@viewer.on_*`

## At a glance

**Audience**
- Wet lab / beginner: use `@viewer.on_selection` and print the size.
- Computational users: read payload schemas + edge cases.
- Developers: read delivery semantics + `/_cellucid/events` details in {doc}`05_architecture_message_routing_http_vs_postmessage`.

## Supported hooks (Python API)

All hooks are decorators attached to the `viewer` object:

```python
@viewer.on_selection
def handle(event):
    ...
```

Supported decorators:
- `@viewer.on_ready`
- `@viewer.on_selection`
- `@viewer.on_hover`
- `@viewer.on_click`
- `@viewer.on_message` (fires for *all* event types; best for debugging and custom events)

## Event catalog (what the viewer can emit)

### `ready`

Fires when the viewer has finished wiring notebook integration and the UI is ready.

Payload (typical):

```python
{
  "n_cells": 12345,
  "dimensions": 3
}
```

Register:

```python
@viewer.on_ready
def on_ready(event):
    print("Ready:", event)
```

### `selection`

Fires when the user completes a selection step (lasso, KNN, proximity, annotation selection, etc.).

Payload (typical):

```python
{
  "cells": [0, 10, 42, ...],   # 0-based row indices
  "source": "lasso"           # e.g. "lasso", "knn", "proximity", "annotation", "click", ...
}
```

Register:

```python
@viewer.on_selection
def on_selection(event):
    cells = event.get("cells", [])
    print("Selected:", len(cells), "source:", event.get("source"))
```

### `hover`

Fires when the user hovers a cell (debounced/throttled).

Payload (typical):

```python
{
  "cell": 123,                     # int or None when not hovering a cell
  "position": {"x": 0.1, "y": 0.2, "z": -0.3}  # may be null in some contexts
}
```

Register:

```python
@viewer.on_hover
def on_hover(event):
    cell = event.get("cell")
    if cell is None:
        return
    print("Hover:", cell)
```

### `click`

Fires when the user clicks on a cell.

Payload (typical):

```python
{
  "cell": 123,
  "button": 0,     # 0=left, 1=middle, 2=right
  "shift": False,
  "ctrl": False    # ctrl/cmd
}
```

Register:

```python
@viewer.on_click
def on_click(event):
    print("Clicked cell:", event.get("cell"))
```

### `console` (debugging)

Best-effort forwarding of frontend warnings/errors to Python (used by `viewer.debug_connection()`).

Payload (typical):

```python
{
  "level": "warn" | "error",
  "message": "...",
  "ts": "2026-01-01T12:34:56.789Z",
  # optional: filename/lineno/colno for window.onerror
}
```

Consume via `@viewer.on_message` if you want:

```python
@viewer.on_message
def debug_all(event):
    if event.get("event") == "console":
        print(event)
```

### `session_bundle` (no-download session capture)

Emitted when the viewer uploads a session bundle for `viewer.get_session_bundle()`.

Payload (typical):

```python
{
  "requestId": "...",
  "status": "ok" | "error",
  "bytes": 123456,     # on success
  "path": "/tmp/....cellucid-session",  # on success (server temp path)
  "error": "..."       # on error
}
```

See:
- {doc}`10_session_bundles_get_session_bundle`

### `pong` and `debug_snapshot` (diagnostics)

Emitted in response to `viewer.debug_connection()` probes.

These are primarily for debugging, but you can also hook them if you want:

```python
@viewer.on_message
def on_any(event):
    if event.get("event") in {"pong", "debug_snapshot"}:
        print(event)
```

### Custom event types

The frontend can send arbitrary events. Any event with at least:

```json
{
  "type": "<something>",
  "viewerId": "<id>"
}
```

will be routed to:
- `@viewer.on_message` (always), and
- `viewer.wait_for_event("<something>")` (if you use the synchronous API).

## How `@viewer.on_message` differs from other hooks

For any event type `X`, the message hook receives:

```python
{"event": "X", **payload_without_type_and_viewerId}
```

Whereas `@viewer.on_selection` receives only the selection payload.

This makes `@viewer.on_message` ideal for:
- debugging (“show me everything”)
- handling custom events without changing the Python package

## Delivery semantics (important in real notebooks)

### Transport
- Events are **HTTP POST** requests from the iframe to `/_cellucid/events`.
- The frontend sends events as “fire-and-forget” (errors are not surfaced to the UI by default).

### Size limits
The Python server currently enforces a **1MB request limit** for `/_cellucid/events`.

Implications:
- Very large selections (hundreds of thousands of indices) can exceed the limit and will be rejected.
- In that situation, prefer higher-level workflows (e.g. exporting highlight groups via session bundles) or reduce event payload size.

### Hover frequency
Hover is throttled/debounced to avoid flooding Python.
Treat hover callbacks as **best-effort**, and keep them fast.

### Threading
Hook callbacks can run on the server’s request-handling thread.

Rules of thumb:
- keep callbacks quick,
- catch exceptions,
- do heavy work outside the callback (see {doc}`08_writing_robust_callbacks`).

## Troubleshooting

- “No events arrive” → {doc}`14_troubleshooting_hooks`
- Run:
  ```python
  viewer.debug_connection()
  ```

## Next steps

- Writing robust callbacks: {doc}`08_writing_robust_callbacks`
- Viewer state + waiting: {doc}`09_viewer_state_and_wait_for_event`
- Full schemas and endpoints: {doc}`16_reference`
