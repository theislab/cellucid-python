# Python → frontend commands

This page documents how your notebook Python code can **control** the embedded Cellucid viewer.

Under the hood:
- the viewer opens one authenticated WebSocket to the data server (the viewer
  channel, {doc}`05_architecture_viewer_channel`);
- Python sends each command on it as one JSON text frame, and a command that
  names cells (`highlight`, `setVisibility`) as one binary frame;
- a command sent before the viewer has connected waits and is sent, in order,
  when it connects; nothing is sent twice.

If you’re new, start with the quickstart: {doc}`02_quickstart_minimal_roundtrip`.

## At a glance

**Audience**
- Wet lab / beginner: use the “Fast path recipes”.
- Computational users: read “Edge cases” + “Performance”.
- Developers: read “Message schema” + {doc}`05_architecture_viewer_channel`.

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
A command sent before the viewer has connected waits in Python and is sent as
soon as the viewer connects; the viewer applies it once its dataset is open.
`viewer.send_message(...)` returns `True` when the command went out at once and
`False` when it is waiting. After `viewer.stop()` it raises `RuntimeError`.
To act only once the viewer shows the data:

```python
viewer.wait_for_ready(timeout=60)
```
```

## Knowing whether a command applied

Commands are **one-way**. `viewer.set_color_by("CD8A")` returns `None` as soon
as the command is sent (or queued for the viewer); it does not wait for the
viewer, and its return value says nothing about whether the viewer applied it.

When the viewer rejects a command — an unknown field, an out-of-range index, a
malformed color — it reports the failure back over the event channel as a
`command_error` event.

Cellucid prints it:

```text
[cellucid] viewer.set_color_by(...) did not take effect in the viewer: Jupyter set-color-by command requires exactly one field "CD8A"
[cellucid] The viewer is unchanged. Handle this in code with @viewer.on_command_error.
```

To act on it in code, register a hook, which also silences that notice:

```python
rejected = []

@viewer.on_command_error
def record(event):
    rejected.append(event)
```

Record rather than raise: hooks run on the viewer channel's thread, so raising
there only logs the error and drops the event instead of failing your cell.

Two limits are worth knowing:

- The report is **asynchronous**. It arrives after the cell that sent the
  command has finished, so the failing call itself still returns `None`.
- It names the **command type**, not the individual call, because commands
  carry no request id.

See {doc}`07_frontend_to_python_events` for the exact record.

## Command catalog (public Python API)

All of these are methods on the `viewer` object returned by `show_anndata(...)`.

| Python method | Purpose | Frontend message type |
|---|---|---|
| `viewer.send_message(message)` | Low-level escape hatch; returns `True` (sent) or `False` (queued) | `message["type"]` |
| `viewer.highlight_cells(cell_indices, color="#ff0000")` | Highlight specific cells by index | `"highlight"` |
| `viewer.clear_highlights()` | Remove highlights | `"clearHighlights"` |
| `viewer.set_color_by(field)` | Change active color-by field | `"setColorBy"` |
| `viewer.set_visibility(cell_indices=None, visible=True)` | Hide/show specific cells (or all) | `"setVisibility"` |
| `viewer.reset_view()` | Reset camera | `"resetCamera"` |
| `viewer.set_data_mapping(**fields)` | Change what the viewer plots and reads, and wait for it to apply | `"setDataMapping"` |

Diagnostics / utilities (also useful for debugging):

| Python method | Purpose | Frontend message type |
|---|---|---|
| `viewer.debug_connection()` | End-to-end connectivity report | `"ping"`, `"debug_snapshot"` (internal) |
| `viewer.get_session_bundle()` | Request a session bundle upload | `"requestSessionBundle"` |
| `viewer.stop()` | Freeze the view, stop the server, unregister hooks, and report any failure | `"freeze"` |

## Message schema (what actually crosses the boundary)

### Commands carry no credentials

The viewer channel is authenticated once, when the viewer opens it, and
belongs to that one viewer, so a command is only its message-specific fields
below. `viewer.send_message({...})` refuses a message that names `viewerId`
or `viewerToken`.

### Commands that name cells are binary frames

`highlight` and `setVisibility` travel as one binary frame each: a
little-endian `uint32` header length, the UTF-8 JSON header (padded with
spaces to a multiple of 4 bytes, so the cells start on a 4-byte boundary and
the viewer reads them in place), then the cells as little-endian `uint32`,
strictly ascending. A frame that breaks these rules closes the channel
(1008), as a malformed event does in the other direction. You give the cells
as a list of ints or a one-dimensional NumPy integer array, in any order, each
once; Python sorts them once and writes them straight into the frame.

`viewer.send_message()` takes these commands in the same form as the others,
with their cells under `"cells"`:

```python
viewer.send_message({"type": "highlight", "cells": [0, 10, 42], "color": "#ff00ff"})
```

### `highlight`

Header, then the cells:

```json
{ "type": "highlight", "color": "#ff00ff" }
```

Notes:
- `cells` (in `send_message`, or `cell_indices`) are 0-based integer indices,
  at least one.
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

Header, then the cells (at least one):

```json
{ "type": "setVisibility", "visible": false, "all": false }
```

Notes:
- `viewer.set_visibility(None, visible=...)` targets every cell: the header
  says `"all": true` and the frame carries no cells.
- An empty list is refused; pass `None` for every cell.
- If you need reproducible filtering, prefer the UI’s filter system; cell-level visibility is best treated as an interactive convenience.

### `resetCamera`

```json
{ "type": "resetCamera" }
```

### `requestSessionBundle` (internal, used by `viewer.get_session_bundle()`)

```json
{ "type": "requestSessionBundle", "requestId": "..." }
```

### `setDataMapping` (used by `viewer.set_data_mapping()`)

```json
{ "type": "setDataMapping", "requestId": "...", "mapping": { "kind": "anndata", "...": "..." } }
```

`mapping` is a complete data mapping, validated by
`cellucid.data_mapping.validate_data_mapping` before it is sent; the
canonical copy is what crosses the wire. Unlike the one-way commands above,
this one is answered: the viewer replies with a `data_mapping` event of origin
`notebook` carrying the same `requestId`, or with a `command_error` naming
`setDataMapping` when it refuses the mapping.

`viewer.set_data_mapping(**fields)` replaces the named top-level fields of the
mapping the viewer last reported (`embeddings`, `vectorFields`, `geneNames`,
`knn`, `expression`, `log1p`, `highlyVariable`, `latent`; never `kind`, which
is always `"anndata"`),
validates the result, sends it, and returns the mapping the viewer reports
back. It raises `RuntimeError` before sending when the viewer has not reported
a mapping yet, `TypeError` for an unknown field, `ValueError` for an invalid
value, `RuntimeError` with the viewer's reason when the viewer refuses it, and
`TimeoutError` when no answer arrives within `timeout` seconds (30 by default).

```python
viewer.set_data_mapping(
    embeddings={
        "enabled": [{"key": "X_pca", "columns": [0, 1]}],
        "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [0, 1]}, "3d": None},
    },
)
```

### `ping` / `debug_snapshot` (internal, used by `viewer.debug_connection()`)

```json
{ "type": "ping", "requestId": "..." }
```

```json
{ "type": "debug_snapshot", "requestId": "..." }
```

## Performance notes

- **Huge cell lists** travel as binary frames, 4 bytes a cell: pass a NumPy
  array (`np.flatnonzero(mask)`) rather than a list for millions of cells, so
  the cells are checked and written in vectorized passes.
- **High-frequency commands** (sending a highlight every mousemove) are not recommended; prefer debouncing and batching.
- If you want to drive many interactions, consider using hooks events to *decide* what to do, but send commands only when needed.

## Edge cases and footguns

- **Indexing is positional**: indices refer to row positions, not stable cell IDs.
- **Out-of-range indices** violate the current command contract; validate
  indices against the dataset cell count before sending them.
- **Viewer not displayed yet**: commands wait and are sent when it connects.
- **Viewer stopped**: `viewer.send_message(...)` raises `RuntimeError`.
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
3. Check whether the viewer rejected the command rather than ignoring it:
   ```python
   viewer.debug_connection()["recent_events"]["command_errors"]
   ```
4. If ping/pong succeeds and no command was rejected, but a specific command
   is ignored:
   - Re-establish the exact configured source generation:
     ```python
     viewer.clear_web_cache()
     ```
   - Re-run the cell to create a fresh viewer.
5. Use the full guide: {doc}`14_troubleshooting_hooks`

## Next steps

- Viewer → Python events: {doc}`07_frontend_to_python_events`
- Robust callback patterns: {doc}`08_writing_robust_callbacks`
- Full reference: {doc}`16_reference`
