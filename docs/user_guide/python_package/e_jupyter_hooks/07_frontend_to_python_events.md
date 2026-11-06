# Frontend → Python events

The embedded viewer sends closed event records to Python, all on the viewer
channel ({doc}`05_architecture_viewer_channel`): most as one JSON text frame
each; a selection, an action's answer, a result and a session as binary
frames. This page covers the interaction, diagnostic, action, session and
history events; the results a person makes (highlights, fields, gene sets) are
{doc}`../g_api_reference_coverage/api/round_trip`, and every record is in
{doc}`16_reference`. The
channel was authenticated with `viewerId` and `viewerToken` when the viewer
opened it and belongs to that viewer, so no event carries either; Python
dispatches each record to the viewer that owns the channel.

## Hook entry points

Use the event-specific decorators for interaction events:

```python
@viewer.on_selection
def handle_selection(event):
    print(event["cells"], event["source"])

@viewer.on_hover
def handle_hover(event):
    print(event["cell"], event["position"])

@viewer.on_click
def handle_click(event):
    print(event["cell"], event["button"])

@viewer.on_ready
def handle_ready(event):
    print(event["n_cells"], event["dimensions"])

@viewer.on_command_error
def handle_command_error(event):
    print(event["command"], event["requestId"], event["reason"])

@viewer.on_history_change
def handle_history_change(event):
    print(event["kind"], event["head"])
```

`@viewer.on_message` receives every accepted event in a common envelope:

```python
@viewer.on_message
def handle_current_event(event):
    print(event["event"], event)
```

The envelope contains `event` plus the exact payload fields below.

## Closed event catalog

Unknown event types and every record with missing or undeclared fields are
rejected before any hook runs.
Integers are native, non-negative JSON integers; booleans are native JSON
booleans.

### `selection`

A selection can hold millions of cells, so it does not travel as JSON. The
viewer sends it as one binary frame:

```text
<uint32 little-endian header length><{"type":"selection","source":"lasso"}><cells>
```

The cells are the selected row indices, four bytes each (little-endian
unsigned 32-bit), in strictly ascending order, at least one cell; `[0, 10,
42]` is the 12 bytes `00 00 00 00 0a 00 00 00 2a 00 00 00`. The server reads
them in place, without a second copy, and refuses a body that is empty, not a
whole number of 4-byte values, or not strictly ascending by closing the
channel (1008) before any hook runs. A selection of any size is taken; the
kernel's memory is the only bound. `source` is one exact token without
whitespace. A `selection` sent as a text frame is refused the same way.

Hook payload:

```python
{"cells": [0, 10, 42], "source": "lasso"}
```

`cells` is a `list[int]` of unique row indices in ascending order.

### `hover`

Text frame:

```json
{
  "type": "hover",
  "cell": 42,
  "position": {"x": 0.1, "y": 0.2, "z": -0.3}
}
```

`cell` and `position` are both `null` when no cell is hovered. A non-null
position contains exactly finite `x`, `y`, and `z` numbers.

### `click`

Text frame:

```json
{
  "type": "click",
  "cell": 42,
  "button": 0,
  "shift": false,
  "ctrl": true
}
```

`button` is exactly `0`, `1`, or `2`.

### `ready`

Text frame:

```json
{
  "type": "ready",
  "n_cells": 12345,
  "dimensions": 2
}
```

`dimensions` is exactly `1`, `2`, or `3`.

### `pong`

`pong` is the response to the internal `ping` diagnostic:

```python
{"requestId": "exact-request-id", "t": 1722000000000}
```

### `debug_snapshot`

`debug_snapshot` reports the page's environment for a browser-tab viewer's
`viewer.debug_connection()` (a notebook widget's report asks for none):

```python
{
    "requestId": "exact-request-id",
    "ts": "2026-07-26T12:00:00.000Z",
    "locationHref": "https://viewer.example/?viewer=true",
    "origin": "https://viewer.example",
    "serverUrl": "https://notebook.example/cellucid",
    "connected": True,
    "userAgent": "Browser identification",
}
```

`userAgent` may be `None`; every other field is required.

### `action_result`

The viewer's answer to one action ({doc}`06_python_to_frontend_commands`): one
binary frame, header `{"type": "action_result", "requestId": …, "step": …}`,
then the action's outcome as UTF-8 JSON -- a read's result can be as large as
it asks for (a picture of a view, the values of many cells), so it never
travels as a text frame. Python decodes the body and builds this record:

```python
{
    "requestId": "exact-request-id",
    "step": {"id": "s4", "label": "Coloured by cell_type"},   # or None: a read, or nothing changed
    "outcome": None,                                          # a change's outcome, or a read's result
}
```

It completes the `ActionRequest` of that `requestId`: `request.result()`
returns `outcome` and `request.step` is `step`. The body must be UTF-8 JSON
with finite numbers; anything else closes the channel (1008).

### `session_bundle`

The viewer sends the session `viewer.request_session()` asked for as one
binary frame (header `{"type": "session_bundle", "requestId": …}`, then the
`.cellucid-session` bytes of format 4: the whole history). Python streams it
to a temporary file, reads and checks it whole
({class}`~cellucid.CellucidSession`), consumes the one pending request it
answers, and builds this record:

```python
{
    "requestId": "exact-request-id",
    "status": "ok",
    "bytes": 123456,
    "path": "/temporary/path/session.cellucid-session",
}
```

An invalid file is removed and reported as a `command_error` naming
`requestSessionBundle` and the request, so the request's `result()` raises
with the reason instead of waiting for its timeout. A session nobody is
waiting for any more (its request timed out or was cancelled) is read and
dropped; the channel stays open.

### `command_error`

A command this notebook sent -- an action (`viewer.color_by(...)`,
`viewer.act(...)`), a session request -- was refused by the viewer instead of
carried out.

```{note}
The viewer sends this event only to a `cellucid` that lists `command_error`
among the events of `GET /_cellucid/protocol` (every current one does), so the
package never receives an event it would reject.
```

Text frame:

```json
{
  "type": "command_error",
  "command": "action",
  "requestId": "exact-request-id",
  "reason": "There is no field \"cell_typ\"; fields like it: cell_type."
}
```

Hook payload:

```python
{
    "command": "action",
    "requestId": "exact-request-id",
    "reason": "There is no field \"cell_typ\"; fields like it: cell_type.",
}
```

`command` must be one of the command types in
{doc}`06_python_to_frontend_commands` (`action`, `load_session`,
`requestSessionBundle`, `freeze`, `ping`, `debug_snapshot`); a record naming anything else is
rejected, so no browser-supplied text reaches your notebook as a command name.
`requestId` is the failed command's request -- an action's
`request.request_id` -- or `None` for a command that names none (`freeze`).
`reason` is one line: the viewer shortens it to 1000 characters, and Python
accepts at most 4000.

The event completes the request it names: the action's `request.result()`
raises `ActionError` with `reason`.

**If no `command_error` hook is registered, Cellucid prints the failure to
`sys.stderr`, naming the method that sent it:**

```text
[cellucid] viewer.color_by(...) did not take effect in the viewer: There is no field "cell_typ"; fields like it: cell_type.
[cellucid] The viewer is unchanged. Handle this in code with @viewer.on_command_error.
```

Registering a hook takes the report over and silences that notice, so the
failure is reported once, where you chose:

```python
rejected = []

@viewer.on_command_error
def record(event):
    rejected.append(event)
```

Do not raise from this hook. Like every hook it runs on the viewer channel's
thread, so raising only logs the error and drops the event instead of failing
your cell (see {doc}`08_writing_robust_callbacks`). To fail on a refusal,
collect the request instead, at a point you choose:

```python
request = viewer.color_by(field={"source": "obs", "key": "cell_type"})
# next cell
request.result()      # raises ActionError with the viewer's reason
```

`viewer.debug_connection()["recent_events"]["command_errors"]` lists the most
recent refusals after the fact.

### `data_mapping`

The viewer reports its *data mapping* (always `kind: "anndata"`) — the
enabled embeddings and each dimension's default, the bound vector fields (each
`{"key": <obsm key>, "embedding": {"key", "columns"}}`), the neighbour graph,
the gene names, the expression matrix and the latent space.
It sends one after the dataset loads and another whenever the mapping changes,
whoever changed it -- a person, a session, or this notebook's `mapping_set`.
The first view never waits for gene names, so the first report is sent once
the gene names of the matrix the dataset opens with are read (right after
the first view); its `geneNames` is always `"index"` or `"column:<name>"`.

```{note}
The viewer sends this event only to a `cellucid` that lists `data_mapping`
among the events of `GET /_cellucid/protocol`, so the package never receives an
event it would reject.
```

Text frame:

```json
{
  "type": "data_mapping",
  "origin": "default",
  "mapping": {
    "kind": "anndata",
    "embeddings": {
      "enabled": [{"key": "X_umap", "columns": [0, 1]}],
      "defaults": {"1d": null, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": null}
    },
    "knn": null,
    "geneNames": "index",
    "expression": "X",
    "log1p": {},
    "highlyVariable": null,
    "latent": null,
    "vectorFields": []
  }
}
```

Hook payload: the record without `type` — `origin` and `mapping`.

- `origin` is one of `default`, `best_guess`, `session`, `user`, and
  `notebook` (this notebook's `viewer.mapping_set(...)`, the `mapping.set`
  action, whose own answer is its `action_result`).
- `mapping` must pass the same validation sessions use
  (`cellucid.data_mapping.validate_data_mapping`); an invalid one rejects
  the whole record before any hook runs. The validated canonical copy is what
  the hook receives and what `viewer.data_mapping` returns.

```python
@viewer.on_data_mapping
def handle(event):
    print(event["origin"], event["mapping"]["embeddings"])

viewer.data_mapping          # the latest mapping, or None before the first report
viewer.data_mapping_origin   # why it is what it is
```

### `history_changed`

The viewer tells every change of its history: a step recorded (by the
person, an agent or this notebook), the step shown moved (undo, redo, a jump,
another branch), the history's structure edited (a rename, a note, a
bookmark, a delete, a restore), a tracker started, used, renamed or removed,
or a session loaded.

Text frame:

```json
{
  "type": "history_changed",
  "kind": "step",
  "tracker": "t1",
  "head": "s5",
  "step": "s5",
  "branch": "b1"
}
```

Hook payload: the record without `type`.

- `kind` is one of `step`, `moved`, `edited`, `tracker` and `loaded`.
- `tracker` is the id of the tracker in use; `head` the step shown.
- `step` is the step recorded for `kind: "step"`, else `None`.
- `branch` is the id of the active branch, or `None` when the viewer names
  none.

```python
@viewer.on_history_change
def handle(event):
    if event["kind"] == "step":
        print("recorded", event["step"])

viewer.state.history          # the latest record, or None before the first
viewer.history.last_change    # the same, through LiveHistory
```

### `session_loaded`

The viewer's answer to `viewer.load_session(path)`, once the file's trackers
are in its history.

Text frame:

```json
{
  "type": "session_loaded",
  "requestId": "exact-request-id",
  "tracker": {"id": "t3", "name": "Annotation"},
  "head": "s7"
}
```

It completes the `LoadSessionRequest` of that `requestId`:
`request.result()` returns `{"tracker": ..., "head": ...}`. A refusal is a
`command_error` naming `load_session` and the request, and `result()` raises
`RuntimeError` with its reason. There is no decorator for it;
`viewer.register_hook("session_loaded", callback)` and `@viewer.on_message`
see it.

## Delivery and failure semantics

- A text frame holds at most 1 MiB of UTF-8 JSON; `selection`, `result`,
  `action_result` and `session_bundle` arrive only as binary frames.
- The channel was opened by a registered viewer with its exact token.
- The event name, full field set, and every payload value are validated before
  hook dispatch; a record that fails closes the channel (1008) and dispatches
  nothing.
- Events reach the hooks one at a time, in the order the viewer sent them.
- `@viewer.on_message` handlers run before the event-specific handlers.
- Each event is recorded once, in `viewer.state`, and wakes
  `viewer.wait_for_event(...)` before any hook runs.
- Every hook runs, even after one raises; the failures are logged (one
  `ExceptionGroup` when several fail) and the channel stays open. A failing
  hook never costs the notebook the event.

Hook callbacks run on the viewer channel's thread, one event at a time, so
the next event waits for the current hook. Keep them short and move expensive
analysis out of the callback path.

## Troubleshooting

If no event arrives:

1. Run `viewer.debug_connection()` in a new cell (a notebook widget's events
   arrive between cells, never while a cell runs).
2. Check `channel_connected`, `connected_outputs` and `recent_events`; a
   browser-tab viewer also reports `server_health`, `frontend_roundtrip` and
   `frontend_debug_snapshot`.
3. For a browser-tab viewer, check the browser network panel (filter `WS`)
   for the `/_cellucid/ws` connection and its frames.
4. Confirm the notebook still owns the displayed viewer; stop stale viewers
   with `viewer.stop()`.

See {doc}`14_troubleshooting_hooks` for the full diagnostic sequence.
