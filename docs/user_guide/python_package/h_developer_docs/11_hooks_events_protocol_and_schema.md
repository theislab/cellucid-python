# Hooks/events protocol and schema

This page specifies the **bidirectional communication protocol** between:

- the embedded Cellucid web UI (in the browser), and
- Python code in `cellucid-python` (in the notebook/kernel).

It is the “developer contract” behind:
- `@viewer.on_selection` and friends (frontend → Python),
- `viewer.highlight_cells(...)` and friends (Python → frontend).

If you only need user-facing usage examples, start with:
- {doc}`../e_jupyter_hooks/index`

---

## One channel, both directions

Transport: one WebSocket per viewer, which the viewer page opens at
`<server URL>/_cellucid/ws` ({doc}`../e_jupyter_hooks/05_architecture_viewer_channel`).

Opening it:
- `GET /_cellucid/protocol` names the channel; the viewer refuses a server
  without one and asks for a package upgrade;
- the upgrade passes the `Host` check (`421`) and the `Origin` check (`403`);
- the first frame is `{"type": "hello", "viewerId", "viewerToken"}`; the
  server compares the token in constant time and answers
  `{"type": "welcome"}`, or closes with `4401`;
- from then on the socket is that one viewer's, so no frame carries
  credentials, and a newer channel of the same viewer replaces it (`4409`).

### Python → frontend (commands)

- one JSON text frame per command, exactly the keys of its type
  (`_COMMAND_MESSAGE_FIELDS` in `jupyter/_wire.py`), validated and serialized
  by `_require_frontend_message`
- a command that names cells (`highlight`, `setVisibility`) is one binary
  frame instead: a little-endian `uint32` header length, the JSON header
  (`_BINARY_COMMAND_HEADER_FIELDS`: `{"type": "highlight", "color"}` or
  `{"type": "setVisibility", "visible", "all"}`) padded with spaces to a
  multiple of 4 bytes, then the cells as little-endian `uint32`, strictly
  ascending (none for `setVisibility` with `"all": true`). The cells are
  checked and sorted once (`_require_cell_indices`, vectorized for a NumPy
  array) and written straight into the frame; the padding puts them on a
  4-byte boundary, where the viewer reads them in place
  (`decodeViewerBinaryCommand` in `data/viewer-protocol.js`), checked as
  strictly as the server checks a selection
- sent by `cellucid._viewer_channel.send_viewer_command`: on the open socket,
  or into the viewer's outbox until it connects, flushed in order on
  `welcome`; a frame that failed to go out whole stays queued, so nothing is
  sent twice
- used for “do something now” commands (highlight, set color-by, reset camera,
  set the data mapping, request a session bundle)

### Frontend → Python (events)

- JSON text frames, at most 1 MiB each, for `hover`, `click`, `ready`, `pong`,
  `debug_snapshot`, `command_error`, `data_mapping`, `color_by`,
  `result_removed` and `results`
- binary frames for `selection`, `result` and `session_bundle`: a
  little-endian `uint32` header length, a UTF-8 JSON header
  (`{"type": "selection", "source"}`, a result's header, or
  `{"type": "session_bundle", "requestId"}`), then the body
- a selection's body is the selected row indices as little-endian `uint32`,
  strictly ascending, one or more; the server reads it into one buffer
  (unmasked in place, no second copy), checks it
  (`_require_selection_cell_array`), and the viewer's validator decodes it
  into the `list[int]` its hooks receive
- a session bundle's body is streamed to a temporary file, validated, and
  turned into the server-built `session_bundle` event

A selection is binary because a JSON array of indices costs about 8 bytes per
cell to encode, parse and validate; four bytes per cell, sorted, lets both
sides check uniqueness in one pass without a hash set.

Delivery:
- each event is validated (`_require_inbound_jupyter_event`) and handed to the
  viewer the channel belongs to, on the channel's thread, one at a time and in
  order;
- a frame that breaks the rules closes the channel with `1008` before any hook
  runs; a hook that raises is logged and the channel stays open.

Implementation:
- `cellucid-python/src/cellucid/_viewer_channel.py` (`accept_viewer_channel`,
  `send_viewer_command`, `register_event_callback`, `route_event`) over
  `_websocket.py` (RFC 6455 framing)
- `cellucid-python/src/cellucid/jupyter/_viewer.py`
  (`_register_viewer_for_messages`, `AnnDataViewer._handle_frontend_message`)
- `cellucid/assets/js/data/viewer-channel.js` and `viewer-protocol.js` (web)
- the shared rules: `cellucid/tests/fixtures/viewer-protocol.json`, kept byte
  for byte in cellucid-python (`tests/fixtures/viewer_protocol.json`) and
  cellucid-r, each copy SHA-pinned

### Event validation in Python

The viewer accepts exactly:

- `selection`
- `hover`
- `click`
- `ready`
- `pong`
- `debug_snapshot`
- `command_error`
- `data_mapping`
- the server-generated `session_bundle`

Each type owns one exact required field set. Unknown types, missing fields,
undeclared fields (a `viewerId` among them) and invalid values raise before
hook dispatch. Adding an event therefore requires an explicit coordinated
schema change in both repositories, **in that order** — see “Adding a new
event type” below.

---

## Canonical event schemas (what Python code should expect)

Each event is a native JSON object with an exact non-empty `type` and the
type's fields, and nothing else; `AnnDataViewer` validates the complete record
before callbacks receive it without `type`.

Below are the expected shapes for core events.

### `ready`

Meaning:
- the viewer finished initial dataset load and is interactive.

Required payload keys:
- `n_cells`: non-negative native integer
- `dimensions`: exactly `1`, `2`, or `3`

### `selection`

Meaning:
- the user selected a set of cells (lasso, click+shift, neighbor selection, annotation, etc.).

Expected keys:
- `cells`: list[int] (0-based indices in dataset row order)
- `source`: exact non-empty selection-method token

### `hover`

Meaning:
- the user is hovering a single cell.

Expected keys:
- `cell`: non-negative int or `null`
- `position`: exact finite `{x, y, z}` object or `null`

### `click`

Meaning:
- the user clicked a single cell.

Expected keys:
- `cell`: non-negative int
- `button`: `0`, `1`, or `2`
- `shift`: bool
- `ctrl`: bool

### `session_bundle` (Python-side event)

Meaning:
- the server received a `.cellucid-session` bundle frame and is reporting status back to Python.

Expected keys:
- `status`: exactly `"ok"`
- `bytes`: non-negative int
- `path`: non-empty temporary file path
- `requestId`: str (correlates with `requestSessionBundle`)

Note:
- this event is generated by the server, not directly by the frontend;
- an invalid bundle never dispatches an error-shaped `session_bundle`: it is
  reported as a `command_error` naming `requestSessionBundle`.

### `command_error`

Meaning:
- a command Python sent was rejected by the viewer instead of applied.

Expected keys:
- `command`: exactly one of the ten command types this protocol defines
- `reason`: one line, at most 500 characters

Commands carry no `requestId`, so this event correlates with a command *type*
and not with an individual call. Adding a `requestId` to the commands is not a
compatible change: the viewer validates each command's key set exactly, so a
new field would make every command fail on an already-deployed web build.

`AnnDataViewer` prints this event to `sys.stderr` when no `command_error` hook is
registered, because the command methods are one-way — `set_color_by` has
already returned `None` by the time the rejection arrives on the channel's thread,
and nothing can make it raise without blocking every command on a browser round
trip.

### `data_mapping`

Meaning:
- the viewer resolved or changed its data mapping (which arrays it plots and
  reads).

Expected keys:
- `origin`: one of `default`, `best_guess`, `session`, `user`, `notebook`
- `requestId`: a string exactly when `origin` is `notebook` (the event answers
  that `setDataMapping` command); `null` otherwise
- `mapping`: a complete data mapping, validated by
  `cellucid.data_mapping.validate_data_mapping` and stored as its canonical copy

`setDataMapping` is the one command that carries a `requestId`: it is new, so
no deployed build validates a key set for it, and `viewer.set_data_mapping()`
needs the answer to tell an applied mapping from a refused one. Its key set is
exactly `type`, `requestId`, `mapping`.

---

## Hook registration semantics (Python)

Python provides a hook registry with two useful patterns:

### Decorator style

```python
@viewer.on_selection
def handle(event):
    print(event["cells"])
```

### “Catch-all” message hook

```python
@viewer.on_message
def handle_any(event):
    # event includes "event": <event_name> plus the payload keys
    print(event)
```

Internals:
- `HookRegistry.trigger` always calls all `message` handlers for any event.

“Catch-all” means *every declared event*, not *any event*. There is no
custom or app-specific event channel: `_require_inbound_jupyter_event` resolves
the event `type` against `_INBOUND_EVENT_FIELDS` and raises on anything else, so
an undeclared type never reaches a hook. Widening a hook is not how a new event
arrives — the checklist below is, and step 2 is the only place an event type
comes into existence.

---

## Adding a new event type (developer checklist)

If you add a new event type:

1) Define one closed field/type schema.
2) Add the schema to `_INBOUND_EVENT_FIELDS` in
   `cellucid-python/src/cellucid/jupyter/_wire.py`, and **release it**. That
   mapping is what `_require_inbound_jupyter_event` validates against and what
   `GET /_cellucid/protocol` publishes, so one entry both accepts the event and
   announces that it is accepted.
3) Only then add the frontend producer and its exact validation, gated on
   `GET /_cellucid/protocol` naming the type in its `events` list. The gate is
   positive evidence rather than an assumption: a build that emits an event the
   notebook does not accept produces a worse failure than not emitting at all.
4) Add the type to the shared `viewer-protocol.json` (the web app's
   `tests/fixtures/viewer-protocol.json`, copied byte for byte to
   cellucid-python's `tests/fixtures/viewer_protocol.json` and to cellucid-r,
   each copy SHA-pinned), and the R package's validator with it.
5) Add positive preservation tests and negative missing/undeclared-field tests.
6) Update docs:
   - add the schema to this page,
   - add user-facing examples in {doc}`../e_jupyter_hooks/index` if appropriate.

### Why the reader ships first

`_require_inbound_jupyter_event` raises `ValueError` on an event type it does
not know, and the channel treats that as a protocol violation: it closes with
`1008`, and the viewer freezes and reports it. A producer that ships ahead of
its validator therefore converts a silent problem into a broken viewer.

The gate in step 3 is the hard part, because the notebook and the web build are
upgraded independently — the Python package fetches the web UI from
`https://www.cellucid.com`, so a notebook running an older `cellucid` can be
driving a newer web build the same day it deploys.

So the viewer asks before it sends: every new event a build could emit needs
a channel the notebook answers with what it accepts. That channel is
{ref}`GET /_cellucid/protocol <python-protocol-capability-endpoint>`, which
the viewer reads before it opens its WebSocket. Because the route is answered
from the same declaration the validator reads, “the route lists it” is
exactly equivalent to “the validator accepts it”. Inferring the capability
from the health `version` string instead would be a proxy rather than
evidence, and would depend on release numbering that is decided after the
code is written.

The document names the channel and holds three lists, and no version number,
so that it can never be sniffed the way a version would be:

```json
{
  "channel": "_cellucid/ws",
  "commands": ["clearHighlights", "debug_snapshot", "freeze", "highlight",
               "ping", "requestSessionBundle", "resetCamera", "setColorBy",
               "setDataMapping", "setField", "setGeneSet", "setVisibility"],
  "binaryCommands": ["highlight", "setField", "setVisibility"],
  "events": ["click", "color_by", "command_error", "data_mapping",
             "debug_snapshot", "hover", "pong", "ready", "result",
             "result_removed", "results", "selection", "session_bundle"]
}
```

`binaryCommands` names the commands sent as binary frames; a viewer refuses a
host whose list disagrees with how it reads a command it knows, or that lists
none at all (a `cellucid` from before binary commands), and asks for an
upgrade. `events` is every type `_require_inbound_jupyter_event` accepts, and
`commands` is every command this installation may send — which is also exactly the set of
names a `command_error` may report as having failed, so a build whose own
command list has grown ahead of the notebook can avoid naming one the notebook
would reject.

The web build fetches this route during `ViewerHostDataSource.initialize()`
and sends an event only when the notebook lists its type. A notebook too old
to name a channel (no `channel` key, or no route at all) is refused at once:
the viewer asks for a package upgrade instead of opening anything.

The gate covers new event types, not a new transport: a host without the
channel is refused with the upgrade message.

---

## Edge cases and safety considerations

### Huge selections

Selections can be large (millions of indices). On the wire they are 4 bytes
per cell; in the hook payload `event["cells"]` is a Python `list[int]`, about
36 bytes per cell.

Developer guidance:
- treat `event["cells"]` as untrusted and potentially huge,
- avoid copying it multiple times,
- consider downsampling for logging and plots.

### Viewer lifecycle and “stale events”

Each viewer has its own channel:
- if you re-run a cell and create a new viewer, the old viewer may still be
  alive, with its own page and channel,
- ensure you stop old viewers (`viewer.stop()`) to avoid confusion; that
  freezes their pages and closes their channels.

### Server exposure

If the server is reachable by others on a network, the per-viewer token and
the `Origin` check still protect the viewer channel, but the service itself
remains an exposed attack surface.

Recommendation:
- do not expose servers publicly; prefer SSH tunnels or trusted internal networks.

---

## Troubleshooting

### Symptom: “I never receive `ready`”

Likely causes:
- dataset load failed in the frontend (check browser console),
- a cell is still running (a notebook widget's messages arrive between cells),
- the viewer's channel did not open (`channel_connected` is false; for a
  browser tab, a proxy that does not pass WebSockets or a server the page
  cannot reach).

Fix:
- run `viewer.debug_connection()` and follow {doc}`12_debugging_playbook`.

### Symptom: “Events arrive, but `cells` indices are wrong”

Cause:
- indices are positional; you changed row order between “what the UI is showing” and “what Python is indexing”.

Fix:
- keep the same AnnData row order as the served dataset,
- or use stable identifiers and explicit mapping for cross-dataset transfer.
