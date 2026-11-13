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

The viewer protocol is **schema 5** (`VIEWER_PROTOCOL_SCHEMA` in the web app's
`data/viewer-protocol.js`; schema 5 added the agent messages). Ten commands,
exactly the keys of their type
(`_COMMAND_MESSAGE_FIELDS` and `_BINARY_COMMAND_HEADER_FIELDS` in
`jupyter/_wire.py`), validated and serialized by `_require_frontend_message`:

| Command | Keys | Answered by |
|---|---|---|
| `action` | `type`, `requestId`, `action`, `params`, `turn`, `turnLabel` | `action_result` or `command_error` (same `requestId`) |
| `load_session` (binary) | `type`, `requestId`, `name`; the body is the session file | `session_loaded` or `command_error` (same `requestId`) |
| `requestSessionBundle` | `type`, `requestId` | `session_bundle` or `command_error` (same `requestId`) |
| `ping` | `type`, `requestId` | `pong` |
| `debug_snapshot` | `type`, `requestId` | `debug_snapshot` |
| `freeze` | `type` | nothing (the host stopped) |
| `agent_list` | `type`, `requestId`, `connection` | `agent_tools` (binary: the MCP tools, the same whatever the person's settings) |
| `agent_call` | `type`, `requestId`, `tool`, `args`, `agent`, `connection` | exactly one `agent_result` (binary: the MCP tool result, an error included), with `agent_asking` / `agent_answered` while the person is asked |
| `agent_cancel` | `type`, `requestId` | the call's `agent_result` (stopped: a change rolled back, but a playing camera path stopped where it is with its step kept; a question withdrawn) |
| `agent_answer` | `type`, `requestId`, `allowed` | the person answered the card's question in the agent's app (Allow once or Deny) |

The four `agent_*` commands carry a local agent connection (`cellucid mcp`,
`connection` `{id, client}`): the page runs each call as the agent under the
person's permissions and answers with the agent events; a package delivers
those events to the viewer's agent port (`_agent_link.py`), never to the
notebook's hooks ({doc}`../../web_app/p_developer_docs/31_local_agent_connection`).

- **`action` is the one command that changes or reads the viewer's state.**
  `action` names one action of the catalog by its dotted type (`color.by`,
  `highlight.cells`, `history.get`, …) and `params` are its parameters; the
  viewer checks them with the action's schema and dispatches the action with
  origin `host` and the host's name (`python`, `jupyter`, `r`), so it takes
  the path a person's control takes and is one step of the history, labeled
  as the notebook's. Schema 2's state commands -- `setColorBy`,
  `setVisibility`, `highlight`, `clearHighlights`, `resetCamera`,
  `setDataMapping`, `setGeneSet`, `setField` -- are gone: they are the
  actions `color.by`, `cells.visibility.set`, `highlight.cells`,
  `highlight.clear-all`, `camera.reset`, `mapping.set`, `gene-set.save` and
  `field.push`.
- `turn` and `turnLabel` name the notebook turn an action belongs to: both
  exact text (the turn's id and what it was asked to do) or both `null`
  outside a turn. `with viewer.turn("QC"):` (Python, ids `python-<12 hex>`)
  and `cellucid_turn(viewer, "QC", {...})` (R, ids `r-<12 hex>`) set them on
  every action sent inside; the viewer passes them to the dispatcher, so the
  steps are one turn of the history, shown as one group in the State panel
  and undone together, as an agent's turn is. Turns do not nest.
- an `action` whose `params` hold arrays is one binary frame instead: a
  little-endian `uint32` header length, the JSON header
  (`_BINARY_COMMAND_HEADER_FIELDS`: `type`, `requestId`, `action`, `params`
  without the arrays, `turn`, `turnLabel` and `arrays`), padded with spaces
  to a multiple of 4 bytes, then each array's elements, little-endian, in the order `arrays`
  lists them. Each `arrays` entry is `{path, type, length}`: `path` the keys
  (and list indices, for a list of objects) leading to the array in `params`,
  ending in a key `params` does not hold; `type` one of `uint8`, `uint16`,
  `uint32`, `float32`, `float64`. Python takes every NumPy array out of the
  parameters in the order met (`_split_action_params`) and writes each
  straight into the frame (`_binary_command_frame`); the viewer puts each back
  at its path, read in place when it starts on a boundary of its element size
  (`decodeViewerBinaryCommand` in `data/viewer-protocol.js`), checked as
  strictly as the server checks a selection. Cells (`x-cells` in the action's
  schema) are `uint32` original positions; one value per cell (`x-values`)
  any value type.
- `load_session` pushes a `.cellucid-session` file into the viewer
  (`viewer.load_session(path)`, R `cellucid_load_session()`): one binary
  frame whose JSON header `{type, requestId, name}` (`name` is the file's
  name, without a folder) is padded the same way, and whose body is the file
  (`_load_session_frame`). The viewer loads it as the person's Load does, its
  trackers beside the person's own, while it shows the host's dataset, and
  answers `session_loaded`, or `command_error` with why not (another
  dataset's session, a damaged file).
- sent by `cellucid._viewer_channel.send_viewer_command`: on the open socket,
  or into the viewer's outbox until it connects, flushed in order on
  `welcome`; a frame that failed to go out whole stays queued, so nothing is
  sent twice

### The action catalog, one document in three repositories

The actions are defined once, in the web app (`assets/js/app/actions/`), and
the packages are generated from them: the web repository's
`scripts/generate-action-catalog.mjs` writes `catalogDocument` of every
definition as canonical JSON (sorted keys) to its
`tests/fixtures/action-catalog.json` (its home), and
`node scripts/sync-shared-files.mjs` there copies it to this package's
`src/cellucid/action_catalog.json` and to cellucid-r's
`inst/action_catalog.json`, byte for byte; each package's suite compares its
copy with the web repository's file (`tests/test_action_catalog.py`,
`tests/testthat/test-actions.R`; see
{ref}`the shared files <shared-files>`).
`scripts/generate_actions.py` then writes `jupyter/_actions_generated.py`
(one method per action, typed, with the action's title, summary and
parameter descriptions as its docstring); `tests/test_action_catalog.py`
fails while that module differs from what the script writes, so an action
added in the web app cannot be missing from Python (or R) unnoticed.
`_actions.py` checks and converts every action's parameters by the
catalog, the same way for every action.

### Frontend → Python (events)

- JSON text frames, at most 1 MiB each, for `hover`, `click`, `ready`, `pong`,
  `debug_snapshot`, `command_error`, `data_mapping`, `color_by`,
  `result_removed`, `results`, `history_changed` and `session_loaded`
- binary frames for `selection`, `result`, `action_result` and
  `session_bundle`: a little-endian `uint32` header length, a UTF-8 JSON
  header (`{"type": "selection", "source"}`, a result's header,
  `{"type": "action_result", "requestId", "step", "unrecorded"}`, or
  `{"type": "session_bundle", "requestId"}`), then the body
- an `action_result`'s body is the action's outcome as UTF-8 JSON (`null`
  for nothing to report, a read's result otherwise): a binary frame, because
  a read's result is as large as it asks for (a picture of a view, the values
  of many cells) and a text event holds at most 1 MiB; `step` is the history
  step the action made, `{id, label}`, or `null`, and `unrecorded` why a
  change made none (`null` beside a step and for a read)
- a selection's body is the selected row indices as little-endian `uint32`,
  strictly ascending, one or more; the server reads it into one buffer
  (unmasked in place, no second copy), checks it
  (`_require_selection_cell_array`), and the viewer's validator decodes it
  into the `list[int]` its hooks receive
- a session's body is streamed to a temporary file, read and checked whole
  (`cellucid.read_session`: session format 4), and turned into the
  server-built `session_bundle` event

Events the viewer keeps sending on its own -- `color_by`, `data_mapping`,
`result` / `result_removed` / `results`, `selection`, `hover`, `click`,
`history_changed` -- mirror the viewer's state into `viewer.state`,
`viewer.results` and `viewer.history` as it changes, whoever changed it.
`history_changed` tells where the history stands after each change (the web
app's `app/viewer-host/history-events.js`); the history itself is read as a
session (`viewer.request_session()`). An `action_result` answers the request
that sent its action, so `data_mapping` and `result` carry no `requestId`.

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
- `cellucid/assets/js/app/viewer-host/command-handler.js` (web): dispatches a
  host's `action` as origin `host`
- the shared rules: `cellucid/tests/fixtures/viewer-protocol.json` (schema 5,
  with one valid frame of every kind, binary actions and a binary
  `load_session` as exact bytes), kept byte for byte in cellucid-python
  (`tests/fixtures/viewer_protocol.json`) and cellucid-r, each copy compared
  with the web repository's file by its own suite; this package's encoder
  makes exactly the bytes of every example binary command
- the session files a host receives and pushes: session format 4, tested in
  the three repositories against one golden history the web app's own model
  makes (the web repository's `scripts/generate-history-golden.mjs`: three
  trackers, one replayed and one removed, one paused, four branches, a
  notebook turn and an agent turn, a group, a squashed step of what changed
  while recording was paused, notes, a bookmark, a deleted step and its
  journal entry, a stale result, each step's parameters version), written as
  `tests/fixtures/history-golden.cellucid-session` and
  `history-golden.expected.json` there, as
  `tests/fixtures/history_golden/history_golden.cellucid-session` and
  `history_golden.expected.json` here and in cellucid-r's
  `tests/testthat/fixtures/history_golden/`; the expected file holds the
  reading and the apply plans every reader must find; the web suite
  (`tests/history-golden-fixture.test.mjs`) holds the files to what the
  generator makes, and this package (`tests/test_sessions.py`) and cellucid-r
  (`tests/testthat/test-sessions.R`) compare their copies with the web
  repository's, byte for byte

### Event validation in Python

Python accepts exactly (`_INBOUND_EVENT_FIELDS`):

- `selection`
- `hover`
- `click`
- `ready`
- `pong`
- `debug_snapshot`
- `action_result`
- `command_error`
- `data_mapping`
- `color_by`
- `result`, `result_removed` and `results`
- `history_changed`
- `session_loaded`
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

### `action_result` (binary)

Meaning:
- the viewer's answer to one `action`.

Expected keys (after Python decoded the body, `_wire._decode_action_result`):
- `requestId`: the action's request
- `step`: `{id, label}` of the history step it made, or `null` (a read, or a
  change that made none)
- `unrecorded`: why a change made no step -- `"unchanged"` (it was so
  already), `"retracted"` (it put back the step just before it, which was
  taken out), `"off"` (the History settings leave this kind of change out),
  `"paused"` (recording is paused), `"elsewhere"` (another tracker was shown
  meanwhile), `"no-history"` (the dataset's history has not started) -- or
  `null` beside a step and for an action that is no change (a read, a move of
  the history); any other value, or a reason beside a step, is refused
  (`_wire._ACTION_UNRECORDED`, the shared fixture's `action.unrecorded`)
- `outcome`: the body's JSON: a change's outcome (counts, names, ids) or a
  read's result; finite numbers only

It completes the `ActionRequest` registered under `requestId`
(`AnnDataViewer._action_waiters`); an answer nobody waits for any more is
dropped.

### `session_bundle` (Python-side event)

Meaning:
- the server received a `.cellucid-session` file (format 4) and is reporting
  status back to Python.

Expected keys:
- `status`: exactly `"ok"`
- `bytes`: non-negative int
- `path`: non-empty temporary file path
- `requestId`: str (correlates with `requestSessionBundle`)

Note:
- this event is generated by the server, not directly by the frontend;
- an invalid file never dispatches an error-shaped `session_bundle`: it is
  reported as a `command_error` naming `requestSessionBundle` and its
  `requestId`.

### `command_error`

Meaning:
- a command Python sent was refused by the viewer instead of carried out.

Expected keys:
- `command`: exactly one of the six command types this protocol defines
- `requestId`: the failed command's `requestId`, or `null` for one that names
  none (`freeze`, or a frame the viewer could not read as far as its id)
- `reason`: one line; the viewer sends at most 1000 characters, Python
  accepts at most 4000

The event completes the request it names: an action's `result()` raises
`ActionError` with `reason`. `AnnDataViewer` also prints it to `sys.stderr`
when no `command_error` hook is registered, naming the method that sent the
action (`viewer.color_by(...)`), because nothing has to collect a request.

### `data_mapping`

Meaning:
- the viewer resolved or changed its data mapping (which arrays it plots and
  reads), whoever changed it.

Expected keys:
- `origin`: one of `default`, `best_guess`, `session`, `user`, `notebook`
  (`notebook`: this host's `mapping.set` action)
- `mapping`: a complete data mapping, validated by
  `cellucid.data_mapping.validate_data_mapping` and stored as its canonical copy

### `history_changed`

Meaning:
- the viewer's history changed: a step recorded (by the person, an agent or
  a host), the step shown moved (undo, redo, a jump, another branch), the
  history's structure edited (a rename, a note, a bookmark, a delete, a
  restore), a tracker started, used, renamed or removed, or a session loaded.
  A conversation kept with a tracker is not told.

Expected keys:
- `kind`: one of `step`, `moved`, `edited`, `tracker`, `loaded`
- `tracker`: the id of the tracker in use
- `head`: the step shown
- `step`: the step the change made or concerns, or `null`
- `branch`: the id of the active branch, or `null`

`viewer.state.history` keeps the latest; `viewer.history` reads it
(`last_change`, `tracker`, `head`, `branch`) and `@viewer.on_history_change`
runs on each. R: `cellucid_last_history_change(viewer)` and
`cellucid_on(viewer, "history_changed", handler)`.

### `session_loaded`

Meaning:
- the viewer loaded the session file of one `load_session`.

Expected keys:
- `requestId`: the `load_session` it answers
- `tracker`: `{id, name}` of the tracker shown
- `head`: its step shown

It completes the `LoadSessionRequest` registered under `requestId`
(`request.result()` returns `{"tracker": {"id", "name"}, "head"}`); a
refused load is a `command_error` naming `load_session`.

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
4) Add the type to the viewer protocol (the web app's
   `assets/js/data/viewer-protocol.js`) and run
   `node scripts/sync-shared-files.mjs` in the web repository: it writes
   `tests/fixtures/viewer-protocol.json` there and copies it byte for byte to
   cellucid-python's `tests/fixtures/viewer_protocol.json` and to cellucid-r.
   Add the R package's validator with it.
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
  "commands": ["action", "agent_answer", "agent_call", "agent_cancel", "agent_list",
               "debug_snapshot", "freeze", "load_session", "ping", "requestSessionBundle"],
  "binaryCommands": ["action", "load_session"],
  "events": ["action_result", "agent_answered", "agent_asking", "agent_result", "agent_tools",
             "click", "color_by", "command_error", "data_mapping", "debug_snapshot",
             "history_changed", "hover", "pong", "ready", "result", "result_removed",
             "results", "selection", "session_bundle", "session_loaded"]
}
```

`binaryCommands` names the commands that may be sent as binary frames (`action`, `load_session`); a viewer refuses a
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

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-10-26 | The capability document's key set is fixed (`channel`, `commands`, `binaryCommands`, `events`). `binaryCommands` was the one key added, when commands that name cells became binary frames | A viewer validates payloads against exact key sets, so a later capability arrives inside these lists; with the added key an older viewer refuses the host rather than misread its commands | A new top-level key for each capability |
| by 2026-11-12 | The state commands of the earlier schema (`setColorBy` and the others) became actions and are no longer the protocol's: a `command_error` naming one is refused as any unknown command is, in Python and in R, and a `data_mapping` event names no request | A change of the viewer's state goes through the action catalog (above) | Two command sets for one change |
