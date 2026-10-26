# Extension point: add a new hook event or command

This page is a step-by-step guide for extending Cellucid’s **bidirectional interaction layer**:

- add a new **frontend → Python** event (hook), or
- add a new **Python → frontend** command.

Prerequisites:
- understand notebook embedding: {doc}`10_jupyter_embedding_architecture`
- read the protocol spec: {doc}`11_hooks_events_protocol_and_schema`

---

## Part A — Add a new frontend → Python event

### Step A1 — Define the event name and payload schema

Pick:
- `type`: a short, stable string (e.g. `"annotation_created"`)

Define:
- required keys (minimum needed to be useful)
- schema-optional keys and their exact meaning
- payload size expectations (avoid megabytes)

### Step A2 — Implement frontend emission (web app)

In the web app (`cellucid/` repo), emit:

```js
// pseudo-code: one JSON text frame on the viewer channel
channel.sendEvent({
  type: "your_event_name",
  ...payload
});
```

Add the type to `VIEWER_EVENT_FIELDS` in `assets/js/data/viewer-protocol.js`,
regenerate `tests/fixtures/viewer-protocol.json`, and copy it to
cellucid-python and cellucid-r with the new digest.

### Step A3 — Add the Python schema

Add the new type and its complete required-field set to `_INBOUND_EVENT_FIELDS`
in `cellucid-python/src/cellucid/jupyter/_wire.py`, and implement its exact
value checks in `_require_inbound_jupyter_event`. Until this exists, Python
rejects the event before hook dispatch.

That one entry does two jobs: it is what the validator accepts, and it is what
`GET /_cellucid/protocol` publishes. Nothing else needs updating for the web
build to discover the new capability, and nothing else may be updated instead —
a list kept anywhere but here can announce an event the validator would reject.

**Release the Python side before the emitter.** The web build must read
`/_cellucid/protocol` and find the type in its `events` list before it sends
one; a notebook running an older `cellucid` treats an unknown event type as a
protocol violation and closes the channel with `1008`, and the viewer freezes,
which is worse than not emitting. The reasoning is in “Why the reader ships first” in
{doc}`11_hooks_events_protocol_and_schema`, and the route itself is
{ref}`documented here <python-protocol-capability-endpoint>`.

If the event needs a convenience decorator, add a property on `AnnDataViewer`
similar to `on_selection`.

### Step A4 — Update Python event handling if needed

Files:
- `cellucid-python/src/cellucid/jupyter/_wire.py`
- `cellucid-python/src/cellucid/jupyter/_hooks.py`

Required changes:
- extend the closed inbound validator in `_wire.py`
- extend `ViewerState` in `_hooks.py` only if the event owns a latest-state
  snapshot

### Step A5 — Document + test

Docs:
- add schema to {doc}`11_hooks_events_protocol_and_schema`
- add usage example to {doc}`../e_jupyter_hooks/index` if user-facing

Tests:
- create a minimal viewer instance (or test the hook registry in isolation)
- call `_handle_frontend_message({...})` and assert:
  - hooks fire
  - state updates
  - missing and undeclared fields fail before either action

---

## Part B — Add a new Python → frontend command

### Step B1 — Decide command name and parameters

Pick:
- `type`: a stable command string (e.g. `"setThreshold"`)

Define:
- the complete required parameter set
- exact value constraints
- behavior when called before `ready`

### Step B2 — Implement frontend handler (web app)

The frontend runs each command it receives on the viewer channel
(`ViewerHostDataSource._handleCommand` in `assets/js/data/viewer-host-source.js`).

Important:
- commands carry no credentials: the channel was authenticated when it opened
- the frontend checks each command's exact key set and reports a failure as a
  `command_error`
- add the type to `VIEWER_COMMAND_FIELDS` in `viewer-protocol.js` and to
  `_COMMAND_MESSAGE_FIELDS` in `_wire.py` (and `.viewer_command_fields` in
  cellucid-r's `R/viewer-protocol.R`), and regenerate the shared
  `viewer-protocol.json`
- a command that names cells is a binary frame instead: its header goes in
  `VIEWER_BINARY_COMMAND_HEADER_FIELDS` / `_BINARY_COMMAND_HEADER_FIELDS` /
  `.viewer_binary_command_fields`, its cells in the body
  (`decodeViewerBinaryCommand` reads it), never a JSON list of cells

### Step B3 — Implement a Python wrapper method (optional but recommended)

Files:
- `cellucid-python/src/cellucid/jupyter/_wire.py`
- `cellucid-python/src/cellucid/jupyter/_viewer.py`

Add the command to `_require_frontend_message` in `_wire.py`, then add a
method on `AnnDataViewer`, e.g.:

- `def set_threshold(self, value: float): self.send_message({"type": "setThreshold", "value": value})`

### Step B4 — Handle notebook lifecycle edge cases

Commands already wait for a viewer that has not connected yet, and arrive in
order. Decide:
- whether the command needs the dataset open (the viewer runs commands only
  once its handlers are in place),
- what it means when the viewer refuses it (a `command_error`).

At minimum:
- keep error messages actionable.

### Step B5 — Document + test

Docs:
- add the command to the “Python → frontend commands” documentation

Tests:
- Python-side tests can only validate that the correct message dict is produced.
  The actual UI effect must be tested in the web app.

---

## Troubleshooting

### Symptom: “My new event never arrives in Python”

Checklist:
1) confirm the event frame appears on the `/_cellucid/ws` connection
   (DevTools Network tab, filter `WS`)
2) confirm `GET /_cellucid/protocol` lists the event (the viewer sends only
   listed types)
3) confirm the Python viewer is alive and hasn’t been stopped/GC’d

Use:
- `viewer.debug_connection()` (see {doc}`12_debugging_playbook`)

### Symptom: “My new command is sent but UI doesn’t change”

Checklist:
1) confirm the frontend handles the command type
2) confirm the message type matches exactly
3) check `viewer.debug_connection()["recent_events"]["command_errors"]` for a
   refusal
4) confirm the viewer is connected (`channel_connected`)
