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
- optional keys (safe to add later)
- payload size expectations (avoid megabytes)

### Step A2 — Implement frontend emission (web app)

In the web app (`cellucid/` repo), emit:

```js
// pseudo-code
POST /_cellucid/events
{
  type: "your_event_name",
  viewerId: "<viewer id>",
  ...payload
}
```

### Step A3 — Decide if Python needs explicit mapping

Python already maps:
- `"selection"`, `"hover"`, `"click"`, `"ready"`

For any other `type`:
- the Python hook name becomes the raw `type`.

If you want a convenience decorator `@viewer.on_<event>`:
- add a new property on `BaseViewer` similar to `on_selection`.

### Step A4 — Update Python event handling if needed

File:
- `cellucid-python/src/cellucid/jupyter.py`

Possible changes:
- extend `event_map` in `BaseViewer._handle_frontend_message`
- extend `ViewerState` if you want it stored as a “latest snapshot”

### Step A5 — Document + test

Docs:
- add schema to {doc}`11_hooks_events_protocol_and_schema`
- add usage example to {doc}`../e_jupyter_hooks/index` if user-facing

Tests:
- create a minimal viewer instance (or test the hook registry in isolation)
- call `_handle_frontend_message({...})` and assert:
  - hooks fire
  - state updates

---

## Part B — Add a new Python → frontend command

### Step B1 — Decide command name and parameters

Pick:
- `type`: a stable command string (e.g. `"setThreshold"`)

Define:
- required parameters
- default values (if any)
- behavior when called before `ready`

### Step B2 — Implement frontend handler (web app)

The frontend must listen for postMessage commands and execute them.

Important:
- commands include `viewerId` + `viewerToken`
- the frontend should ignore commands missing/invalid token (basic robustness)

### Step B3 — Implement a Python wrapper method (optional but recommended)

File:
- `cellucid-python/src/cellucid/jupyter.py`

Add a method on the viewer, e.g.:

- `def set_threshold(self, value: float): self.send_message({"type": "setThreshold", "value": value})`

### Step B4 — Handle notebook lifecycle edge cases

Decide:
- should the command queue until `ready`?
- should it error if viewer is not displayed?
- what happens if the iframe is gone?

At minimum:
- keep error messages actionable (“call viewer.display() first”).

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
1) confirm the browser is POSTing `/_cellucid/events` (DevTools Network tab)
2) confirm the payload includes the correct `viewerId`
3) confirm the Python viewer is alive and hasn’t been stopped/GC’d

Use:
- `viewer.debug_connection()` (see {doc}`12_debugging_playbook`)

### Symptom: “My new command is sent but UI doesn’t change”

Checklist:
1) confirm the JS postMessage handler exists in the frontend
2) confirm the message type matches exactly
3) confirm the frontend is not rejecting the command due to missing/invalid `viewerToken`
4) confirm timing (viewer not ready yet)
