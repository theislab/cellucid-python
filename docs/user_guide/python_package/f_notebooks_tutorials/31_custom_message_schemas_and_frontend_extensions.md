# Custom message schemas and frontend extensions

This tutorial is for developers and “power notebook” users who want to go beyond:
- basic selection hooks
- basic highlighting/color-by commands

It explains how Cellucid’s Python↔frontend messaging works, and what “custom messages” actually means in practice.

```{important}
If you are using the standard hosted Cellucid web app UI (served via the hosted-asset proxy), you can only use message types that the UI already understands.

To *add new message types*, you must modify the Cellucid web app (repo: `cellucid`) or contribute upstream.
```

---

## At a glance

**Audience**
- Developers (primary)
- Advanced computational users (who need custom integrations)

**Time**
- ~20–45 minutes

**Prerequisites**
- A notebook environment
- A working viewer (`show(...)` or `show_anndata(...)`)
- Comfort reading JSON payloads and debugging request/response flows

---

## The two directions (and why they are not symmetric)

### Python → frontend (commands)

Python sends messages into the iframe via `postMessage(...)`.

In Python, the low-level API is:

```python
viewer.send_message({"type": "...", "any": "payload"})
```

Cellucid also provides convenience wrappers:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
viewer.clear_highlights()
viewer.set_color_by("leiden")
viewer.set_visibility([1, 2, 3], visible=False)
viewer.reset_view()
```

### Frontend → Python (events)

The iframe sends events back to Python via HTTP POST to:
- `POST /_cellucid/events`

Cellucid routes events to the correct Python `viewer` object using `viewerId`.

---

## Built-in “request/response” patterns you can reuse

Some message types follow a request/response pattern using `requestId`.

### Pattern: ping → pong (connectivity probe)

```python
import secrets

req_id = secrets.token_hex(8)
viewer.send_message({"type": "ping", "requestId": req_id})
pong = viewer.wait_for_event("pong", timeout=5, predicate=lambda e: e.get("requestId") == req_id)
pong
```

### Pattern: debug_snapshot (frontend environment introspection)

```python
import secrets

req_id = secrets.token_hex(8)
viewer.send_message({"type": "debug_snapshot", "requestId": req_id})
snap = viewer.wait_for_event("debug_snapshot", timeout=5, predicate=lambda e: e.get("requestId") == req_id)
snap
```

```{tip}
If you are debugging “it works locally but not on JupyterHub”, start with:

`viewer.debug_connection()`
```

---

## Receiving custom events in Python

### Option A: catch everything with `on_message`

```python
@viewer.on_message
def _all_messages(event):
    print(event)
```

This receives:
- known events (selection/hover/click/ready)
- debug events (pong, debug_snapshot, console)
- any future/custom types the frontend posts to `/_cellucid/events`

### Option B: register a handler for a specific custom event type

If the frontend posts an event with `type: "my_custom_event"`, you can register:

```python
def handle_custom(event):
    print("custom:", event)

viewer.register_hook("my_custom_event", handle_custom)
```

```{note}
The `@viewer.on_selection` / `@viewer.on_ready` decorators are just convenience wrappers. `register_hook(...)` works for any event name.
```

---

## Designing a custom message schema (recommended conventions)

If you control the frontend and want to add new message types, follow these conventions:

1) **Namespace your `type`**
   - Example: `x_my_lab_plugin:request_stats`
2) **Use `requestId` for request/response**
   - This makes concurrent requests safe.
3) **Keep payload JSON-serializable**
   - avoid numpy arrays; send small summaries or references
4) **Treat inputs as untrusted**
   - validate shapes, bounds, and expected keys

Example request message:

```python
import secrets

request_id = secrets.token_hex(8)
viewer.send_message({
    "type": "x_example:request_summary",
    "requestId": request_id,
    "params": {"top_n": 10},
})
```

Example response event (posted by the frontend to `/_cellucid/events`):

```json
{
  "viewerId": "<viewerId>",
  "type": "x_example:summary",
  "requestId": "<requestId>",
  "status": "ok",
  "result": {"top_n": 10, "items": ["..."]}
}
```

In Python you would wait for it:

```python
resp = viewer.wait_for_event(
    "x_example:summary",
    timeout=30,
    predicate=lambda e: e.get("requestId") == request_id,
)
resp
```

---

## Frontend extension checklist (high level)

If you want to add a new message type end-to-end, you need:

1) A frontend handler for `postMessage` commands
2) A way to send responses/events back to Python via `POST /_cellucid/events`
3) A stable message schema + versioning/compat policy

For the detailed protocol and design rationale, see:
- `cellucid/markdown/HOOKS_DEVELOPMENT.md`

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “My custom message does nothing”

Likely causes:
- the standard web app UI doesn’t recognize your `type`
- you sent the message before the viewer displayed

Fix:
- confirm viewer is displayed (`show(...)` auto-displays in notebooks)
- verify your message type is supported by the UI you are using
- use `@viewer.on_message` to see what the frontend sends back (if anything)

### Symptom: “I receive some events but not my custom ones”

Likely causes:
- your frontend never POSTs to `/_cellucid/events`
- requestId mismatches (you’re waiting for the wrong response)

Fix:
- use `viewer.debug_connection()` to check server health + recent events
- log/print incoming raw messages with `on_message`

---

## Next steps

- {doc}`23_programmatic_highlighting_and_selection_callbacks` (robust hooks in practice)
- {doc}`32_session_persistence_and_restoring_analysis_artifacts` (request/response pattern via session bundles)
