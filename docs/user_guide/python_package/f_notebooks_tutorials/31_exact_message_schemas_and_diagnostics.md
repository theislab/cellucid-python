# Exact message schemas and diagnostics

This tutorial shows how to inspect the closed Python↔viewer protocol without
inventing message types or accepting undeclared fields.

## Prerequisites

- A displayed viewer from `show(...)` or `show_anndata(...)`
- A viewer that has reached `ready`
- Familiarity with Python dictionaries

## The two authenticated directions

Python commands enter the iframe through `postMessage`. The embedding layer
adds the exact `viewerId` and `viewerToken`; user code supplies only the
command-specific record.

Viewer events return through `POST /_cellucid/events`. The server authenticates
the routing pair, removes the token, and gives the closed event record to the
Python viewer.

Neither direction accepts arbitrary message types.

## Inspect accepted events

`@viewer.on_message` sees every accepted event after routing fields have been
removed:

```python
received = []

@viewer.on_message
def record_current_event(event):
    received.append(event)
    print(event)
```

The `event` discriminator is one of:

- `ready`
- `selection`
- `hover`
- `click`
- `pong`
- `debug_snapshot`
- `session_bundle`
- `command_error`

The remaining fields are the exact payload for that discriminator. See
{doc}`../e_jupyter_hooks/07_frontend_to_python_events` for every field and
value constraint.

```{note}
`command_error` arrives only from the web build that emits `command_error`.
Older builds — including the one currently published — apply or refuse a
command silently, so this tutorial will not observe one.
```

## Exercise the request/response diagnostics

### Ping → pong

```python
import secrets

request_id = secrets.token_hex(8)
viewer.send_message({"type": "ping", "requestId": request_id})
pong = viewer.wait_for_event(
    "pong",
    timeout=5,
    predicate=lambda event: event["requestId"] == request_id,
)
assert pong["requestId"] == request_id
assert isinstance(pong["t"], int)
```

### Debug snapshot

```python
request_id = secrets.token_hex(8)
viewer.send_message({"type": "debug_snapshot", "requestId": request_id})
snapshot = viewer.wait_for_event(
    "debug_snapshot",
    timeout=5,
    predicate=lambda event: event["requestId"] == request_id,
)
assert snapshot["connected"] is True
print(snapshot["locationHref"])
print(snapshot["serverUrl"])
```

### Full diagnostic report

```python
report = viewer.debug_connection(timeout=5)
print(report["server_health"])
print(report["dataset_identity_probes"])
print(report["frontend_roundtrip"])
print(report["frontend_debug_snapshot"])
print(report["recent_events"])
```

The report uses only the two diagnostic exchanges above plus direct server and
verified-web-generation probes.

## Use the supported command helpers

Prefer the typed helpers for user actions:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
viewer.clear_highlights()
viewer.set_color_by("cell_type")
viewer.set_visibility([1, 2, 3], visible=False)
viewer.reset_view()
```

The low-level `send_message(...)` entry point is still closed. It accepts the
current command records used by those helpers and the internal diagnostics; an
unknown type or an undeclared field raises before notebook output is
published.

## Failure checks

### Unknown command

```python
viewer.send_message({"type": "my_custom_command"})
```

This raises `ValueError`. Cellucid does not publish or reinterpret the record.

### Undeclared command field

```python
viewer.send_message(
    {
        "type": "ping",
        "requestId": "probe",
        "extra": True,
    }
)
```

This raises before the command reaches the iframe.

### Undeclared event field

An authenticated event with an extra property fails before `@viewer.on_message`
or an event-specific hook runs. The failed event is not written to
`viewer.state` and does not wake `wait_for_event(...)`.

## Adding a new protocol variant

A new message is a coordinated contract change, not a runtime extension:

1. Define one exact field and value schema.
2. Implement and validate the web producer or consumer.
3. Add the matching Python command or event validator.
4. Add positive payload-preservation tests in both repositories.
5. Add negative tests for unknown, missing, and undeclared fields.
6. Update the protocol reference and user-facing example in the same change.

Until all six steps are present, the current builds reject the new type.

## Next steps

- {doc}`23_programmatic_highlighting_and_selection_callbacks`
- {doc}`32_session_persistence_and_restoring_analysis_artifacts`
- {doc}`../e_jupyter_hooks/06_python_to_frontend_commands`
- {doc}`../e_jupyter_hooks/07_frontend_to_python_events`
