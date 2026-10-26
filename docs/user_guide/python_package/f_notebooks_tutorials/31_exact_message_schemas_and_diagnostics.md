# Exact message schemas and diagnostics

This tutorial shows how to inspect the closed Python↔viewer protocol without
inventing message types or accepting undeclared fields.

## Prerequisites

- A displayed viewer from `show_anndata(...)`
- A viewer that has reached `ready`
- Familiarity with Python dictionaries

## The two authenticated directions

Both directions travel on the viewer channel, one WebSocket the viewer opens
to the server and authenticates once with its exact `viewerId` and
`viewerToken`. Python's commands are the command-specific record only (a
highlight or a visibility change travels as one binary frame of cell indices);
the viewer's events come back as JSON text frames, and selections as binary
cell indices. The server checks each closed event record and gives it to the
Python viewer that owns the channel.

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
- `data_mapping`

The remaining fields are the exact payload for that discriminator. See
{doc}`../e_jupyter_hooks/07_frontend_to_python_events` for every field and
value constraint.

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
print(report["anndata_store_probes"])  # the AnnData store the viewer reads
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
unknown type or an undeclared field raises before anything is sent.

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

This raises `TypeError` before the command reaches the viewer.

### Undeclared event field

An authenticated event with an extra property fails before `@viewer.on_message`
or an event-specific hook runs. The failed event is not written to
`viewer.state` and does not wake `wait_for_event(...)`, and the server closes
that viewer's channel (WebSocket close code 1008, policy violation).

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
