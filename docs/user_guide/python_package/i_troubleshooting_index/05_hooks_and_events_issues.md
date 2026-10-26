# Hooks and events issues

## A callback never runs

Confirm that the callback is registered on the same live viewer object shown in
the notebook, and that `viewer.debug_connection()["channel_connected"]` is
`True`. Reproduce one event, then inspect the `/_cellucid/ws` WebSocket in the
browser's Network panel (filter `WS`); preserve its close code and reason if it
closed, and any exception logged by `cellucid.viewer_channel`.

## A command targets the wrong viewer

Send it from the viewer object that owns the displayed frame. This is
especially important when a notebook contains multiple viewers or a cell was
re-run: each viewer has its own channel.

## Session-bundle request rejected

A session bundle answers exactly one pending `get_session_bundle()` request and
is single-use; one that arrives after the request stopped waiting is dropped.
Start a new request (with a longer `timeout=`) instead.

## A command returned without changing the viewer

Commands are one-way: `viewer.set_color_by(...)` returns `None` whether the
viewer applied it or refused it. A refusal comes back as a `command_error`
event, printed to `sys.stderr` and listed by
`viewer.debug_connection()["recent_events"]["command_errors"]`. Handle it in
code with `@viewer.on_command_error`.

## Event payload rejected

Treat event and command schemas as exact contracts. Compare the emitted event
name and payload fields with {doc}`../e_jupyter_hooks/16_reference`.

For a complete symptom catalog, see
{doc}`../e_jupyter_hooks/14_troubleshooting_hooks`.
