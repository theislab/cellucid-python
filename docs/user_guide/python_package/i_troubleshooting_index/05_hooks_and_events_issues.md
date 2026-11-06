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

## Session request rejected

A session answers exactly one pending request and is single-use; one that
arrives after the request was cancelled or stopped waiting is dropped.
Start a new request instead: `request = viewer.request_session()`, and
`request.result()` -- in a notebook in the next cell, since the widget hears
the viewer only between cells.

## A session file pushed into the viewer did not load

`viewer.load_session(path)` returns its request at once. Collect it --
`request.result()`, in a notebook in a later cell -- and it raises
`RuntimeError` with the viewer's reason when the viewer did not load the
file: a session saved on another dataset, a damaged file, or a viewer that
shows a dataset this notebook did not open. A missing path raises
`FileNotFoundError` before anything is sent.

## A turn raised

`with viewer.turn(label):` needs a label (what the turn does, in a few
words), and turns do not nest: opening one inside another raises
`RuntimeError` naming the turn already open. Close the first, or put the
actions in one turn.

## `viewer.history` says nothing yet

`viewer.history.last_change` is `None` until the viewer tells its first
change of the history (`history_changed`); in a notebook that arrives
between cells. Read it in a later cell, or run code on each change with
`@viewer.on_history_change`.

## An action returned without changing the viewer

An action returns its request at once: `viewer.color_by(...)` does not wait
for the viewer. Collect it -- `request.result()`, in a notebook in a later
cell -- and it raises `ActionError` with the viewer's reason when the viewer
refused it. A refusal also comes back as a `command_error` event, printed to
`sys.stderr` and listed by
`viewer.debug_connection()["recent_events"]["command_errors"]`. Handle it in
code with `@viewer.on_command_error`.

## Event payload rejected

Treat event and command schemas as exact contracts. Compare the emitted event
name and payload fields with {doc}`../e_jupyter_hooks/16_reference`.

For a complete symptom catalog, see
{doc}`../e_jupyter_hooks/14_troubleshooting_hooks`.
