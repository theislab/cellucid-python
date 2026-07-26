# Hooks and events issues

## A callback never runs

Confirm that the callback is registered on the same live viewer object shown in
the notebook. Reproduce one event, then inspect requests to
`/_cellucid/events`; preserve the first non-success response and its body.

## A command targets the wrong viewer

Use the viewer ID returned by the object that owns the displayed frame. This is
especially important when a notebook contains multiple viewers or a cell was
re-run.

## Session-bundle request rejected

Session-bundle uploads require the exact registered viewer/request pair and are
single-use. Start a new request instead of replaying an expired or consumed one.

## Event payload rejected

Treat event and command schemas as exact contracts. Compare the emitted event
name and payload fields with {doc}`../e_jupyter_hooks/16_reference`.

For a complete symptom catalog, see
{doc}`../e_jupyter_hooks/14_troubleshooting_hooks`.
