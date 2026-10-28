# Viewer state and `wait_for_event`

Cellucid offers two complementary ways to consume viewer → Python information:

1. **Reactive**: callbacks like `@viewer.on_selection`
2. **Pull / synchronous**: inspect `viewer.state` and block with `viewer.wait_for_event(...)`

This page focuses on the second style, which is often easier for notebook workflows where you want deterministic “run this cell after the next selection” behavior.

## `viewer.state`: latest-event snapshot

`viewer.state` is a small thread-safe container with:
- `ready`, `selection`, `hover`, `click`, `command_error`, `data_mapping`: the
  *latest payload* for those event types (or `None`)
- `last_event_type`, `last_event`, `last_updated_at`: generic “what happened most recently”

Example:

```python
print(viewer.state.ready)       # None until the dataset has opened

print("Last selection:", viewer.state.selection)
print("Last hover:", viewer.state.hover)
print("Last click:", viewer.state.click)
print("Last rejected command:", viewer.state.command_error)
```

When to use `viewer.state`:
- you want “the latest thing the user did”,
- you want to poll/inspect without blocking,
- you want to debug whether events are arriving at all.

## `viewer.wait_for_event(...)`: block until the *next* event arrives

```{important}
A notebook widget (what `show_anndata()` returns) hears the viewer only
**between cells**: nothing the viewer says can arrive while the cell that
waits is running. So in a notebook `wait_for_event()` and `wait_for_ready()`
(unless the viewer is ready already) raise at once and say what to do
instead: register a hook, or read `viewer.state` in a later cell. They block
as described here for a browser-tab viewer ({class}`~cellucid.AnnDataViewer`),
and on a thread of your own in a notebook (the viewer answers it as soon as
the cell ends).
```

Signature (simplified):

```text
viewer.wait_for_event(event: str, timeout: float | None = 30.0, *, predicate=None) -> dict
```

Behavior:
- It waits for the *next* event of that type after the call begins.
- If `timeout` is `None`, it waits forever.
- If `predicate` is provided, it filters payloads until one matches.

### Example: block until the next selection

```python
event = viewer.wait_for_event("selection", timeout=None)
cells = event["cells"]
print("Selected:", len(cells))
```

### Example: wait until the user changes the data mapping

```python
event = viewer.wait_for_event(
    "data_mapping",
    timeout=None,
    predicate=lambda e: e["origin"] == "user",
)
```

### Convenience: `viewer.wait_for_ready(...)`

```python
viewer.wait_for_ready(timeout=60)
```

If the viewer has already sent a `ready` event, this returns immediately using `viewer.state.ready`.

## How it works (short, accurate mental model)

Internally each `wait_for_event(...)` call registers before it waits and is
handed its event the moment the event is recorded, so a burst of other events
(thousands of hovers) cannot make it miss the one it waits for.

Implications:
- `wait_for_event(...)` is reliable for “next event” workflows.
- It is not an event log: it sees only events that arrive after the call
  begins. (`viewer.debug_connection()` lists the most recent events.)

## Edge cases and footguns

### “I called `wait_for_event("selection")` but I already selected cells”

That’s expected: `wait_for_event(...)` waits for the next event after the call begins.

Options:
- inspect `viewer.state.selection` for the last selection, or
- trigger a new selection in the UI.

### Timeouts

If `timeout` is too small and you don’t interact quickly, you’ll get `TimeoutError`.
Use `timeout=None` when you want an interactive “wait until the user does it” cell.

### Multiple viewers / stale outputs

If you have multiple viewers, make sure you call `wait_for_event` on the correct one.
If you restarted the kernel but didn’t re-run the viewer cell, the output is stale
(it keeps its last view and says the notebook stopped it).

Fix: re-run the cell that creates `viewer`.

## Troubleshooting

If waiting never returns:
- run `viewer.debug_connection()` to confirm events are arriving
- see {doc}`14_troubleshooting_hooks`

## Next steps

- Event payloads: {doc}`07_frontend_to_python_events`
- Session bundles: {doc}`10_session_bundles_request_session_bundle`
