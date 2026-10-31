# Writing robust callbacks

Notebook hooks are powerful, but they can also be the fastest way to create “weird notebook states” if callbacks are slow, crash, or mutate data unexpectedly.

This page gives battle-tested patterns for writing callbacks that behave well for:
- wet-lab users (copy/paste; predictable behavior),
- computational users (large data; heavy analysis),
- developers (threading, debugging, and failure isolation).

## Key rule (read this once)

Hook callbacks run on the viewer channel's thread, one event at a time, in the
order the viewer sent them.

Practical meaning:
- keep callbacks fast and defensive,
- avoid doing heavy analysis directly inside the callback,
- treat hook payloads as **untrusted** and validate keys and types.

## Fast path: safe callback template (copy/paste)

```python
import logging

logger = logging.getLogger("cellucid.hooks")

@viewer.on_selection
def on_selection(event):
    cells = event["cells"]
    if type(cells) is not list or any(type(index) is not int for index in cells):
        raise TypeError("selection cells must be a list of native integers")
    logger.info("Selection: %d cells", len(cells))
```

```{note}
A callback that raises does not take the channel down: the server logs the
exception (logger `cellucid.viewer_channel`), the hooks registered after it
still run for that event, and the next event is delivered. The event itself
was recorded in `viewer.state` and handed to `wait_for_event(...)` before any
hook ran, so a failing callback never costs the notebook the event. Catch only
an error that the callback can handle completely; otherwise let it surface in
the log.
```

## Practical patterns

### Pattern A: “Do the minimum in the callback, do the work elsewhere”

Use a queue:

```python
from queue import Queue

selection_queue: Queue[list[int]] = Queue()

@viewer.on_selection
def enqueue_selection(event):
    cells = event["cells"]
    if cells:
        selection_queue.put(cells)
```

Then in another cell (or later in the notebook), consume:

```python
cells = selection_queue.get()  # blocks
subset = adata[cells].copy()
subset
```

Why this is robust:
- callbacks stay fast,
- heavy work happens in normal notebook execution order,
- you control debouncing and cancellation.

### Pattern B: debouncing selection events

Some UI workflows may emit multiple selection events during an interaction.
You can debounce in Python:

```python
import time

last_at = 0.0

@viewer.on_selection
def on_selection_debounced(event):
    global last_at
    now = time.time()
    if now - last_at < 0.25:  # 250ms
        return
    last_at = now
    print("Selection:", len(event["cells"]))
```

### Pattern C: treat hover as a debounced state stream

Hover is throttled and can drop events under load. Avoid heavy work in hover hooks:

```python
@viewer.on_hover
def on_hover(event):
    cell = event["cell"]
    if cell is None:
        return
    # Good: lightweight UI feedback / logging
    # Bad: expensive plotting, model inference, large AnnData slicing
```

### Pattern D: read the latest result in the next cell

If you want a predictable “do X with what I selected” flow, skip the callback:
select in the viewer, then run a cell that reads it.

```python
selection = viewer.selection      # the last selection, or None
cells = selection.cells
```

This avoids callback threading concerns entirely. (`viewer.wait_for_event(...)`
blocks until the next event instead; it works for a browser-tab viewer, but a
notebook widget hears the viewer only between cells, so there it raises at
once.)

See {doc}`09_viewer_state_and_wait_for_event`.

## Deep path: threading and AnnData pitfalls

### Your `adata` and thread safety

`show_anndata(...)` does not open your data in Python while the viewer runs:
the viewer reads the file or folder itself, and the chunks of an in-memory
AnnData or MuData are made from a snapshot of its arrays when the viewer asks
for them (its expression matrices laid out gene-major in the background, the
arrays they are read from read-only until then). The `adata` your callbacks touch is
therefore only yours — but callbacks run on the viewer channel's thread, not
your notebook's, and
if you opened that object backed (`anndata.read_h5ad(..., backed="r")`),
concurrent access from several threads can be unsafe depending on your
h5py/HDF5 stack.

Recommendations:
- avoid slicing backed AnnData inside hook callbacks,
- if you must, serialize access with a lock,
- or copy the needed arrays eagerly.

### Don’t block the channel thread

If your callback blocks for seconds:
- events queue up behind it (they are delivered one at a time, in order),
- a callback that waits for a later event of the same viewer
  (`viewer.wait_for_event(...)`, a session request's `result(...)`) waits until
  its timeout, because that event is delivered only after the callback returns,
- the UI may feel “laggy” or “stuck” from the user’s perspective.

Move long work out of the callback (queue pattern above).

## Logging: make failures visible

To see Cellucid hook/log output:

```python
import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger("cellucid.jupyter").setLevel(logging.INFO)
logging.getLogger("cellucid.server").setLevel(logging.INFO)
logging.getLogger("cellucid.viewer_channel").setLevel(logging.INFO)
```

If you are debugging deeply, use `DEBUG`.

## Edge cases checklist

- `event["cells"]` holds at least one cell; clearing a selection sends no event.
- `event["cell"]` can be `None` (not hovering).
- Use the documented fields for the exact event type and validate their native
  JSON-derived types before analysis.
- Very large selections arrive whole (a binary frame, 4 bytes a cell); the
  callback receives a Python list, so a list of millions of ints is itself
  large (see {doc}`07_frontend_to_python_events`).

## Troubleshooting

Symptoms and fixes:

- “My callback never runs” → check connectivity; run `viewer.debug_connection()`; see {doc}`14_troubleshooting_hooks`.
- “My callback fails” → read the logged callback exception (logger
  `cellucid.viewer_channel`) and fix it; the failure is not swallowed.
- “Notebook becomes sluggish” → your callback is doing too much; move work out of it.

## Next steps

- Viewer state + `wait_for_event`: {doc}`09_viewer_state_and_wait_for_event`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`
