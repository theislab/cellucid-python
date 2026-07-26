# Writing robust callbacks

Notebook hooks are powerful, but they can also be the fastest way to create “weird notebook states” if callbacks are slow, crash, or mutate data unexpectedly.

This page gives battle-tested patterns for writing callbacks that behave well for:
- wet-lab users (copy/paste; predictable behavior),
- computational users (large data; heavy analysis),
- developers (threading, debugging, and failure isolation).

## Key rule (read this once)

Hook callbacks can run on the server’s request-handling thread.

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
Callback exceptions propagate through event delivery. The server logs the
exception and answers that event request with HTTP 500. Catch only an error
that the callback can handle completely; otherwise let it surface.
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

### Pattern D: prefer `wait_for_event` for deterministic workflows

If you want a predictable “do X after the next selection” flow, it can be clearer to block:

```python
viewer.wait_for_ready(timeout=60)
event = viewer.wait_for_event("selection", timeout=None)
cells = event["cells"]
```

This avoids callback threading concerns entirely.

See {doc}`09_viewer_state_and_wait_for_event`.

## Deep path: threading and AnnData pitfalls

### Backed `.h5ad` and thread safety

When you run
`show_anndata("data.h5ad", dataset_name="My dataset", dataset_id="my-dataset")`,
AnnData is opened read-only-backed.
Depending on your stack (h5py/hdf5), concurrent access from multiple threads can be unsafe.

Recommendations:
- avoid slicing backed AnnData inside hook callbacks,
- if you must, serialize access with a lock,
- or copy the needed arrays eagerly.

### Don’t block the server thread

If your callback blocks for seconds:
- events can queue up behind it,
- the UI may feel “laggy” or “stuck” from the user’s perspective.

Move long work out of the callback (queue pattern above).

## Logging: make failures visible

To see Cellucid hook/log output:

```python
import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger("cellucid.jupyter").setLevel(logging.INFO)
logging.getLogger("cellucid.server").setLevel(logging.INFO)
```

If you are debugging deeply, use `DEBUG`.

## Edge cases checklist

- `event["cells"]` can be empty (user cleared selection).
- `event["cell"]` can be `None` (not hovering).
- Use the documented fields for the exact event type and validate their native
  JSON-derived types before analysis.
- Very large selections can hit the `/_cellucid/events` size limit (see {doc}`07_frontend_to_python_events`).

## Troubleshooting

Symptoms and fixes:

- “My callback never runs” → check connectivity; run `viewer.debug_connection()`; see {doc}`14_troubleshooting_hooks`.
- “My callback request returns HTTP 500” → read the logged callback exception
  and fix it; the failure is not swallowed.
- “Notebook becomes sluggish” → your callback is doing too much; move work out of it.

## Next steps

- Viewer state + `wait_for_event`: {doc}`09_viewer_state_and_wait_for_event`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`
