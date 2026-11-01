# Viewer classes (Jupyter objects you control from Python)

```{eval-rst}
.. currentmodule:: cellucid
```

Two classes share one viewer API:

- {class}`~cellucid.jupyter.Widget` (also `cellucid.Widget`): the notebook
  widget {func}`~cellucid.show_anndata` returns, shown in a cell's output and
  reached through the notebook's own connection;
- {class}`~cellucid.AnnDataViewer`: the same API for a viewer in a browser tab
  of its own, served from this machine (a script, an SSH tunnel).

They are designed for two-way interaction:
- **Viewer → Python:** hooks/events like `@viewer.on_selection`
- **Python → Viewer:** commands like `viewer.highlight_cells(...)`

If you prefer a minimal interface, start with {doc}`jupyter`.

---

## Fast path

### A) Any AnnData, in a notebook

```python
import cellucid

viewer = cellucid.Widget(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer  # the cell's output shows it (show_anndata() shows it for you)
```

### B) A Cellucid-ready AnnData (written by `prepare()`), in a browser tab

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer("my_study.h5ad")  # names itself (uns['cellucid'])
print(viewer.viewer_url)  # open it in a browser; in a notebook a link is shown
```

---

## Practical path (what you can do with a viewer)

### 1) Register hooks (viewer → Python)

```python
@viewer.on_selection
def on_selection(event):
    cells = event.get("cells", [])
    print("Selected", len(cells))

@viewer.on_ready
def on_ready(event):
    print("Viewer ready:", event)
```

### 2) Drive the UI from Python (Python → viewer)

```python
viewer.set_color_by("cell_type")
viewer.highlight_cells([0, 1, 2], color="#00ff00")
viewer.reset_view()
```

These calls return `None` without waiting for the viewer. If the viewer refuses
one, it reports a `command_error` back to Python, which Cellucid prints to
`sys.stderr` unless you take it over.

```python
@viewer.on_command_error
def on_command_error(event):
    print("Rejected:", event["command"], event["reason"])
```

### 3) Robust scripts: wait for readiness

A browser-tab viewer answers while your code runs, so it can wait:

```python
viewer.wait_for_ready(timeout=30)
```

(A notebook widget hears the viewer only between cells; there, use
`@viewer.on_ready` or read `viewer.state.ready` in a later cell — see
{doc}`jupyter`.)

### 4) Pull results, or capture a session bundle

```python
summary = viewer.pull(adata)                 # highlights, labels, scores -> adata.obs
request = viewer.request_session_bundle()
bundle = request.result(timeout=60)          # in a notebook: in the next cell
```

See {doc}`round_trip` for the round trip, and {doc}`sessions` for applying a
bundle back to AnnData.


### 5) Read and change the data mapping

```python
viewer.set_data_mapping(
    embeddings={
        "enabled": [{"key": "X_tsne", "columns": [0, 1]}],
        "defaults": {"1d": None, "2d": {"key": "X_tsne", "columns": [0, 1]}, "3d": None},
    },
)
```

Nothing waits. The viewer reports its mapping once the dataset loads (and
again whenever it changes): `set_data_mapping()` merges the fields into it,
checks the result (an invalid value raises, nothing sent) and sends it -- at
once, or, called before the viewer reported a mapping, as soon as it does.
The applied mapping is `viewer.data_mapping` from a later cell; a mapping the
viewer refuses arrives as a `command_error`.

```python
print(viewer.data_mapping["embeddings"])   # a later cell
```

---

## Lifecycle and cleanup (important)

Under the hood, a viewer starts a local server on a background thread of the
kernel: on port 8765, or the next free port when 8765 is taken (a printed note
says so).

Best practices:
- If you create many viewers in a long-running notebook, call `viewer.stop()` when done.
- If you repeatedly re-run a cell, you may end up with multiple servers unless the old viewer is stopped.

---

## API reference

### `Widget`

```{eval-rst}
.. autoclass:: cellucid.jupyter.Widget
   :members: wait_for_ready, display, debug_connection, stop
```

`Widget` has every method of `AnnDataViewer` below except the ones about its
server and web build (`viewer_url`, `server_url`, `clear_web_cache`,
`ensure_web_ui_cached`).

### `AnnDataViewer`

```{eval-rst}
.. autoclass:: AnnDataViewer
   :members:
```

### `SessionBundleRequest`

```{eval-rst}
.. autoclass:: cellucid.jupyter.SessionBundleRequest
   :members:
```

---

## Edge cases (do not skip)

### Not running in a notebook
- `show_anndata()` refuses outside a notebook and names the browser-tab
  alternatives; an `AnnDataViewer`'s `display()` prints its address.

### Re-running cells
- Each `AnnDataViewer` starts a local server on a port, and each `Widget`
  keeps its dataset (for an in-memory object, the store it serves and the
  files a cell-major matrix is laid out in) until it is stopped. Call
  `viewer.stop()` on viewers you no longer use.

### Waiting in a notebook cell
- A `Widget` hears the viewer only between cells: its waiting calls raise at
  once in the cell that would wait, and say what to do instead
  ({doc}`../../e_jupyter_hooks/03_supported_environments_matrix`).

### Callbacks and exceptions
- Hook callbacks run in registration order (`on_message` callbacks first), on
  the viewer channel's thread. An exception a callback raises is logged with
  its traceback (logger `cellucid.viewer_channel`) and never reaches the
  notebook. The event is recorded in `viewer.state` and handed to
  `wait_for_event` before any callback runs, and the callbacks after a failing
  one still run (several failures are logged together, as one
  `ExceptionGroup`). The channel stays open.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “`viewer.display()` prints a URL instead of showing a viewer”
Fix:
- That is an `AnnDataViewer`, which runs in a browser tab: open the URL. To
  show the viewer in a notebook cell, use `show_anndata()` or `cellucid.Widget`.

### Symptom: “Hooks never fire”
Fix:
- Check the viewer is ready: `viewer.state.ready` (in a later cell, for a
  widget), or `viewer.wait_for_ready()` for a browser-tab viewer.
- Temporarily register `@viewer.on_message` to see raw events and confirm communication.

### Symptom: “I have too many running servers”
Fix:
- Call `viewer.stop()` on old viewers.
- Restart the kernel if you’ve lost track of old viewer instances.

---

## See also

- {doc}`jupyter` for the function-level entry point (`show_anndata`)
- {doc}`server` for serving in a browser tab
