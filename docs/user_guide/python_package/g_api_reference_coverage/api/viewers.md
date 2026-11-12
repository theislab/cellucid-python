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
- **Python → Viewer:** actions -- every action of the viewer is a method,
  like `viewer.highlight_cells(...)` -- each a step of the viewer's history

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

Every action of the viewer -- everything a person can do, and everything an
agent reads -- is a method, generated from the action catalog
({ref}`the full list <generated-action-methods>`):

```python
viewer.color_by(field={"source": "obs", "key": "cell_type"})
viewer.highlight_cells(cells=[0, 1, 2], color="#00ff00")
viewer.camera_reset(view="live")
fields = viewer.fields_list()            # a read
```

Each returns an {class}`~cellucid.jupyter.ActionRequest` without waiting for
the viewer; `request.result()` (in a notebook, in a later cell) is the
viewer's answer -- a read's result, or a change's outcome -- and raises
{class}`~cellucid.jupyter.ActionError` with the viewer's reason when it did
not apply the action. A change that made no history step (it was so already,
recording is paused, ...) is answered with why: `request.unrecorded`, and
`result()` warns {class}`~cellucid.jupyter.ActionNotRecordedWarning`. A
refusal nobody collects is printed to `sys.stderr` unless you take it over:

```python
@viewer.on_command_error
def on_command_error(event):
    print("Refused:", event["requestId"], event["reason"])
```

`viewer.act(type, **params)` sends any action named in data
(`viewer.act("color.by", field=...)`); {func}`cellucid.jupyter.action_catalog`
lists them all ({doc}`../../e_jupyter_hooks/06_python_to_frontend_commands`).

### 3) Robust scripts: wait for readiness

A browser-tab viewer answers while your code runs, so it can wait:

```python
viewer.wait_for_ready(timeout=30)
```

(A notebook widget hears the viewer only between cells; there, use
`@viewer.on_ready` or read `viewer.state.ready` in a later cell — see
{doc}`jupyter`.)

### 4) Pull results, or capture a session

```python
summary = viewer.pull(adata)                 # highlights, labels, scores -> adata.obs
request = viewer.request_session()
session = request.result(timeout=60)         # in a notebook: in the next cell
plan = session.apply(adata)                  # what applying the step shown would change
```

See {doc}`round_trip` for the round trip, and {doc}`sessions` for applying a
session back to AnnData.

### 5) The history: turns, changes, sessions loaded

```python
with viewer.turn("Mark the T cells"):            # one turn of the history, undone together
    viewer.color_by(field={"source": "obs", "key": "cell_type"})
    viewer.highlight_cells(cells=adata.obs["cell_type"] == "T cell", label="T cells")

@viewer.on_history_change
def on_history_change(event):
    print(event["kind"], event["head"])

viewer.history.last_change                       # the latest change, or None
viewer.history_undo()                            # the history's actions are methods too
loading = viewer.load_session("analysis.cellucid-session")
```

A turn's actions carry its id (`python-` and 12 hexadecimal digits) and
label; turns do not nest. `viewer.history` is a
{class}`~cellucid.jupyter.LiveHistory`: the latest `history_changed` event
and a hook on each; the history is changed through the generated `history_*`
methods. `viewer.load_session(...)` adds a session file's trackers beside the
person's own and returns a {class}`~cellucid.jupyter.LoadSessionRequest`
(`loading.result()` in a later cell:
`{"tracker": {"id", "name"}, "head"}`).

### 6) Read and change the data mapping

```python
viewer.mapping_set(
    embeddings={
        "enabled": [{"key": "X_tsne", "columns": [0, 1]}],
        "defaults": {"1d": None, "2d": {"key": "X_tsne", "columns": [0, 1]}, "3d": None},
    },
)
```

Nothing waits. `mapping_set()` (the `mapping.set` action) changes the groups
it names and leaves the others as they are (`knn=None` is "no neighbor
graph"; leaving `knn` out keeps the one shown). The viewer reports its
mapping once the dataset loads and again whenever it changes: the applied
mapping is `viewer.data_mapping` from a later cell; a mapping the viewer
refuses fails the request with the reason.

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

(generated-action-methods)=
### Action methods (generated)

Every action of the catalog, as a method of both classes. They are
generated from the catalog the package ships (`scripts/generate_actions.py`,
never written by hand), so this list is exactly the viewer's actions.

```{eval-rst}
.. autoclass:: cellucid.jupyter.ViewerActions
   :members:
```

### `ActionRequest`, `ActionError`, `ActionNotRecordedWarning` and `OMIT`

```{eval-rst}
.. autoclass:: cellucid.jupyter.ActionRequest
   :members:

.. autoexception:: cellucid.jupyter.ActionError

.. autoexception:: cellucid.jupyter.ActionNotRecordedWarning

.. autodata:: cellucid.jupyter.OMIT
   :annotation:

.. autofunction:: cellucid.jupyter.action_catalog
```

### `SessionRequest`

```{eval-rst}
.. autoclass:: cellucid.jupyter.SessionRequest
   :members:
```

### `LiveHistory` and `LoadSessionRequest`

```{eval-rst}
.. autoclass:: cellucid.jupyter.LiveHistory
   :members:

.. autoclass:: cellucid.jupyter.LoadSessionRequest
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
