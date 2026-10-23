# Viewer classes (Jupyter objects you control from Python)

```{eval-rst}
.. currentmodule:: cellucid
```

The viewer class is the **stateful** notebook object returned by {func}`~cellucid.show_anndata`.

They are designed for two-way interaction:
- **Viewer → Python:** hooks/events like `@viewer.on_selection`
- **Python → Viewer:** commands like `viewer.highlight_cells(...)`

If you prefer a minimal interface, start with {doc}`jupyter`.

---

## Fast path

### A) Any AnnData

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer.display()
```

### B) A Cellucid-ready AnnData (written by `prepare()`)

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer("my_study.h5ad")  # names itself (uns['cellucid'])
viewer.display()
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

### 3) Robust notebooks: wait for readiness

```python
viewer.wait_for_ready(timeout=30)
```

### 4) Capture a session bundle (advanced)

```python
bundle = viewer.get_session_bundle(timeout=60)
```

See {doc}`sessions` for applying the bundle back to AnnData.

### 5) Read and change the data mapping

```python
viewer.wait_for_event("data_mapping")
print(viewer.data_mapping["embeddings"])
viewer.set_data_mapping(
    embeddings={"1d": None, "2d": {"key": "X_tsne", "columns": [0, 1]}, "3d": None},
)
```

The viewer reports its mapping once the dataset loads (and again whenever it
changes); until then `viewer.data_mapping` is `None` and `set_data_mapping()`
raises before sending anything.

---

## Lifecycle and cleanup (important)

Under the hood, a viewer starts a local server process/thread.

Best practices:
- If you create many viewers in a long-running notebook, call `viewer.stop()` when done.
- If you repeatedly re-run a cell, you may end up with multiple servers unless the old viewer is stopped.

---

## API reference

### `AnnDataViewer`

```{eval-rst}
.. autoclass:: AnnDataViewer
   :members:
```

---

## Edge cases (do not skip)

### Not running in a notebook
- If you instantiate a viewer in a plain Python script, `display()` will print a URL instead of embedding an iframe.

### Re-running cells creates multiple servers
- Each viewer starts a local server on a port.
- If you re-run the cell repeatedly without stopping old viewers, you may end up with multiple servers and confusing behavior.

### Callbacks and exceptions
- Hook callbacks run in registration order. A callback exception propagates
  through event delivery and is reported by the event request.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “`viewer.display()` prints a URL instead of showing an embedded viewer”
Fix:
- You are likely not in a notebook context; open the printed URL in a browser tab, or run inside Jupyter/VSCode/Colab.

### Symptom: “Hooks never fire”
Fix:
- Wait for the viewer to be ready: `viewer.wait_for_ready()`.
- Temporarily register `@viewer.on_message` to see raw events and confirm communication.

### Symptom: “I have too many running servers”
Fix:
- Call `viewer.stop()` on old viewers.
- Restart the kernel if you’ve lost track of old viewer instances.

---

## See also

- {doc}`jupyter` for the function-level entry point (`show_anndata`)
- {doc}`server` for serving in a browser tab
