# Viewer classes (Jupyter objects you control from Python)

```{eval-rst}
.. currentmodule:: cellucid
```

Viewer classes are the **stateful** notebook objects returned by {func}`~cellucid.show` and {func}`~cellucid.show_anndata`.

They are designed for two-way interaction:
- **Viewer → Python:** hooks/events like `@viewer.on_selection`
- **Python → Viewer:** commands like `viewer.highlight_cells(...)`

If you prefer a minimal interface, start with {doc}`jupyter`.

---

## Fast path

### A) Exported data (recommended for speed)

```python
from cellucid import CellucidViewer

viewer = CellucidViewer("./my_export")
viewer.display()
```

### B) AnnData mode (convenient for exploration)

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(adata)  # or "data.h5ad" / "data.zarr"
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

### 3) Robust notebooks: wait for readiness

```python
viewer.wait_for_ready(timeout=30)
```

### 4) Capture a session bundle (advanced)

```python
bundle = viewer.get_session_bundle(timeout=60)
```

See {doc}`sessions` for applying the bundle back to AnnData.

---

## Lifecycle and cleanup (important)

Under the hood, a viewer starts a local server process/thread.

Best practices:
- If you create many viewers in a long-running notebook, call `viewer.stop()` when done.
- If you repeatedly re-run a cell, you may end up with multiple servers unless the old viewer is stopped.

---

## API reference

### `CellucidViewer`

```{eval-rst}
.. autoclass:: CellucidViewer
   :members:
   :show-inheritance:
```

### `AnnDataViewer`

```{eval-rst}
.. autoclass:: AnnDataViewer
   :members:
   :show-inheritance:
```

---

## Edge cases (do not skip)

### Not running in a notebook
- If you instantiate a viewer in a plain Python script, `display()` will print a URL instead of embedding an iframe.

### Re-running cells creates multiple servers
- Each viewer starts a local server on a port.
- If you re-run the cell repeatedly without stopping old viewers, you may end up with multiple servers and confusing behavior.

### Callbacks and exceptions
- Hook callbacks are best-effort; exceptions are logged rather than raised to the notebook output.

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

- {doc}`jupyter` for function-level entry points (`show`, `show_anndata`)
- {doc}`server` for serving in a browser tab
