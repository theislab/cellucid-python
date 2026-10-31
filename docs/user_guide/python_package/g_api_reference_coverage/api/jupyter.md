# Jupyter (notebook embedding + hooks)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents the notebook-facing APIs:
- **Notebook view:** {func}`~cellucid.show_anndata` (an AnnData or MuData object, or any file or folder Cellucid reads: `.h5ad`, a Zarr store, `.h5mu`, `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds`/`.rda`, or a file written by {func}`~cellucid.prepare`)
- **Interactive control:** the widget it returns ({class}`~cellucid.jupyter.Widget`) supports hooks + Python→frontend commands

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: copy/paste the “Fast path” and use the troubleshooting section when something looks wrong.
- Computational: use “Practical path” + “Edge cases” to avoid performance and data-shape pitfalls.
- Power user: use “Deep path” to understand how the server + hooks + session capture work.

**Prerequisites**
- `pip install cellucid`
- If you use AnnData: install `anndata` (and typically `numpy`, `scipy`, `pandas`).
- A notebook that shows widgets: JupyterLab, Jupyter Notebook 7, VS Code, Google Colab, marimo or Databricks ({doc}`../../e_jupyter_hooks/03_supported_environments_matrix`).

---

## Fast path (beginner-friendly)

### A) “I have an AnnData”

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

What you should see:
- A Cellucid viewer in the notebook output, showing your embedding, once the
  cell has finished (the viewer reads the data through the kernel between
  cells).

### B) “I have a file written by `prepare()`”

```python
from cellucid import show_anndata

viewer = show_anndata("my_study.h5ad")  # it names itself (uns['cellucid'])
```

---

## Practical path (step-by-step + common workflows)

### 1) Notebook display vs browser tab

Cellucid always runs in a browser context (the viewer UI), but you choose **where it appears**:
- **Notebook output**: {func}`~cellucid.show_anndata` (this page), a widget.
- **Browser tab**: {class}`~cellucid.AnnDataViewer` (the same API), or
  {func}`~cellucid.serve_anndata` / `cellucid serve …` (see {doc}`server`).

### 2) Opening what you have, or a prepared file

Open the object or file you have when:
- you’re iterating quickly in analysis notebooks,
- you want the viewer to read the `.h5ad` or Zarr store you already have, in
  place.

Open a file written by {func}`~cellucid.prepare` when:
- you want reproducibility (the file is a concrete, self-naming artifact),
- you want speed (a gene-major matrix: every gene is one small read),
- you want to share the dataset or host it.

### 3) A minimal “hooks” example (viewer → Python)

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)

@viewer.on_selection
def on_selection(event):
    # event["cells"] is a list[int] of selected indices (row indices into adata)
    print("Selected", len(event.get("cells", [])), "cells")
```

### 4) A minimal “control” example (Python → viewer)

```python
# Highlight a few cells (indices into the dataset)
viewer.highlight_cells([0, 10, 20], color="#ff0000")

# Color points by an obs field (must exist in the dataset/AnnData)
viewer.set_color_by("cell_type")

# Hide a subset (or pass None to apply to all)
viewer.set_visibility([1, 2, 3], visible=False)

# Reset the camera
viewer.reset_view()
```

### 5) Running code once the viewer is ready (robust notebooks)

The viewer loads after the cell that showed it finishes, so register a hook:

```python
@viewer.on_ready
def ready(event):
    print("Viewer ready:", event)
```

or read `viewer.state.ready` in a later cell. (`viewer.wait_for_ready()`
returns at once when it is ready, and raises in a cell that would have to
wait: the viewer is heard only between cells.)

### 6) Capturing a `.cellucid-session` bundle without downloading (advanced)

This is the “no browser download” workflow, in two cells: ask, then collect.

```python
request = viewer.request_session_bundle()
```

```python
bundle = request.result()
print(bundle.list_chunk_ids())
```

See {doc}`sessions` for what a session bundle contains and how to apply it back onto an AnnData.

---

## Deep path (how notebook embedding works)

### What happens when you call `show_anndata(...)`

1. The dataset and the launch choices are checked; an in-memory AnnData is
   planned as a store in priority order, made as the viewer reads it (a MuData
   as a MuData store), its expression matrices laid out gene-major in the
   background. Nothing binds a port.
2. A {class}`~cellucid.jupyter.Widget` is shown. Its output frames the
   Cellucid app, which your browser loads from `web_source_url`
   (`https://www.cellucid.com` by default).
3. The app asks the widget for its credentials and a private channel, reads
   the dataset through it (the kernel answers each request from the dataset)
   and opens the viewer channel through it: commands go to the viewer and
   interaction events come back.
4. The widget routes those events to your hook callbacks (e.g.,
   `@viewer.on_selection`), between cells.

Architecture: {doc}`../../e_jupyter_hooks/05_architecture_viewer_channel`.

---

## API reference

### Functions

#### `show_anndata`

```{eval-rst}
.. autofunction:: show_anndata
```

### Classes (see dedicated page)

- {class}`~cellucid.jupyter.Widget` and {class}`~cellucid.AnnDataViewer` are
  documented on {doc}`viewers`.

---

## Edge cases (do not skip)

### Data size and performance
- A path is never loaded into the kernel's memory: the viewer reads an
  `.h5ad`, a Zarr store, an `.h5mu`, a `.loom` or a 10x file in place with HTTP
  range requests (only what a view needs), and an R object or a MatrixMarket
  folder once from the start.
- An in-memory AnnData is served as a Zarr store with its cells in priority
  order, each chunk made from the object when the viewer asks for it, so a
  subsample transfers only its cells; a cell-major sparse `X` is laid out
  gene-major once, in the background, into memory-mapped files removed when
  the viewer stops or Python exits. A MuData is written to a temporary
  Cellucid-ready `.h5mu` when the viewer starts.

### “My embedding is not 2D/3D”
- Cellucid supports 1D/2D/3D embeddings; if the embedding is missing or has an unexpected shape, you’ll see missing plot / errors.

### Duplicate gene identifiers
- The viewer makes repeated gene names unique (`MATR3`, `MATR3-1`), so every
  gene stays reachable; name genes by another `var` column with
  `gene_id_column` when the index holds identifiers you do not search by.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Viewer does not appear in the notebook output”

Likely causes:
- The notebook cannot show widgets: the output shows text (“This output is a
  Cellucid viewer (a notebook widget)…”) instead.
- You are not in a notebook (a plain script): `show_anndata()` raises and
  names the browser-tab alternatives.

Fix:
- Install `anywidget` where the notebook *server* runs
  (`pip install anywidget`) and reload the page, or use a notebook that shows
  widgets ({doc}`../../e_jupyter_hooks/03_supported_environments_matrix`).
- Outside a notebook use {class}`~cellucid.AnnDataViewer` or `cellucid serve`.

---

### Symptom: “The output says the page has not answered as a Cellucid viewer”

Likely causes:
- The browser cannot load the app from `web_source_url` (offline, a blocked
  network), or that address is not a Cellucid build that can be embedded.

Fix:
- Check that `https://www.cellucid.com` opens in the same browser, or host a
  copy and pass `web_source_url=`.

---

### Symptom: “The output says it is waiting for the notebook kernel”

What’s happening:
- A cell is still running; the kernel answers the viewer between cells.

Fix:
- Let the cell finish (or interrupt it); the viewer connects at once.

---

### Symptom: “`wait_for_ready()` / `request.result()` raised at once”

What’s happening:
- A notebook widget hears the viewer only between cells, so a call that waits
  for the viewer cannot be answered in the cell that waits; the message says
  what to do instead.

Fix:
- Use a hook or `viewer.state` in a later cell; for a session bundle,
  `request = viewer.request_session_bundle()` in one cell and
  `request.result()` in the next.

---

### Symptom: “I registered `@viewer.on_selection` but nothing prints”

Likely causes:
- The viewer is not ready yet (it loads once the cell that showed it ends).
- You selected nothing (or selection tools differ from expectation).
- A cell is still running: hooks run when the kernel is idle.

How to confirm:
- In a new cell, check `viewer.state.ready`.

Fix:
- Ensure the viewer is loaded, then select cells (lasso/shift-click/etc.).

---

### Symptom: “`viewer.set_color_by(...)` does nothing”

Likely causes:
- The field name does not exist in `obs` (or has a different spelling/case).

How to confirm:
- Print `adata.obs.columns` and confirm the field exists.
- Register `@viewer.on_command_error`: a refused command reports why.

Fix:
- Use the exact field name.
- Prefer a file written by {func}`~cellucid.prepare` for stable field naming across sessions.

---

### Symptom: “`viewer.highlight_cells(...)` does nothing”

Likely causes:
- Indices are out of range (wrong dataset / wrong cell order).

How to confirm:
- Try highlighting a small index like `[0]`.
- Register `@viewer.on_command_error`: a refused command reports why.

Fix:
- Ensure indices correspond to the dataset loaded in the viewer.
- If you subset/reorder AnnData, re-run the viewer from the same object.

---

### Symptom: “Selections map to the wrong cells in my AnnData”

Likely causes:
- The viewer is connected to a dataset whose row order does not match your `adata`.

How to confirm:
- Select a single cell and print identifying metadata from `adata.obs.iloc[cell_index]`.
- Compare with what you expect in the UI.

Fix:
- Ensure you apply selection indices to the *same* AnnData used to create/serve the viewer.
- Avoid reordering `adata` after starting the viewer; restart the viewer after reindexing.

---

### Symptom: “My hook callback crashes silently”

What’s happening:
- Hooks run on the viewer channel's thread, not in your cell. An exception a
  hook raises is logged with its traceback (logger `cellucid.viewer_channel`,
  "A viewer event hook failed"), never raised in the notebook. The event is
  still recorded in `viewer.state` and handed to `wait_for_event` (it is
  recorded before any hook runs), and the other hooks still run for it. The
  channel stays open.

How to confirm:
- Enable debug logging:
  ```python
  import logging
  logging.basicConfig(level=logging.DEBUG)
  ```

Fix:
- Add defensive checks inside callbacks (validate keys, lengths, types).
- Start with `@viewer.on_message` to inspect raw event payloads.

---

### Symptom: “I don’t know what’s wrong / I need a diagnostic dump”

Fix:
- Run:
  ```python
  report = viewer.debug_connection()
  report
  ```
- For a notebook widget the report includes where the app is loaded from,
  the connected outputs, whether the channel is open, the readiness state and
  recent events (it waits for nothing). For a browser-tab `AnnDataViewer` it
  also probes the server (health, info, every listed store), the web build
  cache, and pings the viewer over the channel.

---

## See also

- {doc}`server` for browser-tab workflows and SSH tunneling
- {doc}`prepare` for Cellucid-ready AnnData files (`prepare(...)`)
- {doc}`sessions` for `.cellucid-session` bundles (capture + apply)
