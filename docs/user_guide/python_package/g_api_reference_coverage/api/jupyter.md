# Jupyter (notebook embedding + hooks)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents the notebook-facing APIs:
- **Fastest notebook preview:** {func}`~cellucid.show_anndata` (serve AnnData directly)
- **Reproducible notebook view:** {func}`~cellucid.show` (view an exported dataset directory)
- **Interactive control:** viewer objects ({class}`~cellucid.AnnDataViewer`, {class}`~cellucid.CellucidViewer`) support hooks + Python→frontend commands

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: copy/paste the “Fast path” and use the troubleshooting section when something looks wrong.
- Computational: use “Practical path” + “Edge cases” to avoid performance and data-shape pitfalls.
- Power user: use “Deep path” to understand how the server + hooks + session capture work.

**Prerequisites**
- `pip install cellucid`
- If you use AnnData: install `anndata` (and typically `numpy`, `scipy`, `pandas`).
- A notebook environment: classic Jupyter, JupyterLab, VSCode notebooks, or Google Colab.

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
- A Cellucid viewer embedded in the notebook output, showing your embedding.
- Progress while the complete configured web generation is fetched, verified,
  and published for this startup (unless output is quiet).

### B) “I already exported a dataset folder”

```python
from cellucid import show

viewer = show("./my_export")
```

---

## Practical path (step-by-step + common workflows)

### 1) Notebook display vs browser tab

Cellucid always runs in a browser context (the viewer UI), but you choose **where it appears**:
- **Notebook embedding**: viewer appears inline (this page).
- **Browser tab**: use {func}`~cellucid.serve` / {func}`~cellucid.serve_anndata` / `cellucid serve …` (see {doc}`server`).

### 2) Choosing `show_anndata` vs `show`

Use {func}`~cellucid.show_anndata` when:
- you’re iterating quickly in analysis notebooks,
- you don’t want to write an export folder yet,
- you can accept slower loads and less shareable artifacts.

Use {func}`~cellucid.show` when:
- you want reproducibility (the export folder is a concrete artifact),
- you want speed (pre-quantized, pre-packed binaries),
- you want to share the dataset directory or host it.

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

### 5) Waiting for the viewer to be ready (robust notebooks)

If you need to run code only after the UI finished loading:

```python
viewer.wait_for_ready(timeout=30)
```

### 6) Capturing a `.cellucid-session` bundle without downloading (advanced)

This is the “no browser download” workflow:

```python
bundle = viewer.get_session_bundle(timeout=60)
print(bundle.list_chunk_ids())
```

See {doc}`sessions` for what a session bundle contains and how to apply it back onto an AnnData.

---

## Deep path (how notebook embedding works)

### What starts when you call `show(...)` / `show_anndata(...)`

1. A **local HTTP server** starts (on `127.0.0.1:<port>`).
   - Exported data: {class}`~cellucid.CellucidServer`
   - AnnData mode: {class}`~cellucid.AnnDataServer`
2. The exact verified viewer generation is served from that same origin.
3. The viewer UI posts interaction events back to the local server at `/_cellucid/events`.
4. The Python viewer routes those events to your hook callbacks (e.g., `@viewer.on_selection`).

### Environment caveat: when the browser cannot reach the kernel’s localhost

Some notebook setups run the kernel remotely. In those cases, your browser may not be able to reach `http://127.0.0.1:<port>` directly.

If that happens, pass the exact browser-reachable HTTP(S) base URL as
`client_server_url=...` when constructing the viewer.

---

## API reference

### Functions

#### `show_anndata`

```{eval-rst}
.. autofunction:: show_anndata
```

#### `show`

```{eval-rst}
.. autofunction:: show
```

### Classes (see dedicated page)

- {class}`~cellucid.AnnDataViewer` and {class}`~cellucid.CellucidViewer` are documented on {doc}`viewers`.

---

## Edge cases (do not skip)

### Data size and performance
- Very large `adata.X` can be slow in AnnData mode; consider exporting with {func}`~cellucid.prepare`.
- A `.h5ad` path is always opened read-only-backed; a `.zarr` path is loaded
  eagerly.

### “My embedding is not 2D/3D”
- Cellucid supports 1D/2D/3D embeddings; if the embedding is missing or has an unexpected shape, you’ll see missing plot / errors.

### Duplicate gene identifiers
- Direct AnnData viewing rejects duplicate gene IDs and portable filename
  collisions during adapter construction; provide a unique identifier through
  `var.index` or `gene_id_column` (see {doc}`export` and {doc}`adapters`).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Viewer does not appear in the notebook output”

Likely causes:
- You are not actually running inside a notebook (e.g., plain Python script).
- The notebook output is blocked (content security policy / iframe blocked).

How to confirm:
- Print `viewer.viewer_url` and try opening it in a new browser tab.

Fix:
- Use a notebook environment (Jupyter, JupyterLab, VSCode notebooks, Colab).
- If iframes are blocked, open `viewer.viewer_url` manually.

---

### Symptom: remote notebook iframe is blank or unreachable

Likely causes:
- The browser cannot reach the kernel’s localhost server (remote kernel or strict origin policy).

Fix:
- expose the server through your notebook deployment, and
- pass that exact base URL as `client_server_url`.

---

### Symptom: “Port already in use”

Fix:
- Choose a different port (pass `port=...` to the viewer class, or stop the process using the port).

---

### Symptom: “I registered `@viewer.on_selection` but nothing prints”

Likely causes:
- The viewer is not ready yet.
- You selected nothing (or selection tools differ from expectation).

How to confirm:
- Call `viewer.wait_for_ready()` and then try selecting again.

Fix:
- Ensure the viewer is loaded, then select cells (lasso/shift-click/etc.).

---

### Symptom: “The iframe is blank / white / stuck loading”

Likely causes:
- Viewer construction raised because the configured generation could not be
  established.
- The server isn’t reachable from the browser (remote kernel / HTTPS notebook constraints).

How to confirm:
- Open `viewer.viewer_url` in a new browser tab (sometimes notebook output hides useful errors).

Fix:
- Ensure the configured source is reachable at startup and the selected
  generation directory is writable. A previously published generation is not
  substituted after source failure.
- If the browser cannot reach localhost, expose the server and pass its exact
  browser-reachable base URL as `client_server_url`.

---

### Symptom: web-generation startup failure

Likely causes:
- The configured source generation could not be fetched, verified, or
  atomically published.

How to confirm:
- Run `viewer.debug_connection()` and inspect the `web_ui` section.

Fix:
- Ensure the runtime can reach the configured source at viewer startup.
- Pass a writable `web_cache_dir` and correct the exact reported source or
  inventory failure.

---

### Symptom: “`viewer.set_color_by(...)` does nothing”

Likely causes:
- The field name does not exist in `obs` (or has a different spelling/case).
- The viewer is not ready yet.

How to confirm:
- Wait for readiness: `viewer.wait_for_ready()`.
- In AnnData mode, print `adata.obs.columns` and confirm the field exists.

Fix:
- Use the exact field name.
- Prefer exporting with {func}`~cellucid.prepare` for stable field naming across sessions.

---

### Symptom: “`viewer.highlight_cells(...)` does nothing”

Likely causes:
- The viewer is not ready yet.
- Indices are out of range (wrong dataset / wrong cell order).

How to confirm:
- `viewer.wait_for_ready()`
- Try highlighting a small index like `[0]`.

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
- Hook callbacks are wrapped; exceptions are logged rather than crashing the notebook cell.

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
- The report includes:
  - server health/info probes,
  - detected notebook context (Jupyter/Colab/VSCode),
  - viewer UI cache status,
  - recent events and readiness state.

---

## See also

- {doc}`server` for browser-tab workflows and SSH tunneling
- {doc}`export` for reproducible exports (`prepare(...)`)
- {doc}`sessions` for `.cellucid-session` bundles (capture + apply)
