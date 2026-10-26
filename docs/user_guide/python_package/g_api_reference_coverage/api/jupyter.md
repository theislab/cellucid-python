# Jupyter (notebook embedding + hooks)

```{eval-rst}
.. currentmodule:: cellucid
```

This page documents the notebook-facing APIs:
- **Notebook view:** {func}`~cellucid.show_anndata` (an AnnData or MuData object, or any file or folder Cellucid reads: `.h5ad`, a Zarr store, `.h5mu`, `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds`/`.rda`, or a file written by {func}`~cellucid.prepare`)
- **Interactive control:** the viewer object ({class}`~cellucid.AnnDataViewer`) supports hooks + Python→frontend commands

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: copy/paste the “Fast path” and use the troubleshooting section when something looks wrong.
- Computational: use “Practical path” + “Edge cases” to avoid performance and data-shape pitfalls.
- Power user: use “Deep path” to understand how the server + hooks + session capture work.

**Prerequisites**
- `pip install cellucid`
- If you use AnnData: install `anndata` (and typically `numpy`, `scipy`, `pandas`).
- A notebook environment: classic Jupyter or JupyterLab. VS Code's notebook webview and Google Colab's proxy have not been tested by the project ({doc}`../../e_jupyter_hooks/03_supported_environments_matrix`).

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

### B) “I have a file written by `prepare()`”

```python
from cellucid import show_anndata

viewer = show_anndata("my_study.h5ad")  # it names itself (uns['cellucid'])
```

---

## Practical path (step-by-step + common workflows)

### 1) Notebook display vs browser tab

Cellucid always runs in a browser context (the viewer UI), but you choose **where it appears**:
- **Notebook embedding**: viewer appears inline (this page).
- **Browser tab**: use {func}`~cellucid.serve_anndata` / `cellucid serve …` (see {doc}`server`).

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

### What starts when you call `show_anndata(...)`

1. A **local HTTP server** ({class}`~cellucid.AnnDataServer`) starts on
   `127.0.0.1:8765`, or, when another viewer or program holds 8765, on the next
   free port, with a printed note (Google sign-in works only on 8765).
   `AnnDataViewer(port=...)` uses the port given or fails.
2. The exact verified viewer generation is served from that same origin.
3. The viewer UI opens its channel, one WebSocket to that server
   (`/_cellucid/ws`), authenticated with its viewer token; commands go to the
   viewer and interaction events come back on it.
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

### Classes (see dedicated page)

- {class}`~cellucid.AnnDataViewer` is documented on {doc}`viewers`.

---

## Edge cases (do not skip)

### Data size and performance
- A path is never loaded into the kernel's memory: the viewer reads an
  `.h5ad`, a Zarr store, an `.h5mu`, a `.loom` or a 10x file in place with HTTP
  range requests (only what a view needs), and an R object or a MatrixMarket
  folder once from the start.
- An in-memory AnnData is written once to a temporary `.h5ad` (sparse matrices
  column-major, uncompressed), and a MuData to a temporary Cellucid-ready
  `.h5mu`, when the viewer starts; the file is removed when the viewer stops or
  Python exits. For a large object, passing the path it was read from skips
  that write.

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
- You are not actually running inside a notebook (e.g., plain Python script).
- The notebook output is blocked (content security policy / iframe blocked).

How to confirm:
- Print `viewer.viewer_url` and try opening it in a new browser tab.

Fix:
- Use a notebook environment (classic Jupyter or JupyterLab; VS Code and Colab are untested).
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
- {func}`~cellucid.show_anndata` and `AnnDataViewer()` without `port=` never
  fail on a taken port: they take the next free one. This error means a port
  you passed to `AnnDataViewer(port=...)` is taken: pass another, or stop the
  process using it.

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
- Print `adata.obs.columns` and confirm the field exists.

Fix:
- Use the exact field name.
- Prefer a file written by {func}`~cellucid.prepare` for stable field naming across sessions.

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
- The report includes:
  - server health/info probes,
  - `anndata_store_probes`, keyed by every exact listed dataset id,
  - detected notebook context (Jupyter/Colab/VSCode),
  - viewer UI cache status,
  - recent events and readiness state.

---

## See also

- {doc}`server` for browser-tab workflows and SSH tunneling
- {doc}`prepare` for Cellucid-ready AnnData files (`prepare(...)`)
- {doc}`sessions` for `.cellucid-session` bundles (capture + apply)
