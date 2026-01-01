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

viewer = show_anndata(adata)  # or "data.h5ad" / "data.zarr"
```

What you should see:
- A Cellucid viewer embedded in the notebook output, showing your embedding.
- (Optionally) a one-time “viewer UI cache” download progress message on first run.

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

viewer = show_anndata(adata)

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

## Screenshot placeholders (optional but recommended)

<!-- SCREENSHOT PLACEHOLDER
ID: jupyter-embedded-viewer-success
Where it appears: Jupyter (notebook embedding + hooks) → Fast path
Capture:
  - UI location: a notebook cell output showing the embedded Cellucid iframe
  - State prerequisites: dataset loaded successfully (points visible)
  - Action to reach state: run `viewer = show_anndata(adata)` and wait until the plot renders
Crop:
  - Include: the notebook output area + enough of the viewer UI to orient the reader (left sidebar + plot)
  - Exclude: browser tabs, personal file paths, kernel names that reveal identity
Redact:
  - Remove: dataset/sample IDs if private, usernames/emails, local filesystem paths
Annotations:
  - Callouts:
    - (1) viewer panel / sidebar
    - (2) embedding plot area
    - (3) a visible field selector or legend (if present)
  - Style: 2–3 callouts max, consistent color
Alt text:
  - Cellucid viewer embedded in a notebook cell, showing an embedding plot and a sidebar.
Caption:
  - The Cellucid viewer running inline in a notebook after calling `show_anndata(adata)`; the sidebar controls and the embedding plot confirm the viewer is interactive.
-->
```{figure} ../../../../_static/screenshots/placeholder-screenshot.svg
:alt: Cellucid viewer embedded in a notebook cell, showing an embedding plot and a sidebar.
:width: 100%

The Cellucid viewer running inline in a notebook after calling `show_anndata(adata)`; the sidebar controls and the embedding plot confirm the viewer is interactive.
```

---

## Deep path (how notebook embedding works)

### What starts when you call `show(...)` / `show_anndata(...)`

1. A **local HTTP server** starts (on `127.0.0.1:<port>`).
   - Exported data: {class}`~cellucid.CellucidServer`
   - AnnData mode: {class}`~cellucid.AnnDataServer`
2. The viewer UI is loaded from the same origin via a **hosted-asset proxy** (to avoid mixed-content and cross-origin issues).
3. The viewer UI posts interaction events back to the local server at `/_cellucid/events`.
4. The Python viewer routes those events to your hook callbacks (e.g., `@viewer.on_selection`).

### Environment caveat: when the browser cannot reach the kernel’s localhost

Some notebook setups run the kernel remotely. In those cases, your browser may not be able to reach `http://127.0.0.1:<port>` directly.

If that happens, set:
- `CELLUCID_CLIENT_SERVER_URL` to a browser-reachable URL for the running Cellucid server (HTTPS if required).

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
- If your `.h5ad` is huge, prefer serving in **backed/lazy** mode (default) rather than loading fully into RAM.

### “My embedding is not 2D/3D”
- Cellucid supports 1D/2D/3D embeddings; if the embedding is missing or has an unexpected shape, you’ll see missing plot / errors.

### Duplicate gene identifiers
- If gene IDs are duplicated, gene lookup behavior can be ambiguous; ensure you have a stable gene ID strategy (see {doc}`export` and {doc}`adapters`).

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
- If iframes are blocked, open the printed URL manually.

---

### Symptom: “Notebook proxy required” / the iframe shows a proxy warning

Likely causes:
- The browser cannot reach the kernel’s localhost server (remote kernel or strict origin policy).

Fix options:
- Recommended: install/enable `jupyter-server-proxy` (environment dependent).
- Or set `CELLUCID_CLIENT_SERVER_URL` to a browser-reachable URL for the server.

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
- The viewer UI failed to load (offline, blocked network, or cache missing).
- The server isn’t reachable from the browser (remote kernel / HTTPS notebook constraints).

How to confirm:
- Open `viewer.viewer_url` in a new browser tab (sometimes notebook output hides useful errors).

Fix:
- If you are offline, run once while online (to populate the viewer UI cache), then retry.
- If the browser cannot reach localhost, set `CELLUCID_CLIENT_SERVER_URL` (or use a notebook proxy).

---

### Symptom: “Cellucid: Viewer UI unavailable” page

Likely causes:
- The hosted viewer assets could not be fetched and the local cache is empty/stale.

How to confirm:
- Run `viewer.debug_connection()` and inspect the `web_ui` section.

Fix:
- Ensure the runtime can reach the hosted viewer assets once (HTTPS).
- Set `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent, writable directory.

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
