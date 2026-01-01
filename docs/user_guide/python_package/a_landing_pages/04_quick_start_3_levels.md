# Quick start (3 levels)

**Audience:** everyone (choose your depth)  
**Time:** 5–45+ minutes  
**What you’ll get:** a working Cellucid viewer (in a notebook or browser) + a path to scaling/sharing

This page gives you three depth levels. Stop as soon as you have what you need:

- **Level 1**: see an AnnData immediately (minimal setup)
- **Level 2**: create an export folder (fast, shareable, reproducible)
- **Level 3**: treat the viewer as an interactive UI in your notebook (hooks + commands)

```{note}
If you prefer a notebook-style, highly verbose walkthrough (with many edge cases and troubleshooting), start with {doc}`../f_notebooks_tutorials/index`.
```

---

## Level 1 — Fast path: `show_anndata(adata)`

This is the fastest way to prove “Cellucid works” and to start exploring a dataset.

### Prerequisites

- A notebook environment (Jupyter / JupyterLab / VSCode notebooks / Colab)
- An AnnData object **with at least one embedding** (typically `adata.obsm["X_umap"]`)

```{tip}
If you don’t want notebooks at all, skip to Level 2 and use the CLI: `cellucid serve ...`.
```

### Step 1: Load (or create) an AnnData

**Option A: load a real dataset**

```python
import anndata as ad

adata = ad.read_h5ad("my_dataset.h5ad")
```

**Option B: create a tiny toy dataset (recommended for a first-run smoke test)**

```python
import numpy as np
import pandas as pd
import anndata as ad

X = np.array([[0, 1],
              [2, 3],
              [4, 5]], dtype=np.float32)  # cells × genes

obs = pd.DataFrame(
    {"cluster": ["A", "A", "B"], "score": [0.1, 0.2, 0.3]},
    index=[f"cell_{i}" for i in range(X.shape[0])],
)
var = pd.DataFrame(index=["G1", "G2"])

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.obsm["X_umap"] = np.array([[0, 0],
                                 [1, 0],
                                 [0, 1]], dtype=np.float32)
```

### Step 2: Confirm you have an embedding

```python
list(adata.obsm.keys())
```

If you don’t see `X_umap` (or another embedding), Cellucid won’t know where to place your cells in 2D/3D space. Compute an embedding first (typically via Scanpy or your existing pipeline), then retry.

### Step 3: Show the viewer

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=600)
```

### Step 4 (optional): Confirm hooks are working

```python
@viewer.on_ready
def _ready(ev):
    print("Viewer ready:", ev)
```

### Step 5 (optional): Programmatic control

```python
viewer.set_color_by("cluster")
viewer.highlight_cells([0, 2], color="#00cc66")
```

### What this does (mental model)

`show_anndata(...)`:
- starts a small local HTTP server (AnnData → “virtual files”),
- embeds the Cellucid web app as an iframe,
- routes UI events back to Python via `POST /_cellucid/events`.

### What success looks like

- An embedded viewer appears in the notebook output cell.
- You see points (cells) in the embedding.
- Coloring by `cluster` changes point colors.

<!-- SCREENSHOT PLACEHOLDER
ID: python-quickstart-level1-success
Where it appears: Quick start → Level 1 → What success looks like
Capture:
  - UI location: notebook output cell with embedded Cellucid viewer
  - State prerequisites: dataset loaded; a categorical field selected for color-by (e.g., cluster)
  - Action to reach state: run the Level 1 code; choose “cluster” in the color-by UI
Crop:
  - Include: enough of the notebook cell to show it’s embedded + the main viewer canvas + the color-by control
  - Exclude: personal file paths, notebook file names if sensitive
Redact:
  - Remove: sample IDs/patient IDs if present
Alt text:
  - Embedded Cellucid viewer in a notebook showing cells colored by a categorical field.
Caption:
  - Level 1 success state: the viewer is embedded and responding to UI controls (color-by).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Embedded Cellucid viewer in a notebook showing cells colored by a categorical field.
:width: 100%

Level 1 success state: the viewer is embedded and responding to UI controls (color-by).
```

### Common pitfalls (Level 1)

- **Not in a notebook**: `show_anndata(...)` will print a URL instead of embedding.
- **No embedding**: you must provide `X_umap` (or equivalent) in `adata.obsm`.
- **HTTPS notebooks / remote kernels**: see {doc}`03_compatibility_matrix_must_be_explicit`.
- **First-run offline**: the viewer UI assets may not be cached yet; see {doc}`02_installation`.

---

## Level 2 — Practical path: `prepare(...); show("./export")`

This is the “real workflow” for large datasets and collaboration: you create an **export folder** that loads quickly and can be archived/shared.

### Prerequisites (minimum required inputs)

`prepare(...)` currently requires:

- `latent_space` (shape `(n_cells, n_latent_dims)`)  
  Used for outlier quantile calculation (category summaries).  
  Good choices: PCA, scVI latent space, ScanVI, etc.
- `obs` (DataFrame with `n_cells` rows)
- At least one embedding array: `X_umap_1d` and/or `X_umap_2d` and/or `X_umap_3d`

Optional but common:
- `var` + `gene_expression` (for gene coloring/search)
- `connectivities` (for neighbor-based interactions)
- `vector_fields` (velocity / drift / displacement overlays)

```{warning}
If you pass an AnnData’s matrices into `prepare(...)`, **row order is identity**.  
Every array you pass must be aligned so that row `i` always refers to the same cell across embeddings, `obs`, latent space, expression, connectivities, and vector fields.
```

### Step 1: Choose an export directory

Pick a location you can archive/share later:

```python
from pathlib import Path

EXPORT_DIR = Path("exports/pbmc_demo")
```

### Step 2: Extract inputs from AnnData

This is a “starter” extraction pattern. Your pipeline may differ.

```python
import numpy as np

X_umap_2d = np.asarray(adata.obsm["X_umap"], dtype=np.float32)

# Pick a latent space:
if "X_pca" in adata.obsm:
    latent = np.asarray(adata.obsm["X_pca"], dtype=np.float32)
else:
    # Fallback: using UMAP as latent is allowed, but not ideal.
    latent = X_umap_2d
```

### Step 3: Export (recommended defaults for most users)

```python
from cellucid import prepare

prepare(
    # Required
    latent_space=latent,
    obs=adata.obs,
    X_umap_2d=X_umap_2d,

    # Optional but common (enables gene coloring/search)
    var=adata.var,
    gene_expression=adata.X,

    # Output
    out_dir=EXPORT_DIR,
    dataset_name="PBMC demo",

    # Performance / size (recommended starting point)
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
    obs_categorical_dtype="auto",

    # Overwrite behavior
    force=True,
)
```

### Step 4: View the export (pick one)

**Option A: embed in a notebook**

```python
from cellucid import show

viewer = show(EXPORT_DIR, height=700)
```

**Option B: open in a browser via CLI (recommended for non-notebook users)**

```bash
cellucid serve exports/pbmc_demo
```

**Option C: open the folder in the hosted web app**

If you prefer not to run a local server, you can open the web app and load an export folder via the browser file picker. For click-by-click instructions, see {doc}`../../web_app/b_data_loading/index`.

### What you get on disk (high-level)

An export folder is intentionally “boring”:
- a few JSON manifests (what files exist, how to interpret them),
- plus binary files for points/fields/genes.

Common top-level files/folders:
- `dataset_identity.json`
- `obs_manifest.json` and `obs/`
- `var_manifest.json` and `var/` (only if you export gene expression)
- `points_2d.bin` / `points_3d.bin` (depending on what you provided)
- `connectivity_manifest.json` and `connectivity/` (optional)
- `vectors/` (optional)

### What success looks like

- The export folder contains the expected manifests and points files.
- The viewer loads quickly and you can:
  - color by `obs` fields,
  - (optionally) search genes and color by expression.

### Common pitfalls (Level 2)

- **Missing `latent_space`**: `prepare(...)` will error (latent is required in the current implementation).
- **Wrong shapes**: `obs` length, embedding row counts, expression matrix shape must all match.
- **Expression orientation**: `gene_expression` must be `(n_cells, n_genes)` with `len(var) == n_genes`.
- **Stale exports**: if you re-export into the same folder, use `force=True` or pick a fresh directory.

---

## Level 3 — Deep path: custom server + hooks + message handling

This level is for notebook users who want the viewer as a live UI:
- select cells in the browser → run Python analysis on that subset,
- compute results → send highlights/colors back to the viewer.

### Step 1: Create a viewer (exported data recommended)

```python
from cellucid import show

viewer = show(EXPORT_DIR, height=700)
```

You can also do this in AnnData direct mode:

```python
from cellucid import show_anndata
viewer = show_anndata(adata)
```

### Step 2: Register hooks (frontend → Python)

```python
@viewer.on_ready
def _on_ready(ev):
    print("ready:", ev)

@viewer.on_selection
def _on_selection(ev):
    cells = ev.get("cells", [])
    print("selected:", len(cells), "cells")
```

### Step 3: Send commands (Python → frontend)

```python
# highlight first 10 selected cells in green
@viewer.on_selection
def _highlight(ev):
    cells = (ev.get("cells") or [])[:10]
    viewer.highlight_cells(cells, color="#00cc66")
```

### Step 4: Make callbacks robust (recommended patterns)

- Always guard against empty selections.
- Wrap analysis code in `try/except` so one exception doesn’t “silently kill” your workflow.
- If your analysis is slow, consider:
  - running it asynchronously,
  - debouncing rapid selection events,
  - or adding an explicit “run analysis” button (future pattern).

### Step 5: Debug connectivity when things feel “haunted”

Cellucid provides a structured report:

```python
report = viewer.debug_connection()
report
```

Look for:
- `server_health` / `server_info` (server reachable?)
- `client_server_url` (is the iframe using the right URL?)
- `web_ui.cache` (do you have cached UI assets?)
- `frontend_console` (any iframe-side errors forwarded back?)

### Step 6: Clean up (avoid port leaks)

```python
viewer.stop()
```

If you created multiple viewers:

```python
from cellucid.jupyter import cleanup_all
cleanup_all()
```

### What success looks like

- Selecting cells triggers your Python callbacks.
- Python can send a highlight that visibly updates the viewer.

---

## Edge cases (read this before scaling up)

- **Huge datasets**: start with exports + quantization; AnnData-direct mode will be slower.
- **Missing embeddings**: no embedding → no viewer; compute one first.
- **NaN/Inf values**: embeddings, latent space, expression, or obs fields containing NaN/Inf can cause surprising rendering or export errors.
- **Category explosions**: categorical fields with very high cardinality can become unusable in the UI; consider grouping or filtering.
- **Remote kernels**: you may need port forwarding + `CELLUCID_CLIENT_SERVER_URL` (see {doc}`03_compatibility_matrix_must_be_explicit`).
- **Offline first run**: prefetch/caching matters (see {doc}`02_installation`).

## Troubleshooting (quick start)

Organize by symptom. Each entry points to the most likely root causes.

### Symptom: “Nothing shows up in my notebook”

**Likely causes**
- You are not in a notebook environment.
- The notebook is served from HTTPS and blocks HTTP loopback iframes.
- The viewer UI assets were not available (offline / firewall).

**How to confirm**
- If you see a printed URL instead of an embedded viewer, you’re not in a notebook.
- Run:
  ```python
  viewer.debug_connection()
  ```

**Fix**
- Use the browser workflow (`cellucid serve ...`) OR install `jupyter-server-proxy` for managed notebooks.
- Prefetch UI assets once while online.

---

### Symptom: `prepare(...)` fails with “latent_space is required…”

**Likely cause**
- You did not pass `latent_space=...`.

**Fix**
- Use PCA/scVI/etc as latent space. If you truly don’t have one yet, you can pass `X_umap_2d` as a fallback latent space (not ideal but works for export).

---

### Symptom: export succeeds, but the viewer can’t find fields / genes

**Likely causes**
- You didn’t export `gene_expression` + `var`.
- You exported only a subset of obs columns (advanced; see Data Preparation API).

**How to confirm**
- Check that `var_manifest.json` exists (for gene expression).
- Check `obs_manifest.json` for field listings.

**Fix**
- Re-export with `var=adata.var` and `gene_expression=adata.X`.

---

### Symptom: hooks don’t fire / selection doesn’t reach Python

**Likely causes**
- Viewer is not fully loaded yet (wait for `on_ready`).
- You have multiple viewers and are listening on the wrong one.
- Remote/HTTPS constraints prevent the iframe from reaching the server.

**Fix**
- Add an `on_message` logger, call `viewer.debug_connection()`, and stop old viewers (`viewer.stop()`).

## Next steps

- Detailed export inputs/edge cases: {doc}`../c_data_preparation_api/index`
- All viewing modes (CLI/server/notebook): {doc}`../d_viewing_apis/index`
- Hooks system deep dive: {doc}`../e_jupyter_hooks/index`
- Notebook-style tutorials: {doc}`../f_notebooks_tutorials/index`
- Symptom-based troubleshooting index: {doc}`../i_troubleshooting_index/index`
