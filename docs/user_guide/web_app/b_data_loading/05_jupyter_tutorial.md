# Jupyter Integration (Notebook Embedding)

This tutorial shows how to load data into Cellucid **from inside a notebook** (JupyterLab, classic Jupyter, VSCode notebooks, etc.).

You will learn:
- how to display a **pre-exported** dataset with `show()` (fastest + most reproducible)
- how to display an **AnnData / `.h5ad` / `.zarr`** directly with `show_anndata()` (most convenient for analysis)
- how to work with **vector fields** (velocity/drift overlays) in notebooks
- how the integration works under the hood (so you can debug it)
- how to drive the viewer from Python (highlight/color/visibility) and react to UI events (hooks)

If you are not in a notebook environment, start with {doc}`04_server_tutorial`.

## At A Glance

**Audience**
- Wet lab / beginner: copy/paste the “Minimal cells” sections and focus on “What success looks like”.
- Computational users: focus on backed mode, dataset sizes, vector fields, and cleanup.
- Power users: focus on remote/HPC workflows, hooks, and debugging endpoints.

**Time**
- Minimal working embed: ~5 minutes
- Full read (hooks + troubleshooting): ~20–30 minutes

**Prerequisites**
- A Jupyter environment (classic notebook, JupyterLab, or VSCode notebooks)
- `pip install cellucid`

```{important}
**Network requirement (important):** Cellucid serves the viewer UI via a **hosted-asset proxy**.

- On first run (or after the website changes), the Python server downloads `https://www.cellucid.com/index.html` + `/assets/*` and caches them on disk.
- Cache invalidation is driven by the web app’s `<meta name="cellucid-web-build-id" ...>` stamp.
- Notebook embeds then load the UI from your local Cellucid server (`http://127.0.0.1:<port>/`), avoiding HTTPS→HTTP mixed-content blocking.
- If you are offline and no cached UI is available, the iframe will show a “Cellucid viewer UI could not be loaded” page with next steps.
- In Jupyter mode, Cellucid intentionally does **not** fall back to demo data if it cannot reach the Python-side server; it shows a connectivity error so you don’t accidentally analyze the wrong dataset.

Configure the cache location with `CELLUCID_WEB_PROXY_CACHE_DIR`.
```

## Minimal Cells (Copy/Paste)

If you just want it to work, copy/paste one of these flows and then come back for details.

### Minimal: show an in-memory `AnnData`

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=600)
viewer  # (optional) display again in some notebook UIs
```

### Minimal: show a `.h5ad` / `.zarr` (recommended for large datasets)

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", height=600)  # backed mode by default
# viewer = show_anndata("data.zarr", height=600)
```

### Minimal: show a pre-exported dataset directory

```python
from cellucid import show

viewer = show("./exports/pbmc_demo", height=600)
```

```{note}
If you don’t have an export directory yet, create one with `cellucid.prepare(...)`.
See {doc}`02_local_demo_tutorial` for a complete export workflow.
```

## How It Works (Mental Model)

When you call `show(...)` or `show_anndata(...)`, Cellucid does two things:

1) Starts a **local data server** (usually `127.0.0.1:<some_port>`)
   - This server reads your data and exposes a small HTTP API (e.g. `/points_3d.bin`, `/obs_manifest.json`, `/dataset_identity.json`).
   - The server is intentionally **localhost-bound** in Jupyter mode for safety (it is not meant to be public).

2) Displays an **iframe** in your notebook pointing at the **same server**:

   ```text
   http://127.0.0.1:<port>/?jupyter=true&viewerId=<id>&viewerToken=<token>
   ```

The viewer UI and the dataset API share the same origin, so the viewer loads data from relative paths (no mixed-content).

### Why this matters

- If your notebook kernel is **local**, `127.0.0.1:<port>` is your laptop and everything “just works”.
- In **Google Colab**, the kernel runs on a remote VM; Cellucid uses Colab’s built-in HTTPS port proxy so the embed still works (your `viewer.viewer_url` won’t look like `127.0.0.1` in that case).
- If your notebook kernel is **remote (HPC / cloud)**, `127.0.0.1` in the viewer URL refers to **your laptop**, not the remote machine.
  - In that case you must use SSH port forwarding (see {ref}`remote-hpc`).

### Debugging endpoints you can open in a browser

Once `viewer` exists:

```python
print(viewer.server_url)  # e.g. http://127.0.0.1:8765
print(viewer.viewer_url)  # the embedded viewer URL (usually same as server_url + query params)
```

Then try:
- `http://127.0.0.1:<port>/_cellucid/health` (server alive?)
- `http://127.0.0.1:<port>/dataset_identity.json` (dataset id + vector fields metadata)

## Choose `show()` vs `show_anndata()`

| Function | Best for | What you pass | Performance |
|---|---|---|---|
| `show(data_dir)` | Fast, reproducible viewing | a **pre-exported** directory | Best |
| `show_anndata(data)` | Convenience in analysis workflows | `AnnData` / `.h5ad` / `.zarr` | Good (but slower than exports) |

If you’re preparing a dataset for collaborators or repeated viewing, prefer `prepare()` + `show()`.

## Option #12 — `show()` (Pre-exported Dataset)

### Step 1 — Create an export (one-time)

Use `cellucid.prepare(...)` to create an export directory.

For a complete workflow (including GitHub sharing), see {doc}`02_local_demo_tutorial`.

### Step 2 — Show it in a notebook

```python
from cellucid import show

viewer = show("./exports/pbmc_demo", height=600)
```

### What success looks like

- An interactive Cellucid viewer appears directly in the output cell.
- You can pan/zoom/select cells.
- If your export includes vector fields, the overlay can be enabled (see {ref}`vector-fields`).

### Advanced: choose a fixed port (useful for SSH tunneling)

The convenience function `show(...)` auto-picks a port. If you need a fixed port:

```python
from cellucid import CellucidViewer

viewer = CellucidViewer("./exports/pbmc_demo", port=8765, height=600)
viewer.display()
```

```{note}
Fixed ports are especially useful for remote notebooks, because you can pre-configure an SSH tunnel to `127.0.0.1:8765`.
```

```python
from cellucid import show

# Example (replace with your path):
# viewer = show("./exports/pbmc_demo", height=600)
```

## Option #13/#14 — `show_anndata()` (AnnData / `.h5ad` / `.zarr`)

`show_anndata()` is the fastest way to get started in an analysis notebook.

It supports:
- an in-memory `AnnData`
- a `.h5ad` file (opened in *backed* mode by default for lazy loading)
- a `.zarr` directory (chunked storage)

### Minimal examples

```python
from cellucid import show_anndata

viewer = show_anndata(adata)
viewer = show_anndata("data.h5ad")
viewer = show_anndata("data.zarr")
```

### Useful kwargs

`show_anndata(..., **kwargs)` forwards extra options to the AnnData adapter (and to the viewer constructor when relevant):

- `port`: fix the local server port (critical for remote/HPC; see {ref}`remote-hpc`)
- `latent_key`: choose the latent space in `obsm` (auto-detected if omitted)
- `gene_id_column`: which `var` column to use as gene IDs (default: index)
- `normalize_embeddings`: normalize coordinates to `[-1, 1]` (default: True)
- `dataset_name`: label shown in the UI
- `dataset_id`: set a stable dataset identity string (important for sessions; see {doc}`06_dataset_identity_why_it_matters`)

```{tip}
If you plan to use Cellucid sessions or share “the same dataset” across runs, set a stable `dataset_id` instead of relying on auto-generated IDs.
```

```python
from cellucid import show_anndata

# In-memory AnnData
# viewer = show_anndata(adata, height=600)

# File-backed (recommended for large datasets)
# viewer = show_anndata("/path/to/data.h5ad", height=600)
# viewer = show_anndata("/path/to/data.zarr", height=600)

# With options
# viewer = show_anndata(
#     "data.h5ad",
#     height=700,
#     latent_key="X_pca",
#     gene_id_column="gene_symbols",
#     dataset_name="PBMC demo",
# )
```

(vector-fields)=
## Vector Fields in Notebooks (Velocity/Drift Overlay)

Vector fields are optional per-cell displacement vectors (e.g. RNA velocity, drift, directed transitions) that Cellucid can render as an animated overlay on top of your embedding.

This matters for data loading because:
- the overlay is only available if the server advertises `vector_fields` in `dataset_identity.json`
- vector fields must be aligned to the **same cells** and **same embedding basis/dimension**

### Quick checklist (AnnData)

To make a vector field appear when using `show_anndata(...)`, you need:

1) A UMAP embedding in `adata.obsm` (`X_umap_2d` / `X_umap_3d` or compatible `X_umap`)
2) A vector field in `adata.obsm` with a **Cellucid-compatible key**
3) The vector array shape must be `(n_cells, dim)` where `dim` matches the embedding (2 or 3)

### Naming convention (UMAP basis)

Cellucid detects vector fields in `adata.obsm` using keys like:

- Explicit (preferred):
  - `velocity_umap_2d` (shape `(n_cells, 2)`)
  - `velocity_umap_3d` (shape `(n_cells, 3)`)
  - `T_fwd_umap_2d` (shape `(n_cells, 2)`)
- Implicit (allowed, but explicit is clash-safe):
  - `velocity_umap` (shape `(n_cells, 2)` or `(n_cells, 3)`)

For full expectations (including exported-folder layout), see {doc}`07_folder_file_format_expectations_high_level_link_to_spec`.

### Minimal example: attach a 2D vector field

```python
import numpy as np

n = adata.n_obs
adata.obsm["velocity_umap_2d"] = np.zeros((n, 2), dtype=np.float32)  # replace with real vectors

viewer = show_anndata(adata)
```

### Example: compute a drift field from a transition matrix (CellRank-style)

If you have a transition matrix `T` and UMAP coordinates, Cellucid ships helpers:

```python
from cellucid import add_transition_drift_to_obsm

# Adds e.g. "T_fwd_umap_2d" into adata.obsm (key depends on dim/basis)
out_key = add_transition_drift_to_obsm(adata, T, basis="umap", field_prefix="T_fwd")
print("Wrote:", out_key)

viewer = show_anndata(adata)
```

### Verify that the server is advertising vector fields

After `viewer = show_anndata(...)`:

```python
import json, urllib.request

with urllib.request.urlopen(viewer.server_url + "/dataset_identity.json") as f:
    ident = json.load(f)

print("Has vector_fields?", "vector_fields" in ident)
```

If the overlay is missing in the UI and `vector_fields` is absent:
- double-check your `adata.obsm` keys and shapes
- confirm you have a matching embedding dimension (2D vs 3D)
- see {doc}`08_troubleshooting_data_loading` and {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

## Programmatic Control (Python → Viewer)

Once you have a `viewer`, you can drive UI state from Python.

Common actions:

```python
# Highlight a few cells (indices are 0-based row indices into your dataset)
viewer.highlight_cells([0, 10, 42], color="#ff0000")

# Color by an obs column (must exist in the dataset)
viewer.set_color_by("cell_type")

# Hide some cells (or set visible=True to show them again)
viewer.set_visibility([0, 10, 42], visible=False)

# Reset camera
viewer.reset_view()
```

```{important}
Cell indices refer to the **row order** Cellucid is serving.

- For `show_anndata(adata)`, this is the current `adata` row order.
- If you subset/shuffle `adata` in Python, the indices will change.
```

## Reacting to the UI (Hooks: Viewer → Python)

The viewer can send events back to Python so your notebook can react to selection/hover/click.

Supported hooks:
- `@viewer.on_ready`
- `@viewer.on_selection`
- `@viewer.on_hover`
- `@viewer.on_click`
- `@viewer.on_message` (raw debugging)

### Minimal: print selections

```python
@viewer.on_selection
def handle_selection(event):
    print("Selected:", len(event["cells"]))
```

### Practical: analyze the selected cells

```python
@viewer.on_selection
def analyze(event):
    cells = event["cells"]
    subset = adata[cells].copy()
    print(subset)
    # e.g. run Scanpy plots or downstream analysis on `subset`
```

### Debug: print all messages

```python
@viewer.on_message
def debug(event):
    print(event)
```

For more about hooks, see {doc}`../../python_package/e_jupyter_hooks/index`.

## Pulling State into Python (No-Download Sessions)

Hooks are great when you want reactive code. For “pull-style” workflows, Cellucid also exposes:

### `viewer.state` (live snapshot)

`viewer.state` is a small, thread-safe snapshot of the latest events:

```python
viewer.wait_for_ready(timeout=60)
print(viewer.state.selection)  # last selection event (or None)
print(viewer.state.hover)      # last hover event (or None)
print(viewer.state.click)      # last click event (or None)
```

### Session bundle (durable saved state → AnnData)

In Jupyter, you can request the current `.cellucid-session` **as a Python object**:

```python
viewer.wait_for_ready(timeout=60)
bundle = viewer.get_session_bundle(timeout=60)

# Apply to AnnData (adds obs/var columns; stores metadata in adata.uns["cellucid"])
adata2 = bundle.apply_to_anndata(adata, inplace=False)
```

Convenience one-liner:

```python
adata2 = viewer.apply_session_to_anndata(adata, inplace=False)
```

```{important}
Session application is currently **index-based** (cell identity is the row position).

Only apply a session to an `AnnData` whose row order matches the dataset that produced the session.
```

### Debugging: `viewer.debug_connection()`

If hooks/session capture seem “stuck”, run:

```python
report = viewer.debug_connection()
report
```

This checks server endpoints (`/_cellucid/health`, `/_cellucid/info`), performs a ping/pong roundtrip, and includes recent frontend console warnings/errors forwarded to Python.
It also includes a frontend “debug snapshot” (the iframe’s `location.href`, origin, and user agent), which is useful in proxied notebook environments.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-jupyter-embedded-viewer
Suggested filename: data_loading/09_jupyter-embedded-viewer.png
Where it appears: Data Loading → 05_jupyter_tutorial → first successful display
Capture:
  - UI location: Jupyter output cell
  - State prerequisites: viewer successfully displayed
  - Action to reach state: run `viewer = show_anndata(...)` or `viewer = show(...)`
Crop:
  - Include: the notebook cell output + a bit of the notebook context (so readers know it's embedded)
  - Exclude: personal notebook filenames, dataset paths, user names
Redact:
  - Remove: any private dataset names shown in the viewer
Annotations:
  - Callouts: #1 embedded viewer, #2 notebook cell that created it
Alt text:
  - Jupyter notebook cell output containing an embedded Cellucid viewer.
Caption:
  - Explain that the viewer is interactive inside the notebook.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for an embedded Cellucid viewer inside a notebook.
:width: 100%

In Jupyter, Cellucid appears as an interactive iframe inside the notebook output.
```

## Cleanup (Do This If You Re-run Cells Often)

Each viewer starts a local server in the background.

- If you create many viewers and never stop them, you can accumulate background servers.
- If you re-run a notebook cell repeatedly, you may see port increments (8765, 8766, 8767, …).

### Recommended pattern

```python
viewer = show_anndata(...)
try:
    # ... use it ...
    pass
finally:
    viewer.stop()  # stop server + cleanup
```

```{note}
`viewer.stop()` freezes the viewer frame (best-effort) so the notebook output stays visually identical to the last interactive state, but becomes non-interactive.
```

### Stop everything created in this kernel

```python
from cellucid.jupyter import cleanup_all

cleanup_all()
```

```python
# Stop a single viewer
# viewer.stop()

# Or stop all viewers created in this session
# from cellucid.jupyter import cleanup_all
# cleanup_all()
```

(remote-hpc)=
## Remote / HPC Notebooks (SSH Tunneling Guide)

If your kernel runs on a remote machine (HPC/JupyterHub/cloud VM) but your browser is on your laptop, you must ensure the browser can reach the remote kernel’s Cellucid server.

The robust solution is SSH local port forwarding:

1) Pick a port you will use for the Cellucid data server (example: `8765`).
2) Start the viewer with that port in the notebook:

   ```python
   from cellucid import show_anndata
   viewer = show_anndata("data.h5ad", port=8765, height=600)
   print(viewer.viewer_url)
   ```

3) On your laptop, create an SSH tunnel that forwards that same local port to the remote machine:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now, when your browser loads `http://127.0.0.1:8765/?jupyter=true&...`,
it hits your laptop’s `127.0.0.1:8765`, which SSH forwards to the remote kernel’s Cellucid server.

```{important}
If you do **not** set a fixed port, Cellucid will pick the first available port (8765, 8766, ...).
That makes remote tunneling awkward because you need to update your SSH forwarding every time.
```

### Common remote variants

- **You already SSH-tunnel your Jupyter server**: you can forward *both* ports in one command (example ports shown):

  ```bash
  ssh -N \\
    -L 8888:127.0.0.1:8888 \\
    -L 8765:127.0.0.1:8765 \\
    <user>@<remote-host>
  ```

- **VSCode Remote / Remote-SSH**: use VSCode port forwarding for the Cellucid port as well (the viewer still needs it on `localhost`).
- **JupyterHub**: you typically still need a localhost-reachable port; prefer using a fixed port and ask your admin if extra forwarding is needed.

## Common Edge Cases

- **No internet access**: the first time the viewer UI loads (or after a web UI update), the server must fetch `https://www.cellucid.com/index.html` + `/assets/*`. If the cache is empty, the iframe shows a “what to do” page; run once while online (or use a persistent `CELLUCID_WEB_PROXY_CACHE_DIR`).
- **Notebook blocks iframes** (security policy): you may need to open `viewer.viewer_url` in a new browser tab.
- **Port exhaustion**: if many ports are in use, Cellucid may fail to find a free one.
- **Corporate proxies / ad blockers**: can block cross-origin requests or event POSTs (hooks).
- **Huge in-memory `AnnData`**: can exhaust kernel RAM; prefer `.h5ad` backed mode or `.zarr`.

(troubleshooting)=
## Troubleshooting (Massive)

This section is intentionally redundant and explicit: it is designed for “I need to fix this now”.

### Symptom: “The viewer doesn’t appear (blank output cell)”

**Likely causes (ordered)**
1) You are not actually running in a Jupyter environment (e.g. plain Python script).
2) The notebook blocks iframes (security policy).
3) The viewer UI assets could not be loaded (hosted-asset proxy blocked and no cached copy available).

**How to confirm**
- Print the viewer URL:

  ```python
  print(viewer.viewer_url)
  ```

- Open that URL in a normal browser tab.

**Fix**
- Ensure you are using Jupyter/JupyterLab/VSCode notebooks.
- If you have outbound HTTPS, confirm the kernel/runtime can reach `https://www.cellucid.com` (used by the hosted-asset proxy).
- If you are offline, run once while online to populate the hosted-asset proxy cache (and consider setting `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent path).
- If iframes are blocked, open the URL manually in a new tab.

### Symptom: “The viewer loads, but it says it cannot connect / everything is empty”

**Likely causes**
- The local data server is not reachable from your browser.
  - common in remote/HPC notebooks without tunneling
  - also happens if you used a port that is blocked locally

**How to confirm**
- Open `viewer.server_url + "/_cellucid/health"` in a browser.
- Check that you get JSON back (status ok).

**Fix**
- Local notebook: restart the kernel and re-run, then try again.
- Remote/HPC: follow {ref}`remote-hpc` (set a fixed port + SSH forward it).

### Symptom: “Port already in use / it keeps picking new ports”

**Likely causes**
- Old viewers still running from earlier cells.
- Some other process is using the default port range.

**Fix**
- Call `viewer.stop()` when done.
- Or run `cleanup_all()`.
- Restart the kernel if needed.
- If you need a stable port, set `port=...` on `show_anndata(...)` (or use `CellucidViewer(..., port=...)`).

### Symptom: “`show_anndata` says no UMAP embeddings”

**Likely causes**
- Your AnnData lacks `obsm['X_umap_2d']` / `X_umap_3d` (or compatible `X_umap`).

**How to confirm**
- `print(adata.obsm.keys())`

**Fix**
- Compute UMAP and store it under one of the supported keys.

### Symptom: “Vector field overlay is missing (no toggle / no fields)”

**Likely causes**
- Your vectors are not stored under a detected key (e.g. missing `_umap` suffix).
- Shape mismatch: vectors are `(n_cells, 2)` but you are viewing in 3D (or vice versa).
- You have vectors, but they don’t match `adata.n_obs` (filtered/shuffled mismatch).

**How to confirm**
- Inspect `dataset_identity.json`:

  ```python
  import json, urllib.request
  with urllib.request.urlopen(viewer.server_url + "/dataset_identity.json") as f:
      ident = json.load(f)
  print(ident.get("vector_fields"))
  ```

**Fix**
- Follow the naming rules in {ref}`vector-fields` (and {doc}`07_folder_file_format_expectations_high_level_link_to_spec`).
- If you have both 2D and 3D points, make sure you provide the matching vector dimension.

### Symptom: “Gene search returns nothing / wrong gene IDs”

**Likely causes**
- Your `var_names` are Ensembl IDs but you’re searching symbols (or vice versa).
- Your gene IDs are in a different `var` column.

**Fix**
- Pass `gene_id_column="..."` to `show_anndata()`.
- Or export with `prepare(var_gene_id_column="...")` and use `show()`.

### Symptom: “Hooks don’t fire (selection events never reach Python)”

**Likely causes**
- Requests from the viewer to `/_cellucid/events` are blocked by a proxy/ad blocker.
- The local server is unreachable from the browser (hooks need HTTP POST).

**How to confirm**
- Register a raw handler:

  ```python
  @viewer.on_message
  def debug(event):
      print(event)
  ```

- Open browser devtools → Network and look for requests to `/_cellucid/events`.

**Fix**
- Ensure the data server is reachable (`/_cellucid/health` works).
- Temporarily disable extensions/ad blockers for `cellucid.com`.

### Symptom: “It’s slow compared to pre-exported data”

**Explanation**
- Direct AnnData mode is designed for convenience, not maximum speed.

**Fix**
- Export once with `prepare()` and use `show()`.
- For large `.h5ad`, prefer `.zarr` when feasible.

## Next Steps

- Stable dataset identity + sessions: {doc}`06_dataset_identity_why_it_matters`
- File/key expectations (including vector fields): {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- General loading troubleshooting: {doc}`08_troubleshooting_data_loading`
- Vector field overlay usage: {doc}`../i_vector_field_velocity/index`
- Hooks deep dive: {doc}`../../python_package/e_jupyter_hooks/index`
- If you want to share datasets publicly: export + GitHub workflow ({doc}`02_local_demo_tutorial`)
