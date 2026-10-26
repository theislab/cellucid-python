# Jupyter Integration (Notebook Embedding)

This tutorial shows how to load data into Cellucid **from inside a notebook** (JupyterLab, classic Jupyter, VSCode notebooks, etc.).

You will learn:
- how to display an **AnnData or MuData object, an `.h5ad`, a Zarr store, or any
  other file or folder Cellucid reads** with `show_anndata()`
- how to work with **vector fields** (velocity/drift overlays) in notebooks
- how the integration works under the hood (so you can debug it)
- how to drive the viewer from Python (highlight/color/visibility) and react to UI events (hooks)

If you are not in a notebook environment, start with {doc}`04_server_tutorial`.

:::{tip}
Prefer to begin from a real, runnable notebook? The
{doc}`Pancreas Jupyter walkthrough
<../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery>`
links the checked notebook, its exact launch command, and real captures of the
embedded viewer and connection report.
:::

## At A Glance

**Audience**
- Wet lab / beginner: copy/paste the “Minimal cells” sections and focus on “What success looks like”.
- Computational users: focus on how the viewer reads a store in place, dataset
  sizes, vector fields, and cleanup.
- Power users: focus on remote/HPC workflows, hooks, and debugging endpoints.

**Time**
- Minimal working embed: ~5 minutes
- Full read (hooks + troubleshooting): ~20–30 minutes

**Prerequisites**
- A Jupyter environment (classic notebook, JupyterLab, or VSCode notebooks)
- `pip install cellucid`

```{important}
**Network requirement (important):** Cellucid serves the viewer UI from one
verified local generation.

- Each viewer/server startup establishes the complete source generation declared
  by `cellucid-web-assets.json`, including `index.html`, root browser metadata,
  and `/assets/*`.
- In notebooks that step is silent, like the rest of the build: a progress bar
  appears only when you establish the generation yourself with
  `viewer.ensure_web_ui_cached()`.
- Generation identity and every asset byte are verified against the inventory.
- Notebook embeds load the UI from the exact `client_server_url=` supplied by
  the caller, or from direct loopback when it is omitted.
- A source or validation failure stops startup and is reported directly; no
  stale UI generation or unrelated dataset is substituted.

Configure the generation location with `web_cache_dir=...`.
Clear the selected generation with `cellucid.clear_web_cache()` or
`viewer.clear_web_cache()`. Manually establish the source generation with
`viewer.ensure_web_ui_cached()` (usually not needed). To verify an existing
generation without network access or mutation, call
`viewer.ensure_web_ui_cached(force=False)`.
```

## Minimal Cells (Copy/Paste)

If you just want it to work, copy/paste one of these flows and then come back for details.

### Minimal: show an in-memory `AnnData`

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer  # (optional) display again in some notebook UIs
```

### Minimal: show a `.h5ad` file or a Zarr store, read in place

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A Zarr store folder or a zipped Zarr store (.zip) takes the same arguments.
```

### Minimal: show a Cellucid-ready AnnData

A file written by `cellucid.prepare(...)` names itself in `uns['cellucid']`,
so it needs no `dataset_name` or `dataset_id`:

```python
from cellucid import show_anndata

viewer = show_anndata("pbmc_demo.zarr", height=600)
```

### What success looks like

```{figure} ../../../_static/screenshots/jupyter/show-anndata-driven-from-python.png
:alt: The viewer cellucid.show_anndata() serves, opened at its viewer URL, after Python called viewer.set_color_by("clusters") and viewer.highlight_cells() on the Ductal cells: the Pancreas sample coloured by clusters with the Ductal cells highlighted.
:width: 1440px

The viewer `show_anndata()` serves, opened at its `viewer.viewer_url` in a
browser tab: in a notebook the same page sits *inside* the output cell, and
scrolls with the notebook. Here Python has already coloured it by `clusters`
and highlighted the Ductal cells; the Dataset box reads *Zarr v2 folder ·
Python*.
```

Three things confirm the call worked:

1. an interactive viewer appears in the output cell (not a blank area);
2. `viewer.wait_for_ready(...)` returns rather than timing out;
3. the cell and dimension counts it reports match the data you passed.

## How It Works (Mental Model)

When you call `show_anndata(...)`, Cellucid does two things:

1) Starts a **local data server** (usually `127.0.0.1:<some_port>`)
   - It computes nothing: it serves the store's bytes under
     `/_cellucid/anndata/<dataset_id>.h5ad` (or `.zarr.zip`, or `.zarr/<key>`
     for a Zarr store folder, or the file's own ending for another format),
     and the viewer's own HDF5/Zarr reader reads them in place with HTTP range
     requests, only what a view needs. An in-memory `AnnData` (`MuData`) is
     first written once to a temporary `.h5ad` (`.h5mu`), removed when the
     viewer stops.
   - The server is intentionally **localhost-bound** in Jupyter mode for safety (it is not meant to be public).

2) Displays an **iframe** in your notebook pointing at the **same server**:

   - Local notebooks often use direct loopback:
     ```text
     http://127.0.0.1:<port>/?viewer=true&viewerId=<id>&viewerToken=<token>
     ```
   - For HTTPS/remote notebooks, the caller may expose the selected port
     through an HTTPS proxy and pass that exact browser base:
     ```text
     https://<notebook-origin>/<base>/proxy/<port>/?viewer=true&viewerId=<id>&viewerToken=<token>
     ```

The viewer UI and the dataset API share the same origin, so the viewer loads data from relative paths (no mixed-content).

### Why this matters

- If your notebook kernel is **local**, `127.0.0.1:<port>` is your laptop and everything “just works”.
- In **Google Colab**, the kernel runs on a remote VM. Obtain Colab's HTTPS
  proxy base for a fixed port, pass it as `client_server_url=`, and declare its
  host name as `allowed_hosts=[...]`.
- For a **remote/HTTPS Jupyter server** (common on JupyterHub), configure a
  browser-reachable route, pass its exact base as `client_server_url=`, and
  declare the route's host name as `allowed_hosts=[...]`.
- `allowed_hosts=` is required behind any proxy: the proxy forwards the
  browser's `Host` header unchanged, and Cellucid refuses an authority that does
  not name it with `421 Misdirected Request`. That refusal is what stops DNS
  rebinding, and a proxied request carries no signal that distinguishes it from
  a rebound one, so the proxy has to be declared.
- Cellucid does not discover Jupyter Server Proxy or call Colab's proxy API.
- If your kernel is remote but your notebook frontend cannot use a server proxy (e.g. unusual VSCode/webview setups), you may still need SSH port forwarding (see {ref}`remote-hpc`).

### Debugging endpoints you can open in a browser

Once `viewer` exists:

```python
print(viewer.server_url)  # e.g. http://127.0.0.1:8765
print(viewer.viewer_url)  # the embedded viewer URL (usually same as server_url + query params)
```

Then try:
- `http://127.0.0.1:<port>/_cellucid/health` (server alive?)
- `http://127.0.0.1:<port>/_cellucid/datasets` (the dataset id and name,
  where the store is, and the launch choices you set)

## `show_anndata()` (AnnData / `.h5ad` / Zarr store)

`show_anndata()` is the fastest way to get started in an analysis notebook.

It supports:
- an in-memory `AnnData`, written once to a temporary `.h5ad` (sparse
  matrices column-major, uncompressed) that is removed when the viewer stops,
  or an in-memory `MuData`, written once to a temporary `.h5mu`
- a `.h5ad` file, served as it is on disk
- a Zarr store folder (Zarr v2 with `.zgroup`, or v3 with `zarr.json`) or a
  zipped Zarr store (`.zip`, such as `.zarr.zip`), served as it is on disk
- any other file or folder Cellucid reads (`.h5mu` or MuData Zarr, `.loom`, a
  10x `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds` / `.rda`),
  served as it is on disk ({doc}`07_file_format_expectations`)

Nothing is copied or converted: the viewer reads the store itself, with range
requests. The one exception is a store holding an array stored with a codec
the browser cannot decode (an HDF5 filter such as Blosc, Zstandard or szip; a
Zarr codec such as bz2, lzma or zfpy): it is read once with `anndata` and its
codec plugins and served as a decodable temporary `.h5ad`.

### Minimal examples

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer.stop()

viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### Exact parameters

`show_anndata(...)` has a closed signature. It accepts `data` and `height`,
followed by these exact keyword-only parameters:

- `client_server_url`: exact browser-reachable base URL for remote notebooks
- `allowed_hosts`: extra exact `Host` names the data server answers to; required
  when `client_server_url` points at a reverse proxy. One bare host name per
  entry — no port, no scheme, no wildcard
- `web_source_url`: exact origin publishing the web asset inventory
- `web_cache_dir`: directory where the verified generation is published
- `latent_key`: the `obsm` latent space (scVI, PCA) the outlier filter
  measures in; `None` takes the store's own default, else no latent space (the
  filter then waits for one chosen in **Data mapping**)
- `gene_id_column`: `None` uses `var.index` (`raw.var.index` with
  `expression_source="raw"`); any non-blank string names that exact column
  (default: `None`). Empty values take the index name and repeated names are
  made unique (`MATR3`, `MATR3-1`)
- `obs_keys`: exact `obs` columns the viewer shows, in this order; `None`
  (the default) shows every column
- `normalize_embeddings`: whether the viewer scales embeddings into its
  coordinate range (default `None`: the viewer's default, which scales them)
- `centroid_min_points`: the fewest cells a category needs to get outlier
  quantiles (default `None`: 10)
- `embedding_keys`: the `obsm` array a dimension shows by default, such as
  `{"2d": "X_tsne"}` or `{"2d": ("X_pca", [0, 1])}` with 0-based columns
  (default: `None`, which uses `X_umap_1d`/`X_umap_2d`/`X_umap_3d` or a bare
  `X_umap`, then a deterministic best guess)
- `knn_key`: the `obsp` matrix the viewer draws as the neighbor graph, such
  as `"connectivities"` (default: `None`: `obsp['connectivities']` when
  present). The graph is read only when a view draws it
- `expression_source`: where expression is read: `"X"`, `"raw"`, or
  `"layer:<name>"` (default `None`: `X`, or without `X` `raw`, then a layer,
  reported as guessed)
- `dataset_name`: non-empty, unpadded label without control characters;
  human-readable Unicode is preserved. Defaults to the store's
  `uns['cellucid']['name']`; required when it has none
- `dataset_id`: portable 1–180 character ASCII identity (letter or digit
  first; then letters, digits, `.`, `_`, or `-`; no trailing `.` or reserved
  Windows device name), important for sessions; see
  {doc}`06_dataset_identity_why_it_matters`. Defaults to the store's
  `uns['cellucid']['id']`; required when it has none
- `vector_field_default`: the vector field the overlay opens on: the
  `<field>` of an `obsm['<field>_<basis>']` array drawn over `obsm['X_<basis>']`
  (`"velocity"` for `velocity_umap`)
- `modality`: the modality a MuData dataset opens on (default `None`: its
  default modality); see {doc}`19_multimodal_mudata`

Every option you name is checked against the store's metadata before the
server starts — a key that does not exist raises at once — and passed to the
viewer as a launch default. An option left at `None` is taken from the store's
own `uns['cellucid']['defaults']`, else from the viewer's default rule (see
{doc}`13_data_mapping`).

The convenience function does not accept `port`. Use
`AnnDataViewer(..., port=<fixed port>)` when a tunnel or proxy requires a fixed
port.

```{tip}
Use the same exact `dataset_id` whenever you reopen the same dataset generation.
```

```python
from cellucid import show_anndata

# In-memory AnnData
# viewer = show_anndata(
#     adata,
#     height=600,
#     dataset_name="My study",
#     dataset_id="my-study-v1",
# )

# An .h5ad file, read in place
# viewer = show_anndata(
#     "/path/to/data.h5ad",
#     height=600,
#     dataset_name="My study",
#     dataset_id="my-study-v1",
# )

# A Zarr store folder (or a zipped Zarr store), read in place
# viewer = show_anndata(
#     "/path/to/data.zarr",
#     height=600,
#     dataset_name="My study",
#     dataset_id="my-study-v1",
# )

# With options
# viewer = show_anndata(
#     "data.h5ad",
#     height=700,
#     latent_key="X_pca",
#     gene_id_column="gene_symbols",
#     dataset_name="PBMC demo",
#     dataset_id="pbmc-demo",
# )
```

### Which `obsm` key becomes the embedding

`X_umap_1d` / `X_umap_2d` / `X_umap_3d` each decide their own dimension. An
object declaring none of them is shown from a bare `X_umap` at the dimension its
column count states, so a Scanpy object needs no rename. Full rules, including
what an unreadable width reports:
{doc}`../../python_package/d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.

### The neighbor graph is read when it is drawn

The viewer shows `obsp['connectivities']` as the neighbor graph when the store
has one, and reads it only when a view draws connectivity edges or the KNN drag
tool is used. To open on another graph, name it:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

Naming a matrix the store does not have raises before the server starts:

```text
knn_key names obsp['connectivities'], which this store does not have. obsp holds: []. Compute the neighbour graph with sc.pp.neighbors(adata), or leave knn_key unset.
```

A graph is checked when it is read: it must be symmetric, so the `distances`
matrix `sc.pp.neighbors()` also writes is refused. Any other `obsp` graph can be
chosen later under **Session › Data mapping › Neighbour graph**
({doc}`13_data_mapping`).

### What the constructor does

A notebook viewer suppresses server output, so the four-step report a terminal
prints — `[1/4] Detecting format` through `[4/4] Starting server` (see
{doc}`04_server_tutorial`) — never appears in a cell. The constructor does the
same work, in this order:

1. **Format** — the path is classified by the input format registry (an
   `.h5ad` file, a Zarr store folder, a zipped Zarr store, an `.h5mu`, …); an
   in-memory `AnnData` is written once to a temporary `.h5ad` (a `MuData` to an
   `.h5mu`), as is a store with an array the browser cannot decode.
2. **Metadata** — the cell and gene counts and the names of `obs`, `var`,
   `obsm`, `obsp` and the layers are read; no array is read.
3. **Launch choices** — every option you named is checked against those names.
4. **Viewer generation and bind** — the configured web generation is
   established and the port is bound.

Everything else happens in the browser, as each view needs it.

(vector-fields)=
## Vector Fields in Notebooks (Velocity/Drift Overlay)

Vector fields are optional per-cell displacement vectors (e.g. RNA velocity, drift, directed transitions) that Cellucid can render as an animated overlay on top of your embedding.

This matters for data loading because:
- vector fields must be aligned to the **same cells** and **same embedding basis/dimension**
- for `show_anndata()`, every `obsm` array named `<field>_<basis>` after an
  embedding `X_<basis>` with the same number of columns is offered as the
  field `<field>`, and is read only when the overlay shows it

### Quick checklist (AnnData)

To make a vector field appear when using `show_anndata(...)`, you need:

1) An embedding in `adata.obsm` named `X_<basis>` (`X_umap`, `X_umap_3d`,
   `X_pca`, `X_branches`, ...) that a view shows
2) A vector field in `adata.obsm` named `<field>_<basis>` after it (scVelo's
   convention: `velocity_umap` for `X_umap`)
3) The vector array shape must be `(n_cells, n_columns)`, with exactly as many
   columns as the embedding
4) If the store declares more than one field, the overlay opens on the one
   you name with `vector_field_default="velocity"` (or, for the
   forward-drift example, `vector_field_default="T_fwd"`)

### Naming convention (scVelo's `<field>_<basis>`)

Cellucid detects vector fields in `adata.obsm` by the key `<field>_<basis>`
beside the embedding `X_<basis>`, for any basis:

- `velocity_umap` (shape `(n_cells, 2)`, over a 2-column `X_umap`)
- `velocity_umap_3d` (shape `(n_cells, 3)`, over `X_umap_3d`)
- `T_fwd_umap` (shape `(n_cells, 2)`, over `X_umap`)
- `velocity_branches` (over `X_branches`, as many columns as it)

The keys of one field on every embedding are one field: `velocity_umap` and
`velocity_umap_3d` are both listed as `velocity`. A `velocity_umap_2d` beside a
bare `X_umap` (no `X_umap_2d`) is not detected; name it `velocity_umap`.

For full expectations (files, catalogs and `uns['cellucid']`), see {doc}`07_file_format_expectations`.

### Minimal example: attach a 2D vector field

```python
import numpy as np

# A complete synthetic counter-clockwise field for an interface check.
# Do not interpret this constructed field as biological velocity.
coordinates = np.asarray(adata.obsm["X_umap"], dtype=np.float32)  # 2 columns
offsets = coordinates - coordinates.mean(axis=0, keepdims=True)
adata.obsm["velocity_umap"] = np.column_stack(
    (-offsets[:, 1], offsets[:, 0])
).astype(np.float32)

viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### Example: compute a drift field from a transition matrix (CellRank-style)

If you have a transition matrix `T` and UMAP coordinates, Cellucid ships helpers:

```python
from cellucid import add_transition_drift_to_obsm

# Adds "T_fwd_umap" (<field_prefix>_<basis>) into adata.obsm, as wide as X_umap
out_key = add_transition_drift_to_obsm(
    adata,
    T,
    basis="umap",
    field_prefix="T_fwd",
    normalize_rows=False,
)
print("Wrote:", out_key)

viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### Verify that the viewer bound the vector fields

After `viewer = show_anndata(...)`, the data mapping the viewer reports lists
every bound field:

```python
viewer.wait_for_event("data_mapping")
print(viewer.data_mapping["vectorFields"])
```

If the overlay is missing in the UI and the list is empty:
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

- For a viewer constructed directly from `adata`, this is the current
  `adata` row order.
- If you subset/shuffle `adata` in Python, the indices will change.
```

### Read and change the data mapping

The *data mapping* says which stored arrays the viewer plots and reads: the
embedding of each dimension slot, the bound vector fields, and for an AnnData
source the neighbour graph, the gene names, the expression matrix and the
latent space. The viewer reports it after the dataset loads, and
`viewer.data_mapping` holds the latest report:

```python
viewer.wait_for_event("data_mapping")
viewer.data_mapping["embeddings"]
# {'1d': None, '2d': {'key': 'X_umap', 'columns': [0, 1]}, '3d': None}
viewer.data_mapping_origin   # "default", "best_guess", "session", "user" or "notebook"
```

`viewer.set_data_mapping(...)` replaces the top-level fields you name, checks
the result before sending it, and waits until the viewer has applied it:

```python
viewer.set_data_mapping(
    embeddings={
        "enabled": [{"key": "X_pca", "columns": [0, 1]}],
        "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [0, 1]}, "3d": None},
    },
)
```

A field this kind of mapping does not have raises `TypeError`, an invalid
value raises `ValueError`, and a mapping the viewer refuses — an array that
fails its load-time check, for example — raises `RuntimeError` with the
viewer's reason. The viewer reads `GET /_cellucid/protocol` when it connects
and reports mappings (and refused commands) only to a `cellucid` that lists
those events there, as this release does. A viewer older than this release
never reports a mapping, so `set_data_mapping()` raises because there is
nothing to change yet.

## Reacting to the UI (Hooks: Viewer → Python)

The viewer can send events back to Python so your notebook can react to selection/hover/click, and to a command the viewer refused.

Supported hooks:
- `@viewer.on_ready`
- `@viewer.on_selection`
- `@viewer.on_hover`
- `@viewer.on_click`
- `@viewer.on_command_error` (a command this notebook sent was rejected by the
  viewer)
- `@viewer.on_data_mapping` (the viewer resolved or changed which arrays it plots)
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
adata2 = bundle.apply_to_anndata(
    adata,
    expected_dataset_id="my-study-v1",
    inplace=False,
)
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

This checks server endpoints (`/_cellucid/health`, `/_cellucid/info`,
`/_cellucid/datasets`), probes the served AnnData store with a `HEAD` and a
one-byte range request (`anndata_store_probes`), performs a
ping/pong roundtrip, and reports recent accepted-event counts.
It also includes a frontend “debug snapshot” (the iframe’s `location.href`, origin, and user agent), which is useful in proxied notebook environments.

An abridged healthy report (the real one also carries the viewer and server
URLs with the viewer's token, the web cache's absolute path, the browser's
user agent and timestamps; `<version>` is your package version):

```text
{'server_running': True,
 'channel_connected': True,
 'server_health': {'status': 'ok', 'type': 'anndata', 'version': '<version>',
                   'format': 'zarr', 'n_cells': 3696, 'n_genes': 3753},
 'anndata_store_probes': {'pancreas': {'head': {'status': 200, ...},
                                       'range': {'status': 206,
                                                 'content_range': 'bytes 0-0/19089'}}},
 'recent_events': {'count': 2, 'by_type': {'data_mapping': 1, 'ready': 1},
                   'command_errors': []},
 'frontend_roundtrip': {'ok': True, 'pong': {...}},
 'frontend_debug_snapshot': {'ok': True, 'snapshot': {'connected': True, ...}},
 ...}
```

`viewer.debug_connection()` returns a dictionary; it prints nothing by itself.
Read it top to bottom: `server_health` proves the server is alive,
`anndata_store_probes` proves it is serving the dataset you think it is, and
`frontend_roundtrip` proves the browser can talk back to Python. A hook problem
that survives all three is a proxy or extension problem, not a Cellucid one.

:::{note}
`debug_connection()` is a diagnostic, and it degrades rather than raising. Any
probe that fails is recorded as a `*_error` string beside the others instead of
stopping the report, so read the keys rather than assuming an exception would
have told you.
:::

## Cleanup (Do This If You Re-run Cells Often)

Each viewer starts a local server in the background.

- If you create many viewers and never stop them, you can accumulate background servers.
- A viewer without `port=` takes 8765, or the next free port when another
  viewer holds it (a printed note says so); repeated cells can therefore leave
  several servers on 8765, 8766, …. Google sign-in (Google Drive data,
  community annotation) works only on port 8765.

### Recommended pattern

```python
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
try:
    viewer.wait_for_ready(timeout=60)
finally:
    viewer.stop()  # stop server + cleanup
```

```{note}
`viewer.stop()` requests an exact freeze, stops the server, unregisters hooks,
and reports any freeze or shutdown failure. After a successful stop, the last
frame remains visible but is non-interactive.
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
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       "data.h5ad",
       port=8765,
       height=600,
       dataset_name="My study",
       dataset_id="my-study-v1",
   )
   print(viewer.viewer_url)
   ```

3) On your laptop, create an SSH tunnel that forwards that same local port to the remote machine:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now, when your browser loads `http://127.0.0.1:8765/?viewer=true&...`,
it hits your laptop’s `127.0.0.1:8765`, which SSH forwards to the remote kernel’s Cellucid server.

```{important}
If you do **not** set a fixed port on `AnnDataViewer`, Cellucid takes 8765, or
the next free port when another viewer holds it. A tunnel forwards one fixed
port, so name it with `port=`: a viewer given a port uses that port or fails.
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

### Kernel on an HPC compute node (two-hop forward)

The tunnel above assumes your `ssh` reaches the machine running the kernel. On
many clusters it does not: only the login node accepts SSH, and your kernel runs
on a compute node inside a job. A kernel there has exactly the same two-hop
problem as `cellucid serve`, and the same shape of solution — the forward is
opened *by the login node*, which then connects onward to the compute node:

```bash
ssh -N -L 8765:compute-node:8765 you@login-node
```

The forward target `compute-node:8765` is resolved by the login node, so the
second connection arrives at the compute node from another machine. A server
bound to `127.0.0.1` refuses it: loopback accepts its own machine only. That is
what the wildcard bind is for, and `AnnDataViewer` does not offer one — it
always binds loopback, for safety in notebooks. So on a compute node, start a
server from the kernel instead of an embedded viewer, and open it in its own
browser tab:

```python
from cellucid import serve_anndata

serve_anndata(
    "data.h5ad",
    port=8765,
    host="0.0.0.0",
    open_browser=False,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

The address you open stays `http://127.0.0.1:8765/?anndata=true`: the tunnel's
near end is on your own machine, whatever the server bound to.

The full recipe, including `ProxyJump` for clusters that *do* admit SSH to the
compute node, is in
{doc}`../../python_package/d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

## Common Edge Cases

- **No source access at startup**: viewer/server startup establishes the current
  source generation and therefore fails directly if the configured source is
  unreachable. `viewer.ensure_web_ui_cached(force=False)` can verify a
  previously established generation without network access, but it does not
  change startup into an offline substitution path.
- **Notebook blocks iframes** (security policy): you may need to open `viewer.viewer_url` in a new browser tab.
- **Socket allocation failure**: the operating system can reject a new listener
  when process or network resource limits are exhausted.
- **Corporate proxies**: can block cross-origin requests or drop the viewer's WebSocket (hooks and commands); the viewer then says it cannot open its channel.
- **Huge in-memory `AnnData`**: the object is written once to a temporary
  `.h5ad`, which needs disk space about the size of the object. Pass the path
  of a saved `.h5ad` or Zarr store instead, and it is read in place.

(troubleshooting)=
## Troubleshooting (Massive)

This section is intentionally redundant and explicit: it is designed for “I need to fix this now”.

### Symptom: “The viewer doesn’t appear (blank output cell)”

**Likely causes (ordered)**
1) You are not actually running in a Jupyter environment (e.g. plain Python script).
2) The notebook blocks iframes (security policy).
3) The exact viewer generation could not be established from its configured source.

**How to confirm**
- Print the viewer URL:

  ```python
  print(viewer.viewer_url)
  ```

- Open that URL in a normal browser tab.

**Fix**
- Ensure you are using Jupyter/JupyterLab/VSCode notebooks.
- Confirm the kernel/runtime can reach the configured viewer source at startup.
- Correct the reported inventory/object verification error or pass a writable
  `web_cache_dir`.
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

### Symptom: “Port already in use / the port is different”

**Likely causes**
- Old viewers still running from earlier cells.
- Some other process is using the explicit port you requested.

**Fix**
- Call `viewer.stop()` when done.
- Or run `cleanup_all()`.
- Restart the kernel if needed.
- If you need a stable AnnData port, construct
  `AnnDataViewer(..., port=..., dataset_name=..., dataset_id=...)`.

### Symptom: “The viewer says Cellucid found nothing to plot”

**Likely causes**
- No `obsm` array can be plotted: each one is a table, sparse, not one row per
  cell, or non-finite in the columns that would be drawn.

**How to confirm**
- `print(adata.obsm.keys())`
- Read the message the viewer shows: it lists every `obsm` entry and why it
  was not used.

**Fix**
- Compute UMAP (`sc.tl.umap(adata)` writes a two-column `X_umap`, which is read
  as 2D as written), or
- name a usable array explicitly:
  `show_anndata(adata, embedding_keys={"2d": "X_tsne"}, dataset_name=..., dataset_id=...)`.
- If you show a `.h5ad` or Zarr store path rather than an in-memory object,
  re-save the file after a repair.

### Symptom: “The KNN graph is there, but there is no edge overlay”

**Likely causes**
- The store's graph is not named `connectivities`, so the viewer did not pick
  it, or the store's `uns['cellucid']['defaults']` names no graph.
- The graph was refused when it was read (not symmetric, such as `distances`).

**How to confirm**

```python
viewer.wait_for_event("data_mapping")
print(viewer.data_mapping["knn"])
```

`None` means no neighbor graph is chosen.

**Fix**
- Choose the graph under **Session › Data mapping › Neighbour graph**, or
  recreate the viewer with `knn_key="<obsp key>"`.

### Symptom: “The cell runs for a long time before the viewer appears”

**Likely causes**
- You passed an in-memory `AnnData`, and it is being written to a temporary
  `.h5ad`.
- The store holds an array stored with a codec the browser cannot decode, and
  it is being read once and written to a decodable temporary `.h5ad`.
- The viewer generation is being established from its configured source.

**Fix**
- Pass the path of a saved `.h5ad` or Zarr store, which is read in place.
- Save the store without the undecodable codec (gzip, or no compression).

### Symptom: “Vector field overlay is missing (no toggle / no fields)”

**Likely causes**
- Your vectors are not stored under a detected key (e.g. missing `_umap` suffix).
- Shape mismatch: vectors are `(n_cells, 2)` but you are viewing in 3D (or vice versa).
- You have vectors, but they don’t match `adata.n_obs` (filtered/shuffled mismatch).

**How to confirm**
- Inspect the bound fields the viewer reports:

  ```python
  viewer.wait_for_event("data_mapping")
  print(viewer.data_mapping["vectorFields"])
  ```

**Fix**
- Follow the naming rules in {ref}`vector-fields` (and {doc}`07_file_format_expectations`).
- If you have both 2D and 3D points, make sure you provide the matching vector dimension.

### Symptom: “Gene search returns nothing / the wrong gene names”

**Likely causes**
- Your genes were named from `var_names`, which holds a vocabulary you are not searching in.
- The names you want are in a different `var` column.

Each gene carries exactly one name — whichever column you select — and Cellucid
applies no alias mapping, so it is worth selecting the column your readers will
type.

**Fix**
- Pass `gene_id_column="..."` to `show_anndata()`.
- Or write the file with `prepare(var_gene_id_column="...")`, which records it.

### Symptom: “Hooks don’t fire (selection events never reach Python)”

**Likely causes**
- The viewer's channel (`/_cellucid/ws`, a WebSocket) is blocked by a proxy that does not pass WebSockets; the viewer says it cannot open its channel.
- The local server is unreachable from the browser.

**How to confirm**
- Register a raw handler:

  ```python
  @viewer.on_message
  def debug(event):
      print(event)
  ```

- Check `viewer.debug_connection()["channel_connected"]`.
- Open browser devtools → Network, filter `WS`, and look at the `/_cellucid/ws` connection: each selection, hover or click adds a frame.

**Fix**
- Ensure the data server is reachable (`/_cellucid/health` works).
- Behind a proxy, make it pass WebSockets (Jupyter Server Proxy does), or use SSH port forwarding.

### Symptom: “The first gene takes noticeably longer than the next ones”

**Explanation**
- The viewer reads a gene from the store itself. For a row-major (CSR) sparse
  matrix it first builds a gene index over the stored entries, once per
  matrix; every later gene reads only its own entries. With the kernel on the
  same machine the index starts as soon as the dataset opens, so a gene asked
  for within the first seconds waits for it; behind a remote JupyterHub (a
  slower link, as the browser measures it) the index starts with the first
  gene.

**Fix**
- Nothing is needed. To make the first gene as fast as the rest, save the
  matrix column-major (CSC); an in-memory `AnnData` is written that way
  already.

## Next Steps

- Stable dataset identity + sessions: {doc}`06_dataset_identity_why_it_matters`
- File/key expectations (including vector fields): {doc}`07_file_format_expectations`
- General loading troubleshooting: {doc}`08_troubleshooting_data_loading`
- Vector field overlay usage: {doc}`../i_vector_field_velocity/index`
- Hooks deep dive: {doc}`../../python_package/e_jupyter_hooks/index`
- Runnable Pancreas notebook and verified Jupyter captures:
  {doc}`../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery`
- If you want to share datasets publicly: a dataset catalog on GitHub or a static host ({doc}`02_local_demo_tutorial`)
