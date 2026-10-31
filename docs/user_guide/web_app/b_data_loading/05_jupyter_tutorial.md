# Jupyter Integration (Notebook Embedding)

This tutorial shows how to load data into Cellucid **from inside a notebook**
(JupyterLab, Jupyter Notebook 7, VS Code, Google Colab, marimo, Databricks).

You will learn:
- how to display an **AnnData or MuData object, an `.h5ad`, a Zarr store, or any
  other file or folder Cellucid reads** with `show_anndata()`
- how to work with **vector fields** (velocity/drift overlays) in notebooks
- how the integration works under the hood (so you can debug it)
- how to drive the viewer from Python (highlight/color/visibility), react to UI
  events (hooks), and pull what you made in the viewer back onto your AnnData

If you are not in a notebook environment, start with {doc}`04_server_tutorial`.
From R, see {doc}`../../r_package/d_viewing_loading/08_view_in_jupyter_notebooks`.

:::{tip}
Prefer to begin from a real, runnable notebook? The
{doc}`Pancreas Jupyter walkthrough
<../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery>`
links the checked notebook and its exact launch command.
:::

## At A Glance

**Audience**
- Wet lab / beginner: copy/paste the “Minimal cells” sections and focus on “What success looks like”.
- Computational users: focus on how the viewer reads a store in place, dataset
  sizes, vector fields, and pulling results back.
- Power users: focus on what happens between cells, hooks, and debugging.

**Time**
- Minimal working embed: ~5 minutes
- Full read (hooks + troubleshooting): ~20–30 minutes

**Prerequisites**
- A notebook that shows widgets: JupyterLab (tested by the project), Jupyter
  Notebook 7, VS Code, Google Colab, marimo or Databricks (expected to work,
  not tested; see
  {doc}`../../python_package/e_jupyter_hooks/03_supported_environments_matrix`)
- `pip install cellucid` (it installs `anywidget`; JupyterLab and Notebook 7
  need it where the *notebook server* runs, the same environment in most
  setups)

```{important}
**Network requirement:** the viewer is the Cellucid web app, which **your
browser** loads from `https://www.cellucid.com` into the cell's output (or
from a copy you host, named with `web_source_url=`). The kernel needs no
network access and opens no port: the app reads your data from the kernel
through the notebook's own connection, so nothing is uploaded anywhere. The
widget's page is a viewer page, which sends no usage analytics.
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
```

In marimo, end the cell with `viewer` (marimo shows the widget a cell
evaluates to); every other notebook shows it at once.

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

### Google Colab

```ipython
!pip -q install cellucid
```

then the same `show_anndata(...)` cell. Colab renders anywidgets; no port,
proxy or `google.colab` call is involved.

### What success looks like

```{figure} ../../../_static/screenshots/jupyter/show-anndata-driven-from-python.png
:alt: The Cellucid viewer driven from Python, after viewer.set_color_by("clusters") and viewer.highlight_cells() on the Ductal cells: the Pancreas sample coloured by clusters with the Ductal cells highlighted.
:width: 1440px

The viewer after Python called `viewer.set_color_by("clusters")` and
highlighted the Ductal cells. This capture was taken in a browser tab; in a
notebook the same app sits *inside* the output cell, and scrolls with the
notebook.
```

Three things confirm the call worked:

1. an interactive viewer appears in the output cell (not an empty area);
2. in a later cell, `viewer.state.ready` is no longer `None`;
3. the cell and dimension counts it reports match the data you passed.

## How It Works (Mental Model)

When you call `show_anndata(...)`, Cellucid:

1) **Checks the dataset and your launch choices.** It reads only metadata (the
   cell and gene counts and the names of `obs`, `var`, `obsm`, `obsp` and the
   layers) and checks every option you named against it. An in-memory
   `AnnData` or `MuData` is served as a store with its cells in priority
   order, made from it as the viewer reads it (no copy of it is written as a
   file).
2) **Shows a widget** (`cellucid.Widget`, an [anywidget](https://anywidget.dev))
   in the cell's output. The widget frames the app at
   `https://www.cellucid.com/?viewer=embedded` and hands it, in a message, the
   viewer's id and token and a private channel.
3) **Answers the viewer's reads.** The viewer's own HDF5/Zarr reader reads the
   dataset in place with range requests, only what a view needs; the widget
   relays each request to the kernel, which answers with the bytes on disk,
   over the notebook's own connection (a Jupyter comm).

No server is started and no port is opened, so the same cell works on a
laptop, on a remote JupyterHub, in VS Code Remote and in Colab.

### Between cells

The kernel handles the widget's messages **between cells**, as for every
notebook widget:

- the viewer loads once the cell that showed it finishes;
- hooks (`@viewer.on_selection`, …) run, and `viewer.state` updates, whenever
  the kernel is idle;
- a call that waits for the viewer's answer (`wait_for_ready()`,
  `wait_for_event()`, a session request's `result()`) cannot be answered
  inside the cell that waits, so it raises at once and says what to do instead;
- commands (`highlight_cells`, `set_color_by`, …) are sent at once, or
  delivered in order as soon as the viewer connects.

### A browser tab instead

`cellucid.AnnDataViewer(...)` is the same API for a viewer in a browser tab of
its own: it starts a local server (port 8765, or the next free one) that
serves the app and the data, and gives you `viewer.viewer_url`. Use it from a
plain script or through an SSH tunnel (see {ref}`remote-hpc`).

## `show_anndata()` (AnnData / `.h5ad` / Zarr store)

`show_anndata()` is the fastest way to get started in an analysis notebook.

It supports:
- an in-memory `AnnData` or `MuData`, served as a store with its cells in
  priority order (its arrays copied into memory-mapped files removed when the
  viewer stops: the object as it was shown, its expression matrices laid out
  gene-major in the background while the viewer opens, read-only until then;
  each chunk made from that copy when the viewer asks for it)
- a `.h5ad` file or a Zarr store folder (Zarr v2 with `.zgroup`, or v3 with
  `zarr.json`) or zipped Zarr store (`.zip`, such as `.zarr.zip`): one
  `prepare()` wrote read as it is on disk, any other served in priority order,
  read from the file
- a `.h5mu` or MuData Zarr store, a `.loom` file or a 10x `outs/` folder:
  one `prepare()` wrote (a MuData) read as it is, any other served in priority
  order as the browser presents the file itself (a MuData Zarr store for a
  MuData, or for a 10x matrix of several feature types), read from the file
- a MatrixMarket folder or an `.rds` / `.rda`, read as it is on disk, once,
  by the browser ({doc}`07_file_format_expectations`)

Nothing is converted ahead: the viewer reads the store itself, with range
requests, and a subsample it opens (Session › Open data) reads only its cells,
the first rows of every array in priority order. A store holding an array
stored with a codec the browser cannot decode (an HDF5 filter such as
szip or scale-offset; a Zarr codec such as Snappy or zfpy) is decoded on the
kernel, chunk by chunk as the viewer asks.

### Exact parameters

`show_anndata(...)` has a closed signature. It accepts `data` and `height`,
followed by these exact keyword-only parameters:

- `web_source_url`: where the browser loads the Cellucid app from
  (`https://www.cellucid.com` by default; the address of a copy you host,
  without a trailing slash)
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
widget is shown — a key that does not exist raises at once — and passed to the
viewer as a launch default. An option left at `None` is taken from the store's
own `uns['cellucid']['defaults']`, else from the viewer's default rule (see
{doc}`13_data_mapping`).

```{tip}
Use the same exact `dataset_id` whenever you reopen the same dataset generation.
```

```python
from cellucid import show_anndata

# With options
viewer = show_anndata(
    "data.h5ad",
    height=700,
    latent_key="X_pca",
    gene_id_column="gene_symbols",
    dataset_name="PBMC demo",
    dataset_id="pbmc-demo",
)
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

Naming a matrix the store does not have raises before the viewer is shown:

```text
knn_key names obsp['connectivities'], which this store does not have. obsp holds: []. Compute the neighbour graph with sc.pp.neighbors(adata), or leave knn_key unset.
```

A graph is checked when it is read: it must be symmetric, so the `distances`
matrix `sc.pp.neighbors()` also writes is refused. Any other `obsp` graph can be
chosen later under **Session › Data mapping › Neighbour graph**
({doc}`13_data_mapping`).

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

Once the viewer shows the data, the data mapping it reports lists every bound
field. In a later cell:

```python
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

# Show a column you computed as a field, and a gene set
viewer.push(adata.obs["leiden"], key="leiden_notebook")
viewer.push_gene_set("T cell markers", ["CD3E", "CD3D", "CD2"])
```

```{important}
Cell indices refer to the **row order** of the dataset Cellucid shows.

- For a viewer constructed directly from `adata`, this is the current
  `adata` row order.
- If you subset/shuffle `adata` in Python, the indices will change.
```

### Read and change the data mapping

The *data mapping* says which stored arrays the viewer plots and reads: the
embedding of each dimension slot, the bound vector fields, and for an AnnData
source the neighbour graph, the gene names, the expression matrix and the
latent space. The viewer reports it after the dataset loads, and
`viewer.data_mapping` holds the latest report (read it in a later cell):

```python
viewer.data_mapping["embeddings"]
# {'1d': None, '2d': {'key': 'X_umap', 'columns': [0, 1]}, '3d': None}
viewer.data_mapping_origin   # "default", "best_guess", "session", "user" or "notebook"
```

To change it, name the arrays at launch (`embedding_keys=`, `knn_key=`,
`expression_source=`, …) or use **Session › Data mapping** in the viewer.
`viewer.set_data_mapping(...)` replaces the top-level fields you name, checks
the result and sends it without waiting -- at once, or, in a cell run before
the viewer reported its mapping, as soon as it does -- in a notebook widget
and a browser-tab viewer alike:

```python
viewer.set_data_mapping(
    embeddings={
        "enabled": [{"key": "X_pca", "columns": [0, 1]}],
        "defaults": {"1d": None, "2d": {"key": "X_pca", "columns": [0, 1]}, "3d": None},
    },
)
```

A field this kind of mapping does not have raises `TypeError` and an invalid
value raises `ValueError`; the applied mapping is `viewer.data_mapping` in a
later cell, and a mapping the viewer refuses — an array that fails its
load-time check, for example — arrives as a `command_error` with the
viewer's reason.

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

Hooks run when the kernel is idle (between cells).

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

## Pulling Results into Python

Hooks are great when you want reactive code. For “pull-style” workflows, read
what the viewer holds in a later cell.

### Everything you made, in one call

Every highlight group, label field, cluster, score and gene set you made in
the viewer is sent to Python as you make it. One call writes them all onto
your AnnData:

```python
summary = viewer.pull(adata)
summary.columns        # obs columns written: cellucid_highlight__<id>, your labels, scores
adata.uns["cellucid"]["gene_sets"]
```

Related: `viewer.selection` (the last selection, with obs names),
`viewer.subset(adata, "my group")` (a highlight group as a new AnnData),
`viewer.gene_sets()`. Full reference:
{doc}`../../python_package/g_api_reference_coverage/api/round_trip`.

### `viewer.state` (live snapshot)

`viewer.state` is a small, thread-safe snapshot of the latest events:

```python
print(viewer.state.ready)      # the ready event (or None before the data opened)
print(viewer.state.selection)  # last selection event (or None)
print(viewer.state.hover)      # last hover event (or None)
print(viewer.state.click)      # last click event (or None)
```

### Session bundle (durable saved state)

In a notebook, ask for the current `.cellucid-session` in one cell:

```python
request = viewer.request_session_bundle()
```

and collect it as a Python object in the next:

```python
bundle = request.result()
bundle.save("my-study.cellucid-session")

# Apply to AnnData (adds obs columns; stores metadata in adata.uns["cellucid"])
adata2 = bundle.apply_to_anndata(
    adata,
    expected_dataset_id="my-study-v1",
    inplace=False,
)
```

```{important}
Session application names cells by **original position** (each cell's row in
the order the dataset's cells were given, whatever order a file stores them in).

Only apply a session to the `AnnData` that produced it, in that original order
(a Cellucid-ready file read with anndata: `cellucid.restore_original_order(adata)` first).
```

### Debugging: `viewer.debug_connection()`

If hooks or results seem “stuck”, run in a new cell:

```python
report = viewer.debug_connection()
report
```

It waits for nothing and reports where the browser loads the app from
(`web_source_url`), which outputs are connected (`connected_outputs`),
whether the viewer opened its channel (`channel_connected`), the dataset it
serves, the viewer's `ready` state, and the recent events by type with the
commands the viewer refused (`recent_events["command_errors"]`).

`viewer.debug_connection()` returns a dictionary; it prints nothing by itself.

## Cleanup (Do This If You Re-run Cells Often)

Each viewer holds its dataset (for an in-memory object, the store it serves
and the files a cell-major `X` was laid out in) until it is stopped.

```python
viewer.stop()  # freeze the output on its last view; release what it served
```

```{note}
`viewer.stop()` requests an exact freeze, stops answering the viewer,
unregisters hooks, and reports any failure. After a successful stop, the last
frame remains visible but is non-interactive.
```

### Stop everything created in this kernel

```python
from cellucid.jupyter import cleanup_all

cleanup_all()
```

(remote-hpc)=
## Remote / HPC Notebooks

A notebook whose kernel runs elsewhere (a JupyterHub, a notebook on an HPC
compute node, VS Code Remote-SSH, Colab, a cloud VM) needs **nothing extra**:
`show_anndata()` reaches the kernel through the notebook's own connection,
whichever way your browser reaches the notebook.

You need SSH forwarding only for a **browser-tab** viewer (`AnnDataViewer`,
`serve_anndata`, `cellucid serve`) started on a remote machine:

1) Start the viewer with a fixed port:

   ```python
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       "data.h5ad",
       port=8765,
       dataset_name="My study",
       dataset_id="my-study-v1",
   )
   print(viewer.viewer_url)
   ```

2) On your laptop, forward that port to the remote machine:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now, when your browser loads `http://127.0.0.1:8765/?viewer=true&...`,
it hits your laptop’s `127.0.0.1:8765`, which SSH forwards to the remote server.

```{important}
Without `port=`, `AnnDataViewer` takes 8765, or the next free port when another
viewer holds it. A tunnel forwards one fixed port, so name it with `port=`: a
viewer given a port uses that port or fails.
```

On an HPC compute node reached only through a login node, the forward is
opened *by the login node* (`ssh -N -L 8765:compute-node:8765 you@login-node`),
so the server must accept a connection from another machine; `AnnDataViewer`
always binds loopback, so start a server from the kernel instead:

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

The address you open stays `http://127.0.0.1:8765/?anndata=true`. The full
recipe is in
{doc}`../../python_package/d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

## Common Edge Cases

- **The browser cannot reach `https://www.cellucid.com`** (offline laptop,
  blocked address): the output says the app did not load. Host a copy of the
  app and pass `web_source_url=`.
- **The notebook does not show widgets**: install `anywidget` where the
  notebook server runs, or use a front end that shows widgets.
- **Showing one widget twice** (ending another cell with `viewer`): the newest
  output takes the viewer over; the older one says it was opened elsewhere.
- **Kernel restarts**: the output keeps its last view and says the notebook
  stopped it. Re-run the cell.
- **Huge in-memory `AnnData`**: its arrays are copied once into memory-mapped
  files (disk for one uncompressed copy of them, never a second copy in
  memory; a sparse `X` laid out gene-major), so later edits to the object do
  not show. Open a subsample (Session › Open data) to read only
  some of the cells.

(troubleshooting)=
## Troubleshooting

### Symptom: “The viewer doesn’t appear (empty output cell)”

**Likely causes (ordered)**
1) The notebook front end does not show widgets (`anywidget` missing where
   the notebook server runs).
2) You are not in a notebook: `show_anndata()` raises and names the
   browser-tab alternatives.
3) In marimo, the widget is not the cell's output.

**Fix**
- Install `anywidget` in the notebook server's environment and reload the page.
- In marimo, end the cell with `viewer`.
- Outside a notebook, use `AnnDataViewer(...)` or `cellucid serve`.

### Symptom: “The output says the app did not load”

**Likely cause:** the browser cannot load `https://www.cellucid.com`.

**Fix:** allow it, or host a copy of the app and pass `web_source_url=`.

### Symptom: “A call raised: … a notebook widget hears the viewer only between cells”

`wait_for_ready()`, `wait_for_event()` and a session request's `result()`
wait for the viewer's answer, which cannot arrive while the cell runs. Read
`viewer.state` or `viewer.results` in a later cell, register a hook, or for a
session collect `request.result()` (or
`viewer.apply_session_to_anndata(adata, request)`) in the next cell.

### Symptom: “The viewer says Cellucid found nothing to plot”

**Likely causes**
- No `obsm` array can be plotted: each dense one has another row count than
  the cells, places no cell (every row non-finite), or holds a value past
  float32. (An `obsm` with no dense numeric array at all is refused as having
  no cell embedding: {ref}`no-cell-embedding`.)

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

**How to confirm** (in a later cell)

```python
print(viewer.data_mapping["knn"])
```

`None` means no neighbor graph is chosen.

**Fix**
- Choose the graph under **Session › Data mapping › Neighbour graph**, or
  recreate the viewer with `knn_key="<obsp key>"`.

### Symptom: “The cell runs for a long time before the viewer appears”

**Likely causes**
- You passed an in-memory `MuData`, and it is being written to a temporary
  `.h5mu`.
- You passed an in-memory `AnnData` with a large cell-major `X`, and the first
  gene waits for its gene-major layout (made in the background from the
  start).

**Fix**
- Pass the path of a saved `.h5ad` or Zarr store, which is read in place.
- Save the store without the undecodable codec (gzip, or no compression).

### Symptom: “Vector field overlay is missing (no toggle / no fields)”

**Likely causes**
- Your vectors are not stored under a detected key (e.g. missing `_umap` suffix).
- Shape mismatch: vectors are `(n_cells, 2)` but you are viewing in 3D (or vice versa).
- You have vectors, but they don’t match `adata.n_obs` (filtered/shuffled mismatch).

**How to confirm** (in a later cell)

```python
print(viewer.data_mapping["vectorFields"])
```

**Fix**
- Follow the naming rules in {ref}`vector-fields` (and {doc}`07_file_format_expectations`).
- If you have both 2D and 3D points, make sure you provide the matching vector dimension.

### Symptom: “Gene search returns nothing / the wrong gene names”

**Likely causes**
- Your genes were named from `var_names`, which holds a vocabulary you are not searching in.
- The names you want are in a different `var` column.

**Fix**
- Pass `gene_id_column="..."` to `show_anndata()`.
- Or write the file with `prepare(var_gene_id_column="...")`, which records it.

### Symptom: “Hooks don’t fire (selection events never reach Python)”

**Likely causes**
- A cell is still running: the widget's events arrive when the kernel is idle.
- The selection was not confirmed.
- The output belongs to an older viewer (re-run cells, a kernel restart).

**How to confirm**
- Register a raw handler:

  ```python
  @viewer.on_message
  def debug(event):
      print(event)
  ```

- Check `viewer.debug_connection()["channel_connected"]` in a new cell.

**Fix**
- Let the running cell finish; re-run the cell that creates the viewer for a
  stale output.

### Symptom: “The first gene takes noticeably longer than the next ones”

**Explanation**
- The viewer reads a gene from the store itself. For a row-major (CSR) sparse
  matrix it first builds a gene index over the stored entries, once per
  matrix; every later gene reads only its own entries.

**Fix**
- Nothing is needed. To make the first gene as fast as the rest, save the
  matrix column-major (CSC); an in-memory `AnnData` is served that way
  already.

## Next Steps

- Stable dataset identity + sessions: {doc}`06_dataset_identity_why_it_matters`
- File/key expectations (including vector fields): {doc}`07_file_format_expectations`
- General loading troubleshooting: {doc}`08_troubleshooting_data_loading`
- Vector field overlay usage: {doc}`../i_vector_field_velocity/index`
- Hooks deep dive: {doc}`../../python_package/e_jupyter_hooks/index`
- Runnable Pancreas notebook:
  {doc}`../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery`
- If you want to share datasets publicly: a dataset catalog on GitHub or a static host ({doc}`02_local_demo_tutorial`)
