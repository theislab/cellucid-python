# Python: `serve_anndata()` quickstart

This page is for running Cellucid servers **from Python code** (scripts, pipelines, interactive Python).

If you’re in a notebook and want an embedded viewer, go to {doc}`06_jupyter_show_anndata_quickstart` instead.

## At a glance

**Audience**
- Computational users: scripted workflows and reproducibility.
- Developers: embedding Cellucid into larger systems.

**Time**
- Minimal server: ~5 minutes

**Prerequisites**
- `pip install cellucid`

## Quick start

Use this when you have a `.h5ad`, a Zarr store (folder or `.zip`), any other
file or folder Cellucid reads ({doc}`02_loading_routes`), or an in-memory
`AnnData` or `MuData`. The viewer reads the store itself, in place, with range
requests where the format allows it; the server serves its bytes. An in-memory
object is written once to a temporary `.h5ad` (`.h5mu` for a MuData), removed
when the server stops.

```python
from cellucid import serve_anndata

serve_anndata(
    "/path/to/data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A Zarr store (folder or .zip) or an in-memory AnnData uses the same call shape.
```

This starts a local server and opens the viewer in your default browser
(unless `open_browser=False`). A file written by `cellucid.prepare()` declares
its own name and id in `uns['cellucid']`, so `dataset_name` and `dataset_id`
are needed only for other AnnData.

```{note}
`serve_anndata(...)` blocks while serving and returns its closed
`AnnDataServer` only after serving has ended. Use the class below with
`start_background()` when you need a running server object or programmatic
lifecycle control.
```

## What a start prints

Four numbered steps — `Detecting format`, `Reading metadata`,
`Checking launch choices`, `Starting server` — then the banner with the Viewer
URL. No step reads an array.

Full transcript: {doc}`../../web_app/b_data_loading/04_server_tutorial`.

## Controlling host/port and browser opening

`serve_anndata` and `AnnDataServer` accept:

- `port` (default `8765`, bound exactly: a taken port raises; `0` asks the
  operating system for a free one)
- `host` (default `127.0.0.1`)
- `open_browser` (default `True` for `serve_anndata`, `False` for
  `AnnDataServer`)
- `quiet` (default `False`)

Example:

```python
from cellucid import serve_anndata

serve_anndata(
    "data.h5ad",
    port=9000,
    host="127.0.0.1",
    open_browser=False,
    quiet=False,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

### Binding to every interface (`host="0.0.0.0"`)

Use it when the client cannot reach this machine directly — typically an HPC
compute node whose tunnel terminates on a login node. Two things to know:

- `url` and `viewer_url` report the **loopback** origin, because the wildcard is
  not an address anything can open. `network_urls` carries the machine's own name
  for clients elsewhere on the network, and the banner prints both.
- Through a tunnel you still open `http://127.0.0.1:<port>/…`: the near end of
  the tunnel is your own machine.

A non-loopback bind has no authentication of any kind. See
{doc}`13_security_privacy_cors_and_networking` and
{doc}`17_hpc_slurm_and_compute_node_serving`.

## Non-blocking servers (recommended for scripts/tools)

### `AnnDataServer`

```python
from cellucid import AnnDataServer

server = AnnDataServer(
    "data.h5ad",
    open_browser=False,
    quiet=False,
    latent_key="X_pca", # optional; see below
    dataset_name="My study",
    dataset_id="my-study-v1",
)
server.start_background()

print(server.url)         # e.g. http://127.0.0.1:8765
print(server.viewer_url)  # e.g. http://127.0.0.1:8765/?anndata=true

# ... do other work ...

server.stop()
```

Use `server.viewer_url` when opening the application. Its exact `anndata=true`
marker identifies this user-served launch before the first paint, so the
ordinary sample-catalog welcome overlay is not shown, and the viewer opens the
one dataset the server lists at `/_cellucid/datasets`. The data and viewer
remain on the same server origin.

## Options you can pass

Every option is optional. Each one you set is checked against the store's
metadata before the server binds (a name the store does not have is an error;
for a format Python does not read — an R object, Loom, 10x, MatrixMarket — the
viewer checks it when it opens the data) and handed to the viewer as a launch
default. One you leave unset is taken from
the store's own `uns['cellucid']['defaults']`, else the viewer's default rule
(see {doc}`08_anndata_mode_show_anndata_and_serve_anndata`).

- `dataset_name`: non-empty, unpadded human-readable identity without control
  characters; Unicode is preserved. Defaults to `uns['cellucid']['name']`;
  required when the store has none
- `dataset_id`: portable 1–180 character ASCII identity (letter or digit first;
  then letters, digits, `.`, `_`, or `-`; no trailing `.` or reserved Windows
  device name). Defaults to `uns['cellucid']['id']`; required when the store
  has none
- `latent_key`: the `obsm` latent space (scVI, PCA) the outlier filter
  measures in, never a 2D/3D layout. The viewer reads the outlier quantiles
  the store holds for it, or computes them on first use
  ({doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`)
- `normalize_embeddings`: whether the viewer scales embeddings into its
  coordinate range
- `centroid_min_points`: the fewest cells a category needs to get outlier
  quantiles (10 when None)
- `obs_keys`: exact `adata.obs` columns the viewer shows, in this order; every
  column is shown when omitted. Use it to leave out a column Cellucid cannot
  show, such as a `datetime64` collection date
- `vector_field_default`: the vector field the overlay opens on, by its field
  id (`velocity` for `obsm['velocity_umap']`)
- `embedding_keys`: the `obsm` array of a dimension slot, such as
  `{"2d": "X_tsne"}` or `{"2d": ("X_pca", [0, 1])}` with 0-based columns; see
  below
- `knn_key`: the `obsp` matrix the viewer draws as the neighbor graph, such as
  `"connectivities"`; read only when a view draws it
- `gene_id_column`: the column of the gene axis that names genes (default: its
  index); empty values take the index name and repeated names are made unique
  (`MATR3`, `MATR3-1`)
- `expression_source`: `"X"`, `"raw"`, or `"layer:<name>"` (default: `X`, or
  without `X` the guess described in
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata`)
- `modality`: for a multimodal dataset (a MuData), the modality the viewer
  opens on (one of `mdata.mod`); when None, the first gene-expression modality
  (`rna`, `gex`, …), else the first. Refused for data with one feature space

### Which `obsm` key becomes the embedding

A slot named in `embedding_keys` shows that array. Every other slot follows the
default rule: `X_umap_1d` / `X_umap_2d` / `X_umap_3d` each fill their own
dimension, and with none of them present, a plain `X_umap` of 1, 2, or 3 columns
fills the slot of its own column count, so a Scanpy object is shown as written.
Only when neither fills any slot does the viewer make a deterministic best guess
from the other `obsm` arrays, and marks it as guessed.

```python
server = cellucid.serve_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    embedding_keys={"2d": ("X_pca", [0, 1])},
)
```

### The neighbor graph: `knn_key=`

```python
server = cellucid.serve_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

Without `knn_key` the viewer uses `obsp['connectivities']` when the store has
it. Either way the graph is read only when a view draws it, so it costs nothing
at startup. Naming a graph the store does not carry stops the server before it
binds; one that is not symmetric (such as `distances`) is refused by the viewer
when it reads it.

Both rules in full:
{doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## Edge cases (high-signal)

- **Remote machine**: prefer an SSH tunnel over binding to `0.0.0.0` (see {doc}`12_remote_servers_ssh_tunneling_and_cloud`).
- **Huge `AnnData` in memory**: an in-memory object is written once to a
  temporary `.h5ad` before the server binds; a path is served in place with no
  copy, so prefer a path when the object is already on disk.
- **Gene IDs not found**: set `gene_id_column` if `var.index` is not what you search by.

## Troubleshooting

- Server starts but UI shows an error page: {doc}`15_troubleshooting_viewing`
- Data requirements (UMAP keys, vectors): {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
