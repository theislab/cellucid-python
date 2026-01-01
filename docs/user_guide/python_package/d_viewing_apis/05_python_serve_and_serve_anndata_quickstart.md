# Python: `serve()` and `serve_anndata()` quickstart

This page is for running Cellucid servers **from Python code** (scripts, pipelines, interactive Python).

If you’re in a notebook and want an embedded viewer, go to {doc}`06_jupyter_show_and_show_anndata_quickstart` instead.

## At a glance

**Audience**
- Computational users: scripted workflows and reproducibility.
- Developers: embedding Cellucid into larger systems.

**Time**
- Minimal server: ~5 minutes

**Prerequisites**
- `pip install cellucid`

## Quick start: serve an exported directory

Use this when you already ran `cellucid.prepare(...)` and have an export folder.

```python
from cellucid import serve

serve("/path/to/export_dir")  # blocks until Ctrl+C
```

This starts a local server and opens the viewer in your default browser (unless you disable it).

## Quick start: serve AnnData directly

Use this when you have a `.h5ad`, `.zarr`, or an in-memory `AnnData`.

```python
from cellucid import serve_anndata

serve_anndata("/path/to/data.h5ad")  # backed/lazy by default
# serve_anndata("/path/to/data.zarr")
```

```{note}
`serve_anndata(...)` returns an `AnnDataServer` instance. The convenience `serve(...)` function is blocking and does not return a server object.
If you want lifecycle control (start in background, stop programmatically), use the classes below.
```

## Controlling host/port and browser opening

Both functions accept:

- `port` (default `8765`)
- `host` (default `127.0.0.1`)
- `open_browser` (default `True` in the convenience functions)
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
)
```

## Non-blocking servers (recommended for scripts/tools)

### Exported data: `CellucidServer`

```python
from cellucid import CellucidServer

server = CellucidServer("/path/to/export_dir", open_browser=False, quiet=False)
server.start_background()

print(server.url)        # e.g. http://127.0.0.1:8765
print(server.viewer_url) # e.g. http://127.0.0.1:8765/

# ... do other work ...

server.stop()
```

### AnnData: `AnnDataServer`

```python
from cellucid import AnnDataServer

server = AnnDataServer(
    "data.h5ad",
    open_browser=False,
    quiet=False,
    backed=True,        # default; lazy loading for h5ad
    latent_key="X_pca", # optional; see below
)
server.start_background()

print(server.viewer_url)  # includes ?anndata=true

# ... do other work ...

server.stop()
```

## AnnData-specific options you can pass

These are forwarded to the `AnnDataAdapter` and affect how your AnnData is interpreted:

- `backed` (h5ad only): `True`/`"r"` for lazy; `False` for in-memory
- `latent_key`: which `obsm` key to treat as latent space (used for some derived quantities)
- `gene_id_column`: which `var` column to use as gene IDs (default `"index"`)
- `normalize_embeddings`: normalize UMAP coordinates to `[-1, 1]` (default `True`)
- `dataset_name`, `dataset_id`: override identity shown/used by the viewer

If you’re not sure, the defaults are usually correct.

## Edge cases (high-signal)

- **Remote machine**: prefer an SSH tunnel over binding to `0.0.0.0` (see {doc}`12_remote_servers_ssh_tunneling_and_cloud`).
- **Huge `AnnData` in memory**: `serve_anndata(adata)` may duplicate/copy data paths and use a lot of RAM; prefer `.h5ad` backed mode or `.zarr`.
- **Gene IDs not found**: set `gene_id_column` if `var.index` is not what you search by.

## Troubleshooting

- Server starts but UI shows an error page: {doc}`15_troubleshooting_viewing`
- Data requirements for AnnData mode (UMAP keys, vectors): {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
