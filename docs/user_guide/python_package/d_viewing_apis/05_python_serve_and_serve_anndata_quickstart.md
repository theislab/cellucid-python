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

serve_anndata(
    "/path/to/data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# Zarr uses the same call shape and is loaded eagerly by anndata.read_zarr.
```

```{note}
Both convenience functions block while serving. `serve(...)` returns `None`;
`serve_anndata(...)` returns its closed `AnnDataServer` only after serving has
ended. Use the classes below with `start_background()` when you need a running
server object or programmatic lifecycle control.
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
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

## Non-blocking servers (recommended for scripts/tools)

### Exported data: `CellucidServer`

```python
from cellucid import CellucidServer

server = CellucidServer("/path/to/export_dir", open_browser=False, quiet=False)
server.start_background()

print(server.url)        # e.g. http://127.0.0.1:8765
print(server.viewer_url) # e.g. http://127.0.0.1:8765/?source=remote

# ... do other work ...

server.stop()
```

Use `server.viewer_url` when opening the application. Its exact
`source=remote` marker identifies this user-served prepared-data launch before
the first paint, so the ordinary sample-catalog welcome overlay is not shown.
The data and viewer remain on the same server origin. `server.viewer_url` opens
the served catalog exposed by `/_cellucid/datasets`; it does not embed an
arbitrary dataset. If that catalog contains exactly one dataset, the viewer
selects that sole entry. If it contains multiple datasets, the viewer requires
an exact dataset-id selection and never chooses the first entry.

### AnnData: `AnnDataServer`

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

print(server.viewer_url)  # includes ?anndata=true

# ... do other work ...

server.stop()
```

## AnnData-specific options you can pass

These are forwarded to the `AnnDataAdapter` and affect how your AnnData is interpreted:

- `.h5ad` paths are always opened with `backed="r"`; `.zarr` paths are always
  loaded eagerly with `anndata.read_zarr`
- `latent_key`: which `obsm` key to treat as latent space (used for some derived quantities)
- `gene_id_column`: `None` uses `var.index`; any non-blank string names that
  exact `var` column (default `None`)
- `normalize_embeddings`: normalize UMAP coordinates to `[-1, 1]` (default `True`)
- `dataset_name`: required non-empty, unpadded human-readable identity without
  control characters; Unicode is preserved
- `dataset_id`: required portable 1–180 character ASCII identity (letter or
  digit first; then letters, digits, `.`, `_`, or `-`; no trailing `.` or
  reserved Windows device name)
- `vector_field_default`: exact default field ID; required when direct AnnData
  declares more than one vector field

All omitted optional values retain the signature defaults; identity is never
derived, and a multi-field vector declaration has no implicit default.

## Edge cases (high-signal)

- **Remote machine**: prefer an SSH tunnel over binding to `0.0.0.0` (see {doc}`12_remote_servers_ssh_tunneling_and_cloud`).
- **Huge `AnnData` in memory**: direct `serve_anndata(...)` may use substantial RAM;
  prefer a `.h5ad` path for read-only-backed access. A `.zarr` path is loaded
  eagerly.
- **Gene IDs not found**: set `gene_id_column` if `var.index` is not what you search by.

## Troubleshooting

- Server starts but UI shows an error page: {doc}`15_troubleshooting_viewing`
- Data requirements for AnnData mode (UMAP keys, vectors): {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
