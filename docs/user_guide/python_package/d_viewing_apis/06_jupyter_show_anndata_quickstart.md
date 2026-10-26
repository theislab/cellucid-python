# Jupyter: `show_anndata()` quickstart

This page shows how to embed Cellucid in a notebook output cell and interact with it from Python.

Use this when you want:
- a tight analysis loop (select cells → compute something → highlight back),
- reproducible exploration inside a notebook,
- and/or to use hooks/events (selection, hover, click).

If you are not in a notebook environment, start with {doc}`04_cli_cellucid_serve_quickstart`.

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste “Minimal cells”, stop at “What success looks like”.
- Computational: focus on how `show_anndata` reads `.h5ad` and Zarr stores in
  place, and remote notebook caveats.
- Power users: focus on HTTPS notebooks, proxies, and `viewer.debug_connection()`.

**Time**
- Minimal embed: ~5 minutes
- Full read (edge cases + troubleshooting): ~20–30 minutes

**Prerequisites**
- `pip install cellucid`
- A notebook environment: classic Jupyter or JupyterLab. VS Code notebooks and
  Colab have not been tested by the project
  ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`)

## Read once: exact web-generation startup (network requirement)

```{important}
Notebook embeds load the viewer UI from the same local server that serves your dataset.

At startup Cellucid establishes the complete source generation declared by
`cellucid-web-assets.json`, including `index.html`, root browser metadata, and
`/assets/*`. Every declared byte is verified before the server publishes the UI.

Startup requires access to the configured source. It never substitutes a stale
cached generation after a source failure.
Configure the generation directory with `web_cache_dir=...`.
```

## Minimal cells (copy/paste)

### A) Show an in-memory `AnnData`

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

### B) Show a file or folder, read by the viewer itself

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A Zarr store (folder or .zip) uses the same call shape, and so do .h5mu,
# .loom, a 10x .h5 or outs/ folder, a MatrixMarket folder and .rds/.rda files.
```

### C) Show a MuData

```python
viewer = show_anndata(mdata, dataset_name="CITE-seq", dataset_id="cite", modality="rna")
```

The object is written once to a temporary Cellucid-ready `.h5mu` (one
modality at a time); `modality=` names the modality the viewer opens on.

### D) Show a Cellucid-ready AnnData (fastest; names itself)

A file written by `cellucid.prepare()` declares its own name and id in
`uns['cellucid']`, so no identity is passed:

```python
from cellucid import show_anndata

viewer = show_anndata("pbmc_cellucid.h5ad", height=600)
```

## What success looks like

- You see an interactive viewer *inside the notebook output cell*.
- You can pan/zoom and select cells.
- The terminal/kernel stays running (the viewer needs the local server).

## Path or object (how to choose)

| You pass | Best for | Tradeoffs |
|---|---|---|
| a path (`.h5ad`, Zarr store, or any other format Cellucid reads) | large data; files you will reopen | read in place, nothing copied |
| an in-memory `AnnData` or `MuData` | the object you are analysing now | written once to a temporary `.h5ad` (`.h5mu`), removed when the viewer stops |
| a file written by `prepare()` | the fastest reads; sharing | written once by `prepare()` |

## How the embed works (so you can debug it)

When you call `show_anndata(...)`, Cellucid:

1) starts a local data server on `http://127.0.0.1:8765` — or, when another
   viewer or program holds 8765, on the next free port, with a printed note
   that Google sign-in (Google Drive data, community annotation) works only on
   8765 — and
2) renders an iframe pointing at the same server, with query params like:

```text
http://127.0.0.1:<port>/?viewer=true&viewerId=<id>&viewerToken=<token>&anndata=true
```

`anndata=true` makes the viewer read the dataset from that server itself,
with range requests where the format allows it.

The viewer then opens its channel, one WebSocket to the same server, on which
Python's commands and the viewer's events (hooks) travel:

```text
ws://127.0.0.1:<port>/_cellucid/ws
```

## HTTPS / remote notebooks: browser reachability

Some notebook environments (JupyterHub, cloud notebooks) serve the notebook over HTTPS or from a remote origin.
In those cases, a direct `http://127.0.0.1:<port>` iframe may be blocked or unreachable.

For a remote notebook, pass the one browser-reachable base URL explicitly with
`client_server_url=...`. Cellucid otherwise uses the local server URL exactly.

Deep dive + fixes: {doc}`10_notebook_widget_mode_advanced`.

## Programmatic control (quick preview)

Commands travel on the viewer's channel. One sent before the viewer has
connected (the cell output is still loading) waits and is delivered, in
order, the moment it connects:

```python
viewer.highlight_cells([0, 1, 2], color="#ff0000")
viewer.set_color_by("cell_type")
viewer.set_visibility([0, 1, 2], visible=False)
viewer.reset_view()
```

Full hooks/events docs: {doc}`../e_jupyter_hooks/index`.

## Cleanup (important for long notebook sessions)

When you’re done:

```python
viewer.stop()
```

This stops the server, closes its port, and removes the temporary `.h5ad`
(`.h5mu`) an in-memory `AnnData` (`MuData`) was written to.

More lifecycle details (ports, multiple viewers): {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`.

## Troubleshooting

Start here:
- {doc}`15_troubleshooting_viewing`

If you’re in a notebook, also run:

```python
viewer.debug_connection()
```
