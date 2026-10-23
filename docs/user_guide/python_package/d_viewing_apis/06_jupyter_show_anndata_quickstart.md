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
- A notebook environment (classic, JupyterLab, VSCode notebooks, Colab)

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

### B) Show a `.h5ad` or a Zarr store, read in place

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A Zarr store (folder or .zip) uses the same call shape.
```

### C) Show a Cellucid-ready AnnData (fastest; names itself)

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
| a path (`.h5ad`, Zarr store) | large data; files you will reopen | read in place, nothing copied |
| an in-memory `AnnData` | the object you are analysing now | written once to a temporary `.h5ad`, removed when the viewer stops |
| a file written by `prepare()` | the fastest reads; sharing | written once by `prepare()` |

## How the embed works (so you can debug it)

When you call `show_anndata(...)`, Cellucid:

1) starts a local data server (usually `http://127.0.0.1:<port>`), and
2) renders an iframe pointing at the same server, with query params like:

```text
http://127.0.0.1:<port>/?jupyter=true&viewerId=<id>&viewerToken=<token>&anndata=true
```

`anndata=true` makes the viewer read the AnnData store from that server
itself, with range requests.

Frontend → Python events (hooks) are delivered via HTTP POST to:

```text
http://127.0.0.1:<port>/_cellucid/events
```

## HTTPS / remote notebooks: browser reachability

Some notebook environments (JupyterHub, cloud notebooks) serve the notebook over HTTPS or from a remote origin.
In those cases, a direct `http://127.0.0.1:<port>` iframe may be blocked or unreachable.

For a remote notebook, pass the one browser-reachable base URL explicitly with
`client_server_url=...`. Cellucid otherwise uses the local server URL exactly.

Deep dive + fixes: {doc}`10_notebook_widget_mode_advanced`.

## Programmatic control (quick preview)

These work after the viewer is displayed:

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

This stops the server, closes its port, and removes the temporary `.h5ad` an
in-memory `AnnData` was written to.

More lifecycle details (ports, multiple viewers): {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`.

## Troubleshooting

Start here:
- {doc}`15_troubleshooting_viewing`

If you’re in a notebook, also run:

```python
viewer.debug_connection()
```
