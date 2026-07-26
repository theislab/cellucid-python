# Jupyter: `show()` and `show_anndata()` quickstart

This page shows how to embed Cellucid in a notebook output cell and interact with it from Python.

Use this when you want:
- a tight analysis loop (select cells → compute something → highlight back),
- reproducible exploration inside a notebook,
- and/or to use hooks/events (selection, hover, click).

If you are not in a notebook environment, start with {doc}`04_cli_cellucid_serve_quickstart`.

## At a glance

**Audience**
- Wet lab / non-technical: copy/paste “Minimal cells”, stop at “What success looks like”.
- Computational: focus on read-only-backed `.h5ad`, eagerly loaded `.zarr`,
  and remote notebook caveats.
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

### B) Show a `.h5ad` read-only-backed or an eagerly loaded `.zarr`

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A .zarr path uses the same call shape and is loaded eagerly.
```

### C) Show a pre-exported directory (fastest + most reproducible)

```python
from cellucid import show

viewer = show("./export_dir", height=600)
```

## What success looks like

- You see an interactive viewer *inside the notebook output cell*.
- You can pan/zoom and select cells.
- The terminal/kernel stays running (the viewer needs the local server).

## `show()` vs `show_anndata()` (how to choose)

| Function | You pass | Best for | Tradeoffs |
|---|---|---|---|
| `show(export_dir)` | export folder | fastest + reproducible + sharing | requires export step |
| `show_anndata(data, dataset_name=..., dataset_id=...)` | `AnnData` / `.h5ad` / `.zarr` | convenience in analysis | slower than exports; requires exact embedding keys |

## How the embed works (so you can debug it)

When you call `show(...)` or `show_anndata(...)`, Cellucid:

1) starts a local data server (usually `http://127.0.0.1:<port>`), and
2) renders an iframe pointing at the same server, with query params like:

```text
http://127.0.0.1:<port>/?jupyter=true&viewerId=<id>&viewerToken=<token>
```

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

This stops the server and frees file handles/memory (especially important for `.h5ad` backed mode).

More lifecycle details (ports, multiple viewers): {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`.

## Troubleshooting

Start here:
- {doc}`15_troubleshooting_viewing`

If you’re in a notebook, also run:

```python
viewer.debug_connection()
```
