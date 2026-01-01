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
- Computational: focus on `.h5ad` backed mode, `.zarr`, and remote notebook caveats.
- Power users: focus on HTTPS notebooks, proxies, and `viewer.debug_connection()`.

**Time**
- Minimal embed: ~5 minutes
- Full read (edge cases + troubleshooting): ~20–30 minutes

**Prerequisites**
- `pip install cellucid`
- A notebook environment (classic, JupyterLab, VSCode notebooks, Colab)

## Read once: hosted-asset proxy + caching (network requirement)

```{important}
Notebook embeds load the viewer UI from the same local server that serves your dataset.

If the UI assets are not already cached, Cellucid will download `index.html` + `/assets/*`
from `https://www.cellucid.com` and cache them on disk (one-time per web build).

Offline use is supported *after* you have a cached copy.
Configure the cache directory with `CELLUCID_WEB_PROXY_CACHE_DIR`.
```

## Minimal cells (copy/paste)

### A) Show an in-memory `AnnData`

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=600)
viewer  # (optional) display again in some notebook UIs
```

### B) Show a `.h5ad` / `.zarr` (recommended for large datasets)

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", height=600)  # backed (lazy) by default
# viewer = show_anndata("data.zarr", height=600)
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

<!-- SCREENSHOT PLACEHOLDER
ID: notebook-embed-success
Suggested filename: web_app/notebook_01_embed-success.png
Where it appears: Python Package Guide → Viewing APIs → Jupyter quickstart → What success looks like
Capture:
  - UI location: notebook output area
  - State prerequisites: viewer rendered; dataset loaded; side panel visible
  - Action to reach state: run `viewer = show_anndata(...)`
Crop:
  - Include: the iframe area + a small amount of notebook context (cell input + output)
  - Exclude: private dataset names, browser tabs, personal info
Annotations:
  - Callouts: (1) the embedded viewer frame, (2) the dataset name / point count, (3) the URL bar is optional
Alt text:
  - Notebook cell output showing an embedded Cellucid viewer.
Caption:
  - A successful notebook embed: the viewer appears as an iframe and is backed by a local Cellucid server.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a successful Cellucid notebook embed.
:width: 100%

A successful Cellucid embed inside a notebook output cell.
```

## `show()` vs `show_anndata()` (how to choose)

| Function | You pass | Best for | Tradeoffs |
|---|---|---|---|
| `show(export_dir)` | export folder | fastest + reproducible + sharing | requires export step |
| `show_anndata(data)` | `AnnData` / `.h5ad` / `.zarr` | convenience in analysis | slower than exports; requires UMAP keys |

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

## HTTPS / remote notebooks: the “proxy required” case

Some notebook environments (JupyterHub, cloud notebooks) serve the notebook over HTTPS or from a remote origin.
In those cases, a direct `http://127.0.0.1:<port>` iframe may be blocked or unreachable.

Cellucid will try to use **jupyter-server-proxy** automatically (recommended).

If it cannot, you may see an in-iframe error message telling you to install/enable `jupyter-server-proxy`
or set `CELLUCID_CLIENT_SERVER_URL`.

<!-- SCREENSHOT PLACEHOLDER
ID: notebook-proxy-required-error
Suggested filename: web_app/notebook_02_proxy-required.png
Where it appears: Python Package Guide → Viewing APIs → Jupyter quickstart → HTTPS / remote notebooks
Capture:
  - UI location: embedded iframe
  - State prerequisites: notebook served over HTTPS; jupyter-server-proxy missing/disabled
  - Action to reach state: run `show_anndata(...)` in a JupyterHub/HTTPS notebook without jupyter-server-proxy
Crop:
  - Include: the error message in the iframe
Alt text:
  - Embedded Cellucid iframe showing a “notebook proxy required” message.
Caption:
  - When the notebook is served from HTTPS/remote origins, the embed may require jupyter-server-proxy (or an explicit client server URL) to avoid mixed-content and reach the kernel’s server.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Cellucid “proxy required” iframe message.
:width: 100%

Notebook embeds may require `jupyter-server-proxy` in HTTPS/remote environments.
```

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
