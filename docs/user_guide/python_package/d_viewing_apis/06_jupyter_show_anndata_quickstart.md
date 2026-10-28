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
- Power users: focus on what happens while a cell runs, and `viewer.debug_connection()`.

**Time**
- Minimal embed: ~5 minutes
- Full read (edge cases + troubleshooting): ~20–30 minutes

**Prerequisites**
- `pip install cellucid` (it installs `anywidget`)
- A notebook that shows widgets: JupyterLab, Jupyter Notebook 7, VS Code,
  Google Colab, marimo or Databricks
  ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`)

## Read once: where the viewer comes from

```{important}
`show_anndata(...)` shows a notebook widget (an anywidget) in the cell's
output. The widget frames the Cellucid app, which **your browser** loads from
`https://www.cellucid.com` (or the copy you name with `web_source_url=`), and
connects it to the kernel through the notebook's own connection: the viewer
reads your dataset from the kernel through the widget, so no port is opened
and nothing is uploaded anywhere.

The kernel handles the widget's messages between cells: the viewer loads once
the cell that showed it finishes.
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
```

In marimo, make the widget the cell's output (end the cell with `viewer`);
the other notebooks show it at once.

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
- The kernel stays running (the viewer reads the dataset through it).

## Path or object (how to choose)

| You pass | Best for | Tradeoffs |
|---|---|---|
| a path (`.h5ad`, Zarr store, or any other format Cellucid reads) | large data; files you will reopen | read in place, nothing copied (an AnnData file not in priority order is served in priority order, from the file) |
| an in-memory `AnnData` or `MuData` | the object you are analysing now | an AnnData served in priority order, made from it as the viewer reads it; a MuData written once to a temporary `.h5mu`, removed when the viewer stops |
| a file written by `prepare()` | the fastest reads; sharing | written once by `prepare()` |

## How the embed works (so you can debug it)

When you call `show_anndata(...)`, Cellucid:

1) checks the dataset and the launch choices you passed (an in-memory AnnData
   is served in priority order, made from it as the viewer reads it; a MuData
   is written once to a temporary `.h5mu`), and
2) shows a {class}`~cellucid.jupyter.Widget` whose output frames the app at
   `https://www.cellucid.com/?viewer=embedded`. The page asks the widget for
   its credentials and a private channel (a `MessagePort`), then reads the
   dataset's bytes and opens the viewer channel through it; the widget relays
   both to the kernel over the notebook's own connection.

No server is started and no port is opened, so the same cell works on a
laptop, on a remote JupyterHub, in VS Code and in Colab. Details:
{doc}`10_notebook_widget_mode_advanced` and
{doc}`../e_jupyter_hooks/05_architecture_viewer_channel`.

A viewer in a browser tab of its own instead (a plain script, an SSH tunnel)
is {class}`~cellucid.AnnDataViewer`: the same API, served from this machine.

## Programmatic control (quick preview)

Commands travel on the viewer's channel. One sent before the viewer has
connected (the viewer loads once the cell finishes) waits and is delivered,
in order, the moment it connects:

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

This freezes the output on its last view, stops answering it, and releases what
it served: the store of an in-memory `AnnData` with the files its cell-major
`X` was laid out in, or the temporary `.h5mu` of a `MuData`.

More lifecycle details (ports, multiple viewers): {doc}`11_viewer_lifecycle_cleanup_ports_and_multiple_viewers`.

## Troubleshooting

Start here:
- {doc}`15_troubleshooting_viewing`

In a notebook, also run (in a new cell):

```python
viewer.debug_connection()
```
