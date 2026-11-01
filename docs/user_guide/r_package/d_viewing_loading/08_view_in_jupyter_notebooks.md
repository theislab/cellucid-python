# View an Object in a Jupyter Notebook (R kernel)

**Audience:** anyone using R in JupyterLab, Jupyter Notebook 7, VS Code
notebooks or Google Colab (IRkernel)
**Time:** 2 minutes
**What you'll do:** show a Seurat or SingleCellExperiment object, or a file,
in the Cellucid viewer inside the notebook cell's output.

## The one line

```r
library(cellucid)
viewer <- cellucid_view(seu)
```

In an IRkernel session `cellucid_view()` shows the viewer in the cell's
output. It takes the same data as everywhere else: a Seurat,
SingleCellExperiment or SummarizedExperiment object, or the path of any file
the web app opens.

## What you need

- The R packages **IRkernel** and **IRdisplay** (`install.packages(c("IRkernel", "IRdisplay"))`).
- The Python package **anywidget** where the notebook *server* runs
  (`pip install anywidget`): the viewer is a notebook widget, the same one
  cellucid-python's `show_anndata()` shows, and the notebook's widget support
  draws it. VS Code and Colab bring their own.
- A browser that can load `https://www.cellucid.com` (the app itself; your
  data does not go there). Name a copy you host with
  `web_source_url = "https://…"`.

No port, proxy or server address is involved: the viewer reaches this R
session through the notebook's own connection, so the same cell works on a
laptop, on a remote JupyterHub or Binder, and over HTTPS.

## What happens

1. The object is read into the same Cellucid-ready AnnData
   `cellucid_prepare()` writes, held in memory; nothing is written to disk.
2. The cell's output frames the Cellucid app. The app asks the output for the
   viewer's credentials (they never appear in an address), then reads the
   dataset and opens the viewer channel through it; R answers every request
   from the object, exactly as its server would.
3. Commands (`cellucid_highlight()`, `cellucid_color_by()`, …) go to the
   viewer; events (selections, hovers, …) come back to your hooks
   ({doc}`06_hooks_selection_and_sessions`).

## Between cells

The R kernel hands the viewer's messages to R between cells, on R's one
thread. So:

- the viewer loads once the cell that showed it finishes, and while a later
  cell runs it waits for data it has not read yet;
- hooks registered with `cellucid_on()` run, and `cellucid_selection()` is
  current, whenever the kernel is idle — which is when you are using the
  viewer;
- a call that waits for the viewer (`cellucid_wait_ready()`, or
  `cellucid_session_result()` before the viewer answered) cannot be answered
  in the cell that waits: it stops at
  once and says to run it again in a later cell, or to register a hook. A
  session is asked for with `cellucid_request_session()` and collected
  (`cellucid_session_result()`) or applied (`cellucid_apply_session(object,
  request)`) in a later cell. Commands never wait: `cellucid_push()` and
  `cellucid_set_data_mapping()` go as soon as the viewer can take them.

```r
cellucid_on(viewer, "selection", function(event) {
  message(length(event$cells), " cells selected")
})
```

```r
# a later cell, after selecting in the viewer
sel <- cellucid_selection(viewer)
seu$picked <- colnames(seu) %in% sel$cell_names
```

## A cell run again

A cell run again replaces its output, so the viewer it showed before is
stopped: the object's matrices it held, its memoised arrays and its chunk
cache are let go (an object reassigned after `NormalizeData()` would
otherwise stay in memory beside the new one), as cellucid-python stops the
widget a re-run cell replaced. A viewer that another cell shows too stays.
The cell is known by the id JupyterLab, Notebook 7 and VS Code send with each
run; where none is sent, viewing a dataset again stops the earlier viewer of
the same dataset id, and says so.

## Why not a server here

Outside a notebook, `cellucid_view()` and `cellucid_serve()` run a small web
server in the background ({doc}`04_view_in_rstudio_positron`). The R kernel
never runs R's event loop between cells, so such a server would answer
nothing until a cell runs: in a notebook `cellucid_serve()` and
`pane = "browser"` stop with that explanation instead of opening a page that
never loads.

## Where it works

| Notebook | Shows the viewer |
|---|---|
| JupyterLab 4 | yes (tested by the project, in a real browser) |
| Jupyter Notebook 7, VS Code | expected (same widget protocol), not tested by the project |
| Google Colab (R runtime) | expected, not tested by the project (see below) |

### Google Colab

Colab shows a third-party widget only when the widget's display names
Colab's widget manager; for Python, `google.colab.output.enable_custom_widget_manager()`
(which anywidget calls) makes the kernel add that to every widget display.
In Colab's R runtime `cellucid_view()` adds the same metadata itself
(detected from Colab's `COLAB_RELEASE_TAG` environment variable), and the
manager loads the widget from jsDelivr and talks to R over the kernel's
comms, exactly as for a Python widget. The project has checked the messages
R sends, not Colab itself. If no viewer appears: write the data with
`cellucid_prepare(data = seu, out_path = "seu.zarr.zip", format = "zarr.zip")`,
download it, and open it at www.cellucid.com (**Open… › File…**), or use
Python in a Colab notebook (`cellucid.show_anndata("seu.zarr.zip")`).

The full matrix, for Python too: {doc}`../../python_package/e_jupyter_hooks/03_supported_environments_matrix`.
