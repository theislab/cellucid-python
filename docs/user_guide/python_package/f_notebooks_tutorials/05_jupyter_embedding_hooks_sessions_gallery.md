# Jupyter embedding, hooks, and sessions

This runnable notebook connects the embedded viewer, Python callbacks, and
session round-tripping with the real Pancreas AnnData checked into the
repository.

If you want the conceptual docs (schemas, environment matrix, troubleshooting), start with:
- {doc}`../e_jupyter_hooks/index`

If you want the “hands-on but explained” tutorials, start with:
- {doc}`23_programmatic_highlighting_and_selection_callbacks`
- {doc}`32_session_persistence_and_restoring_analysis_artifacts`

## Runnable notebook

{doc}`jupyter_pancreas` is a single, short lifecycle:

- locate and inspect the checked-in 2,531-cell Pancreas H5AD,
- open it with the identity-explicit
  `show_anndata(DATA_PATH, ..., dataset_name=..., dataset_id=...)` path,
- color by the real `clusters` observation field,
- highlight cells from Python,
- summarize confirmed viewer selections in Python,
- pull the highlights and fields made in the viewer onto an in-memory
  AnnData with `viewer.pull()`,
- ask for a session bundle in one cell, collect it in the next and apply it to
  an identity-matched AnnData, and
- stop the viewer cleanly.

From a fresh clone:

```bash
python -m pip install -e . jupyterlab
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/jupyter_pancreas.ipynb
```

```{figure} ../../../_static/screenshots/jupyter/show-anndata-driven-from-python.png
:alt: The Cellucid viewer after Python called viewer.set_color_by("clusters") and viewer.highlight_cells(): the Pancreas sample coloured by clusters with the highlighted cells.
:width: 1440px

The viewer the notebook drives, captured in a browser tab: the same app the
notebook's output frames, after Python ran `viewer.set_color_by("clusters")`
and `viewer.highlight_cells(...)`.
```

## Quick troubleshooting if hooks do not fire

1. Run the cells one at a time: the viewer appears when the cell that showed
   it finishes, and a notebook hands the viewer's messages to the kernel
   between cells, so hooks run while the kernel is idle.
2. Confirm a selection into a highlight group; moving a lasso without
   confirming it does not emit the selection event used by the example.
3. Run `viewer.debug_connection()` in a new cell and check that
   `connected_outputs` names the output and `channel_connected` is `True`.
4. No output at all in JupyterLab or Notebook 7: install `anywidget` where the
   notebook server runs and reload the page
   ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).

A healthy report, as the notebook prints it:

```text
Transport: notebook widget
App loaded from: https://www.cellucid.com
Connected outputs: ['<output id>']
Channel connected: True
Dataset: {'id': 'pancreas-endocrinogenesis-day15-5', 'format': 'h5ad'}
Viewer ready: {'n_cells': 2531, 'n_obs': 2531, 'dimensions': 2}
```

The notebook's output is connected, the viewer's channel is open, the viewer
reads the checked-in `.h5ad` in place, and it opened all 2,531 cells. A final
`viewer.stop()` makes cleanup explicit.

```{note}
The notebook loads the app from www.cellucid.com (the default
`web_source_url`) while passing the exact checked-in data path, dataset name,
and stable dataset ID.
```

```{tip}
The fastest plumbing check is `viewer.debug_connection()`. It waits for
nothing, so it works in any cell, and reports where the app comes from, the
connected outputs, the channel, what the viewer last sent and the commands it
refused.
```

More notebooks, each with an R twin: {doc}`06_cellucid_in_your_pipeline_python_and_r`.

```{toctree}
:maxdepth: 1
:hidden:

jupyter_pancreas
```
