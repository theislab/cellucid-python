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
- capture a session bundle and apply it to an identity-matched AnnData, and
- stop the local viewer server cleanly.

From a fresh clone:

```bash
python -m pip install -e . jupyterlab
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/jupyter_pancreas.ipynb
```

```{figure} ../../../_static/screenshots/jupyter/pancreas-notebook-embed.png
:alt: JupyterLab showing the checked-in 2,531-cell Pancreas H5AD in the embedded Cellucid viewer, colored by its real clusters field, with the successful ready report below.
:width: 1440px

A real JupyterLab session after
`show_anndata(DATA_PATH, height=650, dataset_name=..., dataset_id=...)`
loading. The embedded viewer is rendering all 2,531 Pancreas cells, colored by
the real `clusters` observation field; the output below confirms readiness,
dimensional mode, and dataset size. The notebook remains above the viewer, so
Python analysis and browser interaction stay in one reproducible workflow.
```

## Quick troubleshooting if hooks do not fire

1. Wait for `viewer.wait_for_ready(...)` to return before sending commands.
2. Confirm a selection into a highlight group; moving a lasso without
   confirming it does not emit the selection event used by the example.
3. Run `viewer.debug_connection()` and check that `server_health.status` is
   `ok` and that the dataset identity probe matches.
4. For a remote kernel, expose the server through your notebook deployment or
   an SSH tunnel and pass its exact base URL as `client_server_url=`.

```{figure} ../../../_static/screenshots/jupyter/pancreas-debug-connection.png
:alt: JupyterLab showing the successful direct-AnnData Pancreas connection report, exact web build, and explicit viewer cleanup cell.
:width: 1440px

The same real notebook after its connection check: the local server is running,
health is `ok`, the Pancreas dataset identity is verified, the frontend
round-trip succeeded, and the exact web build is recorded. The final
`viewer.stop()` cell makes cleanup explicit.
```

```{note}
Screenshot provenance: for this documentation capture only, the notebook used
the web app from the same checked-out repository at
`http://127.0.0.1:4173`. The downloadable notebook deliberately keeps the
normal released-build resolution path while passing the exact checked-in data
path, dataset name, and stable dataset ID.
```

```{tip}
The fastest plumbing check is `viewer.debug_connection()`. It reports the
server, exact web generation, dataset identity, message bridge, and forwarded
frontend errors in one object.
```

```{toctree}
:maxdepth: 1
:hidden:

jupyter_pancreas
```
