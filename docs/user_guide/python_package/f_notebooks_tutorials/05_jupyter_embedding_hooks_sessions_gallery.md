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

```{figure} ../../../_static/screenshots/jupyter/show-anndata-driven-from-python.png
:alt: The viewer cellucid.show_anndata() serves, opened at its viewer URL, after Python called viewer.set_color_by("clusters") and viewer.highlight_cells() on the Ductal cells: the Pancreas sample coloured by clusters with the Ductal cells highlighted.
:width: 1440px

What the notebook's output frame holds: the viewer `show_anndata()` serves for
the Pancreas sample, here opened at its `viewer.viewer_url` in a browser tab,
after the notebook ran `viewer.set_color_by("clusters")` and
`viewer.highlight_cells(...)` on the Ductal cells. The Dataset box reads
*Zarr v2 folder · Python*: the data comes from the kernel's server, not from
a website.
```

## Quick troubleshooting if hooks do not fire

1. Wait for `viewer.wait_for_ready(...)` to return before sending commands.
2. Confirm a selection into a highlight group; moving a lasso without
   confirming it does not emit the selection event used by the example.
3. Run `viewer.debug_connection()` and check that `server_health.status` is
   `ok` and that the served store answers its probe (`anndata_store_probes`:
   a `206` range read).
4. For a remote kernel, expose the server through your notebook deployment or
   an SSH tunnel and pass its exact base URL as `client_server_url=`.

A healthy report, abridged (the real one also carries the viewer and server
URLs with the viewer's token, the web cache's absolute path, the browser's
user agent and timestamps; `<build>` and `<version>` are your web build and
package version):

```text
{'server_running': True,
 'channel_connected': True,
 'web_ui': {'cache': {'index_html_build_id': '<build>',
                      'verified_build_id': '<build>', ...}, ...},
 'state': {'ready': {'n_cells': 3696, 'dimensions': 3}, ...},
 'server_health': {'status': 'ok', 'type': 'anndata', 'version': '<version>',
                   'format': 'zarr', 'n_cells': 3696, 'n_genes': 3753},
 'anndata_store_probes': {'pancreas': {'head': {'status': 200, ...},
                                       'range': {'status': 206,
                                                 'content_range': 'bytes 0-0/19089'}}},
 'recent_events': {'count': 2, 'by_type': {'data_mapping': 1, 'ready': 1},
                   'command_errors': []},
 'frontend_roundtrip': {'ok': True, 'pong': {...}},
 'frontend_debug_snapshot': {'ok': True, 'snapshot': {'connected': True, ...}},
 ...}
```

The server is running, `server_health.status` is `ok` and names the served
dataset's size, the store answers a one-byte range read with `206`, the
frontend round-trip succeeded, and `web_ui` records the exact web build the
viewer was checked against. A final `viewer.stop()` makes cleanup explicit.

```{note}
The downloadable notebook keeps the normal released-build resolution path (the
web app comes from www.cellucid.com, cached) while passing the exact
checked-in data path, dataset name, and stable dataset ID.
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
