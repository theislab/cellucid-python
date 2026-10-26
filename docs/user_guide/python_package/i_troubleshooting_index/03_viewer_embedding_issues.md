# Viewer embedding issues

## The notebook cell is empty

Keep the returned viewer object alive, confirm the cell produced HTML output,
and inspect the browser console for the first error. Then verify that the
viewer URL shown by the Python object is reachable from the browser that renders
the notebook.

## The viewer frame cannot reach the local server

Remote notebooks, JupyterHub, and HTTPS notebook origins need a route the
browser can reach, such as `jupyter-server-proxy` or an SSH tunnel; pass its
exact base as `client_server_url=` (and a proxy's host as `allowed_hosts=`).
The browser—not the Python kernel—must be able to reach the displayed URL, and
any proxy in between must pass WebSockets for commands and hooks
({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).

## The UI loads but the dataset does not

Inspect the first failing dataset request and its response body. A viewer-asset
problem and a dataset-server problem are separate paths and should be diagnosed
separately.

See {doc}`../d_viewing_apis/06_jupyter_show_anndata_quickstart`,
{doc}`../d_viewing_apis/10_notebook_widget_mode_advanced`, and
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.
