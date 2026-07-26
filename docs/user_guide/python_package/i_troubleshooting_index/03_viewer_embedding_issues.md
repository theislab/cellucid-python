# Viewer embedding issues

## The notebook cell is empty

Keep the returned viewer object alive, confirm the cell produced HTML output,
and inspect the browser console for the first error. Then verify that the
viewer URL shown by the Python object is reachable from the browser that renders
the notebook.

## The viewer frame cannot reach the local server

Remote notebooks, JupyterHub, and HTTPS notebook origins may need
`jupyter-server-proxy`. The browser—not the Python kernel—must be able to reach
the displayed URL.

## The UI loads but the dataset does not

Inspect the first failing dataset request and its response body. A viewer-asset
problem and a dataset-server problem are separate paths and should be diagnosed
separately.

See {doc}`../d_viewing_apis/06_jupyter_show_and_show_anndata_quickstart`,
{doc}`../d_viewing_apis/10_notebook_widget_mode_advanced`, and
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.
