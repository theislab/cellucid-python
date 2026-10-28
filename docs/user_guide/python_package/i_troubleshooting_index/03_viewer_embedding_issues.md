# Viewer embedding issues

## The notebook cell is empty

`show_anndata()` shows a notebook widget. An empty output means the notebook
front end does not show it: install `anywidget` in the environment that runs
the *notebook server* (installing `cellucid` installs it in the kernel's
environment, which some setups keep separate) and reload the page, or use a
front end that shows widgets — JupyterLab, Notebook 7, VS Code, Colab, marimo
or Databricks ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).
Outside a notebook `show_anndata()` raises and names the browser-tab
alternatives.

## The output says the app did not load

The widget frames the Cellucid app, which the **browser** loads from
`https://www.cellucid.com`. If the browser cannot reach it (offline laptop,
blocked address), allow it, or host a copy of the app and pass its address as
`web_source_url=`. The kernel needs no network access and no port.

## The viewer shows the data, but Python hears nothing

The kernel handles the widget's messages between cells: hooks fire and
`viewer.state` updates once the running cell ends. A call that waits for the
viewer in the same cell raises at once and says what to do instead. Run
`viewer.debug_connection()` in a new cell and check `channel_connected`.

## A browser-tab viewer cannot reach its server

`AnnDataViewer`, `serve_anndata` and `cellucid serve` run a server on the
kernel's machine. From a remote kernel, forward the port over SSH, or pass a
reverse proxy's exact base as `client_server_url=` and its host as
`allowed_hosts=`; the proxy must pass WebSockets for commands and hooks
({doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`).

## The UI loads but the dataset does not

Inspect the first failing dataset request and its response body. A viewer-asset
problem and a dataset problem are separate paths and should be diagnosed
separately.

See {doc}`../d_viewing_apis/06_jupyter_show_anndata_quickstart`,
{doc}`../d_viewing_apis/10_notebook_widget_mode_advanced`, and
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.
