# Supported environments matrix

This page answers: “Will the notebook viewer and its hooks work in *my* environment?”

Short answer: **yes, wherever the notebook shows widgets**, and the browser
can load the Cellucid app from `https://www.cellucid.com` (or the copy you
name with `web_source_url=`).

{func}`~cellucid.show_anndata` returns a {class}`~cellucid.jupyter.Widget`, an
[anywidget](https://anywidget.dev) shown in the cell's output. The widget
frames the Cellucid app and connects it to the kernel through the notebook's
own connection (a Jupyter comm): the viewer's channel (commands and events)
and every byte the viewer reads of the dataset travel on it
({doc}`05_architecture_viewer_channel`). No port is opened and no address of
the kernel's machine is involved, so a remote kernel, an HTTPS hub or a cloud
notebook needs nothing extra.

## At a glance

**Audience**
- Wet lab / beginner: find your environment row, then use the quickstart.
- Computational users: read “What the widget needs” and “While a cell runs”.
- Developers/IT: read the network and security notes.

**Start here**
- Quickstart sanity check: {doc}`02_quickstart_minimal_roundtrip`
- Security model: {doc}`12_security_model`

## Matrix

Legend:
- ✅ tested by the project, end to end, in a real browser
- ☑️ expected to work (it shows anywidgets), not tested by the project
- ❌ not a notebook (use a browser tab)

| Environment | Python (`show_anndata`, `cellucid.Widget`) | R (`cellucid_view()`) |
|---|---|---|
| JupyterLab 4 (local or remote, HTTP or HTTPS) | ✅ | ✅ (IRkernel) |
| Jupyter Notebook 7 | ☑️ (same front end as JupyterLab) | ☑️ |
| VS Code notebooks (local or remote kernel) | ☑️ | ☑️ |
| Google Colab | ☑️ | ☑️ (R names Colab's widget manager itself, as anywidget does for Python) |
| marimo | ☑️ (make the widget the cell's output) | — |
| Databricks notebooks | ☑️ | — |
| JupyterHub, Binder, cloud notebooks | ☑️ (no proxy or port setup) | ☑️ |
| A plain Python script | ❌ `show_anndata()` says so: use {class}`~cellucid.AnnDataViewer` (a browser tab), `serve_anndata()` or `cellucid serve` | ❌ use `cellucid_view()` in R at a console, or RStudio/Positron |

What “tested” means: the project's browser suite runs a real JupyterLab,
executes a cell, and checks that the viewer renders, that a lasso drawn in it
reaches the kernel as a selection, and that a command from the kernel reaches
the viewer — from Python and from R
(`cellucid/tests/browser/notebook-widget-jupyterlab.spec.mjs`). The app there
is loaded from another origin than JupyterLab's, as `https://www.cellucid.com`
is in practice. The other front ends render anywidgets with the same widget
protocol, so the same code runs there, but the project has not run them; if
one does not work for you, `viewer.debug_connection()` says how far it got.

## What the widget needs

1. **A notebook that shows widgets.** JupyterLab and Notebook 7 need the
   `anywidget` package in the environment that runs the *notebook server*
   (`pip install anywidget`; installing `cellucid` installs it in the kernel's
   environment, which is the same one in most setups). VS Code, Colab, marimo
   and Databricks bring their own widget support.
2. **The browser can load the app** from `web_source_url`
   (`https://www.cellucid.com` by default). The browser loads it, not the
   kernel: a compute node without internet is fine, a laptop without internet
   is not. An institution can host a copy of the app and pass its address as
   `web_source_url=`. Without the app the output says that it cannot load it.

Nothing else: no port, no `jupyter-server-proxy`, no `client_server_url`, no
`allowed_hosts`.

## While a cell runs

The kernel handles a widget's messages between cells, as for every notebook
widget. So:

- the viewer loads once the cell that showed it finishes, and while a later
  cell runs it waits for data it has not read yet (what it shows stays
  interactive);
- hooks (`@viewer.on_selection`, …) run, and `viewer.state` updates, whenever
  the kernel is idle — which is when you are interacting with the viewer;
- a call that waits for the viewer's answer (`wait_for_ready()`,
  `wait_for_event()`, a session request's `result()`) cannot be answered
  inside the cell that waits (commands such as `push()` and `set_data_mapping()` never
  wait), so it raises at once and says what to
  do instead: register a hook, read `viewer.state` in a later cell, or — for a
  session bundle — ask with `request = viewer.request_session_bundle()` and
  collect `request.result()` in the next cell
  ({doc}`10_session_bundles_request_session_bundle`).

Commands (`highlight_cells`, `set_color_by`, …) are sent at once.

## A browser tab instead

{class}`~cellucid.AnnDataViewer` is the same API in a browser tab of its own:
it serves the app and the data from this machine and the page talks back over
a WebSocket. Use it from a plain script, or through an SSH tunnel
({doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`). Its
waiting calls are answered while a cell runs, since its server runs on its own
threads.

## Edge cases

- **Showing a widget twice** (displaying `viewer` in another cell): the newest
  output takes over the viewer; the older one says it was opened again
  elsewhere.
- **Kernel restarts**: the output keeps its last view and says the notebook
  stopped it. Re-run the cell.
- **Reopening a saved notebook**: the output needs its kernel; re-run the cell.
- **Corporate networks** that block `https://www.cellucid.com`: the output says
  the app did not load; host a copy and pass `web_source_url=`.

## Troubleshooting

- Run (in a new cell):
  ```python
  viewer.debug_connection()
  ```
- Then follow: {doc}`14_troubleshooting_hooks`

## Next steps

- Minimal working loop: {doc}`02_quickstart_minimal_roundtrip`
- Architecture details: {doc}`05_architecture_viewer_channel`
