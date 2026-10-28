# Notebook / widget mode (advanced)

This page explains the notebook widget behind `show_anndata(...)`, including:
- what the widget is, and why it needs no port or proxy,
- where the app comes from, and what happens offline,
- what happens while a cell runs, and the calls that cannot wait in one,
- how to choose the neighbor graph, which the viewer reads only when a view
  draws it,
- and how to debug the browser ↔ kernel connection.

If you want the minimal notebook quickstart, go to {doc}`06_jupyter_show_anndata_quickstart`.

## At a glance

**Audience**
- JupyterHub, Colab, VS Code and cloud notebook users
- anyone whose notebook output shows no viewer

## The widget

{func}`~cellucid.show_anndata` creates a {class}`~cellucid.jupyter.Widget` and
shows it; `cellucid.Widget(...)` creates one without showing it (end a cell
with it to show it). A widget is an [anywidget](https://anywidget.dev): its
front end runs in the cell's output, and talks to the kernel over the
notebook's own connection (a Jupyter comm), which every notebook already has
between the browser and the kernel, wherever the kernel runs.

The output frames the Cellucid app at `<web_source_url>/?viewer=embedded`.
The app then:

1. tells the output it is ready; the output accepts that only from its own
   frame and the app's origin, and answers with the viewer's id and token and
   a private `MessagePort`, posted to that origin only (the token never
   appears in an address, a widget trait a notebook may save, or a log);
2. reads the dataset's documents and bytes (by `Range`, like any remote file)
   through the port; the output relays them to the kernel, which answers from
   the dataset ({meth}`~cellucid.AnnDataServer.answer_request`);
3. opens the viewer channel through the port: the same `hello`/`welcome`,
   commands and events as over a WebSocket
   ({doc}`../e_jupyter_hooks/05_architecture_viewer_channel`).

Requests that arrive together go to the kernel in one message, and bodies
travel as binary buffers.

## Where the app comes from

The browser loads the app from `web_source_url`, `https://www.cellucid.com`
by default: the build the site publishes, as every visitor gets it. The
protocol between the app and the kernel is negotiated when the viewer opens
(`_cellucid/protocol`), so an app newer than the installed package uses only
what the package understands. An institution can host a copy of the app and
pass its address:

```python
viewer = show_anndata(adata, dataset_name="Example", dataset_id="example",
                      web_source_url="https://cellucid.example.org")
```

Offline: the *browser* needs to reach `web_source_url`; the kernel does not
(a compute node without internet is fine). A browser that cannot load the app
shows a message in the output after the page failed to answer. Without any
internet on the machine you browse from, host a copy of the app on your
network, or open the viewer in a browser tab with
{class}`~cellucid.AnnDataViewer`, which serves a verified copy kept on disk
({doc}`09_server_mode_advanced`).

## While a cell runs

The kernel handles a widget's messages between cells, on the thread that
runs cells. Hence:

- the viewer loads once the cell that showed it finishes; a long cell after
  `show_anndata(...)` in the same cell delays it;
- while any later cell runs, the viewer waits for bytes it has not read yet
  (a new gene, a new embedding); what it shows stays interactive;
- hooks run, and `viewer.state` updates, whenever the kernel is idle;
- `wait_for_ready()`, `wait_for_event()` and a session request's `result()`
  (and `apply_session_to_anndata(adata, request)`) wait for an answer that
  cannot arrive while their own cell runs: they raise at once with what to do
  instead (`set_data_mapping()` and `push()` never wait: they go as soon as
  the viewer can take them). A thread of your own may wait (the viewer answers it as soon
  as the cell ends). A session bundle is asked in one cell and collected in
  the next:

```python
request = viewer.request_session_bundle()
```

```python
bundle = request.result()
adata = bundle.apply_to_anndata(adata, expected_dataset_id="example")
```

Commands (`highlight_cells`, `set_color_by`, …) are sent at once.

## A notebook viewer starts quietly

`show_anndata(...)` prepares the dataset with output suppressed, so the
numbered report a terminal server prints never appears. Nothing reads an
array there: the dataset's metadata is read and the launch choices you passed
are checked; the viewer then reads what each view needs. Writing a temporary
`.h5ad` (`.h5mu`) for an in-memory `AnnData` (`MuData`), once, is what can
take a while. To see the report, run the same object through a terminal
server once — `serve_anndata(..., quiet=False)` prints the four steps.

## The neighbor graph

The viewer reads the neighbor graph only when a view draws it. Without
`knn_key` it uses `obsp['connectivities']` when the store has it; `knn_key`
names another `obsp` matrix. The parameter is keyword-only and accepted by
{func}`~cellucid.show_anndata`, {class}`~cellucid.jupyter.Widget` and
{class}`~cellucid.AnnDataViewer`:

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

```{important}
Naming a matrix the object's `obsp` does not hold raises during construction
rather than producing a graph-less viewer. Compute it first with
`sc.pp.neighbors(adata)`, or leave the parameter unset.
```

The rules are in {doc}`08_anndata_mode_show_anndata_and_serve_anndata`.

## Debugging notebook connectivity (high-signal)

In a new cell, run:

```python
viewer.debug_connection()
```

It waits for nothing, and reports:
- `web_source_url`: where the app is loaded from
- `connected_outputs`: the outputs whose page has a channel open
- `channel_connected`: whether the viewer channel is open
- `state`: whether the viewer said it is ready, and its last event
- `recent_events`: accepted events by type, and the commands the viewer
  refused (`command_errors`)

## Troubleshooting

- The output shows text instead of a viewer (“This output is a Cellucid
  viewer (a notebook widget)…”) → the notebook cannot show widgets: install
  `anywidget` where the notebook *server* runs (`pip install anywidget`) and
  reload the page.
- The output says the page “has not answered as a Cellucid viewer” → the
  browser cannot load `web_source_url` (offline, a blocked network), or it is
  not a Cellucid build that can be embedded.
- The output says it is waiting for the notebook kernel → a cell is still
  running; the viewer connects when it finishes.
- A waiting call raised “… hears the viewer only between cells” → run it in a
  later cell, use a hook, or (session bundles) `request_session_bundle()`.
- The edge overlay is missing although `adata.obsp['connectivities']` exists →
  the graph is read only when a view draws edges; turn the edges on in the
  viewer, or choose the graph under **Session › Data mapping**.
- Construction raised “knn_key names obsp['connectivities'],
  which this store does not have.” → run `sc.pp.neighbors(adata)` first, or
  construct without the parameter.
- Construction raised “Cellucid found nothing to plot” → no `obsm` array can be
  plotted. The message lists every `obsm` entry and why it was not used; name a
  usable one with `embedding_keys={"2d": "<key>"}`.
  {doc}`08_anndata_mode_show_anndata_and_serve_anndata` has the complete rules.
