# Jupyter embedding architecture

This page explains how Cellucid shows the web app inside notebooks and how the
bidirectional communication works. JupyterLab is tested end to end by the
project; Jupyter Notebook 7, VS Code, Colab, marimo and Databricks render
anywidgets with the same widget protocol and are expected to work, but are not
tested ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).

If you only need the user-facing quickstart, start with:
- {doc}`../d_viewing_apis/index`
- {doc}`../e_jupyter_hooks/index`

---

## The mental model (one sentence)

A notebook viewer is:

> **an anywidget whose output frames the Cellucid app, relaying the app's reads and its one viewer channel to the kernel over the notebook's own connection** (commands one way, events the other; no port, no server).

A browser-tab viewer (`AnnDataViewer`) is the same API over a local server: the
page and the data over HTTP, the channel over one WebSocket.

---

## Components and responsibilities

### Viewer classes

Implementation:
- `cellucid-python/src/cellucid/jupyter/`

Key classes:
- `AnnDataViewer` (`jupyter/_viewer.py`): the viewer API — hooks, state,
  commands, results (`pull`, `push`, …), session capture — and a browser-tab
  viewer that serves the data and the app via `AnnDataServer`.
- `Widget` (`jupyter/_widget.py`, front end `jupyter/static/widget.js`): a
  subclass of both `AnnDataViewer` and `anywidget.AnyWidget`. It prepares the
  same `AnnDataServer` but never binds it: requests are answered by
  `AnnDataServer.answer_request`, and the viewer channel runs over
  `_CommSocket`, a socket made of widget messages.

Convenience function:
- `show_anndata(...)` → `Widget`, displayed in the cell's output.

### Identity and routing

Each viewer has:
- `viewerId`: names the viewer when its page opens the channel
- `viewerToken`: the secret that page proves it holds, once, in the channel's
  first frame

A widget hands both to its page in a message; a browser tab carries them in
its address. A browser-tab viewer also runs its own server on its own port:
with `port=None` it binds `127.0.0.1:8765` (the one local port Google sign-in
accepts) or, when another viewer holds it, the next free port, and prints a
note saying so; a port given explicitly is used or refused.

---

## How the widget connects

1. The widget's front end creates a frame at
   `<web_source_url>/?viewer=embedded` (`https://www.cellucid.com` by default).
   The address carries no credentials.
2. The page posts `ready` to its parent. The widget accepts it only from its
   own frame and from the app's origin, and answers `connect` with a
   `MessagePort` and the viewer's id and token, posted to that exact origin.
3. Everything after travels on that port: the page's requests
   (`{t: 'fetch', method, path, range}`, answered with status, headers and
   bytes) and the viewer channel's frames. The widget relays each to the
   kernel as an anywidget message (binary buffers for bodies; requests that
   arrive together go as one message) and back. At most two batches of
   requests are with the kernel at once (one being answered, the next
   ready): the rest wait in the widget, where a request the page withdraws
   (`{t: 'cancel', id}`, sent when the read that asked it is cancelled, from
   the reader worker through the page) is dropped before the kernel reads
   for it. A request the kernel already holds is answered, and its answer
   dropped. A request the kernel leaves unanswered for long (a cell running)
   is waited for and said ("still waiting"), never asked again: a copy would
   only queue behind it. The front end is shared byte for byte with
   cellucid-r, so R's widget behaves the same.
4. In the kernel, `Widget._on_page_message` answers requests from the
   dataset and feeds channel frames to `serve_viewer_channel`, the same code a
   WebSocket runs.

The web side is `cellucid/assets/js/data/viewer-embed.js` (the transport) and
`formats/io/host-fetch.js` (requests to `https://host.cellucid.invalid/…`,
from the page and its reader workers). The full rules are in
{doc}`../e_jupyter_hooks/05_architecture_viewer_channel`.

### Between cells

The kernel handles a widget's messages on its own thread between cells. So
the viewer loads after the cell that showed it ends, hooks run while the
kernel is idle, and a call that would wait for the viewer on that thread
(`wait_for_ready`, `wait_for_event`, a session request's `result`) raises at once
(`Widget._require_answerable`) with what to do instead. Commands never wait:
`set_data_mapping` and `push` hold what they cannot send yet and send it on
the `data_mapping` or `ready` event. `request_session_bundle()` asks now; its
`result()` -- or `apply_session_to_anndata(adata, request)` --
collects in a later cell.

---

## How a browser tab's URL is constructed

The viewer URL is built from the **client-reachable server URL** plus query params:

- `viewer=true`
- `viewerId=<id>`
- `viewerToken=<token>`
- `anndata=true`

Key method:
- `AnnDataViewer._get_client_server_url()`: the exact `client_server_url=`
  argument, when supplied, or the bound data-server URL otherwise. It does not
  infer or probe a proxy (see {doc}`06_configuration_env_vars_and_logging`).

---

## The viewer channel

### Python → frontend: commands

When you call:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
```

Python validates the command (`_wire._require_frontend_message`) and sends it
on the channel (`cellucid._viewer_channel.send_viewer_command`): a command
that names cells (`highlight`, `setVisibility`, `setField`) as one binary
frame (the cells as little-endian `uint32`, strictly ascending), every other
command as one JSON text frame. Before the viewer has connected, the command
waits in the viewer's outbox and is sent, in order, when it connects.

Commands are used for:
- highlighting cells
- changing color-by
- visibility toggles
- reset camera
- pushing a field (`setField`) or a gene set (`setGeneSet`)
- “freeze” before shutdown
- requesting a session bundle
- changing the data mapping, pings and debug snapshots

### Frontend → Python: events (hooks)

The viewer sends events on the same channel:
- JSON text frames (hover, click, ready, pong, debug snapshot, command
  errors, data mapping, color-by, result removals and the results inventory)
- binary frames for selections (little-endian `uint32` cell indices, strictly
  ascending), results (cells, then values) and session bundles

The channel belongs to one viewer, so events name none; Python hands each to
that viewer, one at a time and in order, on the channel's thread.

In Python, `AnnDataViewer._handle_frontend_message`, in this order:
- validates the event (`_wire._require_inbound_jupyter_event`; a bad one
  closes the channel with `1008`),
- records the event into `viewer.state` (latest snapshot) and the results
  into `viewer.results`, hands it to the `wait_for_event(...)` calls waiting
  for it, and keeps the most recent events for `debug_connection()`,
- triggers hook callbacks registered via:
  - `@viewer.on_selection`, `@viewer.on_hover`, `@viewer.on_click`,
    `@viewer.on_ready`, `@viewer.on_command_error`, `@viewer.on_data_mapping`, and
    `@viewer.on_message` (`on_message` callbacks first).

The event is recorded before any hook runs. A hook that raises is logged
(`cellucid.viewer_channel`, "A viewer event hook failed"; several failures as
one `ExceptionGroup`) and the channel stays open; the other hooks still run,
and the `viewer.state` update and the `wait_for_event` hand-off have already
happened.

Protocol details are documented here:
{doc}`11_hooks_events_protocol_and_schema`.

---

## Results and session bundles (the round trip)

Every result the viewer holds (highlight groups, user-defined fields, scores,
gene sets, the visible and the opened cells) arrives as a `result` event when
it is made, and all of them whenever the viewer opens its dataset, so
`viewer.pull(adata)` in a later cell writes them without waiting
({doc}`../g_api_reference_coverage/api/round_trip`).

A whole session is captured on request:

1) Python sends a command on the channel:
   - `{type: "requestSessionBundle", requestId: "..."}`
2) The frontend serializes a `.cellucid-session` bundle in-memory.
3) The frontend sends the bytes back on the channel as one binary frame
   (header `{type: "session_bundle", requestId: "..."}`, then the bundle);
   the widget carries it in pieces.
4) Python:
   - checks a pending request is waiting for it,
   - streams the frame to a temp file and validates the bundle,
   - completes the request with a `session_bundle` event (includes path), or a
     `command_error` for an invalid bundle.
5) `request.result()` returns a `CellucidSessionBundle` handle.

This design keeps:
- memory bounded (streamed to disk),
- notebook UX simple (no manual file download/upload),
- and allows reuse of the same session format as the web app.

---

## Troubleshooting

### Symptom: “The notebook output is empty, or says the app did not load”

Likely causes (ordered):
1) the front end does not render widgets (`anywidget` missing where the
   notebook server runs),
2) the browser cannot load `web_source_url` (the widget says the page has not
   answered as a Cellucid viewer),
3) the app at `web_source_url` is too old to be embedded.

How to confirm:
- run `viewer.debug_connection()` in a new cell and read `connected_outputs`,
  `channel_connected` and `state`.

Fix:
- install `anywidget` in the notebook server's environment and reload,
- allow `https://www.cellucid.com`, or host a current copy and pass
  `web_source_url=`.

### Symptom: “Hooks worked once, then stopped”

Common causes:
- a cell is running (messages wait until the kernel is idle),
- the same viewer was shown in a second output or page, which took the
  channel over (the first copy froze and says the viewer is open elsewhere),
- the viewer was garbage-collected (lost references),
- you re-ran a cell and created a new viewer, but are watching the old output.

Fix:
- keep a reference to the viewer object,
- call `viewer.stop()` when done,
- create a new viewer after notebook reloads.

### Symptom: “A browser-tab viewer's page is blank”

Likely causes: the exact source generation could not be established, the
server is not reachable from the browser (remote kernel), or it died. Read
`server_health`, `anndata_store_probes`, `viewer_index_probe` and
`frontend_roundtrip` in `viewer.debug_connection()`; forward the port or pass
a proxy's base as `client_server_url=`; see
{doc}`09_server_mode_architecture_endpoints_and_security`.
