# Jupyter embedding architecture

This page explains how Cellucid embeds the web app inside Jupyter notebooks and how the bidirectional communication works. Classic Jupyter and JupyterLab are tested; VS Code notebooks and Colab have not been tested by the project ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).

If you only need the user-facing quickstart, start with:
- {doc}`../d_viewing_apis/index`
- {doc}`../e_jupyter_hooks/index`

---

## The mental model (one sentence)

A notebook viewer is:

> **a local Python server + an iframe pointed at that server + one WebSocket from the iframe to that server** (the viewer channel: commands one way, events the other).

---

## Components and responsibilities

### Viewer classes

Implementation:
- `cellucid-python/src/cellucid/jupyter/`

Key class:
- `AnnDataViewer` (`jupyter/_viewer.py`): serves the data via `AnnDataServer`
  (any format Cellucid reads, or an AnnData or MuData object), and owns the
  notebook embedding and hook/event handling

Convenience function:
- `show_anndata(...)` → `AnnDataViewer`

### Identity and routing

Each viewer has:
- `viewerId`: names the viewer when its page opens the channel
- `viewerToken`: the secret that page proves it holds, once, in the channel's
  first frame

Each viewer also runs its own server on its own port: with `port=None` it
binds `127.0.0.1:8765` (the one local port Google sign-in accepts) or, when
another viewer holds it, the next free port, and prints a note saying so; a
port given explicitly is used or refused.

---

## How the iframe URL is constructed

The viewer URL is built from the **client-reachable server URL** plus query params:

- `viewer=true`
- `viewerId=<id>`
- `viewerToken=<token>`
- `anndata=true`

Key method:
- `AnnDataViewer._get_client_server_url()`

### Why this is complicated

Notebook environments differ:

- Local classic/JupyterLab: browser can usually reach `http://127.0.0.1:<port>`.
- Remote kernels (JupyterHub, SSH to server, etc.): browser cannot reach kernel loopback.
- HTTPS notebooks: browser blocks HTTP loopback as mixed content.
- Colab exposes a special HTTPS port proxy
  (`google.colab.kernel.proxyPort`), which the caller must invoke and pass
  explicitly.

Cellucid selects one URL:

- the exact `client_server_url=` argument, when supplied, or
- the bound data-server URL otherwise.

It does not infer or probe a notebook proxy. Remote/HTTPS callers must expose
the port and pass the exact browser-reachable HTTP(S) base (see
{doc}`06_configuration_env_vars_and_logging`).

---

## The viewer channel

The page opens one WebSocket to `<server URL>/_cellucid/ws` (the same address
and path prefix it was loaded from, so a notebook proxy that carries the page
carries the channel too, provided it passes WebSockets) and says
`{"type": "hello", "viewerId", "viewerToken"}`; the server answers
`{"type": "welcome"}`. The full rules are in
{doc}`../e_jupyter_hooks/05_architecture_viewer_channel`.

### Python → frontend: commands

When you call:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
```

Python validates the command (`_wire._require_frontend_message`) and sends it
on the channel (`cellucid._viewer_channel.send_viewer_command`): a highlight
or a visibility change, which name cells, as one binary frame (the cells as
little-endian `uint32`, strictly ascending), every other command as one JSON
text frame. Before the viewer has
connected, the command waits in the viewer's outbox and is sent, in order,
when it connects. The notebook output holds no script, only the iframe.

Commands are used for:
- highlighting cells
- changing color-by
- visibility toggles
- reset camera
- “freeze” before shutdown
- requesting a session bundle
- changing the data mapping, pings and debug snapshots

### Frontend → Python: events (hooks)

The viewer sends events on the same channel:
- JSON text frames (hover, click, ready, pong, debug snapshot, command
  errors, data mapping)
- binary frames for selections (little-endian `uint32` cell indices, strictly
  ascending) and session bundles

The channel belongs to one viewer, so events name none; the server hands each
to that viewer, one at a time and in order, on the channel's thread.

In Python, `AnnDataViewer._handle_frontend_message`, in this order:
- validates the event (`_wire._require_inbound_jupyter_event`; a bad one
  closes the channel with `1008`),
- records the event into `viewer.state` (latest snapshot), hands it to the
  `wait_for_event(...)` calls waiting for it, and keeps the most recent events
  for `debug_connection()`,
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

## Session bundle capture (“no download” notebook workflow)

Notebook users often want to:
- interactively create highlights/annotations in the UI,
- then pull that state back into Python for downstream analysis.

Flow:

1) Python sends a command on the channel:
   - `{type: "requestSessionBundle", requestId: "..."}`
2) The frontend serializes a `.cellucid-session` bundle in-memory.
3) The frontend sends the bytes back on the channel as one binary frame
   (header `{type: "session_bundle", requestId: "..."}`, then the bundle).
4) The server:
   - checks a pending request is waiting for it,
   - streams the frame to a temp file and validates the bundle,
   - routes a `session_bundle` event back to Python (includes path), or a
     `command_error` for an invalid bundle.
5) Python returns a `CellucidSessionBundle` handle to the caller.

This design keeps:
- memory bounded (streamed to disk),
- notebook UX simple (no manual file download/upload),
- and allows reuse of the same session format as the web app.

---

## When things go wrong: browser cannot reach the server

If the notebook is served from HTTPS or the kernel is remote, direct loopback
may produce a blank iframe, a mixed-content error, or a failed network request.
Cellucid does not render a special proxy-selection page. When the page loads
but its channel cannot open, the viewer says so and names a proxy that does not
pass WebSockets as the cause.

---

## Troubleshooting

### Symptom: “The viewer iframe is blank”

Likely causes (ordered):
1) the exact source generation could not be established,
2) mixed-content/proxy URL mismatch,
3) the server is not reachable from the browser (remote kernel),
4) the server died (port conflict, exception).

How to confirm:
- run `viewer.debug_connection()` and read:
  - `server_health`, `anndata_store_probes`, `viewer_index_probe`,
    `channel_connected`, `frontend_roundtrip`.

Fix:
- expose the port through a proxy/tunnel and pass its base as
  `client_server_url=`,
- ensure the exact UI source generation can be established (see
  {doc}`09_server_mode_architecture_endpoints_and_security`),
- check server logs for exceptions.

### Symptom: “Hooks worked once, then stopped”

Common causes:
- the notebook output cell was cleared/collapsed and the iframe was destroyed
  (commands then wait until a page of the viewer connects again),
- the same viewer was shown in a second output or page, which took the
  channel over (the first copy froze and says the viewer is open elsewhere),
- the viewer was garbage-collected (lost references),
- you re-ran a cell and created a new viewer/server, but are watching the old iframe.

Fix:
- keep a reference to the viewer object,
- call `viewer.stop()` when done,
- create a new viewer after notebook reloads.
