# Architecture: the viewer channel

This page explains how a notebook viewer and Python talk to each other. If you
only want something that works, start with {doc}`02_quickstart_minimal_roundtrip`.

## Design goals

1. One protocol for both directions, over two transports: a notebook widget's
   comm (what {func}`~cellucid.show_anndata` and
   {class}`~cellucid.jupyter.Widget` use) and a WebSocket to a local server (a
   browser-tab {class}`~cellucid.AnnDataViewer`, the R package outside a
   notebook). Everything above the transport is one implementation.
2. No port and no proxy in a notebook: the widget reaches the kernel through
   the notebook's own connection, so a remote kernel, an HTTPS hub, Colab,
   VS Code and Databricks need no setup.
3. Data and events stay between the browser and the kernel: never a cloud
   relay. The app itself is loaded from `https://www.cellucid.com` (or the
   copy named by `web_source_url`); the dataset never goes there.
4. Credentials never appear in an address a proxy could log: the widget hands
   them to the page in a message.
5. Failures are diagnosable (`viewer.debug_connection()`).

## System diagram: the notebook widget

```text
┌──────────────────────────── Notebook kernel (Python or R) ───────────────────┐
│  cellucid.Widget (anywidget) / cellucid_view() in IRkernel                    │
│    "fetch" {method, path, Range}  → the dataset's documents and bytes         │
│       (AnnDataServer.answer_request / the R server's routes, no socket)       │
│    "open"/"text"/"binary"/"close" → the viewer channel, the same code a       │
│       WebSocket runs (serve_viewer_channel / .on_ws_open)                     │
└───────────────────────────────────────────▲──────────────────────────────────┘
                                            │ the notebook's connection (comm
                                            │ messages, binary buffers; the
                                            │ page's frames to IRkernel as
                                            │ base64), between cells
┌───────────────────────────────────────────▼──────────────────────────────────┐
│  Widget front end (widget.js, in the cell's output)                          │
│    <iframe src="https://www.cellucid.com/?viewer=embedded">                  │
│    accepts "ready" only from that frame and origin, answers "connect"         │
│    {viewerId, viewerToken} + a MessagePort, posted to that origin only        │
└───────────────────────────────────────────▲──────────────────────────────────┘
                                            │ MessagePort (no other window can post to it)
┌───────────────────────────────────────────▼──────────────────────────────────┐
│  The Cellucid web app (ViewerHostDataSource, data/viewer-embed.js)            │
│    viewer channel over a WebSocket-shaped PortSocket                          │
│    requests to https://host.cellucid.invalid/… answered over the port, in     │
│    the page and in its reader workers (formats/io/host-fetch.js)              │
└──────────────────────────────────────────────────────────────────────────────┘
```

The widget's messages name the output they belong to, so one widget shown in
two outputs is two connections; as with a reloaded page, the newer one takes
the viewer over and the older one is told so (4409). Requests that arrive
together go to the kernel as one message, at once (the page's reader bounds
how many it has under way, as with a web host); the Python kernel answers
each on a thread of its own, as its HTTP server does, so a request that takes
long to make (a chunk of a matrix the layout has not reached yet, made from
its progress and the object) holds back no other, and R answers them in turn. An R kernel has at
most two such messages at once (the widget's `_fetch_messages` trait): the
requests asked while they are answered gather into the next message, since
IRkernel spends about 1.7 ms of its own on every message it reads. Bodies travel as
binary buffers, from Python and R alike (IRkernel's own send writes none, so
the R package makes those messages with the kernel's reply and wire
functions);
the page's binary frames to an R kernel travel as base64, because IRkernel's
reader drops the buffers of a message it receives. The kernel handles them
between cells: see {doc}`03_supported_environments_matrix` for what that
means while a cell runs.

One limit of the notebook server shapes the answers: Jupyter Server stops
forwarding a kernel's messages to the browser (it drops them, it does not
delay them) while they exceed 1,000 a second over a 3-second window
(`ServerApp.iopub_msg_rate_limit`, `rate_limit_window`; comm messages
count, their bytes do not). So answers ready together share a message (at
most 8 MiB of bodies each), and the Python kernel sends at most one message
of answers every 2 ms, half that rate; R sends one message (or a few, by
size) for each message of requests.

## System diagram: a browser tab

```text
┌──────────────────────────── Python or R ─────────────────────────────────────┐
│  AnnDataViewer ── AnnDataServer (one per viewer)                              │
│     GET /                         the viewer's web build                      │
│     GET /_cellucid/datasets       the dataset and its launch choices          │
│     GET /_cellucid/anndata/...    the data's bytes (Range, sendfile)          │
│     GET /_cellucid/protocol       where the channel is, what it carries       │
│     GET /_cellucid/ws  ──────────  the viewer channel (one WebSocket)        │
└───────────────────────────────────────────┬──────────────────────────────────┘
                                            │  ws:// (wss:// behind HTTPS)
┌───────────────────────────────────────────▼──────────────────────────────────┐
│  Browser tab: <server>/?viewer=true&viewerId=…&viewerToken=…&anndata=true     │
│  The Cellucid web app (ViewerHostDataSource) opens the channel:            │
│    → hello {viewerId, viewerToken}      ← welcome                             │
│    ← commands (text and binary frames) → events (text and binary frames)     │
└──────────────────────────────────────────────────────────────────────────────┘
```

The page reaches the server at the address it was loaded from, path prefix
included, so a proxy that carries the page carries the channel too. The proxy
must pass WebSockets.

## Opening the channel

The steps are the same over the widget; there the protocol document is a
`fetch` answered by the kernel, the socket is a port socket, and the `Host`
and `Origin` checks of step 2 are the widget's own: it accepts the page's
`ready` only from its frame and the app's origin, and posts the credentials to
that origin only.

1. The viewer reads `GET /_cellucid/protocol`:
   `{"channel": "_cellucid/ws", "commands": [...], "binaryCommands": [...], "events": [...]}`.
   A server without `channel`, or without `binaryCommands` (it sends the
   cells of a highlight as JSON text), is an older `cellucid`; the viewer says
   to upgrade the package.
2. It opens a WebSocket at `<serverUrl>/_cellucid/ws`. The server checks the
   `Host` header (421 otherwise, as for every route) and the `Origin` header:
   the page must be this server's own host, a loopback name, or a name in
   `allowed_hosts` (403 otherwise).
3. The first frame is `{"type": "hello", "viewerId": …, "viewerToken": …}`.
   The server compares the token in constant time and answers
   `{"type": "welcome"}`, or closes with 4401.
4. From then on the socket is that viewer's; no later frame carries
   credentials.

## Commands (Python → viewer)

- One JSON text frame each, with exactly the keys of its type
  (`_wire._COMMAND_MESSAGE_FIELDS`); see {doc}`06_python_to_frontend_commands`.
- A command that names cells (`highlight`, `setVisibility`) is one binary
  frame: its header (`_wire._BINARY_COMMAND_HEADER_FIELDS`), padded to a
  multiple of 4 bytes, then the cells as little-endian `uint32`, strictly
  ascending, which the viewer reads in place. Two million cells are 8 MB,
  where the same cells as JSON text were 14 MB.
- A command sent before the viewer is connected (the cell output still
  loading, the page reconnecting) waits in the viewer's outbox and is sent,
  in order, the moment it connects. Nothing is sent twice.
- `viewer.send_message(...)` returns `True` when the command was sent now,
  `False` when it is waiting, and raises after `viewer.stop()`.

## Events (viewer → Python)

- Text frames for `hover`, `click`, `ready`, `pong`, `debug_snapshot`,
  `command_error`, `data_mapping`, `color_by`, `result_removed` and
  `results`: exactly the keys of their type, and no viewer id.
- Binary frames for `selection`, `result` and `session_bundle`: a
  little-endian uint32 header length, a UTF-8 JSON header
  (`{"type": "selection", "source": …}`, a result's header, or
  `{"type": "session_bundle", "requestId": …}`), then the body. A result's
  body is its cells and then its values ({doc}`../g_api_reference_coverage/api/round_trip`). A selection's
  body is the cell indices as little-endian uint32, strictly ascending; a
  session bundle's is the `.cellucid-session` file, streamed to a temporary
  file as it arrives.
- The viewer sends only the event types `/_cellucid/protocol` lists.
- Each event is checked (`_wire._require_inbound_jupyter_event`) and then
  delivered to the viewer's hooks on the socket's thread, one at a time and
  in order. A frame that breaks the rules closes the channel (1008) before
  any hook runs. A hook that raises is logged (logger
  `cellucid.viewer_channel`, "A viewer event hook failed") and the channel
  stays open. The event was recorded in `viewer.state` and handed to
  `wait_for_event` before any hook ran, and the hooks after the failing one
  still run.

## Lifetime

| Close code | Meaning | What the viewer does |
| --- | --- | --- |
| 1001 | the viewer or its server was stopped (`viewer.stop()`, which also runs when Python exits; `server.stop()`) | freezes the last view and says so |
| 4401 | unknown viewer or wrong token (a page left over from before a kernel restart) | freezes the view; run the cell again |
| 4409 | the same viewer opened a newer channel (its page reloaded, or the output is shown twice) | freezes this copy and says the viewer is open elsewhere; it does not take the channel back |
| 1008 | a frame broke the protocol | freezes the view and names the reason |
| anything else | the network or a proxy dropped the socket | reopens once, at once; waiting events are sent in order (a newer hover replaces a waiting one) |

- Every close the server makes is orderly (RFC 6455 §7.1.1): it sends the
  close frame, shuts its sending side, and reads and drops what the page
  still sends until the page closes its side, for at most 5 s, before it
  closes the socket. The page reads the code and reason it was sent, on
  every platform; a socket closed with unread input is reset instead, and
  on Windows the reset would discard the close frame (the page would see
  1006).
- The server sends a ping every 20 s, which the browser answers, so an idle
  channel survives proxies that close silent connections.
- Nothing is polled. If a reopen fails, the viewer makes one health request:
  no answer means Python stopped; an answer means the server runs but the
  WebSocket cannot reach it, and the message names the proxy as the cause.

## Viewer identity and routing

Each viewer has its own `viewerId` and a secret `viewerToken`: a widget hands
them to its page in a message, a browser tab carries them in its address.
Python keeps a process-wide registry
`viewerId → (token, callback, outbox, socket)` (`cellucid._viewer_channel`).
After a kernel restart the registry is empty, so an old page's `hello` is
refused (4401) and the old view freezes: run the cell again.

## A browser tab's address

A browser-tab viewer's page is
`http://127.0.0.1:<port>/?viewer=true&viewerId=…&viewerToken=…&anndata=true`
(`viewer.viewer_url`). When the browser cannot reach that address directly,
forward the port over SSH, or pass the browser-reachable address of a reverse
proxy with `client_server_url=` (such as
`https://<hub>/user/<name>/proxy/<port>`) and name the proxy's host in
`allowed_hosts=`. The channel then runs at
`wss://<hub>/user/<name>/proxy/<port>/_cellucid/ws`. A notebook widget needs
none of this.

## The web build

A widget's page is the app at `web_source_url` (`https://www.cellucid.com` by
default), loaded by the browser. A browser tab loads the app from its local
server, which establishes a verified generation from `web_source_url` before
it binds: every file is hash-checked on disk and published atomically.
Controls: `web_cache_dir=`, `viewer.ensure_web_ui_cached(force=True)`,
`viewer.clear_web_cache()`.

## Where to look in the code

- Python channel: `cellucid-python/src/cellucid/_viewer_channel.py` (handshake,
  outbox, events, session bundles) over `_websocket.py` (RFC 6455 framing).
- Python checks: `cellucid-python/src/cellucid/jupyter/_wire.py`.
- Python viewer: `cellucid-python/src/cellucid/jupyter/_viewer.py`.
- Web channel: `cellucid/assets/js/data/viewer-channel.js` and
  `viewer-protocol.js`; the bridge: `data/viewer-host-source.js`.
- The rules both sides follow: `cellucid/tests/fixtures/viewer-protocol.json`
  (cellucid-python and cellucid-r keep byte-identical copies as
  `tests/fixtures/viewer_protocol.json` and
  `tests/testthat/fixtures/viewer_protocol.json`).

## Next steps

- Commands: {doc}`06_python_to_frontend_commands`
- Events: {doc}`07_frontend_to_python_events`
- Security: {doc}`12_security_model` and {doc}`13_security_cors_origins_and_mixed_content`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
