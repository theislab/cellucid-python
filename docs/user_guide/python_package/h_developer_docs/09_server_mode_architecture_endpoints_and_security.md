# Server mode architecture, endpoints, and security

This page documents the Python HTTP server that powers Cellucid: it serves
*one dataset* — a file or folder in any format Cellucid reads, or an in-memory
AnnData or MuData — whose bytes the web app reads itself.

It also covers the web-asset delivery, the viewer channel (hooks), and
security caveats.

---

## The server (`AnnDataServer`)

Every entry point — `cellucid serve`, `serve_anndata`, `show_anndata` — runs
it, with nothing converted, for any format of the input format registry
(`input_formats.json`, shared byte for byte with the web app and cellucid-r):

| `kind` | What | Served as |
|---|---|---|
| `h5ad` | `.h5ad` | one file |
| `zarr` | AnnData Zarr (v2 or v3) | a folder key by key, or one `.zarr.zip` |
| `h5mu` | MuData `.h5mu` | one file |
| `mudata-zarr` | MuData Zarr | a folder key by key, or one `.zarr.zip` |
| `loom` | `.loom` | one file |
| `tenx-h5` | 10x Cell Ranger `.h5` | one file |
| `tenx-outs` | a Cell Ranger `outs/` folder | a folder of files, listed |
| `mtx` | a MatrixMarket folder (`matrix.mtx[.gz]`, or GEO's `GSM…_matrix.mtx.gz` naming) | a folder of files, listed |
| `rds`, `rda` | R objects `.rds`, `.rda`/`.RData`, gzipped or not | one file |

A file is taken by its name and is then what its first bytes and HDF5 root
say (an HDF5 root marked `anndata` is an `h5ad`, one marked `mudata` an
`h5mu`); a folder is a Zarr store by its root metadata (a MuData by its
root's `encoding-type`), else what its files make it. An AnnData file not in
priority order, and an in-memory `AnnData`, are served as a Zarr store in
priority order made as the browser reads it (`_served_store.PriorityStore`); a
`MuData` object is written once to a temporary Cellucid-ready `.h5mu`
(through `prepare`, one modality at a time).

Implementation:
- `cellucid-python/src/cellucid/anndata_server.py` (routes, launch choices)
- `cellucid-python/src/cellucid/_anndata_files.py` (classifying the path,
  reading AnnData and MuData metadata, the temporary `.h5mu`)
- `cellucid-python/src/cellucid/_served_store.py` (the store in priority
  order: chunks made on request from an object or a file)
- `cellucid-python/src/cellucid/_viewer_channel.py` and `_websocket.py` (the
  viewer channel)

The server is a threaded HTTP/1.1 server with keep-alive connections, because
the browser reads one file with many range requests.

Its shared building blocks:
- CORS rules
- exact web-generation establishment and delivery
- the `Host` check and request draining

Shared code:
- `cellucid-python/src/cellucid/_server_base.py`

---

## Common endpoints

These endpoints are always served:

### `GET /`

Serves the exact verified Cellucid web generation.

### `GET /_cellucid/health`

Health check for debugging:

```bash
curl -s http://127.0.0.1:8765/_cellucid/health | python -m json.tool
```

### `GET /_cellucid/info`

Server metadata (version, format, counts, the bound address, and the program
serving the viewer):

```bash
curl -s http://127.0.0.1:8765/_cellucid/info | python -m json.tool
```

(python-protocol-capability-endpoint)=
### `GET /_cellucid/protocol`

Where the viewer channel is, and what this installation accepts on it:

```bash
curl -s http://127.0.0.1:8765/_cellucid/protocol | python -m json.tool
```

```json
{
  "channel": "_cellucid/ws",
  "commands": ["clearHighlights", "debug_snapshot", "freeze", "highlight",
               "ping", "requestSessionBundle", "resetCamera", "setColorBy",
               "setDataMapping", "setField", "setGeneSet", "setVisibility"],
  "binaryCommands": ["highlight", "setField", "setVisibility"],
  "events": ["click", "color_by", "command_error", "data_mapping",
             "debug_snapshot", "hover", "pong", "ready", "result",
             "result_removed", "results", "selection", "session_bundle"]
}
```

- `channel` — where the viewer opens its WebSocket, relative to the server's
  base address. A viewer refuses a server without it and asks for an upgrade
  of the package.
- `events` — every event type `_require_inbound_jupyter_event` accepts on the
  channel. A viewer sends only the types named here.
- `commands` — every command this installation may send to the viewer, which
  is also exactly the set of names a `command_error` may report as having
  failed.
- `binaryCommands` — the commands it sends as binary frames (those that name
  cells). A viewer refuses a server without it and asks for an upgrade.

The lists grow; there is deliberately no version number to compare against,
so a viewer tests for the capability it needs rather than inferring it from a
release number. Read
{doc}`11_hooks_events_protocol_and_schema` for why the route exists at all.

The body is a constant of the installed package — no viewer id, no token, no
dataset — so this route needs no credentials, exactly like `/_cellucid/health`
and `/_cellucid/info`. It is answered from the same declaration the validator
reads, so it cannot announce a capability the notebook would then reject.

Implementation:
- `_protocol_capability_document` in `src/cellucid/jupyter/_wire.py`

### `GET /_cellucid/ws` (the viewer channel)

One WebSocket per notebook viewer, carrying Python's commands to the viewer
and the viewer's events back
({doc}`../e_jupyter_hooks/05_architecture_viewer_channel`):

- the upgrade passes the `Host` check (`421`) and an `Origin` check (`403`):
  the server's own host, a loopback name, or a host in `allowed_hosts`;
- the first frame is `hello` with the viewer's id and token, answered
  `welcome`, or closed with `4401`;
- commands are JSON text frames, and binary frames for `highlight`,
  `setField` and `setVisibility` (uint32 cell indices, read in place by the
  viewer; a pushed column's values after them), queued
  in the viewer's outbox until it connects and flushed in order;
- events are JSON text frames (at most 1 MiB) and binary frames for
  `selection` (uint32 cell indices, read in place), `result` (a result's
  cells and values) and `session_bundle` (streamed to a temporary file);
- the server pings every 20 s; `viewer.stop()` or the server stopping closes
  the channel with `1001`, a newer channel of the same viewer replaces it with
  `4409`, and a frame that breaks the rules closes it with `1008` (`1002`,
  `1007` or `1009` for a frame that breaks the WebSocket framing, is not
  UTF-8, or is a text frame over 1 MiB) before any hook runs;
- every close is the orderly one of RFC 6455 §7.1.1: close frame, then the
  write side shut down, then what the page still sends read and dropped until
  it closes its side (at most `CLOSE_DRAIN_SECONDS`, 5 s), and only then the
  socket closed. A thread waiting to read notices a close made from another
  thread within 0.25 s. Closing with input unread would reset the connection,
  and on Windows the reset discards the close frame before the page reads it.

Implementation:
- `accept_viewer_channel` in `src/cellucid/_viewer_channel.py`, over the
  RFC 6455 framing in `src/cellucid/_websocket.py`

The server answers no `POST`: every one is refused with `405`.

---

## AnnData server details

The server computes nothing per request. The web app's own reader reads the
data; the server names the dataset, serves its bytes, and states the launch
choices.

### Startup

`AnnDataServer(...)` runs three numbered steps before `start()` adds the
fourth:

1. `[1/4] Detecting format` — classify the path (the table above) and resolve
   the dataset identity: `dataset_name`/`dataset_id` as given, else
   `uns['cellucid']['name']`/`['id']` of an AnnData or MuData; data with
   neither is refused (an R object, a Loom, 10x or MatrixMarket dataset always
   needs both named). An in-memory AnnData, and an AnnData file not in
   priority order (or in it, but holding an array whose codec the browser
   cannot decode), is planned as a Zarr store in priority order: `obs`, `var`,
   `uns`, `varm`, `varp` written in memory, every per-cell array's chunks made
   on request, uncompressed; an in-memory MuData is written to a temporary
   Cellucid-ready `.h5mu`. Each is released when the server stops or Python
   exits. Another format holding an array the browser cannot decode is
   refused, naming the arrays.
2. `[2/4] Reading metadata` — for an AnnData or MuData, cell and gene counts
   and the names of the stored elements (a MuData's shared cells, its
   modalities, and the default modality's genes); no array is read. Other
   formats are read only by the browser.
3. `[3/4] Checking launch choices` — every choice named (`embedding_keys`,
   `knn_key`, `expression_source`, `gene_id_column`, `latent_key`,
   `vector_field_default`, `obs_keys`, `normalize_embeddings`,
   `centroid_min_points`, `modality`) must have its exact form, and, for an
   AnnData or MuData, exist in the store with a usable shape (`modality` only
   for a MuData); for the other formats the viewer checks them when it opens
   the data. For a Zarr folder without consolidated metadata, the server
   builds it in memory here.
4. `[4/4] Starting server` — establish the viewer generation and bind.

### `GET /_cellucid/datasets`

One dataset, its store, and the launch choices that were set:

```json
{
  "datasets": [
    {
      "id": "pbmc-demo",
      "name": "PBMC demo",
      "anndata": {
        "file": {"url": "_cellucid/anndata/pbmc-demo.h5ad", "kind": "h5ad",
                 "layout": "file", "size": 2799976,
                 "etag": "\"2ab9e8-18e27b3d58754f73-18e27b3d58754f73-1a058bb9\""},
        "defaults": {"embedding2d": {"key": "X_umap", "columns": [1, 0]},
                     "knn": "connectivities"}
      }
    }
  ]
}
```

- `file.url` is relative to the server's base URL (the server cannot know the
  address a browser reaches it by, through a tunnel or a proxy). `kind` is the
  registry's format id (the table above); `layout` is `"file"`, `"zip"` or
  `"folder"`; `size` and `etag` are given for a file and absent for a
  folder. Both are read from the file when the list is asked for: `etag` is
  the strong validator its bytes are answered with (below), and the viewer
  takes it as the file's version -- what it keeps of the file (its byte
  cache, a kept gene index, a converted store) is filed under it and kept
  across visits, a rewritten file is a new version, and a response with
  another `ETag` than the one listed is refused rather than mixed in. A file
  that is gone is answered `404` (`<name> is no longer readable`).
- A folder of files (`tenx-outs`, `mtx`) also lists them, because a browser
  cannot list a folder: `"files": [{"path": "GSM1_matrix.mtx.gz", "size":
  1234, "etag": "\"4d2-…\""}, …]`, every regular file below it with its
  exact size and validator (hidden files left out). A Zarr store folder lists nothing; its consolidated metadata
  names its keys.
- `defaults` holds only the launch choices given at startup. A choice left out
  is taken from the store's own `uns['cellucid']['defaults']`, else the viewer's
  default rule ({doc}`../../web_app/b_data_loading/13_data_mapping`).

### `GET /_cellucid/health` and `GET /_cellucid/info`

`/_cellucid/health` answers
`{"status": "ok", "type": "anndata", "version", "format", "n_cells", "n_genes"}`,
with `format` the registry's format id. `n_cells` and `n_genes` are `null` for
a format Python does not read (an R object, a Loom, 10x or MatrixMarket
dataset): the browser learns them when it opens the data. `/_cellucid/info`
carries the same `version`, `type`, `format`, `n_cells` and `n_genes`, plus
`host` and `port`, and `runtime`: `"python"` for `cellucid serve` and a viewer
outside a notebook, `"jupyter"` for a viewer in a notebook kernel (the R
package's viewer says `"r"`). The viewer names its host by it -- "Python",
"Jupyter" or "R" in the Dataset box -- with the names the shared protocol
file declares (`viewer-protocol.json`, `runtimes`).

### The store: `/_cellucid/anndata/<dataset_id>…`

The route is named by the dataset id, so the local file name is never
published:

- `/_cellucid/anndata/<dataset_id><ending>` — one file, under the ending its
  format names it by (`.h5ad`, `.h5mu`, `.loom`, `.h5`, `.rds.gz`, …);
- `/_cellucid/anndata/<dataset_id>.zarr.zip` — a zipped Zarr store;
- `/_cellucid/anndata/<dataset_id>.zarr/<key>` — a Zarr store folder, key by
  key. An absent key answers `404`, which a Zarr reader takes as "absent" (an
  unwritten chunk is its fill value); a key may not leave the store;
- `/_cellucid/anndata/<dataset_id>/<path>` — one file of a folder of files;
  an absent or escaping path answers `404`.

Every answer is exact: `GET` and `HEAD`, a single `Range` answered with `206`
and `Content-Range`, an unsatisfiable or multiple range refused with `416`,
`Accept-Ranges: bytes`, and an `ETag`. CORS exposes `Accept-Ranges`,
`Content-Encoding`, `Content-Length`, `Content-Range` and `ETag`, so the
viewer's reader can see them. The viewer sends nothing but one `Range`
header, a CORS-safelisted request, so no preflight is made.

What a request costs is kept small, because a browser reads a file with
thousands of them:

- a file's bytes go from the page cache to the socket with `sendfile`
  (`socket.sendfile`), never read into Python, and are never compressed;
- `TCP_NODELAY` on every connection, so a small answer is not held back
  waiting for the client to acknowledge its headers;
- the listen backlog is the platform's maximum, so the six connections a page
  opens at once never wait for a retry.

The browser may keep what it read: store bytes are sent with
`Cache-Control: no-cache` and a strong `ETag` (the file's size, modification
time, status-change time and inode, read from the open file — on Windows the
status-change time is NTFS's change time, read from the file's handle, not
`st_ctime`, which is the creation time there: a file rewritten in place, even
with its modification time put back, or replaced by a copy moved over it gets
a new one; only two writes of the same size within the file system's timestamp
resolution — nanoseconds on APFS, ext4 and NTFS — could share one), so a reload revalidates them (`If-None-Match`, answered
`304 Not Modified` with no body) instead of downloading them again. A file that
changed has a new `ETag`: a range asked with `If-Range` naming the old one is
answered with the whole current file, and the viewer refuses to mix two
versions of a file it has open.

```bash
curl -sI -H 'Range: bytes=0-7' http://127.0.0.1:8765/_cellucid/anndata/pbmc-demo.h5ad
```

Implementation:
- `AnnDataRequestHandler` in `src/cellucid/anndata_server.py`; the web app's
  reader of this document is `assets/js/data/served-anndata.js`.

---

## Verified web generation (why the UI loads from the server)

The server establishes the configured exact source generation:

- at every `start()`, before it binds, the Python server reads the inventory
  `cellucid-web-assets.json` from the web source (`https://www.cellucid.com`
  by default) and downloads every object it declares,
- byte-verifies them and publishes the complete generation atomically in the
  web cache directory,
- and then serves `/index.html` and `/assets/...` from that verified
  generation, from the same origin as the dataset routes (an object the
  inventory does not declare is a `404`).

Why this exists:
- the page and the data share one origin, so a browser tab needs no CORS and
  meets no mixed content (a notebook widget does not use it: its page is the
  app at `web_source_url`, and its reads go through the widget).

Generation-directory configuration:
- `--web-cache-dir` or `web_cache_dir=` (see
  {doc}`06_configuration_env_vars_and_logging`)

If source download, verification, or atomic publication fails, server startup
raises and does not bind. It never serves a different cached generation.

---

## Security model (practical, not theoretical)

### Default stance: local/trusted

By default the server binds to `127.0.0.1` (local only). This is the recommended mode for private datasets.

The bind address alone is not the boundary: a page in another tab can re-point
its own DNS name at `127.0.0.1` (DNS rebinding) and the browser then treats this
server as same-origin. Cellucid therefore validates the `Host` header on every
request ahead of routing, accepting only one well-formed authority that names
`localhost` or a loopback IP literal on the bound port, and answering anything
else with `421 Misdirected Request`.

Behind a reverse proxy (a hub's proxy route) the forwarded `Host` is
the proxy's own, which is byte-for-byte what a rebound page sends, so the proxy
must be declared rather than detected:
`--allowed-host hub.example.org` or `allowed_hosts=["hub.example.org"]`. Entries
are bare host names with no port, scheme or wildcard, and are matched on any
port. See
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

If you bind to `0.0.0.0` or otherwise expose the server:
- you are making your dataset accessible to anyone who can reach that host/port,
- the `Host` check is off unless you declare `allowed_hosts` (a non-loopback
  bind is an explicit request for exposure whose names the server cannot know),
- while the viewer channel stays closed to anyone without a viewer's token.

### CORS is not full authentication

CORS helps browsers decide whether JS from one origin can read responses from another.
It is not a firewall and does not stop:
- direct HTTP clients (`curl`),
- same-origin attackers,
- or a user with network access.

### The viewer channel uses exact per-viewer authentication

CORS does not apply to WebSockets, so the channel's upgrade checks the page's
`Origin` itself (the server's own host, a loopback name, or a host in
`allowed_hosts`; `403` otherwise), after the `Host` check. Its first frame must
then name a registered viewer and its exact token, compared in constant time;
an unknown viewer or a wrong token closes the channel with `4401` before
anything passes. Every later frame is checked against its closed record.

### The capability route carries no credentials, because it carries no secret

`GET /_cellucid/protocol` answers anyone who can reach the port, exactly as
`/_cellucid/health` and `/_cellucid/info` do. Its body is the same for every
viewer, every dataset and every request — it is a property of the installed
`cellucid` version, not of your data — so there is nothing in it a token would
be protecting. It is also read-only: no state is created, consumed or
authenticated. The `Host` check and the CORS rules above apply to it like any
other route.

### A session bundle answers one request, once

A `session_bundle` frame is accepted only on an authenticated channel and only
for a request Python registered (`register_session_bundle_request`) and is
still waiting for. The server streams it to a temporary file, validates the
complete bundle structure, consumes the request exactly once, and only then
delivers it; an invalid bundle is removed and reported as a `command_error`.
A bundle nobody is waiting for is read and dropped.

For a deeper threat-model discussion:
{doc}`17_security_privacy_and_networking`.

---

## Troubleshooting

### Symptom: “I can open `/` but dataset files 404”

- read `file.url` from `/_cellucid/datasets` and request exactly that route
  (see `AnnDataRequestHandler`); any other path under `/_cellucid/anndata/` is a
  `404`.

```bash
curl -s http://127.0.0.1:8765/_cellucid/datasets | python -m json.tool
```

### Symptom: “Hooks don’t fire”

Confirm:
- the viewer's channel is open (`viewer.debug_connection()["channel_connected"]`,
  or the `/_cellucid/ws` connection in the browser's Network panel, filter `WS`),
- Python registered the callback (viewer is still alive).

Start with: {doc}`11_hooks_events_protocol_and_schema` and {doc}`12_debugging_playbook`.
