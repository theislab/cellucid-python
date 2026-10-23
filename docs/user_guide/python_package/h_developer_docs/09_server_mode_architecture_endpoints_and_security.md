# Server mode architecture, endpoints, and security

This page documents the Python HTTP server that powers Cellucid: it serves
*one AnnData store* (`.h5ad` or Zarr) whose bytes the web app reads itself.

It also covers the web-asset delivery, event endpoints (hooks), and security caveats.

---

## The server (`AnnDataServer`)

Every entry point — `cellucid serve`, `serve_anndata`, `show_anndata` — runs
it, for an `.h5ad`, a Zarr store folder (v2 or v3), a zipped Zarr store, or an
in-memory `AnnData`, with nothing converted.

Implementation:
- `cellucid-python/src/cellucid/anndata_server.py` (routes, launch choices)
- `cellucid-python/src/cellucid/_anndata_files.py` (classifying the store,
  reading its metadata, the temporary `.h5ad`)

The server is a threaded HTTP/1.1 server with keep-alive connections, because
the browser reads one file with many range requests.

Its shared building blocks:
- CORS rules
- exact web-generation establishment and delivery
- hooks event endpoint (`/_cellucid/events`)
- selection delivery endpoint (`/_cellucid/selection`)
- session bundle upload endpoint (`/_cellucid/session_bundle`)

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

Server metadata (version, mode, counts):

```bash
curl -s http://127.0.0.1:8765/_cellucid/info | python -m json.tool
```

(python-protocol-capability-endpoint)=
### `GET /_cellucid/protocol`

The wire capabilities this installation accepts:

```bash
curl -s http://127.0.0.1:8765/_cellucid/protocol | python -m json.tool
```

```json
{
  "commands": ["clearHighlights", "debug_snapshot", "freeze", "highlight",
               "ping", "requestSessionBundle", "resetCamera", "setColorBy",
               "setDataMapping", "setVisibility"],
  "events": ["click", "command_error", "data_mapping", "debug_snapshot",
             "hover", "pong", "ready", "selection", "session_bundle"]
}
```

- `events` — every event type `_require_inbound_jupyter_event` accepts on
  `POST /_cellucid/events`. A viewer emits only the types named here.
- `commands` — every command this installation may send into the iframe, which
  is also exactly the set of names a `command_error` may report as having
  failed.

Both are lists that grow; the key set does not, and there is deliberately no
version number to compare against. A viewer validates payloads against exact
key sets, so a capability announced as a new key would break the builds this
route exists to protect, and a version string would be a proxy for the
capability rather than evidence of it. Read
{doc}`11_hooks_events_protocol_and_schema` for why the route exists at all.

The body is a constant of the installed package — no viewer id, no token, no
dataset — so this route needs no credentials, exactly like `/_cellucid/health`
and `/_cellucid/info`. It is answered from the same declaration the validator
reads, so it cannot announce a capability the notebook would then reject.

Implementation:
- `_protocol_capability_document` in `src/cellucid/jupyter/_wire.py`

### `POST /_cellucid/events`

Frontend → Python hooks/events channel:

- browser posts JSON events here (at most 1 MiB)
- server routes by `viewerId` to the correct notebook viewer callback
- `selection` and `session_bundle` records are refused with `400`: only their
  own routes may produce them

Implementation:
- `CORSMixin.handle_event_post` in `src/cellucid/_server_base.py`

### `POST /_cellucid/selection?viewerId=...&viewerToken=...&source=...`

Frontend → Python selection channel:

- viewer posts the selected row indices as little-endian `uint32`, strictly
  ascending, one or more (`Content-Type: application/octet-stream`; no size
  limit beyond the format's own, every `uint32` index once, and the server's
  memory: a body it cannot hold is refused with `413`)
- server authenticates before reading the body, refuses a malformed body
  with `400`, and routes one `selection` event back to Python

Implementation:
- `CORSMixin.handle_selection_post` in `src/cellucid/_server_base.py`

### `POST /_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`

Frontend → Python session capture channel:

- viewer uploads raw `.cellucid-session` bytes here
- server validates + streams to a temp file
- server routes a `session_bundle` event back to Python

Implementation:
- `CORSMixin.handle_session_bundle_post` in `src/cellucid/_server_base.py`

---

## AnnData server details

The AnnData server computes nothing per request. The web app's own HDF5/Zarr
reader reads the store in place; the server names the dataset, serves the
store's bytes, and states the launch choices.

### Startup

`AnnDataServer(...)` runs three numbered steps before `start()` adds the
fourth:

1. `[1/4] Detecting format` — classify the path (`.h5ad` file, zipped Zarr
   store, Zarr v2/v3 store folder) and resolve the dataset identity:
   `dataset_name`/`dataset_id` as given, else `uns['cellucid']['name']`/`['id']`;
   a store with neither is refused. An in-memory object, and a store holding an
   array whose codec the browser cannot decode, is written once to a temporary
   `.h5ad` (sparse `X`, layers and `raw.X` column-major, uncompressed), removed
   when the server stops or Python exits.
2. `[2/4] Reading metadata` — cell and gene counts and the names of the stored
   elements; no array is read.
3. `[3/4] Checking launch choices` — every choice named (`embedding_keys`,
   `knn_key`, `expression_source`, `gene_id_column`, `latent_key`,
   `vector_field_default`, `obs_keys`, `normalize_embeddings`,
   `centroid_min_points`) must exist in the store
   and have a usable shape. For a Zarr folder without consolidated metadata, the
   server builds it in memory here.
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
                 "layout": "file", "size": 2799976},
        "defaults": {"embedding2d": {"key": "X_umap", "columns": [1, 0]},
                     "knn": "connectivities"}
      }
    }
  ]
}
```

- `file.url` is relative to the server's base URL (the server cannot know the
  address a browser reaches it by, through a tunnel or a proxy). `kind` is
  `"h5ad"` or `"zarr"`; `layout` is `"file"`, `"zip"` or `"folder"`; `size` is
  given for a file and absent for a folder.
- `defaults` holds only the launch choices given at startup. A choice left out
  is taken from the store's own `uns['cellucid']['defaults']`, else the viewer's
  default rule ({doc}`../../web_app/b_data_loading/13_data_mapping`).

### `GET /_cellucid/health` and `GET /_cellucid/info`

`/_cellucid/health` answers
`{"status": "ok", "type": "anndata", "version", "format", "n_cells", "n_genes"}`,
with `format` `"h5ad"` or `"zarr"`. `/_cellucid/info` carries the same
`version`, `type`, `format`, `n_cells` and `n_genes`, plus `host` and `port`.

### The store: `/_cellucid/anndata/<dataset_id>…`

The route is named by the dataset id, so the local file name is never
published:

- `/_cellucid/anndata/<dataset_id>.h5ad` — an `.h5ad` file;
- `/_cellucid/anndata/<dataset_id>.zarr.zip` — a zipped Zarr store;
- `/_cellucid/anndata/<dataset_id>.zarr/<key>` — a Zarr store folder, key by
  key. An absent key answers `404`, which a Zarr reader takes as "absent" (an
  unwritten chunk is its fill value); a key may not leave the store.

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
time and inode), so a reload revalidates them (`If-None-Match`, answered
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

- when the browser requests `/index.html` or `/assets/...`,
- the Python server downloads those files from `https://www.cellucid.com`,
- caches them on disk,
- and serves them from the same origin as the dataset server.

Why this exists:
- avoids HTTPS→HTTP mixed content in notebooks,
- avoids cross-origin issues in embedded iframes.

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

Behind a reverse proxy such as `jupyter-server-proxy` the forwarded `Host` is
the proxy's own, which is byte-for-byte what a rebound page sends, so the proxy
must be declared rather than detected:
`--allowed-host hub.example.org` or `allowed_hosts=["hub.example.org"]`. Entries
are bare host names with no port, scheme or wildcard, and are matched on any
port. See
{doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

If you bind to `0.0.0.0` or otherwise expose the server:
- you are making your dataset accessible to anyone who can reach that host/port,
- and you should assume events and session uploads can be triggered by others on the network.

### CORS is not full authentication

CORS helps browsers decide whether JS from one origin can read responses from another.
It is not a firewall and does not stop:
- direct HTTP clients (`curl`),
- same-origin attackers,
- or a user with network access.

### Hooks event endpoint uses exact per-viewer authentication

`POST /_cellucid/events` requires a JSON object with non-empty `viewerId`,
`viewerToken`, and `type`. The server authenticates the token before delivering
the event exactly once. Unknown viewers return `404`; missing or incorrect
credentials return `403`.

`POST /_cellucid/selection` takes the same credentials in its query and checks
them before it reads the body, so an unauthenticated caller cannot make the
server read up to its 512 MiB limit.

### The capability route carries no credentials, because it carries no secret

`GET /_cellucid/protocol` answers anyone who can reach the port, exactly as
`/_cellucid/health` and `/_cellucid/info` do. Its body is the same for every
viewer, every dataset and every request — it is a property of the installed
`cellucid` version, not of your data — so there is nothing in it a token would
be protecting. It is also read-only: no state is created, consumed or
authenticated. The `Host` check and the CORS rules above apply to it like any
other route.

### Session bundle upload is token-bound and one-shot

`/_cellucid/session_bundle` requires:
- the exact `viewerId`, `viewerToken`, and `requestId` query fields,
- a matching pending request registered by Python
  (`register_session_bundle_request`), and
- `Content-Type: application/octet-stream`.

It also:
- authenticates the viewer token,
- consumes the pending request exactly once after full validation,
- validates the complete bundle structure,
- enforces a size limit,
- writes to disk in a streaming way to avoid memory blowups.

Still: do not expose it publicly for sensitive datasets.

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
- the viewer is posting to `/_cellucid/events` (check browser network tab),
- the event includes the correct `viewerId`, `viewerToken`, and `type`,
- Python registered the callback (viewer is still alive).

Start with: {doc}`11_hooks_events_protocol_and_schema` and {doc}`12_debugging_playbook`.
