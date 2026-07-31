# Server mode architecture, endpoints, and security

This page documents the Python HTTP servers that power Cellucid:

- serving *export folders* (static files), and
- serving *AnnData directly* (dynamic “virtual files”).

It also covers the web-asset proxy, event endpoints (hooks), and security caveats.

---

## The two servers (and why both exist)

### 1) Exported dataset server (`CellucidServer`)

Use when:
- you already ran `prepare(...)`,
- you want fastest load times,
- you want to share/host a stable export folder.

Implementation:
- `cellucid-python/src/cellucid/server.py`

### 2) AnnData server (`AnnDataServer`)

Use when:
- you want convenience during iteration,
- you want read-only backed access to `.h5ad`, or an eager load from `.zarr`,
- you don’t want to write an export folder yet.

Implementation:
- `cellucid-python/src/cellucid/anndata_server.py`
- `cellucid-python/src/cellucid/anndata_adapter.py` (the “virtual export” layer)

Both servers share:
- CORS rules
- exact web-generation establishment and delivery
- hooks event endpoint (`/_cellucid/events`)
- session bundle upload endpoint (`/_cellucid/session_bundle`)

Shared code:
- `cellucid-python/src/cellucid/_server_base.py`

---

## Common endpoints

These endpoints exist in both server modes:

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

### `POST /_cellucid/events`

Frontend → Python hooks/events channel:

- browser posts JSON events here
- server routes by `viewerId` to the correct notebook viewer callback

Implementation:
- `CORSMixin.handle_event_post` in `src/cellucid/_server_base.py`

### `POST /_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`

Frontend → Python session capture channel:

- viewer uploads raw `.cellucid-session` bytes here
- server validates + streams to a temp file
- server routes a `session_bundle` event back to Python

Implementation:
- `CORSMixin.handle_session_bundle_post` in `src/cellucid/_server_base.py`

---

## Exported server details

### Dataset listing: `GET /_cellucid/datasets`

The exported server can serve:
- a single dataset directory, or
- a directory containing multiple dataset subdirectories.

`/_cellucid/datasets` lists what’s available (used for multi-dataset browsing).

Detection rules:
- a dataset is a directory containing a valid version-2
  `dataset_identity.json`, an object-valued `obs_manifest.json`, and at least
  one non-empty supported `points_1d`, `points_2d`, or `points_3d` binary;
- once any immediate child is a dataset candidate, every immediate
  subdirectory must be one complete current dataset. Stray directories reject
  the root;
- `CellucidServer.viewer_url` opens the catalog. A sole entry may auto-open;
  multiple entries require an exact id and are never resolved by list order.

### Static file serving

Exported datasets are served as static files.

That means:
- there is no on-the-fly recomputation,
- and `.gz` files are served as raw gzip bytes (the viewer decompresses them).

---

## AnnData server details

The AnnData server implements a set of HTTP routes that mimic the export format.

Examples:

- `/dataset_identity.json`
- `/obs_manifest.json`
- `/var_manifest.json`
- `/points_3d.bin`
- `/obs/<field>.values.f32`
- `/var/<gene>.values.f32`

Key behaviors:

- **loading behavior**: H5AD is read-only backed; Zarr is loaded eagerly
- **gzip behavior**: an exact unsuffixed binary route may return
  `Content-Encoding: gzip` when the request advertises gzip support; `.gz` route
  variants return `404`
- **one rooted route contract**: data paths such as `/points_3d.bin` are served
  exactly at the root; dataset-prefixed variants return `404`.

---

## Verified web generation (why the UI loads from the server)

Both servers establish the configured exact source generation:

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

Exported mode:
- confirm you are serving the directory that actually contains `obs_manifest.json` and points files.

AnnData mode:
- confirm the requested path matches server routes (see `AnnDataRequestHandler`).

Use:

```bash
curl -I http://127.0.0.1:8765/obs_manifest.json
curl -I http://127.0.0.1:8765/points_3d.bin
```

### Symptom: “Hooks don’t fire”

Confirm:
- the viewer is posting to `/_cellucid/events` (check browser network tab),
- the event includes the correct `viewerId`, `viewerToken`, and `type`,
- Python registered the callback (viewer is still alive).

Start with: {doc}`11_hooks_events_protocol_and_schema` and {doc}`12_debugging_playbook`.
