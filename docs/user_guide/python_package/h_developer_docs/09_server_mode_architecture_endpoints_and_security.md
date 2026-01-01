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
- you want lazy loading from `.h5ad` backed mode or `.zarr`,
- you don’t want to write an export folder yet.

Implementation:
- `cellucid-python/src/cellucid/anndata_server.py`
- `cellucid-python/src/cellucid/anndata_adapter.py` (the “virtual export” layer)

Both servers share:
- CORS rules
- hosted-asset proxy (web UI caching)
- hooks event endpoint (`/_cellucid/events`)
- session bundle upload endpoint (`/_cellucid/session_bundle`)

Shared code:
- `cellucid-python/src/cellucid/_server_base.py`

---

## Common endpoints

These endpoints exist in both server modes:

### `GET /`

Serves the Cellucid web UI (typically via hosted-asset proxy).

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

### `POST /_cellucid/session_bundle?viewerId=...&requestId=...`

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
- a dataset is a directory containing `obs_manifest.json` and at least one `points_*d.bin(.gz)`

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
- `/points_3d.bin` or `/points_3d.bin.gz`
- `/obs/<field>.values.f32` (and `.gz`)
- `/var/<gene>.values.f32` (and `.gz`)

Key behaviors:

- **lazy loading** for backed `.h5ad` and `.zarr` when possible
- **gzip behavior**:
  - `.gz` suffix returns raw gzip bytes (export-folder semantics)
  - requests without `.gz` may return `Content-Encoding: gzip` if the browser supports it
- **dataset id prefix stripping**: the frontend may request paths like `<dataset_id>/points_3d.bin`; the server strips the prefix.

---

## Hosted-asset proxy (why the UI loads “from the server”)

Both servers are typically configured to run in hosted-asset proxy mode:

- when the browser requests `/index.html` or `/assets/...`,
- the Python server downloads those files from `https://www.cellucid.com`,
- caches them on disk,
- and serves them from the same origin as the dataset server.

Why this exists:
- avoids HTTPS→HTTP mixed content in notebooks,
- avoids cross-origin issues in embedded iframes,
- enables offline reuse if the cache is pre-populated.

Cache configuration:
- `CELLUCID_WEB_PROXY_CACHE_DIR` (see {doc}`06_configuration_env_vars_and_logging`)

If the proxy fails and no cached UI exists, the server returns a helpful HTML error page at `/index.html`.

<!-- SCREENSHOT PLACEHOLDER
ID: hosted-asset-proxy-error-page
Suggested filename: developer/hosted-asset-proxy-error.png
Where it appears: Python Package → Developer Docs → Server mode architecture
Capture:
  - UI location: browser tab opened at http://127.0.0.1:<port>/
  - State prerequisites: disconnect network AND clear cache dir (or use a fresh machine)
  - Action to reach state: start server, open viewer URL
Crop:
  - Include: the error message text and the “what to do” bullet list
  - Exclude: browser bookmarks/avatars, private paths
Alt text:
  - Browser page showing “Cellucid viewer UI could not be loaded” with steps to fix caching/network access.
Caption:
  - When the hosted UI cannot be fetched and no cached copy exists, the server shows an explicit error page explaining how to populate the cache.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the hosted-asset proxy error page.
:width: 100%

Hosted-asset proxy failure mode: the server explains how to populate the viewer UI cache for offline use.
```

---

## Security model (practical, not theoretical)

### Default stance: local/trusted

By default the server binds to `127.0.0.1` (local only). This is the recommended mode for private datasets.

If you bind to `0.0.0.0` or otherwise expose the server:
- you are making your dataset accessible to anyone who can reach that host/port,
- and you should assume events and session uploads can be triggered by others on the network.

### CORS is not full authentication

CORS helps browsers decide whether JS from one origin can read responses from another.
It is not a firewall and does not stop:
- direct HTTP clients (`curl`),
- same-origin attackers,
- or a user with network access.

### Hooks event endpoint is viewerId-routed but not authenticated

`POST /_cellucid/events` requires a `viewerId`, but does not currently enforce `viewerToken`.
Treat it as trusted-local plumbing.

### Session bundle upload has a “pending request” guard

`/_cellucid/session_bundle` requires:
- a `viewerId` and `requestId`,
- and a prior registration by Python (`register_session_bundle_request`).

It also:
- validates the magic header,
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
- the event includes the correct `viewerId`,
- Python registered the callback (viewer is still alive).

Start with: {doc}`11_hooks_events_protocol_and_schema` and {doc}`12_debugging_playbook`.
