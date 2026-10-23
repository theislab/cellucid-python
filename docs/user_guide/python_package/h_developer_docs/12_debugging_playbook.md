# Debugging playbook

This page is a **step-by-step debugging runbook** for Cellucid Python developers and power users.

It is optimized for two outcomes:

1) you can fix the issue yourself, or
2) you can file a bug report that a maintainer can act on without guessing.

---

## Step 0 — Identify which mode you are in

Everything downstream depends on this.

Pick the closest match:

### Mode A: `prepare()` failed

- `cellucid.prepare(...)` raised: read its message, which names the argument
  and the fix ({doc}`../c_data_preparation_api/11_troubleshooting_prepare`).

### Mode B: AnnData server

- you are serving `.h5ad`/`.zarr` (prepared or not) or an in-memory AnnData.
- you start the viewer with
  `cellucid serve ./data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`
  or `cellucid.show_anndata(...)`.

### Mode C: Notebook embedding

- you are in Jupyter/VSCode/Colab.
- you see an iframe (or a blank output where you expected an iframe).

---

## Step 1 — Confirm the server is alive

Open a terminal (or use `curl` from wherever the server is running) and check:

```bash
curl -s http://127.0.0.1:8765/_cellucid/health
curl -s http://127.0.0.1:8765/_cellucid/info
```

If these fail:
- your server is not running, crashed, or the port is wrong.

If the port differs, find it in:
- the CLI banner output, or
- `viewer.server_url` (notebook), or
- `server.url` (Python).

---

## Step 2 — Confirm the viewer UI assets can load

Start the viewer/server. Source-generation errors occur before the HTTP server
binds and identify the failing inventory, object, verification, or filesystem
operation.

Fix options:
1) ensure access to the configured `web_source_url`,
2) correct the exact source response or inventory contract, and
3) pass a writable `web_cache_dir`.

See: {doc}`09_server_mode_architecture_endpoints_and_security` and {doc}`06_configuration_env_vars_and_logging`.

---

## Step 3 — Confirm dataset endpoints return data

### The store route

The server serves one store; the viewer reads it with range requests. Read the
store's route from `/_cellucid/datasets`, then probe it the way the viewer
does — a readable store answers a `Range` request with `206` and
`Content-Range`:

```bash
curl -s http://127.0.0.1:<port>/_cellucid/datasets | python -m json.tool
curl -sI -H 'Range: bytes=0-7' http://127.0.0.1:<port>/_cellucid/anndata/<dataset_id>.h5ad
```

For a Zarr store folder the route ends in `.zarr/`, and each key is its own
path (`…/<dataset_id>.zarr/.zmetadata`, or `zarr.json` for Zarr v3); an
absent key is a `404`. `viewer.debug_connection()` runs the same probe and
reports it under `anndata_store_probes`.

---

## Step 4 — Use browser DevTools (network + console)

This is the fastest way to distinguish:
- server failures (404/500),
- format issues (an array the reader cannot decode, a stale Zarr metadata document),
- and frontend/UI issues.

### 4.1 Network tab checklist

In browser DevTools → Network:

1) filter for `_cellucid`:
   - `/_cellucid/health` should be 200
   - `/_cellucid/info` should be 200
2) filter for `_cellucid/anndata/`:
   - range requests should answer `206` with `Content-Range`
3) look for 404s:
   - for a Zarr store folder, a 404 on a chunk key is normal (an unwritten
     chunk is its fill value); a 404 on `.zmetadata` or `zarr.json` means the
     route is wrong
4) look for CORS errors:
   - “CORS blocked” usually means origin mismatch or a proxy issue

### 4.2 Console tab checklist

In DevTools → Console:

- copy/paste the first relevant error (don’t paraphrase)
- look for:
  - “Failed to fetch”
  - “CORS blocked”
  - refusals naming `uns['cellucid']` (a key or value the viewer does not accept)

If you suspect web-app internal bugs, also consult:
- {doc}`../../web_app/p_developer_docs/13_debugging_playbook`

---

## Step 5 — Notebook-specific debugging (`viewer.debug_connection()`)

If you have a `viewer` object (Jupyter/VSCode/Colab), run:

```python
report = viewer.debug_connection()
report
```

The report includes:
- server health/info probes,
- `anndata_store_probes`, keyed by every dataset the server lists, recording
  each store's range-read answer or failure,
- whether the viewer is displayed,
- web UI cache status + build id,
- a Python→frontend ping/pong roundtrip,
- a frontend “debug snapshot” (URL/origin/userAgent as seen by the iframe),
- and recent accepted-event counts grouped by exact event type.

If you are filing a bug report, paste the report as JSON:

```python
import json
print(json.dumps(report, indent=2, default=str))
```

---

## Common failure patterns (fast diagnosis)

### Pattern: “Viewer UI loads, but data is missing”

- the store route answers something other than `206` to a `Range` request
  (probe it as above)
- the viewer's reader refused part of the store (read the notification and
  the browser console)

Start with:
- {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`
- {doc}`07_prepare_writer_architecture`

### Pattern: “Hooks don’t fire”

Confirm in browser Network tab:
- `POST /_cellucid/events` is happening and returns 200.

If it is not happening:
- frontend may not be configured for notebook mode,
- iframe may be blocked by proxy/mixed-content issues.

If it is happening:
- viewerId may mismatch (stale iframe vs new viewer),
- Python viewer may have been stopped/garbage-collected.

Start with:
- {doc}`11_hooks_events_protocol_and_schema`
- {doc}`10_jupyter_embedding_architecture`

### Pattern: “Everything works locally, fails on a remote server”

Likely causes:
- binding host incorrectly (`127.0.0.1` vs `0.0.0.0`)
- missing SSH tunnel
- corporate firewall
- HTTPS notebook mixed-content

Fix:
- prefer SSH tunneling over public binding,
- or provide a stable reverse proxy and pass its exact base as
  `client_server_url=`.

---

## What to include in a high-quality bug report

Minimum:

- mode (`prepare()` vs AnnData server vs notebook)
- `cellucid` version + Python version + OS
- reproduction steps (code or click-by-click)
- expected vs actual
- server URL + what `/_cellucid/health` returns
- browser console error(s)

Notebook:
- attach `viewer.debug_connection()` report.
