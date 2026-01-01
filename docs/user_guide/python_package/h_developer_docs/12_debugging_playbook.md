# Debugging playbook

This page is a **step-by-step debugging runbook** for Cellucid Python developers and power users.

It is optimized for two outcomes:

1) you can fix the issue yourself, or
2) you can file a bug report that a maintainer can act on without guessing.

---

## Step 0 — Identify which mode you are in

Everything downstream depends on this.

Pick the closest match:

### Mode A: Export folder

- you ran `cellucid.prepare(...)` and have a directory export.
- you start the viewer with `cellucid serve ./export_dir` or `cellucid.show("./export_dir")`.

### Mode B: AnnData server

- you are serving `.h5ad`/`.zarr` or an in-memory AnnData.
- you start the viewer with `cellucid serve ./data.h5ad` or `cellucid.show_anndata(...)`.

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

Open:
- `http://127.0.0.1:<port>/`

If you see an error page “viewer UI could not be loaded”:
- the hosted-asset proxy could not fetch the UI and no cached copy exists.

Fix options:
1) ensure network access to `https://www.cellucid.com`,
2) run once while online to populate the cache,
3) set `CELLUCID_WEB_PROXY_CACHE_DIR` to a writable/persistent directory.

See: {doc}`09_server_mode_architecture_endpoints_and_security` and {doc}`06_configuration_env_vars_and_logging`.

---

## Step 3 — Confirm dataset endpoints return data

### Export folder mode

Check that the directory contains:

- `dataset_identity.json`
- `obs_manifest.json`
- at least one `points_*d.bin` or `points_*d.bin.gz`

Then probe:

```bash
curl -I http://127.0.0.1:<port>/dataset_identity.json
curl -I http://127.0.0.1:<port>/obs_manifest.json
curl -I http://127.0.0.1:<port>/points_3d.bin
curl -I http://127.0.0.1:<port>/points_3d.bin.gz
```

(One of the two points paths should exist depending on compression.)

### AnnData server mode

Probe:

```bash
curl -I http://127.0.0.1:<port>/dataset_identity.json
curl -I http://127.0.0.1:<port>/obs_manifest.json
curl -I http://127.0.0.1:<port>/points_3d.bin
```

If you want to test gene/obs endpoints, pick a known key and try:

```bash
curl -I http://127.0.0.1:<port>/var/FOXP1.values.f32
curl -I http://127.0.0.1:<port>/obs/cell_type.codes.u8
```

---

## Step 4 — Use browser DevTools (network + console)

This is the fastest way to distinguish:
- server failures (404/500),
- format/schema issues (manifest parse errors),
- and frontend/UI issues.

### 4.1 Network tab checklist

In browser DevTools → Network:

1) filter for `_cellucid`:
   - `/health` should be 200
   - `/info` should be 200
2) filter for `obs_manifest.json` and `var_manifest.json`:
   - should be 200 and have JSON content
3) look for 404s:
   - missing `.gz` files often means a partial export / compression mismatch
4) look for CORS errors:
   - “CORS blocked” usually means origin mismatch or a proxy issue

<!-- SCREENSHOT PLACEHOLDER
ID: devtools-network-failing-request
Suggested filename: developer/devtools-network-failing-request.png
Where it appears: Python Package → Developer Docs → Debugging playbook
Capture:
  - UI location: browser DevTools Network tab
  - State prerequisites: a failing load (404/500/CORS) reproduced
  - Action to reach state: reload viewer, reproduce the failure
Crop:
  - Include: request URL, status code, and the failing response/headers panel
  - Exclude: unrelated tabs, personal extensions, account avatars
Alt text:
  - Browser DevTools Network tab showing a failed request to a Cellucid dataset endpoint with status and headers visible.
Caption:
  - The Network tab reveals whether a failure is a missing file (404), a server exception (500), or an origin/CORS problem.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for browser DevTools Network debugging.
:width: 100%

Use the browser Network tab to separate “server/format” problems from “UI” problems.
```

### 4.2 Console tab checklist

In DevTools → Console:

- copy/paste the first relevant error (don’t paraphrase)
- look for:
  - “Failed to fetch”
  - “CORS blocked”
  - JSON parse errors (manifest schema mismatch)

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
- whether the viewer is displayed,
- web UI cache status + build id,
- a Python→frontend ping/pong roundtrip,
- a frontend “debug snapshot” (URL/origin/userAgent as seen by the iframe),
- and recent frontend console warnings forwarded to Python.

If you are filing a bug report, paste the report as JSON:

```python
import json
print(json.dumps(report, indent=2, default=str))
```

<!-- SCREENSHOT PLACEHOLDER
ID: notebook-debug-connection-report
Suggested filename: developer/notebook-debug-connection-report.png
Where it appears: Python Package → Developer Docs → Debugging playbook
Capture:
  - UI location: notebook output cell showing the printed JSON report
  - State prerequisites: viewer created (even if broken)
  - Action to reach state: run `print(json.dumps(viewer.debug_connection(), indent=2, default=str))`
Crop:
  - Include: the top of the report (viewer_url/server_url) and one failing section (e.g. viewer_index_probe_error)
  - Exclude: tokens (viewerToken) if present, private hostnames/paths if sensitive
Alt text:
  - Notebook output showing a structured JSON connectivity report produced by viewer.debug_connection.
Caption:
  - `viewer.debug_connection()` produces a single structured report that captures the most common notebook failure modes (proxy, cache, server reachability, hook roundtrip).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for viewer.debug_connection output.
:width: 100%

Notebook debugging: `viewer.debug_connection()` consolidates server reachability, cache status, and frontend roundtrip checks into one report.
```

---

## Common failure patterns (fast diagnosis)

### Pattern: “Viewer UI loads, but data is missing”

Export folder mode:
- manifests might reference files that don’t exist (partial export)
- safe-key mismatch (field names vs filenames)
- `.gz` mismatch (exported compressed, but viewer requests uncompressed or vice versa)

AnnData server mode:
- route for that path may not exist (wrong filename/ext)
- adapter could be raising internally (check server logs)

Start with:
- {doc}`08_export_format_spec_and_invariants`
- {doc}`07_prepare_export_pipeline_architecture`

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
- or provide a stable reverse proxy and set `CELLUCID_CLIENT_SERVER_URL`.

---

## What to include in a high-quality bug report

Minimum:

- mode (export folder vs AnnData vs notebook)
- `cellucid` version + Python version + OS
- reproduction steps (code or click-by-click)
- expected vs actual
- server URL + what `/_cellucid/health` returns
- browser console error(s)

Notebook:
- attach `viewer.debug_connection()` report.
