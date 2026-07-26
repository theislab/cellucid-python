# Debugging mental model (where to look first)

**Audience:** everyone (beginners: follow the checklist; experts: jump to the relevant symptom)  
**Time:** 10–30 minutes  
**Goal:** have a systematic, low-guesswork workflow for debugging viewer/server/hooks/session issues.

---

## The debugging rule: separate “UI problems” from “server problems”

Most Cellucid failures fall into one of these buckets:

1) **The exact viewer generation cannot be established or served**
2) **The browser can’t reach the Python server** (port/network/tunnel/proxy)
3) **The dataset can’t be interpreted** (missing embeddings/manifests)
4) **Hooks/events aren’t being delivered** (routing/viewerId mismatch)
5) **Session capture/apply fails** (not ready, requestId gating, exact dataset fingerprint mismatch)

The fastest way to debug is to test these buckets in order.

---

## Step 0 (do this first): `viewer.debug_connection()`

If you have a `viewer` object, run:

```python
report = viewer.debug_connection()
report
```

This returns a structured dict that checks:
- server probes (`/_cellucid/health`, `/_cellucid/info`, `/_cellucid/datasets`)
- whether the verified viewer UI (`/index.html`) is being served
- a Python → frontend roundtrip (ping/pong)
- a frontend debug snapshot (iframe URL/origin/user agent)
- recent frontend console warnings/errors forwarded to Python

### How to read the report (high signal fields)

- `server_health` / `server_health_error`
  - If this fails: the server isn’t reachable (port/network problem).
- `viewer_index_probe` / `viewer_index_probe_error`
  - If this fails: the active verified generation is not being served.
- `frontend_roundtrip`
  - If `ok=False`: Python→frontend messaging is broken or the iframe isn’t alive.
- `frontend_debug_snapshot`
  - Confirms the live iframe URL, origin, server URL, connection state, parent
    origin, and browser identification.
- `recent_events`
  - Confirms which accepted event types reached Python.
- `client_server_url`
  - This must be the exact HTTP(S) base that the browser can reach.

```{tip}
When asking for help, paste the `report` (redacting private paths/dataset names) instead of describing symptoms loosely.
```

---

## Step 1: confirm the server is alive (no browser required)

In Python, check:

```python
print(viewer.server_url)
```

Then probe the health endpoint:

```python
import urllib.request, json
json.loads(urllib.request.urlopen(viewer.server_url + "/_cellucid/health").read())
```

Expected:
- a JSON payload like `{"status": "ok", ...}`

If this fails:
- the server may not be running,
- the port may be blocked/in use,
- you might be on a remote kernel whose loopback isn’t reachable from your browser (see Step 3).

---

## Step 2: confirm the viewer UI is being served

Cellucid notebook/server mode serves one exact verified web generation from the
Python server.

Check:
- `viewer.server_url + "/index.html"`

Viewer construction establishes that generation before binding the server and
raises on failure. If the published generation is removed while a server is
running, asset requests return an explicit service error. See
{doc}`06_privacy_security_and_offline_vs_online`.

---

## Step 3: confirm the browser can reach the server (remote notebooks)

The most common “it works locally but not on JupyterHub/HPC” issue:
- the browser cannot reach the kernel’s loopback server (`127.0.0.1:<port>`), or
- the notebook page is HTTPS and blocks HTTP loopback (mixed content).

Fix patterns:

### Option A (recommended): `jupyter-server-proxy`

If installed and enabled, the notebook can proxy the viewer through the
notebook server origin. Obtain the proxy’s browser base and pass it explicitly
as `client_server_url=`; Cellucid does not detect or select it.

### Option B: SSH port forwarding (robust, works everywhere)

If your kernel runs on a remote machine:

1) Start the viewer on a fixed port:

   ```python
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       "data.h5ad",
       port=8765,
       dataset_name="Example",
       dataset_id="example",
       client_server_url="http://127.0.0.1:8765",
   )
   print(viewer.viewer_url)
   ```

2) On your laptop:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now the browser hits your laptop’s localhost port, which tunnels to the remote kernel.

See the detailed guide: {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`

---

## Step 4: debug hooks/events (frontend → Python)

### Minimal diagnostic hook

```python
@viewer.on_message
def debug(event):
    print(event)
```

If you see nothing when interacting in the UI:
- the browser is not POSTing to `/_cellucid/events`, or
- the `viewerId` in the browser doesn’t match the viewer object that registered the hooks.

### Browser-side confirmation

Open browser devtools → Network and look for:
- `POST /_cellucid/events`

If there are no POSTs:
- the UI may not be wired to the server URL you think it is,
- or the server URL is unreachable.

If POSTs exist but Python sees nothing:
- the viewerId is mismatched (old iframe / old kernel / multiple viewers).

Mitigation:
- recreate the viewer,
- avoid copying stale viewer URLs across sessions.

---

## Step 5: debug Python → frontend commands

If you send commands like:

```python
viewer.highlight_cells([1, 2, 3])
viewer.set_color_by("cell_type")
```

and nothing happens:
- the iframe may not be alive,
- the notebook output may have been cleared/collapsed,
- the viewer may not be ready yet.

Use `viewer.debug_connection()` and check:
- `frontend_roundtrip`

If `frontend_roundtrip.ok` is false:
- call `viewer.display()` again,
- or recreate the viewer.

---

## Step 6: debug session capture (`get_session_bundle`)

Common failure modes:

### Timeout waiting for bundle

Likely causes:
- viewer never reached ready,
- iframe not displayed,
- browser can’t POST to `/_cellucid/session_bundle` due to proxy/mixed-content issues.

How to confirm:
- `viewer.debug_connection()` (look at `state.ready`, `frontend_roundtrip`,
  `frontend_debug_snapshot`, and `recent_events`)
- browser network tab for
  `POST /_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`

Fix:
- call `viewer.wait_for_ready(timeout=60)` first,
- ensure the viewer is displayed,
- fix remote notebook connectivity (Step 3).

### “No pending session bundle request”

Meaning:
- an upload arrived without a matching requestId (stale iframe, TTL expired, wrong viewerId).

Fix:
- request again from the current viewer instance.

For more: {doc}`05_sessions_to_anndata_bridge`

---

## Step 7: debug “the dataset looks wrong” (identity + ordering)

If selection indices don’t match the cells you expect in Python, the usual culprit is:
- dataset row order changed,
- you are analyzing a different AnnData object than the one the viewer is showing.

Checks:
- confirm which AnnData was passed to `show_anndata(...)`,
- if applying sessions, inspect `bundle.dataset_fingerprint` and `expected_dataset_id`.

See: {doc}`04_dataset_identity_and_reproducibility`

---

## Troubleshooting index (symptom → fix)

### Symptom: “Blank viewer area”

Likely causes:
- the configured source generation could not be fetched, verified, or
  published during viewer construction,
- mixed-content blocked (HTTPS notebook),
- server not reachable.

Fix:
- run `viewer.debug_connection()` and follow the first failing probe.

### Symptom: “Viewer loads but hooks don’t work”

Likely causes:
- no POSTs to `/_cellucid/events`,
- viewerId mismatch.

Fix:
- check browser Network tab,
- recreate the viewer.

### Symptom: “`apply_cellucid_session_to_anndata(...)` raised an error”

Likely causes:
- the required `expected_dataset_id`, cell count, or gene count does not match
  the bundle fingerprint,
- an output column already exists,
- the session does not contain a requested highlight or user-defined-field
  chunk.

Fix:
- read the raised exception; fingerprint mismatches and column collisions are
  rejected before mutation,
- pass the exact dataset ID for the AnnData object,
- inspect `bundle.list_chunk_ids()` before selecting which current chunks to
  apply.

---

## Next steps

- Privacy/offline model: {doc}`06_privacy_security_and_offline_vs_online`
- Hooks docs: {doc}`../e_jupyter_hooks/index`
