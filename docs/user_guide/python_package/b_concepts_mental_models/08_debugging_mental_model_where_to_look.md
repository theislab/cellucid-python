# Debugging mental model (where to look first)

**Audience:** everyone (beginners: follow the checklist; experts: jump to the relevant symptom)  
**Time:** 10–30 minutes  
**Goal:** have a systematic, low-guesswork workflow for debugging viewer/server/hooks/session issues.

---

## The debugging rule: separate “UI problems” from “server problems”

Most Cellucid failures fall into one of these buckets:

1) **The browser can’t load the viewer UI** (offline/cache/mixed-content/proxy)
2) **The browser can’t reach the Python server** (port/network/tunnel/proxy)
3) **The dataset can’t be interpreted** (missing embeddings/manifests)
4) **Hooks/events aren’t being delivered** (routing/viewerId mismatch)
5) **Session capture/apply fails** (not ready, requestId gating, mismatch policy)

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
- whether the viewer UI (`/index.html`) is being served (proxy + cache)
- a Python → frontend roundtrip (ping/pong)
- a frontend debug snapshot (iframe URL/origin/user agent)
- recent frontend console warnings/errors forwarded to Python

### How to read the report (high signal fields)

- `server_health` / `server_health_error`
  - If this fails: the server isn’t reachable (port/network problem).
- `viewer_index_probe` / `viewer_index_probe_error`
  - If this fails: the UI assets aren’t being served (offline/cache/network).
- `frontend_roundtrip`
  - If `ok=False`: Python→frontend messaging is broken or the iframe isn’t alive.
- `frontend_console`
  - If populated: often contains the real error (missing assets, fetch failures, blocked requests).
- `jupyter_server_proxy.installed`
  - If false in a remote/HTTPS notebook: mixed-content issues are likely.

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

Cellucid notebook/server mode serves the UI from the Python server (hosted-asset proxy + cache).

Check:
- `viewer.server_url + "/index.html"`

If the UI cannot be served and no cached copy exists, you’ll see an error page explaining how to populate the cache (see {doc}`06_privacy_security_and_offline_vs_online`).

---

## Step 3: confirm the browser can reach the server (remote notebooks)

The most common “it works locally but not on JupyterHub/HPC” issue:
- the browser cannot reach the kernel’s loopback server (`127.0.0.1:<port>`), or
- the notebook page is HTTPS and blocks HTTP loopback (mixed content).

Fix patterns:

### Option A (recommended): `jupyter-server-proxy`

If installed, the notebook can proxy the viewer through the notebook server origin.

The embed code will try this automatically when it detects HTTPS/remote origins.

### Option B: SSH port forwarding (robust, works everywhere)

If your kernel runs on a remote machine:

1) Start the viewer on a fixed port:

   ```python
   viewer = show_anndata("data.h5ad", port=8765)
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
- `viewer.debug_connection()` (look at `state.ready`, `frontend_roundtrip`, and `frontend_console`)
- browser network tab for `POST /_cellucid/session_bundle?...`

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

## Screenshot placeholder: browser devtools network POSTs to `/_cellucid/events`

This screenshot is extremely helpful for non-developers when debugging hooks.

<!-- SCREENSHOT PLACEHOLDER
ID: debugging-network-events-post
Suggested filename: data_loading/devtools-network-events-post.png
Where it appears: Python Package → Concepts → Debugging → Hooks/events
Capture:
  - UI location: browser devtools (Network tab)
  - State prerequisites: a viewer running; interact to trigger selection/hover
  - Action to reach state: open devtools, filter for `_cellucid/events`, do a lasso selection
Crop:
  - Include: the POST request row(s) + the request URL
  - Exclude: unrelated tabs/bookmarks, personal profile icons
Redact:
  - Remove: dataset identifiers if shown
Annotations:
  - Callouts: #1 POST /_cellucid/events request, #2 status code (200 vs failing)
Alt text:
  - Browser developer tools showing POST requests to the Cellucid events endpoint.
Caption:
  - Hooks are delivered as HTTP POSTs from the viewer to the local Python server; if these requests are missing or failing, Python callbacks cannot fire.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for browser devtools showing POST requests to /_cellucid/events.
:width: 100%

Hooks are delivered as HTTP POSTs from the viewer to the local Python server; if these requests are missing or failing, Python callbacks cannot fire.
```

---

## Troubleshooting index (symptom → fix)

### Symptom: “Blank viewer area”

Likely causes:
- UI assets not cached and network blocked,
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

### Symptom: “`apply_session_to_anndata` did nothing”

Likely causes:
- mismatch policy skipped dataset-dependent chunks,
- session doesn’t contain highlights/fields.

Fix:
- use `return_summary=True`,
- inspect `bundle.list_chunk_ids()`.

---

## Next steps

- Privacy/offline model: {doc}`06_privacy_security_and_offline_vs_online`
- Hooks docs: {doc}`../e_jupyter_hooks/index`
