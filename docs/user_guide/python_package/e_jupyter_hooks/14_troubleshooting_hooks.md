# Troubleshooting (hooks)

This page is a symptom-based guide for when notebook embedding + hooks don’t behave as expected.

Rule of thumb: **most failures are connectivity/proxy issues**, not Python hook code.

## First step (always): collect a debug report

Run:

```python
report = viewer.debug_connection()
report
```

This report includes:
- server probes (`/_cellucid/health`, `/_cellucid/info`, `/_cellucid/datasets`)
- `dataset_identity_probes`, keyed by every exact declared dataset id/path
- ping/pong roundtrip (Python → iframe → Python)
- a frontend “debug snapshot” (iframe URL/origin/user agent)
- recent accepted-event counts grouped by exact event type

If you include this report when asking for help, debugging is usually 10× faster.

---

## Symptom: iframe shows “Cellucid viewer UI could not be loaded”

### Symptom
The iframe loads an error page saying the viewer UI could not be loaded.

### Likely causes (ordered)
1. No outbound network access to the configured web source.
2. The source inventory or one declared asset failed exact validation.
3. The selected generation directory is not writable.

### How to confirm
- In the debug report, check:
  - `report["web_ui"]["cache"]`
  - `report["viewer_index_probe_error"]` (browser-facing index probe failure)
- Try opening in a browser:
  - `<viewer.server_url>/index.html`

### Fix
- Ensure the Python runtime can reach the configured `web_source_url`.
- Correct the exact inventory, object, MIME, length, hash, or build error shown.
- Pass a writable `web_cache_dir` when the default temporary location is not
  suitable.

### Prevention
- Before a shared session, exercise the exact startup source with
  `viewer.ensure_web_ui_cached(force=True)`.
- Remember that the real viewer startup establishes the source generation
  again; the check does not provide an offline mode.

---

## Symptom: remote or HTTPS notebook iframe is blank or blocked

### Symptom
The iframe is blank, browser developer tools report mixed content, or the
browser cannot reach the kernel-side loopback URL.

### Likely causes
1. Your notebook frontend is served over HTTPS (or is remote), so direct `http://127.0.0.1:<port>` is blocked/unreachable.
2. No browser-reachable proxy or tunnel has been configured for the selected
   port.

### How to confirm
- In browser devtools console/network, look for mixed-content blocking.
- In the debug report, inspect `client_server_url`, `viewer_url`,
  `viewer_index_probe`, and any corresponding error keys.

### Fix (recommended order)
1. Configure an HTTPS proxy route for the fixed Cellucid port, or use SSH port
   forwarding when the browser is local.
2. Pass that exact browser-reachable HTTP(S) base as `client_server_url=...`.

### Prevention
- Validate the proxy/tunnel URL before a shared remote-notebook workflow.
- The Cellucid package already depends on `jupyter-server-proxy`, but Cellucid
  does not detect, configure, or select a proxy URL.

---

## Symptom: viewer loads but **no events arrive** (`on_selection` never fires)

### Symptom
The viewer UI is interactive, but Python callbacks do not fire.

### Likely causes (ordered)
1. The iframe cannot reach the server’s `/_cellucid/events` endpoint (proxy/mixed-content/firewall).
2. Requests to `/_cellucid/events` are blocked by an ad blocker / corporate proxy.
3. You are interacting with a stale iframe (kernel restarted; viewerId no longer registered).
4. The selection payload is too large (hits the 1MB request limit).

### How to confirm
1. In Python, run:
   ```python
   viewer.debug_connection()
   ```
2. In browser devtools → Network:
   - look for requests to `/_cellucid/events`
   - confirm status codes (should be 200 JSON)
3. Confirm you are using the right viewer object:
   - print `viewer.viewer_url` and compare with the iframe URL (if visible).

### Fix
- Connectivity/proxy issues:
  - follow {doc}`13_security_cors_origins_and_mixed_content`
- Stale iframe:
  - re-run the cell that creates the viewer
- Payload too large:
  - avoid selecting huge fractions of the dataset in one step
  - prefer session bundle capture for exporting large highlight membership sets

### Prevention
- For large datasets, avoid “select everything”; build highlight groups in smaller steps.
- Keep notebooks reproducible: re-run the viewer creation cell after kernel restart.

---

## Symptom: events arrive, but indices don’t match my AnnData

### Symptom
Python receives indices, but they refer to “the wrong cells” when you index `adata[cells]`.

### Likely causes
1. You reordered/subset `adata` after creating the viewer.
2. You are using a different `adata` object than the one you passed to `show_anndata`.

### How to confirm
- Print basic sanity checks:
  ```python
  event = viewer.wait_for_event("selection", timeout=None)
  cells = event["cells"]
  print(min(cells), max(cells), len(cells))
  print("adata.n_obs:", adata.n_obs)
  ```
- If `max(cells) >= adata.n_obs`, you are definitely mismatched.

### Fix
- Only interpret indices against the exact dataset the viewer is serving.
- If you must reorder/subset, create a new viewer from the modified AnnData.

### Prevention
- Treat the viewer’s dataset as immutable while you’re using hooks.

---

## Symptom: Python commands don’t affect the UI (highlight/color/reset do nothing)

### Symptom
Calling `viewer.highlight_cells(...)` / `viewer.set_color_by(...)` / `viewer.reset_view()` has no visible effect.

### Likely causes (ordered)
1. Viewer not displayed (commands are sent via notebook JS injection).
2. The iframe is stale or not the one tied to this Python viewer.
3. Connectivity/proxy issues prevent the iframe from fully initializing Jupyter mode.
4. The frontend and Python package did not come from the same exact current
   generation.

### How to confirm
- If `viewer._displayed` is False (notebooks): you didn’t display the viewer.
- Run:
  ```python
  viewer.wait_for_ready(timeout=60)
  viewer.debug_connection()
  ```

### Fix
1. Display the viewer:
   ```python
   viewer.display()
   ```
2. If the iframe is stale, stop that viewer:
   ```python
   viewer.stop()
   ```
   Then construct and display one current viewer.
3. Use {doc}`13_security_cors_origins_and_mixed_content` to configure the exact
   browser-reachable URL.

### Prevention
- Always `viewer.wait_for_ready()` before sending commands in scripts/notebooks.

---

## Symptom: session bundle capture fails or hangs

### Symptom
`viewer.get_session_bundle()` times out or raises.

### Likely causes
1. Iframe can’t reach `/_cellucid/session_bundle` (connectivity/proxy).
2. Pending request expired; server responds “No pending session bundle request”.
3. Bundle exceeds size cap (512MB).

### How to confirm
- Browser devtools → Network:
  - look for
    `/_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`
  - check status code and response text
- In Python, inspect the debug report for recent `session_bundle` events.

### Fix
- Resolve any reported connectivity failure. An expired request is invalid;
  one subsequent `viewer.get_session_bundle()` call creates a new exact
  request ID.
- Reduce session size (fewer highlight groups / less state).
- Fix connectivity issues (proxy/mixed content).

---

## Symptom: “Port already in use” / the viewer port changed

### Symptom
An explicitly requested port is busy, or an automatically assigned viewer port
does not match a tunnel configured earlier.

### Likely causes
- You created multiple viewers and didn’t stop them.
- You re-ran cells repeatedly without cleanup.

### Fix
- Stop the viewer when done:
  ```python
  viewer.stop()
  ```
- For remote/HPC workflows, choose a fixed port:
  ```python
  from cellucid import AnnDataViewer

  viewer = AnnDataViewer(
      "data.h5ad",
      port=8765,
      dataset_name="My study",
      dataset_id="my-study-v1",
  )
  ```

### Prevention
- Use a `try/finally` pattern to stop viewers.

---

## Debugging checklist (copy/paste)

```python
print("server_url:", viewer.server_url)
print("viewer_url:", viewer.viewer_url)

viewer.wait_for_ready(timeout=60)
report = viewer.debug_connection()
report
```

In a browser:
- open `<server_url>/_cellucid/health`
- open devtools → Network:
  - filter for `/_cellucid/events` and `/_cellucid/session_bundle`

---

## Next steps

- Reference (schemas + env vars): {doc}`16_reference`
