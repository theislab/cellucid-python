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
- ping/pong roundtrip (Python → iframe → Python)
- a frontend “debug snapshot” (iframe URL/origin/user agent)
- recent frontend warnings/errors forwarded to Python (`console` events)

If you include this report when asking for help, debugging is usually 10× faster.

---

## Symptom: iframe shows “Cellucid viewer UI could not be loaded”

### Symptom
The iframe loads an error page saying the viewer UI could not be loaded.

### Likely causes (ordered)
1. No outbound network access to `https://www.cellucid.com` (first run).
2. You are offline and the UI cache is empty.
3. Cache directory is not writable/persistent.

### How to confirm
- In the debug report, check:
  - `report["web_ui"]["cache"]`
  - `report["viewer_index_probe_error"]` (proxy fetch failure)
- Try opening in a browser:
  - `<viewer.server_url>/index.html`

### Fix
- If you can go online once:
  - re-run the cell; Cellucid will prefetch UI assets with a progress bar.
- If you need a persistent cache directory:
  - set `CELLUCID_WEB_PROXY_CACHE_DIR` to a writable location.
- If you suspect a bad cache:
  - `viewer.clear_web_cache()` and retry.

### Prevention
- On shared/HPC environments, set `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent path.
- Prefetch once in a “setup notebook”:
  - `viewer.ensure_web_ui_cached()`

<!-- SCREENSHOT PLACEHOLDER
ID: jupyter-hooks-ui-unavailable-error
Suggested filename: jupyter_hooks/10_ui-unavailable-error-page.png
Where it appears: Python Package → Jupyter Hooks → Troubleshooting → UI unavailable
Capture:
  - UI location: iframe inside notebook output
  - State prerequisites: run a viewer in an offline environment with no cached UI
  - Action to reach state: clear cache, disable internet, run show/show_anndata
Crop:
  - Include: the full error message text (so users can match it)
  - Exclude: private paths
Alt text:
  - Embedded viewer showing an error page stating the Cellucid viewer UI could not be loaded.
Caption:
  - When the notebook cannot download or find cached viewer UI assets, the iframe shows a clear error page with next steps.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the “viewer UI unavailable” error page in the notebook iframe.
:width: 100%

If UI assets cannot be fetched or found in the cache, the iframe shows an error page with next steps.
```

---

## Symptom: iframe shows “notebook proxy required”

### Symptom
The iframe shows a page saying a notebook proxy is required (or similar messaging).

### Likely causes
1. Your notebook frontend is served over HTTPS (or is remote), so direct `http://127.0.0.1:<port>` is blocked/unreachable.
2. `jupyter-server-proxy` is not installed/enabled.

### How to confirm
- In browser devtools console/network, look for mixed-content blocking.
- In the debug report, check `report["jupyter_server_proxy"]`.

### Fix (recommended order)
1. Install/enable `jupyter-server-proxy` on the Jupyter server.
2. If your kernel is remote and your browser is local: use SSH port forwarding.
3. Advanced: set `CELLUCID_CLIENT_SERVER_URL` to a browser-reachable HTTPS URL for the server.

### Prevention
- On JupyterHub, treat `jupyter-server-proxy` as part of the standard environment.

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
4. You are offline with a cached UI build that doesn’t support the command (rare, but possible).

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
2. Re-run the viewer creation cell.
3. Clear UI cache and retry:
   ```python
   viewer.clear_web_cache()
   ```
4. Use {doc}`13_security_cors_origins_and_mixed_content` to fix proxy/mixed-content failures.

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
  - look for `/_cellucid/session_bundle?viewerId=...&requestId=...`
  - check status code and response text
- In Python, inspect the debug report for recent `session_bundle` events.

### Fix
- Retry `viewer.get_session_bundle()` (generates a fresh requestId).
- Reduce session size (fewer highlight groups / less state).
- Fix connectivity issues (proxy/mixed content).

---

## Symptom: “Port already in use” / ports keep incrementing

### Symptom
You see the viewer using ports like 8765, 8766, 8767… and remote workflows break.

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
  viewer = show_anndata("data.h5ad", port=8765)
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
- Screenshot/diagram checklist for docs: {doc}`15_screenshots_and_diagrams`

