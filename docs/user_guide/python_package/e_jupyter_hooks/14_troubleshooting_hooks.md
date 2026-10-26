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
- `anndata_store_probes`, keyed by every exact listed dataset id (a `HEAD` and
  a one-byte `Range` read of the served AnnData store, which a readable store
  answers with `206` and `Content-Range`)
- ping/pong roundtrip (Python → iframe → Python)
- a frontend “debug snapshot” (iframe URL/origin/user agent)
- recent accepted-event counts grouped by exact event type

If you include this report when asking for help, debugging is usually 10× faster.

---

## Symptom: the viewer UI cannot be established or served

### Symptom
`show_anndata(...)` / `AnnDataViewer(...)` raises while establishing the
viewer's web build, `viewer.display()` raises while verifying it, or an iframe
that was already shown reloads to the plain text “Internal server error”.

### Likely causes (ordered)
1. No outbound network access to the configured web source (the constructor
   raises: every viewer start establishes the build again).
2. The source inventory or one declared asset failed exact validation (the
   constructor raises).
3. The selected generation directory is not writable (the constructor raises).
4. The generation was removed after the viewer started (for example
   `viewer.clear_web_cache()`): `display()` refuses it, and the server can no
   longer read a verified page to answer with (the iframe text above; the
   server log names the missing cache directory).

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
1. The viewer's channel (`/_cellucid/ws`) could not open: a proxy in between
   does not pass WebSockets, or refuses the page's origin. The viewer says so
   in a notice ("…cannot open its channel…").
2. You are interacting with a stale iframe (kernel restarted); it froze and
   says Python stopped.
3. A hook raised: the exception is in the log (`cellucid.viewer_channel`).
   The event was still recorded (`viewer.state`), and the other hooks ran.
4. The kernel's machine could not hold a very large selection. A selection
   costs 4 bytes a cell on the wire and about 36 bytes a cell as the
   `list[int]` hooks receive.

### How to confirm
1. In Python, run:
   ```python
   viewer.debug_connection()
   ```
   `channel_connected` must be `True`.
2. In browser devtools → Network, filter `WS`:
   - the `/_cellucid/ws` connection should be `101 Switching Protocols`;
   - its frames show each event as you select, hover or click.
3. Confirm you are using the right viewer object:
   - print `viewer.viewer_url` and compare with the iframe URL (if visible).

### Fix
- Connectivity/proxy issues:
  - follow {doc}`13_security_cors_origins_and_mixed_content`
- Stale iframe:
  - re-run the cell that creates the viewer
- Proxy without WebSockets:
  - enable WebSocket proxying on the proxy, or open the viewer through SSH
    port forwarding instead
- Selection too large for the kernel's memory:
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
1. The viewer **rejected** the command — the field, index or colour was not one
   it could apply.
2. Viewer not displayed: commands wait in Python until a viewer page
   connects, so nothing happens until one does.
3. The iframe is stale or not the one tied to this Python viewer.
4. Connectivity/proxy issues prevent the viewer channel (`/_cellucid/ws`)
   from opening.
5. The frontend and Python package did not come from the same exact current
   generation.

### How to confirm
- A rejected command prints to `sys.stderr` as soon as the viewer reports it:
  ```text
  [cellucid] viewer.set_color_by(...) did not take effect in the viewer: ...
  ```
  If you missed it, list the rejections after the fact:
  ```python
  viewer.debug_connection()["recent_events"]["command_errors"]
  ```
  A non-empty list means the command arrived and was refused — fix the
  argument, not the connection.
- If `viewer.debug_connection()["displayed"]` is False (notebooks): you
  didn’t display the viewer; `channel_connected` says whether a viewer page is
  connected now.
- Run:
  ```python
  viewer.wait_for_ready(timeout=60)
  viewer.debug_connection()
  ```

### Fix
1. If the command was rejected, correct the argument the reason names — most
   often a field ID the loaded dataset does not have.
2. Display the viewer:
   ```python
   viewer.display()
   ```
3. If the iframe is stale, stop that viewer:
   ```python
   viewer.stop()
   ```
   Then construct and display one current viewer.
4. Use {doc}`13_security_cors_origins_and_mixed_content` to configure the exact
   browser-reachable URL.

### Prevention
- Always `viewer.wait_for_ready()` before sending commands in scripts/notebooks.
- In a script, collect rejections and assert on them at a point you control,
  rather than reading a printed line:
  ```python
  rejected = []

  @viewer.on_command_error
  def record(event):
      rejected.append(event)

  ...
  assert not rejected, rejected
  ```

---

## Symptom: session bundle capture fails or hangs

### Symptom
`viewer.get_session_bundle()` times out or raises.

### Likely causes
1. The viewer's channel is not connected (connectivity/proxy).
2. The request timed out before the bundle arrived; the log says “Dropped a
   session bundle nobody is waiting for”.
3. The temporary directory has less free space than the bundle;
   `get_session_bundle()` raises “The session bundle could not be saved: …”.

### How to confirm
- Browser devtools → Network, filter `WS`: the `/_cellucid/ws` connection
  shows the `requestSessionBundle` command and the binary frame that answers it.
- In Python, inspect the debug report for `channel_connected` and recent
  `session_bundle` and `command_error` events.

### Fix
- Resolve any reported connectivity failure. An expired request is invalid;
  one subsequent `viewer.get_session_bundle()` call (with a longer
  `timeout=`) creates a new exact request ID.
- Reduce session size (fewer highlight groups / less state).
- Fix connectivity issues (proxy/mixed content).

---

## Symptom: “Port already in use” / the viewer port changed

### Symptom
An explicitly requested port is busy, or an automatically assigned viewer port
does not match a tunnel configured earlier.

### Likely causes
- Without `port=`, a notebook viewer takes 8765, or the next free port when
  another viewer (or program) holds it, and prints a note saying so (Google
  sign-in works only on 8765). A port you pass is used or refused, never
  replaced.
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
  - filter `WS` and select the `/_cellucid/ws` connection to see its frames

---

## Next steps

- Reference (schemas + endpoints): {doc}`16_reference`
