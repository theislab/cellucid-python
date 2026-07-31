# Compatibility matrix (must be explicit)

**Audience:** notebook users, IT-restricted environments, power users  
**Time:** 10–15 minutes  
**Goal:** know what works where (and what to do when it doesn’t)

This page is intentionally explicit because many “it doesn’t work” reports are actually **environment + browser security** issues (HTTPS notebooks, remote kernels, blocked network, etc.).

## At-a-glance matrix

### Environment support (notebook embedding + hooks)

Legend:
- `OK`: works out of the box in common setups
- `Setup`: works with setup / common caveats
- `No`: not supported (use browser workflow instead)

| Environment | Viewer embeds (iframe) | Python → frontend commands | Frontend → Python hooks/events | Notes |
|---|---:|---:|---:|---|
| Classic Jupyter Notebook (local) | OK | OK | OK | Usually simplest; notebook is often `http://localhost`. |
| JupyterLab (local) | OK | OK | OK | If JupyterLab is served from **HTTPS**, expose the Cellucid port through an HTTPS route and pass that route explicitly. |
| VSCode notebooks (local) | OK | OK | OK | If VSCode runs in Remote/SSH containers, treat it like “remote kernel”. |
| Google Colab | Setup | OK | OK | Obtain Colab's HTTPS proxy base, pass it as `client_server_url=`, and declare its host as `allowed_hosts=[...]`; Cellucid does not call the proxy API. |
| JupyterHub / remote notebooks (HTTPS) | Setup | OK | OK | Configure an HTTPS route for the selected port, pass its exact base as `client_server_url=`, and declare the route's host as `allowed_hosts=[...]`. |

### Network requirements

Cellucid has two different network surfaces:
- the **viewer UI** (HTML/JS/CSS and other declared static assets) — fetched
  from the configured web source and established as one verified generation
  when the Python server starts.
- your **data** — served from your local/remote Python server (or loaded via browser file picker).

| Feature | Source access required? | Notes |
|---|---:|---|
| `cellucid serve exports/...` with its viewer | Yes, at startup | The server downloads, stages, byte-verifies, and atomically publishes the complete source generation before binding its HTTP server. |
| `show(...)` / `show_anndata(...)` | Yes, at startup | The notebook viewer uses the same establishment contract. |
| `cellucid serve ... --no-web-ui` | No | Explicit data-endpoint-only mode; no viewer is served. |
| Hosted web app (`cellucid.com`) | Yes | The browser loads the website directly. |
| `ensure_web_ui_cached(force=True)` | Yes | Establishes the complete source generation now. |
| `ensure_web_ui_cached(force=False)` | No | Verifies an existing generation only; it does not alter normal viewer startup. |
| Hooks/events after a viewer has started | No additional source request | Messages use `/_cellucid/events` on the Python server. |

## Practical notes per environment

### Classic Jupyter Notebook / JupyterLab (local)

What usually works:
- `show_anndata(adata, dataset_name="Example", dataset_id="example")` embeds an iframe in the output cell.
- Hooks and commands work bidirectionally.

Common caveat:
- If your notebook page is served from **HTTPS** (rare locally but common on managed systems), the browser can block an `http://127.0.0.1:<port>` iframe (“mixed content”).

Fixes:
- Configure an HTTPS route for the Cellucid port, pass its exact
  browser-reachable base as `client_server_url=`, and pass the route's host
  name as `allowed_hosts=["hub.example.org"]`. The Python package already
  depends on `jupyter-server-proxy`; installation alone does not select or
  configure a route. The second argument is required because the proxy forwards
  the browser's `Host` header verbatim and Cellucid refuses an undeclared
  authority with `421 Misdirected Request`
  (see {doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`).
- Or skip embedding and use the browser workflow (`cellucid serve ...`) and open the URL directly.

### VSCode notebooks

Local VSCode works like local Jupyter: embedding + hooks typically work.

Remote VSCode (SSH / containers) is different:
- Your **kernel** runs remotely.
- Your **browser/UI** is local.
- A server bound to `127.0.0.1` on the remote machine is **not** reachable from your local browser unless you forward ports.

Recommended pattern (SSH/remote):
1) Forward the Cellucid server port to your machine (VSCode “Ports” tab or SSH `-L`).
2) Pass the exact URL the **browser** should use when you create the viewer:
   ```python
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       adata,
       port=8765,
       dataset_name="Example",
       dataset_id="example",
       client_server_url="http://127.0.0.1:8765",
   )
   ```

### Google Colab

Colab runs the kernel on a remote VM. The caller must obtain Colab's HTTPS
proxy base for one fixed port and pass it explicitly:

Typical pattern:
1) Install in the Colab runtime:
   ```ipython
   !pip -q install cellucid
   ```
2) Load data in the runtime filesystem (uploaded, mounted Drive, etc.)
3) Create the viewer on that same fixed port:
   ```python
   from urllib.parse import urlparse

   from cellucid import AnnDataViewer
   from google.colab.output import eval_js

   port = 8765
   browser_base = eval_js(f"google.colab.kernel.proxyPort({port})")
   viewer = AnnDataViewer(
       adata,
       port=port,
       dataset_name="Example",
       dataset_id="example",
       client_server_url=browser_base,
       allowed_hosts=[urlparse(browser_base).hostname],
   )
   ```

   Colab's proxy also forwards the browser's `Host`, so its host name must be
   declared exactly as for any other reverse proxy.

Cellucid does not call `google.colab.kernel.proxyPort(...)` itself.

If it fails:
- Run `viewer.debug_connection()` and look at:
  - `client_server_url`
  - `allowed_hosts`
  - `server_health`
  - `dataset_identity_probes`
  - `frontend_roundtrip` and `frontend_debug_snapshot`

## Security / browser constraints (what’s going on under the hood)

The two issues that matter most:

### 1) HTTPS → HTTP mixed content (iframe blocking)

If your notebook UI is served from `https://...`, the browser may refuse to load `http://127.0.0.1:<port>` inside an iframe.

Cellucid does not probe, infer, or select a proxy. It uses the exact
`client_server_url=` supplied by the caller, or the bound loopback server URL
when the argument is omitted. A blank iframe, mixed-content error, or failed
network request means the caller must expose the selected port and pass its
browser-reachable HTTP(S) base.

### 2) Exact UI generation availability

Cellucid’s Python server:

1. fetches `cellucid-web-assets.json` from the configured source,
2. downloads every object in that exact inventory into a staging directory,
3. verifies the complete file set, byte lengths, SHA-256 hashes, MIME types,
   and build identity, and
4. atomically publishes the verified generation before starting the server.

If the source or any declared object is unavailable or invalid, startup stops
and the prior cache remains untouched. Supply a reachable
`--web-source-url` and a writable `--web-cache-dir`; an older cache is never
served as a substitute.

## Edge cases (things that look like “bugs”)

- **Ad blockers / strict corporate CSP**: can break the hosted web UI fetch or iframe behavior.
- **Port collisions**: an explicitly requested busy port fails. In Python,
  `port=0` asks the operating system to allocate one port; use the actual URL
  reported by the viewer/server object.
- **Multiple viewers in one notebook**: hooks are routed by `viewerId`. If you re-run cells many times, stop old viewers (`viewer.stop()`) or restart the kernel.
- **Remote access**: binding `--host 0.0.0.0` exposes your server to the network; do this only if you understand the security implications.

## Troubleshooting (compatibility)

### Symptom: viewer iframe is blank / white

**Likely causes**
- Browser blocked the iframe (mixed content).
- The exact viewer generation could not be established from its configured source.
- The browser can’t reach the data server (remote kernel without port forwarding).

**How to confirm**
- In Python:
  ```python
  report = viewer.debug_connection()
  report
  ```
- In the browser: open DevTools → Console.

**Fix**
- Configure an HTTPS proxy route or use port forwarding, then pass that exact
  browser-reachable base as `client_server_url=`. Behind a proxy, also declare
  the proxy's host name as `allowed_hosts=[...]`, otherwise every request is
  refused with `421 Misdirected Request`.
- Confirm source access with `ensure_web_ui_cached(force=True)` and correct the
  reported source, inventory, object, or cache-directory failure. A previous
  generation is intentionally not substituted.

---

### Symptom: hooks don’t fire (`on_selection`, `on_click`, …)

**Likely causes**
- Viewer isn’t fully loaded yet (wait for `on_ready`).
- Old viewer is still registered (multiple viewers, stale kernel state).
- Browser can’t reach `/_cellucid/events` on the server.

**How to confirm**
- Add a raw message logger:
  ```python
  @viewer.on_message
  def _log(ev):
      print(ev)
  ```
- Check `viewer.debug_connection()` and confirm `server_health` is OK.

**Fix**
- Restart the kernel, or call `viewer.stop()` on old viewers.
- Ensure you’re not blocked by remote/HTTPS constraints.

## Next steps

- End-to-end quick starts: {doc}`04_quick_start_3_levels`
- Bidirectional hooks deep dive: {doc}`../e_jupyter_hooks/index`
