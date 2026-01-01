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
| JupyterLab (local) | OK | OK | OK | If your JupyterLab is served from **HTTPS**, you may need `jupyter-server-proxy`. |
| VSCode notebooks (local) | OK | OK | OK | If VSCode runs in Remote/SSH containers, treat it like “remote kernel”. |
| Google Colab | Setup | OK | OK | Uses Colab port proxy; first run downloads viewer UI assets. |
| JupyterHub / remote notebooks (HTTPS) | Setup | OK | OK | Mixed-content is common; `jupyter-server-proxy` is strongly recommended. |

### Offline vs online requirements (what “offline” really means)

Cellucid has two different “web things”:
- the **viewer UI** (HTML/JS/CSS) — by default fetched from `https://www.cellucid.com` and cached locally when you use the Python server embed/proxy.
- your **data** — served from your local/remote Python server (or loaded via browser file picker).

| Feature | Works offline (after caching UI)? | Requires reaching `cellucid.com` right now? | Notes |
|---|---:|---:|---|
| `cellucid serve exports/...` (browser) | Yes | No | UI is served from your server *if it has a cached copy* of the web assets. |
| `show(...)` / `show_anndata(...)` (notebook) | Yes | No | Same as above. First-ever run may require internet to fetch UI assets. |
| Open the hosted web app in your browser (`cellucid.com`) | No | Yes | This is a normal website; you need internet to load it. |
| Prefetch web UI cache (`ensure_web_ui_cached`) | No | Yes | Prefetch is how you make later offline use smoother. |
| Hooks/events (selection/hover/click) | Yes | No | Events go to `/_cellucid/events` on your Python server. |

```{tip}
If you need an offline demo/teaching setup, run the UI prefetch once while online (see {doc}`02_installation`), set a persistent cache directory, then you can use `cellucid serve ...` without internet.
```

## Practical notes per environment

### Classic Jupyter Notebook / JupyterLab (local)

What usually works:
- `show_anndata(adata)` embeds an iframe in the output cell.
- Hooks and commands work bidirectionally.

Common caveat:
- If your notebook page is served from **HTTPS** (rare locally but common on managed systems), the browser can block an `http://127.0.0.1:<port>` iframe (“mixed content”).

Fixes:
- Install `jupyter-server-proxy` (recommended) and restart the notebook server:
  ```bash
  pip install jupyter-server-proxy
  ```
- Or skip embedding and use the browser workflow (`cellucid serve ...`) and open the URL directly.

### VSCode notebooks

Local VSCode works like local Jupyter: embedding + hooks typically work.

Remote VSCode (SSH / containers) is different:
- Your **kernel** runs remotely.
- Your **browser/UI** is local.
- A server bound to `127.0.0.1` on the remote machine is **not** reachable from your local browser unless you forward ports.

Recommended pattern (SSH/remote):
1) Forward the Cellucid server port to your machine (VSCode “Ports” tab or SSH `-L`).
2) Tell Cellucid what URL the **browser** should use:
   ```python
   import os
   os.environ["CELLUCID_CLIENT_SERVER_URL"] = "http://127.0.0.1:8765"  # your *local* forwarded port
   ```
3) Then create the viewer.

### Google Colab

Colab runs the kernel on a remote VM. Cellucid uses Colab’s HTTPS port proxy so the browser can reach the server.

Typical pattern:
1) Install in the Colab runtime:
   ```ipython
   !pip -q install cellucid
   ```
2) Load data in the runtime filesystem (uploaded, mounted Drive, etc.)
3) Call `show_anndata(...)`

If it fails:
- Run `viewer.debug_connection()` and look at:
  - `client_server_url`
  - `server_health`
  - any `frontend_console` messages

## Security / browser constraints (what’s going on under the hood)

The two issues that matter most:

### 1) HTTPS → HTTP mixed content (iframe blocking)

If your notebook UI is served from `https://...`, the browser may refuse to load `http://127.0.0.1:<port>` inside an iframe.

Cellucid’s embed logic tries to avoid this by:
- using a notebook-provided proxy URL when `jupyter-server-proxy` is available, or
- using environment-specific HTTPS proxies (e.g., Colab).

If you see an iframe message about needing a notebook proxy, that’s what it means.

<!-- SCREENSHOT PLACEHOLDER
ID: python-compat-mixed-content-proxy-required
Where it appears: Compatibility matrix → Security → mixed content
Capture:
  - UI location: embedded Cellucid iframe inside a notebook output cell
  - State prerequisites: notebook served from HTTPS/remote origin; `jupyter-server-proxy` not installed or not working
  - Action to reach state: run `show_anndata(...)` and observe the iframe error panel
Crop:
  - Include: the iframe error content and the surrounding notebook cell output
  - Exclude: private notebook names, hostnames, user names
Redact:
  - Remove: any internal URLs/domains if sensitive
Alt text:
  - Embedded Cellucid viewer showing a “notebook proxy required” message.
Caption:
  - When the notebook page is served from HTTPS, the browser may block an HTTP loopback iframe; installing `jupyter-server-proxy` fixes this in most managed notebook environments.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Embedded Cellucid viewer showing a “notebook proxy required” message.
:width: 100%

When the notebook page is served from HTTPS, the browser may block an HTTP loopback iframe; installing `jupyter-server-proxy` fixes this in most managed notebook environments.
```

### 2) UI asset availability (first-run download + caching)

Cellucid’s Python server runs in **hosted-asset proxy** mode:
- it downloads the viewer UI from `https://www.cellucid.com`,
- caches it on disk,
- and serves it from the same origin as your dataset server.

If you’re offline on your first run (or behind a firewall), you may see a “viewer UI unavailable” error page.

Fixes:
- prefetch UI assets once while online (recommended),
- set `CELLUCID_WEB_PROXY_CACHE_DIR` to a persistent, writable location.

## Edge cases (things that look like “bugs”)

- **Ad blockers / strict corporate CSP**: can break the hosted web UI fetch or iframe behavior.
- **Port collisions**: if 8765 is busy, servers may choose a nearby port; always copy the printed URL.
- **Multiple viewers in one notebook**: hooks are routed by `viewerId`. If you re-run cells many times, stop old viewers (`viewer.stop()`) or restart the kernel.
- **Remote access**: binding `--host 0.0.0.0` exposes your server to the network; do this only if you understand the security implications.

## Troubleshooting (compatibility)

### Symptom: viewer iframe is blank / white

**Likely causes**
- Browser blocked the iframe (mixed content).
- The viewer UI assets were not available (offline / firewall).
- The browser can’t reach the data server (remote kernel without port forwarding).

**How to confirm**
- In Python:
  ```python
  report = viewer.debug_connection()
  report
  ```
- In the browser: open DevTools → Console.

**Fix**
- Install `jupyter-server-proxy` for HTTPS/remote notebooks.
- Use port forwarding + `CELLUCID_CLIENT_SERVER_URL` for remote kernels.
- Prefetch UI assets once while online.

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
