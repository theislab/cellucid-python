# Supported environments matrix

This page answers: “Will notebook embedding + hooks work in *my* environment?”

Short answer: **yes, usually**, as long as the **browser can reach the Cellucid server URL** used by the embedded iframe.

The hooks system is designed to avoid brittle notebook-specific machinery:
- **Python → viewer** commands use `postMessage` into the iframe.
- **Viewer → Python** events use **HTTP POST** to `/_cellucid/events` on the local server.

## At a glance

**Audience**
- Wet lab / beginner: check your environment row, then use the quickstart.
- Computational users: pay attention to HTTPS/remote and to “browser can reach the server”.
- Developers/IT: read the proxy + security notes and the network requirements.

**Start here**
- Quickstart sanity check: {doc}`02_quickstart_minimal_roundtrip`
- Mixed content / HTTPS details: {doc}`13_security_cors_origins_and_mixed_content`

## Matrix

Legend:
- ✅ should work with defaults
- ⚠️ works, but you may need one extra step
- ❌ not supported (or not a notebook environment)

| Environment | Embed viewer | Python → viewer commands | Viewer → Python events | Session bundle “no-download” | Notes |
|---|---:|---:|---:|---:|---|
| Classic Jupyter (local) | ✅ | ✅ | ✅ | ✅ | Direct loopback (`http://127.0.0.1:<port>`) usually works. |
| JupyterLab (local) | ✅ | ✅ | ✅ | ✅ | Same as classic; base URL/proxy path handling is supported. |
| VSCode notebooks (local kernel) | ✅ | ✅ | ✅ | ✅ | Works as long as the VSCode webview can reach the server URL. |
| Google Colab | ✅ | ✅ | ✅ | ✅ | Uses Colab’s HTTPS port proxy (kernel is remote). Cache is ephemeral. |
| JupyterHub / hosted Jupyter (HTTPS) | ⚠️ | ✅ | ✅ | ✅ | Usually needs **Jupyter Server Proxy** to avoid HTTPS→HTTP mixed-content blocking. |
| Remote kernel over SSH (browser on laptop) | ⚠️ | ✅ | ✅ | ✅ | Use **SSH port forwarding** or set `CELLUCID_CLIENT_SERVER_URL`. |
| “Notebook” served from file:// (rare) | ⚠️ | ✅ | ✅ | ✅ | Proxy selection logic avoids file:// issues; you may need manual URLs. |
| Not in Jupyter (plain Python script) | ❌ | ⚠️ | ⚠️ | ⚠️ | You can still run a server (`cellucid serve ...`), but there’s no notebook iframe or hooks UI. |

```{important}
The biggest real-world failure mode is **connectivity**, not hooks logic:

> The viewer iframe is loaded in a browser tab, but the Python server is not reachable from that browser origin.

When this happens, the UI may show errors like “Failed to fetch”, selection events won’t arrive, and highlighting commands won’t do anything.
```

## Network requirement: viewer UI asset download (hosted-asset proxy)

In notebooks, Cellucid serves the viewer UI via a **hosted-asset proxy**:

- On first run (or after the website build changes), the Python side downloads the viewer UI assets from:
  - `https://www.cellucid.com/index.html`
  - `https://www.cellucid.com/assets/*`
- Assets are cached on disk (configure with `CELLUCID_WEB_PROXY_CACHE_DIR`).
- If you are offline but the cache exists, embedding still works.

See:
- {doc}`13_security_cors_origins_and_mixed_content` (why this exists + mixed content)
- {doc}`14_troubleshooting_hooks` (“Viewer UI could not be loaded”)

## Environment-specific notes (with actionable fixes)

### JupyterHub / HTTPS notebooks: “notebook proxy required”

If your notebook is served over HTTPS and your Cellucid server is plain HTTP on loopback, browsers may block it as mixed content.

Cellucid tries to auto-fix this by embedding through:

```text
https://<notebook-origin>/<base>/proxy/<port>/?jupyter=true&viewerId=...&viewerToken=...
```

That requires `jupyter-server-proxy` to be installed/enabled on the notebook server.

If you see a “notebook proxy required” message inside the iframe:
- install `jupyter-server-proxy` on the Jupyter server, or
- set `CELLUCID_CLIENT_SERVER_URL` to a browser-reachable HTTPS URL for the Cellucid server (advanced), or
- use SSH port forwarding (if your notebook runs on a remote host).

### Google Colab

Colab runs the kernel on a remote VM. Cellucid uses:

```python
from google.colab.output import eval_js
eval_js("google.colab.kernel.proxyPort(<port>)")
```

to obtain an HTTPS URL the browser can reach.

Practical implications:
- The viewer URL will *not* be `127.0.0.1` in Colab.
- The UI cache is not persistent between Colab sessions, so first-run downloads are common.

### Remote / HPC kernels (browser on laptop)

If your notebook kernel runs on a remote machine but you view the notebook in a browser on your laptop, you must bridge the server port.

Most robust: SSH local port forwarding (example uses port 8765):

1. In the notebook (remote):
   ```python
   viewer = show_anndata("data.h5ad", port=8765)
   print(viewer.viewer_url)
   ```
2. On your laptop:
   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now `http://127.0.0.1:8765/...` in your browser forwards to the remote kernel.

## Edge cases

- **Ad blockers / corporate proxies** can block POSTs to `/_cellucid/events` (hooks) or block downloads from `https://www.cellucid.com`.
- **Kernel restarts** can orphan an iframe tab; re-run the cell to create a fresh viewer/server pair.
- **Multiple viewers** in one notebook are supported (routing uses `viewerId`), but can be confusing during debugging—use `viewer.debug_connection()`.

## Troubleshooting

- Run:
  ```python
  viewer.debug_connection()
  ```
- Then follow: {doc}`14_troubleshooting_hooks`

## Next steps

- Minimal working loop: {doc}`02_quickstart_minimal_roundtrip`
- Architecture details: {doc}`05_architecture_message_routing_http_vs_postmessage`

