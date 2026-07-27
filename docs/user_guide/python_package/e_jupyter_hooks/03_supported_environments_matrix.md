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
| JupyterLab (local) | ✅ | ✅ | ✅ | ✅ | Same as classic when the browser can reach loopback. |
| VSCode notebooks (local kernel) | ✅ | ✅ | ✅ | ✅ | Works as long as the VSCode webview can reach the server URL. |
| Google Colab | ⚠️ | ✅ | ✅ | ✅ | Obtain Colab’s browser-reachable HTTPS proxy base and pass it as `client_server_url=`. |
| JupyterHub / hosted Jupyter (HTTPS) | ⚠️ | ✅ | ✅ | ✅ | Expose the Cellucid port through an HTTPS proxy and pass its exact base as `client_server_url=`. |
| Remote kernel over SSH (browser on laptop) | ⚠️ | ✅ | ✅ | ✅ | Forward the port, then pass the forwarded browser base as `client_server_url=`. |
| “Notebook” served from `file://` | ⚠️ | ✅ | ✅ | ✅ | Pass an exact HTTP(S) `client_server_url`; `file://` is not an accepted server URL. |
| Not in Jupyter (plain Python script) | ❌ | ❌ | ❌ | ❌ | Use `cellucid serve ...`; notebook display and hook channels require a notebook frontend. |

```{important}
The biggest real-world failure mode is **connectivity**, not hooks logic:

> The viewer iframe is loaded in a browser tab, but the Python server is not reachable from that browser origin.

When this happens, the UI may show errors like “Failed to fetch”, selection events won’t arrive, and highlighting commands won’t do anything.
```

## Network requirement: exact viewer generation

In notebooks, Cellucid serves one verified local viewer generation:

- At each viewer/server startup, Python establishes the complete generation
  declared by `https://www.cellucid.com/cellucid-web-assets.json`.
- The generation includes `index.html`, root browser metadata, and all declared
  assets; every byte is verified before publication.
- The selected generation lives on disk (configure with `web_cache_dir=`).
- Startup requires access to the configured source and never substitutes a
  previous generation after source failure.

See:
- {doc}`13_security_cors_origins_and_mixed_content` (why this exists + mixed content)
- {doc}`14_troubleshooting_hooks` (“Viewer UI could not be loaded”)

## Environment-specific notes (with actionable fixes)

### JupyterHub / HTTPS notebooks: browser cannot reach loopback

If your notebook is served over HTTPS and your Cellucid server is plain HTTP on loopback, browsers may block it as mixed content.

Expose the Cellucid port through an HTTPS proxy, then construct its exact base
URL, for example:

```text
https://<notebook-origin>/<base>/proxy/<port>/?jupyter=true&viewerId=...&viewerToken=...
```

With Jupyter Server Proxy, the notebook administrator must install and enable
that extension's server route and the caller must pass the resulting browser
base as `client_server_url=...`. The Cellucid package already depends on
`jupyter-server-proxy`; package installation alone does not expose a route or
select its URL.

If the direct iframe is blocked:

- expose the port through an HTTPS reverse proxy and pass its exact base as
  `client_server_url=...`, or
- use SSH port forwarding and pass the local forwarded base.

### Google Colab

Colab runs the kernel on a remote VM. Obtain its proxy URL in the notebook:

```python
from google.colab.output import eval_js
eval_js("google.colab.kernel.proxyPort(<port>)")
```

Then pass the returned base to `show_anndata(..., client_server_url=base)`.
Cellucid does not invoke Colab’s proxy API implicitly.

Practical implications:
- The viewer URL must be the explicit HTTPS proxy URL, not remote loopback.
- Each viewer startup still establishes the configured exact web generation.

### Remote / HPC kernels (browser on laptop)

If your notebook kernel runs on a remote machine but you view the notebook in a browser on your laptop, you must bridge the server port.

Most robust: SSH local port forwarding (example uses port 8765):

1. In the notebook (remote):
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
2. On your laptop:
   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now the exact printed
`http://127.0.0.1:8765/?jupyter=true&viewerId=...&viewerToken=...` URL in your
browser forwards to the remote kernel.

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
