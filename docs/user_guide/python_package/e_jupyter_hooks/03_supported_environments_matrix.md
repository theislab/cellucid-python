# Supported environments matrix

This page answers: “Will notebook embedding + hooks work in *my* environment?”

Short answer: **yes, usually**, as long as the **browser can reach the Cellucid server URL** used by the embedded iframe.

The hooks system uses no notebook-specific machinery. The notebook output is a
plain iframe (no script), and the viewer in it opens **one WebSocket** to the
same server it was loaded from (`/_cellucid/ws`); commands and events both
travel on it ({doc}`05_architecture_viewer_channel`). So every front end that
can show the iframe behaves the same, and there are exactly two requirements:

1. the browser can load the iframe's URL (the Cellucid server, directly or
   through a proxy), and
2. anything between the browser and the server passes WebSockets
   (`jupyter-server-proxy` does; check reverse proxies in front of it).

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

Commands, events and the session bundle “no-download” workflow all travel on
the one channel, so they work exactly when the viewer channel opens; the table
says what that takes in each setting.

| Environment | What the browser must reach | Extra step |
|---|---|---|
| Classic Jupyter, JupyterLab (local) | `http://127.0.0.1:<port>` directly | none (✅) |
| VS Code notebooks (local kernel), untested | the loopback server from VS Code's webview | none, when the webview can load loopback URLs |
| Google Colab, untested | Colab's HTTPS proxy for the port, with WebSockets | pass the proxy base as `client_server_url=` and its host as `allowed_hosts=` (⚠️) |
| JupyterHub / hosted Jupyter (HTTPS) | an HTTPS proxy route (e.g. `jupyter-server-proxy`) that passes WebSockets | pass its base as `client_server_url=` and its host as `allowed_hosts=` (⚠️) |
| Remote kernel over SSH (browser on laptop) | the forwarded port | forward the port and pass the forwarded base as `client_server_url=` (⚠️) |
| “Notebook” served from `file://` | an HTTP(S) server URL | pass an exact HTTP(S) `client_server_url`; `file://` is not an accepted server URL (⚠️) |
| Not in Jupyter (plain Python script) | — | use `cellucid serve ...`; notebook display and hooks need a notebook front end (❌) |

What has been checked, and what has not: the channel is tested end to end
against a real browser for a page loaded straight from the server (classic
Jupyter and JupyterLab on the same machine, a browser tab, an SSH-forwarded
port). **VS Code's notebook webview and Colab's proxy have not been tested by
the project.** In VS Code with a local kernel the iframe's page is itself
served from `http://127.0.0.1:<port>`, so its WebSocket is a same-origin
loopback connection that no proxy touches; whether a given VS Code build lets
the webview load a loopback iframe at all is what decides it. Colab's
`proxyPort` URL must also carry the `/_cellucid/ws` upgrade, which is
Colab's behaviour, not Cellucid's. In either case `viewer.debug_connection()`
says exactly which part fails (`channel_connected`, `frontend_roundtrip`).

```{important}
The biggest real-world failure mode is **connectivity**, not hooks logic:

> The viewer iframe is loaded in a browser tab, but the Python server is not reachable from that browser origin.

When the page loads but the WebSocket cannot open (a proxy that does not pass
WebSockets), the viewer says so: “The notebook or R session is running, but this viewer
cannot open its channel (…). A proxy between this page and the host must pass
WebSockets.” Commands sent meanwhile wait in Python and arrive once a viewer
connects.
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
- {doc}`14_troubleshooting_hooks` (“the viewer UI cannot be established or served”)

## Environment-specific notes (with actionable fixes)

### JupyterHub / HTTPS notebooks: browser cannot reach loopback

If your notebook is served over HTTPS and your Cellucid server is plain HTTP on loopback, browsers may block it as mixed content.

Expose the Cellucid port through an HTTPS proxy, then construct its exact base
URL, for example:

```text
https://<notebook-origin>/<base>/proxy/<port>/?viewer=true&viewerId=...&viewerToken=...&anndata=true
```

With Jupyter Server Proxy, the notebook administrator must install and enable
that extension's server route, and the caller passes the resulting browser
base as `client_server_url=...` and the hub's host name as
`allowed_hosts=[...]` (the server answers only to `Host` names and page
origins it knows). The viewer channel then runs at
`wss://<notebook-origin>/<base>/proxy/<port>/_cellucid/ws`; Jupyter Server
Proxy passes WebSockets, and any reverse proxy in front of the hub must too.
The Cellucid package already depends on `jupyter-server-proxy`; package
installation alone does not expose a route or select its URL.

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

Then pass the returned base, without its trailing `/`, to
`show_anndata(..., client_server_url=base, allowed_hosts=[<proxy host>])`
(the host name of that base). The proxy URL names one port, and
`show_anndata` serves on 8765 or, when 8765 is taken, the next free port; to
be sure the two agree, construct `AnnDataViewer(data, port=<port>,
client_server_url=base, allowed_hosts=[...])`, which uses that port or fails.
This setup has not been tested by the project. Colab's proxy must pass the
viewer channel's WebSocket for commands and hooks to work; if the viewer
reports that it cannot open its channel, use `cellucid serve` in a browser tab
instead.
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
`http://127.0.0.1:8765/?viewer=true&viewerId=...&viewerToken=...` URL in your
browser forwards to the remote kernel.

## Edge cases

- **Corporate proxies** can drop WebSocket connections (the viewer then names the proxy as the cause) or block downloads from `https://www.cellucid.com`. An idle channel is kept open by a ping from the server every 20 s.
- **Kernel restarts** stop the server; the old iframe freezes its last view and says Python stopped. Re-run the cell to create a fresh viewer/server pair.
- **Multiple viewers** in one notebook are supported (routing uses `viewerId`), but can be confusing during debugging—use `viewer.debug_connection()`.

## Troubleshooting

- Run:
  ```python
  viewer.debug_connection()
  ```
- Then follow: {doc}`14_troubleshooting_hooks`

## Next steps

- Minimal working loop: {doc}`02_quickstart_minimal_roundtrip`
- Architecture details: {doc}`05_architecture_viewer_channel`
