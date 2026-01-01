# Architecture: message routing (HTTP vs postMessage)

This page is the “deep path” explanation of how notebook embedding + hooks are implemented.

If you just want something that works, start with:
- {doc}`02_quickstart_minimal_roundtrip`

## Design goals

Cellucid’s notebook integration is designed to:

1. Work across notebook frontends (classic, JupyterLab, VSCode notebooks, Colab).
2. Avoid fragile notebook-specific transport mechanisms.
3. Make failures diagnosable (clear error pages + `viewer.debug_connection()`).
4. Keep data local by default (localhost server, not a cloud relay).

## System diagram

```text
┌─────────────────────────────── Notebook Kernel ───────────────────────────────┐
│                                                                               │
│  Python: cellucid-python                                                      │
│   - BaseViewer / CellucidViewer / AnnDataViewer                               │
│   - CellucidServer / AnnDataServer                                            │
│   - /_cellucid/events (hooks)                                                 │
│   - /_cellucid/session_bundle (uploads)                                       │
│                                                                               │
└───────────────────────────────────┬───────────────────────────────────────────┘
                                    │ HTTP (browser → server)
                                    │
┌───────────────────────────────────▼───────────────────────────────────────────┐
│                           Browser (notebook frontend)                         │
│                                                                               │
│  Output cell contains an iframe:                                               │
│    <iframe src=".../?jupyter=true&viewerId=...&viewerToken=...">               │
│                                                                               │
│  The iframe runs the Cellucid web app UI (served by the same server origin).  │
│                                                                               │
└───────────────┬───────────────────────────────────────────────────────┬───────┘
                │ postMessage (Python → iframe)                          │
                │ (commands)                                             │
                │                                                       │ HTTP POST
                ▼                                                       │ (events)
┌───────────────────────────────────────────────────────────────────────▼───────┐
│                           Iframe (Cellucid UI)                                 │
│                                                                               │
│  JupyterBridgeDataSource (web app)                                             │
│   - receives postMessage commands (requires viewerToken)                       │
│   - sends events via POST /_cellucid/events                                     │
│                                                                               │
└───────────────────────────────────────────────────────────────────────────────┘
```

## Why there are two channels

### Channel A: Python → iframe (commands) via `postMessage`

Why `postMessage`?
- it’s the universal way to talk to an iframe,
- it does not require any server-to-browser persistent socket,
- it works in restricted notebook frontends (Colab, VSCode webviews).

Security:
- notebook frontends can proxy/transform origins,
- so Cellucid authenticates commands using a per-viewer `viewerToken` rather than strict origin allowlisting.

Where to look:
- Python send path: `cellucid-python/src/cellucid/jupyter.py` (`BaseViewer.send_message`)
- JS receive path: `cellucid/assets/js/data/jupyter-source.js` (`_handleMessage`)

### Channel B: iframe → Python (events) via HTTP POST

Why HTTP POST?
- it works everywhere (no notebook-specific comm channels),
- the Python server already exists (serving the dataset), so it can receive events too,
- it’s robust under proxies and iframes (compared to WebSockets in some contexts).

Where to look:
- JS send path: `cellucid/assets/js/data/jupyter-source.js` (`_postEventToPython`)
- Python receive path: `cellucid-python/src/cellucid/_server_base.py` (`handle_event_post`)
- Python routing path: `cellucid-python/src/cellucid/_server_base.py` (`route_event`)
- Python hook dispatch: `cellucid-python/src/cellucid/jupyter.py` (`_handle_frontend_message`, `HookRegistry`)

## How the iframe URL is chosen (mixed-content avoidance)

The core constraint:
- many notebooks are served on **HTTPS** (JupyterHub, hosted Jupyter),
- browsers often block if an HTTPS page tries to load or fetch from an **HTTP loopback** server (`http://127.0.0.1:<port>`) as *active mixed content*.

Cellucid’s solution:
1. Serve the viewer UI from the same origin as the dataset server (hosted-asset proxy).
2. When the notebook is HTTPS/remote, prefer a notebook proxy URL (Jupyter Server Proxy).

In the notebook output HTML (`BaseViewer._generate_viewer_html`):
- it builds a `directSrc` (often loopback),
- it also tries a `proxySrc` of the form:

```text
https://<notebook-origin>/<base>/proxy/<port>/?jupyter=true&viewerId=...&viewerToken=...
```

and probes `/_cellucid/health` before choosing it.

For Colab:
- Python computes an HTTPS proxy URL via `google.colab.kernel.proxyPort(port)` (best effort).

Override (advanced):
- `CELLUCID_CLIENT_SERVER_URL` forces the browser-facing base URL for the server.

## Hosted-asset proxy (web UI caching)

Notebook embeds load the Cellucid UI from the local server origin, but those HTML/JS/CSS assets are fetched from `https://www.cellucid.com` when not packaged locally.

Behavior:
- The server proxies `index.html` and `/assets/*` from `CELLUCID_WEB_URL`.
- Files are cached to disk (default: temp dir).
- Cache invalidation uses the `<meta name="cellucid-web-build-id" content="...">` stamp in `index.html`.

Controls:
- Cache directory: `CELLUCID_WEB_PROXY_CACHE_DIR=/path/to/cache`
- Prefetch in notebooks (progress bar): `viewer.ensure_web_ui_cached()`
- Clear cache: `viewer.clear_web_cache()` or `cellucid.clear_web_cache()`

## Viewer identity: `viewerId` and routing

Multiple viewers in one kernel are supported.

How routing works:
- Each Python viewer instance generates a unique `viewerId`.
- The server keeps a process-global registry mapping `viewerId → callback`.
- Incoming events include `viewerId` and are delivered to the correct viewer’s `_handle_frontend_message`.

This is why “stale iframes” happen:
- if you restart the kernel, the registry resets,
- an old iframe may still be open and sending events to a server that no longer knows its `viewerId`.

Fix: re-run the cell to create a new viewer.

## Event dispatch inside Python

The flow for a selection event is:

```text
iframe -> POST /_cellucid/events {type:"selection", viewerId:"...", cells:[...], source:"lasso"}
server -> route_event(viewerId, event)
viewer._handle_frontend_message(event)
  - maps type → event name
  - strips internal fields (type, viewerId)
  - updates viewer.state + wakes waiters
  - triggers HookRegistry callbacks (with exception guards)
```

## Threading / reentrancy notes

Events arrive on the HTTP server’s request-handling thread.

Implications:
- keep callbacks fast,
- avoid long-running work directly in the hook handler,
- be careful when touching backed AnnData from multiple threads.

See {doc}`08_writing_robust_callbacks`.

## Debugging primitives built into the architecture

Use `viewer.debug_connection()` to collect:
- server health/info probes
- ping/pong roundtrip
- a frontend “debug snapshot” (`location.href`, origin, user agent)
- recent forwarded frontend warnings/errors (`console` events)

This is the fastest way to determine whether your failure is:
- a notebook proxy/mixed-content problem,
- a stale iframe / wrong viewerId problem,
- or a real application bug.

## Next steps

- Commands: {doc}`06_python_to_frontend_commands`
- Events: {doc}`07_frontend_to_python_events`
- Security: {doc}`12_security_model` and {doc}`13_security_cors_origins_and_mixed_content`
- Troubleshooting: {doc}`14_troubleshooting_hooks`

