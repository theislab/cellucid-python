# Glossary and mental model

Hooks are much easier to debug if everyone uses the same words for the same things. This page defines the key terms and gives you a “one diagram” mental model.

## Glossary (terms used in this chapter)

### Cellucid (web app)
The browser-based viewer UI. In notebooks, this UI is embedded in an iframe.

### cellucid-python (Python package)
The helper package you use from Python/CLI. It can:
- serve datasets (exported or AnnData-backed),
- embed the Cellucid UI in notebooks,
- send commands to the UI and receive events back (hooks),
- capture session bundles and apply them to AnnData.

### Viewer / `viewer`
The Python object returned by `show(...)` or `show_anndata(...)`.

It owns:
- a running background server (`viewer.server_url`)
- an iframe URL (`viewer.viewer_url`)
- a hooks registry (`viewer.on_selection`, etc.)
- a small latest-state snapshot (`viewer.state`)

### Data server
A small local HTTP server started by `viewer`. It serves:
- dataset files (exported mode) or dynamic endpoints (AnnData mode)
- the viewer UI (via hosted-asset proxy)
- special endpoints like:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/events` (viewer → Python hooks)
  - `/_cellucid/session_bundle` (session bundle uploads)

### Iframe
An embedded browser frame inside the notebook output. It runs the Cellucid web app UI.

### Command (Python → viewer)
A message sent from Python to the iframe using `postMessage`. Examples:
- highlight cells
- set color-by
- reset camera

### Event (viewer → Python)
A message sent from the iframe to Python using HTTP POST to `/_cellucid/events`. Examples:
- selection
- hover
- click
- ready

### `viewerId`
A unique ID used to route events and commands to the correct viewer instance (important when you have multiple viewers in one notebook).

### `viewerToken`
A per-viewer secret token used to authenticate `postMessage` commands sent into the iframe.

Why it exists:
- notebook environments can proxy/transform iframe origins (Colab, VSCode webviews),
- so the integration does not rely on strict origin allowlisting for `postMessage`,
- instead it requires the correct token.

### Origin, CORS, mixed content
Browser security concepts that explain many failures:
- **Origin**: scheme + host + port (e.g. `https://example.com:443`)
- **CORS**: rules about which origins can make requests to which servers
- **Mixed content**: browsers blocking `https://...` pages that try to load `http://...` resources

See {doc}`13_security_cors_origins_and_mixed_content`.

### Hosted-asset proxy (web UI cache)
In notebook mode, the Python server can download the viewer UI assets from `https://www.cellucid.com` and serve them from the same origin as the dataset server (cached on disk).

This avoids:
- mixed-content failures (HTTPS notebook embedding HTTP resources),
- and many cross-origin issues.

Controls:
- Cache location: `CELLUCID_WEB_PROXY_CACHE_DIR`
- Clear cache: `cellucid.clear_web_cache()` or `viewer.clear_web_cache()`

## Mental model (one diagram)

```text
┌──────────────────────────────────────────────────────────────────────┐
│                            Your notebook                              │
│                                                                      │
│  Python kernel                                                       │
│   ├─ viewer = show_anndata(...)                                      │
│   ├─ starts local HTTP server: http://127.0.0.1:<port>               │
│   ├─ embeds iframe: viewer.viewer_url                                │
│   └─ registers hooks: @viewer.on_selection / @viewer.on_hover / ...  │
│                                                                      │
│  Browser (notebook frontend)                                         │
│   └─ iframe runs the Cellucid web app UI                             │
│                                                                      │
│  Two communication channels:                                         │
│   (A) Python → iframe: postMessage (commands)                        │
│   (B) iframe → Python: HTTP POST /_cellucid/events (events)          │
└──────────────────────────────────────────────────────────────────────┘
```

If “hooks don’t work”, it’s almost always because one side can’t reach the other:
- iframe can’t reach the server URL (proxy / mixed content / firewall)
- Python is sending commands before the iframe is displayed
- `viewerId` mismatch (stale iframe after kernel restart)

## Common misconceptions (and the correct model)

### “Hooks are a Jupyter feature”
They aren’t. Hooks are just **HTTP POST requests** from the iframe to the server, and `postMessage` commands from Python to the iframe.

### “If I can see the viewer, hooks must work”
Not necessarily:
- you can see the UI, but the iframe may not be able to reach `/_cellucid/events` (blocked by proxy/ad-blocker),
- or the dataset loads but event POSTs are blocked.

### “Cell indices are stable cell IDs”
In the hooks API, indices are **row positions** (0-based). If you reorder/subset the dataset in Python, indices change.

### “I can treat hover events like clicks”
Hover is debounced and high frequency; treat it as “best effort” and keep callbacks light.

## Next steps

- Architecture deep dive: {doc}`05_architecture_message_routing_http_vs_postmessage`
- Commands and events:
  - {doc}`06_python_to_frontend_commands`
  - {doc}`07_frontend_to_python_events`
- Troubleshooting: {doc}`14_troubleshooting_hooks`

