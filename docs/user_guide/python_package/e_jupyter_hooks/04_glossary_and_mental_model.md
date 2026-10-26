# Glossary and mental model

Hooks are much easier to debug if everyone uses the same words for the same things. This page defines the key terms and gives you a “one diagram” mental model.

## Glossary (terms used in this chapter)

### Cellucid (web app)
The browser-based viewer UI. In notebooks, this UI is embedded in an iframe.

### cellucid-python (Python package)
The helper package you use from Python/CLI. It can:
- serve datasets (an AnnData store the viewer reads itself),
- embed the Cellucid UI in notebooks,
- send commands to the UI and receive events back (hooks),
- capture session bundles and apply them to AnnData.

### Viewer / `viewer`
The Python object returned by `show_anndata(...)`.

It owns:
- a running background server (`viewer.server_url`)
- an iframe URL (`viewer.viewer_url`)
- a hooks registry (`viewer.on_selection`, etc.)
- a small latest-state snapshot (`viewer.state`)

### Data server
A small local HTTP server started by `viewer`. It serves:
- the AnnData store's bytes, read by the viewer with range requests
- the exact verified viewer generation
- special endpoints like:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/protocol` (where the viewer channel is, and what it carries)
  - `/_cellucid/ws` (the viewer channel, a WebSocket)

### Viewer channel
The one WebSocket the viewer opens to the data server (`/_cellucid/ws`). Both
commands and events travel on it. See {doc}`05_architecture_viewer_channel`.

### Iframe
An embedded browser frame inside the notebook output. It runs the Cellucid web app UI.

### Command (Python → viewer)
A message Python sends to the viewer on the viewer channel. Examples:
- highlight cells
- set color-by
- reset camera

### Event (viewer → Python)
A message the viewer sends to Python on the viewer channel. Examples:
- selection
- hover
- click
- ready

### `viewerId`
A unique ID naming one viewer instance (important when you have multiple viewers in one notebook): the viewer says it when it opens its channel, and the server binds the channel to that viewer.

### `viewerToken`
A per-viewer secret token in the iframe's address. The viewer sends it once,
when it opens its channel, and the server compares it in constant time before
the channel carries anything. A page without the token cannot open a channel,
so no other page can send events as the viewer or receive its commands.

### Origin, CORS, mixed content
Browser security concepts that explain many failures:
- **Origin**: scheme + host + port (e.g. `https://example.com:443`)
- **CORS**: rules about which origins can make requests to which servers
- **Mixed content**: browsers blocking `https://...` pages that try to load `http://...` resources

See {doc}`13_security_cors_origins_and_mixed_content`.

### Verified local web generation
In notebook mode, the Python server establishes the viewer UI assets declared by `https://www.cellucid.com` (verified, cached on disk) at every viewer start and serves them from the same origin as the dataset server.

This avoids:
- mixed-content failures (HTTPS notebook embedding HTTP resources),
- and many cross-origin issues.

Controls:
- Generation location: `web_cache_dir=...`
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
│  One viewer channel (WebSocket /_cellucid/ws):                       │
│   Python → viewer: commands      viewer → Python: events             │
└──────────────────────────────────────────────────────────────────────┘
```

If “hooks don’t work”, it’s almost always because one side can’t reach the other:
- iframe can’t reach the server URL (proxy / mixed content / firewall)
- a proxy in between does not pass WebSockets (the viewer says so)
- a stale iframe after a kernel restart (it freezes and says Python stopped)

Commands sent before the iframe is displayed are not lost: they wait in
Python and arrive, in order, when the viewer connects.

## Common misconceptions (and the correct model)

### “Hooks are a Jupyter feature”
They aren’t. Hooks are messages on one WebSocket between the viewer and the local server, which works in any front end that can show the iframe.

### “If I can see the viewer, hooks must work”
Not necessarily:
- you can see the UI (loaded over HTTP), but a proxy may refuse the WebSocket;
  the viewer then says it cannot open its channel and names the proxy.

### “Cell indices are stable cell IDs”
In the hooks API, indices are **row positions** (0-based). If you reorder/subset the dataset in Python, indices change.

### “I can treat hover events like clicks”
Hover is a debounced, high-frequency state stream rather than a discrete click
record. Keep callbacks light and do not assume that every intermediate pointer
position produces an event.

## Next steps

- Architecture deep dive: {doc}`05_architecture_viewer_channel`
- Commands and events:
  - {doc}`06_python_to_frontend_commands`
  - {doc}`07_frontend_to_python_events`
- Troubleshooting: {doc}`14_troubleshooting_hooks`
