# Glossary and mental model

Hooks are much easier to debug if everyone uses the same words for the same things. This page defines the key terms and gives you a “one diagram” mental model.

## Glossary (terms used in this chapter)

### Cellucid (web app)
The browser-based viewer UI. In notebooks, it runs inside a widget in the cell's output.

### cellucid-python (Python package)
The helper package you use from Python/CLI. It can:
- serve datasets (an AnnData store the viewer reads itself),
- show the Cellucid UI in notebooks (a widget),
- send commands to the UI and receive events back (hooks),
- capture session bundles and apply them to AnnData.

### Viewer / `viewer`
The Python object returned by `show_anndata(...)`: a `cellucid.Widget`.

It owns:
- the dataset it answers the viewer's reads from
- a hooks registry (`viewer.on_selection`, etc.)
- a small latest-state snapshot (`viewer.state`) and the viewer's results
  (`viewer.results`, `viewer.pull(adata)`)

A browser-tab viewer (`AnnDataViewer`) owns the same, plus a running
background server (`viewer.server_url`) and the tab's address
(`viewer.viewer_url`).

### Widget
The notebook output that shows the viewer: an anywidget whose front end
frames the Cellucid app and relays between it and the kernel over the
notebook's own connection (a Jupyter comm). No port, no server.

### Data server
A small local HTTP server, started by a browser-tab viewer (`AnnDataViewer`,
`serve_anndata`, `cellucid serve`). It serves:
- the AnnData store's bytes, read by the viewer with range requests
- the exact verified viewer generation
- special endpoints like:
  - `/_cellucid/health`
  - `/_cellucid/info`
  - `/_cellucid/protocol` (where the viewer channel is, and what it carries)
  - `/_cellucid/ws` (the viewer channel, a WebSocket)

### Viewer channel
The one channel the viewer opens to Python: over the widget's comm in a
notebook, a WebSocket to the data server (`/_cellucid/ws`) in a browser tab.
Both commands and events travel on it. See {doc}`05_architecture_viewer_channel`.

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
A unique ID naming one viewer instance (important when you have multiple viewers in one notebook): the viewer says it when it opens its channel, and Python binds the channel to that viewer.

### `viewerToken`
A per-viewer secret token: a widget hands it to its page in a message, a
browser tab carries it in its address. The viewer sends it once, when it
opens its channel, and Python compares it in constant time before
the channel carries anything. A page without the token cannot open a channel,
so no other page can send events as the viewer or receive its commands.

### Origin, CORS, mixed content
Browser security concepts that explain many failures:
- **Origin**: scheme + host + port (e.g. `https://example.com:443`)
- **CORS**: rules about which origins can make requests to which servers
- **Mixed content**: browsers blocking `https://...` pages that try to load `http://...` resources

See {doc}`13_security_cors_origins_and_mixed_content`.

### Verified local web generation
A browser-tab server establishes the viewer UI assets declared by `https://www.cellucid.com` (verified, cached on disk) at every start and serves them from the same origin as the dataset. A notebook widget keeps none: the browser loads the app from `web_source_url`.

Controls (browser-tab viewers):
- Generation location: `web_cache_dir=...`
- Clear cache: `cellucid.clear_web_cache()` or `viewer.clear_web_cache()`

## Mental model (one diagram)

```text
┌──────────────────────────────────────────────────────────────────────┐
│                            Your notebook                              │
│                                                                      │
│  Python kernel                                                       │
│   ├─ viewer = show_anndata(...)   (a cellucid.Widget)                │
│   ├─ answers the viewer's reads of the dataset                       │
│   └─ registers hooks: @viewer.on_selection / @viewer.on_hover / ...  │
│                                                                      │
│  Browser (notebook frontend)                                         │
│   └─ the widget's output runs the Cellucid web app                   │
│      (loaded from https://www.cellucid.com)                          │
│                                                                      │
│  The notebook's own connection carries one viewer channel:           │
│   Python → viewer: commands      viewer → Python: events             │
└──────────────────────────────────────────────────────────────────────┘
```

If “hooks don’t work”, it’s almost always timing or a stale output:
- a cell is still running (the widget's events arrive between cells)
- the browser could not load the app (the output says so)
- a stale output after a kernel restart (it keeps its last view and says the
  notebook stopped it)

Commands sent before the viewer is shown are not lost: they wait in
Python and arrive, in order, when the viewer connects.

## Common misconceptions (and the correct model)

### “Hooks are a Jupyter feature”
Not quite. Hooks are messages on one viewer channel. In a notebook the
channel rides on the widget's comm (so any front end that shows widgets
works); a browser tab carries the same channel on a WebSocket.

### “If I can see the viewer, hooks must work”
Not necessarily:
- you can see the UI while a cell runs, but its events reach Python only
  when the cell ends;
- in a browser tab behind a proxy, the proxy may refuse the WebSocket; the
  viewer then says it cannot open its channel and names the proxy.

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
