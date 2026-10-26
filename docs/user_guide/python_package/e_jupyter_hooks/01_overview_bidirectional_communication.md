# Overview: bidirectional communication (Python ↔ viewer)

Cellucid’s notebook integration is intentionally **simple and robust**:

- **Python** starts a small **local HTTP server** that serves your dataset.
- The notebook displays the Cellucid **web app UI** in an **iframe**.
- The viewer opens **one channel** to that server, a WebSocket at `/_cellucid/ws`, and both directions travel on it:
  1. **Python → viewer**: commands (e.g., highlight, set color-by).
  2. **Viewer → Python**: events (e.g., selection, hover, click).

This design needs no Jupyter extension and no code in the notebook output (the output is a plain iframe), so it works the same in every notebook front end that can show the iframe, as long as the **browser can reach the server URL** the iframe is using (and, behind a proxy, the proxy passes WebSockets).

## At a glance

**Audience**
- Wet lab / beginner: read the “Fast path” sections and use the copy/paste snippets.
- Computational users: read “Practical path” + the environment matrix + troubleshooting.
- Developers: read “Deep path” + reference (schemas + endpoints).

**Time**
- First working “select → highlight” loop: ~5–10 minutes
- Full read: ~25–45 minutes (more if you try every troubleshooting recipe)

**Prerequisites**
- `pip install cellucid`
- A notebook environment (classic Jupyter or JupyterLab; VS Code and Colab are untested, see {doc}`03_supported_environments_matrix`)
- A dataset you can view via `show_anndata(...)`

## Fast path (what you can do)

Once you have a `viewer` (from `show_anndata(...)`), you can:

- **Send commands (Python → viewer)**
  - highlight cells by index (`viewer.highlight_cells(...)`)
  - clear highlights (`viewer.clear_highlights()`)
  - color by an obs field (`viewer.set_color_by(...)`)
  - hide/show subsets (`viewer.set_visibility(...)`)
  - reset camera (`viewer.reset_view()`)

- **Receive events (viewer → Python)**
  - `@viewer.on_ready`: viewer is ready for interaction
  - `@viewer.on_selection`: user selects cells
  - `@viewer.on_hover`: user hovers cells (debounced)
  - `@viewer.on_click`: user clicks a cell
  - `@viewer.on_command_error`: a command you sent was rejected by the viewer
  - `@viewer.on_data_mapping`: the viewer resolved or changed its data mapping
  - `@viewer.on_message`: one envelope stream for all nine current event types

- **Pull durable state (session bundles)**
  - `bundle = viewer.get_session_bundle()` (no browser download required)
  - `adata2 = viewer.apply_session_to_anndata(adata)` (apply highlights + user-defined fields back onto AnnData)

If you want the minimal working loop, start with {doc}`02_quickstart_minimal_roundtrip`.

## Practical path (how it works, conceptually)

### One viewer = one server + one iframe

When you run:

```python
from cellucid import show_anndata
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

Cellucid:
1. Starts a local server (usually `http://127.0.0.1:<port>`).
2. Embeds the viewer as an iframe pointing at that same server:
   - `.../?viewer=true&viewerId=<id>&viewerToken=<token>&anndata=true`

The key idea is **same-origin**: the UI and the dataset endpoints share the same origin, which makes the integration much less fragile (no mixed-content, fewer CORS surprises).

### The viewer channel

When the viewer page loads, it opens a WebSocket to `/_cellucid/ws` on the
same server and says hello with its `viewerId` and `viewerToken`. The server
checks the token and answers `welcome`; from then on the socket belongs to
that one viewer.

- When you call `viewer.highlight_cells(...)`, Python sends the command on the
  channel. A command sent before the viewer has connected waits and is sent,
  in order, as soon as it does.
- When a user selects, hovers or clicks, the viewer sends an event on the
  channel (a selection as binary cell indices), and Python triggers your
  registered hooks.

See {doc}`05_architecture_viewer_channel` for the frames and the rules.

## Deep path (implementation pointers)

If you want to trace the implementation:

- Python: `cellucid-python/src/cellucid/jupyter/` (viewer classes + hook registry)
- Python: `cellucid-python/src/cellucid/_viewer_channel.py` (the channel:
  handshake, command outbox, events, session bundles)
- Python: `cellucid-python/src/cellucid/_server_base.py` (verified web-generation delivery)
- Web app: `cellucid/assets/js/data/viewer-channel.js` (the channel) and
  `data/viewer-host-source.js` (ViewerHostDataSource; commands and events)
- Web app: `cellucid/assets/js/app/main.js` (wires session bundle capture + ready event)

## Important notebook constraint: HTTPS / mixed content

If your notebook is served on **HTTPS** (common on JupyterHub), browsers may block iframes that try to load an **HTTP loopback** server directly.

Cellucid uses exactly the `client_server_url=` argument when one is supplied;
otherwise it uses the bound server URL. It does not infer or probe a proxy.
For Jupyter Server Proxy, Colab, or another reverse proxy, obtain the exact
browser-reachable HTTP(S) base URL and pass it to `show_anndata(...)` or the
viewer constructor. See
{doc}`13_security_cors_origins_and_mixed_content` and
{doc}`03_supported_environments_matrix`.

## Edge cases (high-level)

- Multiple viewers: each viewer has a unique `viewerId` and its own channel.
- Large selections: sending huge index lists back and forth can be slow.
- Callbacks run on the viewer channel's thread, not the notebook's: be careful
  with thread-safety when they touch shared objects such as your `adata`.
- Source-restricted environments: viewer/server startup fails directly when it
  cannot establish the configured exact UI generation.

## Troubleshooting (jump table)

- Viewer doesn’t load / blank iframe → {doc}`14_troubleshooting_hooks`
- Events don’t arrive → {doc}`14_troubleshooting_hooks` (“No events arrive”)
- Commands don’t affect the UI → {doc}`06_python_to_frontend_commands`
- HTTPS notebook issues → {doc}`13_security_cors_origins_and_mixed_content`

## Next steps

- Minimal working loop: {doc}`02_quickstart_minimal_roundtrip`
- Commands (Python → viewer): {doc}`06_python_to_frontend_commands`
- Events (viewer → Python): {doc}`07_frontend_to_python_events`
- Robust callback patterns: {doc}`08_writing_robust_callbacks`
- Full reference (schemas + connectivity arguments): {doc}`16_reference`
