# Overview: bidirectional communication (Python ↔ viewer)

Cellucid’s notebook integration is intentionally **simple and robust**:

- **Python** starts a small **local HTTP server** that serves your dataset.
- The notebook displays the Cellucid **web app UI** in an **iframe**.
- You get **two channels** for communication:
  1. **Python → viewer**: commands via **`postMessage`** (e.g., highlight, set color-by).
  2. **Viewer → Python**: events via **HTTP POST** to `/_cellucid/events` (e.g., selection, hover, click).

This design avoids “special Jupyter extensions” and works across notebook frontends (classic Jupyter, JupyterLab, VSCode notebooks, Colab) as long as the **browser can reach the server URL** the iframe is using.

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
- A notebook environment (classic, JupyterLab, VSCode, or Colab)
- A dataset you can view via `show_anndata(...)` or `show(...)`

## Fast path (what you can do)

Once you have a `viewer` (from `show_anndata(...)` or `show(...)`), you can:

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
  - `@viewer.on_message`: raw debugging stream (all event types)

- **Pull durable state (session bundles)**
  - `bundle = viewer.get_session_bundle()` (no browser download required)
  - `adata2 = viewer.apply_session_to_anndata(adata)` (apply highlights + user-defined fields back onto AnnData)

If you want the minimal working loop, start with {doc}`02_quickstart_minimal_roundtrip`.

## Practical path (how it works, conceptually)

### One viewer = one server + one iframe

When you run:

```python
from cellucid import show_anndata
viewer = show_anndata(adata)
```

Cellucid:
1. Starts a local server (usually `http://127.0.0.1:<port>`).
2. Embeds the viewer as an iframe pointing at that same server:
   - `.../?jupyter=true&viewerId=<id>&viewerToken=<token>[&anndata=true]`

The key idea is **same-origin**: the UI and the dataset endpoints share the same origin, which makes the integration much less fragile (no mixed-content, fewer CORS surprises).

### Channel A: Python → viewer (commands)

When you call `viewer.highlight_cells(...)`, Python injects a small JavaScript snippet into the notebook output that uses `postMessage(...)` to send a command into the iframe.

Commands are authenticated with a per-viewer token (`viewerToken`) to avoid accepting arbitrary parent-frame messages.

### Channel B: viewer → Python (events)

When a user selects / hovers / clicks in the UI, the viewer sends an HTTP POST:

```text
POST /_cellucid/events
Content-Type: application/json
{ "type": "selection", "viewerId": "...", ...payload... }
```

The Python server routes that JSON event to the right `viewer` instance using `viewerId`, then triggers your registered hooks.

## Deep path (implementation pointers)

If you want to trace the implementation:

- Python: `cellucid-python/src/cellucid/jupyter.py` (viewer classes + hook registry)
- Python: `cellucid-python/src/cellucid/_server_base.py` (`/_cellucid/events`, session bundle upload, hosted-asset proxy)
- Web app: `cellucid/assets/js/data/jupyter-source.js` (JupyterBridgeDataSource; message handling + event POST)
- Web app: `cellucid/assets/js/app/main.js` (wires session bundle capture + ready event)

## Important notebook constraint: HTTPS / mixed content

If your notebook is served on **HTTPS** (common on JupyterHub), browsers may block iframes that try to load an **HTTP loopback** server directly.

Cellucid handles this by preferring a notebook proxy URL when available (see {doc}`13_security_cors_origins_and_mixed_content` and {doc}`03_supported_environments_matrix`):
- **Jupyter Server Proxy** (`.../proxy/<port>/...`) for JupyterHub / many hosted setups
- **Colab port proxy** for Google Colab
- Manual override: `CELLUCID_CLIENT_SERVER_URL` (advanced)

## Edge cases (high-level)

- Multiple viewers: each viewer has a unique `viewerId`; events route by ID.
- Large selections: sending huge index lists back and forth can be slow.
- Backed AnnData (`.h5ad` backed mode): be careful with thread-safety in callbacks.
- Offline environments: first-time UI asset download may fail without a cached copy (see troubleshooting).

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
- Full reference (schemas + env vars): {doc}`16_reference`

