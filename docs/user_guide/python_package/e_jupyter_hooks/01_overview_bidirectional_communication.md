# Overview: bidirectional communication (Python ↔ viewer)

Cellucid’s notebook integration is intentionally **simple and robust**:

- `show_anndata(...)` shows a notebook **widget** (an anywidget) whose output
  holds the Cellucid **web app**, loaded by the browser from
  `https://www.cellucid.com`.
- The widget connects the app to the kernel through the notebook's **own
  connection** (a Jupyter comm): the viewer reads the dataset's bytes from the
  kernel over it, and opens **one channel** on it, on which both directions
  travel:
  1. **Python → viewer**: commands (e.g., highlight, set color-by).
  2. **Viewer → Python**: events (e.g., selection, hover, click).

No port is opened and no server address is involved, so it works the same in
every notebook front end that shows widgets — locally, on a remote JupyterHub,
in VS Code or in Colab. The kernel handles the widget's messages **between
cells**. A browser-tab viewer ({class}`~cellucid.AnnDataViewer`) has the same
API, with the channel on a WebSocket to its local server (`/_cellucid/ws`).

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
- A notebook that shows widgets (JupyterLab is tested; Notebook 7, VS Code, Colab, marimo and Databricks are expected to work, see {doc}`03_supported_environments_matrix`)
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

- **Pull results back (one call)**
  - `viewer.pull(adata)` writes every highlight group, label, cluster, score
    and gene set the viewer holds onto `adata.obs` / `adata.uns`
    ({doc}`../g_api_reference_coverage/api/round_trip`)
  - `viewer.push(adata.obs["leiden"])` and `viewer.push_gene_set(...)` send
    columns and gene sets the other way

- **Capture durable state (session bundles)**
  - `request = viewer.request_session_bundle()`, then `bundle = request.result()`
    in the next cell (no browser download required)

If you want the minimal working loop, start with {doc}`02_quickstart_minimal_roundtrip`.

## Practical path (how it works, conceptually)

### One viewer = one widget (or one browser tab)

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
1. Checks the dataset and the launch choices (an in-memory AnnData is served
   in priority order, made from it as the viewer reads it).
2. Shows a widget whose output frames `https://www.cellucid.com/?viewer=embedded`.
   The page asks the widget for the viewer's id and token and a private
   channel (a `MessagePort`); the address carries no credentials.

A browser-tab viewer (`AnnDataViewer`) instead serves the app and the data
from one local origin and opens
`http://127.0.0.1:<port>/?viewer=true&viewerId=<id>&viewerToken=<token>&anndata=true`.

### The viewer channel

When the viewer page loads, it opens its channel — over the widget in a
notebook, a WebSocket to `/_cellucid/ws` in a browser tab — and says hello
with its `viewerId` and `viewerToken`. Python checks the token and answers
`welcome`; from then on the channel belongs to that one viewer.

- When you call `viewer.highlight_cells(...)`, Python sends the command on the
  channel. A command sent before the viewer has connected waits and is sent,
  in order, as soon as it does.
- When a user selects, hovers or clicks, the viewer sends an event on the
  channel (a selection as binary cell indices), and Python triggers your
  registered hooks.

See {doc}`05_architecture_viewer_channel` for the frames and the rules.

## Deep path (implementation pointers)

If you want to trace the implementation:

- Python: `cellucid-python/src/cellucid/jupyter/` (viewer classes + hook
  registry; `_widget.py` and `static/widget.js` are the notebook widget)
- Python: `cellucid-python/src/cellucid/_viewer_channel.py` (the channel:
  handshake, command outbox, events, session bundles)
- Python: `cellucid-python/src/cellucid/_server_base.py` (verified web-generation delivery)
- Web app: `cellucid/assets/js/data/viewer-embed.js` (the embedded page's
  transport), `data/viewer-channel.js` (the channel) and
  `data/viewer-host-source.js` (ViewerHostDataSource; commands and events)
- Web app: `cellucid/assets/js/app/main.js` (wires session bundle capture + ready event)

## Important notebook constraint: between cells

The kernel handles a widget's messages between cells. Hooks fire and
`viewer.state` updates while the kernel is idle; a call that waits for the
viewer in the cell that asks (`wait_for_ready()`, `wait_for_event()`, a
session request's `result()`) raises at once and says what to do instead; commands
(`push()`, `set_data_mapping()`, …) never wait. See {doc}`03_supported_environments_matrix`.

## Edge cases (high-level)

- Multiple viewers: each viewer has a unique `viewerId` and its own channel.
- Large selections: sending huge index lists back and forth can be slow.
- Callbacks run on the viewer channel's thread, not the notebook's: be careful
  with thread-safety when they touch shared objects such as your `adata`.
- Restricted networks: the browser must load the app from `web_source_url`
  (`https://www.cellucid.com` by default); the output says when it cannot.

## Troubleshooting (jump table)

- Viewer doesn’t load / empty output → {doc}`14_troubleshooting_hooks`
- Events don’t arrive → {doc}`14_troubleshooting_hooks` (“No events arrive”)
- Commands don’t affect the UI → {doc}`06_python_to_frontend_commands`
- Origins and browser security → {doc}`13_security_cors_origins_and_mixed_content`

## Next steps

- Minimal working loop: {doc}`02_quickstart_minimal_roundtrip`
- Commands (Python → viewer): {doc}`06_python_to_frontend_commands`
- Events (viewer → Python): {doc}`07_frontend_to_python_events`
- Robust callback patterns: {doc}`08_writing_robust_callbacks`
- Full reference (schemas + connectivity arguments): {doc}`16_reference`
