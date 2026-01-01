# Quickstart: minimal round-trip (select → highlight)

Goal: in one notebook, get a working loop where you **select cells in the viewer** and Python **highlights those same cells back** in the viewer.

This is the fastest sanity check that:
- the iframe loads,
- the server is reachable,
- viewer → Python events are arriving,
- Python → viewer commands are being delivered.

## At a glance

**Audience**
- Wet lab / beginner: follow the “Copy/paste” cells and focus on “What success looks like”.
- Computational users: also read “Common pitfalls” and “Troubleshooting”.
- Developers: jump to {doc}`05_architecture_message_routing_http_vs_postmessage` and {doc}`16_reference`.

**Time**: ~5–10 minutes

**Prerequisites**
- A Jupyter environment (classic/JupyterLab/VSCode/Colab)
- `pip install cellucid`
- An `AnnData` object named `adata` (or a path to `.h5ad`/`.zarr`)

```{note}
If you don’t have an `AnnData` handy, you can still follow this guide using a pre-exported directory and `cellucid.show("./export_dir")`.
```

## Step 1 — Create and display a viewer

### Option A (most common): show an AnnData

Run this in a notebook cell:

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=600)
viewer
```

### Option B: show a pre-exported dataset directory

```python
from cellucid import show

viewer = show("./exports/my_dataset", height=600)
viewer
```

<!-- SCREENSHOT PLACEHOLDER
ID: jupyter-hooks-quickstart-viewer-visible
Suggested filename: jupyter_hooks/01_quickstart-viewer-visible.png
Where it appears: Python Package → Jupyter Hooks → Quickstart → Step 1
Capture:
  - UI location: notebook output cell (iframe embed)
  - State prerequisites: the viewer loads successfully (no error page)
  - Action to reach state: run the Step 1 cell and wait for the UI to render
Crop:
  - Include: the embedded viewer + a small amount of notebook context (so beginners recognize it’s “inside the notebook”)
  - Exclude: personal paths, usernames, dataset identifiers if private
Redact:
  - Remove: any sensitive dataset name shown in the top-left title area (if present)
Annotations:
  - Callouts: #1 the viewer iframe, #2 the notebook cell that created `viewer`
Alt text:
  - Notebook cell output showing the Cellucid viewer embedded as an iframe.
Caption:
  - The Cellucid web app UI embedded inside a notebook; once it loads, Python can send commands and receive events.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Cellucid viewer embedded in a notebook output cell.
:width: 100%

Cellucid appears as an interactive iframe inside the notebook output.
```

## Step 2 — Wait until the viewer reports “ready”

Run this in a new cell:

```python
viewer.wait_for_ready(timeout=60)
print("Ready payload:", viewer.state.ready)
```

Why this matters:
- It confirms the viewer finished loading.
- It makes later steps more reliable (you won’t race against initialization).

## Step 3 — Make a selection in the viewer UI

In the embedded viewer:
- use any selection tool (lasso is usually easiest), or
- click cells (depending on your workflow), or
- use the UI’s selection/highlight features.

For a deep dive on selection tools, see the web app guide:
- {doc}`../../web_app/f_highlighting_selection/02_selection_tools_document_each_tool`

## Step 4 — Block until Python receives the selection event

Run this in a new cell. It will **wait** until you make a selection in the viewer:

```python
event = viewer.wait_for_event("selection", timeout=None)  # wait forever
cells = event.get("cells", [])
print(f"Selected {len(cells)} cells (first 10):", cells[:10])
```

```{important}
Cell indices are **0-based row indices** in the dataset that Cellucid is currently serving.

- For `show_anndata(adata)`, indices refer to the current `adata` row order.
- If you reorder/subset `adata` after creating the viewer, indices will no longer match what you expect.
```

## Step 5 — Highlight the same cells back in the viewer

Run:

```python
if cells:
    viewer.highlight_cells(cells, color="#ff00ff")
```

You should see the selection become highlighted (exact visual style depends on the web app UI).

<!-- SCREENSHOT PLACEHOLDER
ID: jupyter-hooks-quickstart-highlight-roundtrip
Suggested filename: jupyter_hooks/02_quickstart-selection-highlight-roundtrip.png
Where it appears: Python Package → Jupyter Hooks → Quickstart → Step 5
Capture:
  - UI location: embedded viewer
  - State prerequisites: a visible selection exists, and `viewer.highlight_cells(...)` was run successfully
  - Action to reach state:
      1) select a distinct group of cells in the viewer (lasso)
      2) run Step 4 to confirm Python received indices
      3) run Step 5 to highlight them
Crop:
  - Include: enough of the viewer UI to show the highlighted points clearly, plus the legend/indicator if it changes
  - Exclude: unrelated notebook cells, personal information
Redact:
  - Remove: private dataset names
Annotations:
  - Callouts: #1 the selected region, #2 highlighted points after roundtrip
Alt text:
  - Embedded viewer showing a group of cells highlighted after Python sends a highlight command.
Caption:
  - Round-trip sanity check: a selection event arrived in Python and the same indices were highlighted back in the viewer.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a selection highlighted back in the viewer from Python.
:width: 100%

After Python receives selected indices, it can send a highlight command back to the viewer.
```

## Optional: “reactive” style with callbacks

If you prefer event-driven code, you can register hooks:

```python
@viewer.on_selection
def handle_selection(event):
    cells = event.get("cells", [])
    print("Selected:", len(cells))
```

```{important}
Hook callbacks run when the viewer POSTs an event to the local server.
In practice, that means callbacks can run in a background server thread.

- Keep callbacks fast and defensive.
- If you do heavy work (Scanpy, model inference), consider queueing the work and handling it outside the callback.

See {doc}`08_writing_robust_callbacks`.
```

## Common pitfalls (read before you debug)

- **Nothing prints in Step 4**: you haven’t made a selection that triggers a `selection` event yet.
- **Selection indices don’t match your `adata`**: your AnnData row order changed after the viewer started.
- **Viewer is blank or shows an error page**: you likely have a notebook proxy / mixed-content / offline UI cache issue (see troubleshooting).

## Troubleshooting

Start here:

```python
report = viewer.debug_connection()
report
```

Then use the symptom-based guide:
- {doc}`14_troubleshooting_hooks`

## Next steps

- Commands catalog (Python → viewer): {doc}`06_python_to_frontend_commands`
- Events catalog (viewer → Python): {doc}`07_frontend_to_python_events`
- Viewer state + waiting: {doc}`09_viewer_state_and_wait_for_event`
- Session bundles (no-download): {doc}`10_session_bundles_get_session_bundle`

