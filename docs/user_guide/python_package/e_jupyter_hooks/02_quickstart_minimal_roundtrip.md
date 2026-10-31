# Quickstart: minimal round-trip (select → highlight)

Goal: in one notebook, get a working loop where you **select cells in the viewer** and Python **highlights those same cells back** in the viewer.

This is the fastest sanity check that:
- the widget shows the viewer and the viewer reads the dataset,
- viewer → Python events are arriving,
- Python → viewer commands are being delivered.

## At a glance

**Audience**
- Wet lab / beginner: follow the “Copy/paste” cells and focus on “What success looks like”.
- Computational users: also read “Common pitfalls” and “Troubleshooting”.
- Developers: jump to {doc}`05_architecture_viewer_channel` and {doc}`16_reference`.

**Time**: ~5–10 minutes

**Prerequisites**
- A notebook that shows widgets (JupyterLab is tested; Notebook 7, VS Code, Colab, marimo and Databricks are expected to work, see {doc}`03_supported_environments_matrix`)
- `pip install cellucid`
- An `AnnData` object named `adata`, or a path to a file or folder Cellucid reads (`.h5ad`, AnnData Zarr, `.h5mu`, `.loom`, 10x, MatrixMarket, `.rds`)

```{note}
If you have a file written by `cellucid.prepare()`, `show_anndata("my_dataset.h5ad", height=600)`
opens it under the name and id it declares.
```

## Step 1 — Create and display a viewer

### Option A (most common): show an AnnData

Run this in a notebook cell:

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer
```

### Option B: show a Cellucid-ready AnnData file

```python
from cellucid import show_anndata

viewer = show_anndata("my_dataset.h5ad", height=600)  # written by cellucid.prepare()
viewer
```

## Step 2 — Check that the viewer reported “ready”

The viewer loads once the cell that showed it finishes. When it shows the
data, run this in a new cell:

```python
print("Ready payload:", viewer.state.ready)
```

`None` means it has not opened the dataset yet; run the cell again in a
moment. (A notebook widget hears the viewer only between cells, so
`viewer.wait_for_ready()` in the same cell would raise at once; a hook,
`@viewer.on_ready`, runs as soon as it is ready.)

## Step 3 — Make a selection in the viewer UI

In the embedded viewer:
- use any selection tool (lasso is usually easiest), or
- click cells (depending on your workflow), or
- use the UI’s selection/highlight features.

For a deep dive on selection tools, see the web app guide:
- {doc}`../../web_app/f_highlighting_selection/02_selection_tools_document_each_tool`

## Step 4 — Read the selection in Python

After you made the selection, run this in a new cell:

```python
selection = viewer.selection          # None before the first selection
cells = selection.cells               # rows of obs (0-based), ascending
print(f"Selected {len(selection)} cells (first 10):", cells[:10])
print(selection.obs_names[:10])       # their obs names
```

The selection reached Python while the kernel was idle, so the next cell sees
it at once.

```{important}
Cell indices are the cells' **original positions** (0-based): a cell's row in
the object or file as it was shown, in its original order.

- For an object shown in its own order, that is its row as it was when shown;
  `selection.obs_names` are the names it had then, whatever is done to it
  after.
- An object read from a Cellucid-ready file in priority order (`anndata.read_zarr`
  of `prepare()`'s output) holds its rows in another order: the viewer still
  names cells by original position, `selection.obs_names` name the right
  cells, and `viewer.pull(adata)` / `viewer.subset(adata)` refuse that object
  until `cellucid.restore_original_order(adata)` puts its cells back.
- To take the cells themselves, use `viewer.subset(adata)` or
  `adata[selection.obs_names]`; an object reordered or subset after it was
  shown is refused by `pull` and `subset` rather than written in the wrong rows.
```

## Step 5 — Highlight the same cells back in the viewer

Run:

```python
if cells:
    viewer.highlight_cells(cells, color="#ff00ff")
```

You should see the selection become highlighted (exact visual style depends on the web app UI).

```{figure} ../../../_static/screenshots/highlighting_selection/highlighting-selected-page.png
:alt: The Highlighting panel with one page renamed "Left lobe" carrying 826 cells, the summary "826 cells highlighted", and one group row labelled "Lasso (826 cells) (826)".
:width: 532px

Highlighted cells are stored on a named page that can be used by analysis and session workflows.
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
Hook callbacks run when the viewer's event reaches Python on its channel, on
that channel's background thread, not in a notebook cell.

- Keep callbacks fast and defensive.
- If you do heavy work (Scanpy, model inference), consider queueing the work and handling it outside the callback.

See {doc}`08_writing_robust_callbacks`.
```

## Common pitfalls (read before you debug)

- **Nothing prints in Step 4**: you haven’t made a selection that triggers a `selection` event yet.
- **Selection indices don’t match your `adata`**: your AnnData row order changed after the viewer started.
- **The output is empty or says the app did not load**: the notebook does not
  show widgets (`anywidget` missing where the notebook server runs), or the
  browser cannot reach `https://www.cellucid.com` (see troubleshooting).

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
- Session bundles (no-download): {doc}`10_session_bundles_request_session_bundle`
