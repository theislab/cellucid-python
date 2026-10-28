# What is a “viewer” object?

**Audience:** everyone (wet lab → power user)  
**Time:** 10–15 minutes  
**Goal:** understand what `viewer = show_anndata(...)` returns, what it controls, and how to clean it up safely.

---

## Fast path (wet lab / non-technical): what you should remember

A **viewer** is the interactive **Cellucid web app UI** running in your browser (or inside a notebook cell).

A **viewer object** is the **Python handle** you use to:
- show that UI,
- send commands to it (highlight, color-by, reset camera, push a column),
- receive events back (selection/hover/click),
- pull what you made in it (highlights, labels, scores, gene sets) back onto your AnnData,
- and stop it when you’re done.

You can think of it like: “a remote control for the web app”.

---

## Practical path (computational users): the viewer types

There are two viewer classes with one API:

- **`cellucid.Widget`** (`cellucid.jupyter.Widget`) — what
  `show_anndata(...)` creates and shows: a notebook widget in the cell's
  output. The browser loads the app from `https://www.cellucid.com`, and the
  app reads the data from the kernel through the notebook's own connection.
  No port, no server.
- **`cellucid.AnnDataViewer`** — the same API for a viewer in a browser tab of
  its own: it starts a local data server and gives you the tab's address
  (`viewer.viewer_url`). Use it from a script, or through an SSH tunnel.

Both view an AnnData or MuData object, or any file or folder Cellucid reads,
by path (`.h5ad`, AnnData or MuData Zarr, `.h5mu`, `.loom`, a 10x `.h5` file
or Cell Ranger `outs/` folder, a MatrixMarket folder, or an R object
`.rds`/`.rda`/`.RData`). A file written by `prepare(...)` names itself, so it
needs no `dataset_name`/`dataset_id`:

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=650,
    dataset_name="Example",
    dataset_id="example",
)
```

```{note}
You can also run Cellucid without notebooks using the CLI (`cellucid serve ...`) or Python (`serve_anndata(...)`).
In that mode you typically **don’t** have a “viewer object”; you just open the printed URL in a browser.
```

---

## Deep path: what a viewer object actually owns

A notebook widget:
- answers the viewer's reads of the dataset (the bytes of one file or folder,
  which the web app reads itself, with range requests) through the widget,
- accepts the viewer's channel on the same connection, on which Python sends
  commands and the web app sends events, selections, results and “session
  bundles” for the no-download workflow.

A browser-tab viewer does the same through a **local HTTP server** (by default
on `127.0.0.1`), which also serves one exact inventory-verified web
generation and accepts the channel as a WebSocket (`/_cellucid/ws`).

Both also own:
- a **viewer ID** (`viewerId`) and **token** (`viewerToken`) the viewer page uses to open its channel,
- hook registries for `on_selection`, `on_hover`, `on_click`, `on_ready`,
  `on_command_error` (a command this notebook sent was refused by the viewer),
  `on_data_mapping`, and `on_message`,
- a small, thread-safe **latest-event snapshot** (`viewer.state`) and every
  result the viewer holds (`viewer.results`, `viewer.pull(adata)`),
- lifecycle methods (`viewer.display()`, `viewer.stop()`).

---

## Lifecycle: create → interact → pull → stop

### 1) Create (and show) the viewer

```python
from cellucid import show_anndata
viewer = show_anndata(
    "data.h5ad",
    height=650,
    dataset_name="Example",
    dataset_id="example",
)
```

In JupyterLab, Notebook 7, VS Code, Colab and Databricks it is shown at once;
in marimo, make `viewer` the cell's output. To show it again in another
output, end a cell with `viewer` (the newest output takes the viewer over).
Outside a notebook `show_anndata()` raises and names the browser-tab
alternatives.

### 2) The viewer loads after the cell ends

The kernel handles the widget's messages **between cells**. So the viewer
loads once the cell that showed it finishes, hooks run while the kernel is
idle, and a call that waits for the viewer in the same cell
(`wait_for_ready()`, `wait_for_event()`, a session request's `result()`)
raises at once and says what to do instead. Commands -- `push()`,
`set_data_mapping()`, highlights -- never wait: they go as soon as the viewer
can take them, and the viewer's answer arrives as live state. In a later
cell:

```python
viewer.state.ready     # None until the dataset has opened
```

### 3) Interact (hooks + commands)

Hooks (frontend → Python):

```python
@viewer.on_selection
def handle_selection(event):
    print("Selected cells:", len(event["cells"]))
```

Commands (Python → frontend) are sent at once, or delivered in order when the
viewer connects:

```python
viewer.set_color_by("cell_type")
viewer.highlight_cells([0, 1, 2], color="#ff00aa")
```

### 4) Pull results back

```python
summary = viewer.pull(adata)   # highlights, labels, scores -> adata.obs
summary.columns
```

See {doc}`05_sessions_to_anndata_bridge` and
{doc}`../g_api_reference_coverage/api/round_trip`.

### 5) Stop (cleanup)

```python
viewer.stop()
```

```{note}
`viewer.stop()` freezes the web app on its last view (the output stays
visible but non-interactive), stops answering it, and removes the temporary
`.h5ad` an in-memory AnnData was written to. A browser-tab viewer also stops
its server.
```

To stop everything created in the current kernel:

```python
from cellucid.jupyter import cleanup_all
cleanup_all()
```

---

## `viewer.state`: the simplest way to “peek” at what just happened

`viewer.state` is a small snapshot that updates when events arrive.

```python
print(viewer.state.ready)       # last ready event (or None)
print(viewer.state.selection)   # last selection event (or None)
print(viewer.state.hover)       # last hover event (or None)
print(viewer.state.click)       # last click event (or None)
print(viewer.state.last_event)  # last event payload (any type)
```

Use `viewer.state` (or `viewer.selection`, `viewer.results`) when you want
“pull” style logic in the next cell. Use hooks when you want “push/reactive”
logic.

---

## Multi-viewer mental model (IDs)

It is normal to have multiple viewers in one notebook.

- Each viewer gets its own `viewerId`/`viewerToken`.
- Each viewer's page opens its own channel and authenticates it with that
  viewer's `viewerId`/`viewerToken`, so its events reach only its viewer object.
- Browser-tab viewers each start their own server. With the default
  `port=None` it serves on 8765 -- the one local port Google sign-in accepts
  -- or, when another viewer holds it, on the next free port, which it prints
  (as the R package's `cellucid_view()` does outside a notebook).

```{important}
Selections/highlights are sent as **cell indices** (row positions). If your dataset changes row order, indices refer to different cells.
For durable state transfer, use `viewer.pull(adata)` or session bundles, and stable dataset identity practices: {doc}`04_dataset_identity_and_reproducibility`.
```

---

## Edge cases and footguns

### “I re-ran a cell and now I have several viewers”

This is expected. Each run creates a new viewer.

Mitigations:
- call `viewer.stop()` when done,
- or call `cleanup_all()` occasionally during exploratory work.

### “I created a viewer but didn’t see anything”

- In a notebook: the front end must show widgets ({doc}`../a_landing_pages/03_compatibility_matrix_must_be_explicit`).
- Outside notebooks: use `AnnDataViewer(...)` and open `viewer.viewer_url` in a browser.

### “Hooks work sometimes, then stop”

Common causes:
- a cell is running (events wait for the kernel to be idle),
- the output was removed or the notebook reloaded (re-run the cell),
- the kernel restarted (the old output keeps its last view and says the
  notebook stopped it).

Start with: {doc}`08_debugging_mental_model_where_to_look`.

---

## Troubleshooting

### Symptom: “The viewer area is empty”

Likely causes (ordered):
1) The notebook front end does not show widgets (`anywidget` missing where
   the notebook server runs).
2) The browser cannot load the app from `https://www.cellucid.com` (the output
   says so).

How to confirm:
- Run `viewer.debug_connection()` in a new cell and inspect
  `connected_outputs`, `channel_connected` and `state.ready`.
- In browser devtools console, look for blocked network requests.

Fix:
- Install `anywidget` where the notebook server runs, and reload.
- Allow `https://www.cellucid.com`, or host a copy and pass `web_source_url=`.

### Symptom: “`viewer.stop()` doesn’t seem to do anything”

Notes:
- `stop()` reports freeze or shutdown failures. After a successful stop, the
  notebook output still shows the last rendered frame.
- The key effect is that hooks stop firing and the viewer is no longer answered.

---

## Next steps

- Learn how the pieces connect: {doc}`02_data_flows`
- Learn hooks in depth: {doc}`../e_jupyter_hooks/index`
- Learn viewing modes: {doc}`../d_viewing_apis/index`
