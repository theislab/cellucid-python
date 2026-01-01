# What is a “viewer” object?

**Audience:** everyone (wet lab → power user)  
**Time:** 10–15 minutes  
**Goal:** understand what `viewer = show(...)` returns, what it controls, and how to clean it up safely.

---

## Fast path (wet lab / non-technical): what you should remember

A **viewer** is the interactive **Cellucid web app UI** running in your browser (or inside a notebook cell).

A **viewer object** is the **Python handle** you use to:
- open that UI,
- send commands to it (highlight, color-by, reset camera),
- receive events back (selection/hover/click),
- and stop the underlying local server when you’re done.

You can think of it like: “a remote control for the web app”.

---

## Practical path (computational users): there are two viewer types

Cellucid has two main notebook viewer classes:

1) **`CellucidViewer`** (export-first; fastest UI for large datasets)  
   Use when you already have an export folder created by `prepare(...)`.

   ```python
   from cellucid import show

   viewer = show("./my_export", height=650)
   viewer.wait_for_ready(timeout=60)
   ```

2) **`AnnDataViewer`** (no-export convenience; great for exploration)  
   Use when you want to view an AnnData object (or `.h5ad` / `.zarr`) directly.

   ```python
   from cellucid import show_anndata

   viewer = show_anndata(adata, height=650)
   viewer.wait_for_ready(timeout=60)
   ```

Under the hood, both are subclasses of `cellucid.jupyter.BaseViewer`.

```{note}
You can also run Cellucid without notebooks using the CLI (`cellucid serve ...`) or Python (`serve(...)` / `serve_anndata(...)`).
In that mode you typically **don’t** have a “viewer object”; you just open the printed URL in a browser.
```

---

## Deep path: what a viewer object actually owns

When you create a viewer, Python starts a **local HTTP server** (by default on `127.0.0.1`):
- it serves the dataset (export folder or AnnData-backed endpoints),
- it serves the web app UI via a “hosted-asset proxy” cache,
- it accepts event POSTs from the web app (`/_cellucid/events`),
- and it can accept “session bundle” uploads (`/_cellucid/session_bundle`) for the no-download workflow.

The viewer object also owns:
- a **viewer ID** (`viewerId`) and **token** (`viewerToken`) used to route/authenticate messages,
- a **URL** you can open manually (`viewer.viewer_url`),
- hook registries for `on_selection`, `on_hover`, `on_click`, `on_ready`, `on_message`,
- a small, thread-safe **latest-event snapshot** (`viewer.state`),
- lifecycle methods (`viewer.display()`, `viewer.stop()`).

---

## Lifecycle: create → display → interact → stop

### 1) Create the viewer

```python
from cellucid import show_anndata
viewer = show_anndata("data.h5ad", height=650)
```

### 2) Display it (if you’re in a notebook)

Most of the time it auto-displays. If you created it with `auto_open=False` or you’re in a non-notebook context:

```python
viewer.display()
```

If you are not in Jupyter, `display()` prints a link:

```text
Not in Jupyter environment. Open manually: http://127.0.0.1:8765/?jupyter=true&...
```

### 3) Wait for readiness (recommended before automation)

```python
viewer.wait_for_ready(timeout=60)
```

This is especially important before:
- requesting session bundles,
- relying on hooks,
- sending programmatic commands you expect to “stick” immediately.

### 4) Interact (hooks + commands)

Hooks (frontend → Python):

```python
@viewer.on_selection
def handle_selection(event):
    print("Selected cells:", len(event["cells"]))
```

Commands (Python → frontend):

```python
viewer.set_color_by("cell_type")
viewer.highlight_cells([0, 1, 2], color="#ff00aa")
```

### 5) Stop (cleanup)

Each viewer starts a background server. Stop it when you’re done:

```python
viewer.stop()
```

```{note}
`viewer.stop()` also tries to “freeze” the web app so the notebook output stays visually stable, but it becomes non-interactive.
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
viewer.wait_for_ready(timeout=60)

print(viewer.state.ready)       # last ready event (or None)
print(viewer.state.selection)   # last selection event (or None)
print(viewer.state.hover)       # last hover event (or None)
print(viewer.state.click)       # last click event (or None)
print(viewer.state.last_event)  # last event payload (any type)
```

Use `viewer.state` when you want “pull” style logic. Use hooks when you want “push/reactive” logic.

---

## Multi-viewer mental model (ports + IDs)

It is normal to have multiple viewers in one notebook.

- Each viewer gets its own `viewerId`/`viewerToken`.
- Each viewer starts its own server (usually a different port: `8765`, `8766`, `8767`, …).
- Events from the browser include `viewerId` so Python can route them correctly.

```{important}
Selections/highlights are sent as **cell indices** (row positions). If your dataset changes row order, indices refer to different cells.
For durable state transfer, use session bundles and stable dataset identity practices: {doc}`04_dataset_identity_and_reproducibility`.
```

---

## Screenshot: what “a viewer in a notebook” looks like

If you want to add a screenshot to orient beginners, capture the *first successful embedded viewer* output.

<!-- SCREENSHOT PLACEHOLDER
ID: python-viewer-object-embedded
Suggested filename: web_app/embedded-viewer-from-python.png
Where it appears: Python Package → Concepts & Mental Models → What is a “viewer” object?
Capture:
  - UI location: Jupyter notebook output cell
  - State prerequisites: any dataset successfully displayed (AnnData or export folder)
  - Action to reach state: run `viewer = show_anndata(...)` or `viewer = show(...)`
Crop:
  - Include: the embedded iframe + a small amount of surrounding notebook context (so it’s obvious this is “in a notebook”)
  - Exclude: personal file paths, usernames, dataset private identifiers
Redact:
  - Remove: private dataset names displayed in the viewer title bar (if any)
Annotations:
  - Callouts: #1 “this is the embedded viewer”, #2 “this cell started the viewer”
Alt text:
  - Jupyter notebook cell output containing an embedded Cellucid viewer.
Caption:
  - The Cellucid web app is embedded as an iframe; the Python `viewer` object controls it and receives events from it.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a Cellucid viewer embedded in a Jupyter notebook.
:width: 100%

Cellucid embedded in a notebook: the UI is the web app; the Python `viewer` object is the handle that controls it.
```

---

## Edge cases and footguns

### “I re-ran a cell and now I have multiple servers”

This is expected. Each run can create a new viewer + new port.

Mitigations:
- call `viewer.stop()` in a `finally:` block,
- or call `cleanup_all()` occasionally during exploratory work.

### “I created a viewer but didn’t see anything”

- In notebooks: call `viewer.display()`.
- Outside notebooks: open `viewer.viewer_url` in a browser.

### “Hooks work sometimes, then stop”

Common causes:
- the notebook output iframe was destroyed (cell cleared, output collapsed, notebook reloaded),
- the server port changed and the browser is still pointing to the old port,
- mixed-content/proxy issues in remote notebook environments.

Start with: {doc}`08_debugging_mental_model_where_to_look`.

---

## Troubleshooting

### Symptom: “The viewer area is blank”

Likely causes (ordered):
1) The viewer UI assets could not be loaded (offline/no cache).
2) Notebook is served from HTTPS, but the viewer tries to load an HTTP loopback URL (mixed content).
3) The local server is not reachable from the browser (remote kernel, missing proxy/tunnel).

How to confirm:
- Run `report = viewer.debug_connection()` and inspect:
  - `server_health`, `viewer_index_probe`, `frontend_roundtrip`, `frontend_console`.
- In browser devtools console, look for blocked network requests.

Fix:
- If offline: run once while online to populate the web UI cache, then retry.
- If remote: use SSH tunneling or `jupyter-server-proxy` (see {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`).

### Symptom: “`viewer.stop()` doesn’t seem to do anything”

Notes:
- `stop()` is best-effort. The notebook output may still show the last rendered frame.
- The key effect is that the underlying server stops and hooks stop firing.

How to confirm:
- Open `viewer.server_url + "/_cellucid/health"` in a browser; it should stop responding after `stop()`.

---

## Next steps

- Learn how the pieces connect: {doc}`02_data_flows`
- Learn hooks in depth: {doc}`../e_jupyter_hooks/index`
- Learn viewing modes: {doc}`../d_viewing_apis/index`
