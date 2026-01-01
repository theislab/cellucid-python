# Programmatic highlighting + selection callbacks

This tutorial shows how to make Cellucid interactive *as part of a Python workflow*:
- user selects cells in the viewer → Python receives indices (hooks/events)
- Python runs analysis → sends highlights or commands back to the viewer

This is the “bridge” between exploratory visualization and reproducible analysis.

---

## At a glance

**Audience**
- Computational users: primary
- Developers: read the “deep path” and debugging sections

**Time**
- Minimal round-trip (select → print → highlight): ~10 minutes
- Full read + troubleshooting: ~30–60 minutes

**Prerequisites**
- A notebook environment
- `pip install cellucid`
- An `AnnData` (or an export folder from `prepare(...)`)

---

## The mental model (one page)

There are two communication channels:

### Viewer → Python (events/hooks)

The viewer sends HTTP POSTs to the local server:
- `POST /_cellucid/events`

Events include:
- `ready`
- `selection`
- `hover`
- `click`
- `message` (catch-all / custom types)

### Python → Viewer (commands)

Python sends `postMessage(...)` commands into the iframe, for example:
- highlight selected cells
- set color-by field
- reset camera

---

## Minimal round-trip (copy/paste)

### 1) Open the viewer

```python
from cellucid import show_anndata

viewer = show_anndata(adata, height=650)
viewer
```

### 2) Wait until it’s ready (optional but recommended)

```python
viewer.wait_for_ready(timeout=30)
```

### 3) Print selection sizes in Python

```python
@viewer.on_selection
def _handle_selection(event):
    cells = event.get("cells", [])
    print(f"Selected {len(cells)} cells")
```

Now select cells in the viewer (lasso, etc.) and confirm Python prints.

### 4) Highlight the selection back in the UI

```python
@viewer.on_selection
def _highlight_selection(event):
    cells = event.get("cells", [])
    viewer.highlight_cells(cells, color="#ff2d55")
```

---

## Practical patterns (what you’ll actually want in real notebooks)

### Pattern A — “Wait for one selection, then continue”

This is a clean workflow for tutorials and reproducible notebooks:

```python
event = viewer.wait_for_event("selection", timeout=120)
cells = event["cells"]
print("Selected:", len(cells))
viewer.highlight_cells(cells, color="#00c853")
```

Why this is useful:
- you can write notebooks that pause until the user selects something
- you avoid “callback soup” for simple interactive flows

### Pattern B — Turn a selection into an `AnnData` subset safely

```python
event = viewer.wait_for_event("selection", timeout=120)
cells = event.get("cells", [])

subset = adata[cells].copy()
subset
```

```{important}
Selection indices refer to the **row order of the dataset at the time you created the viewer**.

If you reorder or subset `adata` after creating `viewer`, your indices can become semantically confusing.
Best practice: treat `adata` as immutable after opening the viewer, or recreate the viewer after major changes.
```

### Pattern C — Robust callbacks (don’t crash your notebook)

Even though Cellucid tries to log callback errors without hard-crashing, you should still guard analysis code:

```python
@viewer.on_selection
def _safe_handler(event):
    try:
        cells = event.get("cells", [])
        if not cells:
            return
        # ... do analysis ...
    except Exception as e:
        print("Selection handler failed:", e)
```

### Pattern D — Use `viewer.state` for “latest selection” workflows

Cellucid stores the latest event payloads in `viewer.state`:

```python
viewer.state.selection
viewer.state.last_event_type
```

This is helpful when:
- you want to poll from another cell (“what was the last selection?”)
- you want to decouple UI interactions from analysis code

---

## Python → viewer commands (what you can do)

These convenience methods send commands into the iframe:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
viewer.clear_highlights()
viewer.set_color_by("leiden")  # obs field name
viewer.set_visibility([1, 2, 3], visible=False)  # hide
viewer.reset_view()
```

If you need low-level control:

```python
viewer.send_message({"type": "highlight", "cells": [1, 2, 3], "color": "#ff0000"})
```

```{note}
The viewer must be displayed before messages can be delivered. The `show(...)` / `show_anndata(...)` helpers display automatically in notebook contexts.
```

---

## Screenshot placeholder (optional: “selection → highlight” success state)

<!-- SCREENSHOT PLACEHOLDER
ID: python-notebooks-selection-highlight-roundtrip
Suggested filename: highlighting_selection/01_selection-highlight-roundtrip.png
Where it appears: User Guide → Python Package → Notebooks/Tutorials → 23_programmatic_highlighting_and_selection_callbacks.md
Capture:
  - UI location: viewer embedded in notebook + Python output
  - State prerequisites: user selects cells; Python receives event; highlight applied back
  - Action to reach state: make a selection; run highlight callback; confirm highlight color changes
Crop:
  - Include: viewer canvas with highlighted cells
  - Include: notebook output showing “Selected N cells”
Redact:
  - Remove: private dataset names/paths
Alt text:
  - Notebook-embedded viewer with a selection highlighted and Python output confirming selection size.
Caption:
  - Hooks let the viewer send selected cell indices to Python; Python can respond by highlighting those cells back in the UI.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a selection-highlight roundtrip.
:width: 100%

Selection event received in Python, and the same cells highlighted back in the viewer.
```

---

## Edge cases (things that confuse even experienced users)

### Multiple viewers in one notebook

Each viewer has its own `viewerId` token, and events are routed to the correct Python object.

Best practices:
- keep a variable name per viewer (`viewer_a`, `viewer_b`)
- call `viewer.stop()` when you’re done to free ports and avoid confusion

### Long-running analysis inside callbacks

If your callback takes seconds/minutes, the UI may feel laggy.

Options:
- do minimal work inside the callback and enqueue to a background thread/process
- use `wait_for_event(...)` and run analysis in the main flow instead

### Hover events can be very frequent

Treat hover as “streaming telemetry”:
- debounce/throttle your handler logic
- do not run heavy analysis on every hover

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “I select cells but Python receives nothing”

Likely causes:
- viewer never reached `ready`
- browser cannot reach the kernel-side server port (remote kernel without proxy/tunnel)
- a proxy/corporate environment strips requests

How to confirm:
```python
viewer.debug_connection()
```

Fix:
- if on remote Jupyter, enable `jupyter-server-proxy` or use SSH tunneling (see {doc}`22_large_dataset_server_mode_and_lazy_gene_expression`)
- verify `/_cellucid/health` is reachable from the browser

### Symptom: “Python receives events but highlight doesn’t show”

Likely causes:
- you sent highlight before the viewer finished initializing
- your cell indices are wrong (reordered dataset)

Fix:
- wait for ready: `viewer.wait_for_ready()`
- print a few indices and sanity-check by clicking cells

### Symptom: “My notebook becomes slow over time”

Likely causes:
- multiple servers/viewers left running
- accumulating callbacks

Fix:
- call `viewer.stop()` for viewers you no longer use
- clear hooks when needed: `viewer.clear_hooks()`

---

## Next steps

- {doc}`32_session_persistence_and_restoring_analysis_artifacts` (save/restore state; reproducibility)
- {doc}`05_jupyter_embedding_hooks_sessions_gallery` (runnable hook/session notebooks)
- Deep docs: {doc}`../e_jupyter_hooks/index`
