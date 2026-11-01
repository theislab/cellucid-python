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
- An `AnnData` (or a Cellucid-ready AnnData file from `prepare(...)`)

---

## The mental model (one page)

Both directions travel on one channel, the viewer channel. In a notebook
(`show_anndata()`, a widget) it is carried by the notebook's own connection;
in a browser tab (`AnnDataViewer`) it is a WebSocket the viewer opens to this
machine (`/_cellucid/ws`). The notebook hands a widget's messages to the
kernel **between cells**, so events reach Python, and hooks run, while the
kernel is idle.

### Viewer → Python (events/hooks)

The viewer sends events on the channel: JSON text frames, and selections as
binary cell indices.

Events include:
- `ready`
- `selection`
- `hover`
- `click`
- `command_error` (a command from Python the viewer refused)
- `data_mapping` (what the viewer plots and reads)

the results (`result`, `result_removed`, `results`: every highlight group,
label field, score and gene set, which `viewer.pull(adata)` writes onto `obs`),
and the replies `pong`, `debug_snapshot` and `session_bundle` to a browser-tab
viewer's `debug_connection()` and to `request_session_bundle()`.

`message` is not one of them. It is a hook that receives every event the
viewer sends and nothing else — there is no custom or app-specific event type. An
inbound event whose `type` this `cellucid` does not declare is rejected before
any hook runs, so a viewer cannot invent one.

### Python → Viewer (commands)

Python sends commands on the channel, for example:
- highlight selected cells
- set color-by field
- reset camera

---

## Minimal round-trip (copy/paste)

### 1) Open the viewer

```python
from cellucid import show_anndata

viewer = show_anndata(
    adata,
    height=650,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer
```

The viewer appears when this cell finishes. Run the next cells one at a time.

### 2) Print selection sizes in Python

```python
@viewer.on_selection
def _handle_selection(event):
    cells = event["cells"]
    print(f"Selected {len(cells)} cells")
```

Now select cells in the viewer (lasso, etc.) and confirm Python prints.

### 3) Highlight the selection back in the UI

```python
@viewer.on_selection
def _highlight_selection(event):
    cells = event["cells"]
    if cells:
        viewer.highlight_cells(cells, color="#ff2d55")
    else:
        viewer.clear_highlights()
```

---

## Practical patterns (what you’ll actually want in real notebooks)

### Pattern A — “Select, then run the next cell”

The cleanest workflow for tutorials and reproducible notebooks needs no
callback: select in the viewer, then read the selection in the next cell.

```python
selection = viewer.selection          # None before the first selection
print("Selected:", len(selection), "cells with the", selection.source)
viewer.highlight_cells(selection.obs_names, color="#00c853")
```

`viewer.selection` carries the cells as rows of `obs` (`selection.cells`)
and, when the kernel can read them, their obs names (`selection.obs_names`).

```{note}
A notebook cell cannot wait for a selection (`viewer.wait_for_event("selection")`):
the widget hears the viewer only between cells, so it raises at once and says
so. In a browser-tab viewer (`AnnDataViewer`), whose server answers on its own
threads, `wait_for_event` does wait.
```

### Pattern B — Turn a selection into an `AnnData` subset safely

```python
subset = viewer.subset(adata)                     # the last selection
group = viewer.subset(adata, "Lasso selection")   # a highlight group, by label or id
```

To keep every highlight group and label field as `obs` columns instead, pull
them all at once: `viewer.pull(adata)` ({doc}`42_pull_everything_back`).

```{important}
Selection indices refer to the **row order of the dataset at the time you created the viewer**.

If you reorder or subset `adata` after creating `viewer`, your indices can become semantically confusing.
Best practice: treat `adata` as immutable after opening the viewer, or recreate the viewer after major changes.
```

### Pattern C — Validate callback input and surface failures

Hooks run on the viewer channel's thread, not in your cell. An exception a hook
raises is logged (logger `cellucid.viewer_channel`, "A viewer event hook
failed") and the channel stays open. The event itself is not lost: it is
recorded in `viewer.state` (and, in a browser-tab viewer, handed to
`wait_for_event(...)`) before any hook runs, and the other hooks still run for it. Validate the exact input, and
record failures where your notebook can see them:

```python
@viewer.on_selection
def _validated_handler(event):
    cells = event["cells"]
    if type(cells) is not list or any(type(index) is not int for index in cells):
        raise TypeError("selection cells must be a list of native integers")
    if cells:
        # ... do analysis; an exception here is logged, not raised in a cell ...
        print(f"Analyzing {len(cells)} cells")
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

These convenience methods send commands to the viewer on its channel:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
viewer.clear_highlights()
viewer.set_color_by("leiden")  # obs field name
viewer.set_visibility([1, 2, 3], visible=False)  # hide
viewer.reset_view()
```

If you need low-level control:

```python
viewer.send_message({"type": "highlight", "cells": [1, 2, 3], "color": "#ff0000", "label": None})
```

```{note}
A message sent before the viewer is displayed waits and is delivered, in order, once the viewer connects; `send_message` returns `False` while it waits and `True` when it was sent at once. The `show_anndata(...)` helper displays automatically in notebook contexts.
```

To show your own columns and gene sets in the viewer, push them:
`viewer.push(adata.obs["leiden"])`, `viewer.push_gene_set("T cells", genes)`
({doc}`43_push_clusters_refine_and_pull`).

---

## Edge cases (things that confuse even experienced users)

### Multiple viewers in one notebook

Each viewer has its own viewer id and token and its own channel, so events
reach the right Python object. (Browser-tab viewers also have their own
server each: the first on port 8765, the next on the next free port.)

Best practices:
- keep a variable name per viewer (`viewer_a`, `viewer_b`)
- call `viewer.stop()` when you’re done

### Long-running analysis inside callbacks

If your callback takes seconds/minutes, the UI may feel laggy.

Options:
- do minimal work inside the callback and enqueue to a background thread/process
- read `viewer.selection` in the next cell and run analysis in the main flow instead

### Hover events can be very frequent

Treat hover as “streaming telemetry”:
- debounce/throttle your handler logic
- do not run heavy analysis on every hover

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “I select cells but Python receives nothing”

Likely causes:
- the kernel is busy: a running cell holds the viewer's messages until it
  ends
- the selection was not confirmed (a lasso that is drawn but not confirmed
  sends nothing)
- the viewer's channel is not open (the app did not load)

How to confirm (in a new cell):
```python
viewer.debug_connection()
```

Fix:
- let the cell finish; hooks run while the kernel is idle
- check `channel_connected` and `connected_outputs` in the report; an empty
  output in JupyterLab or Notebook 7 needs `anywidget` where the notebook
  server runs ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`)

### Symptom: “Python receives events but highlight doesn’t show”

Likely causes:
- your cell indices are wrong (reordered dataset)

Fix:
- name the cells by obs name (`viewer.highlight_cells(adata.obs_names[mask])`)
- print a few indices and sanity-check by clicking cells

### Symptom: “My notebook becomes slow over time”

Likely causes:
- multiple viewers left running
- accumulating callbacks

Fix:
- call `viewer.stop()` for viewers you no longer use
- clear hooks when needed: `viewer.clear_hooks()`

---

## Next steps

- {doc}`32_session_persistence_and_restoring_analysis_artifacts` (save/restore state; reproducibility)
- {doc}`05_jupyter_embedding_hooks_sessions_gallery` (runnable hook/session notebooks)
- {doc}`06_cellucid_in_your_pipeline_python_and_r` (pull, push, gene sets, exports; Python and R)
- Deep docs: {doc}`../e_jupyter_hooks/index`
