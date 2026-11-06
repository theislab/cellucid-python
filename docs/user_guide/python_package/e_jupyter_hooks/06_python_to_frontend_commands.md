# Python → frontend actions

This page documents how your notebook Python code **controls and reads** the
embedded Cellucid viewer.

Everything a person can do in the viewer, and everything an agent reads from
it, is an **action** of one catalogue: colour by a field, highlight cells,
filter, rename a category, push a column, move the camera, read the history,
read a field's categories. The catalogue is defined once in the web app, and
the Python package is generated from it: the viewer has **one method per
action** (`viewer.color_by(...)`, `viewer.highlight_cells(...)`,
`viewer.history_get()`, …) and one generic call, `viewer.act(type, **params)`.
A notebook's action is applied exactly as a person's click is: one step of the
viewer's history, labelled as the notebook's, undone like any other.

Under the hood:
- the viewer opens one authenticated channel to Python (the viewer channel,
  {doc}`05_architecture_viewer_channel`): over the notebook widget's comm, or
  a WebSocket to a browser-tab viewer's server;
- Python sends each action as one `action` command: a JSON text frame, or one
  binary frame when its parameters hold arrays (cells, a column's values);
- the viewer answers each action once, with its outcome (and a read's result)
  or with the one sentence that says why it did not apply it;
- an action sent before the viewer has connected waits and is sent, in order,
  when it connects; nothing is sent twice.

If you’re new, start with the quickstart: {doc}`02_quickstart_minimal_roundtrip`.

## At a glance

**Audience**
- Wet lab / beginner: use the “Fast path recipes”.
- Computational users: read “Parameters” and “Knowing whether an action applied”.
- Developers: read “Message schema” + {doc}`05_architecture_viewer_channel`.

## Fast path recipes (copy/paste)

### Colour by an obs field (or a gene)

```python
viewer.color_by(field={"source": "obs", "key": "cell_type"})
viewer.color_by(field={"source": "var", "key": "CD3E"})
```

### Highlight cells

Cells are rows of `obs` (0-based), a boolean mask, or obs names:

```python
viewer.highlight_cells(cells=[0, 10, 42], color="#ff00ff")
viewer.highlight_cells(cells=adata.obs["leiden"] == "3", label="Leiden 3")
```

### Clear every highlight

```python
viewer.highlight_clear_all()
```

### Reset the camera of the live view

```python
viewer.camera_reset(view="live")
```

### Hide a set of cells, then show every cell again

```python
viewer.cells_visibility_set(cells=[0, 10, 42], visible=False)
viewer.cells_visibility_set(cells=None, visible=True)
```

### Show a column computed in Python

```python
viewer.push(adata.obs["leiden"], key="leiden")
```

### Read something back

A read changes nothing and answers with its result:

```python
request = viewer.fields_list()
# next cell (a notebook widget hears the viewer between cells)
request.result()
```

```{important}
An action sent before the viewer has connected waits in Python and is sent as
soon as the viewer connects; the viewer applies it once its dataset is open.
After `viewer.stop()` an action raises `RuntimeError`. To act only once the
viewer shows the data, register a hook:

```python
@viewer.on_ready
def _ready(event):
    viewer.color_by(field={"source": "obs", "key": "cell_type"})
```
```

## Every action is a method

The methods are generated from the catalogue the package ships
(`cellucid/action_catalog.json`, byte for byte the web app's), so the Python
API names every action the viewer has and nothing it lacks:

```python
from cellucid.jupyter import action_catalog

sorted(action_catalog())            # every action type: 'camera.reset', 'color.by', ...
action_catalog()["color.by"]["summary"]
```

A method's name is its action's type with dots and dashes as underscores
(`filter.category.visibility` → `viewer.filter_category_visibility(...)`);
its docstring is the action's title, summary and the description of each
parameter. `viewer.act("color.by", field=...)` is the same action named in
data. The full list is in {doc}`../g_api_reference_coverage/api/viewers`.

## Parameters

Parameters are keyword arguments and are checked against the action before
anything is sent: a parameter the action does not have, or a required one
that is missing, raises `TypeError` naming the action's parameters; an unknown
action raises `ValueError` with the nearest names. Values are checked by the
viewer, which answers with the reason when it refuses one.

- **Fields, genes and categories are named**, never positions: a field is
  `{"source": "obs", "key": "cell_type"}` (or `"var"` and a gene name as the
  viewer shows it), a category by its label, a highlight page by its name. A
  step therefore replays on another dataset with the same names.
- **Cells** (the `cells` of `highlight_cells`, `cells_visibility_set`,
  `cells_values`, …) are 0-based rows of `obs` -- original positions,
  whatever subsample the viewer opened --, given as a list of ints, a NumPy
  integer array (the fast form for millions of cells), a boolean mask of every
  cell, a boolean `pandas.Series` (its index's cells that are True), or obs
  names. Order is kept; a cell named twice, a negative row or one beyond the
  dataset is refused. They travel as `uint32` in a binary frame.
- **One number per cell** (`field_push(values=...)`) travels as its value
  type: `uint8`/`uint16`/`uint32` codes or `float32`/`float64` as they are,
  booleans as `uint8`, other integers as `float64`.
- Everything else is JSON: NumPy scalars and arrays, tuples and pandas
  indexes are turned into Python values and lists.
- **Optional parameters** default to the action's own default where it is a
  single value (`tool="annotation"`, `all_branches=False`), else to `OMIT`:
  left out, the viewer decides. `None` is sent as JSON `null` (for
  `mapping_set(knn=None)`, "no neighbour graph", unlike leaving `knn` out).
- A parameter the catalogue calls `in` (the analysis actions: which window)
  is written `in_`.

## Knowing whether an action applied

Every action returns an `ActionRequest` at once; nothing waits for the viewer.
The viewer answers each request once:

```python
request = viewer.color_by(field={"source": "obs", "key": "cell_type"})
# next cell
request.done          # True once the viewer answered
request.result()      # the outcome (None for most changes; a read's result)
request.step          # {"id": "s4", "label": "Coloured by cell_type"}: the history step
```

`request.result(timeout=60.0)` waits for the answer (in a browser-tab viewer,
or from a thread of your own). In a notebook widget the answer arrives once
the cell that sent the action has finished, so call it in a later cell; in the
same cell it raises at once and says so. When the viewer did not apply the
action -- a field that does not exist, a category that is not there, a value
out of range -- `result()` raises `cellucid.jupyter.ActionError` with the
viewer's reason (`error.reason`), for example:

```text
ActionError: The viewer did not apply color.by: There is no field "CD8A"; fields like it: CD8, CD8B.
```

Nothing has to collect a request: an action applies whether or not anyone
waits for it. A refusal nobody collects is still reported -- Cellucid prints
it:

```text
[cellucid] viewer.color_by(...) did not take effect in the viewer: There is no field "CD8A"; fields like it: CD8, CD8B.
[cellucid] The viewer is unchanged. Handle this in code with @viewer.on_command_error.
```

To act on refusals in code, register a hook, which also silences that notice:

```python
rejected = []

@viewer.on_command_error
def record(event):
    rejected.append(event)        # {"command": "action", "requestId": ..., "reason": ...}
```

Record rather than raise: hooks run on the viewer channel's thread, so raising
there only logs the error and drops the event instead of failing your cell.
`event["requestId"]` is the failed request's `request.request_id`. See
{doc}`07_frontend_to_python_events` for the exact records.

## Python conversions kept as their own methods

A few methods turn Python objects into an action's parameters:

| Python method | What it does | Action |
|---|---|---|
| `viewer.push(values, key=None)` | Shows a pandas column or an array as a field: labels, clusters and bools as a categorical, numbers as a continuous field, aligned by obs names | `field.push` |
| `viewer.request_session()` | Asks for the viewer's whole history as a session file ({doc}`10_session_bundles_request_session_bundle`) | (the `requestSessionBundle` command) |
| `viewer.apply_session_to_anndata(adata, request, ...)` | Plans applying a step of that session to an AnnData; `plan.apply()` writes it ({doc}`11_session_bundles_apply_session_to_anndata`) | (none: Python) |
| `viewer.load_session(path)` | Adds a saved session's trackers to the viewer's history ({doc}`10_session_bundles_request_session_bundle`) | (the binary `load_session` command) |
| `with viewer.turn(label):` | Makes every action sent inside the block one turn of the history | (`turn` and `turnLabel` on each `action`) |

Diagnostics / utilities:

| Python method | Purpose | Command |
|---|---|---|
| `viewer.act(type, **params)` | Any action, named in data | `action` |
| `viewer.send_message(message)` | Low-level escape hatch; returns `True` (sent) or `False` (queued) | `message["type"]` |
| `viewer.debug_connection()` | End-to-end connectivity report | `ping`, `debug_snapshot` |
| `viewer.stop()` | Freeze the view, stop answering it (a browser-tab viewer: stop its server), unregister hooks, and report any failure | `freeze` |

## The history

Each action a notebook applies is a step of the viewer's history, labelled as
the notebook's (Python or R), beside the person's own steps and an agent's.
`viewer.history_get()` reads it: the trackers, the branches, every step with
who took it and what it did, and the step shown now
(`viewer.history_get(all_branches=True)` for every branch). Undo, redo and
jumps are actions too, as the viewer's own buttons do them:
`viewer.history_undo()` and `viewer.history_redo()` answer `{"moved", "head"}`
(`moved` false when there was no step to move to), `viewer.history_jump(step="s3")`
shows any step `history_get()` lists and answers `{"head"}`. They move the
history and record no step, so `request.step` is `None`. Every other part of
the history is an action of the same family (`viewer.history_branch_switch(...)`,
`viewer.history_note_add(...)`, `viewer.history_bookmark_set(...)`, …), as the
State panel and agents use them.

A **turn** groups a notebook's actions: every action sent inside
`with viewer.turn(label):` carries the turn's id (`python-` and 12
hexadecimal digits, which the block yields) and its label, and the State
panel shows those steps as one group, undone together, as it does an agent's
turn. Turns do not nest: opening one inside another raises `RuntimeError`.

```python
with viewer.turn("Mark the T cells") as turn_id:
    viewer.color_by(field={"source": "obs", "key": "cell_type"})
    viewer.highlight_cells(cells=adata.obs["cell_type"] == "T cell", label="T cells")
```

`viewer.history` follows the history as it changes: `last_change` (the
latest `history_changed` event), `tracker`, `head` and `branch` from it,
`on_change(callback)`, and `read()`, which asks for the history as a session
(`viewer.request_session()`). The whole history travels with the session file
({doc}`10_session_bundles_request_session_bundle`), a saved one is added to
the viewer with `viewer.load_session(path)`, and the state at any of its
steps can be written onto an AnnData
({doc}`11_session_bundles_apply_session_to_anndata`). How the history looks and
works in the viewer: {doc}`../../web_app/m_state_history/index`.

## Message schema (what actually crosses the boundary)

### Commands carry no credentials

The viewer channel is authenticated once, when the viewer opens it, and
belongs to that one viewer, so a command is only its own fields below.
`viewer.send_message({...})` refuses a message that names `viewerId` or
`viewerToken`.

### `action`

```json
{ "type": "action", "requestId": "...", "action": "color.by", "params": { "field": { "source": "obs", "key": "cell_type" } },
  "turn": null, "turnLabel": null }
```

`action` is the action's dotted type; `params` its parameters; `turn` and
`turnLabel` the notebook turn it belongs to (both text inside
`viewer.turn(...)`, both null outside one). The viewer answers with one `action_result` (binary: its outcome as UTF-8 JSON, and the
history step `{id, label}` or `null`) or one `command_error` carrying the same
`requestId`.

### An `action` with arrays is one binary frame

When the parameters hold arrays -- cells, a column's values -- the action
travels as one binary frame: a little-endian `uint32` header length, the
UTF-8 JSON header (padded with spaces to a multiple of 4 bytes, so the body
starts on a 4-byte boundary and the viewer reads the arrays in place), then
each array's elements, little-endian, one after another:

```json
{ "type": "action", "requestId": "...", "action": "highlight.cells",
  "params": { "color": "#e41a1c" }, "turn": null, "turnLabel": null,
  "arrays": [ { "path": ["cells"], "type": "uint32", "length": 3 } ] }
```

`path` names where the array goes in `params` (keys, and list indices for a
list of objects, ending in a key `params` does not hold); `type` is `uint8`,
`uint16`, `uint32`, `float32` or `float64`. A frame that breaks these rules
closes the channel (1008), as a malformed event does in the other direction.
Two million cells are 8 MB written in one pass.

### `requestSessionBundle` (used by `viewer.request_session()`)

```json
{ "type": "requestSessionBundle", "requestId": "..." }
```

Answered by a binary `session_bundle` event (the `.cellucid-session` file of
format 4) or by a `command_error` with the same `requestId`.

### `load_session` (used by `viewer.load_session()`)

One binary frame: the header below, padded the same way, then the
`.cellucid-session` file. `name` is the file's name, without a folder.

```json
{ "type": "load_session", "requestId": "...", "name": "my-session.cellucid-session" }
```

Answered by a `session_loaded` event (`{tracker: {id, name}, head}`) or by a
`command_error` with the same `requestId`.

### `freeze` (used by `viewer.stop()`)

```json
{ "type": "freeze" }
```

### `ping` / `debug_snapshot` (internal, used by `viewer.debug_connection()`)

```json
{ "type": "ping", "requestId": "..." }
```

```json
{ "type": "debug_snapshot", "requestId": "..." }
```

## Performance notes

- **Huge cell lists** travel as binary frames, 4 bytes a cell: pass a NumPy
  array (`np.flatnonzero(mask)`) or the mask itself rather than a list for
  millions of cells, so the cells are checked and written in vectorized
  passes.
- **High-frequency actions** (a highlight every mousemove) are not
  recommended: each is a step of the history. Prefer debouncing and batching.
- If you want to drive many interactions, use hooks to *decide* what to do,
  and send actions only when needed.

## Edge cases and footguns

- **Cells are rows**: they refer to row positions of `obs`, not stable cell
  IDs; pass obs names when the order of your object may differ.
- **Fields and categories are names**: a misspelt name is refused by the
  viewer with what it has ("its categories are …").
- **Viewer not displayed yet**: actions wait and are sent when it connects.
- **Viewer stopped**: an action, and `viewer.send_message(...)`, raise
  `RuntimeError`.
- **Another dataset shown**: while the person shows a dataset the notebook
  did not open in the viewer, the notebook's actions are refused.
- **Multiple viewers**: make sure you are sending actions to the same
  `viewer` instance whose output you are looking at.

## Troubleshooting

If “nothing happens” when sending actions:

1. Confirm the viewer is displayed and ready (in a new cell):
   ```python
   viewer.state.ready      # None until the dataset has opened
   ```
2. Collect the request: `request.result()` raises with the viewer's reason
   when it refused the action.
3. Run:
   ```python
   viewer.debug_connection()
   ```
   - `channel_connected` False: the viewer has not opened its channel (for a
     browser-tab viewer the report also runs a ping/pong).
4. Check whether the viewer refused actions rather than ignoring them:
   ```python
   viewer.debug_connection()["recent_events"]["command_errors"]
   ```
5. If the channel is connected and nothing was refused, but an action is
   ignored: re-run the cell to create a fresh viewer (for a browser-tab
   viewer, `viewer.clear_web_cache()` first re-establishes the app build).
6. Use the full guide: {doc}`14_troubleshooting_hooks`

## Next steps

- Viewer → Python events: {doc}`07_frontend_to_python_events`
- Robust callback patterns: {doc}`08_writing_robust_callbacks`
- Full reference: {doc}`16_reference`
