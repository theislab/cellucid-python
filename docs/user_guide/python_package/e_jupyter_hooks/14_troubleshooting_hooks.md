# Troubleshooting (hooks)

This page is a symptom-based guide for when the notebook viewer and its hooks
don’t behave as expected. It covers the notebook widget (`show_anndata()`)
first, then what is particular to a browser-tab viewer
({class}`~cellucid.AnnDataViewer`, `serve_anndata`, `cellucid serve`).

Rule of thumb: in a notebook, **most surprises come from timing** — the widget
hears the viewer only between cells — not from hook code.

## First step (always): collect a debug report

Run, in a new cell:

```python
report = viewer.debug_connection()
report
```

A notebook widget's report (it waits for nothing, so it works in any cell):

- `transport` (`"notebook widget"`) and `web_source_url` (where the browser
  loads the app from)
- `displayed` and `connected_outputs` (the outputs whose page reached the kernel)
- `channel_connected` (whether the viewer opened its channel)
- `dataset` (the id and format the widget serves)
- `state` (`ready`, the last event type and when it arrived)
- `recent_events` (counts by event type, and the commands the viewer refused)

A browser-tab viewer's report adds its server: `server_health`,
`anndata_store_probes` (a `HEAD` and a one-byte `Range` read of the served
store), a ping/pong round trip (`frontend_roundtrip`) and the page's
`frontend_debug_snapshot`.

If you include this report when asking for help, debugging is usually 10× faster.

---

## Symptom: the output is empty, or says the notebook cannot show widgets

### Likely causes (ordered)
1. `anywidget` is not installed where the **notebook server** runs (JupyterLab,
   Notebook 7). Installing `cellucid` installs it in the kernel's environment,
   which is a different one in some setups.
2. The front end does not render widgets at all.
3. In marimo, the widget is not the cell's output.

### Fix
- `pip install anywidget` in the notebook server's environment and reload the
  page.
- Use a front end that shows widgets ({doc}`03_supported_environments_matrix`).
- In marimo, end the cell with `viewer`.

---

## Symptom: the output says the app did not load

The widget's message reads “The page at https://www.cellucid.com has not
answered as a Cellucid viewer…”.

### Likely causes
1. The browser cannot reach `web_source_url` (offline laptop, blocked
   address, strict content policy).
2. `web_source_url` names a copy of the app that is too old to be embedded.

### Fix
- Allow `https://www.cellucid.com`, or host a current copy of the app and pass
  `show_anndata(..., web_source_url="https://cellucid.example.org")`.

---

## Symptom: a call raises “… a notebook widget hears the viewer only between cells”

`wait_for_ready()`, `wait_for_event()`, a session request's `result()` and
`viewer.apply_session_to_anndata(adata, request)` wait for the viewer's
answer, and in a notebook nothing the viewer says can arrive while the
cell that waits is running. They raise at once and name the alternative:

- readiness: `viewer.state.ready` in a later cell, or `@viewer.on_ready`;
- the next selection: `viewer.selection` in a later cell, or `@viewer.on_selection`;
- results (highlights, labels, scores, gene sets): `viewer.pull(adata)` in a
  later cell ({doc}`../g_api_reference_coverage/api/round_trip`);
- a session bundle: `request = viewer.request_session_bundle()`, then
  `request.result()` in the next cell.

---

## Symptom: viewer loads but **no events arrive** (`on_selection` never fires)

### Likely causes (ordered)
1. A cell is still running: the widget's events are handled while the kernel
   is idle, and they all arrive once it is.
2. You are interacting with a stale output (kernel restarted); it keeps its
   last view and says the notebook stopped it.
3. A hook raised: the exception is in the log (`cellucid.viewer_channel`).
   The event was still recorded (`viewer.state`), and the other hooks ran.
4. The selection was not confirmed: moving a lasso without confirming it
   sends no selection event.
5. The kernel's machine could not hold a very large selection. A selection
   costs 4 bytes a cell on the wire and about 36 bytes a cell as the
   `list[int]` hooks receive.

### How to confirm
- `viewer.debug_connection()`: `channel_connected` must be `True`, and
  `recent_events["by_type"]` counts what arrived.
- Register a raw logger:
  ```python
  @viewer.on_message
  def _log(event):
      print(event)
  ```

### Fix
- Let the cell finish (or interrupt a long one).
- Stale output: re-run the cell that creates the viewer.
- Very large selections: build highlight groups in smaller steps, and bring
  them back with `viewer.pull(adata)`.

---

## Symptom: events arrive, but indices don’t match my AnnData

### Likely causes
1. You reordered/subset `adata` after creating the viewer.
2. You are using a different `adata` object than the one you passed to `show_anndata`.
3. The viewer opened a subsample of the cells: event cells are still rows of
   the full `obs` (original positions), but only the opened cells can appear.

### How to confirm
```python
selection = viewer.selection
print(selection.cells.min(), selection.cells.max(), len(selection))
print("adata.n_obs:", adata.n_obs)
```

If `max >= adata.n_obs`, you are definitely mismatched.

### Fix
- Only interpret indices against the exact dataset the viewer is showing.
- If you must reorder/subset, create a new viewer from the modified AnnData.
- Prefer obs names: `viewer.selection.obs_names`.

---

## Symptom: Python commands don’t affect the UI (highlight/color/reset do nothing)

### Likely causes (ordered)
1. The viewer **rejected** the command — the field, index or colour was not one
   it could apply.
2. The viewer has not connected yet: commands wait in Python and are delivered,
   in order, once it does (the viewer loads after the cell that showed it ends).
3. The output you are looking at belongs to another (or a stopped) viewer.

### How to confirm
- A rejected command prints to `sys.stderr` as soon as the viewer reports it:
  ```text
  [cellucid] viewer.set_color_by(...) did not take effect in the viewer: ...
  ```
  If you missed it, list the rejections after the fact:
  ```python
  viewer.debug_connection()["recent_events"]["command_errors"]
  ```
  A non-empty list means the command arrived and was refused — fix the
  argument, not the connection.
- `viewer.debug_connection()["channel_connected"]` says whether a viewer page
  is connected now.

### Fix
1. If the command was rejected, correct the argument the reason names — most
   often a field ID the loaded dataset does not have.
2. If the output is stale, stop that viewer (`viewer.stop()`), then create and
   show a current one.

### Prevention
- Send commands from a hook that runs when the viewer is ready:
  ```python
  @viewer.on_ready
  def _ready(event):
      viewer.set_color_by("cell_type")
  ```
- In a script, collect rejections and assert on them at a point you control,
  rather than reading a printed line:
  ```python
  rejected = []

  @viewer.on_command_error
  def record(event):
      rejected.append(event)

  ...
  assert not rejected, rejected
  ```

---

## Symptom: session bundle capture fails

### Likely causes
1. In a notebook, it was asked for and collected in the same cell; use
   `request = viewer.request_session_bundle()` and `request.result()` in the
   next cell.
2. The viewer's channel is not connected.
3. The temporary directory has less free space than the bundle;
   `result()` raises “The session bundle could not be saved: …”.

### Fix
- Ask again in one cell and collect in the next.
- To bring results onto AnnData, `viewer.pull(adata)` needs no bundle.

---

## Browser-tab viewers only

### Symptom: the viewer UI cannot be established or served

`AnnDataViewer(...)` / `serve_anndata(...)` raises while establishing the
viewer's web build, `viewer.display()` raises while verifying it, or a page
that was already open reloads to the plain text “Internal server error”.

**Likely causes (ordered)**
1. No outbound network access to the configured web source (every server
   start establishes the build again).
2. The source inventory or one declared asset failed exact validation.
3. The selected generation directory is not writable.
4. The generation was removed after the viewer started (for example
   `viewer.clear_web_cache()`).

**How to confirm**
- In the debug report: `report["web_ui"]["cache"]` and
  `report["viewer_index_probe_error"]`.
- Open `<viewer.server_url>/index.html` in a browser.

**Fix**
- Ensure the Python runtime can reach the configured `web_source_url`.
- Correct the exact inventory, object, MIME, length, hash, or build error shown.
- Pass a writable `web_cache_dir` when the default temporary location is not
  suitable.
- Before a shared session, exercise the source with
  `viewer.ensure_web_ui_cached(force=True)`.

### Symptom: the browser cannot reach the server (remote kernel, HTTPS page)

**Likely causes**
1. The kernel runs on another machine and nothing forwards its port.
2. A reverse proxy forwards the page, but its host is not declared
   (`421 Misdirected Request`) or it does not pass WebSockets (the viewer says
   it cannot open its channel).

**Fix (recommended order)**
1. In a notebook, use `show_anndata()` instead: it needs no address.
2. Use SSH port forwarding with a fixed `port=`
   ({doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`).
3. Behind a reverse proxy, pass its exact browser-reachable base as
   `client_server_url=` and its host as `allowed_hosts=[...]`, and make it pass
   WebSockets ({doc}`13_security_cors_origins_and_mixed_content`).

### Symptom: “Port already in use” / the viewer port changed

- Without `port=`, an `AnnDataViewer` takes 8765, or the next free port when
  another viewer (or program) holds it, and prints a note saying so (Google
  Drive files open only on 8765). A port you pass is used or refused, never
  replaced.
- Stop viewers when done (`viewer.stop()`, or
  `cellucid.jupyter.cleanup_all()`), and choose a fixed port for tunnels:
  ```python
  from cellucid import AnnDataViewer

  viewer = AnnDataViewer(
      "data.h5ad",
      port=8765,
      dataset_name="My study",
      dataset_id="my-study-v1",
  )
  ```

### Debugging checklist (browser tab)

```python
print("server_url:", viewer.server_url)
print("viewer_url:", viewer.viewer_url)
viewer.wait_for_ready(timeout=60)
viewer.debug_connection()
```

In a browser:
- open `<server_url>/_cellucid/health`
- open devtools → Network, filter `WS`, and select the `/_cellucid/ws`
  connection to see its frames

---

## Next steps

- Reference (schemas + endpoints): {doc}`16_reference`
