# Data flows

**Audience:** everyone  
**Time:** 15–20 minutes  
**Goal:** understand where data lives (disk/server/browser) and how the Python package connects to the Cellucid web app.

---

## Fast path: pick a workflow (you can change later)

Cellucid supports multiple “paths” from data → interactive viewer. The main decision is:

- Do you want **maximum convenience** right now? Use **direct AnnData** (`show_anndata(...)`).
- Do you want **maximum performance + shareability**? Use **export-first** (`prepare(...)` then `show(...)` / `serve(...)`).

### Workflow A — Direct AnnData (no export)

```python
from cellucid import show_anndata
viewer = show_anndata("data.h5ad")  # or a zarr path, or an in-memory AnnData
```

Good for:
- exploration in notebooks,
- avoiding an “export step” while iterating.

Tradeoffs:
- slower than a pre-export for very large datasets,
- more moving parts at runtime (server must compute/stream on demand).

### Workflow B — Export-first (recommended for large, reproducible, shareable viewing)

```python
from cellucid import prepare, show

prepare(..., out_dir="./my_export")
viewer = show("./my_export")
```

Good for:
- fast loads,
- consistent results across machines,
- publishing/sharing an export folder as an artifact.

Tradeoffs:
- you spend time once to export,
- you need to manage versions/IDs (covered in {doc}`04_dataset_identity_and_reproducibility`).

---

## The “where does it run?” map

When you use `cellucid-python`, you almost always have three actors:

1) **Python** (your script / notebook kernel)  
2) **A local HTTP server** (started by Python, usually on `127.0.0.1:<port>`)  
3) **The Cellucid web app** (running in your browser or notebook iframe)

The key idea: the viewer is a web app, so it loads data by making HTTP requests to a server.

```text
               (static files) or (dynamic endpoints)
Python ──starts──────────────────────────▶ Local Cellucid HTTP server
  │                                             │
  │ postMessage                                 │ fetch()
  ▼                                             ▼
Notebook iframe / browser tab  ◀────────────  Cellucid web app UI
  ▲
  │  HTTP POST /_cellucid/events
  └─────────────────────────────────────────────── frontend → Python events
```

---

## Data flow by mode (what changes)

### Mode 1: Notebook + `show(...)` (pre-exported folder)

**What you have:** an export folder on disk.

**What happens:**
1) Python starts a `CellucidServer` that serves the export folder over HTTP.
2) The browser loads the Cellucid web app UI from the same server (hosted-asset proxy cache).
3) The web app fetches `points_2d.bin(.gz)` / `obs_manifest.json` / expression binaries, etc.
4) Interactions (selection/hover/click) can be sent back to Python via `/_cellucid/events`.

**Why this is fast:** the export folder is already in a viewer-optimized format (binary, quantized, compressible).

### Mode 2: Notebook + `show_anndata(...)` (direct AnnData)

**What you have:** an AnnData object or an AnnData-backed file (`.h5ad`, `.zarr`).

**What happens:**
1) Python starts an `AnnDataServer` with an `AnnDataAdapter`.
2) The web app asks the server for the same *logical* resources as in export mode (points, obs fields, gene expression),
   but the server produces them on demand from AnnData.
3) Hooks/events work the same way (`/_cellucid/events`).

**Why this is convenient:** you don’t have to export first.

**Why it can be slower:** dynamic conversion + network transfer happen at view time, and gene queries require server work.

### Mode 3: CLI / standalone server (`cellucid serve ...`)

**What you have:** an export folder or AnnData path.

**What happens:**
1) You run a long-lived server process in a terminal.
2) Anyone who can reach the server URL can load the viewer and dataset (depending on host/binding).

This mode is ideal for:
- remote/HPC workflows (with SSH tunnels),
- demos for teammates on the same network,
- separating “viewer runtime” from “analysis notebook”.

### Mode 4: Web app file picker (no Python server)

This is primarily a **web app** workflow: you open the Cellucid web app and load an export folder using the browser’s file picker.

Key mental model:
- files are read locally by the browser,
- nothing is “uploaded” unless you explicitly host/share the folder yourself.

Start here: {doc}`../../web_app/b_data_loading/index`

---

## How hooks/events actually travel (frontend → Python)

If you register a hook like:

```python
@viewer.on_selection
def handle(event):
    print(event["cells"][:10])
```

Here is the concrete path:

1) In the web app, a user selects cells (lasso/click/etc.).
2) The web app POSTs JSON to the server:
   - URL: `http://127.0.0.1:<port>/_cellucid/events`
   - Body includes `viewerId` so the server can route it.
3) The Python server calls the correct viewer object’s internal handler.
4) The viewer object triggers your hook callback.

```{important}
The payload uses **cell indices** (row positions). If the dataset changes row order, indices refer to different cells.
This is a common source of “I selected one thing but analyzed another”.
```

---

## How “no-download sessions” work (frontend → Python → file)

In notebooks, `viewer.get_session_bundle()` is a *pull* workflow:

```python
bundle = viewer.get_session_bundle(timeout=60)
```

What happens:
1) Python sends a `requestSessionBundle` command to the iframe (postMessage).
2) The web app serializes the current session state into bytes (`.cellucid-session`).
3) The web app uploads those bytes to:
   `/_cellucid/session_bundle?viewerId=...&requestId=...`
4) The server streams the upload to a temporary file and notifies Python.
5) Python returns a `CellucidSessionBundle(Path(...))`.

This is covered in depth in: {doc}`05_sessions_to_anndata_bridge`.

---

## Multi-dataset folders: one server, many datasets

If you point a server at a directory containing multiple dataset subfolders, the server can expose a list:

- `GET /_cellucid/datasets` → dataset list + relative paths
- Each dataset folder should contain a `dataset_identity.json` (written by `prepare(...)`)

This matters for:
- demos (“choose a dataset”),
- hosting multiple exports behind one URL,
- stable dataset IDs across collaborators.

See: {doc}`04_dataset_identity_and_reproducibility`.

---

## Screenshot placeholders (optional but helpful)

### Screenshot 1: server banner (terminal)

<!-- SCREENSHOT PLACEHOLDER
ID: python-data-flow-server-banner
Suggested filename: data_loading/python-server-banner.png
Where it appears: Python Package → Concepts & Mental Models → Data flows → CLI / standalone server
Capture:
  - UI location: terminal window
  - State prerequisites: a server started via `cellucid serve ...` or `python -c 'from cellucid import serve; serve(...)'`
  - Action to reach state: run `cellucid serve ./my_export --port 8765 --no-browser`
Crop:
  - Include: the printed Local URL + Viewer URL lines (these anchor the mental model)
  - Exclude: usernames, machine names, private paths
Redact:
  - Remove: dataset names/paths if sensitive
Annotations:
  - Callouts: #1 Local URL (server), #2 Viewer URL (web app entry point)
Alt text:
  - Terminal output showing Cellucid server running with local and viewer URLs.
Caption:
  - The Python server hosts both the dataset and the viewer UI entry point; you open the viewer URL in a browser.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a Cellucid server banner in a terminal.
:width: 100%

The server banner prints both the dataset server URL and the viewer URL you open in a browser.
```

### Screenshot 2: web app loading from a local server

<!-- SCREENSHOT PLACEHOLDER
ID: python-data-flow-viewer-loaded-from-local-server
Suggested filename: data_loading/viewer-loaded-from-localhost.png
Where it appears: Python Package → Concepts & Mental Models → Data flows → Notebook/CLI modes
Capture:
  - UI location: browser (or notebook iframe)
  - State prerequisites: viewer loaded successfully from `http://127.0.0.1:<port>/`
  - Action to reach state: open `viewer.viewer_url` or the CLI-printed viewer URL
Crop:
  - Include: the viewer canvas + any dataset name/point count indicator if present
  - Exclude: browser bookmarks, personal profile icons, private dataset IDs
Redact:
  - Remove: any sensitive dataset identifiers
Annotations:
  - Callouts: #1 viewer URL bar (localhost), #2 dataset loaded indicator
Alt text:
  - Cellucid viewer open in a browser and connected to a local Python server.
Caption:
  - The Cellucid web app UI runs in the browser while the dataset is served by the local Python server.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for Cellucid loaded from a local Python server.
:width: 100%

The viewer UI runs in the browser while the dataset is streamed from the local Python server.
```

---

## Edge cases (common confusion points)

### “Is my data being uploaded to cellucid.com?”

In notebook/server modes, the dataset is served from your Python process. The viewer UI assets may be downloaded once and cached (see {doc}`06_privacy_security_and_offline_vs_online`), but your dataset is not uploaded by default.

### “Why do I see `127.0.0.1:8765` instead of `https://cellucid.com`?”

To avoid mixed-content and cross-origin issues, the Python server serves the viewer UI from the same origin as the dataset (via a hosted-asset proxy + cache).

### “Why do events contain indices instead of cell IDs?”

Indices are fast, compact, and universal across export/AnnData modes. The downside is that they are fragile if you reorder/subset cells after the fact.
Treat index order as part of dataset identity (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Troubleshooting

### Symptom: “The notebook viewer says proxy required / mixed-content blocked”

Likely causes:
- Your notebook page is served from HTTPS (JupyterHub, Colab, remote), but the viewer tries to load `http://127.0.0.1:<port>`.

How to confirm:
- The embedded iframe shows a message about installing `jupyter-server-proxy`, or the browser console shows “mixed content” errors.

Fix options:
1) Install/enable `jupyter-server-proxy` in that environment (recommended).
2) Use SSH port forwarding if your kernel is remote.
3) Set `CELLUCID_CLIENT_SERVER_URL` to an HTTPS-reachable server URL if you have one.

### Symptom: “Hooks never fire”

Likely causes:
- the viewer is not fully loaded,
- the browser cannot reach the server endpoint `/_cellucid/events`,
- you registered hooks on a different viewer instance than the one in the iframe.

How to confirm:
- Run `viewer.debug_connection()` (see {doc}`08_debugging_mental_model_where_to_look`).
- In the browser network tab, look for POST requests to `/_cellucid/events`.

Fix:
- call `viewer.wait_for_ready(timeout=60)` before relying on hooks,
- ensure the viewer URL you opened matches the `viewer.viewer_url` printed by Python.

---

## Next steps

- Understand persistence: {doc}`03_state_persistence_and_scope`
- Understand identity/reproducibility: {doc}`04_dataset_identity_and_reproducibility`
- Understand privacy/offline: {doc}`06_privacy_security_and_offline_vs_online`
