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
viewer = show_anndata(
    "data.h5ad",
    dataset_name="Example",
    dataset_id="example",
)
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

prepare(
    ...,
    out_dir="./my_export",
    dataset_name="Example",
    dataset_id="example",
    obs_categorical_dtype="uint16",
)
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
2) Before binding, Python establishes the exact source web generation; the
   browser then loads it from the same server.
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
   `/_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`
4) The server streams the upload to a temporary file and notifies Python.
5) Python returns a `CellucidSessionBundle(Path(...))`.

This is covered in depth in: {doc}`05_sessions_to_anndata_bridge`.

---

## Multi-dataset folders: one server, many datasets

If you point a server at a directory containing multiple dataset subfolders, the server can expose a list:

- `GET /_cellucid/datasets` → dataset list + relative paths
- Each dataset folder should contain a `dataset_identity.json` (written by `prepare(...)`)
- `server.viewer_url` opens this served catalog without choosing an arbitrary
  entry. Exactly one declared dataset may auto-open; multiple datasets require
  an exact dataset-id selection.

The root is intentionally strict. Once any immediate subdirectory looks like
an exported-dataset candidate, every immediate subdirectory must be a complete
current export. A stray directory rejects the root rather than disappearing
from the catalog.

This matters for:
- demos (“choose a dataset”),
- hosting multiple exports behind one URL,
- stable dataset IDs across collaborators.

See: {doc}`04_dataset_identity_and_reproducibility`.

---

## Edge cases (common confusion points)

### “Is my data being uploaded to cellucid.com?”

In notebook/server modes, the dataset is served from your Python process. Each
viewer-serving startup downloads and verifies the configured web generation
(see {doc}`06_privacy_security_and_offline_vs_online`), but your dataset is not
uploaded by default.

### “Why do I see `127.0.0.1:8765` instead of `https://cellucid.com`?”

To avoid cross-origin issues, the Python server serves the exact verified
viewer generation from the same origin as the dataset.

### “Why do events contain indices instead of cell IDs?”

Indices are fast, compact, and universal across export/AnnData modes. The downside is that they are fragile if you reorder/subset cells after the fact.
Treat index order as part of dataset identity (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Troubleshooting

### Symptom: “The notebook iframe is blank / mixed-content blocked”

Likely causes:
- Your notebook page is served from HTTPS (JupyterHub, Colab, remote), but the viewer tries to load `http://127.0.0.1:<port>`.

How to confirm:
- The browser console shows a mixed-content or connection error, and
  `viewer.viewer_url` is not reachable from the browser.

Fix options:
1) Configure an HTTPS route for one fixed Cellucid port, or use SSH port
   forwarding if the kernel is remote.
2) Pass that browser-reachable server base as `client_server_url=` when
   constructing the viewer.

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
- Understand privacy/network requirements:
  {doc}`06_privacy_security_and_offline_vs_online`
