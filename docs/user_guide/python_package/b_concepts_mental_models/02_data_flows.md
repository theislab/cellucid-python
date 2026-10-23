# Data flows

**Audience:** everyone  
**Time:** 15–20 minutes  
**Goal:** understand where data lives (disk/server/browser) and how the Python package connects to the Cellucid web app.

---

## Fast path: pick a workflow (you can change later)

Every path from data to the viewer reads **AnnData** with the web app's own
reader. The main decision is whether to prepare the file first:

- **Maximum convenience** right now: open the AnnData you have (`show_anndata(...)`).
- **Maximum performance + shareability**: write a Cellucid-ready AnnData once
  (`prepare(...)`), then open that file anywhere.

### Workflow A — Open the AnnData you have

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
- avoiding a preparation step while iterating.

Tradeoffs:
- the first gene of a row-major (CSR) matrix waits for a gene index the
  browser builds once,
- the file carries no name, citation or launch defaults of its own.

### Workflow B — Prepare first (recommended for large, shareable viewing)

```python
from cellucid import prepare, show_anndata

path = prepare(
    "example.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="example",
    dataset_name="Example",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
)
viewer = show_anndata(path)
```

Good for:
- fast loads (a gene-major matrix: every gene is one small read),
- a file that names itself and opens with your chosen defaults,
- publishing/sharing a standard AnnData as an artifact.

Tradeoffs:
- you spend time once to write the file,
- you manage ids across versions (covered in {doc}`04_dataset_identity_and_reproducibility`).

---

## The “where does it run?” map

When you use `cellucid-python`, you almost always have three actors:

1) **Python** (your script / notebook kernel)  
2) **A local HTTP server** (started by Python, usually on `127.0.0.1:<port>`)  
3) **The Cellucid web app** (running in your browser or notebook iframe)

The key idea: the viewer is a web app, so it loads data by making HTTP requests to a server.

```text
               (the AnnData store's bytes)
Python ──starts──────────────────────────▶ Local Cellucid HTTP server
  │                                             │
  │ postMessage                                 │ fetch()
  ▼                                             ▼
Notebook iframe / browser tab  ◀────────────  Cellucid web app UI
  ▲
  │  HTTP POST /_cellucid/events     (JSON events)
  │  HTTP POST /_cellucid/selection  (selections, binary cell indices)
  └─────────────────────────────────────────────── frontend → Python events
```

---

## Data flow by mode (what changes)

### Mode 1: Notebook + `show_anndata(...)`

**What you have:** an AnnData object or an AnnData store on disk (`.h5ad`, a
Zarr store folder, or a zipped Zarr store).

**What happens:**
1) Python starts an `AnnDataServer`. It reads only the store's metadata, checks
   the launch choices you named against it, and serves the store's bytes as
   they are on disk (an in-memory object is first written once to a temporary
   `.h5ad`).
2) The web app's own HDF5/Zarr reader reads the store in place with HTTP range
   requests: only the bytes a view needs are transferred. The server computes
   nothing per request.
3) Hooks/events work the same way (`/_cellucid/events`; selections via `/_cellucid/selection`).

**Why this is convenient:** nothing is converted; any AnnData opens.

**What it costs:** the browser does the reading, so the first gene of a
row-major (CSR) matrix waits for a gene index the browser builds once; later
genes are fast. A file written by `prepare()` is gene-major and needs no index.

### Mode 2: CLI / standalone server (`cellucid serve ...`)

**What you have:** an AnnData path (`.h5ad`, Zarr folder, or `.zip`).

**What happens:**
1) You run a long-lived server process in a terminal.
2) Anyone who can reach the server URL can load the viewer and dataset (depending on host/binding).

This mode is ideal for:
- remote/HPC workflows (with SSH tunnels),
- demos for teammates on the same network,
- separating “viewer runtime” from “analysis notebook”.

### Mode 3: Web app file picker (no Python server)

This is primarily a **web app** workflow: you open the Cellucid web app and load an `.h5ad`, a `.zarr.zip` or a `.zarr` folder with its **Files** buttons.

Key mental model:
- files are read locally by the browser,
- nothing is “uploaded” unless you explicitly host/share the file yourself.

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
2) The web app POSTs the selection to the server:
   - URL: `http://127.0.0.1:<port>/_cellucid/selection?viewerId=...&viewerToken=...&source=...`
   - Body: the selected cell indices as little-endian `uint32`, ascending
     (other events, such as hover and click, are JSON posted to `/_cellucid/events`).
   - The `viewerId` in the query lets the server route it.
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

## Edge cases (common confusion points)

### “Is my data being uploaded to cellucid.com?”

In notebook/server modes, the dataset is served from your Python process. Each
viewer-serving startup downloads and verifies the configured web generation
(see {doc}`06_privacy_security_and_offline_vs_online`), but your dataset is not
uploaded by default.

### “Why do I see `127.0.0.1:8765` instead of `https://www.cellucid.com`?”

To avoid cross-origin issues, the Python server serves the exact verified
viewer generation from the same origin as the dataset.

### “Why do events contain indices instead of cell IDs?”

Indices are fast, compact, and the same in every mode. The downside is that they are fragile if you reorder/subset cells after the fact.
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
- the browser cannot reach the server endpoint `/_cellucid/events` (or
  `/_cellucid/selection`, for selection hooks),
- you registered hooks on a different viewer instance than the one in the iframe.

How to confirm:
- Run `viewer.debug_connection()` (see {doc}`08_debugging_mental_model_where_to_look`).
- In the browser network tab, look for POST requests to `/_cellucid/events`
  (and to `/_cellucid/selection` when you select cells).

Fix:
- call `viewer.wait_for_ready(timeout=60)` before relying on hooks,
- ensure the viewer URL you opened matches the `viewer.viewer_url` printed by Python.

---

## Next steps

- Understand persistence: {doc}`03_state_persistence_and_scope`
- Understand identity/reproducibility: {doc}`04_dataset_identity_and_reproducibility`
- Understand privacy/network requirements:
  {doc}`06_privacy_security_and_offline_vs_online`
