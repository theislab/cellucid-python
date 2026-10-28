# Data flows

**Audience:** everyone  
**Time:** 15–20 minutes  
**Goal:** understand where data lives (disk/server/browser) and how the Python package connects to the Cellucid web app.

---

## Fast path: pick a workflow (you can change later)

Every path from data to the viewer reads the data (AnnData, MuData, Loom,
10x, MatrixMarket or an R object) with the web app's own reader. The main
decision is whether to prepare the file first:

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

When you use `cellucid-python`, you almost always have two actors and one
connection between them:

1) **Python** (your script / notebook kernel), which answers the viewer's
   reads of the dataset's bytes
2) **The Cellucid web app**, running in your browser: inside a notebook
   widget's output (`show_anndata`), or in a browser tab of its own
3) **The connection**: in a notebook, the notebook's own connection (a
   Jupyter comm, no port); for a browser tab, a local HTTP server started by
   Python (usually on `127.0.0.1:<port>`)

The key idea: the viewer is a web app, so it reads the data itself, with
range requests, and Python only hands it bytes.

```text
Notebook widget (show_anndata):

Python kernel ◀── the notebook's own connection (comm) ──▶ widget output
   │  reads (Range) answered                               │ frames the app from
   │  one viewer channel: commands ⇄ events                ▼ https://www.cellucid.com
                                                        Cellucid web app

Browser tab (AnnDataViewer, serve_anndata, cellucid serve):

Python ──starts──▶ Local Cellucid HTTP server ◀── fetch() (Range) ── Cellucid web app
                          ▲                                            (browser tab)
                          └──── one WebSocket, /_cellucid/ws ──────────┘
```

---

## Data flow by mode (what changes)

### Mode 1: Notebook + `show_anndata(...)`

**What you have:** an AnnData or MuData object, or any file or folder Cellucid
reads (`.h5ad`, a Zarr store, `.h5mu`, `.loom`, a 10x `.h5` or `outs/`
folder, a MatrixMarket folder, an `.rds`/`.rda`).

**What happens:**
1) Python prepares the dataset: for AnnData and MuData it reads only the
   metadata and checks the launch choices you named against it; it answers
   reads with the bytes as they are on disk, or -- for an in-memory AnnData,
   and an AnnData file not in priority order -- with the chunks of a store in
   priority order made as they are asked for (a MuData object is written once
   to a temporary `.h5mu`). No server is started.
2) The widget shows the app (loaded by the browser from
   `https://www.cellucid.com`); the app's own reader reads the data itself
   through the widget: in place with range requests where the format allows
   it (only the bytes a view needs are transferred), once from the start for
   an R object or a MatrixMarket folder. Python computes nothing per request.
3) Commands and events travel on the viewer channel, over the same
   connection, handled by the kernel between cells.

**Why this is convenient:** nothing is converted; any AnnData opens.

**What it costs:** the browser does the reading, so the first gene of a
row-major (CSR) matrix waits for a gene index the browser builds once; later
genes are fast. A file written by `prepare()` is gene-major and needs no index.

### Mode 2: CLI / standalone server (`cellucid serve ...`)

**What you have:** a path to any file or folder Cellucid reads (see Mode 1).

**What happens:**
1) You run a long-lived server process in a terminal.
2) Anyone who can reach the server URL can load the viewer and dataset (depending on host/binding).

This mode is ideal for:
- remote/HPC workflows (with SSH tunnels),
- demos for teammates on the same network,
- separating “viewer runtime” from “analysis notebook”.

### Mode 3: Web app file picker (no Python server)

This is primarily a **web app** workflow: you open the Cellucid web app and load a file or folder it reads (an `.h5ad`, a `.zarr.zip` or `.zarr` folder, an `.h5mu`, a `.loom`, a 10x `.h5` or `outs/` folder, a MatrixMarket folder, an `.rds`/`.rda`) with its **Open…** button (or by dropping it on the page).

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
2) The web app sends the selection on its viewer channel (opened with its
   `viewerId` and `viewerToken`: over the widget in a notebook, a WebSocket to
   `/_cellucid/ws` in a browser tab) as one
   binary frame: the selected cell indices as little-endian `uint32`,
   ascending (other events, such as hover and click, are JSON text frames).
3) Python hands it to the viewer object the channel belongs to.
4) The viewer object triggers your hook callback.

```{important}
The payload uses **cell indices** (row positions). If the dataset changes row order, indices refer to different cells.
This is a common source of “I selected one thing but analyzed another”.
```

---

## How results come back (frontend → Python)

Every result the viewer holds — highlight groups, labels, clusters, scores,
gene sets — is sent to Python as it is made, so one call in a later cell
writes them all onto your AnnData:

```python
viewer.pull(adata)
```

A whole session is a *pull* workflow in two cells (a notebook widget hears
the viewer only between cells):

```python
request = viewer.request_session_bundle()
```

```python
bundle = request.result()
```

What happens:
1) Python sends a `requestSessionBundle` command on the viewer channel.
2) The web app serializes the current session state into bytes (`.cellucid-session`).
3) The web app sends those bytes back on the channel as one binary frame.
4) Python streams the frame to a temporary file.
5) `request.result()` returns a `CellucidSessionBundle(Path(...))`.

This is covered in depth in: {doc}`05_sessions_to_anndata_bridge`.

---

## Edge cases (common confusion points)

### “Is my data being uploaded to cellucid.com?”

No. In notebook and server modes the dataset is read from your Python
process by the app running in your own browser. In a notebook the browser
loads the app's code from `https://www.cellucid.com`; a browser-tab server
downloads and verifies the app's code at startup (see
{doc}`06_privacy_security_and_offline_vs_online`). Your dataset is not
uploaded.

### “Why does a browser tab show `127.0.0.1:8765` instead of `https://www.cellucid.com`?”

To avoid cross-origin issues, a browser-tab server serves the exact verified
viewer generation from the same origin as the dataset. A notebook widget shows
the app from `https://www.cellucid.com` and reaches the data through the widget.

### “Why do events contain indices instead of cell IDs?”

Indices are fast, compact, and the same in every mode. The downside is that they are fragile if you reorder/subset cells after the fact.
Treat index order as part of dataset identity (see {doc}`04_dataset_identity_and_reproducibility`).

---

## Troubleshooting

### Symptom: “The notebook output is empty”

Likely causes:
- The notebook front end does not show widgets (`anywidget` missing where the
  notebook server runs).
- The browser cannot load the app from `https://www.cellucid.com` (the output
  says so).

Fix options:
1) Install `anywidget` where the notebook server runs, and reload the page.
2) Allow the address, or host a copy of the app and pass `web_source_url=`.

### Symptom: “A browser-tab viewer cannot be reached”

The kernel runs on another machine: forward its port over SSH, or pass a
reverse proxy's address as `client_server_url=` and its host as
`allowed_hosts=` ({doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`).

### Symptom: “Hooks never fire”

Likely causes:
- the viewer is not fully loaded, or a cell is still running (a widget's
  events arrive between cells),
- you registered hooks on a different viewer instance than the one you see,
- in a browser tab, the channel (`/_cellucid/ws`) could not open (a proxy
  that does not pass WebSockets; the viewer says so).

How to confirm:
- Run `viewer.debug_connection()` in a new cell (see {doc}`08_debugging_mental_model_where_to_look`).

Fix:
- let the cell finish; check `viewer.state.ready` in a later cell,
- re-run the cell that creates the viewer, and register hooks on that object.

---

## Next steps

- Understand persistence: {doc}`03_state_persistence_and_scope`
- Understand identity/reproducibility: {doc}`04_dataset_identity_and_reproducibility`
- Understand privacy/network requirements:
  {doc}`06_privacy_security_and_offline_vs_online`
