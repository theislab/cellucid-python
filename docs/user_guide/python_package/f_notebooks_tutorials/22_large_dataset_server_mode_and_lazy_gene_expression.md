# Large dataset: server mode + lazy gene expression

This tutorial is the “make it work when it’s big and remote” guide.

It explains:
- what “server mode” means in Cellucid (and why it exists)
- how “lazy gene expression” avoids shipping a giant matrix to the browser
- which workflow to pick for large datasets
- what direct-AnnData startup does, as four numbered steps you can read
- what the browser reads, and when, and what that costs at scale
- why a notebook viewer needs nothing on a remote kernel, and how a browser
  tab reaches one (JupyterHub/HPC)
- what to do when you hit performance limits (browser/GPU/memory)

---

## At a glance

**Audience**
- Computational users (primary)
- Power users (remote/HPC notebooks)

**Time**
- ~30–60 minutes (longer if you set up SSH tunneling for the first time)

**Prerequisites**
- `pip install cellucid`
- A dataset:
  - `.h5ad`, or
  - a Zarr store (a folder, or zipped as `.zip`), or
  - a Cellucid-ready AnnData written as in {doc}`21_prepare_a_cellucid_ready_anndata`
  - (MuData, Loom, 10x and MatrixMarket files open the same way:
    {doc}`24_citeseq_mudata_prepare_and_show`,
    {doc}`25_serve_loom_10x_and_matrixmarket`)

---

## Key concepts

### “Server mode”

Cellucid’s viewer runs in the browser, but your dataset may be:
- too large to load into browser memory
- too large to ship as a single download
- on a remote machine the browser cannot open files from

So Cellucid uses an HTTP server that serves:
- the viewer UI assets
- the data: the bytes of one dataset (an `.h5ad` or Zarr store, or any other
  format the web app reads) exactly as they are on disk

The server computes nothing. The web app's own HDF5/Zarr
reader reads the store in place with HTTP range requests, so only the bytes a
view needs cross the wire. A notebook viewer (`show_anndata`, a widget) runs
no server: the same reads travel through the notebook's own connection to the
kernel, and are answered the same way.

### “Lazy gene expression”

The browser does not download the entire expression matrix upfront.

Instead, when you search/select a gene, the viewer reads only that gene's
entries in the store's expression matrix, by range request.

This means:
- time-to-first-view can be fast even for huge `n_genes`
- gene coloring is “pay as you go” (first time you request a gene, you pay I/O)

For a row-major (CSR) matrix — the layout Scanpy writes — a gene's entries are
spread across every row, so the browser builds a gene index once, on the first
gene, and every later gene reads only its own entries. A column-major (CSC)
matrix — the layout `prepare()` writes — needs no index: one gene is one
contiguous range.

### “Read when drawn”

The neighbor graph in `adata.obsp`, and any vector field in `obsm`, is read
only when a view draws it. Nothing about them costs startup time.

---

## Choose your workflow (decision guide)

### Option A (recommended for sharing/repeated viewing): `prepare()` once, then open the file

Use this if:
- you plan to view the dataset many times
- you want a stable artifact to share (one standard AnnData file)

Workflow:
- `prepare("my_dataset.h5ad", format="h5ad", ...)` → a Cellucid-ready AnnData
- `show_anndata("my_dataset.h5ad")` (notebook) or `cellucid serve my_dataset.h5ad` (browser); the file names itself

See:
- {doc}`21_prepare_a_cellucid_ready_anndata`

### Option B (a large `.h5ad` or Zarr store as it is): open it directly

Use this if:
- you want convenience and are iterating quickly
- you have a large `.h5ad` or Zarr store and don’t want to convert it

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=650,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer
```

`dataset_name` and `dataset_id` may be omitted when the file declares them in
`uns['cellucid']['name']` and `uns['cellucid']['id']`.

### Option C (browser-only, no notebook): `cellucid serve /path/to/data`

Use this if:
- you don’t need notebook embedding
- you want a simple “serve and open browser” workflow

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
```

---

## Reading the startup report

Startup prints four numbered steps: `[1/4] Detecting format`,
`[2/4] Reading metadata` (cell and gene counts and the names of the stored
elements — no array is read), `[3/4] Checking launch choices` (each option you
named must exist in the store), and `[4/4] Starting server`.

The full transcript is in
{doc}`../../web_app/b_data_loading/04_server_tutorial`.

## What it costs at scale

Measured on one Mac over localhost, for a 562k-cell × 8,192-gene `.h5ad` (CSR
matrix, gzip-compressed):

| What | Time |
| --- | --- |
| Server startup | 0.16 s |
| obs fields | about 90 ms |
| 2D and 3D UMAP | about 40 ms each |
| Gene index, built once in the background from the moment the dataset opens | about 2 s |
| A gene, once the index is built | a few ms |
| A gene asked for before the index is built | 1.3 s |
| 100 genes at once | about 20 ms |
| 1,000 cell rows | 0.55 s |
| Neighbor graph, when drawn | 0.65 s |

Choosing a graph at launch costs nothing at startup either:

```python
viewer = cellucid.show_anndata(
    "data.h5ad",
    dataset_name="My dataset",
    dataset_id="my-dataset",
    knn_key="connectivities",
)
```

`knn_key` names the `obsp` matrix the viewer draws; it must exist in the store,
and it is read only when a view draws the graph.

Details: {doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.

## Paths and in-memory objects

When you pass a path the browser can read as it is, it is served as it is on
disk — never copied — and the browser reads it in place: an `.h5ad`, a Zarr store folder (v2 or v3), or a
zipped Zarr store (`.zip`) — or a MuData, Loom, 10x or MatrixMarket file or
folder, or an R object, which the browser reads too
({doc}`25_serve_loom_10x_and_matrixmarket`). A file whose cells are not in
priority order (one `prepare()` did not write) is served instead as a store in
priority order made from the file as the viewer reads it, its compressed
arrays decoded once into temporary memory-mapped files removed when the server
stops ({doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`).

An in-memory `AnnData` object is served as a Zarr store with its cells in
priority order, each chunk made when the viewer asks for it from a snapshot of
the object in memory-mapped files removed when the server stops or Python
exits: its expression matrices laid out gene-major in the background after the
call returns (so a gene is one range; the arrays read stay read-only until
then, and later edits to the object do not show). A
subsample the viewer opens is then the first rows of every array and the first
entries of every gene, so only its cells are sent
({doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading` has the
measurements).

```python
import anndata
from cellucid import show_anndata

adata = anndata.read_h5ad("data.h5ad")
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

A store holding an array the browser cannot decode (an HDF5 filter such as
Blosc, Zstandard, szip or scaleoffset; a Zarr codec such as bz2, lzma or zfpy)
is decoded on the kernel, chunk by chunk as the viewer asks, through the same
store; the terminal names those arrays.

### Which `obsm` array gets drawn

A key that names its dimension decides it. Failing that, a plain `X_umap` of 1,
2, or 3 columns is read at the dimension its column count states, so a Scanpy
object is drawn as written. Full rules:
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.

## Remote notebooks (HTTPS / JupyterHub / HPC)

A notebook viewer needs nothing for a remote kernel: `show_anndata(...)` is a
notebook widget, which reaches the kernel through the notebook's own
connection, wherever the kernel runs and whether the notebook is HTTP or
HTTPS. The browser loads the app from `https://www.cellucid.com`; the
dataset's bytes go from the kernel to the viewer through the widget
({doc}`../e_jupyter_hooks/03_supported_environments_matrix`).

A viewer in a browser tab of its own is the other choice for a remote
machine: serve it there and forward its port.

### SSH port forwarding for a browser tab

On the remote machine:
```bash
cellucid serve my_dataset.h5ad --dataset-name "My study" --dataset-id my-study-v1 --host 127.0.0.1 --port 8765 --no-browser
```

On your laptop:
```bash
ssh -L 8765:127.0.0.1:8765 user@remote-host
```

In your browser:
- open the printed Viewer URL:
  `http://127.0.0.1:8765/?anndata=true`

```{warning}
Be careful with `--host 0.0.0.0`. It exposes the server to your network. Use it only if you understand the security implications.
```

### Serving from a compute node

A server on an HPC compute node usually cannot be tunnelled to directly; the
forward terminates on a login node, so the server must listen on every
interface: `cellucid serve ... --host 0.0.0.0` (or
`serve_anndata(..., host="0.0.0.0")`; an `AnnDataViewer` always binds
`127.0.0.1`, and a notebook widget binds nothing). The address you open stays `http://127.0.0.1:<port>/…`. See
{doc}`../d_viewing_apis/17_hpc_slurm_and_compute_node_serving`.

## Performance and scaling checklist

### 1) Browser/GPU limits are real

Very large point clouds can hit:
- WebGL limits (context loss)
- GPU memory pressure

If the viewer becomes unstable, consider:
- downsampling for exploratory viewing
- filtering to a subset of cells first
- using the web app performance guidance:
  - {doc}`../../web_app/n_benchmarking_performance/index`

### 2) Prefer a Cellucid-ready AnnData for speed and stability

`prepare()` writes:
- the expression matrix gene-major, so every gene is one small read
- arrays chunked for range reads, with codecs browsers decode natively
- one stable file that names itself

See:
- {doc}`21_prepare_a_cellucid_ready_anndata`

### 3) Expect “first gene request” latency

Lazy gene expression means:
- the first time you request a gene, it may take noticeable time (disk/network);
  for a row-major (CSR) matrix the browser builds a gene index once: in the
  background as soon as the dataset opens when the link to the server is as
  fast as a disk (the same machine or network), else with the first gene
- subsequent genes read only their own entries

If gene requests are consistently slow:
- your storage is slow (network filesystem)
- the server runs far from the data, or the tunnel has high latency
- write the file with `prepare()` (column-major), so each gene is one contiguous range

### 4) Serve only the obs columns you need

`obs_keys=` (or `--obs-key`, repeated) fixes both the shown set and its order.
Each name must exist in `obs`; the check runs before the server starts.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Viewer loads but gene search is extremely slow”

Likely causes:
- the `.h5ad` sits on slow storage
- a row-major (CSR) matrix: the first gene waits for the gene index
- a slow link between the browser and a remote kernel (a notebook widget
  reads through the notebook's connection)

Fix:
- move data to fast local SSD (if possible)
- write it once with `prepare()` (gene-major) for repeated viewing

### Symptom: “Startup sits for a long time before printing a URL”

Likely causes:
- the viewer generation is being fetched and verified (`[4/4] Starting server`):
  every start downloads the whole declared web build again and checks every
  file, so a slow link to the configured source slows every start
- a gene asked for before its matrix is laid out (an in-memory object's
  expression matrices, and a cell-major `X` of a file, are laid out gene-major
  in the background; the viewer asks again when they are)

How to confirm:
- read the last line printed; each of the four steps names its phase.

Fix:
- check the link to the configured source (`https://www.cellucid.com`, or
  `web_source_url=`); a server started with `serve_anndata(..., serve_web_ui=False)`
  (`cellucid serve --no-web-ui`) serves the data without fetching the viewer build
- serve the file from disk instead of an in-memory object

### Symptom: “The edge overlay is missing, but my object has a KNN graph”

Likely causes:
- the viewer opened with a different (or no) neighbor graph.

How to confirm:
- open **Session › Data mapping › Neighbour graph** and check which `obsp`
  matrix is chosen.

Fix:
- choose the graph there, or restart with `knn_key="connectivities"`
  (`--knn connectivities` from the terminal).

### Symptom: “Cellucid found nothing to plot”

Likely causes:
- no `obsm` array can be plotted: each dense one has another row count than
  the cells, places no cell (every row non-finite), or holds a value past
  float32 in the columns that would be drawn (an `obsm` with no dense numeric
  array at all is refused as having no cell embedding:
  {ref}`no-cell-embedding`)

How to confirm:
- read the message in the viewer. It lists every `obsm` entry and why it was
  not used.

Fix:
- repair or compute an embedding, or name a usable one with
  `embedding_keys={"2d": "<key>"}` (`--embedding-2d <key>` from the terminal),
  then re-save the file if you serve a path rather than an in-memory object.

### Symptom: a browser-tab viewer on a remote machine is blank or unreachable

Likely causes:
- the browser cannot reach the remote server (no tunnel), or the page is
  reached through a proxy that was not declared

Fix:
- forward the port with SSH (above), or
- behind a reverse proxy, pass its exact base URL as `client_server_url` and
  its host as `allowed_hosts=[...]` to `AnnDataViewer` (`--allowed-host` for
  `cellucid serve`); the proxy must pass WebSockets (`/_cellucid/ws`, the
  viewer channel).

A notebook viewer (`show_anndata`) needs none of this: it reaches the kernel
through the notebook's own connection, remote or not.

### Symptom: “WebGL context lost / blank canvas”

Likely causes:
- GPU memory pressure (too many points, too high render settings)
- browser/GPU driver issues

Fix:
- reduce dataset size (downsample/filter)
- close other GPU-heavy tabs
- confirm hardware acceleration and WebGL2 are enabled, then update the
  browser, operating system, and GPU driver when applicable

---

## Next steps

- {doc}`23_programmatic_highlighting_and_selection_callbacks` (hooks + robust callback patterns)
- {doc}`04_real_world_dataset_recipes_gallery` (runnable real Pancreas preparation)
- Web app performance guidance: {doc}`../../web_app/n_benchmarking_performance/index`
