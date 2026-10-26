# Large dataset: server mode + lazy gene expression

This tutorial is the “make it work when it’s big and remote” guide.

It explains:
- what “server mode” means in Cellucid (and why it exists)
- how “lazy gene expression” avoids shipping a giant matrix to the browser
- which workflow to pick for large datasets
- what direct-AnnData startup does, as four numbered steps you can read
- what the browser reads, and when, and what that costs at scale
- how to make notebook embeds work on HTTPS/remote kernels (JupyterHub/HPC)
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
view needs cross the wire.

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

When you pass a path, it is served as it is on disk — never copied — and the
browser reads it in place: an `.h5ad`, a Zarr store folder (v2 or v3), or a
zipped Zarr store (`.zip`) — or a MuData, Loom, 10x or MatrixMarket file or
folder, or an R object, which the browser reads too
({doc}`25_serve_loom_10x_and_matrixmarket`).

An in-memory `AnnData` object is written once to a temporary `.h5ad`, with
sparse matrices (`X`, layers, `raw.X`) column-major (CSC) and uncompressed, and
that file is removed when the server stops or Python exits. On the same
562k-cell object, writing it took 0.8 s (925 MB on disk); after that the first
gene took 8 ms and 100 genes 68 ms, because a CSC matrix needs no gene index.

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
is read once with anndata and its codec plugins and served from a decodable
temporary `.h5ad`; the terminal names those arrays.

### Which `obsm` array gets drawn

A key that names its dimension decides it. Failing that, a plain `X_umap` of 1,
2, or 3 columns is read at the dimension its column count states, so a Scanpy
object is drawn as written. Full rules:
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.

## Remote notebooks (HTTPS / JupyterHub / HPC)

### The problem

If your notebook runs at an HTTPS origin (common on JupyterHub), the browser may not be able to load an `http://127.0.0.1:<port>` iframe directly:
- HTTPS → HTTP is often blocked (mixed content)
- remote kernels mean the browser’s “localhost” is not the kernel machine

### Configure one browser-reachable URL

Expose the Cellucid port through your notebook deployment or reverse proxy,
then pass that exact base URL, and name the proxy's host so the server answers
it (and accepts the viewer channel's WebSocket from its pages):

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    client_server_url="https://notebooks.example.org/user/alice/proxy/8765",
    allowed_hosts=["notebooks.example.org"],
)
```

A notebook viewer serves on port 8765, or on the next free port when another
viewer holds it (a printed note says which); the proxy route names that port.

`client_server_url` must be an absolute HTTP(S) base URL without a trailing
slash, query, fragment, or credentials.

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
`serve_anndata(..., host="0.0.0.0")`; a notebook viewer always binds
`127.0.0.1`). The address you open stays `http://127.0.0.1:<port>/…`. See
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
- remote notebook proxy adding latency

Fix:
- move data to fast local SSD (if possible)
- write it once with `prepare()` (gene-major) for repeated viewing

### Symptom: “Startup sits for a long time before printing a URL”

Likely causes:
- the viewer generation is being fetched and verified (`[4/4] Starting server`):
  every start downloads the whole declared web build again and checks every
  file, so a slow link to the configured source slows every start
- an in-memory object, or a store with arrays the browser cannot decode, is
  being written to a temporary `.h5ad`

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
- no `obsm` array can be plotted: each one is a table, sparse, not one row per
  cell, or non-finite in the columns that would be drawn

How to confirm:
- read the message in the viewer. It lists every `obsm` entry and why it was
  not used.

Fix:
- repair or compute an embedding, or name a usable one with
  `embedding_keys={"2d": "<key>"}` (`--embedding-2d <key>` from the terminal),
  then re-save the file if you serve a path rather than an in-memory object.

### Symptom: remote viewer iframe is blank or unreachable

Likely causes:
- notebook served over HTTPS and the kernel server is only reachable over HTTP loopback

Fix:
- expose the server through the notebook deployment, and
- pass that exact base URL as `client_server_url` and its host as
  `allowed_hosts=[...]`; the proxy must pass WebSockets
  (`/_cellucid/ws`, the viewer channel).

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
