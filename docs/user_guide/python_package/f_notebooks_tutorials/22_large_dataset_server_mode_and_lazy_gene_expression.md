# Large dataset: server mode + lazy gene expression

This tutorial is the “make it work when it’s big and remote” guide.

It explains:
- what “server mode” means in Cellucid (and why it exists)
- how “lazy gene expression” avoids shipping a giant matrix to the browser
- which workflow to pick for large datasets
- what direct-AnnData startup does, as five numbered steps you can read
- why the neighbor graph is never read unless you ask for it, and what that
  saves at scale
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
  - `.h5ad` (recommended for read-only-backed access), or
  - `.zarr` (loaded eagerly), or
  - a pre-exported folder from {doc}`21_prepare_exports_with_quantization_and_compression`

---

## Key concepts

### “Server mode”

Cellucid’s viewer runs in the browser, but your dataset may be:
- too large to load into browser memory
- too large to ship as a single file
- stored in formats that benefit from on-demand access

So Cellucid uses an HTTP server that serves:
- the viewer UI assets
- metadata manifests (obs/var/connectivity)
- per-gene and per-field binary vectors **on demand**

### “Lazy gene expression”

The browser does not download the entire expression matrix upfront.

Instead, when you search/select a gene, the viewer requests something like:
- `var/<gene_id>.values.*.bin` (exported mode)
- or a dynamic endpoint that returns that gene’s vector (AnnData server mode)

This means:
- time-to-first-view can be fast even for huge `n_genes`
- gene coloring is “pay as you go” (first time you request a gene, you pay I/O)

### “Opt-in connectivity”

The neighbor graph in `adata.obsp['connectivities']` is not read unless you ask
for it. Everything else — points, obs fields, gene expression, centroids, vector
fields — is unaffected by that choice.

---

## Choose your workflow (decision guide)

### Option A (recommended for sharing/repeated viewing): export once, then `show()` / `cellucid serve`

Use this if:
- you plan to view the dataset many times
- you want a stable artifact to share (export folder)

Workflow:
- `prepare(...)` → export folder
- `show("./exports/<id>")` (notebook) or `cellucid serve ./exports/<id>` (browser)

See:
- {doc}`21_prepare_exports_with_quantization_and_compression`

### Option B (recommended for large `.h5ad` without exporting yet): direct read-only-backed AnnData viewing

Use this if:
- you want convenience and are iterating quickly
- you have a large `.h5ad` and don’t want to export yet

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

### Option C (browser-only, no notebook): `cellucid serve /path/to/data`

Use this if:
- you don’t need notebook embedding
- you want a simple “serve and open browser” workflow

```bash
# Exported directory:
cellucid serve /path/to/export

# Direct AnnData:
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
```

---

## Reading the startup report

Startup prints five numbered steps. Step 2, `Loading AnnData`, is where a large
object spends nearly all of its time, and it names each phase of the build as it
starts, so a long pause is labelled rather than silent. A prepared export runs
three steps instead.

The full transcript is in
{doc}`../../web_app/b_data_loading/04_server_tutorial`.

## The graph and the overlay are opt-in

This is the setting that matters most at scale. The graph is **not** read unless
you ask for it:

```python
viewer = cellucid.show_anndata(
    "data.h5ad",
    dataset_name="My dataset",
    dataset_id="my-dataset",
    serve_connectivity=True,
)
```

Left off, `adata.obsp['connectivities']` is never touched and the dataset
publishes no connectivity. Asked for, the whole graph is read and validated
before the server binds, because the manifest the viewer fetches first declares
the edge count.

Why it dominates: the cost tracks *stored neighbors*, not cells. A 50-neighbor
graph over millions of cells is hundreds of millions of them, and validating one
needs several times its own memory. Ask for it when you will draw it.

Details and the arithmetic: {doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.

## Large `.h5ad` specifics (backed mode)

When you pass a `.h5ad` path, Cellucid always opens it read-only-backed:
- it does **not** load the entire matrix into RAM
- it reads data as needed (I/O-bound but memory-friendly)

To use an in-memory dataset, load an `AnnData` object yourself and pass that
object to `show_anndata` with the required identity. A `.zarr` path is loaded
eagerly by `anndata.read_zarr`.

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

### Which `obsm` array gets drawn

A key that names its dimension decides it. Failing that, a plain `X_umap` of 1,
2, or 3 columns is read at the dimension its column count states, so a Scanpy
object serves as written and a backed file needs no rewrite. Full rules:
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.

## Remote notebooks (HTTPS / JupyterHub / HPC)

### The problem

If your notebook runs at an HTTPS origin (common on JupyterHub), the browser may not be able to load an `http://127.0.0.1:<port>` iframe directly:
- HTTPS → HTTP is often blocked (mixed content)
- remote kernels mean the browser’s “localhost” is not the kernel machine

### Configure one browser-reachable URL

Expose the Cellucid port through your notebook deployment or reverse proxy,
then pass that exact base URL:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    client_server_url="https://notebooks.example.org/user/alice/proxy/8765",
)
```

`client_server_url` must be an absolute HTTP(S) base URL without a trailing
slash, query, fragment, or credentials.

### SSH port forwarding for a browser tab

On the remote machine:
```bash
cellucid serve exports/my_dataset --host 127.0.0.1 --port 8765 --no-browser
```

On your laptop:
```bash
ssh -L 8765:127.0.0.1:8765 user@remote-host
```

In your browser:
- open the printed Viewer URL:
  `http://127.0.0.1:8765/?source=remote`

```{warning}
Be careful with `--host 0.0.0.0`. It exposes the server to your network. Use it only if you understand the security implications.
```

### Serving from a compute node

A kernel or a server on an HPC compute node usually cannot be tunnelled to
directly; the forward terminates on a login node, which needs
`host="0.0.0.0"`. The address you open stays `http://127.0.0.1:<port>/…`. See
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

### 2) Prefer exports for speed and stability

Exports allow you to:
- quantize and compress values
- cache “what the browser needs” in a stable folder

See:
- {doc}`21_prepare_exports_with_quantization_and_compression`

### 3) Expect “first gene request” latency

Lazy gene expression means:
- the first time you request a gene, it may take noticeable time (disk/network)
- subsequent requests may be faster (server/browser caching)

If gene requests are consistently slow:
- your storage is slow (network filesystem)
- your `.h5ad` is not chunked/optimized for the access pattern
- use a read-only-backed `.h5ad` on fast local storage

### 4) Only ask for the neighbor graph when you will draw it

`serve_connectivity=True` / `--connectivity` is the most expensive switch on the
direct path, and it is paid in full at startup on every session. Leave it off
for exploration, and turn it on for the session that actually inspects edges. If
you need the graph repeatedly, `prepare(connectivities=...)` validates it once
into an export instead of once per launch.

### 5) Serve only the obs columns you need

`obs_keys=` (or `--obs-key`, repeated) fixes both the served set and its order.
Every column is classified during step 2, so naming a short list shortens the
longest step and drops columns Cellucid cannot serve, such as a `datetime64`
collection date.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Viewer loads but gene search is extremely slow”

Likely causes:
- backed `.h5ad` on slow storage
- large `.h5ad` with access pattern that forces many small reads
- remote notebook proxy adding latency

Fix:
- move data to fast local SSD (if possible)
- consider exporting (quantize+compress) for repeated viewing
- move the `.h5ad` to fast local storage

### Symptom: “Startup sits for minutes before printing a URL”

Likely causes:
- you are inside `[2/5] Loading AnnData` on a very large object
- you asked for the neighbor graph, and the edge list is being built

How to confirm:
- read the last line printed. `Obs columns: classifying …` is column
  classification; `Embeddings: resolving obsm keys` is embedding resolution;
  `Connectivity: reading obsp['connectivities']` is the edge list, and for
  50,000,000 or more stored neighbors the `cellucid.connectivity` logger has
  already warned with both counts.

Fix:
- drop `serve_connectivity=True` / `--connectivity` if this session will not
  draw edges
- narrow `obs_keys=` to the columns you actually color by
- export once with `prepare()` for repeated viewing

### Symptom: “The edge overlay is missing, but my object has a KNN graph”

Likely cause:
- the neighbor graph is opt-in, and this run did not ask for it.

How to confirm:
- read the `Connectivity:` line of `[3/5] Analyzing dataset`. It has three
  states:

  ```text
        Connectivity: yes (125,481,922 edges)
        Connectivity: not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)
        Connectivity: no
  ```

  The middle line means you own a graph and did not ask for it. `no` means the
  object really has none.

Fix:
- restart with `serve_connectivity=True`, or `--connectivity` from the terminal.

### Symptom: “No UMAP embedding Cellucid can read was found in adata.obsm”

Likely causes:
- `obsm` carries no `X_umap_1d`, `X_umap_2d`, `X_umap_3d`, or bare `X_umap`
- it carries a bare `X_umap` that is not 1, 2, or 3 columns wide

How to confirm:
- read the message. It lists the `obsm` keys the object does have, names the
  shape of `X_umap` when its width is the reason nothing resolved, and prints an
  exact assignment to run for every other array that could be an embedding.

Fix:
- assign the columns you mean to draw under the key for their own count, then
  re-save the file if you serve a path rather than an in-memory object:

  ```python
  adata.obsm["X_umap_2d"] = umap_coordinates_2d
  ```

### Symptom: remote viewer iframe is blank or unreachable

Likely causes:
- notebook served over HTTPS and the kernel server is only reachable over HTTP loopback

Fix:
- expose the server through the notebook deployment, and
- pass that exact base URL as `client_server_url`.

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
- {doc}`04_real_world_dataset_recipes_gallery` (runnable real Pancreas export)
- Web app performance guidance: {doc}`../../web_app/n_benchmarking_performance/index`
