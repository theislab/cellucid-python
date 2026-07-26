# Large dataset: server mode + lazy gene expression

This tutorial is the “make it work when it’s big and remote” guide.

It explains:
- what “server mode” means in Cellucid (and why it exists)
- how “lazy gene expression” avoids shipping a giant matrix to the browser
- which workflow to pick for large datasets
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

---

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
- open `http://127.0.0.1:8765/`

```{warning}
Be careful with `--host 0.0.0.0`. It exposes the server to your network. Use it only if you understand the security implications.
```

---

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
- {doc}`04_real_world_dataset_recipes_gallery` (dataset-specific export notebooks)
- Web app performance guidance: {doc}`../../web_app/n_benchmarking_performance/index`
