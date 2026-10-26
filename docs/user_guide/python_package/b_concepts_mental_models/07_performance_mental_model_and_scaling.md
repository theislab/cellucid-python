# Performance mental model and scaling

**Audience:** everyone (computational users get the most value)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What costs happen when a file is prepared vs when it is viewed
- Why `prepare(...)` exists (and when it’s worth it)
- Where bottlenecks usually come from (CPU, disk, network, GPU)
- Practical scaling recommendations and edge cases

---

## Mental model (one sentence)

Cellucid performance is best when the file is laid out for the browser **once** (gene-major, chunked, natively decodable) and view-time work is mostly **range reads + GPU rendering**.

---

## One reader, two kinds of file

Every route — `show_anndata`, `serve_anndata`, `cellucid serve`, the web app's
file buttons and links — hands the viewer the data, which its own reader reads
(an AnnData or MuData in place with its HDF5/Zarr reader; a Loom or 10x file
in place too; an R object or a MatrixMarket folder once from the start). For
an AnnData, what differs is how the file is laid out.

### A Cellucid-ready AnnData (`prepare(...)`)

What gets “paid once”, when `prepare()` writes the file:
- converting the expression matrix to gene-major (CSC) float32,
- chunking every array for range reads and compressing it with a codec
  browsers decode natively.

What gets “paid every time you view”: the range reads a view needs, and GPU
rendering. On a 561,947-cell dataset over a simulated 80 ms / 50 MB/s link, the
first gene took 86 ms from the prepared `.h5ad` and 4,110 ms from the same
data as anndata writes it by default
({doc}`../c_data_preparation_api/10_performance_and_storage`).

This is the recommended path for:
- large datasets,
- sharing/hosting,
- repeated viewing.

### Any AnnData (`show_anndata(...)`)

The server computes nothing. It reads the store's metadata once (cell and gene
counts, the names of the stored elements — never an array), checks the launch
choices you named, and serves the store's bytes. The web app's own HDF5/Zarr
reader then reads the store in place with HTTP range requests.

What gets “paid at view time”, in the browser:
- range reads of exactly the arrays a view needs (an embedding, an obs column,
  a gene),
- for a row-major (CSR) expression matrix, a gene index the browser builds once:
  in the background as soon as the dataset opens when the server is on the same
  machine or network, on the first gene over a slow link,
- the neighbor graph, read only when a view draws it.

This is the recommended path for:
- quick exploration,
- iterating on preprocessing before writing a file to share.

Measured on one Mac over localhost, for a 562k-cell × 8,192-gene `.h5ad` (CSR,
gzip): server startup 0.16 s; obs fields about 90 ms; the 2D and 3D UMAPs
about 40 ms each; the one-time gene index builds in the background from the
moment the dataset opens (about 2 s), after which a gene takes a few ms and 100
genes about 20 ms (a gene asked for at once waits for the index, 1.3 s); the
neighbor graph 0.65 s when drawn. An in-memory object is written once to a temporary `.h5ad` with
column-major (CSC) sparse matrices — 0.8 s for the same object — after which
the first gene took 8 ms and 100 genes 68 ms.

Startup reports four numbered steps (`Detecting format`, `Reading metadata`,
`Checking launch choices`, `Starting server`).

For the step-by-step cost breakdown, see
{doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.

---

## The bottleneck map (what can be slow)

### 1) Disk I/O

- The browser's range reads of the `.h5ad` file or Zarr store, served from
  disk as it is.

Symptoms:
- slow first load,
- stutter when switching fields/gene expression.

### 2) CPU (serialization + decompression + preprocessing)

- decompression in the browser (deflate natively; Blosc/LZ4/Zstandard in JavaScript)
- the one-time gene index of a row-major (CSR) matrix
- centroid computations when a field is coloured, and outlier quantiles on the
  outlier slider's first move when the file does not store them

Symptoms:
- high CPU usage during load,
- UI responsiveness drops when scrubbing controls quickly.

### 3) Network transfer (even on localhost)

The viewer loads data via HTTP requests.

For large datasets, the number of requests matters as much as their size: a
single file (`.h5ad`, `.zarr.zip`) lets the reader merge neighbouring chunk
reads, while a `.zarr` folder costs one request per chunk.

### 4) GPU / WebGL

Rendering millions of points is GPU-bound.

Symptoms:
- low FPS,
- “WebGL context lost”,
- smoke/volumetric modes are especially heavy.

Web app performance guidance lives here:
- {doc}`../../web_app/n_benchmarking_performance/index`

---

## Practical recommendations by dataset size (rule of thumb)

```{note}
These are guidelines, not hard limits. Hardware, browser, and field complexity matter.
```

### Small (≤ ~50k cells)

- `show_anndata(...)` on the data you have is usually fine.
- `prepare()` is optional unless you want to share or host it.

### Medium (~50k–500k cells)

- `prepare()` once is recommended for repeated viewing, sharing and hosting.

### Large (≥ ~500k cells)

- Prefer a Cellucid-ready AnnData (`prepare()`) for reliable UX: its
  column-major (CSC) matrix spares the browser the one-time gene index a
  row-major matrix needs.
- Host one file (`.h5ad` or `.zarr.zip`) rather than a `.zarr` folder.
- Consider server mode + SSH tunneling for remote workflows.

Opening the AnnData you have is the right tool while preprocessing is still
moving. Once the object stops changing and is to be shared or hosted, prepare
it.

---

## Performance knobs you can control (Python-side)

### Preparing (`prepare(...)`)

`prepare()` has no compression or precision knobs: the layout is the one
measured fastest for the browser
({doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`).
What you control is what goes in:

- `gene_identifiers=[...]` keeps chosen genes,
- subset cells before the call (`adata[mask]`),
- pass only the embeddings and `obs` columns you need,
- choose the container: a `.zarr` folder for GitHub (a repository or its
  Pages site), an `.h5ad` (or `.zarr.zip`) for everything else, local files
  included ({ref}`which-format`).

### Opening any AnnData (`show_anndata(...)`)

The viewer reads the data where it is, and the costs follow the storage:
- a row-major (CSR) matrix costs one gene index, built in the browser on the
  first gene; a column-major (CSC) matrix reads each gene as one contiguous
  range,
- an in-memory AnnData is written once to a temporary `.h5ad` (sparse
  matrices CSC, uncompressed) and a MuData to a temporary `.h5mu`, removed when
  the server stops or Python exits; it needs disk space for one uncompressed
  copy,
- a store holding an array the browser cannot decode (an HDF5 filter such as
  Blosc or Zstandard, a Zarr codec such as bz2) is read once with anndata and
  served from a decodable temporary `.h5ad`; startup names those arrays.

---

## Edge cases (performance-specific)

### “It’s fast on my laptop but slow on the cluster”

Remote notebooks add complexity:
- proxy/tunnel overhead,
- weaker CPU per user,
- shared disks.

Mitigations:
- use a fixed port and SSH tunnel,
- keep the file on a local SSD of the compute node if possible,
- ensure the configured viewer source and selected generation directory have
  adequate throughput (see
  {doc}`06_privacy_security_and_offline_vs_online`).

### “My prepared file is huge”

Most often:
- you kept every gene (or did not subset cells),
- integer counts above 2^24 made the matrix float64,
- the expression matrix is dense.

Fix:
- keep fewer genes with `gene_identifiers=[...]` if appropriate,
- pass a sparse matrix: `prepare()` stores sparse input as CSC,
- drop rarely used `obs` columns before the call.

### “Gene expression is slow”

Likely causes:
- the first gene of a row-major (CSR) matrix: the browser builds the gene index
  once, and every later gene is fast,
- running on network-mounted storage, or over a high-latency tunnel.

Fix:
- write the file with `prepare()` (column-major) when first-gene latency matters,
- ensure the `.h5ad` is on fast storage, close to the server.

---

## Troubleshooting

### Symptom: “The first load takes forever”

Likely causes:
- the complete configured web generation is being fetched and verified for
  this startup,
- very large embeddings (millions of cells),
- slow disk.

How to confirm:
- Run `viewer.debug_connection()` and check `web_ui.cache` and `viewer_index_probe`.
- Watch network requests in the browser devtools network tab.

Fix:
- diagnose source download or generation-directory throughput,
- write a Cellucid-ready AnnData with `prepare()`,
- move data to faster storage.

### Symptom: “Startup sat there for minutes with nothing happening”

Likely causes:
- the viewer generation is being fetched (every viewer-serving startup
  downloads and verifies the complete generation again),
- an in-memory object, or a store with arrays the browser cannot decode, is
  being written to a temporary `.h5ad`.

How to confirm:
- read the last line printed. Each of the four steps names the phase it is in.

Fix:
- make the configured web source fast to reach from the Python runtime (each
  startup fetches the generation in full; a previous one is never served
  instead),
- serve the `.h5ad` or Zarr store from disk instead of an in-memory object,
- rewrite a store whose arrays use a codec the browser cannot decode with one
  it can (gzip, LZF for HDF5).

### Symptom: “Interactions lag when I scrub sliders”

Likely causes:
- expensive recomputation in the viewer (filters, outlier thresholds),
- CPU-bound rendering options.

Fix:
- use coarser steps (avoid continuous scrubbing),
- reduce simultaneous overlays,
- switch to a lighter render mode.

---

## Next steps

- Debugging checklist: {doc}`08_debugging_mental_model_where_to_look`
- Sizes and costs of a prepared file: {doc}`../c_data_preparation_api/10_performance_and_storage`
