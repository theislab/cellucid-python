# Performance, scaling, and lazy loading

This page explains what makes Cellucid fast (or slow) and how to choose a workflow that scales.

If you are new, you can treat this as “rules of thumb + common failure modes”.

## At a glance

**Audience**
- Computational users working with large datasets
- Anyone hitting slow loads, memory spikes, or browser crashes

## The three biggest performance drivers

1) **Number of cells (`n_cells`)**
   - affects point rendering, selection, and most UI interactions
2) **Number of genes (`n_genes`)**
   - affects gene search lists and the cost of gene expression queries
3) **Data layout**
   - a gene-major (CSC) matrix is one read per gene; a row-major (CSR) one is
     indexed by gene first

## Prepared or not: performance tradeoffs

### A Cellucid-ready AnnData (`prepare()`)

Best for:
- repeated viewing,
- sharing and hosting,
- very large datasets.

Why it’s fast: the expression matrix is gene-major, every array is chunked for
range reads, and the codecs are ones browsers decode natively. On a
561,947-cell dataset the first gene took 86 ms from the prepared `.h5ad`
against 4,110 ms from the same data written by anndata's defaults, over a
simulated 80 ms / 50 MB/s link
({doc}`../c_data_preparation_api/10_performance_and_storage`).

### Any AnnData (best convenience)

Best for:
- interactive notebook exploration,
- quick sanity checks,
- viewing a large `.h5ad` or Zarr store as it is.

How it works: the viewer reads the store itself, in place, with HTTP range
requests. The Python server reads the store's metadata once at startup, then
only serves bytes; it computes nothing per request and holds no matrix in
memory. The browser holds what it has read.

## Lazy loading (what it means in practice)

### The problem

Gene expression is huge (`n_cells × n_genes`). Loading it all into browser memory is often impossible.

### The solution

Every route reads gene expression **on demand**:
- the viewer fetches only the genes you request,
- an embedding, an `obs` column, the neighbour graph, or a vector field is read
  when a view needs it.

### Which modes are lazy?

- `.h5ad` or Zarr store served by Python: ✅ the browser reads the store in
  place with range requests (so are `.h5mu`, MuData Zarr, Loom and 10x files)
- An R object (`.rds`, `.rda`) or a MatrixMarket folder: ⚠️ read once from the
  start, because these formats allow no random access
- In-memory `AnnData` served by Python: ✅ served in priority order, each chunk
  made from the object when the viewer asks for it, so a subsample transfers
  only its cells
- Browser file picker for `.h5ad`: ✅ lazy (read from disk in a worker)

See every route: {doc}`02_loading_routes`.

## What direct AnnData costs, measured

Measured on one Mac, on localhost, against a 562k-cell × 8,192-gene `.h5ad`
whose expression matrix is row-major (CSR) and gzip-compressed:

| Step | Time |
| --- | --- |
| Server startup (`[1/4]` to `[4/4]`, excluding the web build) | 0.16 s |
| `obs` fields | ~90 ms |
| UMAP, 2D or 3D | ~40 ms each |
| Gene index (row-major matrix), in the background from the moment the dataset opens | about 2 s |
| First gene, shown after the index is built | a few ms |
| First gene, asked for at once | 1.3 s (it waits for the index) |
| Every later gene | under 1 ms |
| 100 genes at once | ~20 ms |
| 1,000 cell rows | 0.55 s |
| Neighbour graph, when a view draws it | 0.65 s |

A row-major matrix stores cells, not genes, so the browser builds a gene index
over it once. The browser measures the link to the server while the dataset
opens: from the same machine or network it reads as fast as a local disk, so
the index starts at once, in the background, and is usually ready before the
first gene is chosen. Behind a slow link (a remote JupyterHub over the
internet) it starts with the first gene instead, so opening downloads only what
the first view draws. A column-major (CSC) matrix needs no index.

A file whose cells are not in priority order is served the same way, each
chunk made from the file's arrays when it is asked for. Rows in priority order
lie all over the file, so an array stored in chunks (compressed, or any sparse
matrix anndata wrote) is decoded once, each source chunk the first time a row
in it is asked for, into a memory-mapped file in the temporary folder, and a
served chunk's rows are taken from it with one indexing step; an uncompressed
array stored in one piece is read where it lies. Every chunk of each array,
read once in order, on one Mac (1,000,000 cells × 2,000 genes, 4.5 million
stored values, gzip `.h5ad`; dense: 200,000 × 100), in the server's process
and over HTTP on localhost:

| Served | Before | Now, in process | Now, over HTTP |
|---|---|---|---|
| `X`, cell-major (CSR) file | 7 MB/s | 109 MB/s | 75 MB/s |
| `X`, gene-major (CSC) file | 28 MB/s | 62 MB/s | 46 MB/s |
| `X`, dense file | 126 MB/s | 143 MB/s | 128 MB/s |
| `obsm` (UMAP, 30-column PCA), 1M cells | 0.5 MB/s | 148 MB/s | 124 MB/s |
| `obs` | 2,300 MB/s | 2,300 MB/s | 880 MB/s |
| `X` of an in-memory object | — | 1,970 MB/s | 263 MB/s |

The first pass includes decoding the file's gzip chunks; the served bytes are
the same as before.

An in-memory `AnnData` is served as a store with its cells in priority order,
each chunk made from the object when the viewer asks for it, so a subsample
the viewer opens transfers only its cells. Measured on one Mac over localhost,
with the web app's own reader, for a 1,000,000-cell × 2,000-gene object
(100 million stored values, CSR, `obs` with a categorical and a float column,
a 2-D UMAP and a 30-column PCA), against the temporary `.h5ad` this package
served before:

| Opened | Bytes sent (open, a field, the UMAP, two genes) | Requests |
|---|---|---|
| 1 % of the cells, served in priority order | 0.18 MB | 19 |
| 1 %, before (`blocks`) | 2.1 MB | 168 |
| 1 %, before (`uniform`) | 10.6 MB | 23 |
| 10 %, served in priority order | 1.4 MB | 19 |
| 10 %, before (`blocks` or `uniform`) | 10.6 MB | 23 |
| every cell, served in priority order | 14.2 MB | 39 |
| every cell, before | 10.6 MB | 23 |

Opening every cell costs 4 MB more: the cells' original positions, which the
viewer reads to name every cell by its row in your object. The server starts
in 0.4 s (before: 1.9 s to write the temporary file) and lays the cell-major
`X` out gene-major in about 2 s (measured when the layout ran in the
background; it now runs before the call returns, the object's snapshot:
{doc}`08_anndata_mode_show_anndata_and_serve_anndata`), after which a gene
takes a few milliseconds; the notebook widget answers the same requests with
the same bytes through the kernel.

Your numbers depend on the machine, the disk, and the network between browser
and server.

## Startup reads metadata only

```bash
cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

Startup prints four numbered steps — `[1/4] Detecting format`,
`[2/4] Reading metadata`, `[3/4] Checking launch choices`,
`[4/4] Starting server` — and none of them reads an array: the server reads
cell and gene counts and the names of the stored elements, checks the launch
choices you set against them, and binds. The transcript is in
{doc}`../../web_app/b_data_loading/04_server_tutorial`. A long pause inside
`[4/4] Starting server` is the viewer's web build being fetched and verified.

Inputs that add work:

- an in-memory `AnnData`, or an AnnData file not in priority order: its `obs`,
  `var`, `uns` and the other small elements are written once in memory, in
  priority order, before the bind; an object's arrays are copied once, before
  the bind too (its snapshot: sparse expression matrices laid out gene-major,
  stored as `prepare()` stores them);
- an in-memory `MuData`, written once to a temporary `.h5mu` by `prepare()`;
- a store holding an array stored with a codec the browser cannot decode (an
  HDF5 filter such as Blosc or Zstandard, a Zarr codec such as bz2): each chunk
  is decoded on the kernel when it is asked for. Step 1 names those arrays.

## The neighbour graph is read when drawn

`knn_key` (`--knn`) names the `obsp` graph the viewer opens on; without it the
viewer uses `obsp['connectivities']` when the store has it. Either way the
graph costs nothing at startup: the browser reads it the first time a view
draws edges (0.65 s for the measured dataset above).

## Practical recommendations by dataset size (rule of thumb)

```{note}
Exact cutoffs depend on your machine, browser, and GPU. These are starting points.
```

- **Small (≤50k–100k cells):**
  - Direct AnnData viewing is fine and very convenient.
- **Medium and large (100k–millions):**
  - prefer a path (`show_anndata("data.h5ad", dataset_name="My dataset", dataset_id="my-dataset")`)
    over an in-memory object, so nothing is written first.
  - a column-major (CSC) matrix avoids the first-gene index build.
  - `prepare()` once pays off if you share the dataset or view it often.

## Browser/GPU considerations (web app side)

Cellucid’s rendering load is ultimately in the browser (WebGL).

Common symptoms of GPU pressure:
- “WebGL context lost”
- extreme slowdowns when switching to high-detail modes
- browser tab crashes

For web-app-specific performance tuning (render modes, LOD, etc.), see:
{doc}`../../web_app/n_benchmarking_performance/index`.

## Troubleshooting performance

If performance is unexpectedly bad:

- check which route opened it ({doc}`02_loading_routes`)
- the first gene of a large row-major (CSR) matrix builds the gene index once;
  later genes are fast
- read the numbered startup steps: a long `[1/4]` is an in-memory MuData being
  written to its temporary `.h5mu`
- over an SSH tunnel, run the server on the machine that holds the file
- consider writing a Cellucid-ready AnnData once (`prepare`) and opening that

More symptom-based help: {doc}`15_troubleshooting_viewing`.
