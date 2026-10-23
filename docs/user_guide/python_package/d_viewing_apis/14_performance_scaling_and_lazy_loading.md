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
  place with range requests
- In-memory `AnnData` served by Python: ✅ written once to a temporary `.h5ad`,
  then read like a file
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

An in-memory `AnnData` is written once to a temporary `.h5ad` with sparse
matrices column-major and uncompressed. For the same dataset that write took
0.8 s and 925 MB on disk; the first gene then took 8 ms, and 100 genes 68 ms.

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

Two inputs add work before the bind:

- an in-memory `AnnData`, written once to a temporary `.h5ad`;
- a store holding an array stored with a codec the browser cannot decode (an
  HDF5 filter such as Blosc or Zstandard, a Zarr codec such as bz2), read once
  with `anndata` and its codec plugins and served as a temporary `.h5ad`. Step 1
  names those arrays. Rewriting the store with gzip (or no compression) avoids
  the copy.

## The neighbour graph is read when drawn

`knn_key` (`--knn`) names the `obsp` graph the viewer opens on; without it the
viewer uses `obsp['connectivities']` when the store has it. Either way the
graph costs nothing at startup: the browser reads it the first time a view
draws edges (0.66 s for the measured dataset above).

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
- read the numbered startup steps: a long `[1/4]` with a `Browser cannot decode`
  note is the one-time decodable copy
- over an SSH tunnel, run the server on the machine that holds the file
- consider writing a Cellucid-ready AnnData once (`prepare`) and opening that

More symptom-based help: {doc}`15_troubleshooting_viewing`.
