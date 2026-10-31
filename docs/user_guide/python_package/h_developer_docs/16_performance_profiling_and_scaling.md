# Performance profiling and scaling

This page explains where Cellucid Python performance bottlenecks typically come from and how to profile them.

It covers:
- preparation-time performance (`prepare`)
- server performance
- notebook performance (embedding + hooks)

---

## Performance “north star”: what users experience

Users mostly care about:

- time to first view (TTFV): “how long until I can rotate and color-by?”
- interaction latency: “does selecting/filtering feel instant?”
- memory footprint: “does this crash my laptop/server?”

Developer work should aim to improve one without catastrophically harming the others.

---

## Preparation-time performance (`prepare`)

### Main cost drivers

1) **Converting the expression matrix** to CSC float32 (a full copy of its
   stored values).
2) **Compressing chunks**: both containers compress chunks on every core. An
   `.h5ad`'s numeric chunks are deflated by a thread pool and stored with
   HDF5's direct chunk write (`_write_chunks_in_parallel` in `_store.py`),
   because HDF5's own filters compress one chunk after the other on the
   calling thread.
3) **Symmetry check of a large graph** (`connectivities` against its transpose).

### Practical tuning knobs

- `gene_identifiers=[...]` reduces size and time if you only need a gene panel.
- Subset cells before the call.
- There are no codec or precision knobs: the layout is the measured one
  ({doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`).

### Memory footguns

- Dense gene expression matrices can explode memory.
- Sparse matrices are preferred; `prepare` converts them to CSC for per-gene access.

---

## Server performance

The Python server computes nothing per request: it answers
range requests with the store's bytes from disk, on a threaded HTTP/1.1 server
with keep-alive connections. The reading — and so the performance — is the web
app's readers (`cellucid/assets/js/data/anndata-file-source.js`,
`data/anndata/`, `base-anndata-adapter.js`). Formats read with range requests
(AnnData, MuData, Loom, 10x) transfer only what a view needs; an R object or a
MatrixMarket folder is read once from the start, so its whole size crosses the
link before the first view.

What decides the cost:
- matrix layout: a column-major (CSC) matrix reads one gene as one contiguous
  range; a row-major (CSR) matrix costs a gene index the browser builds once,
  on the first gene,
- storage and latency between the server and the data (network filesystems,
  tunnels),
- an in-memory object, or a file not in priority order, served as a store in
  priority order: its small elements written in memory at startup, every chunk
  made on request (a store with browser-undecodable codecs is decoded on the
  kernel so), an object's arrays copied once at startup into memory-mapped
  files (its snapshot; sparse expression matrices laid out gene-major, stored
  as `prepare()` stores them); an in-memory MuData served the same way, as a
  MuData Zarr store, with no file written.

Measured on one Mac over localhost, for a 562k-cell × 8,192-gene `.h5ad` (CSR,
gzip): startup 0.16 s, obs fields about 90 ms, 2D/3D UMAP about 40 ms each,
the gene index about 2 s in the background from the moment the dataset opens
(the browser measures a localhost link as fast as a disk), then a gene in a few
ms and 100 genes in about 20 ms (1.3 s for a gene asked before the index is
built), 1,000 cell rows 0.55 s, neighbour graph 0.65 s when drawn. The same object in memory: 0.8 s to
write the temporary CSC `.h5ad` (925 MB), then first gene 8 ms and 100 genes
68 ms.

Practical guidance:
- If users complain “gene coloring is slow”, profile the browser's reads in the
  network tab first: count the range requests per gene and check the matrix
  layout.

---

## Notebook embedding performance

Notebook mode adds overhead:

- notebook layout around the widget's output
- every read crosses the notebook's connection (the widget relays it to the
  kernel and coalesces requests that arrive together; while a cell runs, reads
  wait until it ends)
- large hook payloads (huge selections)

Guidance:
- keep hook payloads lightweight when possible,
- avoid sending large “command” payloads frequently (e.g., repeated highlight updates every frame).

---

## Profiling playbook (developer steps)

### Step 1 — Measure first

Before optimizing, capture:
- dataset size (`n_cells`, `n_genes`)
- whether `prepare()` wrote the file, and its container (`.h5ad`, `.zarr`, `.zarr.zip`)
- rough timings (“prepare took 20 minutes”, “gene color-by takes 3 seconds”)

### Step 2 — Profile the Python side

Useful built-in tools:

```bash
python -m cProfile -o profile.out -m your_script_that_calls_prepare
python -c "import pstats; p=pstats.Stats('profile.out'); p.sort_stats('cumtime').print_stats(30)"
```

For memory:
- prefer tooling like `tracemalloc` in targeted code paths (avoid broad guesses).

### Step 3 — Profile HTTP behavior

In browser DevTools Network:
- filter `_cellucid/anndata/` and check which ranges or Zarr keys dominate
  (expression chunks, `obs` columns, `obsm` arrays, the graph),
- check whether requests are repeated (a store answer should revalidate with
  `304`, not download again).

### Step 4 — Decide optimization strategy

Common strategies:
- reduce data size (subset cells or genes)
- reduce recomputation (cache)
- reduce per-request overhead (batch requests, avoid repeated conversions)

---

## Edge cases

### Millions of cells

Expect:
- large point arrays
- heavy browser memory pressure
- a row-major (CSR) matrix: the gene index the browser builds on the first gene

Mitigation:
- recommend writing the file with `prepare()` (gene-major)
- consider downsampling workflows in documentation

### Extremely high category counts

Even when a file is written, categories with tens of thousands of levels can become unusable in the UI.

Mitigation:
- document recommended practice (collapse categories, use higher-level labels).

---

## Troubleshooting

### Symptom: “The prepared file is enormous”

Likely causes:
- every gene kept,
- a dense expression matrix,
- integer counts above 2^24 (stored as float64).

Fix:
- keep a gene panel with `gene_identifiers=[...]`,
- pass a sparse matrix.

### Symptom: “The server uses too much memory”

Likely causes:
- an in-memory object is written to a temporary CSC `.h5ad` or `.h5mu` (a copy of its matrices),
- caches storing large arrays.

Fix:
- serve a path rather than an in-memory object for large datasets,
- reduce gene interactions, or restart kernel to clear caches after exploration.
