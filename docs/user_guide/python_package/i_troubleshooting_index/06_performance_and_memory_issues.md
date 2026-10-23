# Performance and memory issues

First identify the phase that is slow:

1. Python preparation and validation
2. writing a Cellucid-ready AnnData (`prepare()`)
3. direct-AnnData reading, where the browser reads the `.h5ad`/Zarr store in
   place with HTTP range requests
4. network or disk transfer
5. browser parsing and GPU upload
6. interaction, filtering, analysis, or figure export

Phase 3 exists only when you serve an AnnData object or an `.h5ad` / Zarr
path directly; the Python server computes nothing there, and the browser's
reader does the work. It is also where a gene or column holding an infinity is
refused — see {doc}`../d_viewing_apis/15_troubleshooting_viewing`.

Record dataset shape, sparse or dense representation, artifact sizes, browser,
viewport, number of visible views, render mode, and the exact operation being
timed.

## Direct AnnData startup looks like it has hung

Read the last line the terminal printed. A direct AnnData server runs four
numbered steps and reads only metadata before it binds, so each step is short
except two that are labelled:

| Last line | What is running | What to do |
| --- | --- | --- |
| `Serving: a decodable copy, written once to a temporary .h5ad` | The store holds arrays the browser cannot decode (named on the line above); anndata reads it once and writes a copy | Rewrite the store with codecs the browser decodes (gzip or LZF for HDF5) |
| `Format: written once to a temporary .h5ad` | An in-memory object is being written to the temporary `.h5ad` | Serve the file from disk instead |
| `[4/4] Starting server...` | The configured viewer generation is being fetched and verified | Wait once; the cache is reused afterwards |

The transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

## The neighbor graph is slow to draw

The browser reads an `obsp` matrix only when a view
draws the neighbor graph; nothing about it costs startup time. Measured on one
Mac over localhost, the graph of a 562k-cell `.h5ad` took 0.66 s to read when
drawn. The cost tracks *stored neighbors*, not cells: a 50-neighbor graph over
five million cells holds 250,000,000 to 500,000,000 of them.

`prepare(..., connectivities=...)` checks a graph once, when the file is
written: the exact symmetry check compares the matrix against its transpose, so
two further copies of the graph are alive at once. If that runs out of memory,
subset the object first or leave the graph out (`connectivities=None`).

## Repeated sessions over one very large object

Serving an `.h5ad` directly costs a fraction of a second at startup, so
repeated sessions pay little. What a session pays is the browser's reading: for
a row-major (CSR) matrix, the first gene builds a gene index once per session.
Run `prepare(...)` once to write the matrix column-major (CSC), which removes
that cost, especially when the dataset is to be shared or hosted.

## Python memory grows during preparation

Preserve sparse matrices when the API accepts them, avoid unnecessary dense
copies, and verify expression orientation before preparing.

## Python memory grows while serving an AnnData directly

A server for a path reads no array: it serves the file's bytes from disk, so
its memory does not grow with the genes a session visits. Memory is spent:

- by your own `adata`, when you pass an in-memory object; it is written once to
  a temporary `.h5ad` (sparse matrices column-major, uncompressed) and then
  served from disk, so the temporary file needs disk space for one
  uncompressed copy;
- once, at startup, when a store holds arrays the browser cannot decode: the
  whole store is read with anndata to write a decodable copy.

Gene values are held by the browser, not by Python.

## The browser becomes slow

Use the in-app Performance Benchmark and browser performance tools to distinguish
frame-rate pressure from long CPU tasks or slow requests. Additional views and
large pixel dimensions add independent rendering work.

## Figure export is slow or large

Record format, dimensions, DPI, visible point count, and whether all views are
included. Test the smallest scientifically sufficient output first.

See {doc}`../b_concepts_mental_models/07_performance_mental_model_and_scaling`,
{doc}`../c_data_preparation_api/10_performance_and_storage`, and
{doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.
