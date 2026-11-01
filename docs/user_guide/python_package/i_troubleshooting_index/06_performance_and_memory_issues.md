# Performance and memory issues

First identify the phase that is slow:

1. Python preparation and validation
2. writing a Cellucid-ready AnnData (`prepare()`)
3. direct reading, where the browser reads the served file or folder itself:
   in place with HTTP range requests (`.h5ad`, Zarr, `.h5mu`, `.loom`, 10x), or
   once from the start (an `.rds`/`.rda`, a MatrixMarket folder)
4. network or disk transfer
5. browser parsing and GPU upload
6. interaction, filtering, analysis, or figure export

Phase 3 exists only when you serve an object or a path directly
(`cellucid serve`, `serve_anndata`, `show_anndata`); the Python server computes
nothing there, and the browser's reader does the work. It is also where a gene or column holding an infinity is
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
| `Format: a Zarr store in priority order, made as the browser reads it` (an in-memory object) | Its expression matrices are counted, then laid out gene-major in the background; the viewer opens at once, and a gene read before the layout ends is made from the layout's progress (what is placed copied, the rest found by one scan of the object's column indices; slower than after the layout, which is waited for only when it ends sooner) | Wait for the layout, or serve a Cellucid-ready file (`prepare()`), which needs none |
| `[4/4] Starting server...` | The configured viewer generation is being downloaded and verified, every asset, at every server or viewer start | Give the machine fast outbound access to the web source; on a machine without it, serve with `--no-web-ui` |

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

A server serves a file in priority order (one `prepare()` wrote) as its bytes
from disk, and anything else as a store in priority order whose chunks it
makes when they are asked for, from the object or from the file's own arrays
(and a bounded cache of renumbered genes and chunk entries: a 64th of the
machine's memory, 256 MiB on a 16 GB machine, as cellucid-r's served stores
keep; `cache_bytes=` of {func}`~cellucid.serve_anndata`,
{func}`~cellucid.show_anndata` and {class}`~cellucid.AnnDataServer` sets it,
`0` keeping nothing). A file's
array stored in chunks (compressed, or any sparse matrix anndata wrote) is
decoded once, a source chunk at a time as rows in it are first asked for,
into a memory-mapped file in the system temporary folder, so rows in priority
order are taken from it at memory speed: that needs disk space for the
decoded arrays the viewer reads (removed when the viewer stops), and its
pages are file cache the system can reclaim, not Python memory. An
uncompressed array stored in one piece is read where it lies. Memory is spent:

- by your own `adata`, when you pass an in-memory object: the store reads its
  arrays and copies none of them, but a cell-major (CSR) sparse `X` is laid
  out gene-major once, in the background from the start (a layer when it is
  first read), into memory-mapped files in the system temporary folder, which
  need disk space for one uncompressed copy of the matrix and are removed when
  the viewer stops; a MuData is written to a Cellucid-ready `.h5mu` by
  `prepare()`;
- by the obs table and the other small elements (`var`, `uns`, `varm`,
  `varp`), held in memory in priority order.

Gene values are held by the browser, not by Python.

## The browser becomes slow

Use the in-app **Benchmark** panel and browser performance tools to distinguish
frame-rate pressure from long CPU tasks or slow requests. Additional views and
large pixel dimensions add independent rendering work.

## Figure export is slow or large

Record format, dimensions, DPI, visible point count, and whether all views are
included. Test the smallest scientifically sufficient output first.

See {doc}`../b_concepts_mental_models/07_performance_mental_model_and_scaling`,
{doc}`../c_data_preparation_api/10_performance_and_storage`, and
{doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.
