# Performance and memory issues

First identify the phase that is slow:

1. Python preparation and validation
2. writing or reading export artifacts
3. direct-AnnData serving, where one request materializes a dense float32
   column, validates it, and caches it
4. network or disk transfer
5. browser parsing and GPU upload
6. interaction, filtering, analysis, or figure export

Phase 3 exists only when you serve an AnnData object or an `.h5ad` / `.zarr`
path directly; a prepared export has done that work already and reads finished
artifacts instead. It is also the phase that answers `422` rather than serving
when a column is not entirely finite — see
{doc}`../d_viewing_apis/15_troubleshooting_viewing`.

Record dataset shape, sparse or dense representation, artifact sizes, browser,
viewport, number of visible views, render mode, and the exact operation being
timed.

## Direct AnnData startup looks like it has hung

It probably has not. Read the last line the terminal printed. A direct AnnData
server runs five numbered steps and names each part of the work as it starts, so
a silent stretch always has a label above it. A prepared export runs three steps
instead, because it has no adapter to build.

Startup prints five numbered steps; the transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

Match the last line printed against the phase it names:

| Last line | What is running | What to do |
| --- | --- | --- |
| `Connectivity: reading obsp['connectivities']` | The neighbor graph is being validated and deduplicated into edges | Drop `serve_connectivity=True` / `--connectivity` if the session will not draw the graph |
| `Obs columns: classifying N` | One validation pass per continuous obs column | Narrow the served set with `obs_keys=` / `--obs-key` |
| `[4/5] Building manifests...` | One centroid per category per available dimension | Reduce the number of categorical fields served, or export once |
| `[5/5] Starting server...` | The configured viewer generation is being fetched and verified | Wait once; the cache is reused afterwards |

## Serving the neighbor graph is slow or runs out of memory

`adata.obsp['connectivities']` is not read at all unless you pass
`serve_connectivity=True`, or `--connectivity` on the command line. It is off by
default because on a large object it is the longest part of startup, and most
sessions never draw the graph.

The scale is worth knowing before you ask for it. A 50-neighbor graph over five
million cells is 250,000,000 directed neighbor relations; symmetry stores the
reverse of every relation that is not already mutual, so the matrix holds
250,000,000 to 500,000,000 **stored neighbors**. The served edge list is the
strict upper triangle, so exactly half of that: 125,000,000 to 250,000,000
**deduplicated edges**.

The memory goes to the passes, not to the result:

- the exact symmetry check in topology and weight materializes the transpose and
  compares it against the original, so two further copies of a matrix already
  costing about 12 bytes per stored neighbor are alive at once;
- duplicate-coordinate detection sorts one 64-bit key per stored neighbor (2 GB
  at 250,000,000) next to 24 bytes per stored neighbor of row, column, and
  weight arrays (6 GB);
- the ordering pass holds a 64-bit permutation per edge (1 GB at 125,000,000
  edges) plus the reordered arrays alongside the originals (about 3 GB more);
- what survives is 16 bytes per edge — a 4-byte source index, a 4-byte
  destination index, and an 8-byte Float64 weight — so 2 GB at 125,000,000
  edges, served from `connectivity/edges.src.bin`,
  `connectivity/edges.dst.bin`, and `connectivity/edges.weights.f64.bin`.

Above **50,000,000 stored neighbors**, which is exactly a 50-neighbor graph over
one million cells, the `cellucid.connectivity` logger names both counts and the
cost before the validation starts:

```text
Reading a neighbor graph of 250,000,000 stored neighbors over 5,000,000 cells. Validating it is symmetric, deduplicating it into edges, and ordering them takes minutes and several times the graph's own memory at this size. Serve without connectivity to skip this work entirely.
```

The whole graph is read before the server binds its port, because
`connectivity_manifest.json` is among the first documents the viewer fetches and
it declares the edge count and the neighbor maximum. There is nothing to defer.

If the pass runs out of memory, the options are: serve without connectivity,
subset the object, or export once with `prepare(..., connectivities=...)` and
serve the export.

Asking for connectivity when `adata.obsp` holds no `connectivities` matrix is an
error naming `sc.pp.neighbors(adata)`, not a silent absence. `Connectivity: no`
in the step 3 report means the object carries no graph; `Connectivity: not
served (...)` means it carries one that this server was not asked for.

## Repeated sessions over one very large object are slow every time

Serving an `.h5ad` directly rebuilds the adapter and recomputes every centroid
on each startup. Run `prepare(...)` once and serve the resulting export from
then on: the export pays that work a single time, and the prepared server reads
finished artifacts in three steps rather than five. For an object you open more
than a couple of times, this is faster than serving the `.h5ad` every time, and
the difference grows with the number of cells.

## Python memory grows during preparation

Preserve sparse matrices when the API accepts them, avoid unnecessary dense
copies, and verify expression orientation before export.

## Python memory grows while serving an AnnData directly

Serving colours genes on demand, and each gene requested is one dense float32
column — four bytes per cell, so 288 KB on a 72,000-cell object and 73 MB on an
18,142,044-cell one. Those columns are cached so that revisiting a gene is
instant, and the cache is bounded by **bytes**, not by a count of entries:
`GENE_EXPRESSION_CACHE_BYTES` is 256 MiB, and `LRUCache` evicts the
least-recently-used columns until a new one fits.

The unit matters more than the number. A budget of a hundred entries is 29 MB
on a 72,000-cell object and 7.3 GB on an 18,142,044-cell one — unbounded in
practice on exactly the size of object the cache exists to help. Counting the
bytes that actually run out makes the ceiling the same everywhere.

What that budget buys:

- roughly 32 cached columns of a 2,000,000-cell object, or 3 of an
  18,000,000-cell one;
- a column larger than the whole 256 MiB budget — over about 67,000,000 cells —
  is served normally and simply not cached, because a cache that refused to
  serve what it cannot hold would be deciding what the server may publish;
- steady-state resident memory for gene serving that does not grow with how
  many distinct genes a session visits.

If serving still runs out of memory, the object itself is the cost, not the
cache: an `.h5ad` served read-only-backed reads columns from disk, while an
in-memory object and an eagerly loaded Zarr store hold the whole matrix. Run
`prepare(...)` once and serve the export instead.

## The browser becomes slow

Use the in-app Performance Benchmark and browser performance tools to distinguish
frame-rate pressure from long CPU tasks or slow requests. Additional views and
large pixel dimensions add independent rendering work.

## Figure export is slow or large

Record format, dimensions, DPI, visible point count, and whether all views are
included. Test the smallest scientifically sufficient output first.

See {doc}`../b_concepts_mental_models/07_performance_mental_model_and_scaling`,
{doc}`../c_data_preparation_api/10_performance_tuning_guide_prepare_export`, and
{doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.
