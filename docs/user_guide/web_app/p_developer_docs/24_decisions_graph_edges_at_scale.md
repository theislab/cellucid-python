# Decision record: graph edges at scale

This page answers: **"How does the web app turn a neighbour graph (`obsp`)
into the edges a pane draws and the rows the k-NN tool expands through, at
10M–50M cells and on a subsample, and why is it built that way?"** It records
the current design, the decisions behind it with their dates, the
measurements that justify them, what was tried and taken out, and what is
still open. Update this page in the same change as the code it describes.
How the K shown edges reach the GPU is in
{doc}`07_rendering_pipeline_webgl_and_performance_notes`; how a share's rows
are chosen and priced is in {doc}`23_decisions_partial_reading`.

## At a glance

**Audience**
- Developers changing the graph reader, the edge builder,
  the edge layer or the k-NN tool.

**Time**
- 15 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (where the reader worker and the page
  meet) and {doc}`09_data_loading_pipeline_and_caching` (the long form of
  "Graphs and embeddings arrive finished").

:::{note}
This page describes the code as of 2026-10-30.
:::

---

## What it does now

Paths are under `cellucid/assets/js/`.

1. **The graph is read line by line** (`data/anndata/graph-lines.js`). A line
   is a row of a CSR matrix or a column of a CSC one (a symmetric matrix is
   its own transpose). Lines arrive in slabs of about `SLAB_BYTES` of stored
   entries; the next slab is asked for before the current one is handed out,
   so a 10M-cell graph (about 2 GB stored) is never held whole.
   `sparseGraphLines` reads a stored CSR/CSC group, `denseGraphLines` a dense
   `n × n` array by row slab, `placedGraphLines` a MuData modality's graph in
   its own rows.
2. **A subsample's graph is its induced subgraph, read in stored order**
   (`inducedGraphLines`). The open cells' stored rows are sorted once
   (`storedRowOrder`, a radix order from `utils/radix-order.js`), and each
   kept entry's column becomes its row's rank among the open rows, so the
   lines ascend. The reader picks one of two reads
   (`AnnDataReader._graphReadPlan`, a `GraphReadPlan`):
   - `stream`: the graph from the first open row to the last, in contiguous
     slabs, keeping the open rows;
   - `rows`: the open rows' pointers and entries only, through the chunks
     that hold them.

   `stream` is chosen exactly when it downloads and decodes no more bytes
   than `rows` (`_readsNoMoreThanRows`, priced by `firstGeneReadCost` with
   the link's merge gap; if every gap between open rows is shorter than a
   chunk, the answer is known without pricing). Over a link slower than a
   disk, slabs are read ahead, bounded by the decoded-chunk share of memory
   (`_graphReadAhead`).
3. **The edges are built in one checked pass** (`data/connectivity-edges.js`,
   `createEdgeBuilder` via `connectivityEdgesFromLines`). As line `major`
   arrives, its upper entries are its forward edges; each lower entry must
   equal the next due edge of its partner's line. That one pass checks bounds,
   order, finite non-negative weights and **exact** symmetry of topology and
   weight. A diagonal entry (Seurat's self-neighbour) and a stored zero are
   not edges.
4. **What leaves the worker is the forward adjacency** (`ConnectivityEdges`):
   `forwardOffsets` (`Uint32Array`, 4 bytes a cell), `partners` (`Uint32Array`)
   and `strengths` (`Float32Array`, `weight / maxWeight`), 8 bytes an edge.
   The three buffers are transferred to the page, not copied. Nothing is kept
   in draw order.
5. **A subsample's edges are placed on the open cells** (`relabelledEdges`):
   each edge is counted at its smaller cell as its line arrives (the
   builder's `cellOffsets`, `createEdgeBuilder` `relabel`), then prefix sum,
   scatter (typed loops), then sort each cell's few partners. The output is
   bit-identical to building from the cells' own lines. Strengths and
   partners are placed in one walk when the four edge arrays fit what every
   stored line's edges take to build (`wholeBuildBytes`); else strengths are
   placed and released before partners (each partner made its cell in place
   on the first walk), so about 1.5 times the edges' bytes are held while it
   runs (`connectivityEdgeBytes` `relabelCells`, the cells' offsets counted
   from the build on).
6. **Memory is counted before anything is allocated.** `connectivityEdgeBytes`
   gives `build` and `held`. If `build` is more than
   `memoryAllowanceBytes()`, the build is refused at once with a
   `RangeError` (`code: 'TOO_LARGE'`) that names both sizes. While the build
   runs, it is the `'graph build'` part of the working set
   (`data/memory-budget.js`), so every thread's caches make room. On the
   page, the published edges are the `'graph'` part (`app/main.js`).
7. **A pane shows the first K visible edges in a fixed random order**
   (`rendering/connectivity-edge-selection.js`). The order is a keyed
   four-round Feistel bijection of `[0, E)`; the source of an edge is found by
   binary search inside a 128-edge block. The first 4K ranks are kept as a
   shown-order prefix. Drawing is in `rendering/connectivity-edge-layer.js`;
   details are on page 07.
8. **The k-NN tool reads rows, not the whole graph.** `KnnExpansion`
   (`rendering/selection-gesture-queries.js`) asks a `KnnRowSource` for
   65,536 frontier cells at a time (`KNN_ROW_BATCH_CELLS`), with the next
   batch in flight. In `app/main.js`, `knnRowReader` answers from the file
   (`AnnDataReader.getConnectivityRows`, which streams the graph when a large
   frontier's rows lie in every chunk) or, **while a view holds the edges**,
   from them (`forwardEdgeRows`: two passes over every edge). It picks
   whichever is faster by rates measured on that page: milliseconds per row
   from the file, edges per millisecond for a pass. `KNN_SCAN_EDGES_PER_MS`
   (100,000) is used only until the first pass is timed, and the seed's rows
   always come from the file. The tool keeps only a record that the graph was
   checked (`checkedConnectivity`); it keeps none of the edges.

## Decisions and why

| Date | Decision | Reason | Rules out |
|---|---|---|---|
| before 2026-10-29 | Shown edges follow a stateless keyed bijection (Feistel, cycle walking), not a stored shuffle | O(√E) memory at any size; every K-prefix is a uniform sample, the same in every session | A per-edge permutation array (4 bytes an edge) |
| 2026-10-28 | A build that cannot fit is refused before any array is allocated, naming its bytes | At 50M an earlier version allocated gigabytes and then failed at its last array with "Array buffer allocation failed" and no size | Partial builds; a fixed edge-count cap (there is none) |
| 2026-10-29 | Hold only the forward adjacency: 8 B/edge + 4 B/cell | The edges in draw order (12 B/edge) and every cell's neighbour list (8 B/edge) were removed; at 50M × k = 15 the earlier layout needed 10.4 GB to build and 7.4 GB to hold | Any per-edge array beyond partners and strengths on either thread |
| 2026-10-29 | The k-NN tool reads frontier rows from the file in the reader worker | No thread holds every cell's neighbour list; works for every format and source (`.h5ad` CSR/CSC/dense, Zarr, Loom, MuData, `.rds`, served stores, links, providers, subsamples) | A page-side neighbour index |
| 2026-10-29 | A subsample's induced subgraph streams into the builder | On the full Human Lung Cell Atlas (2.3M cells, 102.9M stored entries), 1.3 GB was held before the build began | Gathering the induced subgraph first |
| 2026-10-29 | 100M graphs: symmetry checked by exact per-cell counts plus a keyed 128-bit hash summed order-independently in both directions (designed, not built) | A disk-backed build cannot afford the random access of the exact check; a false pass has probability 2^-128 | Holding 700M edges in RAM at 100M |
| 2026-10-30 | Read a share's graph in **stored** order, by the cheaper of `rows` and `stream`, and build over ranks before relabelling | In open order, a uniform share was scattered through the file, so every slab spanned the whole graph and decoded its chunks again | Reading in open order; an extra builder that sorts lines |
| 2026-10-30 | Keep the ordered, exact builder for shares too | An order-free builder (keyed pair-hash symmetry) was measured 3–4× slower per edge | The hash builder below 100M |
| 2026-10-30 | Answer a k-NN frontier from the edges held when that is faster, by rates measured on the page | A 65,536-cell frontier took 1.7 s from a 10M-cell file but 0.18 s from 20M held edges | A fixed rule or size threshold |
| 2026-10-30 (owner) | A share never costs more than a larger share or every cell, in any measure; a cost of about s/n (share s of n cells) is the target to work towards | Partial reading has to be a speed gain on every source | Share plans priced above the whole |
| 2026-10-30 (owner) | Speed comes before memory, as long as memory is counted and explained | Placing a share's edges holds about 1.5× the edges' bytes for a faster, exact build | Trading speed for an unexplained memory saving |

## Measurements that justify it

All measurements were taken on a 16 GB Apple Silicon laptop shared with
other work (±15 % run to run). "Node" means the reader worker's handler run
in Node on a local file; "browser" means the whole app in Chromium on the
Metal GPU.

**Forward adjacency (2026-10-29/30, synthetic k = 15 graphs).**

| Case | Before | After |
|---|---|---|
| 50M cells, 349.7M edges, Zarr v3 sharded, browser edges step | not built ("Array buffer allocation failed") | 40.0 s, 250k shown; `'graph'` working set 2.9 GB |
| 50M, Node build | refused (9.7 GB build over the allowance) | 36.3 s, payload 3.0 GB |
| 10M `.h5ad` (69.9M edges), browser edges step | 16.4 s | 8.1 s |
| 10M, reader worker peak / page peak (browser) | 6.6 GB / 2.5 GB | 3.0 GB / 1.2 GB |
| 10M, Node payload / build garbage left | 1,479 MB / 1.3 GB | 600 MB / 0.2 GB |
| 2M local, payload / peak reachable buffers | 296 MB / 434 MB | 120 MB / 259 MB |
| Full HLCA (51.4M edges), Node build peak | 1.3 GB held before the build began | 0.99 GB peak |

At 10M, the earlier "reader worker holds 5.1 GB, the model counts 2.1 GB"
turned out to be uncollected build garbage, not a leak. After a forced
collection, the worker holds what the model counts.

**Edges at a share (Node, 2026-10-30; equal edges by both reads, checked
against the induced subgraph computed from the stored arrays).**

| File | Share | Before | After | Every cell |
|---|---|---|---|---|
| 5M CSR `.h5ad`, 70M entries | uniform 2.5M | 17.9–19.9 s | 2.1–2.5 s, 44 MB | 1.9–2.5 s, 44 MB |
| same | uniform 0.5M | 3.3–3.7 s | 1.2 s | |
| same | blocks 2.5M / 50k | 3.1 s / 79 ms | 2.1 s, 17 MB / 64 ms, 0.4 MB | |
| same | uniform 2.5M, peak RSS | 1,460 MB | 750–840 MB | 770 MB |
| 10M CSR `.h5ad` | uniform 5M | 35.2 s, 860 MB | 3.8 s, 66 MB | 3.2–3.3 s, 66 MB |
| 5M Zarr v3 | uniform 2.5M | 21.1 s, 19,319 req | 2.88 s, 2,215 req | 2.81 s, 2,215 req |
| 2M `.h5mu` (MuData graph) | uniform 1M | 4.7 s | 0.93 s | 0.89 s |
| Full HLCA over S3 | blocks 1 % | | 41.9 s, 107 MB | 296 s, 703.5 MB |

A uniform share's open rows lie in every chunk of the graph (a 65,536-entry
chunk holds about 4,700 cells, so even a 1 % share has about 47 open rows in
each chunk). That makes the whole graph's bytes the floor for any uniform
share, and the read streams it. Block shares read about s/n of the bytes.

**k-NN rows.** A frontier of 1–14 cells reads in 0.2–8.5 ms at 10M. One batch
of 65,536 scattered cells read from the file took 1.7 s (10M) and 5.3 s (50M);
answered from 20M held edges it took 0.18 s. In the browser (Chromium), a drag
to degree 10 lands 600 ms after the press at 10M and 722 ms at 50M, and the
page holds nothing more afterwards.

## Tried and reverted

| What | When | Why it was taken out | What would have to change |
|---|---|---|---|
| Order-free edge builder (symmetry by a keyed 128-bit pair hash, lines in any order) | 2026-10-30 | 3–4× slower per edge than the ordered build plus relabel | Only the 100M disk-backed build needs it, where random access is impossible |
| Unbounded slab read-ahead over a slow link | 2026-10-30 | The full HLCA's graph over S3 held 1.58 GB of slabs in flight | Now bounded by the decoded-chunk share (`_graphReadAhead`) |
| Per-line JavaScript array and comparator sort for unsorted lines | 2026-10-30 | 1.5 s of a 2.5M-cell share's 20 s | Replaced by typed scratch sorting (`lineSorter`); lines of a share now ascend anyway |
| Builder keeping its capacity arrays until it returned | 2026-10-30 | Held them beside the placed arrays (HLCA: 1,085 MB of buffers against 778 MB counted) | `createEdgeBuilder` now hands its arrays over and keeps none (not re-measured on HLCA) |

## Known limits and open items

- **Edge count.** Edge positions are `uint32`, so a graph above 2^32 − 1
  edges is refused (`requireEdgeCount`). A MuData modality's or an R
  object's graph above 2^31 − 1 stored entries is refused rather than
  wrapping its int32 pointers.
- **100M cells × k = 15 is designed, not built.** That is 700M edges: 5.6 GB
  of forward adjacency next to the points, more than a 16 GB machine holds.
  The design: stream the build to an origin-private file, check symmetry with
  counts plus the hash above, and select the shown edges in the worker.
- **Half shares on fast links.** A uniform 50 % share streams every cell's
  bytes and then places its edges on its cells after the last line (the
  relabel); every cell's build overlaps decoding. With the counting moved
  into the build, typed walks and one walk where memory allows, a 5M-cell
  local `.h5ad`'s 2.5M share took 2.08–2.11 s against every cell's
  2.03–2.05 s (before: 2.2–2.5 s against 1.95–2.0 s; 10M before: 3.8 s
  against 3.3 s). Placing the edges overlapped with the stream is not
  built.
- **HLCA over S3.** On 2026-10-29, every cell's edges took 172 s as about
  2,300 small requests. With stored-order reads, which make fewer and larger
  requests, they took 252–296 s: S3 delivers per connection, so only about 7
  of 42 connections were busy. With request splitting over every connection,
  a re-run took 219 s. The 172 s target is not yet shown again, and the HLCA
  share runs have not been redone on the final code. See
  {doc}`22_decisions_link_model_and_request_planning`.
- **50M reader worker memory.** After collection it holds 2.34 GB against
  0.6 GB counted. The edges are not the cause (the build keeps 2 MB once the
  payload has moved); the gene index built ahead is. See
  {doc}`27_decisions_memory_model`.
- **50M publication.** Main-thread long tasks of 1.2–2.8 s at the edges step
  (block sources, prefix and first count slices, plus the harness's forced
  collection).
- **Browser proof of edges at a share** (all three engines) has not been
  run; the Node measurements above are the evidence so far.

## Known issues (pre-commit)

Defects, regressions and unverified behaviour known when this version was
committed (2026-10-30), each with where it is recorded above. An item is
removed in the change that fixes or verifies it.

- **HLCA over S3:** every cell's edges took 172 s on 2026-10-29 and 252–296 s on 2026-10-30 (see Known limits: S3's per-connection delivery); not yet re-measured on the final tree.
- **Edges at a share:** the browser proof was not run; the numbers above are from Node.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| Lines, induced subgraph, plan, k-NN rows | `cellucid/assets/js/data/anndata/graph-lines.js`, `cellucid/assets/js/data/anndata/reader.js` (`_connectivityLines`, `_graphReadPlan`, `_graphReadAhead`, `getConnectivityRows`) | `cellucid/tests/graph-lines-contract.test.mjs`, `cellucid/tests/edges-at-a-share-contract.test.mjs`, `cellucid/tests/analysis-subsample-differential.test.mjs` |
| Edge builder, relabel, bytes, held-edge rows | `cellucid/assets/js/data/connectivity-edges.js` | `cellucid/tests/connectivity-extraction-contract.test.mjs`, `cellucid/tests/scale-boundaries-contract.test.mjs` |
| Graph summary published by the adapter | `cellucid/assets/js/data/connectivity-summary-contract.js` | `cellucid/tests/connectivity-summary-contract.test.mjs` |
| Radix order of open rows | `cellucid/assets/js/utils/radix-order.js` | `cellucid/tests/scale-boundaries-contract.test.mjs` |
| Shown-edge selection and drawing | `cellucid/assets/js/rendering/connectivity-edge-selection.js`, `cellucid/assets/js/rendering/connectivity-edge-layer.js` | `cellucid/tests/connectivity-edge-selection-contract.test.mjs`, `cellucid/tests/connectivity-edge-layer-contract.test.mjs`, `cellucid/tests/connectivity-lod-membership-contract.test.mjs`, `cellucid/tests/browser/connectivity-edges.spec.mjs`, `cellucid/tests/browser/connectivity-edge-gather-parity.spec.mjs` |
| k-NN expansion and row source | `cellucid/assets/js/rendering/selection-gesture-queries.js` (`KnnExpansion`, `KnnRowSource`), `cellucid/assets/js/app/main.js` (`knnRowReader`, `KNN_SCAN_EDGES_PER_MS`) | `cellucid/tests/highlight-knn-lifecycle.test.mjs`, `cellucid/tests/data-mapping-knn-contract.test.mjs`, `cellucid/tests/browser/knn-missing-connectivity-warning.spec.mjs` |
| Memory parts (`'graph build'`, `'graph'`) | `cellucid/assets/js/data/memory-budget.js` | `cellucid/tests/connectivity-extraction-contract.test.mjs` |

## Related pages

- {doc}`07_rendering_pipeline_webgl_and_performance_notes`: the edge layer, K
  records on the GPU, the Feistel order in detail.
- {doc}`09_data_loading_pipeline_and_caching`: the long form of graph reads
  at a share.
- {doc}`20_read_path_architecture`: where payloads are made in the worker.
- {doc}`19_subsampling_specification`: which cells a share keeps.
- {doc}`23_decisions_partial_reading`: the share-cost rule and block order.
- {doc}`22_decisions_link_model_and_request_planning`: merge gaps, request
  splitting, S3.
- {doc}`27_decisions_memory_model`: working-set parts and allowances.
- {doc}`25_decisions_exports`: the other long reads of many rows.
- {doc}`01_codebase_map_and_entry_points`.
