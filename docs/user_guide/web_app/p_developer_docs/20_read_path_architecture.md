# Read path architecture

This page answers: **"Which module owns each step between a dataset's bytes
and its pixels, and what must stay true there?"** It is a map, short on
purpose: each step names its module and entry points (all under
`cellucid/assets/js/`) and the invariants that module's header and tests
hold. The long-form description of sources, caches and failures is
{doc}`09_data_loading_pipeline_and_caching`; which cells a subsample keeps is
{doc}`19_subsampling_specification`; drawing is
{doc}`07_rendering_pipeline_webgl_and_performance_notes`. In the code tree,
`assets/js/ARCHITECTURE.md` has the layers, the glossary and a "where to
change X" index, and each folder has a short `README.md`.

## At a glance

**Audience**
- Developers and coding agents changing how data is read, priced, kept or
  drawn.

**Time**
- 15 minutes.

**Prerequisites**
- {doc}`01_codebase_map_and_entry_points`.

---

## The path

```
 page thread                                   reader worker (one per open dataset)
 ───────────                                   ─────────────────────────────────────
 1 source        data/*-source.js,             
                 data/sources/ (provider,      
                 identity, remote-file)        
        │  open(source, {cellSubset})           
        ▼                                      
 6 client        data/anndata-reader-client.js ──postMessage──▶ data/anndata-reader-worker.js
                                                                 │
                                               2 bytes           formats/io/: byte-source → byte-reader
                                                                 │  (+ kept-parts, range-cache, link, hosts)
                                               3 containers      formats/{hdf5,zarr,zip,rds,text}
                                                                 │  → store model (formats/store.js)
                                               4 AnnData         data/anndata/open.js → views/ | transcode/
                                                                 │  → AnnDataReader (reader.js)
                                               5 share planner   subsample-plan.js → views/cell-subset.js
                                                                 │  gene sample, column indexes, graphs
        ◀──── payloads (typed arrays, transferred) ──────────────┘
 7 page state    app/state/ (DataState, field managers, gene-sample.js),
                 app/open-cells-resize.js
 8 renderer      rendering/ (alpha layers, read-cell layer, HighPerfRenderer, LOD)
```

Cross-cutting, on both threads: memory (`data/memory-budget.js`),
cancellation (`utils/abort.js`), progress and logs (`formats/io/activity.js`,
`data/read-activity-log.js`, `data/transfer-log.js`,
`data/opening-progress.js`).

---

## 1. Sources: what is opened, and which version

| Owns | Module | Entry points |
| --- | --- | --- |
| What a file or folder is | `data/input-formats.js` | the format registry (mirrored byte for byte as `input_formats.json` in Python and R, digest-pinned in all three suites) |
| One class per way a dataset arrives | `data/local-user-source.js`, `link-source.js`, `local-demo-source.js`, `remote-source.js`, `viewer-host-source.js`; all read through `anndata-file-source.js` | `DataSourceManager` (`data/data-source-manager.js`) |
| Link providers | `data/sources/providers/<host>.js`, `data/sources/registry.js` | `getProviderRegistry`, `formatLocator` |
| A remote file described without downloading it | `data/sources/remote-file.js` | `describeRemoteFile`, `downloadRemoteFile` (hosts that serve no ranges) |
| The identity and version of a file | `data/sources/identity.js` | `fileIdentityOf`, `persistsAcrossVisits` |
| A link's public form | `data/sources/link-credentials.js` | `publicLink`, `linkIdentity` |

Invariants: only a file whose identity names a version is kept across visits;
credentials never enter an identity, a log, a session or a URL the page
shows; feature state never goes into URL parameters.

## 2. Bytes: requests, the link, kept bytes

| Owns | Module | Entry points |
| --- | --- | --- |
| Where bytes come from | `formats/io/byte-source.js` | `blobByteSource`, `httpByteSource`, `bytesByteSource`; `fetchOnce` (the one place a request to a web host is sent); priors `REMOTE_POLICY`, `LOCAL_POLICY` |
| Turning reads into requests | `formats/io/byte-reader.js` | `createByteReader`; `planMergeGap`, `requestPieceBytes` |
| The measured link | `formats/io/link.js` | `LinkMeter`: `latencyMs` (shortest round trip), `loadedLatencyMs` and `loadedRoundTripMs` (round trip under load), `loadedRequestMs`, `bytesPerMs`, concurrency |
| Each web host | `formats/io/hosts.js` | `createHostBook`, `realmHosts`: answer times, deliveries, rate-limit budgets |
| What a host delivered on earlier visits | `data/host-link-memory.js` | seeds a link before it is measured |
| Bytes kept across visits | `formats/io/kept-parts.js` | `withKeptParts`, journals checked by `page-checks.js` |
| Bytes other readers of the page fetched | `formats/io/range-cache.js` (worker), `data/byte-range-cache.js` (page) | `withRangeCache` |

Order of a read: the byte reader's own memory, then kept parts (this
machine's disk), then the page's range cache, then the host.

Invariants:

- A gap between two ranges is fetched only when it saves time on the link as
  measured now; plans never merge across bytes already kept or held, and a read
  joins a request in flight that carries its bytes.
- Request **shaping** uses the shortest round trip; **prices** of reads use
  the round trip under load (`loadedRoundTripMs`) and no less than the host's
  measured cost per request.
- Until the link is measured, the source's prior stands (a local file is a
  disk; a web host is read ahead only once measured).

## 3. Containers: the store model

| Owns | Module | Entry points |
| --- | --- | --- |
| The model every format is opened as | `formats/store.js` | `Group.load`, `get`; `ArrayNode.read`, `storedBytes`, `readRows`, `chunkIndexReads` |
| HDF5 (lazy tree, chunk indexes, read plans) | `formats/hdf5/` | `openHdf5` |
| Zarr v2/v3 (folders, zip, web hosts) | `formats/zarr/` | `openZarr`, `httpStore`, `directoryStore`, `zipStore` |
| ZIP, R serialization, MatrixMarket | `formats/zip/`, `formats/rds/`, `formats/text/` | `openZipArchive`, `readRFile`, `readMatrixMarket` |
| Decoding off the reader thread | `formats/io/decode-pool.js`, `decoder-choice.js`, `wasm-codecs.js`, `inflate.js` | `decodeChunkBytes`, `warmDecodePool` |

Invariants: a value returned by `read()` belongs to the caller; an unloaded
group throws `StoreNotLoadedError` rather than answering empty;
`storedBytes(selection, {toFetch: true})` counts only bytes a read would
fetch now, and every time estimate prices that.

## 4. AnnData: one reader for every format

| Owns | Module | Entry points |
| --- | --- | --- |
| From a source to a store root | `data/anndata/open.js` | `openDatasetStore` |
| Formats presented as AnnData in place | `data/anndata/views/` (`loom`, `tenx-h5`, `tenx-outs`, `mudata`, `modalities`) | `views/index.js` by input-format id |
| Streams converted once into a kept Zarr store | `data/anndata/transcode/` (`r-object`, `matrix-market`, `tables`; storage in `space.js`) | `transcode/index.js`, `session.js` |
| The reader | `data/anndata/reader.js` | `AnnDataReader` (open, fields, genes, embeddings, graphs, category codings) |
| Gene columns and cell rows by layout, and their prices | `data/anndata/matrix-access.js` | `geneLayoutOf`, `requestPriceOf`, each layout's `readCost`, `columnIndexResidentBytes` |

Invariant: a format-specific path is a view or a transcoder, never a branch
in `reader.js`.

## 5. The share planner: which cells, which genes, at what cost

| Owns | Module | Entry points |
| --- | --- | --- |
| The request an opening sends | `data/cell-subset-policy.js` | `requireCellSubsetRequest`, `setTabCellSubsetRequest`, `takeCellSubsetRequest` (modes `auto`, `all`, `subsample`, `kept`, `prefix`, `chunks`; listed in the `subsample-plan.js` header) |
| Which cells a subsample keeps | `data/subsample.js` | the specification ({doc}`19_subsampling_specification`), `priorityOrder`, `createSubsampleCache` |
| Resolving the request against the store, before reading | `data/anndata/subsample-plan.js` | `resolveCellSubset`, `autoFitCount` (fits the device's memory; an auto-fit opens `blocks`, or `uniform` under the file's seed for a priority-ordered file), `defaultBlockSize` and `blockRung` (the `blocks` size: the storage's block doubled up to what the link reads in a round trip, or the rung of a kept gene index) |
| The gene sample (open cells a cell-major gene is read for) | `data/anndata/subsample-plan.js` | `geneSampleOf`, `geneTimeFit` (the first gene's time over the measured link, `GENE_READ_TARGET_MS`), `geneIndexFit` (its index's memory) |
| The subsample as a smaller AnnData | `data/anndata/views/cell-subset.js` | `cellSubsetView`, `gatheredArray` (kept rows as boxes in source order, in waves), `KeptPositions`, `CellRows` |
| Rows read in store order | `data/anndata/values.js` | `readRowsInStoreOrder`, `readInSourceOrder` |
| Gene indexes of cell-major matrices | `data/sparse-column-index.js`; `data/anndata/prefix-column-index.js` (segments over priority-ordered rows, doubling, carried to the next opening of nested cells); `data/anndata/kept-column-index.js` (kept across visits) | `AnnDataReader` `_segmentedColumns` (rows spread through the store form one segment in store order, or every stored row's index restricted to them: `_restrictsWhole`, `restrictedAccess`, `PassRates`) |
| Neighbour graphs at a share | `data/anndata/graph-lines.js` | `inducedGraphLines`, `readInducedGraph` |
| Progressive first view and growth | `data/anndata/progressive-open.js`, `data/opening-cost.js` | grown in place by `app/open-cells-resize.js` |

Invariants:

- **The share-cost rule.** A share never costs more than a larger share or
  than every cell in what a person meets: time, download, memory, rendering
  and analysis work (the request count is a means, not a cost: decision
  2026-10-30, {doc}`23_decisions_partial_reading`); the ideal is about s/n
  of every cell's cost. A read is priced before it is made, over the link as
  measured, and the cheaper of equivalent routes is taken (a cell subset by
  rows or by genes, an induced graph by rows or streamed whole). In code the
  "never more than every cell" half is the `wholeSampleCost` ceiling of
  `subsample-plan.js`; both halves are asserted by
  `tests/subsample-share-cost-contract.test.mjs` (stored chunks and requests)
  and `tests/obs-strings-share-contract.test.mjs`.
- **Exact cells.** A subsample is exactly the cells its request names, never
  raised; values read through the view are bit-identical to the same rows of
  a whole read; no chunk holding no kept row is fetched.
- **Counted memory.** What grows with cells or values is counted in
  `data/memory-budget.js` (auto-fit, `columnIndex` share, decoded chunks, held
  bytes); a read holds at most two slabs of decoded boxes, never a
  full-length array.
- What a gene costs never sizes the map: a cell-major gene is read for the
  gene sample, and the sample grows on the page (step 7).

## 6. Worker and client

| Owns | Module | Entry points |
| --- | --- | --- |
| The worker's API over one reader | `data/anndata-reader-worker.js` | one handler per client call |
| The page's end | `data/anndata-reader-client.js` | `AnnDataReaderClient`: `open`, `reopenCells`, `getObsField`, `getMatrixColumn`, `getMatrixColumns`, `getMatrixCellSubset`, `readEmbedding`, `geneSampleOf`, `growGeneSample`, `readMatrixRows`, `linkMeasurement` |
| What is read first, and never twice at once | `data/fetch-scheduler.js` | |
| Genes straight to an analysis worker | `data/gene-delivery.js` | `deliverMatrixColumnsSparse`, `deliverMatrixCellSubset` (client) |
| Payloads made in the worker | `data/embedding-payload.js`, `data/connectivity-edges.js`, `data/discovered-categories.js` | coordinates, forward adjacency, codes ready to draw |
| The adapter the app talks to | `data/base-anndata-adapter.js`, `data/data-loaders.js` | |

Invariants: the page never blocks on a file; typed arrays cross by transfer;
failures cross in one wire format (`utils/error-wire.js`); cancellation is
`utils/abort.js` on both sides.

## 7. Page state

| Owns | Module | Entry points |
| --- | --- | --- |
| The state coordinator | `app/state/core/data-state.js` | `DataState` (mixins from `app/state/managers/`) |
| Field loading and residency | `app/state/managers/field/loading.js`, `residency.js` | |
| The gene sample on the page | `app/state/managers/field/gene-sample.js` | `getGeneSample`, `readGenesForMoreCells` (the legend's Read for all cells), `topUpGeneSample` (what is left of the gene time after the first gene reads more cells) |
| Which cells are open, across boundaries | `data/cell-identity.js` | `createCellIdentity` (`toOriginal`, `fromOriginal`, `openValues`) |
| Changing the open cells in place | `app/open-cells-resize.js`, `data/open-cells-work.js` | |

Invariants: inside the app per-cell arrays are in open-cell positions;
everything that leaves (sessions, exports, viewer protocol, Python and R)
names cells by original position, converted only through the cell identity.
A cell not read for a gene holds NaN (no value), never a zero.

## 8. Renderer

| Owns | Module | Entry points |
| --- | --- | --- |
| What a colour's alpha byte means | `rendering/alpha-visibility.js` | `BENEATH_ALPHA_BYTE` (254, no-value grey beneath), `UNREAD_ALPHA_BYTE` (253, not read), 255 coloured |
| How unread, zero and read cells look | `rendering/read-cell-layer.js` | `zeroCellColor`, `unreadCellColor`, `readCellSizeMultiplier` |
| Points, LOD, culling | `rendering/high-perf-renderer.js`, `rendering/high-perf/` | `HighPerfRenderer`; LOD levels are the first cells of the priority order (`data/cell-priority.js`) |

Invariants: every surface that draws cells (viewer, figure export, smoke)
takes its layer colours and alpha meanings from these two modules; unread
cells are drawn first, so every read cell shows above them.

---

## Why it is this way: decision records

Each step's dated decisions, the measurements per engine that justify them,
what was tried and reverted, and the open items:
{doc}`21_decisions_decoders_and_codecs` (step 3 decoding),
{doc}`22_decisions_link_model_and_request_planning` (step 2),
{doc}`23_decisions_partial_reading` (step 5),
{doc}`24_decisions_graph_edges_at_scale` (graphs in steps 5 and 8),
{doc}`25_decisions_exports`,
{doc}`26_decisions_served_store_and_widgets` (served sources in steps 1 and 2),
{doc}`27_decisions_memory_model` (counted memory, every step),
{doc}`28_decisions_rendering_layers` (step 8).

---

## Where the read path is tested

| Step | Tests (in `cellucid/tests/`) |
| --- | --- |
| Bytes | `formats-io-contract`, `host-*`, `http-request-resilience-contract`, `kept-parts-*`, `formats-decode-pool`, `codec-layer-contract`, `wasm-codecs-contract` |
| Containers | `hdf5-*`, `zarr-*`, `rds-*`, `text-matrix-*`, the corpus conformance tests |
| AnnData and share planner | `anndata-*`, `subsample-*`, `gene-sample-contract`, `kept-column-index-contract`, `graph-lines-contract`, `edges-at-a-share-contract`, `obs-strings-share-contract`, `subsample-share-cost-contract` |
| Page and renderer | `high-perf-*`, `alpha-*`, browser specs `gene-sample`, `renderer-smoke`, the `*-parity` specs |

`cellucid/scripts/efficiency-floor/` measures every source × format × share
against the floor the link allows.

---

Next: {doc}`01_codebase_map_and_entry_points` for the rest of the tree.
