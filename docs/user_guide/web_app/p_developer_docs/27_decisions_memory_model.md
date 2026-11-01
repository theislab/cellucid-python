# Decision record: memory model

This page answers: **"How much memory may the web app use, how is what it
holds counted, what happens when memory is short, and why is it built that
way?"** The full description of the model, its shares and its rules is
{ref}`memory-budgets` in {doc}`06_state_datastate_and_events`; this page
keeps the decisions behind it, the measurements, what was replaced, and the
limits still open. Module paths are under `cellucid/assets/js/` unless they
start with `cellucid`. Update this page in the same change as the code it
describes.

## At a glance

**Audience**
- Developers and coding agents adding a cache, a worker, a per-cell array,
  or anything that holds memory for longer than one call.

**Time**
- 10 minutes.

**Prerequisites**
- {ref}`memory-budgets` (the model as it stands).

---

## What it does now

1. **The device** is estimated once, on the page
   (`data/memory-budget.js` `estimatedDeviceMemory`): `navigator.deviceMemory`,
   raised by the JavaScript heap limit × 4 where a report of 8 GiB may be a
   cap, and 8 GiB where the browser reports neither (Firefox, Safari). The
   page hands it to every worker it starts (`connectMemoryWorker`).
2. **The allowance** is half the device, divided among the Cellucid pages of
   one origin (`data/memory-census.js`: Web Locks tell who is open,
   a `BroadcastChannel` what each holds): each page may use what the others
   do not hold, and at least an equal share.
3. **The working set** takes its part first: `CELL_WORKING_SET_BYTES`
   (77 bytes) per open cell, plus named parts (`setWorkingSet`): `'graph'`
   while edges are shown, a worker's `'graph build'`, the page's
   `'cell identity'`, each reader's `'cell rows N'` (an open subsample's
   rows, runs, inverse and ranks), and `'unaccounted'` where the browser
   measures the page (cross-origin isolated pages only, which Cellucid's
   hosts are not).
4. **The caches** divide what is left by `RESIDENCY_SHARES` (summing to
   exactly one). Every cache registers as a holder (`registerResidentHolder`:
   `bytes()` and `trim(budget)`), in whichever thread it lives; workers
   report to the page. When the allowance or the working set changes, every
   share over budget is trimmed (`settleMemory`), cheapest to rebuild first
   (`RELEASE_ORDER`). Decoded chunks of every open store share one budget
   (`formats/io/decoded-budget.js`).
5. **Short of memory.** A refused allocation smaller than the allowance
   (`allocateTypedArray` → `recordAllocationFailure`, any thread) halves
   the allowance everywhere; it doubles back a step at a time on evidence
   (a granted allocation at least as large while the app holds the lowered
   allowance, or the refused size granted when another dataset opens). A
   tab the browser ended mid-opening is found on the next load
   (`app/page-memory.js`, a session-storage mark) and halves it too.
6. **Hidden, frozen, idle.** A hidden page lets every cache but the gene
   index and undo steps go; frozen or in the back/forward cache, the gene
   index too (`app/page-memory.js`). Workers of every pool end after
   `WORKER_IDLE_MS` (10 s) idle (`workers/idle-retirement.js`).
7. **Budgets never refuse.** Each is a residency policy: what is evicted,
   spilled to disk or read again, never what can be opened or shown. No
   budget has a floor or ceiling of its own.

The Python and R servers keep a 64th of the machine's memory for what a
served store makes (`cache_bytes`, `cellucid-python/src/cellucid/_served_store.py`,
`cellucid-r/R/viewer-api.R`).

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-20 | **No fabricated limits.** A size or count bound exists only for a real platform limit (one `ArrayBuffer`'s length, a WebAssembly memory's 2^32 bytes, GPU texture size, a protocol field) or to bound untrusted input against what it declares. Every named limit is classified with its reason in `cellucid/tests/size-limit-allowlist-contract.test.mjs`. | An invented "512 MiB browser limit" refused real files and had spread to about ten places. | Refusing, truncating or hiding data because of its size; memory budgets that decide what opens. |
| 2026-10-26 | 77 bytes a cell in the working set, measured, not modelled. | Measured steady at 100M synthetic 3D cells with levels of detail (renderer and GPU processes, Chromium, M1 Pro): 64.6 bytes for the points, plus a field's codes, one continuous field and its quantile. An earlier 54 left out the spatial tree, draw order and GPU overhead, and fitted 148M cells to a 16 GB machine where 125M already swapped. | A per-cell figure nobody measured. |
| 2026-10-28 | **Count bytes, not the heap.** One ledger of holders in every thread replaces heap-percentage pressure. | `performance.memory.usedJSHeapSize` read 10,000,000 before and after allocating 512 MiB of typed arrays (quantised, rarely refreshed outside cross-origin isolation), counts no `ArrayBuffer`, and exists only on Chromium pages: pressure never arrived. | Thresholds on heap percentage; clock-driven cleanups. |
| 2026-10-28 | **Allowance = half the device; the working set first, the caches divide what it leaves.** | Before, the auto-fit and the caches were each sized against the full allowance: at the auto-fit limit they reached 100 % of the device (Chromium, 16 GiB Mac: 111.6M cells, 8 + 8 GiB). | Two budgets spending the same bytes. |
| 2026-10-28 | **The page estimates the device and hands it to workers.** | A worker has no heap limit to read: a Chromium that caps `deviceMemory` at 8 gave the worker 8 GiB and the page 14. | Per-thread estimates that disagree. |
| 2026-10-28 | **Pages of one origin share the allowance**, liveness by Web Locks. | Two tabs each assumed the whole allowance. A lock is released by the browser when a page is gone (closed, crashed, discarded), with no timer. | Heartbeats and timeouts. |
| 2026-10-28 | **A refusal halves the allowance in every thread; evidence raises it.** | The browser's refusal is the only memory signal every engine gives; most per-cell arrays are allocated in the reader worker, whose refusals the page did not hear before. | A fixed floor or ceiling on any budget. |
| 2026-10-28 | **One idle rule (10 s) for decode, analysis, plot raster and figure-export workers.** | A worker is an engine of its own (tens of MB; the figure-export worker holds WebGL contexts). | Pools that keep workers for the page's life. |
| 2026-10-30 | **Speed and less download before memory, while memory is counted, bounded and explained.** Unexplained memory is a defect. | The product's priority order. Memory above every cell's is acceptable when counted and explained. | Trading time or re-downloads for memory (for example, releasing a subsample's ranks after growth). |
| 2026-10-30 | **An index pass under way counts as held** (`building()` on index access and segments, in the reader's `_heldBytes`), so genes kept for a cell-subset run are never held beside a building gene index past the share. | At 50M, a kept cell subset beside a 50M gene-index pass took the reader worker to 8.7 GB and failed (`Array buffer allocation failed`). | Keeping a run's genes against an unbuilt index's reservation. |
| 2026-10-30 | **Analysis workers' marker contexts are counted at what they hold** (`markerContextBytes`: the codes and the 2-4 B group index each worker builds from them, `heldBytes`). | They were counted at the codes alone (1 B a cell a worker). | Uncounted per-worker indexes. |
| 2026-10-30 | **What a share holds is counted where it is held.** An open subsample's rows and ranks (`'cell rows N'`), the page's cell identity, the view's row pointers (`heldBuffers`, reported under `columnIndex`). | A 50M reader worker held far more than it reported. | Per-cell arrays outside the ledger. |
| 2026-10-31 | **The page's inactive field data is counted from the cached fields only** (`state/managers/field/residency.js`: each inventory indexed by cache key once per fields array and length, then a memory report looks up the cached keys, each view's active field and the user-defined aliases). | The count is asked on every memory report (each worker report, each working-set change); it walked every field of every view, 30,000 genes, each time: 0.5-1.4 s of the page's main thread per markers run at 10M cells (traced). The same answer as the walk, bit for bit (`residency-cached-keys-contract`). | Walking the inventories per report. |

## Measurements that justify it

Apple M1 Pro, 16 GB, Chromium on Metal unless stated.

- **The model replacing heap pressure** (1M cells, 2026-10-28, before →
  after): analysis workers 15 s after a differential expression 11 → 2; a
  tab hidden 15 s holds 191 → 46 MB of array buffers and 1,135 → 921 MB of
  renderer and GPU memory; two tabs' allowance 8,192 MB each → 8,035 and
  8,079 MB; a local dataset left for a sample keeps its reader worker →
  closes it.
- **Differential expression at 10M** (2026-10-28): the tab crashed (renderer
  +2.5 to 3 GB within 2 s); the cause was the reader worker's heap, 4.0 GB of
  about a million one-cell reads in flight, not the analysis workers. With
  every read priced at least one request of the link: Wilcoxon 2.8 s,
  renderer peak 4.5 → 2.4 GB, reader heap 4.0 GB → 4 MB.
- **Neighbour graphs at scale** (2026-10-30): forward adjacency (8 bytes an
  edge + 4 bytes a cell held) replaced edge records. 50M cells, 349.7M edges
  build and show in 40.0 s (refused before: "Array buffer allocation failed"
  after a 10.4 GB build), `graph` working set 2.9 GB; 10M edges reader worker
  peak 6.6 → 3.0 GB, page 2.5 → 1.2 GB. See {doc}`24_decisions_graph_edges_at_scale`.
- **Held against counted** (2026-10-30). Node, reader worker handler, live
  `ArrayBuffer`s after collection against `memoryReport()`:

  | Case | Before | After |
  |---|---|---|
  | 50M Zarr v3, every cell, first k-NN rows | +1,087 MB (two 537 MB tables), 0 counted | +214 MB (array inverse), counted |
  | 10M Zarr v3, 5M uniform share, after opening | 201 MB held, 20 MB counted | 140 MB counted (one `CellRows`) |
  | 10M Zarr v3, every cell, edges + 21 genes, released | 488 MB held, 487 counted | unchanged |

  Browser heap snapshot, 50M Zarr v3 sharded store picked, after edges: the
  reader worker's 2,357 MB were all live `ArrayBuffer`s (no allocator
  slack); uncounted were the open view's row pointers (400 + 400 + 200 MB)
  and a kept index copy's pointers (42 MB), now counted.
- **Page side at 10M** (2026-10-30): a gene crosses from the worker in its
  sparse form when smaller, and inactive genes are kept as that form
  alone; a gene shown on the page 66 to 70 → 6 to 10 ms of main thread
  (Node).

## Tried and reverted

| What | When | Result | Would need |
|---|---|---|---|
| Heap-percentage monitor (75 / 85 / 95 %), the dimension manager's own heap check, a 5-minute periodic cleanup | replaced 2026-10-28 | Blind on every browser (above); the cleanup cleared caches on a clock whatever the memory. | Nothing: the ledger replaces them. |
| Releasing a subsample's rank arrays after the last growth step | 2026-10-30 | Built, then removed: the next Apply had to rank again, seconds at 50M. The ranks stay, counted in `'cell rows N'`. | A case where those bytes are unexplained, which they are not. |
| The inverse of an open subsample's rows always as a hash table | replaced 2026-10-30 | 1.07 GB at 50M cells, uncounted. Now the smaller of the table and a 4-byte-a-row array (200 MB), counted. | — |

## Known limits and open items

Real platform limits the model designs around (details in
{doc}`07_rendering_pipeline_webgl_and_performance_notes`): one `ArrayBuffer`
(Chromium 2,145,386,496 bytes, WebKit 2^32, Firefox 2^34), so larger
per-cell arrays are made in a WebAssembly memory (`utils/large-arrays.js`);
a WebAssembly memory's 2^32 bytes; Firefox and Safari report no device
memory (taken as 8 GiB; a tab ended mid-opening halves it).

Open:

- **50M reader worker, held against counted.** After the pointer fix, a
  168 MB read held in a promise and about 250 MB of arrays with no owner in
  the snapshot were still unattributed at 50M (not reproduced at 10M).
- **Not counted** (stated in {ref}`memory-budgets`): GPU objects beyond the
  working set's measured share (kept views' colour stores, highlights, edge
  geometry), sharded Zarr shard indexes, and the read activity log's
  operations under way.
- **A subset read's short-lived buffers** in the reader (about 1 GB of
  transient `ArrayBuffer`s at 10M during a differential expression) are not
  in the ledger.
- **Firefox during differential expression**: about 870 MB more than
  counted, unattributed (analysis buffers held account for about 55 MB at
  2.28M cells; attributing the rest needs a per-thread measurement).

## Known issues (pre-commit)

Defects, regressions and unverified behaviour known when this version was
committed (2026-10-30), each with where it is recorded above. An item is
removed in the change that fixes or verifies it.

- **Firefox during differential expression:** about 870 MB more than counted, unattributed.
- **50M reader worker:** about 418 MB held but not attributed to a counted holder (a 168 MB read held in a promise and about 250 MB of arrays with no owner in the snapshot).
- **DE at 10M:** about 1 GB of short-lived buffers in the reader, not in the ledger.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| Device, allowance, working set, shares, ledger, refusals, report | `data/memory-budget.js` | `cellucid/tests/memory-budget-contract.test.mjs`, `cellucid/tests/memory-pressure-release-contract.test.mjs` |
| Pages of one origin | `data/memory-census.js` | `cellucid/tests/memory-census-contract.test.mjs` |
| Hidden, frozen, ended mid-opening | `app/page-memory.js` | `cellucid/tests/interrupted-opening-contract.test.mjs`, `memory-pressure-release-contract` |
| Worker lifetimes | `workers/idle-retirement.js` | `cellucid/tests/worker-idle-retirement-contract.test.mjs` |
| Decoded chunks' shared budget | `formats/io/decoded-budget.js`, `data/anndata/decoded-chunks.js` | `cellucid/tests/decoded-chunks.test.mjs` |
| Allocations that report refusals | `data/sparse-column-index.js` (`allocateTypedArray`), `utils/large-arrays.js` | `cellucid/tests/per-cell-allocation-contract.test.mjs` |
| Open rows and ranks counted | `data/anndata/reader.js`, `data/anndata/views/cell-subset.js`, `data/subsample.js` | `cellucid/tests/cell-rows-lookups-contract.test.mjs` |
| Analysis memory | `app/analysis/shared/gene-stream-runner.js`, `app/analysis/compute/worker-pool.js`, `app/analysis/data/data-layer.js` | `cellucid/tests/analysis-gene-reads-bounded-contract.test.mjs` |
| Auto-fit against the allowance | `data/anndata/subsample-plan.js` | `cellucid/tests/scale-auto-fit-budget.test.mjs` |
| Page holders and per-cell array ownership | `app/state/core/data-state.js` (inactive fields, unshown layouts, cell identity), `app/state/managers/field/residency.js`, `app/state/managers/per-cell-array-ownership.js` | `cellucid/tests/view-memory-ownership-contract.test.mjs`, `cellucid/tests/residency-cached-keys-contract.test.mjs`, `cellucid/tests/browser/kept-views-memory.spec.mjs` |
| Named limits and their reasons | (every module) | `cellucid/tests/size-limit-allowlist-contract.test.mjs` |
| Server caches | `cellucid-python/src/cellucid/_served_store.py`, `cellucid-r/R/viewer-api.R` | — |

## Related pages

- {ref}`memory-budgets` in {doc}`06_state_datastate_and_events`: the model,
  its shares and a table of budgets per device.
- {doc}`09_data_loading_pipeline_and_caching`: caches, kept data and
  browser memory.
- {doc}`20_read_path_architecture`: memory as a cross-cutting concern of the
  read path.
- {doc}`07_rendering_pipeline_webgl_and_performance_notes`: GPU memory and
  `ArrayBuffer` limits.
- {doc}`21_decisions_decoders_and_codecs`: decode workers and their memory.
- {doc}`23_decisions_partial_reading`: what a share costs, memory included.
- {doc}`24_decisions_graph_edges_at_scale`: the graph's working set.
- {doc}`26_decisions_served_store_and_widgets`: the servers' caches.
