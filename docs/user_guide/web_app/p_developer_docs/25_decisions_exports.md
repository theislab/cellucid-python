# Decision record: cell exports

This page answers: **"How does the web app write the highlighted cells or a
field out as a CSV or an `.h5ad`, from a file of millions of cells over a
slow link, and why is it built that way?"** It records the current design,
the decisions behind it with their dates, the measurements that justify
them, what was tried and reverted, and what is still open. Update this page
in the same change as the code it describes. Figure export (PNG and SVG of
the views) is a separate system: see {doc}`12_figure_export_architecture`.

## At a glance

**Audience**
- Developers changing the export buttons, the export
  writers, or the reader's many-row reads.

**Time**
- 10 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (reader worker, client, adapter) and
  {doc}`22_decisions_link_model_and_request_planning` (how a read is priced
  and split into requests).

---

## What it does now

Paths are under `cellucid/assets/js/`. The user-facing side is
{doc}`../f_highlighting_selection/03_highlight_ui`.

1. **Three exports** (`app/exports/cell-exports.js`, buttons in
   `app/ui/modules/cell-exports-ui.js`):
   - a page's highlight groups as an obs-column CSV (`highlightsCsv`): one
     `cellucid_highlight__<group id>` column per group, one row per open
     cell;
   - one field as an obs-column CSV (`fieldCsv`);
   - the highlighted cells as an `.h5ad` subset (`highlightsH5ad`): obs
     names, every obs field, each group as a boolean column, the genes, the
     expression rows of the matrix the data mapping reads (CSR, float32),
     and the embeddings the mapping shows.

   An obs-column CSV declares itself with a header row, a kind row
   (`highlight`, `categorical`, `continuous`) and a categories row (each
   `categorical` column's categories in the field's order, a JSON array of
   strings in one CSV cell; empty for the other kinds). The readers build a
   categorical's categories from that row, never from the rows' order of
   first appearance, which lost a trajectory's order (Pancreas `clusters`)
   and every category no cell held. `cellucid/tests/fixtures/obs-columns.csv`
   pins that format for the web app, Python and R. Python's
   `cellucid.pull()` and R's `cellucid_pull()` read all three files back.
2. **Cells are listed in ascending original position**, whatever order the
   file stores them in or the app opened them in. The order comes from a
   radix sort of 32-bit original positions (`app/exports/cell-exports.js`
   `inOriginalOrder` →
   `utils/radix-order.js` `orderByKey`): two linear passes and one
   `Uint32Array`, never a `Map` of cells or a comparator sort.
3. **The open cells stay fixed while an export runs** (`useOpenCells`,
   `data/open-cells-work.js`). The export starts once no growth step or
   Apply is changing the open cells, and none commits until the file is
   written, so every row it reads belongs to the cells it started on.
4. **A CSV is written a batch at a time.** Each batch is 65,536 rows
   (`CSV_BATCH`), with a yield between batches, progress on the card, and
   the signal checked. Obs names are read by range when the cells are
   consecutive; for the listed cells only when they ascend
   (`getObsNames({cells})` → `AnnDataReader.readObsNames({cells})`: their
   rows of the index, in stored order); and, when more than a batch is
   listed out of the stored order (every cell of a priority-ordered file,
   or of a share opened in a subsample's order, by original position), once,
   as the range of the open cells, the batches named from it: batch by batch
   each batch's cells lie throughout the file and every name chunk was read
   again (every cell of a 10M-cell priority-ordered file: 15 GB in 400 s).
   The names of the next two batches are asked for while one is made and
   written (`CSV_READ_AHEAD`): read one after another, each waiting on the
   link for its references and then its strings, the link idled between
   them. Each batch goes to the browser's
   origin-private storage as it is made (`createSpool`: a file for each
   segment of about 32 MB written), or to in-memory Blobs where that storage is unavailable. When
   the site's storage is full (a private window's is small; gene indexes
   kept for a revisit share it), the files of this page's earlier exports
   are removed and the write tried again, and a write that still fails, and
   every one after it, is kept in memory: a full storage never fails an
   export, and no byte written is lost. After the first batch the card
   states the expected file size.
5. **An `.h5ad` export says what it will read before it starts.**
   `highlightsH5adCost` → adapter `getMatrixRowsCostOf` →
   `AnnDataReader.matrixRowsCost` prices the expression rows from the file's
   layout alone, without making a request. It uses the same model as the gene
   sample and *Genes read for* (`firstGeneReadCost`), and returns bytes,
   bytes already in this browser, requests, the time over the link as
   measured (`ms`, null until the link is measured), `linkMs` and
   `decodedBytes`. The card shows one compact line, with the parts under
   *Details*.
6. **The expression rows are one streamed read.** `readMatrixRowsOf` →
   `AnnDataReader.readMatrixRows` reads the rows' pointers once and plans the
   whole read. It then reads parts of about `SLAB_BYTES` of entries in the
   cells' order, asking for each part before the previous one is handed out.
   Each part goes to `onPart`, which appends it to two spools (`data`,
   `indices`). A spool on a full storage frees this page's earlier exports'
   files and then keeps the rest in memory (item 4): no part fails to be
   written for want of room.
   The export holds about two parts at a time, not every cell's rows.
   `getMatrixRows` is the same read with its parts gathered, for callers that
   want all rows at once.
7. **Fields and embeddings.** Each obs field is loaded, written, and then
   released again under the residency budget
   (`releaseInactiveFieldsOverBudget`). Embedding rows are read for the
   exported cells alone, in the export's order (`getObsmColumns({cells})` →
   `AnnDataReader.readObsmColumns({cells})`: their rows in the stored
   order, `gatheredArray`), never every open row between them.
   `formats/hdf5/writer.js` lays out the file as a Blob of parts; a spooled
   dataset stays a Blob and is never copied.
8. **The page does not keep bytes that only an export reads.** A read of
   kind `cells` is `ReadActivity.readOnce` (`formats/io/activity.js`).
   Requests that serve only such reads are not handed to the page's cache of
   fetched bytes (`formats/io/byte-reader.js` `keepInMemory`), where they
   would push out what other reads use. That cache finds ranges by binary
   search (`data/byte-range-cache.js`).
9. **Cancel and failures.** Each export runs under its own signal. *Cancel*
   on its card stops it, and nothing is saved. A read that fails after
   Cancel was pressed counts as canceled. A second click while an export is
   running gets a notice instead of being silently ignored. A failure whose
   message is empty is reported with its class and cause (`failureText`,
   `utils/errors.js`).

## Decisions and why

| Date | Decision | Reason | Rules out |
|---|---|---|---|
| 2026-10-28 | Order export rows by a radix sort of original positions | An export of more than 2^24 cells threw "Map maximum size exceeded"; at 10M cells the `Map` cost 672 MB | A `Map` per export; comparator sorts over cells |
| 2026-10-28 | Write CSVs a batch at a time to origin-private storage, reading names per batch | A CSV held every cell's name and built every row in one task | One string or one names array for the whole file |
| 2026-10-29 | One streamed read of every exported row (`readMatrixRows`) | Six separate 4,096-cell batch reads made six pointer plans and six gather plans, re-read chunk parts across batches, and merged less | Batch reads; the batch `getMatrixRowsOf` (removed) |
| 2026-10-29 | Price the read first, from the layout, and show it on the card | A person should see a long read's size and time before it runs; same model as *Read for all cells* | Starting a multi-minute read without saying so |
| 2026-10-29 | Cancel on the card; nothing is saved after Cancel; a second click is answered | Exports were long operations nobody could stop | Exports without a signal |
| 2026-10-29 | Every failure is shown with words (`failureText`) | Firefox showed an export failure as an empty-message `Error` | Notifications that show `error.message` alone |
| 2026-10-30 | Names of more than a batch listed out of the stored order are read once, as one range | Every cell of a 10M-cell priority-ordered file: each 65,536-cell batch read every name chunk again, 15 GB in 400 s (the committed tree read every name once) | Batch-by-batch listed reads of cells spread over the whole file (the 2026-10-28 rule "never every name at once" holds where the cells ascend) |
| 2026-10-30 | A CSV reads the names of two batches ahead of the one it makes | Every cell of a 10M-cell file took 11–25 s, against 7 s for the committed tree's one read of every name: each batch's names were asked for only once the one before was written | Batches read strictly one after another |
| 2026-10-30 | Embedding rows of an export are read for its listed cells alone | Runs of open rows merged across 4,096-row gaps (`rowRuns`, removed) covered every open row at a uniform share: one cell type of a 50% share of 10M cells read all 5M rows, 173 MB in 5.1 s, more time than every cell's export (0.8 s now) | Runs of open rows between exported cells |
| 2026-10-30 | A spool writes a file for each segment of about 32 MB written; on a full storage it frees this page's earlier exports' files, then keeps the rest in memory | An `.h5ad` export of every cell of a 10M-cell file failed with `QuotaExceededError` after a 148 MB CSV export in the same page (a spool's files stayed until the next page); a file for each 1 MB CSV batch cost a file's creation each: 153 files for every cell of a 10M-cell file, seconds of an export the committed tree made in memory in 5 s | An export failed by storage the browser could free, or by bytes a single file write had already swallowed; a file a write |
| 2026-10-30 | Rows of a gene-major (CSC) matrix are a cell subset of every gene (`readMatrixRows` → `getMatrixCellSubset`), not an index of every stored row | 2,000 cells of a 1% share of a 10M-cell priority-ordered file decoded every chunk of X (244 MB, 10 s, 1.9 GB peak); now each gene's prefix (0.8 s, 0.7 GB) | A transposed index of the whole matrix for an export |
| 2026-10-29 | Public S3 objects are read over seven addresses (42 connections, not 24) | S3 speaks HTTP/1.1, six connections per host; every connection was busy | (see {doc}`22_decisions_link_model_and_request_planning`) |
| 2026-10-30 | Bytes only an export reads skip the page's cache of fetched bytes | They were copied and posted to the page (812 MB for an HLCA export) only to evict what other reads use | Caching write-once bytes |
| by 2026-11-12 | The page's cache of fetched bytes finds a range by binary search over the ranges kept, which are ordered because none holds another (`data/byte-range-cache.js`) | An export of tens of thousands of scattered chunks kept each one in a scan of every range kept: 2.8 s of the page's main thread on a 21,285-cell export of HLCA core | A scan of every kept range per chunk |
| 2026-10-30 | A running read keeps its gap bytes only within the memory model's budget for fetched bytes, oldest first out | A 10M-cell export held 210 of 216 MB of requests until it ended | Unbounded retention by a long read |
| 2026-10-30 (owner) | Speed comes before memory, as long as memory is counted and explained | The stream holds two parts in flight for a busier link | Trading time for an unexplained memory saving |

## Measurements that justify it

**Browser, HLCA core over real S3, 21,285 scattered cells, `.h5ad`
(2026-10-29 before, 2026-10-30 after; s / MB fetched / requests / requests a
second / renderer peak).**

| Engine | Before | After |
|---|---|---|
| Chromium | 479 s / 808 MB / 24,462 / 51 / 2,767 MB | 468 s / 726 MB / 32,685 / 70 / 3,552 MB (end footprint 1,974 MB) |
| Firefox | 463 s / 964 MB / 20,999 / 45 / 1,807 MB | 514 s / 590 MB / 34,873 / 68 / 1,774 MB |

Other runs:
- Firefox, full HLCA, whole journey: the `.h5ad` export now finishes (790 s,
  812 MB, 52,085 requests). Before, it failed at 600 s.
- WebKit, a served local 10M-cell file with a 333,333-cell cluster: 6.9 s,
  763 MB, 1,805 requests.
- The cell-exports spec passes in all three engines.

**Why the time is what it is.** The time is set by request count ×
latency under load ÷ connections, not by stalls. Throughout the run all 42
connections were busy. This machine's path to S3 (us-west-2, from Europe)
completes about 70–88 small requests a second whatever the connection
count, and moves about 3.2 MB/s in bulk. So a plan's floor is the larger of
requests ÷ ~70 a second and bytes ÷ 3.2 MB/s. The plan above (about 35k
requests, 590 MB) has a floor of 498 s: Firefox ran at 1.03× of it and
Chromium at 1.0×.

**Node, the app's reader, real S3 (HLCA core, 21,285 cells).**

| Setup | Time / requests / MB |
|---|---|
| 4 hosts, batches (before) | 493 s / 28,373 / 596 |
| 7 hosts, batches | 383 s / 30,516 |
| 7 hosts, one stream | 271 s / 27,201 / 620 (1.18× its floor) |

**Simulated link, 10M-cell CSR file.** At 150 ms and 5 MB/s, 21,285 cells
took 65.8 s with batches and 57.6 s with the stream (1.20× → 1.08× of the
floor; 1,246 → 263 requests). 333,333 cells took 93.8 s → 62.6 s (1.75× →
1.09×). The stream's peak memory is two parts in flight: 333k cells went
from 310 to 718 MB.

**More than 2^24 cells (Node, 2026-10-28).** Ordering 10M cells in priority
order went from 1.4 s to 0.24 s, and the 672 MB `Map` is gone. 17M cells
went from a `RangeError` to an export that works.

**Deflate (Node 26, Apple M1 Pro).** The Zarr v2
writer's deflate is fflate's, synchronous, with the same bytes in every
engine and in Node. On chunks of 2^14 shuffled `int32` gene indices, fflate
runs at 79 MB/s, Node's zlib at 134 and `CompressionStream` at 93; on
float32 expression values all three run at 23 MB/s, and fflate's output is 5%
smaller. Deflate is the slow part of writing a store, about 25 MB/s of
float values on one core, which is why a copy's chunks are encoded on every
decode worker (`formats/zarr/encode-task.js`).

**Memory traced (Firefox, full HLCA journey).** The renderer footprint grows
from 2,202 MB to 4,438 MB during the `.h5ad` export. Of that, 233 MB is
counted (obs fields read for the export, within their budget). The rest is
allocation churn that the engine collects late: the export decodes 2.7 GB of
whole chunks to keep 0.2 GB of rows. With a collection every 50 ms the live
peak is 420 MB.

## Tried and reverted

| What | When | Why it was reverted | What would have to change |
|---|---|---|---|
| Decode workers decoding into reused per-task buffers and sending back only the part of each chunk an export keeps (HDF5 boxes, Zarr runs) | 2026-10-30 | No read shared a chunk with its neighboring runs, and none was cached, so chunks touched by two runs (in-flight splits, part boundaries) were fetched and decoded twice. HLCA core: Chromium timed out at 780 s with 1,351 MB (754 s / 970 MB before); Firefox 723 s / 801 MB against 751 s / 622 MB; peak renderer memory barely moved | One read per export part that lists all its row ranges (HDF5 and Zarr multi-range reads, the matrix access `rows()`), each chunk decoded once and its kept ranges cut in the worker: a reader redesign |
| Pricing a request at the link's shortest round trip | replaced 2026-10-30 | The export was priced at about half its real time (round trip 230 ms idle against 0.35–0.5 s under load) | Now priced at latency and request rate under load ({doc}`22_decisions_link_model_and_request_planning`) |

## Known limits and open items

- **Possible regression.** A Chromium run of the HLCA-core export later on
  2026-10-30 took 754 s / 970 MB / 28,775 requests, against 468 s that
  morning. It may be the day's link and byte-reader changes or a slower S3
  path; this has not been settled.
- **Request-bound plans.** At about 35k requests the plan is limited by
  requests, not bytes (the earlier plan made 21k requests and 964 MB, with a
  floor near 300 s). Merges are now priced at the measured request rate under
  load (`loadedRequestMs`, `formats/io/link.js`), but the HLCA-core export
  has not been re-run since.
- **Decode churn.** The export decodes whole chunks to keep a small part of
  them (2.7 GB for 0.2 GB on HLCA). The fix is the multi-range read
  described above; it has not been started. Chromium's transient renderer
  peak (3.55 GB, live end footprint 1.97 GB) comes from this.
- **Part boundaries.** A chunk that straddles two parts is fetched twice
  (WebKit 10M run: 128.6 MB of 763 MB).
- **Export reads still pass through the decoded-chunk cache.** A one-pass
  read could skip it.
- **No price for some matrices.** `matrixRowsCost` returns null for a matrix
  read by gene (CSC) or a modality's matrix. The export still runs, without
  a price line.
- **What an `.h5ad` export holds.** Only the expression rows are spooled.
  The cells' obs names (read in one call), their obs columns and their
  embedding rows are held in memory until the file is laid out.
- **A half share's export holds more than every cell's.** A 10M-cell CSR
  file at 50% uniform peaked the browser at 5.3 GB (renderer 4.1 GB), against
  3.8 to 4.6 GB (3.4 GB) for every cell, in the same 13.1 s (Chromium, two
  runs each, 2026-10-30); not repeated since, cause not found.
- **The empty Firefox `Error`** has not been reproduced. It is now shown
  with its class and cause.

(browser-anndata-copy)=
## A copy of the open dataset (`.zarr.zip`)

Built 2026-11-04 (`data/anndata/write/`, its `README.md`). The writer that
saves the dataset the page has open as one Cellucid-ready AnnData
`.zarr.zip`, with changes written in: its output profile is the
specification's {ref}`browser copy <browser-copy-profile>`. The History ⋯
menu's **Save AnnData with this state…** and the action
`history.anndata.save` call it, with the apply plan of the step shown as its
changes ({ref}`history-anndata-save-action`).

**Three layers**, each its own module, so another container or format plugs
in at its own level:

| Layer | Module | Owns |
|---|---|---|
| Container | `formats/zip/writer.js` | A stored ZIP written as a stream, ZIP64 exactly where a field overflows, the same bytes for the same entries. The figure export's batch archive uses it too, so its former 4 GiB / 65,535-entry ceiling is gone. |
| Array | `formats/zarr/format-v2.js` (shared with the transcoders), `formats/zarr/writer.js`, `formats/zarr/encode-task.js` | Zarr v2: the documents, dtypes, `vlen-utf8`; the byte shuffle, zlib level 6 by fflate, float16 bits; a streamed store (metadata after the chunks, the consolidated metadata the last entry, chunks in grid order); chunks encoded on the decode workers (`encode-task.js`) and put back in order. |
| AnnData | `data/anndata/write/` | The element layer the Zarr writers share (`element-writer.js`: the transcoders' kept stores and the copy), the specification's layout by role (`layout.js`), the streamed writer (`ready-writer.js`), the copy (`copy.js`, `gene-major.js`), the checks before a byte is written (`changes.js`, `names.js`), the way to the page (`acknowledged-sink.js`) and to every core (`pool-encoder.js`). |

**The way in.** `BaseAnnDataAdapter.writeAnnDataCopy(changes, {sink, signal,
owner})` → `AnnDataReaderClient.writeAnnDataCopy` → the reader worker. It
answers `{bytes, entries, zip64, chunks, layout: 'anndata' | 'mudata', nObs,
priorityOrder, notes}`. The sink takes one `Blob` per batch
(`COPY_BATCH_BYTES`, 4 MiB; at most `COPY_WINDOW_BYTES`, 16 MiB, unanswered),
a File System Access writable or the cell exports' spool; on a failure or a
cancel the caller (`app/history/anndata-copy-save.js`) aborts its writable,
never closes it, and removes its spool. Progress
is the copy's read activity (kind `cells`, of the job `owner`), the one
channel the page shows reads by. `changes` adds or replaces `obs` columns
(one value per original position, placed at each cell's file row), `var`
columns (by var row; `varModality` in a MuData), `uns` keys and
`uns['cellucid']` keys; `priority` and unknown keys are refused, and so is
any change that writes, renames or replaces
`obs['cellucid_original_position']` or takes its name
(`PRIORITY_COLUMN_REFUSAL`: it holds the rows' original positions of a file
in priority order, which the copy keeps); the whole set is checked
(`changes.js`) before the writer exists.

**What it does.** Every element of the presented store is copied (every
cell, whichever are open, rows in stored order, so a priority-ordered file
stays one), every value with its dtype and bits. Every expression matrix is
written gene-major: a `csc_matrix` array by array, a `csr_matrix`
transposed in passes within a budget (count, values, rows), with the
priority order's prefix counts computed from the rows, equal value for value
to `prepare()`'s; a dense matrix in strips of a band of rows and a group of
genes. A dataset of several modalities is written as a Cellucid-ready MuData.

| Date | Decision | Reason | Rules out |
|---|---|---|---|
| 2026-11-04 | Zarr format 2 | What `prepare()` and `cellucid_prepare()` write by default; the app reads it as fast as v3 | A v3 copy (no reader gains) |
| 2026-11-04 | fflate at level 6, not `CompressionStream` | Measured: 79 against 93 MB/s on indices, 23 against 23 MB/s on float values, fflate's output 5% smaller on values; it runs in a worker of the decode pool with a level the specification names | The native stream (no level, main-thread scheduling) |
| 2026-11-04 | Encode chunks on the decode workers, in order | Profiling put about 38% of a copy in compression (45% in reading the source) | One encoding thread |
| 2026-11-04 | Re-encode every chunk, never pass a compressed chunk through | Passing through needs a raw-chunk read on the Zarr and HDF5 readers and helps only a source already in the specification's chunks | Raw chunk passthrough |
| 2026-11-04 | Deterministic bytes | The same source and changes give the same archive however the work is cut (budget, slab, encoders in flight), so a copy can be pinned and compared | Output that depends on scheduling |
| 2026-11-04 | One element layer for the transcoders and the copy | No second AnnData writer. Consequence: a change to `format-v2.js` or `element-writer.js` changes kept stores: the code fingerprint of converted stores covers them, so a change makes browsers convert again (`tests/kept-code-fingerprints-contract.test.mjs`) | A copy-only writer |
| 2026-11-12 | Every writable stream of a page is opened through `utils/file-writable.js` (`openFileWritable`); no other module calls `createWritable` (`tests/file-writable-contract.test.mjs`) | WebKit writes the whole buffer under a typed-array view handed to `createWritable().write(view)`, not the view (measured: a 1 MB `subarray` of an 8 MB buffer, and a `DataView` of it, each wrote an 8 MB file; Chromium and Firefox wrote 1 MB). A view of part of a buffer is handed over as a Blob of its bytes | A writer that hands a stream a view itself |
| 2026-11-04 | Names exactly as stored; every name a change makes a path follows the portable rule; a name in two Unicode forms refused | An NFD column read back as its fill value (the reader normalized ZIP entry names, the Zarr metadata did not); the reader now compares keys in NFC too (`formats/zarr/open.js`), and the ZIP and folder stores find a key in either form | Normalizing names on write |

**Normalized, as anndata does:** numbers little-endian; a float16 NaN loses
its sign and payload (the reader drops them already); a record array's text
as fixed-width Unicode. **Left out, with a note:** an HDF5 attribute that is
an object reference, an element of `uns` the browser cannot decode, a
compound or opaque array. Anything else unreadable fails the copy in the
reader's words.

**Measured** (Node, `scripts/benchmark-anndata-copy.mjs`, peak = heap +
array buffers during the copy; 2026-11-04):

| Copy | Time | Archive rate | Peak |
|---|---|---|---|
| Suo (561,947 cells, 233 MB), one thread | 49.3 s | 4.7 MB/s | 151 MB |
| Suo, 8 encode workers | 24.7 s | 9.5 MB/s | 129 MB |
| Garcia (219,731 cells) as cell-major, 64 MB budget, 8 workers | 3.7 s | 31 MB/s | 251 MB |
| Garcia, same, one thread | 12.3 s | 9.4 MB/s | 142 MB |

Pooled and single-thread outputs were byte-identical; a transposed `X`
matched the source's gene-major `X` byte for byte.

**Memory** ({doc}`27_decisions_memory_model`): two slabs of the array being
copied (`SLAB_BYTES`), the budget (the memory model's `columnIndex` share: a
transposition's block of genes, or a dense matrix's band of every gene), the
chunks in flight and the batches posted, all counted in the reader worker's
working set under the part `anndata copy <ticket>` while it runs, plus the
archive's central directory (about 70 bytes an entry). Bytes a copy reads
are kept neither in the page's memory nor in kept parts
(`ReadActivity.readOnce`).

**Open:** it was not measured in a browser on a large dataset; prefix counts for a MuData modality whose
matrix is cell-major are not computed; a large dense `X` re-reads source
chunks once per group of genes; `showSaveFilePicker` with a file on disk is
not exercised by a test (the private-storage writable is, in
`history-anndata-copy.spec.mjs`). The highlighted cells' `.h5ad`
export still builds its HDF5 tree itself: an HDF5 subclass of
`AnnDataElementWriter` is the follow-up.

**Kept data removed once.** The shared element layer and an
`older-layouts.js` fix (a table anndata 0.7 kept in `uns` was not converted,
so a copy lost its categories) changed code the fingerprints of converted
stores and gene indexes cover: every user's converted R and MatrixMarket
stores and kept gene indexes are removed once on the first visit of this
version.

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| A copy of the open dataset | `cellucid/assets/js/data/anndata/write/` (`copy.js`, `gene-major.js`, `element-writer.js`, `layout.js`, `ready-writer.js`, `changes.js`, `names.js`, `acknowledged-sink.js`, `pool-encoder.js`), `cellucid/assets/js/formats/zarr/writer.js`, `cellucid/assets/js/formats/zip/writer.js` | `cellucid/tests/anndata-copy-contract.test.mjs`, `cellucid/tests/anndata-copy-cross-language.test.mjs`, `cellucid/tests/browser/anndata-copy.spec.mjs` |
| Writers, order, spools, price | `cellucid/assets/js/app/exports/cell-exports.js`, `cellucid/assets/js/utils/radix-order.js`, `cellucid/assets/js/formats/hdf5/writer.js` | `cellucid/tests/scale-boundaries-contract.test.mjs`, `cellucid/tests/round-trip-contract.test.mjs` |
| Buttons, card, Cancel | `cellucid/assets/js/app/ui/modules/cell-exports-ui.js`, `cellucid/assets/js/data/open-cells-work.js` | `cellucid/tests/cell-export-reads.test.mjs`, `cellucid/tests/open-cells-lease-contract.test.mjs`, `cellucid/tests/browser/cell-exports.spec.mjs` |
| Streamed rows, gathered rows, price, listed names | `cellucid/assets/js/data/anndata/reader/matrix-rows.js` (`readMatrixRows`, `getMatrixRows`, `matrixRowsCost`), `cellucid/assets/js/data/anndata/reader/names.js` (`readObsNames`), `cellucid/assets/js/data/anndata/matrix-access.js` (`pipelined`), `cellucid/assets/js/data/base-anndata-adapter.js` (`readMatrixRowsOf`, `getMatrixRowsCostOf`) | `cellucid/tests/cell-export-reads.test.mjs` |
| Read-once bytes, range cache | `cellucid/assets/js/formats/io/activity.js`, `cellucid/assets/js/formats/io/byte-reader.js`, `cellucid/assets/js/data/byte-range-cache.js` | `cellucid/tests/cell-export-reads.test.mjs` |
| Failure words, cancellation | `cellucid/assets/js/utils/errors.js` (`failureText`), `cellucid/assets/js/utils/abort.js` | `cellucid/tests/cell-export-reads.test.mjs`, `cellucid/tests/cancellation-and-errors-contract.test.mjs` |
| CSV format shared with Python and R | `cellucid/tests/fixtures/obs-columns.csv`, `cellucid-python/src/cellucid/anndata_session/_results.py` (`pull`), `cellucid-r/R/round-trip.R` (`cellucid_pull`) | `cellucid/tests/round-trip-contract.test.mjs`, `cellucid-python/tests/test_round_trip_contract.py`, `cellucid-r/tests/testthat/test-round-trip.R` |

## Related pages

- {doc}`12_figure_export_architecture`: figure export (PNG/SVG), its worker
  and its Cancel.
- {doc}`09_data_loading_pipeline_and_caching`: the reader's many-row reads;
  {doc}`22_decisions_link_model_and_request_planning`: link pricing in long form.
- {doc}`22_decisions_link_model_and_request_planning`: request rate and
  latency under load, S3 addresses, request splitting.
- {doc}`21_decisions_decoders_and_codecs`: the decode pool the churn lives in.
- {doc}`27_decisions_memory_model`: residency budgets and fetched-byte
  budgets.
- {doc}`24_decisions_graph_edges_at_scale`: the other long many-row read.
- {doc}`10_sessions_persistence_and_serialization`: saving the state rather
  than the cells.
- {doc}`20_read_path_architecture`, {doc}`01_codebase_map_and_entry_points`.
