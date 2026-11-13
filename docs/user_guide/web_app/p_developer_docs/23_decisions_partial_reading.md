# Decision record: partial reading

This page answers: **"Why does Cellucid read only part of a large dataset the
way it does, what was measured to justify it, and what was tried and taken
back?"** Partial reading (a share of the cells, a gene read for a sample of
the open cells, kept data reused on a later visit) is the main way the web app
gets fast on large files, judged first on public links. Its "What it does
now" is the one full description of what is read of a large file: gene
indexes, the gene sample, cell subsets, edges at a share, reading the open
file again and progressive opening; {doc}`09_data_loading_pipeline_and_caching`
is the reference that links here. The modules and invariants are mapped in
{doc}`20_read_path_architecture`; which cells a subsample keeps is
{doc}`19_subsampling_specification`. Update it in the same change as the code
it describes.

## At a glance

**Audience**
- Developers changing the share planner, the gene sample, gene indexes, cell
  subset reads, kept gene indexes or progressive opening.

**Time**
- 30 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture`, {doc}`19_subsampling_specification`.

:::{note}
A date written "by 2026-11-12" is the day the entry was first written into
these records, where the day of the decision itself was not recorded.
:::

---

## What it does now

All paths are under `cellucid/assets/js/`.

(partial-reading-gene-sample)=
### Shares and the gene sample

1. **The share-cost rule.** A share of the cells never costs more than a
   larger share of the same request, or than every cell. About s/n of every
   cell's cost is the aim, not a rule. In code, the "never more than every
   cell" half is the ceiling `wholeSampleCost` in
   `data/anndata/subsample-plan.js`: the bytes and requests every cell's own
   gene sample would make in the same time. A share's gene sample
   (`geneTimeFit`, `ceiling`) and the page's top-up
   (`AnnDataReader.geneSampleCost`, `geneTime`; `geneSampleCost`'s
   `geneTime.spentRequests`) stay within it.
2. **Exact cells.** A subsample is exactly the cells its request names. It is
   never made larger. A size under `FEWEST_OPEN_CELLS` (2, `data/subsample.js`)
   is refused with the same reason everywhere: one cell has no coordinate
   variation to draw. Values read through the subsample view are bit-identical
   to the same rows of a whole read (`views/cell-subset.js`, invariant 1).
3. **Auto-fit sizes the map, not the gene.** `{mode: 'auto'}` opens every cell
   unless the per-cell working set (`CELL_WORKING_SET_BYTES`, 77 B a cell) or
   the gene indexes every open cell holds would pass the device's allowance
   (`autoFitCount`). Then it opens the largest count that fits. On a file in
   its own order that count is a `blocks` subsample at the file's rung. On a
   priority-ordered file it is `uniform` under the file's own seed, so the
   cells are the file's first rows (`defaultSubsampleMethod`). A cell-major
   root matrix (CSR, dense in row chunks) is left out of this fit: its genes
   are read for a gene sample, weighed by its own fit (`geneIndexFit`).
4. **The block ladder.** The storage's block (`defaultBlockSize`) is the run of
   cells that spans `BLOCK_CHUNKS` (4) chunks of the dominant cell-major array.
   For an array read byte by byte it is the cells whose bytes fill
   `BYTE_READ_BLOCK_BYTES` (64 KiB). The layout alone sets it, so a file's
   default sample is a property of the file. `blockSizeLadder` doubles that
   block up to the cell count. `blockRung` picks a rung in this order: the one
   a gene index kept in this browser was built under (it looks for a kept
   index under the share's order at each rung, `keptOrder`), otherwise the
   lowest rung where a block fills its round trips, at the higher of the
   measured rate and a seed rate (the host's remembered rate, or this page's
   own downlink). The chosen rung is recorded with the cells (`blockSize`,
   `cellOrder`). The formulas are in {doc}`19_subsampling_specification`.
5. **The gene sample.** For a cell-major root matrix, a gene is read for the
   first `cells` open cells (`geneSampleOf`). `geneSampleOf` answers before
   any gene or row pointer is read (the open cells' view of a CSR matrix is
   described by its shape), so the page knows the sample at open and every
   analysis asks it of the source first. The count is the smaller of two
   fits. `geneTimeFit` takes the most cells within `FIRST_GENE_SHARE` (½) of
   `GENE_READ_TARGET_MS` (20 s, a response-time target, not a platform
   limit) over the measured link, chunk-index reads included; the fit's
   `budgetMs` is the whole 20 s. `geneIndexFit` takes the most whose gene
   index fits the `columnIndex` share of what the open cells leave. The fewer
   is the sample, recorded as `genes` `{cells, of, fit}` (`fit` one of the
   two, with `limit` `geneTime` or `columnIndex`). Open order decides what
   the sample reads:
   - Every cell of a file in its own order opens in the priority order of the
     file's default `blocks` subsample (`cellOrder`), so the first `m` cells
     are whole storage blocks and the sample is a nested prefix.
   - A `uniform`, `stratified` or `spatial` share of such a file opens in that
     same order with its own cells filtered in (`shareInCellOrder`,
     `priorityOrder` of the cell order's `blocks` over the share's cells). Its
     sample is its cells in every cell's first blocks, a subset of the same
     reads. The share's cells are unchanged; the order's rung is recorded
     (`cellOrder`), a request reading the same cells again names it (`order`,
     `requestOfOpenCells`), and the kept gene index files it with the share's
     specification (`openOrderSpecOf`, `orderKey`).
   - A `blocks` share keeps its own order. A share of a priority-ordered file
     opens in file order (`placeInFile`: its cells by file row, a uniform
     sample of them at every prefix, in the file's first rows); a share whose
     every open cell is read anyway keeps its priority order, as every cell
     then keeps file order.
   - Kept runs ending within the sample cost nothing (`breaks`).
   - When the same cells were open before (a share lowered and raised back),
     they reuse the order and sample they were read in (`seenGeneSampleOf`).
     They are not fitted again. Each opening of the same bytes carries the
     cells earlier openings opened, with their cell order and sample as they
     were when replaced (`AnnDataReader._geneSamplesSeen`, noted by
     `noteGeneSampleSeen`, keyed by every cell or by the share's spec): the
     same cells opened again with a request naming no sample take that order
     (a gene index kept under it is found), a person's count or every cell
     read as they were, and a gene-time sample fitted again over the link as
     measured now but floored at the cells it read and never every cell when
     it read fewer (`geneSampleOf`'s `seen`/`seenFitted`). A sample grown by
     the page's top-up keeps the fit that sized it (the page passes `topUp` to
     `growGeneSample`), and a share chosen while every cell was open is held
     to every cell's sample as it was read (`wholeSampleCost` `seen`).
   - A progressive opening's steps carry the fitted sample with its fit
     (`cell-subset-policy.js` `geneSampleRequestOf`: `genes: {cells, fit}`; a
     count alone is a person's or a session's choice, `fit` null), so an
     atlas grown to every cell is topped up too.
   - A gene-major matrix (CSC, a Cellucid-ready file) reads each gene for
     every open cell and has no sample.
6. **Pricing the sample's first gene.** Over a measured link
   `firstGeneReadCost` prices the index pass for the sample's rows from the
   layout alone — the stored chunks those rows lie in (each row taken to hold
   the mean stored values) for the indices (and values, unless read per
   gene) and the pointers, the requests their runs make, and the chunk-index
   reads that find them (`chunkIndexReads` of the store contract: an HDF5
   version 1 chunk B-tree's nodes not read yet, each counted once, a level a
   round trip; a sharded Zarr array's shard indexes), in dependent steps
   (`steps`: the pointers' index levels, their chunks, then the entries'
   index levels and chunks) — and `linkReadMs` turns that into time (a round
   trip per wave of requests the link runs at once in each step, or a round
   trip a step and the bytes at its rate). The round trip is the link's
   under load (`link.js` `loadedLatencyMs`, `loadedRoundTripMs`), not its
   shortest, and a request costs no less than the host takes one under load
   (`loadedRequestMs`, part of `requestCostMs`, so of `requestMs`); how the
   link is measured and why: {doc}`22_decisions_link_model_and_request_planning`.
   Every price of a read's time (the gene sample, Genes read for, an export,
   a progressive first view, an analysis run's read-ahead) takes the loaded
   round trip; the shortest shapes reads. The rows are priced as runs
   (`PricedRows`): a `blocks` order a block at a time
   (`subsampleBlockRuns`), a subsample's rows ordered once.
7. **The sample read, and grown.** The reader indexes a matrix in segments
   up to `genes.cells` only (`prefix-column-index.js`, `readRows`): a gene's
   column is NaN on the other open cells, a sparse column has no entry there,
   and `getMatrixCellSubset` reads a cell past it by its row (its own values,
   never the index's absence of it as a zero). `AnnDataReader.growGeneSample`
   reads it for more cells in segments doubling from the sample
   (`growthSegmentBounds`: the first as many cells as the sample, about the
   first gene's time; `geneSampleCost` prices them first, and when the gene
   index of every open cell is kept nowhere here names the most cells it is
   kept for, `index.fits`, which the legend then offers), keeping what is
   built when stopped; a request names the sample to reopen with (`genes`: a
   count, or `'all'`) and every cell's order (`order`).
8. **The top-up on the page.** After the first gene of a sample the gene
   time sized, `app/state/managers/field/gene-sample.js` (`DataState.topUpGeneSample`)
   spends the rest of `budgetMs` reading genes for more cells, each step
   through `readGenesForMoreCells`, sized from the reads' measured times
   (`cellsWithin`: a fixed time and a time a cell from the last two reads of
   different sizes; after the first alone, only a step at least as large, at
   its rate; the first gene's cells answered by a kept gene index,
   `fit.kept`, not counted), capped by `geneSampleCost`'s `index.fits`, never
   by a price, and for a share by the rule its first gene was fitted by: its
   sample transfers no more over the gene time than every cell's would
   (`geneSampleCost`'s `geneTime` option: `wholeSampleCost`, every cell's
   sample at the gene time over the link measured after the first gene, less
   what the sample transferred and requested so far; `ceiling` its bytes and
   requests, `target` the cells lowered to them), since a scattered share's
   steps each fetch nearly every stored chunk again. On a later visit the
   cells whose gene index was kept (`fit.kept`) count as read at the rate the
   first gene measured on the others, so the gene time buys no more cells in
   all than it would have read. A gene asked for never waits behind a step
   (`yieldGeneSampleTopUp`, called as every gene load starts: the step under
   way stops at once, the reader trimming its unbuilt segments on the abort,
   and the steps go on after the gene within the time left). A step adding
   under a quarter of the cells read is not made. An analysis holds the
   sample still before it runs (`holdGeneSample`). A person's **Genes read
   for** (item 9) takes over: the top-up stops, keeping what it read. The
   sample keeps its `fit` (still the gene time's: the reader's own,
   `growGeneSample`'s `topUp`; a person's count has none).
9. **Genes read for: one control, one path** (`app/ui/modules/gene-read-control.js`,
   in Session › Subsampling below Cells to open). The share of the *open*
   cells a gene is read for, on a slider that snaps to steps: the share read
   now, then 1, 5, 10, 25, 50 and 100% of the open cells above it
   (`geneReadSteps`), the last held to the most cells whose gene index is
   kept here (`AnnDataReader.geneSampleCost` `index.fits`, asked when the
   control is first used, not while the dataset opens: pricing reads the
   sample's row pointers). Its one status line prices the chosen step
   (`geneSampleCost(target)`: bytes and, over a measured link, time) and
   prices it again every 5 s while it is shown (`REPRICE_MS`) and after each
   growth step (the time left), so the estimate follows the link as it
   changes. **Read** prices the step once more and holds it to the index's
   room at that moment, then calls `readGenesForMoreCells`, which grows the
   sample in doubling steps (`growthSegmentBounds`), each step reading only
   the rows it adds: the cells read are kept, never read again. **Stop** (and
   the notification's Cancel) keeps what is read. 100% is what was *Read for
   all cells*; there is no other path. The gene's legend and an analysis's
   gene sample line are links to it (`createGeneReadLink`,
   `createGeneReadNote`, `revealGeneReadControl`, which opens the Session
   panel and its Subsampling section, `ui/modules/session-sections.js`
   `expandSessionSection`). The share is of the open cells, so it follows
   Cells to open: 100% is every open cell, and a share lowered there lowers
   it. A session records the sample read (`genes.cells`,
   `session/cell-subset-record.js`): a growth stopped, or under way when the
   session is saved, records the cells read by then.
10. **Unread cells** hold NaN, never 0, and are drawn faded beneath the read
    cells. The look is decided in {doc}`28_decisions_rendering_layers`.

(partial-reading-gene-indexes)=
### Gene indexes

1. **Gene indexes per matrix.** Switching the expression matrix (X, a layer,
   `raw.X`, in the data mapping) keeps each built gene index while every kept
   index fits the `columnIndex` memory budget, the most recently read first:
   switching back to a matrix already indexed reads nothing. An index still
   being built when another matrix is chosen stops (the genes waiting on it
   are still answered). A sparse index keeps 8 bytes an entry (a uint32 cell
   and a float32 value: 414 MB for Suo's 51.7 million entries); building it
   keeps the first pass's indices for the second in the smallest type that
   holds a gene position (2 bytes for up to 65,536 genes, 103 MB for Suo
   instead of 207 MB), and the values a waiting gene made the first pass
   read, so no slab is decoded twice.
2. **What the first gene costs, said before it is read.** How a matrix will
   be indexed is decided from its metadata before anything is read
   (`sparseIndexPlan`, the store factory's `plan`; `geneReadCost` in
   `data/anndata/matrix-access.js`, which a subsampled opening weighs too),
   and `AnnDataReader.geneReadPlan(element)` reports it: whether the matrix
   is stored cell by cell, the stored bytes the first gene reads, the index
   it builds and where that is kept (memory, disk, or nowhere, with the
   reason), and the measured rate. The gene's loading notification states it
   in one line, with the reasons as fact rows in its details
   (`geneReadPlanFacts`, `data/gene-read-plan.js`; the cause from the layout
   through one helper, `geneLayoutReason`/`geneLayoutFacts`, which the
   auto-fit notice and the index pass's read activity use too), and has
   Cancel: a build asked for by genes stops when the last of them withdraws
   — one round trip of the measured link later, so a gene picked instead of
   the withdrawn one joins the pass under way with its reads in flight — and
   its queued reads are never sent.
3. **A stopped pass goes on; a failed one starts over.** A stopped pass keeps
   the slabs it counted (`_pass1`: counts, kept indices, row pointers, values
   read whole for genes), and the next gene goes on from the first slab not
   counted; pointer windows end at fixed multiples, so the slabs are the ones
   a pass from row 0 makes (pass 2 finds them by number). A failed pass
   starts over. A matrix no store can hold is scanned once per batch of
   waiting genes, each scan its own read activity naming why
   (`onUnindexedScan`).
4. **The rows' pointers.** The kept rows' row pointers a subsample's first
   gene needs are read as one batch through the store's decoded chunks, as a
   read activity of their own (`_readyRows`), and only the rows the gene is
   read for: a gene sample's first cells (`readRowPointers` given `rows`,
   then `rowsMatrix`: the sample's rows as a CSR matrix of their own;
   `_openRowsMatrix`), never every open cell's; a stored matrix opened whole
   reads its first `rows + 1` pointers. The rows' pointers are read in one
   pass over the rows ascending (`CellRows.ascending`, made once), pointer
   ranges up to three pointers apart read and held as one (a range of its own
   holds as many numbers), and no run of the rows in open order is placed in
   them until something reads the rows in open order
   (`RowPointers.runIndexFor`). Pointers read once are kept together and
   copied from by every later read, piece by piece (a range only part of
   which is held reads only what it lacks; a grown sample reads only its
   added rows'), and a read of fewer cells than are open by their rows (a
   cell export, cells past the sample) reads only theirs (`_namedRows`). A
   raise reads the pointer chunks its added rows share with the open ones
   from those decoded chunks, as every gathered read of a subsample does
   (`readKeptRowPointers`, `inducedGraphLines`, a modality's cell map:
   `chunks`).
5. **Values beside the indices.** Over a link, values stored in chunks read
   whole (compressed or filtered) are read beside every slab's indices while
   a gene waits, and kept, when they fit beside the kept indices
   (`sparseIndexPlan` `buysValues`): the first gene was priced with them
   (`firstGeneReadCost`), and no later gene reads a value again; a raw
   array's values are read per gene by range.
6. **The subsample view** (`data/anndata/views/cell-subset.js`):
   - A read is bit-identical to the same rows of a whole read.
   - No chunk without a kept row is fetched.
   - A read holds its output and at most two slabs of decoded boxes.
   - Only the kept rows' CSR row pointers are read, and only those of the rows
     a read covers (`readRowPointers({rows})`).
   - Many rows are read in waves over the stored chunks in ascending order
     (`gatheredArray`, `values.js` `readRowsInStoreOrder`). The whole view of
     a column (obs, codes, `obsm`, dense matrices) is read row by row in the
     open rows' file order (`rowWaves`: the rows ascending and the view row
     of each, `CellRows.byRow`), each row one copy to its place, with no
     pieces listed or ordered. Runs are ordered by source row with a radix
     order (`rankedOrder`), never a comparison sort.
   - Describing the open cells' CSR matrix makes none of its arrays
     (`_requireMatrixElement` reads its shape; the stored matrix was checked
     at opening): its arrays in open order are its kept rows' entry runs in
     open order, which a pass in store order never uses. Row runs hold 32-bit
     positions wherever every one fits (`runsOfPositions`).
7. **Gene indexes in segments.** Rows in priority order are indexed in
   doubling segments (`prefix-column-index.js`; the first is 1/32 of the rows,
   at least a slab). A later opening of the same method and seed takes over
   the segments it shares (`carriedSegments`). Rows spread through the store
   (`rowsInterleaved`) are one segment in **store order** (`reader/gene-segments.js`
   `_segmentedColumns`, `_storeOrderMatrix`; `_inStoreOrder`), answered in
   open order through `fileOrderedAccess`, in the open rows' own file order
   (`CellRows.byRow`, made once, its rows the rows ascending). That segment
   reads each stored chunk once; the gene plan of such a pass makes no matrix
   of the rows in open order (`_segmentProgress`); the kept values and the
   share of index chunks they lie in are one walk forward
   (`keptValueStats`), and pricing a read lists no piece (`planBoxes`
   `pieces: false`). A gene index's pass in store order over a CSR matrix
   whose decoded chunks pass the decoded-chunk budget reads through a window
   of its own instead (`AnnDataReader._scanChunks`: the two chunks of
   `indices` and of `data` a slab shares with the next): each pass reads each
   chunk once, in order, and kept by recency the chunks are gone before a
   later pass comes back to them, so keeping them only held memory.
8. **Every stored row's index restricted to a share.** When the spread rows
   reach to the last open row, lie in every chunk (the whole reads no more
   bytes: `_readsNoMoreThanRows`) and their own pass is more work than every
   cell's (`wholeIsCheaper`), the segment is instead **every stored row's
   index restricted to them** (`_restrictsWhole`, `restrictedAccess`): every
   cell's pass over the stored matrix as stored (the same reads and decoded
   chunks), kept in this browser as every cell's index (either opening takes
   the other's), with no kept row's pointer read (`_readyGeneRows`) and no
   entry gathered. A gene is placed from it on the open rows (each value that
   is not zero by the open cells' lookup; once built, from the gene's entries
   alone). The work is weighed per entry over the index built and each of
   its genes answered once, at rates measured on the source as passes and
   genes run (`PassRates`: the index's own work per entry, `onWork`; the
   view's own -- its kept rows' pointers, rows put in file order, entries
   gathered -- its `ViewMeter`; each path's placement of a gene in open
   order, `placing('whole')` per entry of every stored row's,
   `placing('own')` per entry of the share's), carried to every opening of
   the same store; before any is measured, a share keeping at least two
   thirds of the entries reads the whole. Placing a gene from every stored
   row's entries tests each against a bit of the segment's open rows
   (`openRowBits`, 1.25 MB at 10M stored rows, held in cache) and looks up
   only those of open rows.

(partial-reading-genes-arriving)=
### Genes as they arrive

1. **A subsample's genes fill in as its cells are read.** A subsample opens
   its cells in an order whose every run of first cells is a sample of them
   (its priority order, the file's cell order, or file order:
   `subsample-plan.js`). Its cell-major matrix is indexed in segments of that
   order (`data/anndata/prefix-column-index.js`): the first about a
   thirty-second of the cells (at least one slab of the index's work), each
   next one as large as all before it, one after another, sharing one memory
   budget (each segment's store placed with what the built ones left,
   `withMemoryBudget`). A segment's `indices` and `data` are the subsample's
   gathered entries cut to its rows, read as one gather of the stored arrays
   (`gatheredArray` composes the two sets of runs, `composedRuns`), so a gene
   of k cells holds memory of their entries and the decoded-chunk budget,
   never of the file (`tests/subsample-gene-read-memory.test.mjs`). Their
   builds report into one read activity (`wholeGeneIndexActivity`): it
   expects the stored bytes of every unbuilt segment as it starts and counts
   each segment's progress as its share of their rows, so its bytes and time
   left are the whole gene's.
2. **Partial answers.** `AnnDataReader.getMatrixColumn(…, {onPartial})`
   answers a gene after each segment with its values on the cells indexed so
   far and NaN on the rest; the reader worker posts those as `{event:
   'partial', column}` before the reply (`getMatrixColumn`'s `progressive`
   argument), as the cells read so far in their sparse form when that is
   smaller (below), so a partial answer moves a gene's entries, not a whole
   column (40 MB at 10M cells, 200 MB at 50M); the client hands them to the
   callers that asked (`onPartial`), and the page draws them as colors only
   (`DataState.showPartialVarField`: the gene is not loaded until it is
   whole, so filters, analyses and sessions never see partial values; Cancel
   puts the colors of the field shown back). `geneReadPlan` says what remains
   to read and when the first colors come (`progressive`).
3. **Two forms of a gene.** A view colors by a gene's dense column (float32,
   one value per cell), made in the reader worker and transferred, so the
   page never expands one. Beside it the worker makes, in one pass, its
   sparse form when that is smaller (`data/column-forms.js`: `{length,
   indices, values, readRows}`, cells from `readRows` on not read, a gene
   sample's NaN tail), and the client notes it beside the column
   (`noteColumnForm`): the color kernels
   (`state/managers/continuous-colors.js`) and the range scan read the form,
   a word fill and the gene's entries, instead of two passes over every cell
   on the main thread, and a gene grown or cut in place keeps its form
   (`resizedForm`). A column is never written in place, so its form stays
   exact while it lives. A computation reads the same gene as its non-zero
   entries, `{length, indices: Uint32Array, values: Float32Array}`
   (`getMatrixColumnsSparse`; `DataLayer.readGeneColumn` →
   `DataState.readVarSparseColumn` → the var loader's `readSparseColumn` →
   `getGeneSparseValues`): the dense column's values bit for bit (a duplicate
   coordinate summed in stored order, a zero left out), about 45 times fewer
   bytes for a median Suo gene (0.2% non-zero), and nothing installed in
   field state. The analysis kernels take either form
   (`compute/gene-statistics.js`). Every typed array a worker reply carries
   owns its whole buffer (a view is copied out first), so replies are
   transferred, never cloned. What the page then holds of a gene:
   {doc}`27_decisions_memory_model`.

(partial-reading-cell-subsets)=
### Some genes in some cells

1. **Two routes.** Differential expression between two named groups reads
   one gene-major block of the groups' cells (`getMatrixCellSubset`): each
   gene's entries, with each entry's position in the cells asked for
   (ascending) and its value. It is read by whichever route decodes or copies
   fewer bytes (`readCellSubset` in `matrix-access.js`): the cells' rows (CSR
   or dense: every chunk the rows touch, since a chunk is decoded whole;
   nearby rows are read as one range, entries filtered into the block
   without a per-row expansion) or the genes' columns (a CSC matrix's
   ranges, or a gene index's entries in memory plus what is left of building
   it). Few cells of a large matrix go by cell, many cells, or any once the
   gene index is built, by gene. A subsample view's rows are priced at the
   stored chunks they lie in (`touchedChunkBytes`). An analysis's first pass
   on the first cells of priority order reads those cells' part of each gene
   the same way.
2. **A request is priced.** Every read of the cells' rows is priced at no
   less than a request of the store's link (`requestPriceOf` in
   `matrix-access.js`, from `formats/io/link.js`: its merge gap, and a round
   trip's worth of bytes, or the merge gap at every request the link runs at
   once until it is measured), in three places: rows nearer than the merge
   gap are read as one range even where an array is read in part (an
   uncompressed or contiguous `.h5ad`, whose chunks do not join rows), a
   read in flight counts as at least a request's worth of entries
   (`pipelined`), so no more requests are in flight than the in-flight
   entries hold at that price, and the route's cost (`readCost`) adds a
   request's worth for every range. CSC gene ranges and dense rows are
   priced the same way. Nor does an analysis ask for every gene at once: two
   groups' blocks are read as many genes at a time as the run's memory holds
   ({doc}`11_analysis_architecture`).
3. **CSC by gene.** A CSC matrix read by gene reads consecutive genes as one
   range of the indices and one of the values, several ranges in flight
   while the entries already read are kept, and keeps each gene's entries in
   one pass without a branch on which cells the subset holds; only a gene
   stored with a duplicate coordinate or unsorted rows goes through
   `sparseColumnOf` first. Each range of genes is then kept on a decode
   worker (`keepSubsetRun` in `subset-entries.js`, sent with
   `runOnDecodeWorker`: the range's arrays are transferred, the cells'
   positions are sent to each worker once as a context), while the reader
   thread reads and decodes the next, byte for byte the block the reader
   thread keeps.
4. **Genes an analysis computes never reach the page.** Marker discovery,
   differential expression, the gene-moments pass, the dot plot and UCell
   compute in the analysis workers, and the reader delivers their genes to
   those workers directly (`data/gene-delivery.js`): a run opens a gene route
   (`app/analysis/data/gene-route.js`), one `MessageChannel` per analysis
   worker, the reader holding one end of each and each worker the other; the
   page asks for a delivery (`DataLayer.deliverGeneColumns` /
   `deliverGeneCellSubset` → `DataState.deliverVarSparseColumns` → the var
   loader's `deliverSparseColumns` →
   `AnnDataReaderClient.deliverMatrixColumnsSparse`) and hears back only the
   genes' bytes, then hands each worker a task naming its delivery. A gene
   only the page holds (a user-defined copy) is read here and sent with its
   task.

(partial-reading-edges)=
### Edges at a share

A subsample's neighbor graph (the induced subgraph of the open cells) streams
into its edges a slab at a time, as a whole graph does; it is never held
whole. It is read in the stored graph's row order, never in the order the
cells were opened, by the cheaper of two reads (`AnnDataReader._graphReadPlan`,
priced as a matrix's rows are, `firstGeneReadCost`): the whole graph streamed
once from the first open row to the last, keeping the open rows, whenever
that reads no more bytes than the open rows alone (they lie in every chunk: a
uniform share of a graph in chunks of thousands of cells), else the open
rows' pointers and entries only (`readKeptRowPointers`, `gatheredArray`, a
block share). Each line's column is its row's rank among the open rows in
stored order, so the lines ascend and the edges are built and checked exactly
as a whole graph's, then placed on the open cells (`ConnectivityLines`
`relabel`, `relabelledEdges`): strengths and partners placed in one walk when
the four edge arrays held at once take no more than every stored line's edges
take to build (`wholeBuildBytes`), else strengths placed and let go before
the partners, about 1.5 times the edges' bytes, counted in `'graph build'`.
Each edge is counted at its smaller cell as its line arrives (the builder's
`cellOffsets`, while the walk waits for the next slab), each placement walk
is a typed loop of its own (no closure an edge, no count after the last
line), the open rows' stored order is the open cells' own (`CellRows`
`byRow`), and a stream's row pointers are asked for before the lines' order
and ranks are made. Over a link slower than a disk, slabs are read ahead
toward filling every request the link runs at once (`_inFlightEntries`),
within the decoded chunks' share of memory (`_graphReadAhead`); from a disk,
one ahead. A MuData modality's graph is read the same way in its own rows
(`placedGraphLines`, `readPlacedGraphRows`). The whole graph's build:
{doc}`24_decisions_graph_edges_at_scale`.

(partial-reading-reading-again)=
### Reading the open file again, and progressive opening

1. **The same reader worker, the same store.** Opening another cell subset of
   the open file (the subsample control, `app/ui/modules/subsample-control.js`;
   a session's recorded cells) happens in the same reader worker
   (`AnnDataFileSource.reopenCells` → `AnnDataReaderClient.reopenCells`, a
   new adapter on the same client) and reads through the store already open
   (`anndata-reader-worker.js`, `sourceKey`): its structure, chunk indexes,
   cached metadata blocks, decoded chunks and measured link stay, and only
   the cells change, so a new subset costs no request before its own cells'
   bytes. The bytes are the same when the source names the same version: a
   local file by its name, size and modification time, a remote file by its
   identity's key and version (`ETag`, provider revision) or its URL, size
   and `ETag`; a picked folder (a Zarr store, a Cell Ranger `outs/` or
   MatrixMarket folder) or a set of files by every file's path and that key,
   so a file added, removed or rewritten is another source; a Zarr store
   folder on a web host read at a version naming its content (a GitHub
   folder at a commit: the provider's key and commit) by that key and
   version, whatever the host lets the page see of its headers (GitHub's raw
   host exposes none to a page, so its validators read as absent); any other
   store folder on a web host (a served store, a folder link) by its URL
   while the validators of its consolidated metadata (`ETag` and
   `Last-Modified`, each the host names) are the ones it was opened with,
   asked again with one `HEAD` (`httpStore.unchanged`), or, when they moved,
   while what was read is what the store holds (a site deployed again:
   {ref}`link-model-kept-bytes`). That requires the metadata to move with any
   rewrite of the store: every Zarr writer writes the consolidated metadata
   again when it writes the store, so its write time (`Last-Modified`), or a
   tag made from it (`cellucid serve`'s `ETag`: size, times, inode), changes.
   Another version, a folder whose host names no validator, or a set of
   files one of which is known only by its size opens anew. A growth step of
   a store folder on a web host asks its version beside its own reads
   (`anndata-reader-worker.js` `open`): the step reads through the held store
   while the `HEAD` of `.zmetadata` is out, and is kept only once the folder
   is known unchanged; found changed (or failing while it was), it stops, is
   discarded and the folder is opened anew.
2. **What is carried.** The new reader also takes over
   (`AnnDataReader.open`'s `carry`) what the open one read that does not
   depend on the cells — CSR row pointers, a CSC matrix's column pointers —
   and the gene index of the open cells both share as a prefix: a
   subsample's cells open in priority order (`subsample-plan.js`,
   `priorityOrder`), so a larger `uniform` or `blocks` subsample under the
   same seed begins with the smaller one, and its index (built in segments of
   that order, `anndata/prefix-column-index.js`) is extended for the new
   cells only; a smaller one cuts it and reads nothing. A failed or canceled
   re-reading reads the cells open before again (carried, so it fetches
   nothing) and leaves the dataset on screen as it was.
3. **The open cells grown or cut in place** (`app/open-cells-resize.js`).
   Every reading of the open dataset names its own cells (a cell subset
   request passed to `reopenCells`; the next opening's request of
   `cell-subset-policy.js` is a link's, or the cells a dataset of My datasets
   was last opened with, for the next opening only): Apply during a
   progressive opening's step, and the earlier cells read back after a
   canceled one, never take or overwrite another's. Before reading, the work
   under way on the open cells ends on them (`data/open-cells-work.js`:
   field, gene, embedding and vector field loads, the gene list read after
   the first view, a data-mapping change, cancellable jobs such as an
   analysis run or outlier quantiles; a run's result arrives and is marked
   stale; and the operations that hold the open cells in use until they
   finish, `useOpenCells`: saving a session, the notebook's session request
   and the CSV and `.h5ad` cell exports); the change begins in the same turn
   as it sees no work under way (`beginCellsChangeOnceSettled`); from then
   until the commit, work about to start waits and starts on the cells after
   it (a field picked meanwhile loads on them). When the new opening's first
   cells are all of the open ones, or its cells are all among their first
   (`sharedOpenCells`, `anndata-file-source.js`: `reopenCells` answers
   `resized: {from, to}`), nothing shown is read again: the new adapter does
   not read the first view's embedding (`initialize({firstView: false})`),
   the reader keeps every embedding in the frame the cells were first shown
   in (`EmbeddingPayloads.adoptTransforms`) and the cell-order sums of the
   rows both openings share (`adoptCellOrderSums`, the session's cell check:
   {doc}`10_sessions_persistence_and_serialization`), and the page reads only
   rows `[from, to)` of what it holds — every resident layout
   (`loadCoordinateRows`), every loaded obs field and gene (the loaders'
   `readRows`; `BaseAnnDataAdapter`, "Rows of the open cells"). The reader
   answers `rows` reading only what those rows need where the layout allows:
   an obs column or embedding's row range, a cell-major matrix's index
   segments the rows lie in, a priority-ordered file's gene-major `X` between
   the prefix-count levels around the rows (`_placedColumnRows`,
   `_prefixEntryRows`). Everything the commit needs (the identity and
   priority of the new count, the live view's layout, the length of every
   added row) is checked before it changes anything; an embedding whose added
   rows hold a value that cannot be drawn is let go, as an opening of those
   cells refuses it (when it is the one shown, Apply shows the cells as a new
   reading and a progressive opening stops at the cells open). Then, in one
   synchronous step, the dimension manager's layouts (`resizeCells`) and
   every per-cell structure of `DataState` (`state/managers/open-cells.js`,
   which tabulates each) grow or are cut, and the viewer takes the new live
   points with the camera kept (`resizeLiveCells`). A failure anywhere before
   it reads the earlier cells back with their own request. Kept views are
   retired and made again from their record (the renderer sizes them by the
   count); a selection under way stays when cells are added and is retired
   when cells are cut; the neighbor graph, the vector field overlay and
   outlier quantiles in use are read again on the new cells; analyses mark
   their results stale. What is derived from the file takes the added cells
   as opening them would: a user-defined category made from a column (a copy,
   edits, a consensus) labels them as it labels its other cells of the same
   source category, and a category or range highlight group takes those in
   its category or range that are visible; labels and groups made by hand
   stay as they were. Each per-cell structure is bit for bit what opening at
   the new count directly holds after the same actions
   (`tests/open-cells-resize-state.test.mjs`). A dataset read again with
   other cells gives its newly open cells the same (`cell-subset-carry.js`,
   `takeInNewCells`): a category or range group keeps its type and rule on
   both paths.
4. **Progressive opening** (`anndata/progressive-open.js`). The same step
   grows a progressive opening. Its first view is sized by one rule for every
   source, a picked file as a web host: the measured link (a picked file's is
   its reads'), the stored bytes of the embedding chunks its cells lie in
   (`firstViewCost`), their decoding at the rate chunks are measured to
   decode in the reader (`decodingRate`; before any chunk has been decoded
   there, the first open cell's chunk is read to measure it when a read costs
   no round trip worth pricing, a read the first view makes anyway) and the
   page's measured cost of showing a cell (`data/opening-cost.js`). When the
   first view would take longer than about a second, the app opens the most
   cells predicted within that target (at least √n); it is taken only when it
   arrives at least twice as soon as every cell, then doubles to the cells
   asked for (`cellSubset.growth`):
   - Nested orders grow by prefix.
   - A share opened in an order its smaller sizes are not prefixes of -- one
     spread through a file in its own order, opened in the file's cell order
     (`cellOrder`), or a share of a priority-ordered file opened in file order
     -- grows through the first cells of that order (`prefix` requests,
     `subsample-plan.js` `prefixStep`, carrying the `order` and the fitted
     sample): each step's first cells and gene sample are the share's, so each
     grows the one before in place (`sharedOpenCells`) and the last is the
     share itself, also in place.
   - A subsample of a file in its own order whose first-view embedding is
     chunked much coarser than its blocks grows from few chunks instead
     (`progressive-open.js`, `chunkMajorOrder`, `chunkedPlan`): the same kept
     cells, taken chunk by chunk in a seeded random chunk order (each chunk's
     cells in priority order), each step a `chunks` request
     (`subsample-plan.js`, `chunkStep`) ending at a chunk boundary so each
     chunk is read once; then one `reordered` reading puts the same cells in
     priority order, every per-cell array moved with its cell
     (`DimensionManager.reorderCells`, `DataState.reorderOpenCells`), the
     frame kept. It is planned only over a measured link (any source's) and
     when its first view is predicted at least twice as soon as the priority
     order's (or as every cell's, when that is not progressive), each
     ordering priced by the embedding chunks its first view touches; the
     prefix plan is weighed against it as a priority order is. A blocks
     share, whose first cells in priority order lie in few chunks already,
     keeps priority order; a uniform share of a coarsely chunked file grows
     from few chunks.
   - A share grows only when its cells are all on the map no later than one
     reading of them would put them (`grownLater`, `chunkedPlan`; the
     owner's rule judged on completion, not on the first picture). A growth
     reads the same bytes in more readings, a loaded round trip each (and
     from chunks shows every cell again in its order), so over any priced
     link a share opens whole; it grows only on the page's measured cost of
     showing cells, with no link to pay.
   - Every cell of a file in its own order does not grow. It opens in its
     default subsample's cell order (for a gene sample) only when the sample
     brings the first gene sooner by more than the order costs the opening
     end to end: made and timed at opening, its runs found a block at a time
     (`blockOrderRuns`), and each of the `ORDER_OPENING_PASSES` passes over
     the open cells it adds (its copy in the reply, the page's cell identity,
     the first view's embedding and codes gathered and scattered) priced at
     that time (`subsample-plan.js` `orderOpeningMs`, `orderPays`: the
     sample's priced first gene, `geneTimeFit` `sampleMs`, against every
     cell's `neededMs`); else in file order, its genes read for every cell.
   - The plan is made after the store's structure is read, and opening reads
     no gene name (the var index is read when the names are first asked for,
     after the first view), so nothing else is priced as the first view's.
     Only a `uniform` or `blocks` subsample or a priority-ordered file opened
     whole is priced at all: the chunks of any other are never looked at.
   - A step's cells are recorded with the subsample's specification of their
     count, which does not give them again (`reproduced: false`; a prefix
     step's `part`), so a session saved during the growth opens exactly them
     (in the same order: the request names the cell order).
   - A later visit's first view is priced by each chunk's own bytes to fetch
     (`firstViewCost` per chunk once anything is kept, `storedFor` for
     decoding a kept chunk), and it never opens fewer cells than those whose
     chunks are all kept (`progressivePlan`).
   - A folder on a web host read at a content version (a GitHub folder at a
     commit) is taken over by every reading again without asking its host
     (`sourceKey`), so its growth steps stay in place. A source no later
     reading can tell to be the same bytes (a store folder with no key and
     whose host names no version: no `ETag` or `Last-Modified` on its
     metadata) is never opened progressively: the reader worker drops the
     plan and opens every cell asked for at once (`anndata-reader-worker.js`
     `cellsToOpen`; `open-cells-folder-reuse`).
   - A store whose matrix is read after its map (a MatrixMarket folder,
     `transcode/matrix-market.js`) is read in **one pass for the cells the
     opening grows to** (`AnnDataReader.shareRows`: the first `fitCells` of
     them, `_grownShareRows` planned with the growth), not for each step's.
     Every step reopens the same bytes and joins the pass under way or reads
     the share kept (`MatrixStream.ensure`: a pass serves every caller asking
     as many rows or fewer, stops only once none waits, and more rows asked
     while a smaller pass runs replace it). A completed larger share replaces
     the smaller, whose directory is let go once no reader reads it
     (`letGo`, `space.js` `letGoTranscodedStore`): one share kept a folder
     version. Every other format read once into a store of its own opens
     progressively from that one store; a step reopens the same bytes, which
     the reader worker takes over (`sourceKey`). A file a host serves only
     whole is downloaded once by the page before it opens. The gene index of
     a cell-major matrix (CSR `.h5ad` and Zarr, a 10x `.h5`), kept in this
     browser, grows by the cells each step adds. Gene-major matrices (Loom, a
     CSC `.h5ad`) are read in place gene by gene, with no derived store.
5. **The cost of showing a cell** is what the prediction needs: each
   opening's time from its cells in the page (their GPU stores made next) to
   the first frame that draws them (`recordFirstFrame`, from
   `reloadActiveDatasetInPlace`; a frame that drew none of them, or cells
   replaced first, records nothing). An opening also costs what showing any
   dataset costs, so the measure is the time a cell adds, the least-squares
   slope of the records' time against their cells (through the origin while
   they hold one size, or when the line would give an opening a negative
   cost), each earlier record counting half as much as the next, so no single
   slow moment keeps sizing later openings; kept in the browser under
   `cellucid:opening-cost:v2` as the records' decayed sums. An in-place
   resize is not recorded: the GPU keeps the cells kept and takes only those
   added, while every open cell's colors and filters are made again, so its
   time is neither its added cells' cost on an opening nor its open cells'.
6. **One subsample cache per opening.** Every size the opening resolves --
   the cells asked for, the first view, each step, and the cells asked for
   again at the last -- is selected through one subsample cache the reader
   holds for the opening and carries to each step and to Apply
   (`subsample.js` `createSubsampleCache`): the cells (a `blocks`
   subsample's blocks, a stratum's cells) are ranked in one pass for the
   largest size, each smaller size is a prefix of that rank order and its
   priority order the prefix itself, and a priority-ordered file's uniform
   steps are prefixes of its computed row order, so a step costs its own
   cells, not the file's.
7. **Decoded chunks carried between steps.** A subsample's cells lie across
   the whole file, so every step's added cells lie in chunks earlier steps
   decoded: a subsample's view reads every array decoded a chunk at a time
   through whole decoded chunks the store keeps (`anndata/decoded-chunks.js`,
   carried from each opening of the store to the next, let go with the
   store), within the one decoded-chunk budget the store's own cache shares
   (`formats/io/decoded-budget.js`: together they never hold more than the
   `decodedBlocks` share; chunks stay until something needs their room, the
   allowance or the working set changes, or the page is hidden:
   `trimDecodedChunks`), and a step, or a column read again, copies its rows
   out of them instead of decoding them again. Each read of an `obsm` array
   names the open cells it read and what for in the read activity (`open
   cells 7,435–20,000 (embedding, cells added)`).

### Kept gene indexes and revisits

`kept-column-index.js`: each built run is kept in the origin-private file
system, filed by file identity, matrix, shape and order key (method, seed and
settings, never size). It is checked before use: manifest, file rows, the
file's row pointers read now, and a CRC-32 per gene piece. A store-order run
is used whole or not at all. A run whose values were all read is kept with
them (`sparse` layout), so a revisit reads no value again. A run used is read
as an index built in the opening is held: whole, in one read of each array at
its first gene, when the gene-index memory budget holds it (counted in that
share), then every gene from memory; otherwise the genes asked for together
are read as one range per run of neighboring genes, never one read per gene.
A run whose reading is priced above building it again (`rebuildMs`, at the
read rate measured in the thread) is not used. A revisit's gene sample
(`geneTimeFit`) is the kept prefix (the kept runs ending within the count the
target prices): its first gene reads no new cell, so coming back is never
slower than the first visit. The count priced at the link as measured or
remembered (`pricingLink`) is `growTo`, the page's top-up's first step
(`topUpGeneSample`, in the background, yielding to genes asked for); later
steps are sized by the reads. Over a link with neither (not read ahead,
`linkPriced` false: an opening served wholly from kept data, host speed
forgotten) a new cell has no price: `neededMs` null, no `growTo`, and the
first step is a quarter of the sample. A change to code a kind's fingerprint
covers (`transcode/kept-code-fingerprints.js`) removes that kind of kept data
from every browser once; nothing is raised by hand.

(partial-reading-strings)=
### Strings and categories at a share

- A one-dimensional column of variable-length strings (an obs index, a text
  column) reads the rows asked for alone (store contract `readRows`): of an
  HDF5 column the references of the chunks they lie in, then of each global
  heap collection those rows point into the prefix up to the furthest of
  their objects (the chunk's other references, and its neighbors', name the
  objects before it; HDF5 `HeapStream`), a collection let go once the rows
  pointing into it are decoded; of a Zarr column the chunks they lie in, held
  as stored (bytes and an offset an item) with only the rows asked for
  decoded. A share of such a column decodes and holds its own strings, not
  every string of the chunks it touches.
- The byte reader never fetches bytes already in flight. A running read
  keeps its gap bytes until it ends, within the fetched-bytes budget.
- A categorical obs column whose dictionary has more categories than there
  are open cells, and more than one request of the link reads (its stored
  elements past the link's `maxRequestBytes`; a smaller dictionary is read
  whole, every category listed, the unused ones too), is coded by the
  categories its open cells hold (`AnnDataReader._categoryCoding`,
  `_readPresentCategorical`): its codes are read, then the labels and
  `uns['<key>_colors']` entries at the stored positions they name alone
  (`readCategoriesAt`, the dictionary gathered by those positions through
  `gatheredArray`) -- unless the rows read are enough to name a category in
  every chunk of a chunked dictionary (more than c · ln c rows for c chunks),
  when the dictionary and its colors are read whole beside the codes, in one
  plan, and the present ones picked (no round trip waits for the codes; no
  byte or request more than reading the dictionary); the field's categories
  are those in stored order, its codes positions among them, and the worker
  hands the page `present: {count, stored, colors}` (the adapter's
  `_adoptPresentCategories`, `presentCategoryPalette`: each category's color
  is the file's or Cellucid's at its stored position). The coding is decided
  on the field's first read and kept: a clear of the reader's cache keeps it,
  an opening grown or cut in place carries it (`_takeCarried`), so added rows
  are coded with every category found plus the ones they add, in stored
  order, and the page grows the field by label (`category-growth.js`); the
  page marks such a field `discovered` (a saved category it lacks is one no
  open cell holds, `pending-categories.js`) with `storedCategoryCount` for
  the legend's line.

### The subsampling specification

It is shared byte for byte with Python and R:
{doc}`19_subsampling_specification`.

## Decisions and why

### Shares and the gene sample

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| 2026-10-27 | `blocks` is the default method, aligned to about four chunks of the dominant array | a run of neighbors is a run of storage: a block share reads about its share of the bytes | a uniform default on a file in its own order (10% of HLCA's cells touch 45% of `X`'s chunks) |
| 2026-10-27 | a subsample opens in its priority order | its first `m` cells are a sample at every `m`, so growth and segment indexing append | rereading on growth |
| 2026-10-28 | a share that rounds below two cells is refused with the reason, never made larger | exact cells: sessions, links, Python and R name the same set | silent changes to the cell set |
| 2026-10-28 | large row-stored atlases open **all cells the device can hold**. A gene is read for a stated sample sized to the gene-time target, and the other cells wait until the user chooses *Read for all cells* (its cost shown). Analyses say which cells they used | the map costs bytes a cell; a gene of a CSR file costs each cell's whole row | sizing the map by what a gene costs |
| 2026-10-28 | the gene sample is the first `g` open cells in the opening's priority order. Every cell of a file in its own order opens in its default `blocks` order. A session records `genes` only while a sample is in effect, and one saved with every cell read restores with a fitted sample | the sample is then exactly a `blocks` subsample, and kept indexes of earlier visits count | a silent multi-GB read on restore |
| 2026-10-29 | genes on a sample: unread cells fade, read cells drawn on top and larger ({doc}`28_decisions_rendering_layers`) | unread cells should not look like 0 | the same gray for 0 and not read |
| 2026-10-30 | speed and smaller downloads outrank memory, as long as memory is counted and bounded | the owner's priority | trading time or downloads for memory |
| 2026-10-30 | **the share-cost rule** in every metric (network, memory, rendering, analysis, opening); s/n is the ideal, not a rule | a person choosing a smaller share expects it to be cheaper | any planner path where a share costs more than the whole |
| 2026-10-30 | the rule is judged on the costs a person meets: time, download, memory, rendering and analysis work. The request count is a means | a 1% share's obs names made 575 requests against every cell's 530, but took 30.6 s against 37.9 s and downloaded 75.7 MB against 77.9 MB | spending bytes to save requests (see Tried and reverted) |
| 2026-10-30 | the first gene is priced to half the gene time; the rest is read by measured top-up steps | the opening's link is a lower bound (CELLxGENE atlases priced at 20 s were read in 6-15 s; on Census S3 the opening's resolve measured 0.96–1.1 MB/s where the first gene then measured 11 MB/s) | one fixed sample from an early price |
| 2026-10-30 | block size comes from one ladder: the base from the layout, the rung from what is kept and the link, and a revisit keeps its rung | a size taken from the open-time link stayed at the storage block on first visits, and revisits then chose another size and read their cells again | sizes that do not nest |
| 2026-10-30 | spread shares open in every cell's block order, with their own cells filtered in | in its own order each of a 1% share's first cells lay in a chunk of its own: 71 cells in 10.7 s, against every cell's 1,572,326 in 5.5 s (10M-cell CSR `.h5ad` over an emulated S3 link) | a share's sample costing more than the whole |
| 2026-10-30 | levels of detail rank a cell by its identity under the share's seed | the same cells come first in every order, step and revisit | ranking by open position |
| 2026-10-30 | a share lowered and raised back comes back in its order with the gene sample it was read for (the gene time's, grown by its top-up, or the count a person read), never fitted again | fitted again at that moment's measured link and decode rate, a 50M store's sample of 4.49M became every cell (markers 15.4 → 62.3 s, a reader worker of 8.7 GB) and a 10M store's markers went 7.8 → 80.7 s; reused alone, an opening's sample of 68 cells (priced on the link's first small reads, its top-up cut short by the share change) came back as 68 where the measured link reads 13,290 | re-pricing the same cells on each Apply; reusing a sample without its fit |
| 2026-10-30 | a run that reads genes first stops and settles the gene sample's background top-up (`DataState.holdGeneSample`); a person's *Read for all cells* is left alone | the top-up grew the sample under a running DE, which then failed as invalidated | runs racing the top-up |
| by 2026-11-12 | a progressive opening's steps carry the fitted sample with its fit (`geneSampleRequestOf`) | its last step's request named the count only, and the top-up never ran after a progressive opening | a count without its fit between steps |
| by 2026-11-12 | a share's sample, and the top-up after it, never makes more bytes or requests than every cell's sample in the same time | on Tabula Sapiens a 1% share in the cell order, bounded by bytes alone, read 3.4 MB in 434 requests for 83 cells where every cell's sample made 351 | bounding a share's sample by bytes alone |
| 2026-10-31 | a file on this computer, or a host on it (a loopback origin: `cellucid serve`, Python's and R's served stores), prices the first gene as this machine at rest whatever its link measured: at no less than a disk's rate (`LOCAL_POLICY`, 1 GB/s), at its shortest round trip, no time a request beyond it (`host-link-memory.js` `pricingLink` `thisMachinePricing`; a local file's link snapshot carries its priced rate and `onThisMachine`) | an opening's metadata reads measure a fraction of such a link (a 10M-cell store on a loopback host in Chromium: 45 MB/s from 0.9 MB, Zarr v3 9 MB/s from 51 KB); priced at that, every cell's first gene (read in about 6 s by the code before the change) was priced at 15-37 s, and genes were read for 32-68% of the cells, the first differential expression on them ({doc}`22_decisions_link_model_and_request_planning`) | the gene sample on this computer from the opening's lower bound; pricing a host here as measured once a read slowed by the machine's own load counts as filling the link, or at its round trip under load (5 of 6 openings of 10M-cell stores with the machine busy sampled 19,779-1,947,096 cells, the first DE following the sample); a notebook kernel priced as this computer (it may run across the internet) |
| 2026-10-31 | the gene sample is grown from one control, **Genes read for** in Session › Subsampling (below Cells to open): a percent of the open cells on a slider snapping to steps (the share now, 1, 5, 10, 25, 50, 100%, held to the index's room), its cost said and priced again as the link changes, Read / Stop; 100% is the former *Read for all cells*, and the legend's and analyses' lines only link to it | owner, on HLCA ("Gene read for 3.8K of 2.3M cells", only *Read for all cells*, 9.1 GB, ~11 h): not reading all is fine, but the share must be choosable in the UI; reading changes coloring and analyses alike, so it sits with the reading controls, one place | a second size control (in the legend or an analysis); a count slider (the percent of the open cells reads as Cells to open does) |
| 2026-10-31 | growing the gene sample reads only the rows of the cells it adds, in doubling steps, so it never costs more than reading that share from scratch (the share-cost rule applied to growth) | the cells already read are kept; a step that read them again would make a growth cost more than a fresh read of its size | regrowing from 0 to a new size |
| 2026-11-06 | every rate measured of the machine or the link that a plan weighs is read through one module, `utils/machine.js` (`measureOf`): the link the reader plans on (`plannedLink`), the decoding rate, the order time, the kept-read rate, the index pass rates, the opening cost, the cores. A test or a benchmark holds them (`holdMachine`) on the page, and a thread hands what it holds to every worker it starts (`worker-script.js` `handToWorkers`) | the pins beside each measurement (`pinDecodingRate`, `pinOrderMadeMs`, `holdShownCells`) reached no worker from a page, and the link, the kept-read rate and the pass rates had none, so browser specs planned on whatever a loaded runner measured | a pin of its own beside a measurement; a spec asserting a plan on a measured link |
| — | the open cells a session, page address, Share link and Copy link name are read from the request the opening ends with, not from the step shown | a session saved while Suo grew reopened its first 131,072 cells, and its link named them (`data/cell-subset-policy.js`) | naming the step shown |
| 2026-10-27 | the subsample specification is part of what a link shares: the one feature state that lives in the address (`app/subsample-url.js`) | the specification gives the same cells everywhere, so a link reproduces the opened cells | any other feature state in URL parameters |
| 2026-11-01 | the gene-read cost is that of the matrix genes are actually read from (`defaultExpressionPath`) | pricing `X` while genes are read from `raw.X` opened every cell of a file whose first gene then read the whole of `raw.X` (an anndata 0.6 file on Dropbox: 122 MB for one gene) | pricing the matrix the file calls `X` |
| 2026-11-01 | a file's default sample is a property of the file, so the byte-read block (`BYTE_READ_BLOCK_BYTES`) is a constant, never the link as measured | a merge gap measured from one read's timing gave the same `.h5ad` blocks of 28 to 66 cells from one opening to the next | a block size from the link |
| — | a raise to 100% (Apply, `{mode: 'all'}`) names no sample and is not fitted again over the link as measured then | a 50M-cell store's sample of 4.49M became every cell after a round trip (markers 15.4 → 62.3 s, a reader worker of 8.7 GB); a 2M-cell `.h5ad` over an emulated 40 MB/s link went from 1,690,543 to 352,680 cells in another block order | fitting a raise again at the link as it is |
| — | while the matrix answers for cells newly open (a text matrix still arriving), the gene sample is the opening's, never every cell | the page said null (every cell read) and a gene picked then was refused as having no value in any cell (a browser-test run, Linux Firefox, `text-matrix-stream.spec.mjs`) | a null sample while the matrix is arriving |
| — | a gene asked for before a text-matrix share arrives is read from the share, not from the map's `X` | the map's `X` held no entry, so the gene was NaN on every cell and the page refused it | reading a gene from the map's `X` |
| 2026-11-10 | the Sample dataset listing is asked for once the open dataset has its first view | asked for as the interface was built it ran a millisecond before startup opened the default sample, so it read every sample's file, the one being opened too, against that first view | listing the samples beside the first view |

### Gene indexes and reads at a share

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| 2026-10-30 | a 50% share's slow first gene is fixed by cutting the CPU of the store-order pass, not by opening in file order | file order would make a larger share on a later visit read its kept index again | the file-order design |
| 2026-10-30 | a read on the gene sample checks its row pointers once a run: the most first rows known covered are kept per matrix (pointers held only grow) | markers at 50M read whole columns a few genes at a time (about 3,750 reads), and each read walked the sample's rows twice to ask whether its pointers were held: about 290 s on the reader alone | a per-read coverage walk |
| by 2026-11-12 | a subsample's first gene reads only the rows' pointers it is read for (`readRowPointers` given `rows`) | an 11.4M-cell CELLxGENE atlas's 2,011-cell sample read its whole `indptr` before, 31 MB in 4,379 requests: the writer lays the pointer chunks between the matrix's, so none merge | every open cell's pointers for a sample |
| by 2026-11-12 | a segment's entries are one gather of the stored arrays (`composedRuns`), not reads through the subsample's estimated chunks | read through the estimated chunks, each slab was tens of thousands of reads at once, and a gene of 3,401 cells of HLCA (4 billion stored entries) held gigabytes in the reader worker for 20 MB of entries | a slab of reads per estimated chunk |
| by 2026-11-12 | a store-order pass whose decoded chunks pass the budget reads through its own window (`_scanChunks`) | kept by recency the chunks are gone before a later pass comes back to them, so keeping them only held memory: a 50% uniform share of 5M cells held 132 MB above every cell's | the shared decoded-chunk cache for a scan |
| by 2026-11-12 | describing the open cells' CSR matrix makes none of its arrays | its arrays in open order are 2.5M runs and 140 MB for a 50% share of 5M cells, which a pass in store order never uses | building the view's arrays at describe |
| 2026-10-31 | a radix sort's digit has no more values than the keys sorted (`radix-order.js` `digitPlanOf`, 16 bits at most; the passes and the digit width chosen for the fewest counts walked) | a gene's few hundred entries put in open order (`fileOrderedAccess`) walked 65,536 counters a pass: 6 s of a 30,000-gene DE at 1M | a fixed 16-bit digit |
| 2026-10-31 | every stored row open in another order (a share raised back to every cell in its order) takes every stored row's kept index for the rows past the carried ones (`_restrictsWhole`, kept runs past `start` only), its genes placed by the open rows' runs when there are fewer runs than a gene's entries (`restrictedAccess` `entriesByRun`: no lookup scattered through memory, no sort) | the rows past the share's were indexed again by a pass of their own (10M Zarr: 27 MB read and decoded, marker genes 2.7 → 5.3 s when no gene was shown before) while every cell's index was kept; placed by a lookup of each row and a sort, a markers run's 30,000 genes cost 1.0 s (Node), by runs 0.6 s | a second index of every cell after the round trip |
| 2026-10-31 | a share's choice between every stored row's index restricted to it and its own pass in store order weighs, besides the index's and the view's work, each path's placement of a gene in open order (`wholeIsCheaper`, `PassRates.placing`): per entry over the index built and each gene answered once; unmeasured, the whole for two thirds of the entries or more | every stored row's index places every stored row's entries of each gene, twice a half share's: a 50% uniform share of a 10M-cell store took 1.3-1.7x on it what its own pass took (Chromium; first gene CSR 7.7 / Zarr 6.4 s against 6.0 / 4.5 s, marker genes 11.6 / 12.4 s against 8.0 / 7.4 s; every cell's first gene 6.0 / 4.9 s), the same values | weighing the index pass alone (a gene was taken to cost the same after either) |
| 2026-10-31 | a gene placed from every stored row's entries tests each against a bit of the segment's open rows before looking its open row up (`openRowBits`) | the open rows' table is 40 MB at 10M stored rows and each lookup a cache miss; half of a half share's entries are of rows not open (Node: 18.8 → 13.8 ns an entry) | a lookup of every entry |
| 2026-10-30 | categories read the present codes only when the dictionary is larger than the open cells | a 1% share of a column of 1,029,016 categories sent every label to the page | whole dictionaries at every share |
| by 2026-11-12 | a subsample's neighbor graph is read in the stored graph's row order, by the cheaper of the whole graph streamed and the open rows alone | a uniform share's open order is scattered through the file: read that way, every slab spanned the whole graph and decoded its chunks again (a 2.5M-cell uniform share of a 5M-cell `.h5ad` took 19.9 s for its edges against 2.0 s for every cell's) | reading a share's graph in open order |
| 2026-11-01 | a gather of a gathered array reads through the inner gather's stored chunks and its decoded chunks (kept rows per stored chunk), never through its own estimated chunks | on CELLxGENE's HLCA a gene of 3,401 cells held gigabytes of reads in flight for 20 MB of entries when each of a slab's thousands of estimated chunks was a read of its own | a read per estimated chunk |
| by 2026-11-12 | kept entry runs are held in typed arrays, at most one a row run | lists grown an item at a time held a 50% uniform share's 2.5M runs as boxed numbers | lists of boxed numbers |
| — | analyses read a gene the views do not show as its entries, not as a dense column | a dense column per gene held 20 GB for 50 genes at 100M cells (a scale check) | a dense column per gene |
| — | a run's first read over some cells keeps every gene of the cells; a CSR matrix read by rows in batches does not read each row of the cells again on every batch | a 10M-cell file with 30,000 genes read 470 times about 130 MB | one pass of the cells' rows a batch |
| — | the rows of a sparse `obsm` array are made only when read: describing the data mapping loads `obsm` deep enough to list each entry and makes none of a sparse array's `indptr`, `indices` and `data` | the share's CSR view made them on that load, before its kept rows' pointers were read, so every subsample of such a file failed ("could not be read again with these cells") | building a sparse array's children at describe |
| 2026-11-11 | an obs column of 64-bit integers holding values past 2^53 is read exactly as text in every reader and format and listed as a categorical field, saying why in one line (owner decision); a column of 64-bit integers that are all safe stays numbers (float64) | such values are identifiers, never measurements, and no number holds them exactly | reading identifiers as rounded numbers |
| — | stored outlier quantiles are accepted for a category of any size; under a subsample, stored ones describe the whole dataset and are gathered for the open cells | above 2^24 placed cells float32 no longer tells neighboring ranks apart, and recovering a rank with `round(q * n)` refused quantiles the rule wrote correctly: from about 17M cells in one category every opening read and ranked the whole latent space | a size above which quantiles are refused |
| 2026-11-13 | the gene index's code fingerprint follows the reader's methods through `this` (`scripts/kept-code-fingerprints.mjs`, `compositionEdges`), so `reader/row-pointers.js`, whose pointers a kept run is written with, is covered; the reader's own core is left out (it plans how a matrix is read) | after the reader was divided into modules composed onto one class, the closure followed imports only and missed the module of the row pointers, so a change to how they are made would have kept stale runs | a fingerprint that follows imports only |

### Cell subsets and analyses

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| 2026-10-28 | every read of the cells' rows is priced at no less than a request of the link (`requestPriceOf`) | unpriced, a two-group subset of a 10-million-cell uncompressed `.h5ad` (two categories, a million cells scattered through the file, about 10 entries each) put a read per cell in flight at once: a million requests, each with its promises, plan and pieces on the reader's heap, which grew to its 4 GB limit within ten seconds and took the tab down ({doc}`27_decisions_memory_model`) | a request per row |
| by 2026-11-12 | an analysis's genes go from the reader to the analysis workers directly (`data/gene-delivery.js`), never through the page | relayed through the page, a run's genes had cost its main thread more than the messages: WebKit counts every array a worker sends the page as the page's own allocation, so a run over every gene drove the page's garbage collector, whole-heap collections included, while it ran | relaying genes through the page |
| 2026-10-30 | a cell-subset run (differential expression) is priced over the run, not read by read: the gene stream passes the run's remaining genes (`runGenes`) down to `getMatrixCellSubset`, which prices rows every read, an index pass once, or the cells kept once, and takes the cheapest | read by read, every batch of a 10M run re-decoded the cells' rows (about 64 GB decoded for one run) | a second ask before keeping the cells; charging an unbuilt index's whole pass to every batch |
| 2026-10-31 | a cell subset is read through the first rows of each gene (`_sparseColumnRows`) only where those rows are fewer than the gene's (`_readsFewerColumnRows`: segments ending before the last, or prefix counts below every row); any other read takes the route priced over the run, and a read within a gene sample takes it by the sample's columns (only its rows' pointers are read) | a share raised back to every cell opens in the share's order with its gene index in segments, and two groups of every cell read each batch through the first rows of every gene: whole columns, every gene of the read at once, unpriced (1M: DE 9.0 s against 1.5 s; with no gene index yet, a read of 8,610 genes failed after 3 min 50 s for lack of memory) | the first-rows read for cells that reach the last segment |
| 2026-10-31 | a run's cell subset is checked once (`_checkedCells`: distinct, in range, its largest) and compared on each later read | a 10M DE reads 1.9M cells about 400 times; checked, bounded by the gene sample, by the rows read and by a prefix each read, they were 8.5 of the reader's 9 s | walking the cells again each read |
| 2026-10-31 | once a run's cells are kept (`_keptCellSubset`), a later read of them is a copy and nothing else (`_heldCellSubset`): no row view, pointer check or named rows again | every cell open with a gene sample (a local 10M CSR file: 8.1M sampled), a DE's groups reach past the sample; each of ~470 reads still made a row view of 2M rows and walked its pointers' coverage three times before copying (0.6 s a read): DE 4 s in the code before the change, 180–360 s with it | per-read work on kept cells |
| 2026-10-31 | a gene stream over some cells starts a read once half the memory gate is free, or every task left fits (`gene-stream-runner.js` `roomToRead`) | a read over some cells costs its cells whatever its genes; started as each task left the gate, a 10M DE of two groups (1.9M cells) made ~15,000 reads: 54 s against 4 s when it read whole columns, a task a worker a read | a read a task |
| 2026-10-31 | a subsample view's rows are priced at the stored chunks they lie in (`matrix-access.js` `touchedChunkBytes`: a view chunk touched costs the stored chunk it is part of, decoded whole) | priced at the view's part of each chunk, a 10% uniform share's groups (lying in every stored chunk) read by their rows cost 5× what was priced (10M CSR: 161 MB priced, 750 MB read and decoded) | pricing a view's rows at its own part of each chunk |
| 2026-10-31 | a run keeping a subsample's cells by their rows reads them in store order (`reader/cell-subset-reads.js` `_subsetInStoreOrder`: the cells' file rows ascending through the open rows' store-order matrix, each stored chunk decoded once, the entries put back at the cells' positions) | read in open order, a spread share's rows made each read of the run's pipeline decode nearly every stored chunk again: a 50% uniform share's DE groups (944,727 cells) read 7.5 GB of a 2.6 GB 10M CSR file (39–41 s; now 0.18 GB, 5.7 s; every cell 0.22 GB, 5.2 s) | reading a share's rows in open order |
| 2026-10-31 | a run's genes per task follow the genes' bytes as read (`run-tuner.js`: the densest form until a read is measured, sized again at each read and each shrink of the budget) | held at the densest form (8 B a cell of 1.9M) against what the run's context leaves of the share, a 10M-cell differential expression went one gene a task: 30,000 tasks, deliveries and page messages (traced: 30,040 against 168 before the change), every-cell DE 1.2x slower; marker genes likewise | a task size fixed at the run's start |

### Opening, growth and reading again

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| by 2026-11-12 | the first view reads only the embedding the dataset opens on; every other is read when a view first shows it | Suo before its first render: 5.5 MB in 27 requests, not the 10.9 MB in 139 that reading all three embeddings costs (Measurements) | reading every embedding at open |
| 2026-10-31 | a matrix read after its map (MatrixMarket) is read in one pass per opening, sized to the cells a progressive opening grows to; a pass serves every caller, stops only when none waits, is replaced by a larger one asked meanwhile, and the share it replaces is let go (one kept per folder version) | each growth step (499 → 998 → 1,996 → 3,000 cells) asked for its own share, and a pass reads the stream to its share's last cell, nearly all of it in priority order: the 3,000-cell spec folder read its matrix 4 times and kept 4 shares; a 1M-cell folder 11 times (996 MB of a 90.6 MB `.mtx.gz`), its genes at 24.8-25.8 s against 6.8-7.5 s at once (now 1 pass, 90.6 MB, 6.5-6.8 s: Measurements) (CI WebKit, where the folder opens progressively, failed `text-matrix-stream.spec`) | a pass a step; a smaller pass let finish before a larger one (the rest of one reading and a whole other); a canceled caller stopping the pass others wait for |
| by 2026-11-12 | a Zarr store folder read at a content version (a GitHub folder at a commit) is taken over by the provider's key and version, whatever the host shows of its headers | GitHub's raw host shows a page no validator: every reading of a GitHub folder, a growth step included, was opened anew, and every step read other cells than those shown (Suo: growth stopped at 32,768 of 561,947 cells, and its revisit at 750) | reusing a folder only by its host's validators |
| by 2026-11-12 | a source no later reading can tell to be the same bytes is never opened progressively (`cellsToOpen`) | the Garcia sample on such a host stopped growing at 48,860 of 219,731 cells | a growth that cannot tell its steps read the same bytes |
| by 2026-11-12 | a share grows only when its cells are all on the map no later than one reading would put them (`grownLater`) | a growth reads the same bytes in more readings, a loaded round trip each: a 50% uniform share of 10M cells on a local web host, grown from chunks, had all 5M at 4.5 s, whole at 2.8 s | growing whenever the first view arrives sooner |
| by 2026-11-12 | every cell of a file in its own order opens in its default subsample's cell order only when the gene sample saves more than the order costs end to end (`orderPays`) | a sample of 9.84M of 10M cells over a local web host saved tens of milliseconds while the order cost the opening hundreds and made it grow in two steps (all 10M at 3.4 s against 1.5 s in file order) | an order taken for any gene sample |
| 2026-11-11 | a growth step of a store folder on a web host (its version its metadata's validators) asks that version beside its own reads, not before them: the step reads through the held store while the `HEAD` of `.zmetadata` is out, and is kept only once the folder is known unchanged; found changed (or failing while it was), it stops, is discarded and the folder is opened anew (`anndata-reader-worker.js` `open`) | each step waited one round trip in series for the check (Suo on 3G: 0.3 s × 5 of a 40 s opening, about 4%); nothing read of a changed folder is ever shown, as with the check first: a change before the check's answer is seen by it, one after it would be missed by a check made first too | checking once per opening (a folder rewritten between steps would be mixed into one dataset) |
| 2026-10-28 | anything kept is filed under a version (content name, else size plus modification time or a strong ETag), checked before any kept byte is used | a changed file must never be served from what was kept | keeping by size or name alone |
| 2026-11-10 | a stopped opening stops its reads where they are | a superseded Apply kept reading the first view's embedding for about 15 s on a 1.6 Mbps link, because the reader stays open when other cells are read again and the read joined no caller that could stop it (`data/anndata-file-source.js`) | a read no one can cancel |
| — | a file picked from this computer is not copied to the reader worker a second time when other cells are read again: the worker reopens the source it holds | on hosted Linux WebKit the second copy read as empty ("WebKit encountered an internal error") where the held copy read on | sending the source again |
| — | an Apply in Cells to open that cannot read the cells asked for resolves `{shown: false}`, never rejects | an Apply that rejected was an unhandled rejection the person never saw (a browser-contract run, Linux WebKit) | a rejection no one handles |
| 2026-10-30 | the cell-rows lookups trade memory for speed where the ranks are concerned: the ranks stay for the next size (no second rank pass) and are counted in the memory model (owner decision) | a second rank pass to save the bytes | a second rank pass |
| 2026-10-30 | opening every cell costs no more than it did before gene samples and cell orders: a blocks order finds its runs a block at a time (no walk of 10M positions, twice); an embedding's display bounds are the transform of its stored bounds; its cell-order sums are taken after the payload is handed over | the opening of every cell is the baseline the share-cost rule is held to (the passes table below) | — |
| 2026-11-06 | kept directories are read, removed and rewritten under a lock of their own (shared for reads, exclusive for removal and for the in-place rewrite of the last-used time), taken after the directory's own lock (`transcode/space.js` `readingFiles`, `removingFiles`) | WebKit logs a read of a file removed or rewritten after the read began as a failed load ("WebKitBlobResource error 1"), whatever the caller does with the error: the Stored listing read directories a reader's opening had removed (two hosted macOS WebKit runs of the web contracts, `kept-parts.spec.mjs`), and the last-used time a later visit's reader writes was logged in one run of 120 | a read racing a removal or a rewrite |

### Opening other formats

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| 2026-10-26 | data with no cell embedding is refused, not shown with a computed one: a 10x matrix, a Loom file, a Seurat object, a MatrixMarket folder; the refusal says how to compute one in Python and in R, or to pass coordinates already at hand (`data/anndata/views/refusal.js`) | Cellucid computes none; `prepare()` and `cellucid_prepare()` write one that exists and refuse data without one | computing an embedding in the viewer |
| 2026-10-26 | a Cell Ranger `outs/` folder shows the folder's own analysis (UMAP, t-SNE and PCA projections, clusterings) as embeddings and obs columns (`views/tenx-outs.js`) | the matrix alone has no embedding | — |
| 2026-10-26 | h5Seurat (SeuratDisk's own HDF5 layout of a Seurat object) is recognized and refused with what to do (`views/h5seurat.js`) | SeuratDisk has been unmaintained since 2023, and an `.rds` of the same object opens directly | — |
| 2026-11-01 | when a file holds several single-cell objects the person picks one (`transcode/session.js`) | — | choosing for the person |

## Measurements that justify it

Real atlases were read from CELLxGENE Census over S3 (HTTP/1.1, several
mirror hosts). The machine was an Apple M1 Pro laptop; browser runs used real
GPUs (Chromium on Metal). "Node" means the app's own reader in Node over the
same real network, with six connections per address as a browser has. Suo is
562k cells × 8,192 genes.

| Dataset (cells) | Engine, date | Step | Before → after |
| --- | --- | --- | --- |
| HLCA core (584,944) | Chromium, 2026-10-28 | whole map; gene sample; first gene; later genes | 17.7 s / 7.7 MB (first picture 10.2 s); 2,143 cells; 11.9 s / 14.5 MB; 0.3 s / 0 MB |
| full HLCA (2,282,447) | Chromium, 2026-10-28 → 10-29 | first picture; whole map; first gene | 17.5 → 8.6 s; 33.4 → 13.7 s; 25.7 s / 1,384 req → 11.6 s / 17.7 MB / 950 req |
| mouse embryo (11,441,407) | Chromium, 2026-10-28 → 10-29 | first picture; whole map; first gene | 27.8 → 10.7 s; 47.6 → 22.7 s; 81 s / 5,379 req → 12.1 s / 13.9 MB / 913 req |
| mouse embryo | Node, 2026-10-29 | first gene once only the read rows' pointers are read | 27.2 s / 39.5 MB / 4,872 req / 1.34 GB → 7.3 s / 13 MB / 863 req / 0.72 GB |
| core / full / embryo | Chromium, 2026-10-30 | revisit: first picture; requests to open | 8.4 → 1.6 s, 9.0 → 2.0 s, 9.4 → 2.8 s; 189 / 455 / 921 → 0 (version checks only); later genes 0 |
| full HLCA / TS (1,136,218) / embryo | Node, 2026-10-30 | rung; cells with the gene at ~5 s and ~20 s | 19 → 152, 15 → 30, 20 → 40; embryo 1,465 @ 5.6 s → 1,400 @ 3.4 s and 5,232 → 8,235 @ 20 s; TS 3,125 → 6,792 @ 15 s |
| full HLCA | Node, 2026-10-30 | before the top-up → with it | 2,729-3,826 cells at 6-8 s and no more → 1,767 @ 4.6 s, 11,055 @ 17.5 s |
| HLCA core (Census S3), revisit from a kept index | Node / Chromium | a new cell priced as free → at the remembered rate | Node refitted 1,518 kept cells to 102,271 and read its first gene in 147.8 s, 575 MB, 14,680 requests → Chromium read 5,400-5,800 new cells with it in 6.6-7.9 s where the first visit's took 2.9-3.1 s |
| synthetic 10M x 30,000, `.h5ad` CSR / Zarr v3 sharded, loopback web host | Chromium real GPU, 2026-10-31, same hour against the code before the change | every cell: first gene, cells without a value, first DE top genes | before 6.6-6.9 s / 6.1-6.2 s, 0 cells; with the lower-bound pricing 4.8-5.0 s / 0.8 s for a sample, 2.4-6.7M cells without a value, Zarr DE rank 3 G02402 (before: G24308); after 5.1-5.2 s / 4.4 s, 0 cells, sums exact, DE top 3 G00000, G00001, G24308 as before |

Share against every cell, first gene, measured in the same hour:

| Case | Every cell | Share |
| --- | --- | --- |
| 10M CSR `.h5ad`, emulated S3 (230 ms, 20 MB/s, HTTP/1.1), Node, 1% uniform | 347,913 cells, 8.04 s, 18.7 MB, 193 req | before 53 cells, 11.13 s, 4.4 MB, 280 req; after 3,488 cells (= s/n), 8.16 s, 18.4 MB, 193 req |
| TS, real S3, Node, 1% uniform | 1.0 | time 0.40-0.56, bytes 0.16-0.24, requests 0.50-0.65 |
| full HLCA, real S3, Node, 1% uniform | 1.0 | time 0.78, bytes 0.19, requests 0.82 |
| TS / full HLCA, Chromium, before spread shares used the block order | 8.1 s / 6.1 s | 12.2 s / 17.5 s (the break that design fixes) |
| full HLCA / TS, Firefox, same run | 7.1 s / 10.3 s | 13.5 s / 6.1 s |

Synthetic stores (Node unless named; 20 ms / 50 MB/s HTTP/2 emulated unless
named):

| Case | Result |
| --- | --- |
| 10M CSR `.h5ad`, 10% uniform first gene | 517 MB → 266 MB (= every cell's) once spread rows became one store-order segment |
| 10M, 50% uniform first gene, `genes: 'all'` | CSR 16.8 → 7.50 s (every cell 6.76 s); sharded Zarr v3 10.16 s (7.75); dense rows 4.89 s (3.77); bytes ≤ every cell's on all three |
| 4M CSR `.h5ad`, 50% uniform, first gene, after the pass's CPU cut (pointer pieces, rows ascending, no open-order matrix) | 5.06 → 4.69 s, every cell 4.83 s; on a link faster than decoding the share's own pass took about 1.4x every cell's first gene, through the whole every cell's |
| 5M CSR `.h5ad`, emulated 20 ms link, the reader worker's handler, first gene with `genes: 'all'` | 50% uniform 3.98 s, peak ArrayBuffers 922 MB, RSS 1,582 MB → 3.86 s, 711 MB, 1,383 MB; every cell 3.49 s, 599 MB, 1,119 MB |
| 5M CSR `.h5ad`, local, k = 15, edges of a share read in stored row order | 2.5M uniform share 19.9 → 2.2 s (every cell 1.9–2.5 s), 0.5M 3.7 → 1.2 s; edges bit-identical |
| 10M opening, every cell in a cell order (150 ms / 5 MB/s HTTP/1.1) | first view 344 → 48 requests and open 1,512 → 659 ms, once rows were read in store order |
| 10M growth (`createSubsampleCache`) | priority file growth steps 566 → 46 ms; uniform 50% growth 1,310 → 186 ms; blocks 50 % 827 → 125 ms |
| HLCA-shaped `.h5ad`, 1% blocks of 19 cells, decoded chunks carried between steps | the last step decoded 33 embedding chunks and 2 of `cluster`, not 84 and 46 |
| 10M, a 50% uniform view read row by row in file order | first view 4.7-7.0 s → 0.9-1.0 s, every cell's 0.9 s |
| 50M sharded Zarr v3, Chromium (Metal, loopback) | markers 15.4 s (a 1,800 s timeout before the row-pointer coverage fix), DE 10.2 s, signature 2.2 s, dot plot 1.2 s |
| 50M and 10M Zarr after a share round trip (10% then 100% in place), Chromium | 50M DE 70.3 s (reader 1.8 GB, no allocation failure), markers 62.3 s, signature 4.5 s; 10M DE 83.6 s, markers 80.7 s, signature 1.6 s |
| Node, the reader on a 50M gene sample (9.49M cells), whole columns 8 genes a read | 2,000 genes 22.0 s before; 30,000 genes 3.4 s / 83 MB after the coverage fix (10M: 1,600 genes 19.4 → 2.6 s) |
| Node, 10M CSR `.h5ad`, DE on 1.89M cells over a simulated link | about 13-17 ms a gene (400-500 s for 30,000 genes) → 30,000 genes 14.2 s / 192 MB once priced over the run |
| Node, 2M CSR `.h5ad`, emulated link, share lowered and raised back | 40 MB/s: raise 352,680 cells at rung 15,624 (priced again, 53 s) → 1,690,543 cells at rung 7,812 with the first opening's fit; 15 MB/s: 1,695,490 → 1,086,048 (= the first opening); the raise's first gene reads nothing new |
| 10M rich CSR `.h5ad` DE in Chromium (gene sample 8.35M) | 218.6 s (321 s on a sample, 660 s without one, before); markers 4.5 s |
| 10M-cell revisit, DE from a kept index read gene by gene | 71.7 s against 16.0 s from an index built in memory (now read as one range per run of genes) |
| 3M Zarr, 1% uniform obs names | heap 931 → 8 MB |
| a column of 1,029,016 categories (1.2M-cell CELLxGENE atlas, real S3), 1% uniform | labels sent to the page 1,029,016 → 12,050, worker heap 216 → 94 MB, same 38.9 MB downloaded |
| 2,058,652-cell MatrixMarket stream (GEO, Node) | first view 19.8 s / 106 MB, against at least 7 min for the whole 3.95 GB stream |

Reads of genes and cells (Suo unless named; Node and a local `.h5ad` unless
named):

| Case | Result |
| --- | --- |
| Suo, first view | before its first render, 5.5 MB in 27 requests (its metadata and its 3D embedding), not the 10.9 MB in 139 requests that reading all three embeddings costs; over HTTP/2 at 80 ms and 2 MB/s its first render takes 3.8 s in WebKit (6.9 s reading all three) |
| Some genes in some cells: all 8,192 genes, index built | the two largest cell types (103,736 cells) 8.8 s → 0.69 s (by gene), two types of 16,376 cells 0.92 → 0.30 s, one type of 61 cells 16 → 15 ms (by cell); the run's peak memory 3.6 → 2.0 GB |
| Priced requests: a synthetic uncompressed CSR `.h5ad` of 200 genes and 10 values a cell, Chromium | 10M-cell differential expression finishes in 2.8 s with the reader's heap at 4 MB; at 1M cells 5.8 → 0.8 s (unpriced, the same read held 1 GB of heap for 5.5 s) |
| CSC by gene: every gene of Suo, 51.8 million entries | 26 reads of each array instead of 8,192 |
| CSC by gene: the two largest cell types (118,905 cells), every gene, Chromium | 0.82 s → 0.36 s for the read alone (decoding the 1,580 chunks it touches takes 121 ms on the decode pool; the rest is the reader thread keeping 51.8 million entries and copying them out of their chunks), 1.16 → 0.62 s in the benchmark's sequence of steps |
| Keeping on a decode worker (`keepSubsetRun`): the rebuilt Suo (5,103 genes), the same cells, Chromium | 0.31 s → 0.21 s, against 0.15 s for reading the whole of `X` (2.1× → 1.4×) |
| Two forms of a gene, 10M cells | a gene shown 66–70 ms → 6–10 ms, a partial answer 146 → 6 ms ({doc}`27_decisions_memory_model`) |
| k-NN rows from the page's edges | 20M edges, a 65,536-cell frontier: 0.18 s, against 1.7 s from a 10M-cell file |
| An index built ahead of need beside an extraction, Chromium | Suo's graph after 1,000 cell rows 0.80 → 0.47 s |
| Reading the open file again in the same worker | a CSR gzip `.h5ad` over 80 ms: the 31 requests and 0.8 s of its structure; the lung atlas: its 130 requests and 10–30 s, before its own cells' bytes → none |
| The cost of showing a cell, Chromium | about 25 ms an opening, against 9 ns a cell: 10K cells took 27 ms, 5M 70 ms |
| Ordering every cell for the gene sample, on a browser trace | made in 11 ms, 60-70 ms more |
| Full HLCA, Census S3, Chromium on the GPU, 2026-10-31 | the first gene was read for 6,449 of 2,282,447 cells (0.28%); Genes read for offered 1% (22,824 cells: 70 MB, ~25 s) up to 240,000 cells, the most whose gene index this device keeps (994 MB, ~5 min), not 100%; opening and the first gene downloaded 49 MB |

Growing the gene sample against a fresh read of the same size (Node, the
reader: `growGeneSample` from the sample to 10× it, after its first gene,
against an opening asked for that many cells, `cellSubset.genes`, and its
first gene). Synthetic CSR store of 200,000 cells × 2,000 genes, 40 entries a
row, zlib chunks of 65,536 entries, served over an emulated link of 0.4 MB/s
and 40 ms a request (2026-10-31). The growth reads only the rows it adds:
sample plus growth is the fresh read's bytes (3.89 + 28.91 = 32.80 MB), and
the growth alone is 0.88× of the fresh read in bytes, requests and time. Its
estimate before Read (`geneSampleCost`) priced the growth at 30.9 MB and
77.3 s: 7% and 6.5% above what it read and took.

| Store, link | Sample → grown to | Grown: bytes, requests, time | Fresh read of that size: bytes, requests, time |
| --- | --- | --- | --- |
| 200K cells, every cell (blocks order), 0.4 MB/s | 19,654 → 196,540 | 28.9 MB, 222, 72.6 s | 32.8 MB, 255, 82.1 s |
| 200K cells, 50% uniform share, 0.4 MB/s | 9,875 → 98,750 | 29.2 MB, 224, 73.2 s | 33.1 MB, 257, 82.8 s |

A MatrixMarket folder opened at once and progressively (showing a cell
measured at 2 ms, as on CI's software GL), before and after one pass a
progressive opening (Node, the reader worker in-process, files on this
computer, the origin-private file system on disk; 2026-10-31). The 1M-cell
folder is `scripts/efficiency-floor/scale-matrix.py --format mtx --cells
1000000 --per-cell 20` (a 90.6 MB `matrix.mtx.gz`, 30,000 genes), two
interleaved rounds; the 3,000-cell folder is `text-matrix-stream.spec`'s
(a 282 KB `matrix.mtx`, 40 genes). Matrix bytes include the map's read of
its first lines; requests are every file's reads; kept are the converted
stores (the map and the shares); genes is the time until every open cell
holds its genes.

| Folder, opening | Before: passes, matrix bytes, requests, kept, genes | After |
| --- | --- | --- |
| 1M cells, at once | 1, 90.6 MB, 45, 2, 6.8-7.5 s | 1, 90.6 MB, 45, 2, 6.4-6.7 s |
| 1M cells, progressive (first view 1,000, 10 steps) | 11, 996.5 MB, 295, 12, 24.8-25.8 s | 1, 90.6 MB, 45, 2, 6.5-6.8 s |
| 3,000 cells, at once | 1, 348 KB, 8, 2, 0.58 s | 1, 348 KB, 8, 2, 0.49 s |
| 3,000 cells, progressive (499 → 998 → 1,996 → 3,000) | 4, 1.19 MB, 17, 5, 0.94 s | 1, 348 KB, 8, 2, 0.59 s |

Opened progressively, the folder costs what opening it at once costs, byte
for byte and read for read; before, each step's pass read the whole stream
again (each a share of the first cells in priority order, spread over the
file) and kept its own share. An R object read once into a store of its own
(`seurat5_gzip.rds`, `sce.rds`, `Seurat_3.0.0_pbmc_small.rda`) reads the same
bytes grown in 4 steps as at once, one store kept; a kept gene index grows by
the cells each step adds (`kept-column-index-contract`, nested reuse across
k).

On the test machine, Playwright's WebKit reports every public host as
offline, so public links were walked in Chromium and Firefox only. WebKit
runs used local files.

### Passes over a share's rows

Why each pass over the rows of an opening, a view or a gene read is shaped
as it is. Dates are where the source gave one; the machine is an Apple M1
Pro, Node 26 or Chromium as the row says, and a row that names neither was
not tied to one.

| What | Conditions | Result |
|---|---|---|
| Row pointers of the kept rows only (`readKeptRowPointers`, `rowsMatrix`) | a 1% subsample of CELLxGENE's HLCA (2.3M cells; int64 row pointers in gzip chunks of about 1,250 rows) | the whole row pointers of the three matrices were 19 of the 20 MB the opening fetched, so a subsample reads only its kept rows' pointers, never those of cells it left out or of a matrix it does not read; the first gene of a 2,011-cell sample of an 11.4M-cell atlas read the whole `indptr` (9,146 pointer chunks, which CELLxGENE's writer lays between the matrix's, so none merge): 31 MB in 4,379 requests, where the pointers of the first rows alone need about a hundred of the chunks |
| Row lookup for the whole view of a share of more than an eighth of the axis | a 50% uniform share of 10M cells | a table of the axis (one write and one read a row in file order) about 40 ms, the radix order 132 ms (machine and engine not named) |
| Ordering the runs of a subsample | a 10% share, a million single-row runs | a comparator sort took 1.7 s of its first view, more than opening every cell: runs are sorted as numbers and each placed at its start's rank |
| Ordering run starts | a 50% uniform share of a 10M-cell file (5M single-row runs) | a sort of the values with a binary search placing each run took 2.3 s of the first view, more than every cell's whole first view; a stable radix order over the 32-bit words is linear in the runs |
| Ordering rows ascending | up to the largest row | a bit a row is linear in the rows and the file, where a sort of a priority-ordered selection of 50M rows takes 1.4 s |
| `copyRows`: a loop or a view and a copy | Node 26, float32 | 16 values: 21 against 37 ns; 32 values even; 64 values: 75 against 41 ns, so below a threshold a loop wins |
| Chunks of the rows of a share found once per chunking (`chunksOfRows`), not once a column | a 50% share of 10M cells | 70 ms a column, its pricing as much |
| A 2-D embedding's two values copied as two assignments | 5M rows of a 50% share | a loop a row was a third of the placement |
| A whole view of an array read row by row in the rows' file order (`rowWaves`) | a 50% uniform share of a 10M-cell file | 4.7-7.0 s for the first view through the pass over pieces (5M single-row pieces placed through three permuted arrays a row) against 0.9 s for every cell, the same bytes |
| The rows path only where runs are mostly single rows | every cell opened in a blocks order is 153 runs of 65,536 rows at 10M cells | the rows path took them as 10M rows (0.6 s of an opening) |
| Reading by source order | every cell of a 10M-cell file in blocks order | the embedding read stretch by stretch, each stretch's chunks lying apart, took 318 requests, against 22 to 34 (two counts of the same opening) read in source order (`efficiency-floor.test.mjs`) |
| A copy in file order and a scatter after it | a 50% share of 10M cells | a pass more a column than every cell's read, so each wave is placed in its view row as it is read |
| Pointer reads by a cursor (`PointerCursor`, a walk forward) | a 50% uniform share's 5M rows in open order | two lookups a row took 1.3 s of its first gene |
| `keptChunkShare` by halving, chunk by chunk, where rows are fewer than the chunks' span | a 50% share of 10M rows | a walk of every row took 50 ms of an opening |
| The inverse of the kept rows (`KeptPositions`) | every cell of a 50M-cell priority-ordered store | 1.07 GB as a table (2^27 slots of 8 bytes) against 200 MB as an array, used once an array holds fewer bytes than the table; its graph's first rows read in 2.26 s, then 0.18 s |
| Placing each cell of an order by a bit per cell of the file instead of a search of `kept` (`data/cell-identity.js`) | 10M open cells of 20M, Node | 0.08 s against 3.3 s |
| `keptRanks` (the bits of the open cells) against a search of `kept` a cell | a 50% share of 10M cells put in priority order in place | 2.4 s for the search |
| The cell-order digest keeps its words in locals (`cell-order-digest.js`) | a 10M-cell opening | an array of three cost the digest a third of its time |
| Cell-order sums taken after the payload is handed over, a slice at a time (`embedding-payload.js`) | 10M cells | hashing each cell cost a first view 0.1-0.15 s, and the next request waited 154 ms behind the sums |
| `placedBounds` as a loop a width with its bounds in locals | 10M-cell openings | a call and an array a value cost tens of milliseconds a pass |
| Pointer lookups by bucket (`RowPointers`) | the pointers of 5M kept rows in view order, a 5M-cell induced subgraph | a binary search over millions of ranges misses the cache at nearly every step: 11 s of the build; bucketed (about one range a bucket) a lookup searches a few, and the pointers are read in one pass over the rows ascending, ranges a few pointers apart as one |
| Two forms of a gene (`data/column-forms.js`) | coloring a gene and measuring its range, 10M cells | from the sparse form a word fill and its few entries, about 4 ms; from the column two passes over every cell on the main thread, about 50 ms |
| Read activity reports (`data/read-activity-log.js`) | a run over every gene reports thousands of operations | describing, freezing and announcing each as it arrived cost the page 83-129 ms a run; they are taken in at once and acted on once a frame |
| Opening cost (`data/opening-cost.js`) | Chromium, a 20M-cell store | raising 10M cells to 20M in place took 340-380 ms to its first frame, 34-38 ns an added cell, against an opening's 9 ns a cell |
| An opening's fixed cost | Chromium, first frame by cells | 10K cells 27 ms, 100K 23 ms, 1M 31 ms, 5M 70 ms: about 25 ms and 9 ns a cell, hence the measure is the slope of the least-squares line, not the ratio (which would make a cell of the 10K-cell opening cost 300 times one of the 5M) |
| Decoded chunks kept across the steps of a progressive opening (`anndata/decoded-chunks.js`) | HLCA's 1% blocks share | at its last step it decoded the 182 chunks of `X_umap` again (7.8 MB, mostly page-cache hits, 575 ms) |

### Pricing the gene sample and the first gene

| What | Conditions | Result |
|---|---|---|
| The first gene's time target (about 20 s) | a 20% share of a 2.3M-cell CSR file on a 1.2 MB/s link | 20-25 minutes for its first gene, which looks endless |
| A price made before the link was measured | CELLxGENE atlases over S3 | priced 20 s and read in 6-15 s: it rests on early, small reads that measure a fraction of the link's rate, which is why the gene sample is priced to fill half the gene time (`FIRST_GENE_SHARE`) |
| Rows whose elements lie in chunks already counted are passed over (`rowChunks`) | pricing a 50% uniform share's first gene at each step of the gene sample's fit, 20 ms link | 5M single-row spans visited a step took 2 s of an opening |
| Gap units priced as fetched bytes (`gapUnits`) | a 1% uniform share of a 10M-cell `.h5ad` over a 20 ms link | 1,365 rows priced at 81 MB without them were read at 142 MB |
| The fit prices up to every row about log2(n) times | 50M rows; every cell of 10M in blocks order | priced a row at a time 50M rows were 11 s of an opening; priced as runs of rows, every cell of 10M in blocks order was still 0.5 s, so a block order is priced as runs, never a step a row |
| The fit's counts read only their own cells' rows | a 50% uniform share of 10M cells | passing over every open cell in row order for each of the fit's two dozen counts took 0.47 s of an opening |
| First-gene cost includes the chunk index | an 11.4M-cell CELLxGENE atlas | 311 of the first gene's 539 requests were chunk B-tree nodes, in five dependent steps |
| Values read per gene by range only from a raw data array | a 200,000-cell gzip CSR file at 4 MB/s | a compressed array is read in whole chunks and a gene's entries lie in nearly every chunk of the rows read, so the first gene read 146 MB, twice the indices alone |
| Block size on a first visit | an 11.4M-cell CELLxGENE atlas | 0.8-1 MB/s at its first resolve (its delivery taken for a full link, slow-start connections) and 11.5 MB/s on its first gene, so the rate used is the higher of the measured one and the seed (the host's remembered rate or the page's downlink); an underestimate halves the cells per halving: the embryo atlas in blocks of 20 read 5,232 cells in the 20 s gene time, in blocks of 80 11,543 |
| A revisit keeps its block rung | CELLxGENE atlases reopened | an atlas reopened in blocks of 80 after a first visit in blocks of 20 read its first gene again (15 MB, 311 requests); the HLCA at 76 → 19 read 9.5 MB in 512 requests, more than the first visit |
| Sample representativeness against block size | CELLxGENE; category shares of the sample against the whole, total variation distance | on files whose rows are mixed (the HLCA core and full) the block size changes nothing at a given count of cells, so more cells is pure gain; on a file sorted by donor (an 11.4M-cell mouse embryo atlas) larger blocks at the same count are coarser (donor 0.163 at 20-cell blocks, 0.402 at 160, for 2,000 cells), but at the same count of blocks (the same round trips) they are as good or better: 2,000 cells in blocks of 20 gave cell type 0.070, donor 0.163; 20,000 cells in blocks of 160 gave 0.050, 0.134 |
| The ceiling on a share's sample (`wholeSampleCost`) | a 1% uniform share of a 10M-cell CSR `.h5ad` fitted to the same time | read 165-172 MB for its first gene against every cell's sample 150 MB; in the file's cell order a 1% share of Tabula Sapiens (CELLxGENE, S3) read 3.4 MB in 434 requests for 83 cells where every cell's sample made 351 (its cells reach past every cell's blocks while its bytes are fewer) |
| A read priced at the link's round trip under load, not its shortest (`loadedRoundTripMs`) | an `.h5ad` export of 21,285 HLCA cells over S3 | priced at the shortest round trip (230 ms) it was 198 s and took 271-386 s (0.36-0.5 s a request): a read of many requests at once comes out far below its time ({doc}`22_decisions_link_model_and_request_planning`) |
| A source on this machine priced as the machine at rest | a 10M-cell store on a loopback host, Chromium, openings while the machine is busy (another browser starting, a run ending beside it) | a slow opening in 5 of 6: `slowRead` counted as transfer (100-13,000 B/ms measured as the link's rate), `slowTurns` small reads of which most waited (a round trip under load of 863 ms where 1.5-7 ms is usual), the first gene's many waves of requests priced at 45.7 s; at rest the opening's 0.9 MB of metadata read at 45 MB/s, a lower bound |
| A revisit priced as free | the HLCA core's 585,000-cell revisit (Census S3, kept index of 1,518 cells, nothing measured) | took the gene index's memory fit, 102,271 cells: its first gene read 575 MB in 14,680 requests over 148 s, where the first visit read 9.4 MB in 3.7 s |
| A later visit answers from what it kept | Chromium, the HLCA core | a revisit priced at the remembered rate read 5,400 new cells with its first gene, 6.6 s and 26.3 MB, where the first visit's took 2.9 s: so the sample is the kept cells and the priced count is returned as `growTo` |
| The pass priced over the rows in open order | a 50% share of 5M cells | the runs of every row sorted (3 s) and the kept values counted a row at a time (1.3 s); the gene's plan counts the work by comparison sorts and row-pointer lookups so these stay out |
| A 10% blocks share of a priority-ordered file | a 10M-cell `.h5ad` | read 380 MB of the 385 MB every cell reads, so the file's own order is offered |
| A partly kept file priced chunk by chunk | a GitHub Zarr folder whose first visit kept 65,536 of its 561,947 cells | priced at the mean, every chunk at 88% of its bytes, a revisit opened 750 cells |
| `costOf` walks the share's cells | a 50% share of 10M cells | 60 ms of an opening, which is why a share over a link whose round trips cost anything is decided to open whole before its first view is priced |
| Growing a gene sample | a sample of 2,192 of 2,282,447 cells | the first step grew it by 71,327 rows and took five minutes, so the first step of a growth adds about as many rows as are read (doubling), never a step of a share of every row |

### Gene indexes, kept copies and cell subsets

| What | Conditions | Result |
|---|---|---|
| Kept copies read in long ranges (`KEPT_READ_GAP_BYTES`, 1 MiB; `sparse-column-index.js`) | the origin-private file system, 60,000 reads of a gene's piece of a 10M-cell index, Chromium, Apple M1 Pro | a read costs about a millisecond before its first byte, in which a disk delivers about a mebibyte (70 s for the 60,000 reads); read one gene at a time the 10M-cell revisit's DE took 71.7 s against 16.0 s from the index built in memory |
| Partial chunks: arrays read byte by byte (stored in one piece, or an unfiltered HDF5 chunk) | 1,000 cells of a 1M-cell CSR `.h5ad` over 30 ms / 5 MB/s | 4.5 MB in 12 s, not 121 MB in 25 s |
| Analyses after a share lowered and raised back in place | a 10M-cell Zarr, every cell, Chromium; at 1M | differential expression 9.4 s on the first opening and 35.9 s after the round trip (1M: 1.5 s and 9.0 s); without a gene read before the round trip the run failed after 3 min 50 s (a read of 8,610 genes asked every gene's column at once); the cause was the share's cell order, the gene index in segments and the read priced by route |
| A gene that waits makes pass 1 read every slab's values beside its indices | values stored in chunks read whole, a CELLxGENE atlas | a gene sample's third gene needed 316 requests when each later gene read the chunks its entries lie in |
| `_valuesBeside` | a 40,000-cell file over 150 ms at 5 MB/s | a 1% subsample's first gene read pointers, indices, then values: three round trips where two do; reading every slab with its values, a 40,000-cell CSR `.h5ad`'s first gene over every cell read 35.4 MB where 28.4 MB do (hence a slab is read with its values only when the indices still to read take less than a round trip) |
| A gene waits for every unbuilt segment's read together (`segmentedColumns`) | when their indexes fit the memory budget; a 1% subsample on a 150 ms link | its first gene took 2 × its floor reading one segment at a time |
| Cost of a cell subset by genes against by rows (`sourceWeight`) | two large groups of a 1M-cell CSR `.h5ad` over 5 MB/s | read from the built index in memory, not as 400 MB of rows again |
| A gene's kept entries alone | 100 genes of a 1% subsample, 150 ms at 5 MB/s | reading the values of the kept entries' runs (one round trip more) took 3.2 s against 2.7 s for 13% fewer bytes, so it is read only when the value chunks it would skip cost more than the round trip |
| `subsetCells` orders the cells by row with a two-pass 16-bit radix sort | a million cells | a comparison sort took a quarter of a second |
| Keeping a run's entries on a decode worker | every gene of the rebuilt Suo in its two largest cell types, a local file in Chromium | 305-326 ms with about 150 ms of it keeping entries on the reader thread, 206-215 ms with the workers keeping them, against 150 ms for reading the whole of `X` |
| Nested `uns` entries are left to the reading that needs them | a CELLxGENE file's objects | about half of its objects are the children of nested `uns` entries (the HLCA's: 259), so opening loads an AnnData group to depth 2 |
| A dense expression matrix is not indexed ahead of need (`reader/matrices.js` `prepareMatrix`) | one dataset of 28,317 x 3,136 float32 values, read from disk | its index holds every value: 663 MB resident when built at open, for a dataset whose genes may never be shown; built when the first gene is asked for it took 549 ms |
| A categorical obs field of one category a cell is coded by the categories its open cells hold (`reader/fields.js` `_categoryCoding`) | the barcodes of CELLxGENE's soundlife atlas kept as a categorical | 1.2M categories, 35 MB of dictionary: read whole for the few a share of the cells names |
| 64-bit integer arrays through their 32-bit words (`anndata/values.js`) | a 10x `.h5`'s int64 indices (44 million at 100,000 cells) | 0.5 s through a BigInt per element's comparisons and `Number()`, a few tens of milliseconds through words |

### Kept stores, local files and opening other formats

| What | Conditions | Result |
|---|---|---|
| The kept-directory lock (`transcode/space.js`) | WebKit, a file a worker rewrites through a synchronous access handle while the page reads it | 5 of 2,645 reads refused and logged without the lock, none of 3,001 under it; Chromium and Firefox refuse none either way |
| The `names` attribute of an R matrix's numeric payload is never read (`transcode/r-object.js`) | a 100,000-cell SingleCellExperiment | arithmetic on `@x` in R leaves one name per stored value on it, and decoding those strings cost 7 of the 10.7 s pass |
| Stage size of a transcode payload write (`transcode/entries.js`) | an origin-private file | a megabyte a write at 2.6 times the rate of 64 KB writes in Chromium (670 against 260 MB/s; Firefox 1.7 against 1.2 GB/s) |
| A text table of cells read in one pass without keeping its bytes (`transcode/tables.js`) | a 2M-row table of 30 columns | 0.6 GB less than when its bytes were kept for a second pass |
| Transcoded `X` kept by feature (csc with prefix counts; `transcode/anndata-model.js`) | a 50,000-cell Seurat object (75M entries a layer), Node, Apple M1 Pro | laying `X` out by feature takes 4 s more per matrix (6.7 s with none, 10.9 s with `X`, 14.8 s with `X` and counts), and every later visit reads a gene in 7 ms instead of building a column index first (1.2 s and 0.9 s for the first two genes, 800 MB of memory, on every visit) |
| Dropped folder listing in Chromium (`app/ui/modules/dataset-drop.js`) | real drops, medians of five | through the directory handle faster than through the drop's entries (whose calls the file system answers one at a time however many are asked at once): 22,023 files in 2.7 s against 4.4 s; `suo.zarr`'s 6,272 files in 0.63 s against 1.2 s; so the handle is read where the drop gave one (other engines: {doc}`22_decisions_link_model_and_request_planning`) |
| A `.zarr.zip` stored uncompressed whose stored chunks hold every cell | read as elements | every element past the first chunk's index looked held, the first gene was priced at a few bytes and the atlas opened with no sample |
| The default sample's absent object | `uns['cellucid']['priority']['seed']`, a Zarr chunk holding only the fill value (a web host answers 404), 2026-11-10 | asked for 11 times at start, once per growth step; a web host store now asks once and takes the absence as known on a re-read, another opening of the same version and a reader with no page cache |
| HDF5 global heap collections read for the bytes their objects need (`formats/hdf5/heaps.js`) | CELLxGENE's HLCA over S3, 4 KB collections | padding them to a round trip's worth of the link fetched 3.4 MB of the 20 MB an opening read |
| The HDF5 chunk B-tree read as reads need it (`formats/hdf5/chunk-index.js`) | the parents of the leaves of an 11.4M-cell atlas's `X/indices` (2.26M chunks), S3 | 552 nodes (the leaves' parents: about 35,000 leaves of 2.26M chunks at 64 entries a node) and 13.5 s of requests before the first view when read with the table (562 is every internal node of the tree: those 552, the 9 above them and the root; the comments that gave 552 and 562 were never asserted by a test); a gene sample's 170 blocks need about 170 of them, and the parents of tens of thousands of leaves are why a large tree stops at the first level that is not one wave |
| Chunk-index leaves read ahead as the read's own requests, not as background work | a CSR matrix's first gene over 150 ms at 5 MB/s; a burst of 100 Loom genes | leaves asked for behind the read's chunks came only after them and every slab of a scan waited a round trip more: a 10 MB scan took 3.1 s against a floor of 2.2 s; the burst made 44 leaf requests, then each read's few chunks alone, about 36 requests in seven turns |

## Tried and reverted

| What | When | Why it went | What would change it |
| --- | --- | --- | --- |
| HDF5 string heaps loaded as one union after every piece | 2026-10-28 | 10M at 150 ms: 86.2 → 68.8 s, but 100M at 20 ms HTTP/2 got slower (72.0 → 87.7 s): a barrier on a fast link | — |
| Fixed batches of a round of the link's requests | 2026-10-28 | 100M HTTP/2 87.1 s, against 71.8 s with no batching | replaced 2026-10-29 by the room rule, now in `HeapStream`: load when the link has a free slot, gather while requests queue (10M 150 ms 77.0 s; 20 ms HTTP/2 within 2% of no batching) |
| Never fetching again any byte inside a gap | 2026-10-30 | no duplicate bytes, but more requests (5 → 8 in the HDF5 gene pricing test, about 1,000 more at 10M over HTTP/1.1) | — |
| A request budget per heap run for obs names | 2026-10-30 | 0.1% share 4.95 → 7.69 s and 241 → 376 MB | — |
| Asking held heap bytes at once; a window bounded by file span | 2026-10-30 | every cell 436 → 566 MB; no gain | — |
| A whole heap collection asked up front when its tail was expected to hold a kept name | 2026-10-30 | 15-25% fewer requests, 0.5-16 MB more bytes | knowing a collection's size before its first read |
| A uniform share's gene-index segments merged into one in open order | 2026-10-30 | 1,019 MB against 517 MB | the store-order segment |
| Dropping the decoded-chunk cache during scans | 2026-10-30 | one chunk read again at every slab boundary | — |
| Searching each kept cell among the runs of the stored rows when a CSC gene is placed on a subsample's cells (`cell-maps.js` `keptRowsOf`) | by 2026-11-12 | a 50% uniform share of 10M cells, its rows sorted and each cell searched among 2.5M runs: 1.6 s of a CSC gene's 1.8 s (machine and engine not named) | a radix order of the kept rows (`orderByKey`), each row's cell where the order puts it |
| Walking every cell for every gene when a cell subset is made of sparse gene columns (`cell-subset-blocks.js` `cellSubsetOfSparseColumns`) | by 2026-11-12 | a differential expression of two groups of a 13M-cell gene sample of a 50M-cell file: 5.8 ms a gene whatever the gene held, 174 s for 30,000 genes (machine and engine not named) | galloping each gene's ascending entries through the ascending rows: the entries times log of cells over entries, not one step a cell |
| A single versioned file (a Census `.h5ad`) taken over by version at a growth step | 2026-10-30 | the gene list hung in Node; now done for folders only | the hang traced |
| A block size from the merge gap measured on one read | before 2026-10-30 | the same `.h5ad` gave blocks of 28 to 66 cells from one opening to the next | — (`BYTE_READ_BLOCK_BYTES` is a constant) |

## Known limits and open items

- **CPU on fast links:** a 50% share's first gene reads through every
  stored row's index: at 4M and 5M cells, on a near-zero-latency
  link and at 20 ms / 50 MB/s, 1.00–1.04× every cell's in a loaded machine's
  run-to-run noise (same bytes, same requests, the same pass), where its own
  pass took 1.4×. A 50% share's edges still place them on its cells after
  the last line ({doc}`24_decisions_graph_edges_at_scale`): about 1.03×
  every cell's at 5M. The 10M emulated S3 1% share took 1.015× every
  cell's time on identical reads (its kept rows' pointers are read before
  its pass, a step every cell's pipelined pass does not wait for).
- **A half share's marker genes:** a 50% uniform share of a 10M-cell store
  finds marker genes in 7.1-8.4 s (`.h5ad` CSR) and 7.2 s (Zarr v3 sharded)
  against 3.8-4.1 s and 3.0 s for every cell after its first gene (Chromium,
  same hour, 2026-10-31); its first gene reads in 1.08-1.10x (CSR) and
  0.84-0.88x (Zarr) every cell's, values exact. The rest is each gene's
  entries put in the share's open order, on the reader worker alone: traced,
  4.5 s of its 7.5 s (a radix sort of each gene's entries by open row 2.2 s,
  the lookup of each entry's open row 1.3 s, the segments' copy 0.6 s), work
  every cell's genes, already in open order, never do. Placement over every
  stored row's index cost more (1.3-1.7x, see Decisions), so a half share
  takes its own pass. Removing it needs analyses that take a gene's entries
  in stored order with each worker holding the groups by stored row, or
  spread shares opened in store order: not built.
- **Share round trips (fixed 2026-10-30, browser re-measure pending):** a
  raise in place after a lowered share used to fit the gene sample again and
  opened genes for every cell (50M: 4.49M → every cell; markers after a round
  trip 62-81 s against 15 s without one; the 50M DE allocation failure in
  {doc}`27_decisions_memory_model`). The same cells now come back with their
  sample (see Decisions), verified in Node only; the 10M and 50M browser runs
  of the fixed code are not taken yet. Still fitted over the link as
  measured: a share never opened before (the decode rate measured on small
  reads under-prices it), and the earlier opening's in-memory gene index is
  not carried through the lowered share (the kept index serves it).
- **A gene sample on this computer with every core busy:** a source on this
  computer is priced as the computer at rest (see Decisions: a disk's rate,
  its shortest round trip), and openings right after other runs, or with four
  CPU-bound processes beside the browser, read the first gene for every cell
  (10 of 10, Chromium, 2026-10-31). With ten CPU-bound processes on the ten
  cores, 1 of 10 openings (a 10M-cell `.h5ad` on a loopback host) measured
  even its shortest round trip at 900 ms (1.5-12 ms otherwise), priced every
  cell's first gene at 47.7 s and read it for 1,610,169 of 10M cells. Pricing
  this computer's round trip from a prior (as its rate is from a disk's)
  would remove it; it changes the progressive opening's and the block rung's
  prices too, so it is left for a decision. The same measurement is what
  failed browser specs on loaded runners (2026-11-06): one request of a
  238 KB file answered 370 ms late prices the first gene of 1,500 of its
  cells at 36 round trips, though the reader holds every byte of the file
  from that one request, and reads genes for a sample of them
  (`tests/browser/held-machine.spec.mjs`). The specs hold the link they plan
  on ({doc}`14_testing_ci_and_release_process`, `utils/machine.js`); the
  app's pricing is as it was. Two cheaper remedies than a prior were seen and
  not made: a read whose bytes the reader already holds could be priced as no
  request; and a share's gene index kept under its own order could be looked
  for before a sample is fitted (today only the sample's cell order is, so a
  later visit that reads genes for a sample -- a slower link, a busier
  machine -- indexes the matrix again instead of reading the index an earlier
  visit kept for every cell: shown in Chromium with the decoding rate held
  slow on the later visit only).
- **A later run of a measured analysis kind takes up to about two and a half
  times as long to its final result:** once an analysis kind's cost is
  measured (`progressive-cells.js` `AnalysisCostModel`, kept in memory and in
  this browser's storage), a later run of it (the same page, or a revisit in
  the same browser) runs a first pass on a sample and then the pass on every
  cell. The full pass costs what a single run does, so the time to the final
  result is that plus the sample pass: 10M-cell Zarr store, every cell,
  Chromium, a second run on the same page: differential expression 3.0 s →
  7.3 s, marker genes 3.3 s → 4.2 s. The figures recorded here before
  (markers 16 s and DE 36 s after a share round trip) were mostly not this:
  the harness's `--single-pass` never cleared the costs (it looked the module
  up in the page's resource timing buffer, which a 10M page fills before the
  module loads), and the reads after a round trip were slow of their own
  (fixed: see Decisions, 2026-10-31). To do: the full pass reuses the sample
  pass's work, or the sample pass is skipped when the full pass is priced
  short.
- Per-read O(cells) work remains in DE at 10M (the distinct-cells check and
  the DE cells' list copied to the reader on each read).
- The top-up prices the requests of scattered rows about 2× too low: a TS 1%
  share's top-up made 1,487 requests against 847 priced, 1.19× every cell's.
  Absolute request prices for Zarr folders are low as well.
- A smaller share does not reuse a larger share's kept store-order run on a
  later visit; it indexes again. Applying 1% → 10% of a share in a cell
  order is a new reading, not a growth in place.
- Some limits come from the file's layout. A gene-major (CSC) file reads a
  gene's whole column at any share. String heaps and uniform shares of
  row-chunked arrays reach every chunk, so their bytes stay near the whole.
- Memory at a 50% uniform share stays above every cell's (5M: 711 against
  599 MB peak buffers) because row pointers are held several times. It is
  counted and released.
- A host naming only a tag made from the bytes shows a rewrite of a store
  only when the consolidated metadata's bytes change: one keeping every shape
  and attribute and changing only values is not seen.
- **CSC `.h5ad` over a slow link (40,000 cells, 150 ms, 5 MB/s):** 100 genes
  read in 2.5–3.3 s and 81–122 requests from one run to the next (the plan
  follows the link measured from the opening's few reads); the DE read
  1.8–2.8 s.
- **Spread shares in every cell's block order:** in the browser proofs
  (Chromium and Firefox; Tabula Sapiens and HLCA over S3 at 1% and 10%
  uniform and every cell, two visits each, no page error) a share's first
  picture costs what every cell's does (the opening's metadata) and once
  more: Tabula Sapiens 1 % 12.4 s against 10.9 s for every cell (Chromium),
  HLCA 1 % 17.1 s against 14.7 s (Firefox); single runs over S3, not
  separated from its noise. Revisits draw every share in 1.2–1.9 s. The
  harness's DE read at 1% reported 2–4 cells read (the gene sample's cells in
  the two groups), against 138–406 for every cell. A WebKit crash on the
  second visit of a 1% uniform share is not yet traced.
- **WebKit:** every app run records the page error "Fetch API cannot load
  http: /127.0.0.1:<port>/ due to access control checks" (the harness's page
  port); not traced to the app or the harness.
- **A background top-up on a closed reader** (after Apply) logs "this._root
  is null"; nothing is affected.
- **Kept data removed once:** kept data is stamped with a code fingerprint,
  so everything kept by a version before the fingerprints was removed on its
  first visit.
- Not built yet: 16-bit gene-index positions, an early partial picture within
  the first gene's read, keeping a cut-short conversion or text-matrix pass,
  and genes on the cells read so far during a pass.
- Planned clean break: rename the wire names `{mode: 'kept'}` and `cellOrder`
  and write the shared fixtures again in the web, Python and R repositories
  (`node scripts/sync-shared-files.mjs`).

## Where the code and tests live

| Concern | Module(s) under `cellucid/assets/js/` | Tests in `cellucid/tests/` |
| --- | --- | --- |
| Share-cost ceiling, auto-fit, ladder, gene sample | `data/anndata/subsample-plan.js` | `subsample-share-cost-contract`, `share-store-order-pass-contract`, `subsample-auto-fit`, `gene-sample-contract` |
| Subsampling methods, rank cache | `data/subsample.js` | `subsample-conformance`, `subsample-growth-cache` |
| Requests and modes | `data/cell-subset-policy.js` | `subsample-address-and-control`, `subsample-recorded-set-contract` |
| Subsample view | `data/anndata/views/cell-subset.js` | `subsample-gather-differential`, `open-cells-rows-contract`, `efficiency-floor` |
| Gene index plans, the first gene's cost | `data/anndata/matrix-access.js`, `data/gene-read-plan.js`, `data/sparse-column-index.js` | `anndata-sparse-columns-contract` |
| Segments, store order, every stored row's index restricted to a share (`wholeIsCheaper`, `PassRates`) | `data/anndata/prefix-column-index.js`, `data/anndata/reader/gene-segments.js` | `subsample-share-cost-contract`, `share-restricted-whole-contract`, `anndata-sparse-columns-contract`, `subsample-gene-read-memory` |
| Two forms of a gene, gene delivery | `data/column-forms.js`, `data/gene-delivery.js`, `app/analysis/data/gene-route.js` | `gene-delivery-contract`, `expression-transform` |
| Kept indexes, code fingerprints | `data/anndata/kept-column-index.js`, `data/anndata/transcode/space.js`, `data/anndata/transcode/kept-format.js`, `data/anndata/transcode/kept-journal.js`, `data/anndata/transcode/kept-code-fingerprints.js` | `kept-column-index-contract`, `transcode-cache-contract`, `kept-parts-contract`, `kept-code-fingerprints-contract` |
| Page gene sample, top-up | `app/state/managers/field/gene-sample.js` | `gene-sample-contract`, browser `gene-sample.spec` |
| Genes read for (steps, estimate, Read / Stop, the links to it), the Session sections | `app/ui/modules/gene-read-control.js`, `app/ui/modules/session-sections.js` | `gene-read-control-contract`, browser `gene-sample.spec` |
| Reading again, growth in place, progressive opening | `data/anndata-reader-worker.js`, `data/anndata-file-source.js`, `data/open-cells-work.js`, `app/open-cells-resize.js`, `app/state/managers/open-cells.js`, `data/anndata/progressive-open.js`, `data/anndata/decoded-chunks.js`, `data/opening-cost.js` | `open-cells-resize-state`, `progressive-subsample`, browser `progressive-subsample.spec`, `open-cells-folder-reuse` |
| Strings, categories | `formats/hdf5/nodes.js`, `formats/io/byte-reader.js`, `data/anndata/reader/fields.js`, `data/anndata/reader/names.js` | `obs-strings-share-contract`, `present-categories-contract`, `hdf5-reader-structures` |
| Edges at a share | `data/anndata/graph-lines.js` | `edges-at-a-share-contract` |
| Cell-subset reads (run batching, genes per task, checked and kept cells, after a share round trip), sorting keys | `app/analysis/shared/gene-stream-runner.js`, `app/analysis/shared/run-tuner.js`, `data/anndata/reader/cell-subset-reads.js` (`getMatrixCellSubset`), `data/anndata/matrix-access.js` (`touchedChunkBytes`), `data/anndata/reader/cell-subset-reads.js` (`_subsetInStoreOrder`), `data/anndata/subset-entries.js`, `data/anndata/prefix-column-index.js`, `utils/radix-order.js` | `round-trip-analyses-contract`, `analysis-gene-reads-bounded-contract`, `analysis-run-tuner-contract` |
| Text-matrix streams, one pass a progressive opening | `data/anndata/transcode/matrix-market.js`, `data/anndata-reader-worker.js` (`ensureMatrixShare`, `swapOntoShare`), `data/anndata/reader/gene-sample.js` (`shareRows`) | `text-matrix-stream-contract`, browser `text-matrix-stream.spec` |

`cellucid/scripts/efficiency-floor/share-grid.mjs` measures every downstream
at 1, 10 and 50% shares against every cell, in the whole app.

## Related pages

- {doc}`20_read_path_architecture`, {doc}`19_subsampling_specification`,
  {doc}`09_data_loading_pipeline_and_caching`,
  {doc}`07_rendering_pipeline_webgl_and_performance_notes`,
  {doc}`01_codebase_map_and_entry_points`
- {doc}`22_decisions_link_model_and_request_planning` (how a read is priced),
  {doc}`28_decisions_rendering_layers` (how unread cells look),
  {doc}`24_decisions_graph_edges_at_scale`, {doc}`27_decisions_memory_model`,
  {doc}`21_decisions_decoders_and_codecs`
- User view: {doc}`../b_data_loading/21_subsampled_loading`
