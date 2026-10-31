# Decision record: partial reading

This page answers: **"Why does Cellucid read only part of a large dataset the
way it does, what was measured to justify it, and what was tried and taken
back?"** Partial reading (a share of the cells, a gene read for a sample of
the open cells, kept data reused on a later visit) is the main way the web app
gets fast on large files, judged first on public links. The modules and
invariants are mapped in {doc}`20_read_path_architecture`; which cells a
subsample keeps is {doc}`19_subsampling_specification`. This page holds the
reasons and numbers behind them. Update it in the same change as the code it
describes.

## At a glance

**Audience**
- Developers changing the share planner, the gene sample,
  kept gene indexes or progressive opening.

**Time**
- 20 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture`, {doc}`19_subsampling_specification`.

---

## What it does now

All paths are under `cellucid/assets/js/`.

1. **The share-cost rule.** A share of the cells never costs more than a
   larger share of the same request, or than every cell. About s/n of every
   cell's cost is the aim, not a rule. In code, the "never more than every
   cell" half is the ceiling `wholeSampleCost` in
   `data/anndata/subsample-plan.js`: the bytes and requests every cell's own
   gene sample would make in the same time. A share's gene sample
   (`geneTimeFit`, `ceiling`) and the page's top-up
   (`AnnDataReader.geneSampleCost`, `geneTime`) stay within it.
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
   a gene index kept in this browser was built under, otherwise the lowest
   rung where a block fills its round trips, at the higher of the measured
   rate and a seed rate (the host's remembered rate, or this page's own
   downlink). The chosen rung is recorded with the cells (`blockSize`,
   `cellOrder`). The formulas are in {doc}`19_subsampling_specification`.
5. **The gene sample.** For a cell-major root matrix, a gene is read for the
   first `cells` open cells (`geneSampleOf`). That count is the smaller of two
   fits. `geneTimeFit` fits the first gene into `FIRST_GENE_SHARE` (½) of
   `GENE_READ_TARGET_MS` (20 s) over the measured link, chunk-index reads
   included. `geneIndexFit` fits the gene index into the memory the open cells
   leave. The fit's budget is the whole 20 s. Open order decides what the
   sample reads:
   - Every cell of a file in its own order opens in the priority order of the
     file's default `blocks` subsample (`cellOrder`), so the first `m` cells
     are whole storage blocks.
   - A `uniform`, `stratified` or `spatial` share of such a file opens in that
     same order with its own cells filtered in (`shareInCellOrder`). Its sample
     is therefore a subset of every cell's reads.
   - A `blocks` share keeps its own order. A share of a priority-ordered file
     opens in file order (`placeInFile`).
   - Kept runs ending within the sample cost nothing (`breaks`).
   - When the same cells were open before (a share lowered and raised back),
     they reuse the order and sample they were read in (`seenGeneSampleOf`).
     They are not fitted again. Each opening of the same bytes carries the
     cells earlier openings opened, with their order and sample as they were
     when replaced (`AnnDataReader._geneSamplesSeen`, noted by
     `noteGeneSampleSeen`, keyed by every cell or by the share's spec); a
     sample grown by the page's top-up keeps the fit that sized it (the page
     passes `topUp` to `growGeneSample`), and a share chosen while every cell
     was open is held to every cell's sample as it was read (`wholeSampleCost`
     `seen`).
   - A gene-major matrix (CSC, a Cellucid-ready file) reads each gene for
     every open cell and has no sample.
6. **The top-up on the page.** After the first gene,
   `app/state/managers/field/gene-sample.js` `topUpGeneSample` spends the rest
   of the 20 s reading genes for more cells. It sizes steps from the reads it
   measured (a fixed time plus a time per cell), within the ceiling and the
   index memory (`index.fits`). It stops when a step would add less than a
   quarter of the cells already read. A gene that is asked for stops the step
   at once (`yieldGeneSampleTopUp`). An analysis holds the sample still before
   it runs (`holdGeneSample`). The legend's *Read for all cells*
   (`readGenesForMoreCells`) grows the sample in doubling steps
   (`growthSegmentBounds`), each step reading only the rows it adds.
7. **Unread cells** hold NaN, never 0, and are drawn faded beneath the read
   cells. The look is decided in {doc}`28_decisions_rendering_layers`.
8. **The subsample view** (`data/anndata/views/cell-subset.js`):
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
     pieces listed or ordered: a 50 % uniform share of a 10M-cell `.h5ad`
     took 4.7-7.0 s for its first view against every cell's 0.9 s, and
     takes 0.9-1.0 s. Runs are ordered by source row with a radix order
     (`rankedOrder`), never a comparison sort.
9. **Gene indexes in segments.** Rows in priority order are indexed in
   doubling segments (`prefix-column-index.js`; the first is 1/32 of the rows,
   at least a slab). A later opening of the same method and seed takes over
   the segments it shares (`carriedSegments`). Rows spread through the store
   (`rowsInterleaved`) are one segment in **store order** (`reader.js`
   `_segmentedColumns`, `_storeOrderMatrix`), answered in open order through
   `fileOrderedAccess`. That segment reads each stored chunk once. When those
   rows reach to the last open row, lie in every chunk (the whole reads no
   more bytes: `_readsNoMoreThanRows`) and their own pass is more work than
   every cell's (`wholeIsCheaper`), the segment is instead **every stored
   row's index restricted to them** (`_restrictsWhole`,
   `restrictedAccess`): every cell's pass over the stored matrix (the same
   reads), with no kept row's pointer read and no entry gathered, kept in
   this browser as every cell's (either opening takes the other's). The work
   is weighed per entry over the index built and each of its genes answered
   once, at rates measured on the source as passes and genes run
   (`PassRates`: the index's own work per entry, `onWork`; the view's own,
   its `ViewMeter`; each path's placement of a gene in open order,
   `placing('whole')` per entry of every stored row's, `placing('own')` per
   entry of the share's), carried to every opening of the same store; before
   any is measured, a share keeping at least two thirds of the entries reads
   the whole. Placing a gene from every stored row's entries tests each
   against a bit of the segment's open rows (`openRowBits`, 1.25 MB at 10M
   stored rows, held in cache) and looks up only those of open rows.
10. **Kept gene indexes and revisits** (`kept-column-index.js`). Each built run
    is kept in the origin-private file system, filed by file identity, matrix,
    shape and order key (method, seed and settings, never size). It is checked
    before use: manifest, file rows, the file's row pointers read now, and a
    CRC-32 per gene piece. A store-order run is used whole or not at all. A
    run whose values were all read is kept with them (`sparse` layout), so a
    revisit reads no value again. A run used is read as an index built in
    the opening is held: whole, in one read of each array at its first gene,
    when the gene-index memory budget holds it (counted in that share), then
    every gene from memory; otherwise the genes asked for together are read
    as one range per run of neighbouring genes, never one read per gene (a
    10M-cell revisit's DE read each gene's piece on its own: 71.7 s against
    16.0 s from an index built in memory). A run whose reading is priced
    above building it again (`rebuildMs`, at the read rate measured in the
    thread) is not used. A revisit's gene sample (`geneTimeFit`) is the
    kept prefix (the kept runs ending within the count the target prices):
    its first gene reads no new cell, so coming back is never slower than
    the first visit. The count priced at the link as measured or
    remembered (`pricingLink`) is `growTo`, the page's top-up's first step
    (`topUpGeneSample`, in the background, yielding to genes asked for);
    later steps are sized by the reads. Over a link with neither (not read
    ahead, `linkPriced` false: an opening served wholly from kept data, host
    speed forgotten) a new cell has no price: `neededMs` null, no
    `growTo`, and the first step is a quarter of the sample. Measured on
    HLCA core (Census S3): priced as free, a Node revisit refitted 1,518
    kept cells to 102,271 and read its first gene in 147.8 s, 575 MB,
    14,680 requests; priced at the remembered rate, a Chromium revisit read
    5,400-5,800 new cells with it in 6.6-7.9 s where the first visit's took
    2.9-3.1 s. `transcode/space.js` holds
    `KEPT_INDEX_REVISION`, `KEPT_STORE_REVISION` and `KEPT_BYTES_REVISION`.
    Raising one removes that kind of kept data from every user's browser.
11. **Progressive opening** (`progressive-open.js`). When the first view would
    take longer than about a second, the app opens the most cells predicted
    within that target (at least √n). If that view would not arrive at least
    twice as soon as every cell, the opening is not progressive. Otherwise it
    doubles to the cells asked for, grown in place by `app/open-cells-resize.js`:
    - Nested orders grow by prefix.
    - A share in a cell order, or on a priority-ordered file, grows by `prefix`
      requests (`prefixStep`).
    - A file whose embedding chunks are much coarser than its blocks grows
      chunk by chunk (`chunkMajorOrder`).
    - A share grows only when its cells are all on the map no later than one
      reading of them would put them (`grownLater`, `chunkedPlan`; the
      owner's rule judged on completion, not on the first picture). A
      growth reads the same bytes in more readings, a loaded round trip
      each (and from chunks shows every cell again in its order), so over
      any priced link a share opens whole; it grows only on the page's
      measured cost of showing cells, with no link to pay. Measured, 50 %
      uniform of 10M cells, local web host: grown from chunks all 5M at
      4.5 s, whole 2.8 s.
    - Every cell of a file in its own order does not grow. It opens in its
      default subsample's order (for a gene sample) only when the sample
      brings the first gene sooner by more than the order costs the
      opening end to end: made and timed at opening, every pass over the
      open cells it adds priced at that time (`subsample-plan.js`
      `orderOpeningMs`, `ORDER_OPENING_PASSES`, `orderPays`); else in file
      order, its genes read for every cell.
12. **Incremental growth.** One rank pass per opening: `createSubsampleCache`
    (`data/subsample.js`) is held by the reader through growth and Apply.
    Uniform and blocks steps are prefixes of one ranking.
13. **Strings and categories at a share.**
    - String columns read only a share's rows (store contract `readRows`;
      HDF5 `HeapStream`).
    - The byte reader never fetches bytes already in flight. A running read
      keeps its gap bytes until it ends, within the fetched-bytes budget.
    - A categorical column whose dictionary has more categories than open
      cells, and would cost more than one policy request, reads its codes
      first, then only the present labels and colours (`reader.js`
      `_readPresentCategorical`).
14. **The subsampling specification** is shared byte for byte with Python and
    R: {doc}`19_subsampling_specification`.

## Decisions and why

| Date | Decision | Why | Rules out |
| --- | --- | --- | --- |
| 2026-10-27 | `blocks` is the default method, aligned to about four chunks of the dominant array | a run of neighbours is a run of storage: a block share reads about its share of the bytes | a uniform default on a file in its own order (10 % of HLCA's cells touch 45 % of `X`'s chunks) |
| 2026-10-27 | a subsample opens in its priority order | its first `m` cells are a sample at every `m`, so growth and segment indexing append | rereading on growth |
| 2026-10-28 | a share that rounds below two cells is refused with the reason, never made larger | exact cells: sessions, links, Python and R name the same set | silent changes to the cell set |
| 2026-10-28 | large row-stored atlases open **all cells the device can hold**. A gene is read for a stated sample sized to the gene-time target, and the other cells wait until the user chooses *Read for all cells* (its cost shown). Analyses say which cells they used | the map costs bytes a cell; a gene of a CSR file costs each cell's whole row | sizing the map by what a gene costs |
| 2026-10-28 | the gene sample is the first `g` open cells in the opening's priority order. Every cell of a file in its own order opens in its default `blocks` order. A session records `genes` only while a sample is in effect, and one saved with every cell read restores with a fitted sample | the sample is then exactly a `blocks` subsample, and kept indexes of earlier visits count | a silent multi-GB read on restore |
| 2026-10-28 | anything kept is filed under a version (content name, else size plus modification time or a strong ETag), checked before any kept byte is used | a changed file must never be served from what was kept | keeping by size or name alone |
| 2026-10-29 | genes on a sample: unread cells fade, read cells drawn on top and larger ({doc}`28_decisions_rendering_layers`) | unread cells should not look like 0 | the same grey for 0 and not read |
| 2026-10-30 | speed and smaller downloads outrank memory, as long as memory is counted and bounded | the owner's priority | trading time or downloads for memory |
| 2026-10-30 | **the share-cost rule** in every metric (network, memory, rendering, analysis, opening); s/n is the ideal, not a rule | a person choosing a smaller share expects it to be cheaper | any planner path where a share costs more than the whole |
| 2026-10-30 | the rule is judged on the costs a person meets: time, download, memory, rendering and analysis work. The request count is a means | a 1 % share's obs names made 575 requests against every cell's 530, but took 30.6 s against 37.9 s and downloaded 75.7 MB against 77.9 MB | spending bytes to save requests (see Tried and reverted) |
| 2026-10-30 | the first gene is priced to half the gene time; the rest is read by measured top-up steps | the opening's link is a lower bound (CELLxGENE atlases priced at 20 s were read in 6-15 s) | one fixed sample from an early price |
| 2026-10-30 | block size comes from one ladder: the base from the layout, the rung from what is kept and the link, and a revisit keeps its rung | a size taken from the open-time link stayed at the storage block on first visits, and revisits then chose another size and read their cells again | sizes that do not nest |
| 2026-10-30 | spread shares open in every cell's block order, with their own cells filtered in | in its own order each of a 1 % share's first cells lay in a chunk of its own: 71 cells in 10.7 s, against every cell's 1,572,326 in 5.5 s | a share's sample costing more than the whole |
| 2026-10-30 | a 50 % share's slow first gene is fixed by cutting the CPU of the store-order pass, not by opening in file order | file order would make a larger share on a later visit read its kept index again | the file-order design |
| 2026-10-30 | levels of detail rank a cell by its identity under the share's seed | the same cells come first in every order, step and revisit | ranking by open position |
| 2026-10-30 | a read on the gene sample checks its row pointers once a run: the most first rows known covered are kept per matrix (pointers held only grow) | markers at 50M read whole columns a few genes at a time (about 3,750 reads), and each read walked the sample's rows twice to ask whether its pointers were held: about 290 s on the reader alone | a per-read coverage walk |
| 2026-10-30 | a cell-subset run (differential expression) is priced over the run, not read by read: the gene stream passes the run's remaining genes (`runGenes`) down to `getMatrixCellSubset`, which prices rows every read, an index pass once, or the cells kept once, and takes the cheapest | read by read, every batch of a 10M run re-decoded the cells' rows (about 64 GB decoded for one run) | a second ask before keeping the cells; charging an unbuilt index's whole pass to every batch |
| 2026-10-30 | a run that reads genes first stops and settles the gene sample's background top-up (`DataState.holdGeneSample`); a person's *Read for all cells* is left alone | the top-up grew the sample under a running DE, which then failed as invalidated | runs racing the top-up |
| 2026-10-30 | a share lowered and raised back comes back in its order with the gene sample it was read for (the gene time's, grown by its top-up, or the count a person read), never fitted again | fitted again at that moment's measured link and decode rate, a 50M store's sample of 4.49M became every cell (markers 15.4 → 62.3 s, a reader worker of 8.7 GB) and a 10M store's markers went 7.8 → 80.7 s | re-pricing the same cells on each Apply |
| 2026-10-30 | categories read the present codes only when the dictionary is larger than the open cells | a 1 % share of a column of 1,029,016 categories sent every label to the page | whole dictionaries at every share |
| 2026-10-31 | a cell subset is read through the first rows of each gene (`_sparseColumnRows`) only where those rows are fewer than the gene's (`_readsFewerColumnRows`: segments ending before the last, or prefix counts below every row); any other read takes the route priced over the run, and a read within a gene sample takes it by the sample's columns (only its rows' pointers are read) | a share raised back to every cell opens in the share's order with its gene index in segments, and two groups of every cell read each batch through the first rows of every gene: whole columns, every gene of the read at once, unpriced (1M: DE 9.0 s against 1.5 s; with no gene index yet, a read of 8,610 genes failed after 3 min 50 s for lack of memory) | the first-rows read for cells that reach the last segment |
| 2026-10-31 | a radix sort's digit has no more values than the keys sorted (`radix-order.js` `digitBitsOf`, 16 bits at most) | a gene's few hundred entries put in open order (`fileOrderedAccess`) walked 65,536 counters a pass: 6 s of a 30,000-gene DE at 1M | a fixed 16-bit digit |
| 2026-10-31 | a run's cell subset is checked once (`_checkedCells`: distinct, in range, its largest) and compared on each later read | a 10M DE reads 1.9M cells about 400 times; checked, bounded by the gene sample, by the rows read and by a prefix each read, they were 8.5 of the reader's 9 s | walking the cells again each read |
| 2026-10-31 | once a run's cells are kept (`_keptCellSubset`), a later read of them is a copy and nothing else (`_heldCellSubset`): no row view, pointer check or named rows again | every cell open with a gene sample (a local 10M CSR file: 8.1M sampled), a DE's groups reach past the sample; each of ~470 reads still made a row view of 2M rows and walked its pointers' coverage three times before copying (0.6 s a read): DE 4 s on the pushed tree, 180–360 s here | per-read work on kept cells |
| 2026-10-31 | a gene stream over some cells starts a read once half the memory gate is free, or every task left fits (`gene-stream-runner.js` `roomToRead`) | a read over some cells costs its cells whatever its genes; started as each task left the gate, a 10M DE of two groups (1.9M cells) made ~15,000 reads: 54 s against the pushed tree's 4 s (it read whole columns, a task a worker a read) | a read a task |
| 2026-10-31 | a subsample view's rows are priced at the stored chunks they lie in (`matrix-access.js` `touchedChunkBytes`: a view chunk touched costs the stored chunk it is part of, decoded whole), | priced at the view's part of each chunk, a 10 % uniform share's groups (lying in every stored chunk) read by their rows cost 5× what was priced (10M CSR: 161 MB priced, 750 MB read and decoded) | pricing a view's rows at its own part of each chunk |
| 2026-10-31 | a run keeping a subsample's cells by their rows reads them in store order (`reader.js` `_subsetInStoreOrder`: the cells' file rows ascending through the open rows' store-order matrix, each stored chunk decoded once, the entries put back at the cells' positions) | read in open order, a spread share's rows made each read of the run's pipeline decode nearly every stored chunk again: a 50 % uniform share's DE groups (944,727 cells) read 7.5 GB of a 2.6 GB 10M CSR file (39–41 s; now 0.18 GB, 5.7 s; every cell 0.22 GB, 5.2 s) | reading a share's rows in open order |
| 2026-10-31 | every stored row open in another order (a share raised back to every cell in its order) takes every stored row's kept index for the rows past the carried ones (`_restrictsWhole`, kept runs past `start` only), its genes placed by the open rows' runs when there are fewer runs than a gene's entries (`restrictedAccess` `entriesByRun`: no lookup scattered through memory, no sort) | the rows past the share's were indexed again by a pass of their own (10M Zarr: 27 MB read and decoded, marker genes 2.7 → 5.3 s when no gene was shown before) while every cell's index was kept; placed by a lookup of each row and a sort, a markers run's 30,000 genes cost 1.0 s (Node), by runs 0.6 s | a second index of every cell after the round trip |
| 2026-10-31 | a file on this computer, or a host on it (a loopback origin: `cellucid serve`, Python's and R's served stores), prices the first gene as this machine at rest whatever its link measured: at no less than a disk's rate (`LOCAL_POLICY`, 1 GB/s), at its shortest round trip, no time a request beyond it (`host-link-memory.js` `pricingLink` `thisMachinePricing`; a local file's link snapshot carries its priced rate and `onThisMachine`) | an opening's metadata reads measure a fraction of such a link (a 10M-cell store on a loopback host in Chromium: 45 MB/s from 0.9 MB, Zarr v3 9 MB/s from 51 KB); priced at that, every cell's first gene (read in about 6 s by the last pushed commit) was priced at 15-37 s, and genes were read for 32-68 % of the cells, the first differential expression on them | the gene sample on this computer from the opening's lower bound; pricing a host here as measured once a read slowed by the machine's own load counts as filling the link, or at its round trip under load (5 of 6 openings of 10M-cell stores with the machine busy sampled 19,779-1,947,096 cells, the first DE following the sample); a notebook kernel priced as this computer (it may run across the internet) |
| 2026-10-31 | a share's choice between every stored row's index restricted to it and its own pass in store order weighs, besides the index's and the view's work, each path's placement of a gene in open order (`wholeIsCheaper`, `PassRates.placing`): per entry over the index built and each gene answered once; unmeasured, the whole for two thirds of the entries or more | every stored row's index places every stored row's entries of each gene, twice a half share's: a 50 % uniform share of a 10M-cell store took 1.3-1.7x on it what its own pass took (Chromium; first gene CSR 7.7 / Zarr 6.4 s against 6.0 / 4.5 s, marker genes 11.6 / 12.4 s against 8.0 / 7.4 s; every cell's first gene 6.0 / 4.9 s), the same values | weighing the index pass alone (a gene was taken to cost the same after either) |
| 2026-10-31 | a gene placed from every stored row's entries tests each against a bit of the segment's open rows before looking its open row up (`openRowBits`) | the open rows' table is 40 MB at 10M stored rows and each lookup a cache miss; half of a half share's entries are of rows not open (Node: 18.8 → 13.8 ns an entry) | a lookup of every entry |
| 2026-10-31 | a run's genes per task follow the genes' bytes as read (`run-tuner.js`: the densest form until a read is measured, sized again at each read and each shrink of the budget) | held at the densest form (8 B a cell of 1.9M) against what the run's context leaves of the share, a 10M-cell differential expression went one gene a task: 30,000 tasks, deliveries and page messages (traced: 30,040 against the pushed tree's 168), every-cell DE 1.2x the pushed tree's; marker genes likewise | a task size fixed at the run's start |

## Measurements that justify it

Real atlases were read from CELLxGENE Census over S3 (HTTP/1.1, several
mirror hosts). The machine was an Apple M1 Pro laptop; browser runs used real
GPUs (Chromium on Metal). "Node" means the app's own reader in Node over the
same real network, with six connections per address as a browser has.

| Dataset (cells) | Engine, date | Step | Before → after |
| --- | --- | --- | --- |
| HLCA core (584,944) | Chromium, 2026-10-28 | whole map; gene sample; first gene; later genes | 17.7 s / 7.7 MB (first picture 10.2 s); 2,143 cells; 11.9 s / 14.5 MB; 0.3 s / 0 MB |
| full HLCA (2,282,447) | Chromium, 2026-10-28 → 10-29 | first picture; whole map; first gene | 17.5 → 8.6 s; 33.4 → 13.7 s; 25.7 s / 1,384 req → 11.6 s / 17.7 MB / 950 req |
| mouse embryo (11,441,407) | Chromium, 2026-10-28 → 10-29 | first picture; whole map; first gene | 27.8 → 10.7 s; 47.6 → 22.7 s; 81 s / 5,379 req → 12.1 s / 13.9 MB / 913 req |
| mouse embryo | Node, 2026-10-29 | first gene once only the read rows' pointers are read | 27.2 s / 39.5 MB / 4,872 req / 1.34 GB → 7.3 s / 13 MB / 863 req / 0.72 GB |
| core / full / embryo | Chromium, 2026-10-30 | revisit: first picture; requests to open | 8.4 → 1.6 s, 9.0 → 2.0 s, 9.4 → 2.8 s; 189 / 455 / 921 → 0 (version checks only); later genes 0 |
| full HLCA / TS (1,136,218) / embryo | Node, 2026-10-30 | rung; cells with the gene at ~5 s and ~20 s | 19 → 152, 15 → 30, 20 → 40; embryo 1,465 @ 5.6 s → 1,400 @ 3.4 s and 5,232 → 8,235 @ 20 s; TS 3,125 → 6,792 @ 15 s |
| full HLCA | Node, 2026-10-30 | before the top-up → with it | 2,729-3,826 cells at 6-8 s and no more → 1,767 @ 4.6 s, 11,055 @ 17.5 s |
| synthetic 10M x 30,000, `.h5ad` CSR / Zarr v3 sharded, loopback web host | Chromium real GPU, 2026-10-31, same hour against the last pushed commit | every cell: first gene, cells without a value, first DE top genes | pushed 6.6-6.9 s / 6.1-6.2 s, 0 cells; before this fix (the independent check's runs) 4.8-5.0 s / 0.8 s for a sample, 2.4-6.7M cells without a value, Zarr DE rank 3 G02402 (pushed G24308); after 5.1-5.2 s / 4.4 s, 0 cells, sums exact, DE top 3 G00000, G00001, G24308 as pushed |

Share against every cell, first gene, measured in the same hour:

| Case | Every cell | Share |
| --- | --- | --- |
| 10M CSR `.h5ad`, emulated S3 (230 ms, 20 MB/s, HTTP/1.1), Node, 1 % uniform | 347,913 cells, 8.04 s, 18.7 MB, 193 req | before 53 cells, 11.13 s, 4.4 MB, 280 req; after 3,488 cells (= s/n), 8.16 s, 18.4 MB, 193 req |
| TS, real S3, Node, 1 % uniform | 1.0 | time 0.40-0.56, bytes 0.16-0.24, requests 0.50-0.65 |
| full HLCA, real S3, Node, 1 % uniform | 1.0 | time 0.78, bytes 0.19, requests 0.82 |
| TS / full HLCA, Chromium, before spread shares used the block order | 8.1 s / 6.1 s | 12.2 s / 17.5 s (the break that design fixes) |
| full HLCA / TS, Firefox, same run | 7.1 s / 10.3 s | 13.5 s / 6.1 s |

Synthetic stores (Node unless named; 20 ms / 50 MB/s HTTP/2 emulated):

| Case | Result |
| --- | --- |
| 10M CSR `.h5ad`, 10 % uniform first gene | 517 MB → 266 MB (= every cell's) once spread rows became one store-order segment |
| 10M, 50 % uniform first gene, `genes: 'all'` | CSR 16.8 → 7.50 s (every cell 6.76 s); sharded Zarr v3 10.16 s (7.75); dense rows 4.89 s (3.77); bytes ≤ every cell's on all three |
| 4M CSR `.h5ad`, 50 % uniform, after the pass's CPU cut | 5.06 → 4.69 s, every cell 4.83 s |
| 10M opening, every cell in a cell order (150 ms / 5 MB/s HTTP/1.1) | first view 344 → 48 requests and open 1,512 → 659 ms, once rows were read in store order |
| 10M growth (`createSubsampleCache`) | priority file growth steps 566 → 46 ms; uniform 50 % growth 1,310 → 186 ms; blocks 50 % 827 → 125 ms |
| 50M sharded Zarr v3, Chromium (Metal, loopback) | markers 15.4 s (a 1,800 s timeout before the row-pointer coverage fix), DE 10.2 s, signature 2.2 s, dot plot 1.2 s |
| 50M and 10M Zarr after a share round trip (10 % then 100 % in place), Chromium | 50M DE 70.3 s (reader 1.8 GB, no allocation failure), markers 62.3 s, signature 4.5 s; 10M DE 83.6 s, markers 80.7 s, signature 1.6 s |
| Node, the reader on a 50M gene sample (9.49M cells), whole columns 8 genes a read | 2,000 genes 22.0 s before; 30,000 genes 3.4 s / 83 MB after the coverage fix (10M: 1,600 genes 19.4 → 2.6 s) |
| Node, 10M CSR `.h5ad`, DE on 1.89M cells over a simulated link | about 13-17 ms a gene (400-500 s for 30,000 genes) → 30,000 genes 14.2 s / 192 MB once priced over the run |
| Node, 2M CSR `.h5ad`, emulated link, share lowered and raised back | 40 MB/s: raise 352,680 cells at rung 15,624 (priced again, 53 s) → 1,690,543 cells at rung 7,812 with the first opening's fit; 15 MB/s: 1,695,490 → 1,086,048 (= the first opening); the raise's first gene reads nothing new |
| 10M rich CSR `.h5ad` DE in Chromium (gene sample 8.35M) | 218.6 s (321 s on a sample, 660 s without one, before); markers 4.5 s |
| 3M Zarr, 1 % uniform obs names | heap 931 → 8 MB |
| a column of 1,029,016 categories (1.2M-cell CELLxGENE atlas, real S3), 1 % uniform | labels sent to the page 1,029,016 → 12,050, worker heap 216 → 94 MB, same 38.9 MB downloaded |
| 2,058,652-cell MatrixMarket stream (GEO, Node) | first view 19.8 s / 106 MB, against at least 7 min for the whole 3.95 GB stream |

On the test machine, Playwright's WebKit reports every public host as
offline, so public links were walked in Chromium and Firefox only. WebKit
runs used local files.

## Tried and reverted

| What | When | Why it went | What would change it |
| --- | --- | --- | --- |
| HDF5 string heaps loaded as one union after every piece | 2026-10-28 | 10M at 150 ms: 86.2 → 68.8 s, but 100M at 20 ms HTTP/2 got slower (72.0 → 87.7 s): a barrier on a fast link | — |
| Fixed batches of a round of the link's requests | 2026-10-28 | 100M HTTP/2 87.1 s, against 71.8 s with no batching | replaced 2026-10-29 by the room rule, now in `HeapStream`: load when the link has a free slot, gather while requests queue (10M 150 ms 77.0 s; 20 ms HTTP/2 within 2 % of no batching) |
| Never fetching again any byte inside a gap | 2026-10-30 | no duplicate bytes, but more requests (5 → 8 in the HDF5 gene pricing test, about 1,000 more at 10M over HTTP/1.1) | — |
| A request budget per heap run for obs names | 2026-10-30 | 0.1 % share 4.95 → 7.69 s and 241 → 376 MB | — |
| Asking held heap bytes at once; a window bounded by file span | 2026-10-30 | every cell 436 → 566 MB; no gain | — |
| A whole heap collection asked up front when its tail was expected to hold a kept name | 2026-10-30 | 15-25 % fewer requests, 0.5-16 MB more bytes | knowing a collection's size before its first read |
| A uniform share's gene-index segments merged into one in open order | 2026-10-30 | 1,019 MB against 517 MB | the store-order segment |
| Dropping the decoded-chunk cache during scans | 2026-10-30 | one chunk read again at every slab boundary | — |
| A single versioned file (a Census `.h5ad`) taken over by version at a growth step | 2026-10-30 | the gene list hung in Node; now done for folders only | the hang traced |
| A block size from the merge gap measured on one read | before 2026-10-30 | the same `.h5ad` gave blocks of 28 to 66 cells from one opening to the next | — (`BYTE_READ_BLOCK_BYTES` is a constant) |

## Known limits and open items

- **CPU on fast links:** a 50 % share's first gene reads through every
  stored row's index (item 9): at 4M and 5M cells, on a near-zero-latency
  link and at 20 ms / 50 MB/s, 1.00–1.04× every cell's in a loaded machine's
  run-to-run noise (same bytes, same requests, the same pass), where its own
  pass took 1.4×. A 50 % share's edges still place them on its cells after
  the last line ({doc}`24_decisions_graph_edges_at_scale`): about 1.03×
  every cell's at 5M. The 10M emulated S3 1 % share took 1.015× every
  cell's time on identical reads (its kept rows' pointers are read before
  its pass, a step every cell's pipelined pass does not wait for).
- **Share round trips (fixed 2026-10-30, browser re-measure pending):** a
  raise in place after a lowered share used to fit the gene sample again and
  opened genes for every cell (50M: 4.49M → every cell; markers after a round
  trip 62-81 s against 15 s without one; the 50M DE allocation failure in
  {doc}`27_decisions_memory_model`). The same cells now come back with their
  sample (see Decisions); the 10M and 50M browser runs of the fixed tree are
  not taken yet. Still fitted over the link as measured: a share never opened
  before (the decode rate measured on small reads under-prices it), and the
  earlier opening's in-memory gene index is not carried through the lowered
  share (the kept index serves it).
- Per-read O(cells) work remains in DE at 10M (the distinct-cells check and
  the DE cells' list copied to the reader on each read).
- The top-up prices the requests of scattered rows about 2× too low: a TS 1 %
  share's top-up made 1,487 requests against 847 priced, 1.19× every cell's.
  Absolute request prices for Zarr folders are low as well.
- A smaller share does not reuse a larger share's kept store-order run on a
  later visit; it indexes again. Applying 1 % → 10 % of a share in a cell
  order is a new reading, not a growth in place.
- Some limits come from the file's layout. A gene-major (CSC) file reads a
  gene's whole column at any share. String heaps and uniform shares of
  row-chunked arrays reach every chunk, so their bytes stay near the whole.
- Memory at a 50 % uniform share stays above every cell's (5M: 711 against
  599 MB peak buffers) because row pointers are held several times. It is
  counted and released.
- Not built yet: 16-bit gene-index positions, an early partial picture within
  the first gene's read, keeping a cut-short conversion or text-matrix pass,
  and genes on the cells read so far during a pass.
- Not yet read or traced: browser proofs of spread shares in the block order
  (Chromium and Firefox, TS and HLCA at 1 % and 10 % uniform), and a WebKit
  crash on the second visit of a 1 % uniform share.
- The kept-index code changed after revision 2 (store-order runs, kept values).
  `KEPT_INDEX_REVISION` is raised once when this is committed, which removes
  every user's kept gene indexes once. The revision tests still pin a digest
  of the source rather than of the produced bytes.
- Planned clean break: rename the wire names `{mode: 'kept'}` and `cellOrder`
  and re-pin the fixtures in the web, Python and R repositories.

## Known issues (pre-commit)

Defects, regressions and unverified behaviour known when this version was
committed (2026-10-30), each with where it is recorded above. An item is
removed in the change that fixes or verifies it.

- **A half share's marker genes:** a 50 % uniform share of a 10M-cell store finds marker genes in 7.1-8.4 s (`.h5ad` CSR) and 7.2 s (Zarr v3 sharded) against 3.8-4.1 s and 3.0 s for every cell after its first gene (Chromium, same hour, 2026-10-31); its first gene reads in 1.08-1.10x (CSR) and 0.84-0.88x (Zarr) every cell's, values exact. The rest is each gene's entries put in the share's open order, on the reader worker alone: traced, 4.5 s of its 7.5 s (a radix sort of each gene's entries by open row 2.2 s, the lookup of each entry's open row 1.3 s, the segments' copy 0.6 s), work every cell's genes, already in open order, never do. Placement over every stored row's index cost more (1.3-1.7x, see Decisions), so a half share now takes its own pass. Removing it needs analyses that take a gene's entries in stored order with each worker holding the groups by stored row, or spread shares opened in store order: not built.
- **A gene sample on this computer with every core busy:** a source on this computer is priced as the computer at rest (see Decisions: a disk's rate, its shortest round trip), and openings right after other runs, or with four CPU-bound processes beside the browser, read the first gene for every cell (10 of 10, Chromium, 2026-10-31). With ten CPU-bound processes on the ten cores, 1 of 10 openings (a 10M-cell `.h5ad` on a loopback host) measured even its shortest round trip at 900 ms (1.5-12 ms otherwise), priced every cell's first gene at 47.7 s and read it for 1,610,169 of 10M cells. Pricing this computer's round trip from a prior (as its rate is from a disk's) would remove it; it changes the progressive opening's and the block rung's prices too, so it is left for a decision.
- **Share round trips:** the fix (the same cells come back with their order and gene sample) is verified in Node only; the 10M and 50M browser runs of the fixed tree are not taken.
- **Spread shares in every cell's block order:** in the browser proofs (Chromium and Firefox; Tabula Sapiens and HLCA over S3 at 1 % and 10 % uniform and every cell, two visits each, no page error) a share's first picture costs what every cell's does (the opening's metadata) and once more: Tabula Sapiens 1 % 12.4 s against 10.9 s for every cell (Chromium), HLCA 1 % 17.1 s against 14.7 s (Firefox); single runs over S3, not separated from its noise. Revisits draw every share in 1.2–1.9 s. The harness's DE read at 1 % reported 2–4 cells read (the gene sample's cells in the two groups), against 138–406 for every cell.
- **CSC `.h5ad` over a slow link (40,000 cells, 150 ms, 5 MB/s):** 100 genes read in 3.3 s / 122 requests on the live tree and the checkpoint alike, 2.5 s / 81 requests in one checkpoint run (the plan follows the link measured from the opening's few reads), 2.8 s / 100 requests at the last push; the DE read 1.8 s (2.8 s at the last push).
- **WebKit:** every app run records the page error "Fetch API cannot load http: /127.0.0.1:<port>/ due to access control checks" (the harness's page port); not traced to the app or the harness.
- **A background top-up on a closed reader** (after Apply) logs "this._root is null"; nothing is affected.
- **Top-up pricing:** the top-up prices the requests of scattered rows about 2× too low (Tabula Sapiens 1 %: 1,487 requests against 847 priced, 1.19× every cell's).
- **Kept gene indexes:** `KEPT_INDEX_REVISION` is raised once at this commit, so every user's kept gene indexes are removed once.
- **A later run of a measured analysis kind takes up to about two and a half times as long to its final result:** once an analysis kind's cost is measured (`progressive-cells.js` `AnalysisCostModel`, kept in memory and in this browser's storage), a later run of it (the same page, or a revisit in the same browser) runs a first pass on a sample and then the pass on every cell. The full pass costs what a single run does, so the time to the final result is that plus the sample pass: 10M-cell Zarr store, every cell, Chromium, a second run on the same page: differential expression 3.0 s → 7.3 s, marker genes 3.3 s → 4.2 s. The figures recorded here before (markers 16 s and DE 36 s after a share round trip) were mostly not this: the harness's `--single-pass` never cleared the costs (it looked the module up in the page's resource timing buffer, which a 10M page fills before the module loads), and the reads after a round trip were slow of their own (fixed: see Decisions, 2026-10-31). Fix in phase 4: the full pass reuses the sample pass's work, or the sample pass is skipped when the full pass is priced short.

## Where the code and tests live

| Concern | Module(s) under `cellucid/assets/js/` | Tests in `cellucid/tests/` |
| --- | --- | --- |
| Share-cost ceiling, auto-fit, ladder, gene sample | `data/anndata/subsample-plan.js` | `subsample-share-cost-contract`, `share-store-order-pass-contract`, `subsample-auto-fit`, `gene-sample-contract` |
| Subsampling methods, rank cache | `data/subsample.js` | `subsample-conformance`, `subsample-growth-cache` |
| Requests and modes | `data/cell-subset-policy.js` | `subsample-address-and-control`, `subsample-recorded-set-contract` |
| Subsample view | `data/anndata/views/cell-subset.js` | `subsample-gather-differential`, `open-cells-rows-contract`, `efficiency-floor` |
| Segments, store order, every stored row's index restricted to a share (`wholeIsCheaper`, `PassRates`) | `data/anndata/prefix-column-index.js`, `data/anndata/reader.js` | `subsample-share-cost-contract`, `share-restricted-whole-contract`, `anndata-sparse-columns-contract`, `subsample-gene-read-memory` |
| Kept indexes, revisions | `data/anndata/kept-column-index.js`, `data/anndata/transcode/space.js` | `kept-column-index-contract`, `transcode-cache-contract`, `kept-parts-contract` |
| Page gene sample, top-up | `app/state/managers/field/gene-sample.js` | `gene-sample-contract`, browser `gene-sample.spec` |
| Progressive opening, growth | `data/anndata/progressive-open.js`, `app/open-cells-resize.js`, `data/opening-cost.js` | `progressive-subsample`, browser `progressive-subsample.spec`, `open-cells-folder-reuse` |
| Strings, categories | `formats/hdf5/nodes.js`, `formats/io/byte-reader.js`, `data/anndata/reader.js` | `obs-strings-share-contract`, `present-categories-contract`, `hdf5-reader-structures` |
| Edges at a share | `data/anndata/graph-lines.js` | `edges-at-a-share-contract` |
| Cell-subset reads (run batching, genes per task, checked and kept cells, after a share round trip), sorting keys | `app/analysis/shared/gene-stream-runner.js`, `app/analysis/shared/run-tuner.js`, `data/anndata/reader.js` (`getMatrixCellSubset`), `data/anndata/matrix-access.js` (`touchedChunkBytes`), `data/anndata/reader.js` (`_subsetInStoreOrder`), `data/anndata/prefix-column-index.js`, `utils/radix-order.js` | `round-trip-analyses-contract`, `analysis-gene-reads-bounded-contract`, `analysis-run-tuner-contract` |
| Text-matrix streams | `data/anndata/transcode/matrix-market.js` | `text-matrix-stream-contract`, browser `text-matrix-stream.spec` |

`cellucid/scripts/efficiency-floor/share-grid.mjs` measures every downstream
at 1, 10 and 50 % shares against every cell, in the whole app.

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
