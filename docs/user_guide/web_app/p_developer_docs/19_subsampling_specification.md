# Subsampling specification

This page is the specification of Cellucid's subsampling: which cells a
subsample keeps, given a request. Three implementations follow it and must
agree on every cell:

| Where | Implementation |
|---|---|
| Web app (the reader worker, while a dataset opens) | `cellucid/assets/js/data/subsample.js` |
| Python | {func}`cellucid.subsample_cells` (`cellucid-python/src/cellucid/subsample.py`) |
| R | `cellucid_subsample_cells()` (`cellucid-r/R/subsample.R`, `src/subsample.c`) |

The shared conformance fixture `subsample-conformance.json` (web app
`tests/fixtures/`, Python `tests/fixtures/subsample_conformance.json`, R
`tests/testthat/fixtures/subsample_conformance.json`) is written by
`cellucid/scripts/generate-subsample-conformance.py`, a plain-Python reading of
this page independent of all three, and the three copies are pinned by one
SHA-256 in all three test suites. Changing anything below means regenerating
the fixture once and updating the digest everywhere.

Version: 1.

---

## Terms

- `n` — the cells in the dataset (`n_obs`), `0 ≤ n ≤ 2^32 − 1`.
- A cell's **position** — its row in `obs`, `0 ≤ p < n`, a uint32.
- `k` — the subsample size, an exact count, `1 ≤ k`. `k ≥ n` keeps every cell
  (no subsample). `k = 0` is refused. (The web app *opens* a subsample of at
  least 2 cells, since one cell has no coordinate variation to draw
  ({ref}`fewest-cells-open`); the selection itself, specified here, keeps any
  `k ≥ 1`.)
- `seed` — a uint32. The default is `0`.
- Every result is the **ascending** positions of exactly `k` cells.

All integer arithmetic below is exact. `^` is bitwise XOR, `>>>` an unsigned
shift, `·₃₂` multiplication modulo 2^32.

## The rank of a cell

```
fmix32(h):          # MurmurHash3's finaliser, a bijection on uint32
  h = h ^ (h >>> 16)
  h = h ·₃₂ 0x85ebca6b
  h = h ^ (h >>> 13)
  h = h ·₃₂ 0xc2b2ae35
  h = h ^ (h >>> 16)

seedKey(seed) = fmix32(seed ^ 0x9e3779b9)
rank(seed, p) = fmix32(p ^ seedKey(seed))
```

For one seed `rank` is a bijection of uint32, so **no two positions ever share
a rank**. Wherever cells are ordered by rank, ties (which cannot occur) would
be broken by the smaller position. "The `m` lowest-ranked of a set" is the `m`
members with the smallest ranks.

## Size from a share

An interface offers a share; the size it stands for is computed exactly:

```
share  = the rational a decimal text names: a percentage text over 100
         (`0.35` → 35/10000), or a fraction given as a binary64 number
         taken as its shortest round-trip decimal: the fewest significant
         digits whose text reads back as `x` under correctly rounded
         conversion, and among those the nearest to `x` (JavaScript
         `String(x)`, Python `repr(x)`, R in C with the C library's
         `strtod`)
k      = min(n, max(1, floor((2·n·P + Q) / (2·Q))))   for share = P/Q
```

that is `n × share` rounded half up, in integer arithmetic (no binary64
product). A share must be in (0, 1] (a percentage in (0, 100]). The size is
always stored as `k` with `n`, never as the share. The fixture lists
fractions (`fraction`) and percentages (`percent`) with their counts,
including shares a converter can get wrong at sizes where that changes the
count: 2^-24 of 25,165,824 cells is 2 (its 16-digit text is
`5.960464477539063e-08`; `…062` reads back as another double, since a power
of two's rounding interval is narrower below it), and 16-digit shares that a
string-to-double conversion that is not correctly rounded reads one ulp off
(R's own, on some builds), whose 17-digit text rounds the other way.

## Methods

### `uniform`

The `k` lowest-ranked of all `n` cells.

Nested: with one seed, the cells kept at size `k` are among those kept at any
size `k' > k`. Moving the slider up keeps every cell already shown.

### `blocks` (`blockSize` `B ≥ 1`)

The cells are cut into blocks of `B` neighbours in file order: block `j` is
positions `[j·B, min(n, (j+1)·B))`, `0 ≤ j < ceil(n / B)`. Blocks are taken in
order of `rank(seed, j)` (the block number ranked like a position), whole,
until the next one would exceed `k`; that block contributes its first
`k − taken` cells.

Nested like `uniform`. With `B = 1` a block's number is its cell's position,
so `blocks` keeps exactly `uniform`'s cells.

`B` is part of the specification of a subsample and is recorded with it:
whatever chose it, a session, a link carrying the subsample and Python's or
R's `subsample_cells` reproduce the same cells from the recorded `B`. When a
request leaves `B` out, the web app takes a rung of the file's **ladder**
(`subsample-plan.js`, `blockRung`, `blockSizeLadder`): the storage's block
`B₀` (below) and its doublings `B₀·2^j` up to the cell count, so a larger
block holds whole smaller ones and every rung stays aligned to the storage.
The storage's block comes from the layout alone (`defaultBlockSize`), from
the **dominant per-cell array**. That is the
expression matrix genes are read from by default (`defaultExpressionPath`):
the matrix the file's launch defaults name (`uns['cellucid']['defaults']['expression']`),
else the default expression rule's (`raw.X` when it holds more genes than
`X` or `X` holds scaled values, else `X`, else the layer the rule picks);
written `M` below:

1. a CSR `M` holding values: a chunked index (`M/indices`) gives
   `B = max(1, round(4 × chunk / m))`, `m` the mean stored values per cell
   (about four of its chunks);
2. else a dense `M`: `B = 4 × ` its chunk's rows;
3. else (`M` is CSC, of another encoding, absent, or a CSR matrix with no
   value) the `obsm` array with one row per cell and the most columns:
   `B = 4 × ` its chunk's rows; with no such array, `B = 1`.

The same matrix prices the gene sample (the open cells a gene is read for,
by gene time and the gene index's memory: `geneTimeFit`, `geneIndexFit`,
`geneSampleOf`): a file whose genes come from a cell-major `raw.X` beside a
gene-major `X` reads genes for the sample a gene of `raw.X` fits. Every cell
of a file in its own order, opened with a gene sample, opens in the priority
order of this method's subsample at seed 0 (`cellOrder`), so the gene sample
of `m` cells is the `blocks` subsample of `m` cells. That order, and a
`blocks` subsample requested without `B` (a share typed or chosen in the
Cells to open control, auto-fit), take the same rung (`blockRung`):

- **a revisit keeps its rung**: the rung a gene index kept in this browser
  was built under (`kept-column-index.js`, filed by method, seed and `B`),
  whatever the link measures now; of several kept, the one holding most
  cells. Another rung would open other cells and read them again;
- **a first visit chooses by request economics**: the lowest rung at or past
  `requests × loadedRoundTrip × rate / (connections × stored bytes a cell)`
  (`requests` what a block's first gene reads: the runs of its arrays' chunks
  and the chunk-index nodes that find them). Below it a sample is bound by
  round trips (fewer cells in the gene time); above it by bytes (no more
  cells, in coarser blocks). The opening's own rate is a lower bound (its
  structure reads are a few requests on new connections), so `rate` is the
  higher of it and the seed: the rate the host sustained on earlier visits
  (`host-link-memory.js`) or this browser's downlink measured on the page's
  own assets (`pageDownlinkBytesPerMs`). `B₀` itself without a measured link.

An array read byte by byte instead of by chunks (stored in one piece, or an
uncompressed Zarr array whose chunks are read in part) has no chunk to align
to: `B` is then the cells whose bytes fill 64 KiB of it,
`ceil(65536 / bytes per cell)` (a CSR `M`: `m` × the index's element size;
a dense `M`: genes × element size; an `obsm` array: columns × element size).
`B₀` is never the link: it is a property of the file. The rung depends on
the link and on what this browser keeps, so two browsers may open a request
without `B` in different blocks; the recorded `B` reproduces either, and a
larger share of the same request (same `B`) keeps every cell of a smaller
one.

### Splitting a size over strata

`stratified` and `spatial` split `k` over **strata**, each a set of cells,
numbered `0 … S−1`. Only occupied strata take part, in number order.

1. Each occupied stratum `s` (size `n_s`) first gets `base_s = min(m, n_s)`
   cells, `m` the minimum per stratum. If `Σ base_s > k` the request is
   refused.
2. The rest, `R = k − Σ base_s`, is split over the strata that still have
   room, capacity `c_s = n_s − base_s > 0`, in number order:
   - **proportional** (largest remainder): `q_s = floor(R·c_s / C)`,
     `C = Σ c_s`, with the remainder `(R·c_s) mod C` computed exactly
     (`R·c_s < 2^64`); the `R − Σ q_s` strata with the largest remainders get
     one more, ties to the lower stratum number;
   - **balanced** (water-filling): the largest level `L` with
     `Σ min(c_s, L) ≤ R`; each gets `min(c_s, L)`, and the
     `R − Σ min(c_s, L)` lowest-numbered strata with `c_s > L` one more.
3. Stratum `s` keeps its `base_s + q_s` lowest-ranked cells.

No stratum ever gets more cells than it has.

### `stratified` (a categorical field)

A categorical field's categories are in their stored order; plain values
(not categorical) are made categories in **code-point order** (Python's
pandas `Categorical`, R `sort(method = "radix")`), never a locale's.

Strata are the field's categories in their stored order, codes `0 … C−1`, and
one more stratum, `C`, for every cell without a category (a missing code:
negative, NaN). Categories with no cell take no part. Allocation
`proportional` (default) or `balanced`; a minimum per category `m` (default
0). Recorded: the field's key, the allocation, `m`.

### `spatial` (an embedding)

The cells are placed on a grid of `G` boxes per axis (`gridSize`, default 32)
over one to three columns of an embedding, and `k` is split **balanced** (no
minimum) over the occupied boxes, so sparse regions keep their cells while
dense ones are thinned.

Coordinates are first rounded to **float32** (as the viewer holds them), then
widened exactly to binary64. Over the cells whose every coordinate is finite,
each axis `a` has `low_a` (minimum) and `high_a` (maximum). A cell's box on
axis `a` is

```
w = high_a − low_a
b_a = 0                                   if w = 0 (or no finite cell)
b_a = min(G − 1, floor(((x − low_a) × G) ÷ w))   otherwise
```

where `x − low_a`, `× G` and `÷ w` are each one binary64 operation rounded to
nearest-even, in that order (never fused). The box is
`b_0 + G·b_1 + G²·b_2`. A cell with any non-finite coordinate is in box
`G^d` (`d` axes), after every other. Strata are the occupied boxes in
ascending box number.

Recorded: the embedding's key, its columns, `G`.

## Why `blocks` is sound

`blocks` is cluster sampling: equal-size clusters of neighbouring cells,
chosen uniformly at random (by rank). Every cell has the same chance to be
kept, so proportions and means estimated from the subsample are unbiased.
Only their variance grows with how alike neighbouring cells are: the design
effect is `1 + (m − 1)ρ` for blocks of `m` cells with intra-block
correlation `ρ`. With thousands of small blocks it stays small; a file sorted
by a category (a sample, a cell type) raises `ρ`, and the Session panel shows
the measured largest difference between a category's share of the open cells
and of every cell, so it is visible.

## Priority-ordered files

A file Cellucid writes stores its cells in **priority order**: row `r` is
the cell of the `r`-th lowest `rank(seed, p)` under the file's seed, with
`uns['cellucid']['priority'] = {seed, specVersion: 1, nObs}` (verified
against the file before it is used: the specification's
{ref}`priority order <priority-order>`, "Verified before use") and each row's
original position in `obs['cellucid_original_position']`; a gene-major `X`
may add `varm['cellucid_prefix_counts']` (n_vars × levels, with
`uns['cellucid']['priority']['levels']`): each gene's stored values among the
first `levels[j]` rows. The writers' side of this contract — the level
sequence (each `floor(previous × 3 / 4)` from `n`), both axes of `obsp`, MuData
modalities ordered by their cells' priority — is the Cellucid-ready AnnData
specification's {ref}`priority order <priority-order>`; `specVersion` is this
page's version.

Then the file's seed is the default seed, and a `uniform` subsample of size
`k` under it is file rows `[0, k)`: every per-cell array is read as one
prefix (at most one partly used chunk each); a CSR `X` as entries
`[0, indptr[k])`; a gene of a gene-major `X` as the first
`prefix_counts[gene, j]` entries of its range for the smallest level
`levels[j] ≥ k` (one more entry when the gene has one: every counted row
below `levels[j]` and the next at or past it, else the file is refused);
`obsp` as rows `[0, k)` with columns `< k`. Any other seed or
method is resolved over original positions (per-cell inputs put in original
order through the file's permutation) and read through the general path. The
reader computes the open cells' original positions from the seed (the
positions sorted by rank) instead of reading them, checks them against the
stored column's first chunk of rows (`priority-layout.js`). A file whose
record does not verify (subset, combined or reordered, then saved again with
anndata, which keeps `uns`) is never refused and never misread: it is read
as a file Cellucid did not write — each row its own original position,
`blocks` its default — and the opening says why. Open cells are in priority order (for
that prefix, and for every cell of such a file, file order); every boundary
still names cells by original position (`cell-identity.js`, `openOrder`).

### Where a store in priority order comes from

Not only `prepare()` and `cellucid_prepare()` files: every store a
Cellucid program makes for the browser is in priority order, under seed `0`
(this page's default), so the prefix path above applies to data of any
source layout.

| Data | The store the browser reads |
|---|---|
| A file `prepare()` / `cellucid_prepare()` wrote | the file, as it is |
| An AnnData or MuData object Python serves (`show_anndata`, `cellucid.Widget`, `serve_anndata`, `AnnDataViewer`) | a Zarr v2 store made from the object as the browser reads it (`cellucid._served_store.PriorityStore`; a MuData as a MuData store): `obs`, `uns` and the other small elements written once in memory; the object's arrays snapshotted into memory-mapped files removed with the server: expression matrices (`X`, layers, `raw/X`, sparse or dense) counted as the store is made and laid out gene-major in priority order by the store's background worker, `X` first, by the one layout pass `prepare()` writes with (`prepare_data/_gene_major.py`), stored as `prepare()` stores them (duplicates summed, zeros dropped, float32 unless integer counts need float64, non-finite refused), the arrays read-only while they are read; other dense arrays in priority order; every chunk made from that copy when it is asked for, a chunk of a matrix not laid out yet made from the snapshot at once, exactly as it will be laid out (cellucid-r's served store the same, its layout on a native thread) |
| An AnnData file Python serves that is not in a verified priority order (`.h5ad`, Zarr folder or zip; `cellucid serve`) | the same store, made from the file's own arrays (h5py, zarr): a cell-major `X` laid out gene-major by the worker as above, its values as stored (layers and `raw/X` when first asked for), so a gene is one range; a gene-major one renumbered with `varm['cellucid_prefix_counts']` counted in one pass when first read |
| A MuData file Python serves that is not in a verified priority order (`.h5mu`, MuData Zarr folder or zip) | a MuData Zarr store made the same way: the shared cells in priority order, each modality an AnnData group as above with its rows in the order of its cells (`obsmap` renumbered; a modality's prefix levels are its rows of cells below each global level), so a prefix of the cells is a prefix of every modality |
| A Loom file or a 10x Cell Ranger `outs/` folder Python serves (`cellucid._served_formats`) | the AnnData the browser's layout view presents from the file (`views/loom.js`, `views/tenx-h5.js`, `views/tenx-outs.js`: the same `obs`, `var`, embeddings, graphs and names), as the same store: a Loom matrix read transposed by row spans, a Cell Ranger matrix cell-major with its rows in priority order; several feature types (or genomes) are a MuData of modalities, each a column window of the one matrix. A Cell Ranger `.h5` alone, or a Loom file, with no embedding is refused in the browser's words |
| An R object R serves (`cellucid_view`, `cellucid_serve`, the notebook widget) | the `cellucid_prepare()` plan in memory, in priority order |
| An R object or MatrixMarket folder the browser converts | the converted store (origin-private file system), written in priority order with the same records (`transcode/priority.js`); an uncompressed R file, whose arrays are read in place and never copied, keeps the object's order and is read in blocks |

In each, a cell's original position is the one `prepare()` and
`cellucid_prepare()` give it, the same in every path, by the one rule for a
stored order (`stored_order()` in Python, `.stored_order()` in R,
`readPriorityLayout()` here): `uns['cellucid']['priority']` decides, and it
verifies only when its seed and `specVersion` are valid, its `nObs` is the
cell count, its `levels` are `prefixLevels(nObs)` and the position column,
when there, is the seed's order. An object in a verified stored order (read
back from a file in priority order, its record kept) is put back in its
original order first, and a file in a verified order is served as it is, so
every cell keeps the position the file gave it. Anything else is cells
Cellucid did not write: each row is its own original position (the object's
row order, the file's rows), and a position column with no record that
describes the cells — a column alone, a subset or reordering of such an
object, a column that is not the seed's order — names nothing: it is not
served or written, and a warning says why, in the same words on every path.
So a selection, a highlight or a session lands on the same cells whichever
path wrote or serves the data. A MuData object Python serves is served as a
MuData Zarr store in priority order, from a snapshot, with no file written. The web app's own tests open each
foreign file Python serves both ways — the file through the browser's view,
the served store — and compare them cell for cell
(`tests/served-priority-store.test.mjs`). Files the servers serve as they
are — MatrixMarket folders and R files given by path to Python, which the
browser reads once from the start and converts as above; any file given by
path to R, which has no HDF5 reader and no reader of Zarr chunks, and says
so when it serves one not in priority order — are read by the browser in
storage blocks (`blocks`).

The servers' chunks are uncompressed (read over loopback or a notebook's
connection). The browser reads a run of an uncompressed chunk's elements
— a prefix of rows, a gene's counted entries — as those bytes alone when
the bytes it leaves unread cost more than a request on the link
(`formats/zarr/array.js`), so the bytes a subsample transfers grow with `k`
beyond the chunk layout. Measured on a 1,000,000-cell × 2,000-gene AnnData
object (100 million stored values, CSR) served by Python: see the
performance page of the Python package
({doc}`../../python_package/d_viewing_apis/14_performance_scaling_and_lazy_loading`).

## Priority order: the order a subsample's cells open in

The kept cells (ascending, as every method returns them) are opened in their
**priority order**:

- `blocks`: by `rank(seed, j)` of their block `j = floor(p / B)`, then by
  position within the block — exactly the order `blocks` takes them in;
- every other method: by `rank(seed, p)`.

The first `m` cells of this order are a sample at every `m`: for `uniform`
and `blocks` exactly the subsample of size `m` under the same seed (both are
nested), for `stratified` and `spatial` a uniform sample of the cells they
kept. The fixture lists each case's order (`order`; `orderHead` and
`orderTail` for the large one), which the web app's `priorityOrder` is held
to. It is the web app's own: everything that leaves the app names cells by
original position, and a session saved on a subsample holds its per-cell
state in **ascending original position** (`session-cell-order.js`: cell `s`
of a chunk is `cells.kept[s]`), so Python and R apply it without knowing any
order.

What the web app builds on it (`subsample-plan.js`, `views/cell-subset.js`,
`prefix-column-index.js`):

- **Reads in priority order.** A read of open cells `[0, m)` fetches the
  chunks of the `m` highest-priority cells first; a view's reads are planned
  over the stored rows ascending (two rows share a read under the rule
  below) and each run is copied to its own place, bit for bit.
- **A gene index in segments.** A cell-major matrix (CSR, or dense not
  stored in gene strips) whose open cells are in priority order is indexed
  in segments of rows — the first about a thirty-second of the cells (at
  least one slab of the index's work), each next one as large as all before
  it — built one after another. A gene asked for meanwhile is answered after
  each segment on the cells indexed so far (NaN, shown as missing, on the
  others): every partial answer is a sample of the open cells.
- **Growing and shrinking in place.** Opening the same file again with
  other cells (the subsample control) reads it in the same reader worker
  through the store already open, and the new reading takes over the gene
  index of the open cells both share as a prefix: raising a `uniform` or
  `blocks` subsample under the same seed indexes only the new cells, and
  lowering it reads nothing. Nothing is carried for cells the two openings
  do not share (another method or seed starts the index anew).
- **Kept across visits** (`kept-column-index.js`). Each segment built (and
  the whole index of a matrix opened with every cell) is kept in the
  origin-private file system as one directory of `transcode/space.js`,
  filed under the group: the file's identity, the matrix, its shape and the
  order key (the subsample's method, seed and method settings, never `k`;
  `file` for every cell). It records its first open row, its rows' file rows
  and row pointers, and the store as built: the entry positions of each gene
  (a structure, over a link slower than a disk), each gene's rows and values,
  or a dense column store. A later opening takes every kept run whose file
  rows equal its own at the same open positions — the whole run, or its
  first rows when the opening has fewer (a prefix of every gene's entries) —
  and builds segments only where none is kept. Runs are shared exactly as
  far as two openings' priority orders agree from the start: all of the
  smaller for `uniform` and `blocks` (nested), only the common start for
  `stratified` and the other methods that draw each size apart. Before use
  a run is checked (`findKeptRuns`): manifest revision
  (`KEPT_INDEX_REVISION`), identity digest and sizes; its file rows against
  the opening's; its row pointers against the file's, read now (a mismatch
  removes the run); the CRC-32s of its column starts and per-gene checksum
  tables. The file's version is the one rule of every kept kind: runs are
  filed under the identity this opening established from the file's source
  (a content name where the source has one, else its size with its
  Last-Modified date and strong ETag, `space.js` `writeVersion`; a local
  file's size and modification time), and opening a version first forgets
  what is kept of the file's other versions (`forgetOtherVersions`), so a
  run of another version is never found.
  Each gene's piece of every per-gene array (positions, rows, values)
  carries a CRC-32 checked as it is read; a mismatch removes the run, the
  gene is read from the file and the run is built and kept again (CRC-32,
  not SHA-256: 1.6 µs against 23 µs for a 3.2 KB piece, 45 µs against 93 µs
  for 96 KB, Node 26 on an M1 Pro; it finds accidental damage, not
  tampering). Kept across visits only under an identity that names a
  version (`persistsAcrossVisits`), else for the opening that built it (and
  its readings again with other cells). A kept run makes room
  only from other kept runs (`optional`), and after a `QuotaExceededError`
  removes more and writes again while that frees anything; the index being
  built (`column-store-placement.js` `diskFreeBytes`) and a conversion make
  room from runs before converted stores, least recently used first within
  each. Every run an opening takes is held (a Web Lock) until its reader
  closes, and the writer of a run lets go of it only once the page holds it
  (`{storeHeld}`).
- **Kept work is priced at no download.** The gene sample's fit by the
  first gene's read time (`geneTimeFit`) leaves out of the cost the first
  open cells whose index runs are kept (`kept-column-index.js`
  `keptIndexPrefix`, from the runs' manifests, for the open cells' order:
  a subsample's, or every cell's `cellOrder`, filed as the `blocks`
  subsample whose order it is) and every chunk of the matrix's arrays this
  browser holds already (`heldChunkBands`: chunks of a file's downloaded
  parts, `formats/io/kept-parts.js`, read from this machine's disk). A
  later visit therefore reads genes for more cells within the same target;
  the fit records what was kept (`genes.fit.kept`: `{cells, bytes}`), which
  the details table names (*Kept here*).

## Edge cases

| Case | Result |
|---|---|
| `k = 0` | refused |
| `k ≥ n` | every cell (no subsample) |
| a seed or position outside uint32 | refused |
| an empty category | takes no part |
| a minimum per category that does not fit `k` | refused (fixture `refusals`) |
| `gridSize^d ≥ 2^53` | refused (box numbers would not be exact) |
| a spatial column that is negative or not an integer | refused |
| NaN or infinite coordinates | the last box |
| every coordinate of an axis equal | that axis' box is 0 |
| no cell with finite coordinates | every cell in the last box (`uniform` over them) |

## Cell set digest

The digest names a kept set: 16 lowercase hex digits from two 32-bit lanes
over the ascending positions `p_0 … p_{k−1}` (for "every cell", `0 … n−1`):

```
a = 0x811c9dc5 ^ n
b = 0x9e3779b9 ^ k
for each p:
  a = (a ^ p) ·₃₂ 0x01000193
  x = b ^ p
  b = (rotl32(x, 13) ·₃₂ 5) + 0xe6546b64   (mod 2^32)
digest = hex8(fmix32(a ^ k)) + hex8(fmix32(b ^ n))
```

## Cell set encoding (sessions)

A session stores a kept set as bytes: for each ascending position, the gap
from the previous one (the first from −1) minus one, as an unsigned LEB128
varint. Every cell is no bytes. In the session manifest the bytes are base64,
under `cells.kept`, beside `cells.spec` (the specification: `method`, `nObs`,
`k`, `seed`, and the method's own keys) and `cells.digest`
({doc}`10_sessions_persistence_and_serialization`).

## Where subsampling happens in the web app

The reader worker resolves the request (`assets/js/data/anndata/subsample-plan.js`)
before the dataset's structure is read, reading only what the method needs (a
stratified subsample the field's codes; a spatial one the embedding's
columns). It then reads the dataset through a view of the kept cells
(`assets/js/data/anndata/views/cell-subset.js`): everything above it — the
viewer, filters, analyses, gene reads, exports — sees a dataset of `k` cells.
Cell `i` of the view is original position `openOrder[i]`, the kept cells in
priority order (`assets/js/data/cell-identity.js`); everything that leaves
the app names cells by original position.

| Read path | How it reads the kept cells |
|---|---|
| obs columns, categorical codes, nullable values and masks, obs index | the kept rows' runs, chunk by chunk, in the stored order: a read of many open cells goes over the stored chunks in ascending order, a slab after another (`readInSourceOrder`, `values.js` `readRowsInStoreOrder`), each slab's rows put where they belong (the whole view of a column row by row in the open rows' file order, `cell-subset.js` `rowWaves`: one copy a row, no piece listed), so the requests are the stored array's whatever the open order (every cell of a 10M-cell `.h5ad` in its `blocks` order: the embedding in 24 requests; read a slab of open cells at a time, 318). A column of variable-length strings (the obs index, a text column; HDF5 or Zarr) reads its kept rows alone (`formats/store.js` `readRows`): an HDF5 column the references of the chunks they lie in, then of the global heap only the collections those rows point into, each read up to the furthest object asked for (its other references name the objects before it), asked for as each chunk's references arrive while the link has room for a request, gathered into one plan while requests queue (`hdf5/nodes.js` `HeapStream`, `byte-reader.js` `whenRoom`); a Zarr column the chunks they lie in, each kept as stored (its bytes and an offset an item) and only the rows asked for decoded. A share holds and decodes its own strings, never every string of the chunks it touches, and a heap collection only while a chunk being decoded points into it. The byte reader fetches each byte in flight once: a read of bytes another request carries waits for it, and no plan fetches a gap across them (`byte-reader.js`, "Each byte in flight once") |
| dense `obsm`, dense `X`, layers, `raw/X`, stored outlier quantiles | the same |
| CSR `X` / layers / `raw/X` | the kept rows' `indptr` only, when the matrix is first read (never at opening; another opening of the same store takes over what was read), and of a gene sample only its own rows' (`readRowPointers` given `rows`, `rowsMatrix`); the kept rows' entry ranges only; the gene index is built over them alone, in segments of priority order; over a link, values stored in chunks read whole are read with the first gene beside the indices and kept (`buysValues`) |
| CSC `X` / layers / `raw/X` | a gene's own entries, the kept cells' kept |
| `obsp` graphs | the kept rows' `indptr` and entry ranges; columns not kept dropped; renumbered. The values of the edges kept are read after the indices, and only theirs, when the values of the edges dropped would cost more than one round trip of the link (always for a local file); else every kept row's values are read with the indices. The Dataset box says how much of the open cells' graph that keeps (the other edges lead to cells not open) |
| modalities (MuData, 10x feature types, Seurat assays) | each modality's rows the kept cells measured, its cell map renumbered |

Two kept rows share one read of the stored array only when no chunk lies
wholly between them (a chunk is decoded whole), or, for an array read byte by
byte (stored in one piece, or an unfiltered HDF5 chunk), when the gap costs
less than a request on the link as measured. A read never fetches a chunk
that holds no kept row. `tests/subsample-gather-differential.test.mjs` holds
every path to "read every cell, then take the kept rows", bit for bit, and to
that chunk budget, over randomly drawn stores; `tests/efficiency-floor.test.mjs`
holds every element, at every share, to its floor: each stored chunk that holds
a kept cell fetched once, none other (`scripts/efficiency-floor.mjs` measures
the same floor over real files, links and the HDF5 reader).
