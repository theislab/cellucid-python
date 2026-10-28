# Subsampled loading: open a share of the cells

A dataset of millions of cells can be opened with a share of them: a seeded,
reproducible subsample, read **while the dataset opens**. Every part of
Cellucid then works on those cells alone — the view, colours, filters,
highlights, analyses, gene expression, exports — and only the parts of the
file that hold them are read. Over a network that is less to download; on any
device it is less memory.

## Cells to open

**Session › Open data › Cells to open** holds one slider and the exact count it
stands for ("228,245 of 2,282,447 cells").

1. Move the slider (or type the percentage) below 100 %. The settings of the
   subsample appear: how the cells are chosen, its seed, and what the method
   needs (a field, an embedding).
2. Press **Apply**. A larger or smaller share of the same method and seed
   changes the open cells in place ({ref}`raising-and-lowering-the-share`).
   Other cells are read again, in the same reader and through the parts of
   the file already open. Your
   highlights, label fields, columns a notebook pushed and cells you hid stay
   on the cells that are still open, matched by each cell's original
   position: a highlight none of whose cells is open any more is named and
   dropped, a label field you set by hand leaves the cells it never labelled
   empty, and a pushed column has no value (NaN) on cells it never had one
   for. What you made from the file's columns takes the newly open cells as
   it does when the share is raised in place: a highlight made from a
   category or a range keeps that rule and takes the new cells in it that
   are visible, and a label field made from a column labels them from it.
   Where you were stays
   too: every view's camera and layout, the field each view is coloured by and
   its filters, the sidebar settings, and the open analysis windows. A part
   the new cells cannot take (a filter on a category none of them has) is
   named in a notice; the rest is kept. What the view shows then is exactly
   what opening the dataset at that share directly shows.
3. Datasets you open next in this tab open with the same share and method.
   A method a dataset cannot use (a stratifying field or an embedding it does
   not have) opens the same share by the file's default method instead
   ([How the cells are chosen](#how-the-cells-are-chosen)), and says so.

Back to 100 % and Apply opens every cell: in place for a file Cellucid wrote
(only the cells not open yet are read), read anew for any other file.

The **Dataset** box shows what is open: "228K of 2.3M cells", and under
Details the method and seed ("uniform, seed 0").

## How the cells are chosen

| Method | Keeps | Use it when |
|---|---|---|
| **The file's default** | for a file Cellucid wrote (stored in priority order): uniform under the file's own seed; for any other file: blocks aligned to how it is stored | always the first choice: either way only about the chosen share of the file is read |
| **Blocks** | whole runs of neighbouring cells | a file Cellucid did not write: parts of it that hold only left-out cells are never read. On a file Cellucid wrote, blocks are scattered through it and read more (the panel says so) |
| **Uniform** | every cell equally likely | on a file Cellucid wrote, with its own seed (the default), it reads only its share; with another seed, or on any other file, its cells lie in most of the file's chunks, so it reads more (the panel says how much: "fetches 45 % of the matrix for these 10 % of the cells") |
| **Stratified** | every category of a field in proportion, or the same count of each (*Same count per category*) | rare cell types must stay visible, or each sample should keep its share |
| **Spatial** | cells on a grid over an embedding, dense regions thinned, sparse ones kept | the map's outskirts matter more than its crowded centre |

With **Blocks**, the panel also shows how far its categories are from every
cell's ("largest difference from every cell in cell_type: 0.3 percentage
points").

Every method is **seeded**: the same dataset, method, size and seed always
open the same cells, in this browser, in another one, and in Python and R
(`cellucid.subsample_cells()`, `cellucid_subsample_cells()` give the same
cells). With **Uniform** and **Blocks**, a smaller share keeps a subset of the
cells a larger one keeps, so raising the slider adds cells and keeps the ones
you saw.

(raising-and-lowering-the-share)=
### Raising and lowering the share

With **Uniform** and **Blocks**, a smaller share is the beginning of a
larger one under the same seed: the cells open in that order (their
*priority order*), so the first cells of any share are themselves a sample.
So is every cell of a file Cellucid wrote, opened whole: its first cells are
its uniform share under its own seed. Changing the share then changes the
cells **in place**: nothing shown is read again.

- **Raising** the share reads only the new cells: their positions in every
  embedding loaded, their values of every field and gene loaded, and the gene
  index of a matrix stored cell by cell (most atlases, every CELLxGENE file)
  for them alone. The new cells join the view where they belong, coloured by
  the field shown, under its filters. What you made from the file's columns
  takes them as opening the larger share would: a highlight made by
  clicking a category or brushing a range takes the new cells in it that
  are visible, and a label field made from a column (a copy, merged or
  deleted categories, a consensus) labels them from that column. A highlight
  or label you made by hand keeps its cells, and a selection under way
  stays.
- **Lowering** it reads nothing: the cells that stay keep everything; a
  highlight keeps those of its cells that are still open (one left with none
  is named and removed); a selection under way is cancelled.

On the graphics card the same holds: the cells kept stay where they are
and only the cells added are sent (a lower sends nothing, or at most the
last 2,097,152 cells when the view had already put them in its spatial
order), and every level of detail draws at once, since a level is the first
cells of that same order. The spatial order that culling uses is rebuilt in
the background and taken a few million cells a frame, so no frame stops for
it.

Either way the camera does not move (the cells shown keep their
coordinates: an embedding stays in the frame its first cells were shown
in), each view keeps its field, filters and settings, and analysis windows
stay open, their results marked **Stale · Run** since their cells changed.
An analysis still running, or a field still loading, finishes on the cells
it started on first (its result is then marked **Stale · Run**); a field or
gene picked while the cells change loads on the new cells. An embedding
whose new cells hold a value that cannot be drawn can no longer be shown,
as when opening those cells (said); when it is the one shown, Apply opens
the cells anew.
Statistics follow the cells: a filter or colour range at a field's bound
moves with the bound, one inside it stays (one no open cell lies in is reset,
and said). Kept views are made again from their settings. Outlier filters
in use get their outlier measures again, over the cells now open.

Another method or seed opens a different set of cells, which is read anew
(through the parts of the file already open), as described under Apply above.

### Large datasets open progressively

When reading the cells asked for is measured to take more than about a
second — over a slow link, or for a dataset this device takes long to draw —
Cellucid shows a first view of fewer of them at once, then grows it in place
to the cells asked for, doubling the cells open at each step. Every step is
a sample of those cells (the first cells of their priority order), so the
picture, the colours, filters and analyses are right from the first view;
**Cells to open** says how many are open so far while it grows. The first
view is sized from what reading is measured to cost here: the link's
measured speed and latency, the stored bytes of the parts of the embedding
its cells lie in (a file's parts are the least it can deliver), and how long
this device took to show the cells of the datasets it opened before, from
their being read to their first picture, apart from what showing any
dataset takes (kept in this browser; growing or
cutting the cells shown in place is not counted: it keeps the cells shown on
the graphics card and takes in only those added, so it says nothing of what
opening them costs). With nothing
measured, when every cell is predicted within the second, or when a first
view would not arrive at least twice as soon (its cells lie in most of the
parts every cell does: a blocks share of a file stored in parts larger than
its blocks), a dataset opens at once. A file Cellucid wrote opens
progressively whole or by a **Uniform** or **Blocks** share; any other file
by a **Uniform** or **Blocks** share (a file too large for the device opens
so, fitted). Stratified and spatial shares, every cell of a file Cellucid
did not write, and a session's recorded cells always open at once.

Some files store the embedding in parts much larger than a share's blocks
(the lung atlas: parts of 8,916 cells), so the first cells of a share lie in
every part and even a small first view reads the whole embedding. When
reading the share part by part brings the first view at least twice as
soon, Cellucid opens it so instead: the share stays exactly the one asked
for, its parts are taken in a random order fixed by its seed (so even the
first view is spread over the file), and each step shows the share's cells
in the parts read so far; each part is read once. These steps are not yet a
sample of the share (**Cells to open** says so); once every part is read,
the cells are put in their order in place, and what is open is exactly what
opening the share at once opens.

### Genes fill in as cells are read

The first gene of a matrix stored cell by cell needs every open cell's
entries read once (its notice says how much, before it starts). With a
subsample, the gene's colours appear before that read ends: first on about a
thirty-second of the open cells, then on twice as many, and so on, the cells
not read yet shown in the missing-value colour. Each step is a sample of the
open cells, so the picture is right from the first colours; the notice says
how many cells it is shown on so far. Later genes read only themselves.

### The gene index is kept for later visits

The index a matrix stored cell by cell needs before its first gene is kept
in this browser (its private storage, like the files read whole: the
Session panel's *Stored* line says how much, and **Clear** removes it), under
the file, the matrix, and how the cells were chosen (method, seed and its
settings; never the share). Opening the same file again later — after the
tab closed, in another tab, another day — with the same method and seed reads
its genes at once: the index of the cells it shares with an earlier opening
is read from this browser, and only cells no earlier opening indexed are
indexed now. With **Uniform** and **Blocks** a larger share under the same
seed holds the smaller one's cells first, so it indexes only its new cells,
and a smaller one none. **Stratified** draws each size apart, so two sizes
share only the cells their orders agree on from the start (in a test file of
100,000 cells in five groups, stratified shares of 10,000 and 20,000 cells
share their first 9,678): the rest is indexed again.

It is kept across visits only for a file whose version is known (a local
file's modification time, a CELLxGENE Census release, a link's `ETag` or its
provider's revision or checksum, the `ETag` a viewer host — `cellucid
serve`, a notebook, R — lists a file with), never by its size alone: a file
known only by its size (a link that names no version) keeps it for its tab
only.

Before a kept part is used it is checked:

- against its own record: its rows are this opening's, its sizes and
  checksums are whole;
- against the file's row layout (a matrix stored sparse): the number of
  stored values in each of its cells, read again (a part of a file whose
  rows changed is removed);
- against the file's content, unless the file's version names its content (a
  checksum or content hash: Google Drive, Figshare, Zenodo, Hugging Face,
  Git LFS; a fixed CELLxGENE Census release; a OneDrive content tag). A
  local file's name, size and modification time, a `Last-Modified` date,
  most hosts' `ETag`s and the `ETag` of `cellucid serve` or
  `cellucid_serve()` (the file's size and modification and change times)
  say only when the file was last written: a rewrite of the same size
  within the time they are kept to (a second for `Last-Modified`, two on a
  FAT disk, the file system's own for a served file), or a copy that keeps
  the modification time (`cp -p`, `rsync -t`), leaves them as they were. So
  the part's fingerprint is read again: the values of one cell
  drawn at random in each doubling stretch of its cells (about 17 cells for
  100,000), a few chunks of the matrix, read at once. A file reordered or
  normalised again is found and indexed anew; a change to a few cells only
  can go unnoticed.

Each gene read from a kept part is checked against a checksum written with
it; a part found damaged then is removed, the gene read from the file, and
the index built and kept again. When the browser's storage is full, kept
parts no open tab uses make room first (for the index being built, a file
being converted), least recently used first; a kept part never removes a
converted R or MatrixMarket file to make room for itself.

With **Blocks**, cells next to each other in the file come together. Files
written sample by sample (most atlases) keep whole stretches of a sample: the
subsample is still seeded and reproducible, but less even than Uniform over
samples. Use Stratified on the sample field when that matters.

(in-a-link)=
## In a link

A link can carry the subsample, so everyone who follows it opens the same
cells:

```text
https://www.cellucid.com/?data=cellxgene:9f222629-…&cells=10%&sample=blocks&seed=0
https://www.cellucid.com/?data=gh:me/atlas&cells=250000&sample=stratified&by=cell_type
```

| Parameter | Value |
|---|---|
| `cells` | a share (`10%`, `0.5%`) or an exact count of at least 2 (`250000`; {ref}`one cell cannot be drawn <fewest-cells-open>`); `100%` opens every cell |
| `sample` | `uniform`, `blocks`, `stratified` or `spatial` (default: the file's, as above) |
| `seed` | a whole number from 0 to 4294967295 (default: the file's own seed for a file Cellucid wrote, else 0) |
| `block` | `blocks`: cells per block (default: from how the file stores its cells, {doc}`../p_developer_docs/19_subsampling_specification`) |
| `by`, `alloc`, `min` | `stratified`: the field, categorical or text (a text field's values are its categories, in code-point order); `balanced` for the same count per category; a first minimum per category |
| `embedding`, `grid` | `spatial`: the embedding; boxes per axis |

A link's subsample is used as given, even where the device could hold more
or fewer cells. A malformed value is named in a message, and the dataset opens
as if the link named no subsample. After you Apply a subsample, the address
shows it (as its exact cell count), and **Copy link** under *Cells to open*
copies a link that opens exactly these cells. The link names how the cells
are chosen, never the cells themselves; the Session panel shows the cells'
digest ("Cell set 3f2a…"), so two people can confirm they see the same cells.

## Datasets too large for this device

When opening every cell would not fit this device's memory — the per-cell
data the view holds (77 bytes a cell, as measured with levels of detail on,
against half the device's memory), or, for a cell-major (CSR) matrix, the index that makes
its genes readable — Cellucid opens a subsample sized to fit, by the file's
default method (uniform under the file's own seed for a file Cellucid wrote,
**Blocks** for any other), and says so above **Open all cells**, with how far
its categories are from every cell's. The budget is what the browser reports of the device
(`navigator.deviceMemory` and the JavaScript heap limit, see
{doc}`../n_benchmarking_performance/index`); a dataset that fits always opens
whole. **Open all cells** reads every cell of this dataset anyway; datasets
opened next still open fitted.

## Neighbour graphs

A subsample keeps the edges of the neighbour graph between two open cells;
an edge from an open cell to one left out cannot be drawn. Once the graph is
read, the Dataset box says how much of the open cells' graph that is, beside
its edge count ("10 % of the open cells' edges; the others lead to cells not
open"): with Uniform at 10 %, each neighbour of an open cell is open with a
chance of about a tenth, so about a tenth of the edges remain.

(analyses-on-a-subsample)=
## Analyses, filters and overlays on a subsample

Everything works on the cells open, as if the file held only them:

- **Analyses** (Differential expression, Marker genes, Gene signature,
  Correlation, Composition, Dot plot, Detailed, Quick Insights) compute on
  the open cells: a page, category or group is its open cells, and the
  counts they show are of open cells. A group with fewer than two open
  cells is not tested; a category with no open cell has nothing to test or
  count (Differential expression and Marker genes do not offer it); a gene
  no open cell expresses shows no difference (p = 1) and no correlation
  (r = 0, p = 1). A first answer on a sample
  ({doc}`../h_analysis/02_analysis_ui_overview`) is a sample of the open
  cells, and the final answer is exactly a run on every open cell. Per-cell
  exports name each cell by its row in the whole file.
- **Every count of the open cells says how many of every cell they are**:
  *Showing all 240 of 1,200 cells* (or *Showing 120 of 240 of 1,200 cells*
  with a filter) under Coloring & Filtering, *Counted on the 240 of 1,200
  cells open.* above a categorical legend's counts, *Total: 240 of 1,200
  cells* in the category builder, and the dataset picker's *Atlas (240 of
  1.2K cells)*, besides the Dataset box and **Cells to open**. A count of
  some of the open cells (a category, a highlight, a page) is of the open
  ones.
- **Filters** (categories, value ranges, genes, any of them at once) act on
  the open cells; a range's bounds are the open cells' values. Outlier
  filters use the quantiles a file stores, which describe the whole
  dataset, or compute them on the open cells and say so
  ({ref}`outliers-on-a-subsample`).
- **Highlights, lasso and labels** name open cells; in a notebook and in
  exports they are named by their row in the whole file, and cells not open
  are unknown (NA).
- **Overlays**: the neighbour graph keeps the edges between open cells
  (below), a vector field shows the open cells' vectors, and Data mapping
  offers each neighbour graph on the open cells.
- **When the open cells change**, an analysis result is marked
  **Stale · Run**, and is not saved, exported or sent to a notebook until it
  is run again.

(fewest-cells-open)=
### At least two cells

A subsample opens at least **2 cells**. The limit is not chosen: the viewer
draws an embedding only when its coordinates vary, and one cell has no
coordinate variation (the rule a one-cell file is refused by:
{ref}`embedding-refused`, *“no coordinate variation”*). So *Cells to open* does
not offer a share of one cell (its count line says why and **Apply** stays
off), a link with `cells=1` is refused with that reason, and auto-fit never
fits a dataset to fewer than two cells. A share is never turned into another
count: one that is a single cell of the dataset it opens (a share chosen on
a larger dataset, `cells=0.01%` in a link) opens that one cell, and its
embedding is refused by the same rule, as are two or more open cells that
all share one point.

`cellucid.subsample_cells()` (Python) and `cellucid_subsample_cells()` (R)
select cells, and keep any count from 1: one cell is a selection, just not
one the viewer can draw.

## Sessions

A session saved on a subsample records exactly which cells were open (the
cells themselves, compactly, and how they were chosen). Loading it opens the
dataset with those same cells — reading it again first if it is open with
others — so every highlight and label lands where it was. A session saved on
every cell opens every cell. Applied in Python or R
(`cellucid.apply_cellucid_session_to_anndata()`, `cellucid_apply_session()`),
each cell's state goes to its own row of the whole dataset. The page address
shows a subsample you chose as its specification ({ref}`in-a-link`), never
the cells themselves.

## From a notebook

```python
import cellucid
kept = cellucid.subsample_cells(adata, fraction=0.1, method="uniform", seed=0)
small = adata[kept].copy()          # the cells the viewer opens at 10 %, seed 0
```

```r
kept <- cellucid_subsample_cells(seurat_object, fraction = 0.1, seed = 0)
small <- seurat_object[, kept]
```

A viewer opened from a notebook names every cell by its row in the whole
dataset: selections and events report original positions whatever share is
open.

Python serves what it shows with the cells in priority order, so a share the
viewer opens is the first rows of every array and only its cells are sent,
whatever the source's layout: an in-memory AnnData or MuData, an AnnData file
(`.h5ad`, Zarr), and by path a MuData file (`.h5mu`, MuData Zarr), a Loom
file and a 10x Cell Ranger `outs/` folder, each presented as the browser
presents the file itself. R serves an R object in priority order; a file
given to R by path is served as it is (R has no HDF5 reader), so a share
reads the file's storage blocks, as from a link, and R says so when it
serves it — read it into R and view the object, write it once with
`cellucid_prepare()`, or serve it from Python, to send only the share's
cells.

## What is read, exactly

The specification — the rank of a cell, each method, the rounding of a share
to a count, the float rules of the spatial grid, and how each part of a file
is read for a subsample — is
{doc}`../p_developer_docs/19_subsampling_specification`.
