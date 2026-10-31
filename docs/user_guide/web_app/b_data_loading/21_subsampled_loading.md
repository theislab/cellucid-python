# Subsampled loading: open a share of the cells

A dataset of millions of cells can be opened with a share of them: a seeded,
reproducible subsample, read **while the dataset opens**. Every part of
Cellucid then works on those cells alone — the view, colours, filters,
highlights, analyses, gene expression, exports — and only the parts of the
file that hold them are read. Over a network that is less to download; on any
device it is less memory.

## Cells to open

**Session › Open data › Cells to open** holds one slider and the exact count it
stands for ("228,245 of 2,282,447 cells"). It says *All N cells* only when the
count is every cell: a subsample a few cells short of the whole (an opening
fitted near the device's limit) shows its share with as many digits as keep it
under 100 (99,996 of 100,000 cells is `99.996`), and **Apply** opens exactly
that count.

While a subsample is open, one short line under the settings says how its
cells were chosen ("Uniform, seed 7.", "Fitted to this device: blocks of 19,
seed 0."), or what needs your attention; with every cell open it says, when
they are, that genes are read for a sample of them ("Genes read for 40K of
2.3M cells.", {ref}`genes-on-a-sample-of-the-cells`). **Details**
under it opens a compact table of the facts behind that line: **Open** (k of
n cells), **Chosen by** and **Seed**, **Fitted to** (for an opening fitted to
the device: the memory it needed and what the device has), **Stored order**
and **Matrix read** (what these cells cost on how the file is stored),
**Genes read** and **Sized to** (the cells a gene is read for, and why: the
time a gene takes on this connection, or its gene index's memory, with how
the file is stored; for a MatrixMarket folder, the matrix read after its
cells, which holds genes for the cells it kept so far, none while it is read:
{doc}`17_open_r_objects`), what a gene costs on this file (worked out when Details
first opens), **Cell set** (the cells' digest) and **Categories** (how far a
field's categories are from every cell's). A value with more behind it is
underlined with dots: its tooltip says why.

1. Move the slider (or type the percentage, with a decimal point or a
   decimal comma: *12,5* is 12.5 %) below 100 %. The settings of the
   subsample appear: how the cells are chosen, its seed, and what the method
   needs (a field, an embedding). When the share you chose would read nearly
   the whole file anyway (its cells lie in nearly every stored chunk, and a
   chunk is read whole), the line says so before you apply it, with what the
   cheap method would read, and one button switches to it. On a file in its
   own order that is blocks ("Reads about 100 % of the file for 1.0 % of its
   cells; blocks read about 1 %.", **Use blocks**); on a file Cellucid wrote
   in priority order it is uniform with the file's own seed, since blocks of
   neighbouring cells lie throughout such a file (**Use the file's order**).
2. Press **Apply**. A larger or smaller share of the same method and seed
   changes the open cells in place ({ref}`raising-and-lowering-the-share`).
   Other cells are read again, in the same reader and through the parts of
   the file already open. The control is yours again as soon as the new
   cells are open: what the views showed comes back on them behind it (a
   view coloured by a gene reads that gene for the new cells, which on a
   large atlas takes a while), and another share can be applied meanwhile.
   Your
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
   named in a notice; the rest is kept. A text column whose new cells hold
   labels it had not shown lists them in their code-point places; every
   label it showed keeps its colour and whether it is hidden. A categorical
   with more categories than open cells and than one request reads (one
   category per cell, as barcodes kept as a categorical) lists the categories its open cells hold, and the
   new cells' join them in the file's order, the same way. What the view
   shows then is exactly what opening the dataset at that share directly
   shows, but for one thing: a new label whose own colour a label already
   shown holds takes the next free colour, where a direct opening may have
   given the two the other way round.
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
| **File default** (shown as *File default (uniform)* or *File default (blocks)*) | for a file Cellucid wrote (stored in priority order): uniform under the file's own seed; for any other file: blocks aligned to how it is stored | always the first choice: either way only about the chosen share of the file is read |
| **Blocks** | whole runs of neighbouring cells | a file Cellucid did not write: parts of it that hold only left-out cells are never read. On a file Cellucid wrote, blocks are scattered through it and read more (Details says so: *Stored order*) |
| **Uniform** | every cell equally likely | on a file Cellucid wrote, with its own seed (the default), it reads only its share; with another seed, or on any other file, its cells lie in most of the file's chunks, so it reads more (Details says how much: *Matrix read* "45 % for 10 % of the cells") |
| **Stratified** | every category of a field in proportion, or the same count of each (*Same count per category*) | rare cell types must stay visible, or each sample should keep its share |
| **Spatial** | cells on a grid over an embedding, dense regions thinned, sparse ones kept | the map's outskirts matter more than its crowded centre |

Details also shows how far the open cells' categories are from every cell's
(*Categories* "cell_type: within 0.3 percentage points (B cell)": the largest
difference between a category's share of the open cells and of every cell,
in percentage points), most useful with **Blocks**.

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
cells **in place**: nothing shown is read again. One exception: a share
whose genes are read for a sample of its cells and whose cells lie spread
through a file Cellucid did not write (**Uniform**, **Stratified**,
**Spatial**) opens in the file's own block order instead
({ref}`genes-on-a-sample-of-the-cells`), and another size of it is then
read anew (through the parts of the file already open).

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
the first cells of the order they open in (a sample of them, or, for a
share opened in the file's block order, its cells in the blocks read so
far), so the picture, the colours and filters are right from the first view;
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
its blocks), a dataset opens at once. A share also opens at once whenever
growing it would put all its cells on the map later than reading them at
once (over any link, each step costs a round trip): a share never takes
longer to show whole than every cell does. A later visit never opens fewer cells
at first than this browser keeps the parts of: the parts kept cost nothing
to read, and only the others are priced. A file Cellucid wrote opens
progressively whole or by any share (one that is not its own uniform share
opens in file order: its cells in the file's first rows first, a sample of
them at every step); any other file by a **Uniform** or **Blocks** share (a
file too large for the device opens so, fitted), by any share whose genes
are read for a sample of it (it opens in the file's block order, and each
step is its cells in the blocks read so far) -- each of these only while
growing would not end later than opening at once (above). Every cell of a
file Cellucid did not write opens at once: in its Blocks subsample's order
when its genes are read for a sample of them and that sample brings the
first gene sooner by more than the order costs
({ref}`genes-on-a-sample-of-the-cells`), else in file order with genes read
for every cell. Other stratified and spatial shares and a session's recorded
cells always open at once.

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

(why-one-gene-reads-every-open-cell)=
### Why one gene reads every open cell

How a file stores its expression matrix decides what one gene costs:

| Stored | One gene reads |
| --- | --- |
| cell by cell, sparse (**CSR**: most atlases, every CELLxGENE `.h5ad`) | one gene's values are spread over every cell's row, so the first gene reads the rows of **all** open cells once; later genes reuse the index this builds and read only their own values |
| row by row, dense (in one piece, or in blocks of whole cells) | one gene's values lie in every block, so the first gene reads the blocks of all open cells once; later genes reuse the index |
| dense in strips of genes too wide to keep decoded, whose index fits the memory kept for gene indexes | the first gene reads every strip once and indexes them in memory; later genes read the index |
| gene by gene (**CSC**, dense in gene strips: a Loom file's 64-gene strips, every file `cellucid.prepare()` or `cellucid_prepare()` writes) | only that gene's own values (its strip). A strip larger than the memory kept for decoded chunks is decoded anew for each of its genes: still a fraction of the matrix, never the whole of it |

The first gene's notice is one short line (what it reads, how long at the
speed Cellucid measured, how many cells) with **Cancel**; its **Details** say
how the file stores the matrix, why that makes one gene read every open
cell, what later genes read, and what avoids it: fewer open cells, or a copy
written gene-major (`cellucid.prepare()` in Python, `cellucid_prepare()` in
R). What it reads and how long count only the bytes still to download: what
this browser kept of the file from an earlier visit is read from this
computer, and the Details say how much (*Kept here*).

### Genes fill in as cells are read

The first gene of a matrix stored cell by cell needs the entries of every
open cell it is read for read once (every open cell, or the gene sample's:
its notice says how much, before it starts). Before that, Cellucid reads
where each of those cells' rows starts in the file (their row pointers,
and only theirs), once: a step of its own in the read activity. Over a link,
values stored compressed are read with the first gene, so the genes after it
read nothing more. With a subsample, or a gene sample,
the gene's colours appear before the read ends: first on about a
thirty-second of the cells it is read for, then on twice as many, and so on,
the cells not read yet shown in the missing-value colour. Each step is a sample of the
open cells, so the picture is right from the first colours; the notice says
how many cells are coloured so far and how long the rest of the gene will
take, and the read activity's bytes and time left are the whole gene's,
not one step's. A share whose cells lie spread through the file (a
**Uniform** share of a file in its own order, past a few tens of thousands
of cells) is read in one pass in the file's order instead: its cells lie in
nearly every stored chunk, so each step would read nearly the whole matrix
again, and one pass reads each chunk once (a 10 % share of a 10M-cell file
read its matrix twice over in steps). Its gene's colours then come once, at
the end of that pass, and its notice promises no earlier ones. Later genes
read only themselves.

Picking another gene while the first is still being read does not start
over: the new gene joins the read under way, and nothing read is read
again. A gene cancelled on its own stops the read; the next gene goes on
from where it stopped.

### The gene index is kept for later visits

The index a matrix stored cell by cell needs before its first gene is kept
in this browser (its private storage, like the files read whole: the
Session panel's *Stored* row says how much, and its **Manage…** dialog
lists each kept index and removes the ones you choose), under
the file, the matrix, and how the cells were chosen (method, seed and its
settings; never the share). Over a link, the gene values the first gene
read with it are kept beside it, so the cells it covers need nothing more
from the file. Opening the same file again later — after the tab closed, in
another tab, another day — with the same method and seed reads its genes at
once: the index of the cells it shares with an earlier opening, and their
values, are read from this browser, and only cells no earlier opening
indexed are indexed (and downloaded) now. A kept index is read into memory
at its first gene, in one read, when the memory kept for gene indexes holds
it (as an index built in this opening would be held), so analyses over many
genes run as fast as on the first opening; one larger than that memory is
read from this browser's storage in long ranges of the genes asked for
together. With **Uniform** and **Blocks** a larger share under the same
seed holds the smaller one's cells first, so it indexes only its new cells,
and a smaller one none. A share read in one pass in the file's order is
kept as one part and used whole: a later visit with the same cells reads
its genes at once, a larger share indexes only its new cells (one more
pass), and a smaller one indexes its own cells again (lowering the share in
the same tab reads nothing). **Stratified** draws each size apart, so two sizes
share only the cells their orders agree on from the start (in a test file of
100,000 cells in five groups, stratified shares of 10,000 and 20,000 cells
share their first 9,678): the rest is indexed again.

It is kept across visits only under the file's version, the one rule for
everything Cellucid keeps in this browser: a name of the file's content
where its source has one (a CELLxGENE Census release, a checksum or commit
on Google Drive, Figshare, Zenodo, Hugging Face or GitHub, a OneDrive
content tag), else its exact size with the `Last-Modified` date and the
strong `ETag` its host sends (the `ETag` a viewer host — `cellucid serve`,
a notebook, R — lists a file with; a local file's size and modification
time). Each opening asks the source for the version before anything kept
is used, and what is kept of other versions of the file is removed then. A
file whose source names none of these (no content name, no `Last-Modified`,
no strong `ETag`) keeps nothing for a later opening.

Before a kept part is used it is checked against its own record (its rows
are this opening's, its sizes and checksums are whole) and against the
file's row layout (a matrix stored sparse: the number of stored values in
each of its cells, read again; a part of a file whose rows changed is
removed). What a version cannot tell: a file replaced by other content of
exactly the same size within the same second (`Last-Modified` is kept to
the second), on a host that sends no `ETag` and whose source names no
content, is not noticed. How each source names its version:
{doc}`../p_developer_docs/09_data_loading_pipeline_and_caching`, *When a
dataset changes*.

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
| `block` | `blocks`: cells per block (default: the storage's block or a doubling of it, sized to the connection or kept from an earlier visit; the block used is shown and recorded, {doc}`../p_developer_docs/19_subsampling_specification`) |
| `by`, `alloc`, `min` | `stratified`: the field, categorical or text (a text field's values are its categories, in code-point order); `balanced` for the same count per category; a first minimum per category |
| `embedding`, `grid` | `spatial`: the embedding; boxes per axis |

A link's subsample is used as given, even where the device could hold more
or fewer cells. A malformed value is named in a message, and the dataset opens
as if the link named no subsample. After you Apply a subsample, the address
shows it (as its exact cell count), and **Copy link** under *Cells to open*
copies a link that opens exactly these cells. The link names how the cells
are chosen, never the cells themselves; Details under *Cells to open* shows
the cells' digest (*Cell set* "3f2a…"), so two people can confirm they see the
same cells.

## Datasets too large for this device

When opening every cell would not fit this device's memory — the per-cell
data the view holds (77 bytes a cell, as measured with levels of detail on,
against the memory allowance, half the device's memory), or the gene index of
a matrix every open cell holds (a modality's `X` of a MuData file; in its
share of what those cells leave of the allowance: the cells and the index
divide one allowance, never each the whole of it) — Cellucid opens a
subsample sized to fit, by the file's default method (uniform under the
file's own seed for a file Cellucid wrote, **Blocks** for any other), and
says so ("Fitted to this device: …"), with **Open all cells**; Details shows
the memory it needed and what the device has, and how far its categories are
from every cell's. The budget is what the browser reports of the device
(`navigator.deviceMemory` and the JavaScript heap limit, see
{doc}`../n_benchmarking_performance/index`). What a gene costs never makes
fewer cells open: it sizes the gene sample instead (below).

A dataset that fits always opens whole, and **Open all cells** or a larger
share reads more cells anyway; datasets opened next still open fitted.

(genes-on-a-sample-of-the-cells)=
## Genes on a sample of the cells

The map of an opening and its genes cost per cell on very different scales.
A cell on the map is a few bytes: its position and the codes of its
categories. A gene of a matrix stored cell by cell (CSR, every CELLxGENE
`.h5ad`; or dense row by row) costs every open cell's whole row, read once
({ref}`why <why-one-gene-reads-every-open-cell>`), and its entries kept in
the gene index. So such a dataset opens **every cell the device holds** on
the map, with all its categorical fields, and a gene is read for a
**sample** of the open cells:

- **How many.** As many as the first gene reads in **about 20 seconds** on
  this connection, and as many as its gene index fits in memory beside the
  open cells, whichever is fewer. The first gene is read for the cells its
  read is priced to take **about 10 seconds** for, at the round trip and
  speed Cellucid measured so far — the stored chunks their rows lie in,
  their row pointers, and the requests they take, for the matrix genes are
  read from (`raw.X` rather than `X` when `raw` holds more genes than `X` or
  `X` is scaled; Data mapping › Expression matrix) — and drawn; the rest of
  the 20 seconds then reads genes for more cells, as many as the first
  gene's own read time says fit (the legend counts them; **Stop** keeps what
  is read). The time a read really took sizes it, not a price: on a first
  visit to a host its first, small reads measure a fraction of its speed
  (CELLxGENE atlases over Amazon S3: priced at 20 seconds, read in 6–15),
  and a connection's speed changes. Twenty seconds is a response-time target
  for a requested operation shown with its progress (people keep their
  attention on a task for about ten seconds without feedback, hence the
  first picture within ten, and wait about twice that with a progress line),
  not a limit of the browser or the file. A sample bounded by time is
  taken only when it brings the first gene sooner by more than opening
  every cell in its Blocks order costs the first picture end to end (the
  order is timed as it is made, and every pass over the cells it adds --
  its copy to the page, the page's check, the first view's arrays put in
  its order -- is priced at that time); a sample of nearly every cell is
  not, and every cell's genes are read in file order. A file on this
  computer, or one a server on this computer serves (`cellucid serve`,
  Python's and R's viewers, any web server at `localhost` or `127.0.0.1`),
  is priced at no less than a disk's speed (1 GB/s), whatever its reads
  measure: its few small opening reads measure only a fraction of it, and
  a read slowed by the computer's own work (another program starting
  beside it) measures the computer's load, not the link. A tunnel to a
  remote server (`ssh -L` to a `localhost` port) is priced the same way. So, as on a fast link,
  its genes are read for every open cell, unless the gene index does not
  fit in memory or the first gene would take past ten seconds even at that
  speed; and a matrix stored gene by gene (CSC, dense in gene strips, every file
  `cellucid.prepare()` or `cellucid_prepare()` writes) has **no sample**: a
  gene is read for every open cell. A notebook's viewer is priced as its
  connection measures (its kernel may run across the internet).
- **Which.** The first cells of the open cells' order. Every cell of a
  file Cellucid did not write being open, that is the order of its
  **Blocks** subsample (whole blocks of the file's storage, so no chunk is
  read for nothing; blocks large enough that each one's download outlasts
  its round trips, so the 20 seconds read as many cells as the connection
  carries — measured on CELLxGENE atlases, this keeps the cell types' and
  donors' shares as close to the whole as blocks of the storage's own size
  do). The sample of *m* cells is then exactly the Blocks subsample of *m*
  cells at that block size: a gene index kept on an earlier visit counts — a
  later visit opens the blocks it was kept for, whatever the connection
  measures then — and later visits read genes for more cells within the same
  20 seconds (*Kept here* in the Details names what was kept). Coming back
  is never slower than the first visit: the first gene is drawn at once on
  the cells whose index was kept, downloading nothing, and genes for more
  cells are then read in the background within the 20 seconds (stopping for
  any gene you pick meanwhile). The first of those steps reads the cells the
  connection is priced to deliver — as this visit has measured it, or as
  this browser remembers the host; when neither is known (the opening
  downloaded nothing, and the host's speed was forgotten or never kept) a
  quarter more cells, never all at once as if downloading them were free —
  and the steps after it are sized by what the reads took. A
  **Uniform**, **Stratified** or **Spatial** share of such a file opens in
  the same block order, its cells filtered: its gene sample is its cells in
  the blocks every cell's sample reads first, so it never downloads more,
  nor asks for more, than every cell's would in the same time. In its own
  order each of its first cells lay in a part of the file of its own (a 1 %
  share of a 10-million-cell file read its first gene for 71 cells in 10.7 s
  where every cell's read 1.5 million in 5.5 s); the cells are exactly the share's either way — a session, a
  link, Python's and R's `subsample_cells` give the same ones — and only
  their order changes. A Blocks share keeps its own order (its first cells
  are whole blocks already); one chosen without a block size takes every
  cell's, so its genes cost no more a cell than every cell's. A share of a
  file Cellucid wrote opens in file order (its own uniform share is the
  file's first rows). Every gene is read for the same cells, from the index
  the first one builds. Cells open before in this tab and opened again with
  **Apply** (a share lowered, then raised back to every cell) come back in
  their order: the 20 seconds' sample is sized again for the connection as
  measured then, but never to fewer cells than were read before, nor to
  every cell unless you read them all; the count you read with **Read for
  all cells** comes back as it was. A share chosen while every cell was open
  reads no more than every cell's sample did.
- **The other cells** have no value for the gene: they are drawn faded,
  almost to the background and never in the gray of 0, beneath every cell
  read. The cells read are drawn on top and larger, sized from how dense the
  sample is on screen, so the pattern shows at once; their zeros are gray
  beneath the expressing cells (a gray made from the background: light on a
  light background, mid-dark on a dark one, {doc}`../d_fields_coloring_legends/03_color_by_behavior`).
  The legend names the three states in one
  line, **Faded: not read · gray: 0**, each beside a dot of its colour, and so does an exported figure. The gene's legend says **Gene read for 40K of 2.3M cells** (the
  exact counts in its tooltip: the line fits the narrowest panel), with an (i)
  that says this, and **Cells to open** says *Genes read for …* with **Genes
  read** (exact) and **Sized to** in its Details.
- **Read for all cells** (in the gene's legend) reads the rest. Its cost —
  the download and, at the measured speed, the time — is said beside it
  before you press it. It reads in steps that double the cells read so far
  — the first reads about as many cells as the first gene did, so it takes
  about as long, and each later step twice the one before — and nothing read
  is read again; the gene shown fills in after each step, and its notice
  says how many cells are read, with **Cancel** (the button becomes
  **Stop**, with no cost beside it: the cost said is that of **Read for all
  cells**, which is shown again once the read, or the gene time's own top-up,
  ends). What a stopped read read stays read; a kept view showing the
  gene takes the values read and is drawn from them when you focus it. Once
  every open cell's genes are read, the line goes; later genes then read
  every open cell.
  When the gene index of every open cell fits neither this device's memory
  nor the storage this site may use (the full Human Lung Cell Atlas: about
  15 GB), the action reads genes for as many cells as it does fit, and says
  so (**Read for 640K cells**, the reason in its tooltip); when no more fit,
  it says *no room for more here*.
- **Filters** on a gene judge the cells it was read for: the others stay
  shown (the filter list says *cells not read kept*).
- **Analyses** that read genes (Differential expression, Marker genes, Gene
  signature, Dot plot, Detailed and Correlation on a gene) compute on each
  group's cells among the sample, and say so on their chip (*genes on
  40,000*, with each group's share in its tooltip); a group with fewer than
  two such cells is named with how many it has. Run them again after Read
  for all cells to use every cell. Analyses of categorical and continuous
  fields (Composition, a field's Detailed view) use every open cell.
- **Sessions** record the sample, and a viewer in a notebook hears it with
  the gene it colours by ({ref}`sessions-and-the-gene-sample`).

On the Human Lung Cell Atlas's 2.3 million cells over a home connection this
opens the whole map with its cell types in seconds, and reads genes for a
sample of tens of thousands of cells within 20 seconds; before, the whole
opening was shrunk to the few thousand cells whose first gene read in that
time.

The speed this is measured at is, until the opening's own reads have
measured as much, what the host delivered on earlier visits (an opening's
first requests are small, a round trip each, and measure a host slower
than its gene reads then run: HLCA's opening measured about 640 KB/s,
its gene reads 1.5–3.4 MB/s). It is kept in this browser per host (the
host's address and its speed, nothing of the data; at most 32 hosts;
listed in **Session › Open data › Stored › Manage…** as *Connection speeds
remembered*, whose **Remove** forgets them), scaled down when this opening's round
trip is longer than the one it was measured on (another network), and
never below what this opening measured.

## Neighbour graphs

A subsample keeps the edges of the neighbour graph between two open cells;
an edge from an open cell to one left out cannot be drawn. Once the graph is
read, the Dataset box says how much of the open cells' graph that is, beside
its edge count ("10 % of the open cells' edges; the others lead to cells not
open"): with Uniform at 10 %, each neighbour of an open cell is open with a
chance of about a tenth, so about a tenth of the edges remain.

Reading a share's graph never costs more than reading every cell's: the open
cells' rows are read in the file's order, either alone or, when they lie in
every stored chunk of the graph (a uniform share of a graph stored in chunks
of thousands of cells), by streaming the graph once, whichever downloads
less. A share in blocks reads about its share of the graph; a uniform share
of such a file reads all of it (each chunk holds some open cells), in about
the time every cell's graph takes, or less.

(analyses-on-a-subsample)=
## Analyses, filters and overlays on a subsample

Everything works on the cells open, as if the file held only them:

- **Analyses** (Differential expression, Marker genes, Gene signature,
  Correlation, Composition, Dot plot, Detailed, Quick Insights) compute on
  the open cells: a page, category or group is its open cells, and the
  counts they show are of open cells. A group with fewer than two open
  cells is not tested; a category with no open cell has nothing to test or
  count (Differential expression and Marker genes do not offer it, and
  Composition leaves a cell type with no open cell out of its test); a gene
  no open cell expresses shows no difference (p = 1), and its correlation is
  undefined (*constant in these cells*), never r = 0. A first answer on a
  sample ({doc}`../h_analysis/02_analysis_ui_overview`) is a sample of the
  open cells, and the final answer is exactly a run on every open cell.
  Per-cell exports name each cell by its row in the whole file.
- **Every result says how the open cells were chosen**, as a chip beside it
  (*10% of cells · uniform*). A **uniform** subsample, or a **stratified,
  proportional** one with no minimum per category, is a random sample of the
  file at one rate: its proportions and p-values estimate the file's.
  **Blocks** (whole runs of neighbouring cells: the default for a file
  Cellucid did not write), **stratified, balanced** (or with a minimum per
  category), **spatial**, and the first chunks of a subsample still opening
  are not: a file stored by sample or cell type loses whole samples or types
  in blocks, cells of one block are alike, and the other methods give rare
  categories or sparse regions more of their cells on purpose. Their chip is
  a warning, and its tooltip says what that means for the result; open every
  cell, or a uniform subsample, to analyse the file. **Composition** refuses
  to run when the open cells are stratified at unequal shares by the field
  it compares (its proportions would be the subsample's design).
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
  is run again; so is one that read genes when genes are then read for more
  of the open cells (**Read for all cells**).
- **Genes read for a sample of the open cells**: an analysis that reads
  genes computes on each group's cells among that sample, and its chip says
  so ({ref}`genes-on-a-sample-of-the-cells`).

(fewest-cells-open)=
### At least two cells

A subsample opens at least **2 cells**. The limit is not chosen: the viewer
draws an embedding only when its coordinates vary, and one cell has no
coordinate variation (the rule a one-cell file is refused by:
{ref}`embedding-refused`, *“no coordinate variation”*). So *Cells to open* does
not offer a share of one cell (its count line says why and **Apply** stays
off), a link with `cells=1` is refused with that reason, and auto-fit never
fits a dataset to fewer than two cells. A share is never turned into another
count: one whose exact count on the dataset it opens is fewer than two cells
(a share chosen on a larger dataset, `cells=0.01%` in a link) is refused with
the same reason, naming the share and its count. Two or more open cells that all share one point are refused by the
embedding's rule.

`cellucid.subsample_cells()` (Python) and `cellucid_subsample_cells()` (R)
select cells, and keep any count from 1: one cell is a selection, just not
one the viewer can draw.

## Sessions

A session saved on a subsample records exactly which cells were open (the
cells themselves, compactly, and how they were chosen). Loading it opens the
dataset with those same cells — reading it again first if it is open with
others — so every highlight and label lands where it was. A session saved on
every cell opens every cell.

(sessions-and-the-gene-sample)=
### Sessions and the gene sample

A session saved while genes were read for a sample of the open cells
({ref}`genes-on-a-sample-of-the-cells`) records that sample too: how many of
the open cells, the file's cell order they were open in (every cell's, or a
share's spread through the file), and the digest of those cells. Loading it reads genes for the same cells: when
the dataset opens again it opens with that sample; when it is open already
with the same cells in the same order, it reads genes for more cells if the
session read them for more (with the cost said and **Cancel**); with the
cells in another order it reads them again in the session's order. A session
saved with every open cell's genes read records nothing of it, and opens
with the sample this device and connection fit (**Read for all cells** reads
the rest). A viewer in a notebook hears the sample with the gene it colours
by: the `color_by` event's `cells_read` is how many of the open cells the
gene is read for (the first ones), or null when every open cell's are. Applied in Python or R
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

Both take the object's rows (columns in R) as the cells' original positions.
An object read from a Cellucid-ready file holds its cells in priority order,
not that one, and both refuse it until `cellucid.restore_original_order()` /
`cellucid_restore_original_order()` puts the cells back, as pulling results
and applying a session do: `adata[kept]` would otherwise take other cells.

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
