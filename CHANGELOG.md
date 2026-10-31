# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.6]
<!-- CELLUCID_VERSION -->

Version 0.9.6 is a patch release on 0.9.5. Files prepared, sessions saved and
dataset catalogs written by 0.9.5 are read as before. The viewer protocol is
schema 2: use the Python and R packages of this release with the web app of
this release. The main changes are:

### Highlights

- **Subsampled reading.** Open a seeded share of a dataset's cells; only the
  parts of the file that hold them are read, and every part of Cellucid
  (view, filters, analyses, genes, exports, notebooks) works on those cells.
  The same cells in the web app, Python and R, in a link and in a session.
- **Files Cellucid writes store their cells in priority order**, so a share
  of the cells is the first rows of every array.
- **Datasets too large for the device open fitted, and large ones open
  progressively**, growing in place; a larger or smaller share is applied in
  place, reading only the cells added.
- **One notebook widget for Python and R** (`show_anndata()`, `cellucid.Widget`,
  R `cellucid_view()` in Jupyter): no port, no proxy.
- **Round trip in one call**: `viewer.pull(adata)` / `cellucid_pull()` write
  every highlight, label, score and gene set back onto the object;
  `push()` / `cellucid_push()` send yours in.
- **Up to 100 million cells in 3D** on a 16 GB laptop, with levels of detail
  that are uniform samples of the cells.

### Web app

#### Subsampled reading

- **Cells to open** (Session › Open data): one slider (or a typed share) with
  the exact count it stands for ("228,245 of 2,282,447 cells"). Below 100 %
  the method, its seed and what the method needs appear, with **Apply**.
  Datasets opened next in this tab take the same share; a method a dataset
  cannot use (a stratifying field or embedding it lacks) opens the same share
  by the file's default method, and says so.
- **Four methods, one specification.** `uniform` (every cell equally
  likely), `blocks` (runs of neighbouring cells aligned to how the file
  stores them, so parts holding only left-out cells are never read),
  `stratified` (proportional or the same count per category of a categorical
  or text field, optional minimum per category) and `spatial` (balanced over
  a grid on an embedding, keeping sparse regions). Every method is seeded;
  with `uniform` and `blocks` a smaller share is inside a larger one. The
  specification is published, with a conformance fixture shared by the web
  app, Python and R.
- **The file's default.** A file Cellucid wrote opens a share as `uniform`
  under the file's own seed (the first rows of every array); any other file
  as `blocks`, whose block size is set by the file's layout alone (about four
  chunks of its dominant per-cell array, or 64 KiB of it when it is read byte
  by byte), so the same file and share always open the same cells. The panel
  says what `uniform` would cost on such a file ("fetches 45 % of the matrix
  for these 10 % of the cells") and how far a `blocks` share's categories are
  from every cell's.
- **At least two cells.** One cell has no coordinate variation to draw (the
  rule a one-cell file is refused by): the control does not offer a share of
  one cell, a link asking for `cells=1` is refused with that reason, and
  fitting to the device never goes below two.
- **In a link.** `?data=…&cells=10%&sample=blocks&seed=0` (also `by`,
  `alloc`, `min`, `embedding`, `grid`, `block`) opens exactly that share,
  overriding fitting to the device; `cells=100%` opens every cell. Parsed
  strictly: a malformed value is named and the dataset opens as if none were
  given. A share chosen in the panel is written to the address as its exact
  count, method and seed (never one fitted to the device, never the cells),
  **Copy link** copies one, and the panel shows the cells' digest so two
  people can confirm they see the same cells.
- **Sessions keep their cells.** A session saved on a share records the cells
  (compactly), how they were chosen and their digest; loading it opens
  exactly those cells, reading the dataset again first when needed. Per-cell
  session data names cells by original position, so it applies whatever
  order a file stores its cells in.
- **Counts say "k of n".** The Dataset box ("228K of 2.3M cells"; Details:
  the method and seed), the filter line (*Showing 120 of 240 of 1,200
  cells*), categorical legends, the category builder, the screen-reader
  status line and the dataset picker all say how many of every cell are open.
- **Neighbour graphs on a share** keep the edges between open cells; the
  Dataset box says how much of the open cells' graph that is, and Data
  mapping offers every graph on the open cells.
- **A figure of a share says so**: SVG and PNG provenance name the open
  cells, how they were chosen and their digest.

#### Too large for the device, progressive opening, Apply in place

- **Fitted to the device.** When every cell would not fit this device's
  memory (77 bytes a cell, as measured with levels of detail on, against half
  the device's memory; or, for a matrix stored cell by cell, the index that
  makes its genes readable), a share sized to fit opens by the file's default
  method, with a notice and **Open all cells**. About 104 million cells fit a
  16 GB machine.
- **Progressive opening.** When reading the cells asked for is measured to
  take more than about a second (the link's speed and latency, the stored
  bytes of the embedding's parts the first view needs and their decoding, and
  how long this device took to show cells before), and a first view would
  arrive at least twice as soon, a first view of fewer cells shows at once and
  grows in place, doubling, to the cells asked for. Every step is a sample of
  them. A file stored in parts much coarser than its blocks can grow part by
  part, each part read once, ending exactly as opening the share directly.
  The same rule prices a file picked from this computer.
- **Apply in place.** A larger or smaller `uniform` or `blocks` share under
  the same seed (and every cell of a file Cellucid wrote) changes the open
  cells in place: raising reads only the added cells (their positions, their
  values of every loaded field and gene, the gene index for them alone),
  lowering reads nothing. Camera, every view's field, filters and settings,
  highlights, labels, pushed columns, hand-hidden cells and analysis windows
  stay; highlights and labels made from a column take the new cells as
  opening them would. Each per-cell structure ends bit for bit what opening
  at that share directly holds. On the graphics card only the added cells are
  sent, and lowering cuts the end. Another method or seed reads the dataset
  again through the parts already open, carrying the same state by original
  cell position.
- **Nothing read twice.** Another share of the same file, in the same or a
  new reader, takes what was already fetched and decoded (within the tab's
  memory budget); Zarr folders on a web host, picked folders and file sets
  are reused while their version is unchanged.

#### Reading at the floor

- **A share reads only its cells' parts**, element by element (obs of every
  encoding, embeddings, vector fields, graphs, `X` as CSR, CSC or dense,
  layers, `raw`, stored outlier quantiles, MuData modalities): a permanent
  test holds every element, at every share, to the stored chunks that hold a
  kept cell, each fetched once. Opening 1 % of CELLxGENE's 2.3M-cell Human
  Lung Cell Atlas: first view 3.8 MB in 143 requests at 8.2 s (was 20.3 MB in
  2,645 requests), at the floor computed from the file's layout.
- **The first view never waits for gene names.** The gene list (the var
  index, or the column Data mapping names) is read right after it, and no
  other var column; meanwhile the gene search, gene-based analyses and their
  pickers say *Reading gene names…*, and sessions, exports with genes and
  data-mapping changes wait.
- **Faster reads over a network**: HDF5 chunk-index nodes no longer read
  before a gene (the lung atlas read 12,382 of them first); index nodes of
  one task are read as one burst (all HDF5 formats; 100 Loom genes 11.1 →
  7.8 s over 150 ms / 5 MB/s); every S3 address of a public bucket is used;
  neighbouring reads are merged when that saves time; HTTP/1.1 hosts are
  recognised; sharded Zarr v3 reads all its shard indexes at once, then all
  inner chunks; uncompressed Zarr chunks are read in part; two columns of a
  latent space read only the open cells' rows.
- **Stricter, steadier reads.** A host answering another byte range, a file
  whose size, `ETag` or modification time changed while open, is reported,
  never read into what was read before. Busy hosts (429, 502, 503, 504) are
  waited for with spaced retries and fewer requests at once. A failed gene is
  read again on the next attempt. A shared read is cancelled only when every
  reader waiting on it has withdrawn.
- **Decoders.** Safari reads Zstandard-compressed Zarr about seven times
  faster (55 → ~395 MB/s); a Blosc frame with the delta flag is now decoded
  as c-blosc2 writes it (it was misread).

#### Genes on every path

- **Genes leave the viewer by identifier** (the var index) everywhere: gene
  sets, the gene coloured by, GMT files, DE and marker CSVs, `.h5ad`
  exports. When Data mapping names the genes by a var column (a CELLxGENE
  file's `feature_name`), the name shown goes alongside; named by the var
  index, nothing does.
- **The first gene of a cell-major matrix says what it costs** before it is
  read (bytes and time at the measured speed), with **Cancel**; on a share,
  its colours fill in as cells are indexed, each step a sample of the open
  cells. No silent rescans; errors say what failed.
- **The gene index is kept across visits.** For a matrix stored cell by cell
  (CSR, or dense row-chunked) the gene index is kept in this browser under
  the file, the matrix and how the cells were chosen (never the share):
  opening the same cells again, another day, reads genes at once, and a
  larger `uniform` or `blocks` share indexes only its new cells. Kept only for
  a file whose version is known; checked before use against its record, the
  file's row layout and, unless the version names the content, a fingerprint
  of the file; each gene's piece carries a CRC-32.
- **Fixed:** the first gene of a share of a Cell Ranger `.h5`; a text field
  for stratifying; SingleCellExperiment graph weights (were all 1); a
  MatrixMarket first gene or cell that looked like a header (`IDH1`,
  `cell0`).

#### Files that cannot be read in parts

- **Said up front, kept for next time.** An R object, a MatrixMarket folder
  or any file on a host that serves no parts of files: the opening says it is
  read whole, how much, and that a Cellucid-ready AnnData is read in parts
  instead. What it is read into is kept in this browser's storage across
  visits (for a file with a known version) and opens again in milliseconds,
  reading nothing; converted stores are written in priority order with a
  gene-major `X` (a gene of a 50,000-cell Seurat object: 7 ms). **Session ›
  Open data › Stored** counts what is kept, with **Clear**.
- **Checked as it is read.** Everything kept carries a CRC-32 of every
  16 KiB; damaged bytes fail the read as damaged, the copy is removed and the
  next opening reads the file again.
- **Safari keeps converted files and gene indexes across visits** like
  Chromium and Firefox; only private windows keep them in memory for the tab.

#### Rendering at scale

- **Levels of detail are uniform samples.** A level is the first cells of
  the open cells' priority order (1,000 cells, then a quarter more at each
  level), so a reduced level shows the cloud's density and mix of categories
  as all cells do. Cells already in that order (a file Cellucid wrote, a
  share) draw every level at once (100M 3D: 0.33 s after levels of detail
  were turned on, was about 29 s); the spatial order culling uses is taken a
  few million cells a frame, without stopping a frame.
- **100 million cells in 3D.** Levels of detail at 100M 3D were silently
  unused (a GPU texture Chromium refused); fixed, with a notice if a GPU ever
  refuses. Positions and colours live in segments of 16,777,216 points, so
  no WebGL buffer limits the cell count; an array the browser will not
  allocate in one piece is made in WebAssembly memory (3D up to 357.9M
  cells); the velocity overlay's textures are layered. Memory per cell
  measured at 100M 3D: 54 bytes drawn, about 50 after levels of detail are
  built (was 65); 150M cells in 2D fit 16 GB.
- **Keep view and recolouring above 16.7M cells on Chromium** work (a
  vertex attribute released without a buffer kept stale bindings); released
  kept views free their GPU memory; GPU errors name the call that failed.

#### Analyses on a share

- **First answers on a sample.** Differential expression, Marker genes, Gene
  signature, Correlation, Composition, Dot plot, Detailed and Quick Insights
  show a first pass on a uniform sample of the open cells when a run is
  predicted to take more than about a second, then refine in place to every
  cell; the final result is bit for bit a single run, and only it is saved
  or exported. A first pass reads only its sample's part of each gene where
  the file allows.
- **Analyses follow the open cells.** When the open cells change, results are
  marked **Stale · Run** (its tooltip gives both counts); a stale result is
  not saved in a session, not exported and withdrawn from notebooks until it
  is run again. Differential expression and Marker genes no longer kept the
  earlier cells' groups after an Apply.
- **Stored outlier quantiles hold on a share** (they describe the whole
  dataset); quantiles computed on the open cells say they are an estimate.
- **Sparse statistics.** Differential expression and Marker genes read each
  gene's stored values, never a dense column (100M cells at 10 %: peak
  2.4 GB → 0.6 GB).
- Every analysis is checked against an independent reference on every way
  cells open (each method, progressive, Apply, sessions, served data) and
  every format and source.

#### Round trip and exports

- **Viewer protocol schema 2** (shared fixture re-pinned in all three
  repositories): every cell either side names is an original position; new
  events `result`, `result_removed`, `results`, `color_by`; new commands
  `setField`, `setGeneSet`; `ready` carries `n_obs`. Results carry an epoch,
  so a reloaded page replaces an older page's results.
- **Exports without a notebook**: the highlighted cells as an obs-name CSV
  or an `.h5ad` subset (written in the browser; anndata, h5py and Cellucid
  read it), any obs field as a CSV (**Export** beside **Duplicate**), gene
  sets as GMT. Exports list cells in the dataset's own order; genes by
  identifier, with the names shown beside when a var column names them; DE
  and marker CSVs name genes by identifier.
- **From a notebook**: pushed columns and gene sets show as fields and under
  **From notebook** in Gene signature and Dot plot, and are saved with the
  session.

### Python

- **`cellucid.subsample_cells()`** gives exactly the cells the viewer opens
  for a method, size and seed (any count from 1; the viewer opens at least 2).
- **`prepare()` stores cells in priority order** (`priority_seed=0`;
  `None` keeps the order given) in `.h5ad`, Zarr v2/v3, `.zarr.zip`, AnnData
  and MuData: `uns['cellucid']['priority'] = {seed, specVersion, nObs,
  levels}`, `obs['cellucid_original_position']` and, beside a sparse `X`,
  `varm['cellucid_prefix_counts']`. Chunks are smaller for shares: per-cell
  arrays 2^16 values, `X` data and indices 2^14. Preparing an object read
  from such a file puts it back in its original order first.
- **`cellucid.restore_original_order()`** puts an AnnData or MuData read from
  such a file (also a backed one) back in the order its cells were given. A
  record that no longer describes the cells (a subset, a reordering) is
  never trusted.
- **`show_anndata()` returns a notebook widget** (`cellucid.Widget`, an
  anywidget; `cellucid.Widget(adata)` makes one without showing it). The page
  and the viewer channel travel over the notebook's own connection: no port,
  proxy, `client_server_url` or `allowed_hosts`. Tested end to end in
  JupyterLab. Calls that wait for the viewer raise at once inside the cell
  that would wait, saying what to do. `anywidget` is a dependency;
  `jupyter-server-proxy` is no longer one. `AnnDataViewer` is the
  browser-tab viewer (its unused `height` is gone).
- **Served in priority order.** `show_anndata()`, `serve_anndata()`,
  `AnnDataViewer` and `cellucid serve` serve an in-memory AnnData, and any
  AnnData, MuData, Loom or 10x file `prepare()` did not write, as a Zarr store
  with the cells in priority order, made chunk by chunk as the viewer asks: a
  share transfers only its cells (1 % of a 1M-cell object: 0.18 MB, was
  2.1–10.6 MB). The temporary `.h5ad` is gone. An in-memory object is
  snapshotted before the call returns (memory-mapped, about 2 s per 100M
  values), its matrices stored as `prepare()` stores them. Files are served
  at disk speed (CSR `X` 7 → 109 MB/s). Each served file is listed with a
  strong `ETag`, so the browser keeps what it read across visits and sees a
  rewrite.
- **Round trip.** `viewer.pull(adata)` / `cellucid.pull(adata, viewer or
  paths)` writes every highlight, label, score and gene set onto `obs` /
  `uns` in one call (`<NA>` for cells a subsampled viewer did not open), from
  a live viewer or from the web app's CSV, GMT and session files;
  `viewer.push()`, `push_gene_set()`, `subset()`, `gene_sets()`, `results`,
  `selection` (obs names), `write_gmt()` / `read_gmt()`. Pull and subset
  refuse an object whose cells are in another order or that lacks the
  viewer's cells.
- **One session API.** `request = viewer.request_session_bundle()`, then
  `request.result()` or `viewer.apply_session_to_anndata(adata, request)`.
  `get_session_bundle()` is removed. Sessions saved on a share apply (cells
  not open are `<NA>`); `set_data_mapping()` no longer waits or returns.
- Viewers of a `prepare()` file name cells correctly in Python (obs names in
  original order); count messages read "1 cell".

### Docs

- **Fourteen runnable round-trip notebooks**, seven Python–R twins with
  "Open in Colab" badges: show a dataset in any notebook, pull results back,
  push clusters and gene sets, read the web app's exports, re-analyse an
  exported `.h5ad`, subsampled views; a gallery page lists them side by side.
- New pages: **Subsampled loading** (web app) and the **Subsampling
  specification**; the Cellucid-ready AnnData specification gains **Priority
  order**; the notebook pages describe the widget; the published-data
  tutorials open CELLxGENE's 2.3M-cell Human Lung Cell Atlas.
- Every page was checked against this release's behaviour (served stores,
  kept stores, sessions, analyses on a share, rendering limits).

## [0.9.5]

Version 0.9.5 is a patch release on 0.9.4 that breaks with every earlier
release: no export, session file or dataset catalog written for 0.9.4 is read.
**AnnData is the one data model.** The prepared export format is removed
from the web app, the Python package and the R package, with the export
server, the manifests, quantization and the `exports-v3` layout of the sample
datasets. `prepare()` and `cellucid_prepare()` now write a
**Cellucid-ready AnnData**, a standard AnnData that anndata reads and the web
app opens directly, and take a required `format=` argument. Session files are
**format 2**; a format-1 file is refused. Dataset catalogs (`datasets.json`)
are **schema 2** and list AnnData files. The contract between the writers, the
web app and the dataset repositories is the
[Cellucid-ready AnnData specification](https://docs.cellucid.com/en/latest/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.html)
(`docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md`).
The main changes are:

- **One data model**: `.h5ad`, AnnData Zarr v2 and v3 folders and
  `.zarr.zip`, read by the web app's own HDF5 and Zarr readers, and MuData,
  Loom, 10x Cell Ranger, MatrixMarket folders and Seurat or
  SingleCellExperiment `.rds`/`.rda` files, each opened as AnnData; one input
  format registry names them in the web app, Python and R. The Python and R
  writers produce identical Zarr stores from the same input.
- **Open from a link**: CELLxGENE, GEO, GitHub, Zenodo, Figshare, Hugging
  Face, Dropbox, Google Drive, OneDrive and any web address, read in place with
  range requests that follow the measured link, from one address bar that also
  connects to a server; a private Google Drive or OneDrive file opens for a
  collaborator after one sign-in, with no Cellucid server involved.
- **`cellucid serve` and Jupyter hand the file to the browser**, which reads it
  itself; the server computes nothing per request. They take every format the
  web app reads (AnnData, MuData, Loom, 10x, MatrixMarket, R objects) and
  in-memory AnnData and MuData.
- **One viewer channel**: a notebook viewer and Python talk over one
  authenticated WebSocket; the notebook output is a plain iframe.
- **MuData**: `prepare(data=mdata)` writes a Cellucid-ready MuData, and the
  `modality` launch default names the modality a dataset opens on.
- **Stored outlier quantiles** and an outlier filter that reads or computes
  them on first use, cancellably.
- **Analysis, selection, figure export and edges off the main thread**, every
  analysis cancellable.
- **Sessions** save analysis results, Figure Export settings and the theme,
  in smaller files that restore in short tasks.

### Added

- **Python**: `prepare(modalities=[...])` chooses which modalities of a MuData are written, in that order, as cellucid-r's `cellucid_prepare(modalities = )` does: several are a Cellucid-ready MuData of those modalities (the root `var`/`varmap` hold theirs alone); one is a plain Cellucid-ready AnnData with the root's cells, obs, embeddings, graph and identity and that modality's `X`, `var` and `obsm` (it must have measured every cell; `gene_expression`, `var` and `gene_identifiers` then apply to it).
- **R**: `cellucid_apply_session()` writes a saved session's highlight groups and user-defined fields onto a Seurat object or an experiment, with the rules of `apply_cellucid_session_to_anndata()`; both are held to one web-app-written session file.

- **Web app, Python, R**: The Cellucid-ready AnnData. `uns['cellucid']` names the dataset: `schema` (1), `id` (the id sessions and community annotations key on), `name`, `description`, `citation` (`text`, `url`), `source` (`name`, `url`), the launch `defaults` (`embedding1d`/`2d`/`3d`, `knn`, `expression`, `geneIdColumn`, `latentKey`, `vectorFieldDefault`, `obsKeys`, `normalizeEmbeddings`, `centroidMinPoints`) and `outlierQuantiles`. Every key is optional; an unknown key refuses the file with a message. A server's explicit launch choices override the file's defaults, and the dataset panel shows the citation.
- **Python, R**: The writers store the outlier quantiles. With `latent_key`, `prepare()` and `cellucid_prepare()` compute every categorical obs column's latent-space outlier quantiles (categorical, text, boolean and nullable boolean columns; a column with no category of `centroid_min_points` cells, 10 by default, gets no entry) and store them as `uns['cellucid']['outlierQuantiles'] = {latentKey, minPoints, fields: {<column>: float32[n_obs]}}`, NaN where a cell has none. Python (`cellucid.prepare_data._outliers`), R (`src/latent_outliers.c`) and the web app (`data/latent-outliers.js`) are bit-identical, checked against a numpy conformance fixture. On Suo (562k cells, 14 categorical fields, a 24-column scANVI latent space) they add 25.5 MB and take 3 s in `prepare()`.
- **Web app**: The Session panel's **Link or server** address bar, and the `?data=` page parameter, open data where it lives, with no sign-in: a CELLxGENE dataset (its id or Explore link, read from the CZ CELLxGENE Census copy of the latest release), a GEO series (`GSE…`, its `.h5ad` supplementary file), GitHub (`gh:owner/repo/path`, `github.com` tree, blob and Raw links, `raw.githubusercontent.com`, Git LFS files read from `media.githubusercontent.com`), Zenodo (`zenodo:<record>`, record pages, file links and DOIs; a concept DOI opens the latest version), Hugging Face (`hf:<user>/<repo>[@rev][/path]`, pinned to the commit) and Dropbox file links, or any web address of an `.h5ad`, a `.zarr.zip`, a Zarr folder or a catalog folder. A GitHub link without `@branch` reads the repository's default branch, whatever it is called. Single files are read in place with HTTP range requests in the reader worker, so only what a view needs is downloaded. Each provider is one module behind one contract (`assets/js/data/sources/providers/`).
- **Web app**: The **Folder…** item of the Open… menu, and a folder dropped on the page, open an AnnData `.zarr` store folder, as `prepare()` and `cellucid_prepare()` write it, a Cell Ranger `outs/` folder or a MatrixMarket folder.
- **Web app**: AnnData files written before anndata 0.8 open: anndata 0.7 and 0.6 files and SeuratDisk `Convert()` output (0.7 dataframes with categories under `__categories`, Zarr stores without a root `.zattrs`, 0.6 compound `obs`/`var`/`obsm` with categories in `uns`, `raw.X`/`raw.var`, h5sparse-marked matrices, and every column a SeuratDisk file leaves out of its `column-order`). Each is read as a view of the current layout, so data mapping, annotation and sessions see the same dataset; files written by the real anndata 0.7.8 and 0.6.22 read value-for-value like their current-layout twins.
- **Web app**: A session opens its own dataset. The session manifest records where the dataset came from (`source`: `{kind: 'sample'|'link'|'remote'|'file', …}`, or null for a notebook); **Load State** reads it first and, when another dataset or none is open, opens that one, then restores. Only the fingerprint decides whether the session applies.
- **Web app**: Every analysis can be cancelled. Differential expression, Marker genes (and the heatmap genes it reads afterwards), Correlation, Gene signature, Quick insights, Detailed and the outlier quantiles run as one kind of cancellable job (`app/jobs/cancellable-job.js`): a box shows what the run covers, its progress, bytes read, rate and time left, with **Cancel**. Closing the window or switching dataset cancels the run. A cancelled run keeps nothing, and its reads stop in the reader worker unless another caller still waits for them. Cancelling one analysis worker task no longer ends the whole worker pool.
- **Web app**: The Legend's outlier control is shown for every categorical field and says what it can do: no latent space (with a *Choose latent space* button that opens Data mapping), no category large enough (the fewest cells it needs), not stored (what the first use reads, in MB), reading (progress with *Cancel*) or failed (why; moving the slider tries again). Data mapping marks the stored latent space "outliers stored".
- **Python**: `CellucidSessionBundle` reads session format 2: `dataset_source`, `required_fields`, and `cellucid.session_codecs` (`decode_cell_set` for all five highlight encodings, `digest_cell_set`, `decode_user_defined_codes`, `decode_typed_json`), vectorized with numpy. `apply_cellucid_session_to_anndata` refuses an AnnData missing any obs column or gene the session requires, naming all of them, before any write.
- **Python**: `cellucid serve` gains `--normalize-embeddings` / `--no-normalize-embeddings` and `--centroid-min-points N`, the launch choices `serve_anndata` and `show_anndata` already took.
- **R**: `cellucid_transition_drift(transition_matrix, embedding, normalize_rows)` computes each cell's drift in an embedding from a cell-cell transition matrix (CellRank's, for example), operation for operation as Python's `compute_transition_drift()`. A shared conformance fixture (`transition_drift_conformance.json`) holds both to the same drift: bit for bit where the arithmetic is exact, otherwise within the forward-error bound.
- **Datasets**: The five samples in `cellucid-datasets` are Cellucid-ready AnnData Zarr v2 folders, `anndata/<id>.zarr/`, written by `prepare(format="zarr")` and listed by `anndata/datasets.json` (schema 2, default `suo`, each with a format-2 starting view). They hold the same cells, genes, fields and coordinates as before (every 3D cell-order digest is unchanged), raw counts as float32 instead of 8-bit quantized values, the latent space with its stored outlier quantiles, and Pancreas's velocity under scVelo's names. Sizes: Suo 215.5 MB, Garcia 107.2 MB, He 32.3 MB, Kanemaru 50.2 MB, Pancreas 3.8 MB; no file over 1.6 MB. The app's sample catalog is `https://theislab.github.io/cellucid-datasets/anndata/`.
- **Demo datasets**: `cellucid-demo-custom-datasets` publishes its three datasets as a `catalog/` of Zarr v2 folders with a schema-2 `datasets.json`, and the 2-D dataset as an `.h5ad` and a `.zarr.zip` in `files/`, all written by `prepare()` and rebuilt byte for byte by `generate_datasets.py --check`. The README describes publishing with `.nojekyll` and the link forms.
- **Demo datasets**: `cellucid-demo-custom-datasets` adds a multimodal sample in two forms: `catalog/synthetic-citeseq.zarr/`, a Cellucid-ready MuData Zarr (60 cells, modalities `rna` and `prot`) written by `prepare(data=mdata, format="zarr")`, and `catalog/synthetic-citeseq-seurat.rds`, the same cells as a Seurat v5 object (assays `RNA` and `ADT`, reductions `pca` and `umap`, graph `RNA_snn`, `misc$cellucid`). The generator writes the `.rds` without R, as R serialization format 3 in gzip; R 4.5.3's `serialize()` of the object SeuratObject 5.4.0 builds from the same values is byte for byte the same. Both are listed in the catalog, so its **Sample dataset** list shows an R object.
- **Annotation**: Community annotation keeps every dataset separate. User files are version 2 (the profile once, then one section per dataset id) and every merge in `merges.json` (version 2) names its dataset; schemas `user-v2` and `merges-v2`. Repositories on the version-1 schemas are refused with a message.
- **Docs**: "Host a Dataset on GitHub" (why a Zarr folder is the way to publish on GitHub, `.nojekyll`, measured speed of GitHub Pages and `raw.githubusercontent.com`, GitHub's limits and Git LFS), "Which format?" in `prepare()`, `?cellucid_prepare`, the specification and the preparing guides, and a Memory budgets section in the developer docs.
- **Docs**: Tutorials, each run end to end before publishing: "Tutorials: open published data as it is" (a GEO MatrixMarket sample by link, a Seurat `.rds` on Zenodo by link, a Cell Ranger `outs/` folder from disk); Python "CITE-seq as a MuData: prepare and show" and "Serve Loom, 10x and MatrixMarket files"; R "Explore a Seurat object from R" (`cellucid_view()`, commands, a selection back in R), "Share a SingleCellExperiment on GitHub" and "CITE-seq: a Seurat object with RNA and ADT". The example-link tables of "Open from a link" tag each file Python or R, and list an R object and a GEO MatrixMarket sample.
- **Web app, Python**: Two data-mapping kinds: **Values**, whether each expression matrix (`X`, `raw.X`, each layer) holds `log1p` values and in which base, and **Highly variable genes**, the true/false `var` column that marks them (`raw.var` for `raw.X`). They are chosen under Data mapping › Genes and saved in sessions. The defaults are Scanpy's: with `uns['log1p']`, `X`, every layer and a non-integer `raw.X` are `log1p` in its `base` (none is e), else every matrix is as stored; `var['highly_variable']` when it is true/false, else none. A default read from Scanpy's records is marked guessed. The data mapping (sessions, the Jupyter wire, `cellucid.data_mapping.validate_data_mapping`) gains the keys `log1p` (`{"<matrix>": {"base": null|number}}`) and `highlyVariable`, pinned in the shared conformance fixture. Analyses read them through the gene loader: `log1p({matrix})` and `highlyVariableGenes()`, which reads the column on first use. The sample presets in `cellucid-datasets` carry the two keys.
- **Web app**: A categorical obs column opens in the colours its file stores the way Scanpy keeps them, `uns['<column>_colors']` (one hex colour per category, as `sc.pl.umap(adata, color="leiden")` leaves it), read as Scanpy reads it: a shorter list is not used, a longer one gives its first colours. The colours are read with the column's categories when it is first shown. Without them, or once categories are merged, the colours are Cellucid's palette. A session records a colour only where it differs from the one the field opens with.

- **Web app**: Every opening shows its progress from the first moment: a link (`?data=` or the address bar, every provider), a file, a `.zarr.zip` or a folder, a sample, a catalog dataset, a `cellucid serve` server and a notebook. While no dataset is shown, a panel over the empty viewer names what is opened and the step under way (*Finding the file on NCBI GEO…*, *Reading the file's structure…*, *Reading embedding X_umap…*), with a progress bar, how much has arrived, how fast and how long is left when that is known, and **Cancel**. While another dataset stays on screen, the same appears as a notification with Cancel. A request made again after a failure shows as *Retrying…* with why. Cancel stops the opening where it is (its reader worker and every request end) and keeps the dataset on screen; how an opening ended is kept in the notification center. Every layer that opens data reports into one opening (`data/opening-progress.js`, shown by `app/opening-progress-view.js`), so a link, the file it names and the dataset change that draws it show as one.
- **Web app**: **Loom, 10x Cell Ranger and MuData files open directly**, read in place like an `.h5ad` (only what a view needs), through layout views that present each file as AnnData (`data/anndata/views/`): no conversion step. **Loom** (`.loom`, Loom 2.0.1 and 3.0.0, loompy 2 and 3, velocyto, SCope/pySCENIC, loomR's orientation): `/matrix` and the layers as `X` and `layers` (a transposed view), cell attributes as `obs`, two-dimensional ones (`X_umap`) and coordinate columns (`_X`/`_Y`, `_tSNE1`/`_tSNE2`, `umap_1`/`umap_2`) and SCope's compound embeddings as `obsm`, `col_graphs` as `obsp`, Loom 2's XML character references decoded. **10x Cell Ranger HDF5** (`.h5`, Cell Ranger 1.2 to 9, one group per genome or `matrix/features`): `X` is the stored matrix as it is (CSC over barcodes is CSR over cells), var columns named as `scanpy.read_10x_h5` names them (`gene_ids`, `feature_types`, `genome`). A matrix has no embedding, so a lone `.h5` is refused saying to open its **Cell Ranger `outs/` folder** (Open… › Folder…, or a dropped folder), which opens with Cell Ranger's own `analysis/` projections (`X_umap`, `X_tsne`, `X_pca`) and clusterings (`graphclust`, `kmeans_<k>_clusters`) matched by barcode; `molecule_info.h5` is refused naming the matrix to open. **MuData** (`.h5mu` with its 512-byte user block, and MuData Zarr v2 and v3; mudata 0.2, 0.3 and 0.4): the global `obs`, `obsm` and `obsp`, each modality's obs columns and embeddings as `<modality>:<key>` (mudata 0.4 no longer copies them to the global `obs`), a MuData sharing features instead of cells refused. A Loom file with no embedding is refused with the two lines that compute one in Python and in R. The file picker, drops, links and the format registry list `.loom`, `.h5` and `.h5mu`. A corpus written by loompy 2.0.17 and 3.0.8, as Cell Ranger writes (with scanpy's Cell Ranger 1.2 test files), and by mudata 0.2.4, 0.3.2 and 0.4.1 (`scripts/formats-corpus/`, `tests/fixtures/{loom,tenx,mudata}/`) is read value for value against what loompy, scanpy and mudata read (`tests/formats-corpus-conformance.test.mjs`).
- **Web app, Python, R**: **Modalities**: a dataset of several feature spaces over the same cells (a MuData's modalities, the feature types of a 10x matrix such as Gene Expression and Antibody Capture, the genomes of a Cell Ranger 2 file) opens with its gene-expression modality, and Data mapping › Genes shows one **Modality** select above Expression for it (nothing changes for a dataset of one feature space). Gene names, expression, its `log1p` and the highly variable genes are then the modality's; gene search, colour by gene, Marker genes, Differential expression and Dot plot read it. A modality that measured only some cells reads NaN for the others (not 0), through MuData's `obsmap`. The 10x feature types are column windows of the one matrix and share its gene index. The data mapping (sessions, the Jupyter wire, `cellucid.data_mapping.validate_data_mapping`) gains the key `modality`, present only for a multimodal dataset, so every existing session, preset and fixture is unchanged; the launch defaults gain `modality`. The shared conformance fixture holds the rules (`defaultModality`, the modality a dataset opens with: the first named `rna`, `gex`, `gene expression`, `gene_expression` or `sct`, else the first; `mappingModality`, a mapping's modality checked against the dataset's), pinned in all three repositories; Python adds `default_modality()` and `require_mapping_modality()`. Every multimodal store is presented in MuData's layout (`/mod/<name>`, `mod-order`, `obsmap`), which is the one contract a writer of several modalities follows.
- **Web app**: **Seurat and SingleCellExperiment objects open directly** from `.rds`, `.rda` and `.RData` files: no conversion step and no R. The web app reads R's serialization itself (`assets/js/formats/rds/`: versions 2 and 3; gzip natively, bzip2 and xz with its own decoders, no WebAssembly; every item type R writes, including closures, byte code, environments, references, S4 objects and the ALTREP compact sequences, deferred strings and wrappers), in one pass over the file that keeps only what is shown (`scale.data`, `commands`, `tools`, neighbours, loadings and models are passed over; a Seurat v5 layer, whose name comes after its values, is passed over while the file streams and read in a second pass only if it turns out to be needed), into a store kept in the browser's origin-private file system for the tab's life (memory where there is none), opened as AnnData (`data/anndata/transcode/`). An uncompressed `.rds` (`compress = FALSE`, `readr::write_rds()`) is not copied: its matrices are read in place, from disk or a link. **Seurat** (SeuratObject 3.0 and later, `Assay`, `Assay5` with split layers joined, `SCTAssay`, Signac's `ChromatinAssay`): `meta.data` and `active.ident` as obs, every reduction as `X_<name>`, `graphs` as `obsp`, the default assay's `data` (else `counts`) as `X` with the other as a layer, variable features as `highly_variable`, other assays as modalities (the default assay first), Visium and FOV positions as `spatial`, `misc$cellucid` as `uns['cellucid']`. **SingleCellExperiment** (1.8 and later) and **SummarizedExperiment**: `colData`, `reducedDims` (`X_<name>`, case kept), `logcounts` (else `counts`) as `X`, `rowData`, `altExps` as modalities, `colPairs` as `obsp`, `metadata$cellucid`. A `save()` file of several objects keeps all of them in the one pass and the opening asks which to open (one button per object); choosing another later in the tab reads nothing again. Refused with the R line that fixes it: matrices kept outside the file (BPCells, HDF5Array/DelayedArray, `SaveSeuratRds(move = TRUE)`), Seurat v2 objects, h5Seurat files, objects with no embedding, zstd, ASCII and native-binary saves. On a 100,000-cell Seurat object (2.45 GB gzip, 5.3 GB decompressed, counts and data with 1,500 genes per cell, `scale.data`, PCA, UMAP, an SNN graph) the pass takes 1.3 times plain decompression (16.6 s in Node on an M1 Pro), keeps 3.5 GB on disk and peaks at 277 MB of memory; reopening it in the tab takes milliseconds. The opening shows the pass as one step with bytes, rate, time left and Cancel; a cancelled or failed pass leaves nothing on disk, and a closed tab's stores are removed at the next start.
- **Web app**: **MatrixMarket folders open** with **Open… › Folder…** or a drop: a count matrix (`.mtx`, `.mtx.gz`) with its feature and barcode lists as Cell Ranger 2 and later, STARsolo, kb-python, alevin-fry and GEO uploads lay them out, its orientation read from its size, written once into a cells × features matrix as the entries stream. A matrix that arrives in an order Cellucid cannot keep (a MatrixMarket file in no order, Seurat v5 split layers, an R `dgTMatrix`) is put in order in blocks of 4 Mi entries in browser storage, never whole in memory; one written by feature (kb-python, alevin-fry) is kept as it arrived, as a `csc_matrix`. Several feature types (`features.tsv` column 3: Gene Expression with Peaks or Antibody Capture) open as one modality per type, Gene Expression first, as a 10x `.h5` does. The folder's own cell tables give the embedding: Cell Ranger's `analysis/` projections and clusterings, and any `.csv`/`.tsv` table of the cells (`UMAP_1`/`UMAP_2`-style columns become embeddings, the rest obs fields; a clustering column such as `Cluster`, `seurat_clusters`, `leiden`, `louvain` or `kmeans_10` is categorical even when its labels are numbers); without one the folder is refused with how to compute one. A Cell Ranger `outs/` folder that also holds the `.h5` opens through the `.h5`, in place.
- **Web app**: The Dataset box says what became of an R object or a MatrixMarket folder read once into the browser: one line under its format, *Read once into this browser's storage for this tab*, *Opened from this browser's storage, read earlier in this tab* or *Read once into memory for this tab* (no origin-private file system, as in a Firefox private window), with *· N parts not read* when something was left out, and in Details the rows **Kept** (where, until the tab closes) and **Not read** (one line per part, named by its R path and its size in the file: `assays$RNA@scale.data (15 KB): Cellucid does not use it.`). The R parser records the bytes each dropped part spans and whether it held anything, so the empty `scale.data`, `commands` and `neighbors` every Seurat object carries are not listed; a Seurat v5 `scale.data` layer is listed too. The **Read** row of such a dataset says *Once, from start to end*, and its origin line names the format alone (*R object (.rds) · this computer*).
- **Web app**: The R and MatrixMarket corpora (`tests/fixtures/rds/corpus.zip`): objects written by SeuratObject 4.1.4 and 5.4.0 and SingleCellExperiment 1.32 in every compression and both serialization versions, SeuratObject's own `pbmc_small`, refusal cases, and MatrixMarket folders in each producer's layout (`scripts/rds-corpus/`, `scripts/mtx-corpus/`), read value for value against what R's accessors and scanpy read (`tests/rds-corpus-conformance.test.mjs`); every transcoded store also reads in anndata and mudata (`scripts/rds-corpus/verify_with_python.py`). The parser is checked item type by item type, at every split of its input (`tests/rds-parser-contract.test.mjs`), and the kept stores' lifetime with an origin-private file system on disk (`tests/transcode-space-contract.test.mjs`).
- **Web app**: **Figshare, Google Drive and OneDrive** open from a link (the address bar or `?data=`), read in place with range requests ([Google Drive, OneDrive and private sharing](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/16_private_drive_onedrive.html)). **Figshare**: an item page (figshare.com or an institution's portal), a Figshare DOI or any DOI that resolves to one, or `figshare:<id>[v<n>][/<file>]`; the page asks `ndownloader.figshare.com` for the file's signed S3 link (its bot check refuses requests from a reader worker in some browsers) and the reader reads S3 with it for the link's stated ten seconds, then asks the page for the next with a one-byte request; private links (`figshare.com/s/…`) and bare download links are refused with the reason. **Google Drive**: a file shared with anyone opens with no sign-in through the Drive API and Cellucid's public, referrer-restricted API key, its details and first bytes asked for together (one round trip); a private file shows "Private file on Google Drive. [Sign in to Google]", signs in with the `drive.file` permission only and opens Google's file picker already showing that file, so one click on **Select** grants it. **OneDrive and SharePoint** (personal, work and school accounts): every file opens after **Sign in with Microsoft** (OneDrive serves nothing anonymously; PKCE in the page, no client secret, no library), with `Files.Read`, and `Files.Read.All` asked for only for a colleague's file (which many organisations let only an administrator approve); a share link is redeemed durably, so sessions and `?data=` name the file as `onedrive:<drive>/<item>`; the short-lived download address is renewed through Microsoft Graph once for all reads waiting on it, and a file changed meanwhile is refused. Folder links are refused (share one `.h5ad` or `.zarr.zip`). Sign-ins live in the tab's `sessionStorage` through one store shared with the GitHub sign-in of community annotations, never in a link, a session, the console, analytics or an error; an expired one shows a calm **Sign in again** under the address bar and as a notification, and the reads waiting on it resume. A read behind a sign-in or an expiring link asks the page for its URL and headers once for all parallel reads (`httpByteSource`'s new `access`), and any redirect that states its lifetime (`X-Amz-Expires`) is reused for that lifetime and no longer. Google answers only the origins Cellucid registers (www.cellucid.com and a local `cellucid serve` on port 8765, `provider-config.js` `registeredOrigins`), for public files and sign-ins alike; on any other page Cellucid says where Drive works before asking Google (a local server on another port: start it without `--port`; a notebook viewer on another port: stop the server holding 8765 and show the data again; another host: www.cellucid.com). Microsoft sign-in is refused up front, with the reason, on a plain `http://` page of another computer and inside another site's page (a notebook at `localhost` showing a viewer at `127.0.0.1`, whose storage the browser keeps apart from the sign-in window); a sign-in window found closed offers the step again, and a sign-in finished in it after all still counts until the step ends. The site's privacy policy is `www.cellucid.com/privacy.html`, linked as **Privacy** in the sidebar footer; the public client IDs, API key and project number are in `assets/js/data/sources/provider-config.js` (registered for www.cellucid.com and `http://127.0.0.1:8765` / `http://localhost:8765`), and a copy of the app without them says "not configured on this site yet". Checked against the real services: public Drive files open with no sign-in and are read in place, and Microsoft's sign-in page accepts all three return addresses. Google holds back downloads without a sign-in from a network that made many of them (an HTML "unusual traffic" page, HTTP 403, without CORS headers, while the file's details still answer and signed-in requests go through; seen three times on 2026-10-25, each after some 200–400 MB from one address within minutes, whether in 30 requests or 400, for about 50 minutes): its head read is asked once, and a failure confirmed with one byte (Google sends its page after about seven seconds, so this is told in about eight seconds); a file whose first bytes are held back opens after "Google Drive sends this file only with a Google sign-in for now. [Sign in to Google]" and the picker, and a read held back later says so, in place of the generic "no connection, or no CORS headers" (`registerSourceAccess`'s `unreadable`). A Drive file granted to Cellucid is read with the sign-in while one is held, public or not. `drive.google.com/u/<n>/uc?id=…` links (as gdown writes them) are Drive links, and the address of a OneDrive or SharePoint page that is not a share link is refused with how to copy the file's link. The Dataset box's read line ("5.5 MB read in 43 requests") is shown for a file read at an address with a query string (Google Drive's `…/files/<id>?alt=media`, a Dropbox `?dl=1` link): the box names a file by its address without the query, and the Transfers log now finds its entry by that address. The documentation lists public Google Drive files (Python and R) with what each first view read.
- **Web app**: Google Drive and OneDrive in every browser. Time spent in a sign-in window, a consent screen, Google's file picker or the choice of an R file's object never counts as a stalled download: every wait on the person registers itself (`utils/person-wait.js`), the startup deadline stands still meanwhile and starts over once the person has answered, and the opening goes on by itself after the sign-in; a data link's stall says "while opening the data link". Signed in to Google, a Drive link to a file not yet chosen for Cellucid (the `drive.file` sign-in reads only those) opens Google's picker on it at once and is read with the sign-in; closing the picker reads a public file with the key unless Google holds it back, and a read held back later says to open the file again with the Google account. A Drive download held back after the file opened is told in about 9 s, as at opening: the first read that gets no answer asks Drive, from the file's access, for the file's details and one byte, and the details answering while the byte fails ends the read in Drive's words (`unreadable.confirm()` in `sources/access.js`, asked by the reader worker with `source-unreadable`); a check that does not confirm the hold leaves the read's failure as it was. Signed Drive reads go to one stable address per request in flight (`lane=<n>`, a parameter Drive ignores) under the HTTP cache mode, so Chromium and Firefox preflight each address once: on a Drive-shaped host (0.65 s a request, 0.2 s a preflight) a signed first gene takes 5.3 s, as with no sign-in. The head a provider reads beside a file's details answers every read that starts inside it, so an R object from Drive or OneDrive is not read twice (the test kit's Seurat `.rds`: 9.4 MB in 3 requests). A OneDrive item is asked for whole, since Microsoft Graph leaves the download address out of an answer whose fields are selected (`$select`), even when the selection names it; an answer without what reading needs is said in Graph's terms ("Microsoft Graph described … without a download address"). **Sign in with Microsoft** opens its window empty, on Cellucid's own address, in the click, and sends it to Microsoft once the attempt's PKCE verifier is made, so Safari on a Mac with Microsoft's single sign-on extension (Company Portal) shows Microsoft's ordinary sign-in page for any account instead of handing the window to the extension; a second click brings back a sign-in window still open, and a pop-up blocker's closed window is said as a blocked window. Sign-in failures are said for what they are, with the next step: an organisation that lets only an administrator approve Cellucid (`AADSTS90094`, `90095`, `65001`, naming the organisation), an account it did not assign (`AADSTS50105`), its sign-in policy (`AADSTS53003`), a decline (`AADSTS65004`, which is also how Microsoft reports leaving its "Need admin approval" page, so both are named), a site registration Microsoft or Google does not accept, a Google Workspace administrator's block (`admin_policy_enforced`), and Google's window closed (Google shows "Access blocked" there). Where the same account would fail again the button offers **Use another account**, and the line links the administrator steps in the documentation.
- **Web app**: **Composition** (Analysis → Cells) asks whether cell-type proportions differ between conditions, accounting for samples. It counts the cells of each type per sample in one pass in a worker (a sample spanning conditions counts once per condition) and, with a condition, tests every type across the samples: logit proportions `(count + 0.5) / (total + 1)` with Welch's t-test or one-way ANOVA (after propeller, without its empirical-Bayes moderation), ranks of proportions (Wilcoxon or Kruskal-Wallis) or centred log-ratios, with Benjamini-Hochberg or Bonferroni over the types. Without samples it runs χ² tests (Fisher's exact test for a sparse 2 × 2 table) and says the cells were the replicates. It shows stacked bars per sample or condition, per-sample points per type, a table, the results as CSV and the counts as CSV; the statistics match SciPy on the shared fixture.
- **Web app**: **Dot plot** (Analysis → Genes) draws, for each gene of one or more named gene sets and each group of a categorical field, the mean expression as color and the fraction of cells expressing as dot size, as scanpy's `pl.dotplot` (or as a matrix plot, or as a stacked violin of each gene's distribution per group): scaled 0–1 or to z per gene, a free expression cutoff, the mean of expressing cells only, swapped axes, groups ordered by size or a dendrogram, Reds, Viridis or Blues. Only the listed genes are read, each summed per group in one pass over its non-zero entries in a worker; display changes never read again. Genes the dataset lacks are left out on a chip. It exports the long table (`group, gene, gene_group, mean, fraction_expressing, n_cells`) and opens Gene signature with its sets (**Score as signature**). Means and fractions match scanpy's `dot_color_df` and `dot_size_df` on the shared fixture.
- **Web app**: Marker genes, Differential expression and Dot plot have an **Expression matrix** option: the matrix the data mapping shows (default), or `X`, a layer or `raw` of the same file for that run alone, without changing what is shown. Each matrix's genes are named on its own axis (`raw` has `raw.var`), read through the same `var/{axis}/{names}` routes as a mapped matrix on every way of opening a dataset (a file, a link, `cellucid serve`, a notebook); a gene the matrix lacks is named as not found; Automatic fold change follows that matrix's own `log1p`. The choice is part of the marker cache key, the result line names it, and sessions save it with the options and the result.
- **Web app**: Differential expression has a **Replicates** option: with a sample field the run is a pseudobulk test, each sample's cells of a group one observation (its mean expression, from one pass over each gene's non-zero entries) and Welch's t between the two groups' samples, as `scipy.stats.ttest_ind(equal_var=False)` on the per-sample means; no variance moderation. Each group needs cells in at least 2 samples.
- **Web app**: Clicking a gene in a Marker genes or Differential expression result colours the cells by it through the sidebar's gene selector (**Coloring & filtering → Gene expression**), whose search box then names it.
- **Web app**: Analysis forms share one layout: the essentials, then one collapsed **More options** holding a collapsed group per topic, each counting the options changed from their defaults, with **Reset to defaults**; then **Run**, with the reason above it when it cannot run. Every option has an (i) tip and a default stated in the documentation; the options and which groups were open are saved in sessions. After a run, a change that would change the numbers dims the result with **Stale · Run**; display options (filters, ranking, labels, plot, rescaling) redraw it without computing again. Correlation and Gene signature run on **Run** instead of on every change.
- **Web app**: The site has **terms of use** (`www.cellucid.com/terms.html`: BSD 3-Clause license, no warranty, your data stays yours, other services keep their own terms) and an **accessibility statement** (`accessibility.html`: what works, the WebGL plot's limits, how to report a barrier), reading pages in the privacy policy's design, linked as **Terms** beside **Privacy** and **Accessibility** beside the licence in the sidebar footer; the privacy policy now names Google Analytics' two cookies. `/.well-known/security.txt` (RFC 9116) names the security policy's contact, its `Expires` renewed a year after each build by `scripts/generate-web-asset-inventory.mjs`.

- **Python, web app**: The viewer channel. A notebook viewer opens one WebSocket to its server, `<server>/_cellucid/ws`, named by `GET /_cellucid/protocol` (`{"channel": "_cellucid/ws", "commands": [...], "binaryCommands": [...], "events": [...]}`), and says `hello` with its viewer id and token; the server checks the `Host` and `Origin` headers (the server's own host, a loopback name or an `allowed_hosts` name), compares the token in constant time and answers `welcome`. Commands are JSON text frames, and a command that names cells (`highlight`, `setVisibility`) one binary frame: a `uint32` header length, the JSON header padded to a multiple of 4 bytes, then the cells as little-endian `uint32`, strictly ascending, which the viewer reads in place and checks as strictly as the hosts check its binary events (a host whose protocol document lists no `binaryCommands` is told to upgrade). Commands are queued until the viewer connects and sent in order, never twice; events are JSON text frames, and a selection (uint32 cell indices) and a session bundle (streamed to a temporary file) binary frames. Close codes: 1001 stopped, 4401 unauthorized, 4409 replaced by a newer page of the same viewer, 1008 protocol violation. The server pings every 20 s; the viewer reopens a dropped socket once, at once, and one health request tells a stopped host from a proxy that does not pass WebSockets. The rules are one document, `viewer-protocol.json` (web app `data/viewer-protocol.js`), kept byte for byte in cellucid-python and cellucid-r. Two million selected cells reach a hook in 89 ms.
- **Python**: `cellucid serve`, `serve_anndata()`, `show_anndata()` and `AnnDataViewer` take every format of the input format registry by path: `.h5ad`, AnnData Zarr (folder or `.zip`), `.h5mu`, MuData Zarr, `.loom`, 10x `.h5` files and Cell Ranger `outs/` folders, MatrixMarket folders (with GEO's `GSM…_` prefixes) and `.rds`/`.rda`/`.RData` (gzipped or not); `anndata.file.kind` and the health `format` carry the registry id. A folder of files is served file by file and listed in `/_cellucid/datasets` (`file.files: [{path, size}]`); the health and info `n_cells`/`n_genes` are `null` for a format Python does not read. An in-memory `MuData` is written once to a temporary Cellucid-ready `.h5mu`.
- **Python**: `prepare(data=...)` takes an AnnData or a MuData whole. Every component it holds (obs, X with var, obsm embeddings and vector fields, `obsp['connectivities']`, `uns['cellucid']`) is the default, an explicit argument replaces one component, and the output is byte-identical to passing the components by hand. `embeddings`, `dataset_id` and `dataset_name` are required only without `data`. The caller's arrays are never changed: a matrix that needs canonicalising (sorted indices, summed duplicates, dropped stored zeros) is copied first, and value ranges are checked without whole-matrix temporaries, so peak memory beyond the input is about one converted matrix. A backed AnnData (`backed="r"`) is prepared in two passes over row blocks of its file, never held whole.
- **Python**: `prepare()` writes a Cellucid-ready MuData from a MuData: `format="h5mu"`, or MuData Zarr with `"zarr"` / `"zarr.zip"`. Each modality is a Cellucid-ready AnnData, the shared cells, embeddings and graph are at the root, and modalities are written one at a time to bound memory (30,000 cells, 20,000 genes and 200 proteins: 4.8 s and 316 MB peak to a Zarr folder). New optional extra: `pip install "cellucid[mudata]"`.
- **Python**: New launch default `modality` (`prepare(modality=...)`, `uns['cellucid']['defaults']['modality']`, `serve_anndata`/`show_anndata(modality=...)`, `cellucid serve --modality`): the modality a multimodal dataset opens on.
- **Python**: `apply_cellucid_session_to_anndata` applies a session to a MuData; the added columns go to the global `mdata.obs`, and each modality's arrays are placed among the cells only when read.
- **Web app**: GEO links open MatrixMarket samples: a series' `GSM…_matrix.mtx.gz` with its `_features`/`_genes` and `_barcodes` files is one dataset, opened by its matrix (`geo:GSE…/GSM…_matrix.mtx.gz`), read from the sample's own folder with the sizes `filelist.txt` states. `GSM…` sample accessions are accepted. An incomplete or ambiguous sample is refused naming the file it lacks or the two that compete.
- **Web app**: `.rds.gz`, `.rda.gz` and `.RData.gz` (GEO's gzipped R objects) open from links, providers, drop and the File… picker; the content decides whether a file is gzipped once or twice, and it is decompressed as it streams. A host that sends an R file as `Content-Encoding: gzip` even for a part of it (NCBI's GEO file host does, for `saveRDS()`'s default gzip) cannot be read in parts by any browser: its file is downloaded whole once, decoded by the browser, and the notification says why (GSE327442's 25 MB Seurat object opens in 7–10 s in Chromium, Firefox and WebKit).
- **Web app**: A folder of MatrixMarket files sharing a sample prefix (an extracted GEO `_RAW.tar`) opens with Open… › Folder…; with several samples, Cellucid asks which to open.

### Changed

- **Python, R, web app**: A host program opens its viewer at `?viewer=true&viewerId=…&viewerToken=…` (it was `?jupyter=true`), and `/_cellucid/info` names the program serving it, `runtime`: `"python"` (`cellucid serve`, a viewer outside a notebook), `"jupyter"` (a viewer in a notebook kernel) or `"r"`. The viewer shows it in the Dataset box ("h5ad file · Jupyter", "Seurat object (.rds) · R") and in the dataset list ("R data"); its startup steps name "the viewer's host". The page keys and the runtime names are part of the shared viewer protocol file (`page`, `runtimes`; SHA `642c04e1…` in all three repositories). The old key opens nothing.
- **Python**: A notebook viewer (`show_anndata`, `AnnDataViewer` with `port=None`) serves on 8765, the one local port Google sign-in accepts, or on the next free port when another server holds it, saying so; it used a port the operating system picked, where Google sign-in never worked. cellucid-r's `cellucid_view()` already did this.
- **Python**: `prepare(format="h5ad")` (and `.h5mu`) compresses chunks on every core, written with HDF5's direct chunk write in grid order, so the file is laid out as HDF5's own writer lays it out: 1M cells (150 values each, 1.15 GB of X) in 13.3 s at 4.5 GB peak (63.9 s when HDF5 compresses one chunk after the other on one thread), 2M cells in 22.1 s, 200,000 cells in 2.8 s.
- **Python**: Applying a session copies a continuous field's missing values (NaN) as the viewer copied them, and accepts field names and highlight fields with `:` (ATAC peaks, MuData `<modality>:<column>` columns), as the viewer does; both were refused.
- **Python, development**: The `dev` extra brings the `mudata` extra, so a development install, CI and the release gate run the MuData tests instead of skipping them. On Python 3.11, where pip installs mudata 0.3, the suite (`tests/conftest.py`) sets mudata's `pull_on_update=False`, so every MuData it builds or reads behaves as under mudata 0.4. `tests/test_minimum_python_grammar.py` parses every Python file with the grammar of the oldest supported Python and refuses the f-strings only Python 3.12 accepts.

- **Web app, Python, R**: One input format registry decides what a file or folder is, everywhere (`assets/js/data/input-formats.js`). Each format has a stable id (`h5ad`, `zarr`), the name endings of each layout it arrives in (`file`: `.h5ad`; `zip`: `.zarr.zip`; `folder`: `.zarr`), container signatures and folder markers. The File… picker's `accept` list, a drop, the address bar's routing, the link providers (GitHub, Zenodo, GEO, Hugging Face, any web address), the dataset catalog, the served-store validator, the server health checks and the reader worker all ask it; the rules each of them spelled out before are gone. A name decides whether a file is taken; its content decides what it is, so a misnamed file (an `.h5ad` saved as `.zip`, a Zarr ZIP saved as `.h5ad`) opens as what it is, and an HDF5 file with a user block of up to 2 KiB is recognised by its signature. `tests/fixtures/input-formats.json` is the registry's document; cellucid-python ships the same bytes (`src/cellucid/input_formats.json`), from which `cellucid serve` classifies a path (`classify_anndata_path`) and names its `kind`, and cellucid-r ships them as `inst/input_formats.json`, which `cellucid_view()` and `cellucid_serve()` read; the three repositories pin one SHA-256. Every message that names what Cellucid opens (a link that is not a dataset, a catalog path, a Drive, OneDrive or Dropbox folder, a listing with nothing to open) lists the registry's single-file formats (`openableNamesText`).
- **Web app**: One data source reads every dataset file, `AnnDataFileSource` (`data/anndata-file-source.js`), addressed by one scheme, `anndata://<dataset id>/`; it replaces `H5adDataSource` and `ZarrDataSource` and their `h5ad://` and `zarr://` schemes; the data layer addresses an open file by its AnnData elements (`obs`, `var`, `obsp`, `obsm/{i}/{cols}`, `outliers/{l}/{f}`, …) instead of export file names. Its loaders are `loadFile` (a local or remote file), `loadFolder` and `loadRemoteFolder`; the format is a fact of the open dataset (`format`, `describeStorage()`: `{format, layout, zarrVersion, remote, url, name, bytes}`). A dataset that states no id is still `<format>_<name>` (`h5ad_<file>`, `zarr_<store>`), so every saved session keeps matching. A local file's source line names its format (`AnnData (.h5ad)`, `AnnData Zarr`).
- **Web app**: Opening has one entry, `openDatasetStore({source, expect})` in `data/anndata/open.js`, over one file `{file}`, a picked folder `{files, rootName}` or a web folder `{url, name}`; it opens the container and presents the root through a layout view when one recognises it (`presentStore`, `data/anndata/views/`), or transcodes a single-stream format into an AnnData Zarr store (`data/anndata/transcode/`). The view nodes `older-layouts.js` is built from are shared in `data/anndata/views/view-nodes.js`, with transposed, row-window, column-window and computed arrays for the layout views.
- **Web app**: The usage analytics' local load methods are named by how a dataset arrived, whatever format it holds: `local-user-file`, `local-user-zip` and `local-user-folder` (they were `local-user-h5ad`, `local-user-zarr-zip` and `local-user-zarr-folder`); `data_load` reports the input format it was read as (`data_format`: `h5ad`, `zarr`, `loom`, …) on success. The Dataset box names the format (`Loom (.loom)`, `MuData Zarr v3 folder`, `Cell Ranger outs/ folder`) and, for a multimodal dataset, its modality.
- **Python, breaking**: `prepare(out_path, *, format, embeddings, dataset_id, dataset_name, …)` writes a Cellucid-ready AnnData and returns its path. `format` is required, with no default: `"h5ad"` (one file, for single-file hosts such as Zenodo, Figshare, Hugging Face, S3 or a lab server, local files, `cellucid serve` and Jupyter), `"zarr"` (a `.zarr` folder of small files, the way to publish on GitHub: no file over 100 MB, no LFS) or `"zarr.zip"` (the Zarr store in one file). `out_path` must end in the format's suffix. `zarr_format` (2 or 3; 2 when None) applies to the Zarr formats only. New arguments: `citation_text`, `citation_url`, `latent_key` (an `embeddings` key), `normalize_embeddings`. Every obs and var column is stored; `obs_keys` sets which columns the viewer lists.
- **R, breaking**: `cellucid_prepare(out_path, format, embeddings, dataset_id, dataset_name, …)` writes the same Cellucid-ready AnnData as a Zarr store: `format` is `"zarr"` or `"zarr.zip"`, required, matched as `match.arg()` matches; `zarr_format = 3` writes Zarr v3. There is no `"h5ad"`: the package writes no HDF5 and depends on no HDF5 library.
- **Python, R**: The layout both writers share: `X` gene-major (CSC from sparse input; dense stays dense) in float32, or float64 when the genes kept hold integer counts above 2^24; `obs`/`var` with anndata's own encodings and pandas' category code widths (int8 below 127 categories, int16 below 32,767); `X/indices` as `<i4` unless the axis it indexes has 2^31 entries, `indptr` as `<i8` from 2^31 stored values; values that float32 rounds to zero dropped rather than stored as zeros; every embedding whole in `obsm` (float32, unnormalized); `obsp['connectivities']` as symmetric float32 CSR; chunks of 2^16 values for `X` and 2^18 for arrays read whole; deflate level 6, with byte shuffle on every multi-byte array except the expression values: in Zarr v2 numcodecs' `shuffle` filter and zlib, in Zarr v3 only codecs of the v3 specification, the `blosc` codec (byte shuffle, zlib level 6) on shuffled arrays and gzip on the others (`numcodecs.shuffle` is outside the specification; on He the v3 store is 29.00 MB, against 37.62 MB unshuffled); consolidated metadata. A `.zarr.zip` stores chunks first in grid order and the consolidated metadata last, with fixed timestamps and permissions, so it opens with the read of its end and the same input gives the same bytes. Python's `.h5ad` declares every element before writing values, so a reader gets the whole tree and the values opening needs from its first 512 KiB (431 KB on Suo); the file reads the same in anndata 0.8 to 0.12, scanpy and h5py. R and Python write identical Zarr v2 and v3 stores (every metadata document and every value), checked by cellucid-r's `tests/testthat/test-cross-writer.R`; R writes every NaN and NA as numpy's canonical quiet NaN.
- **Web app, Python, R, breaking**: A vector field binds by scVelo's naming on any embedding: `obsm['<field>_<basis>']` is the field `<field>` over `obsm['X_<basis>']` (`velocity_umap` over `X_umap`, `velocity_branches` over `X_branches`), the longest existing basis winning, and its id is `<field>` (`velocity`), so one field spans every basis. `prepare(vector_fields=...)` and `cellucid_prepare(vector_fields = ...)` write each array under that key, with as many columns as the embedding; arrays of another width, colliding keys and keys the viewer would read otherwise are refused before anything is written. `vector_field_default` names the field and must be drawn on an embedding the viewer shows at open. The rule is pinned in the shared data-mapping conformance fixture (web `parseVectorFieldKey`, Python `cellucid.data_mapping.parse_vector_field_key`, R `.parse_vector_field_key`).
- **Python, breaking**: `add_transition_drift_to_obsm(adata, T, basis=..., field_prefix=...)` computes the drift over `obsm['X_<basis>']` (every column) and stores `obsm['<field_prefix>_<basis>']`; `dim` is removed (use `basis="umap_3d"`). `compute_transition_drift` takes an embedding of any width, and a zero row sum is refused.
- **Python**: `cellucid serve <file>`, `serve_anndata()` and `show_anndata()` hand the AnnData store to the viewer, which reads it in place with range requests, as it reads a file or a link. The server serves the store's bytes at `/_cellucid/anndata/<dataset id>.h5ad` (or `.zarr.zip`, or `.zarr/<key>`) with exact `Range`, `HEAD` and strong `ETag` answers, over HTTP/1.1 with one thread per connection and `sendfile`, and names the dataset, its store and the launch choices in `/_cellucid/datasets`. Name and id come from the store's `uns['cellucid']` (`--dataset-name`/`--dataset-id` are needed only when it has none), and every launch choice is checked against the store before the server binds. An in-memory AnnData, or a store with a codec the browser cannot decode, is written once to a temporary `.h5ad`, one matrix at a time, and removed when the server stops. On a 562k-cell × 8,192-gene `.h5ad`: startup 0.16 s instead of 3.0 s, a UMAP about 40 ms instead of 170–250 ms, 100 genes 0.1–0.6 s instead of 115 s.
- **Python, breaking**: `AnnDataViewer` is the one notebook viewer class (`cellucid.jupyter.BaseViewer` is removed; same methods and signature).
- **Web app, breaking**: `?github=owner/repo/path` is replaced by `?data=gh:owner/repo/path`, and the **GitHub** control by the **Link or server** address bar; the source type `github-repo` is `link`. `?exportsBaseUrl=` is `?sampleCatalog=`, and the `cellucid-exports-base-url` meta element is `cellucid-sample-catalog`.
- **Web app, breaking**: `datasets.json` schema 2: `{"schema": 2, "default": "<id>", "datasets": [{"id", "path", "session"?}]}`, `path` a Zarr folder, an `.h5ad` or a `.zarr.zip` below the catalog folder, `session` (`{path, bytes, sha256}`) the dataset's starting view. The sample list and every folder link read this one format through one loader; listing reads each file's metadata only, and a catalog link without `dataset=` opens its default. A catalog in another schema is refused.
- **Web app, breaking**: Session format 2. Files start with `CELLUCID_SESSION/2\n`; any other file is refused before anything is read ("This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2)."), in the web app, the Jupyter upload and `CellucidSessionBundle` alike. Sessions now save Figure Export settings, the docked analysis panel with analysis results (differential expression, marker genes, gene signature and correlation, as typed columns; floating windows carry theirs), per-view centroid toggles and the app theme, and no longer save derived labels, duplicated cameras, antialiasing or per-cell gene caches. Highlights are stored in the smallest of sorted gaps, a bitset, runs or relative to their category (checked by count and digest on restore); user-defined codes raw, as runs or bit-packed. The fingerprint records the dataset's identity id instead of how it was opened, and the manifest lists every obs field and gene the session uses (`requires`), so a dataset lacking any is refused naming every one. Restores request every field together, rebuild each kept view once and finish in background tabs and off-screen notebook frames. Suo (562k cells), a heavy session of 7 views: file 317 → 193 KB, save 189 → 101 ms, restore 1.89 → about 0.83 s, with no task over 50 ms.
- **Web app, breaking**: A data mapping is always `kind: "anndata"`; the `export` kind is deleted and a mapping or session naming it is refused.
- **Web app**: Own HDF5 and Zarr readers (`assets/js/formats/`) under one AnnData reader (`assets/js/data/anndata/`) replace the vendored h5wasm and the two separate loaders. They read Zarr v2 and v3 (sharded included), from a folder, a `.zarr.zip` or a web host; Zarr codecs gzip, zlib, Blosc, Zstandard, LZ4, shuffle and CRC-32C; HDF5 deflate, shuffle, Fletcher-32 and LZF. A codec or HDF5 filter outside these (Blosc or Zstandard inside an `.h5ad`, bz2 in Zarr) is refused naming it and the dataset, with how to save the file again. A Zarr store on a web host that hides paths starting with `_` or `.` (GitHub Pages without `.nojekyll`) is refused naming `.nojekyll` instead of read with empty cell and gene names. Suo (562k cells × 8,192 genes, CSR `.h5ad`) read locally in Node, against h5wasm: open 210 → 27 ms, first gene 1.9 → 0.85 s, neighbour graph 7.1 → 0.54 s.
- **Web app**: Reading follows the measured link (`formats/io/link.js`): every byte reader measures round trip and bandwidth as it reads, and request size, gap merging, read-ahead and the number of requests in flight follow the measurement and the host's protocol (six to an HTTP/1.1 host; a host that hides its protocol, as GitHub does, starts at HTTP/2's 100 streams until a burst shows otherwise). An `.h5ad`'s size probe carries its first 512 KiB in the same round trip; range reads bypass the HTTP cache so Chromium sends a file's ranges together; a link delivering at least 64 MB/s, a local file and a picked Zarr folder are read ahead like a disk. Suo `.h5ad` over a simulated 80 ms / 50 MB/s link: probe and open 98 ms (1.2× the round-trip bound), dataset ready in the app 1.3 s.
- **Web app**: An `.h5ad` on a web host is opened with 16 KiB metadata pages, since anndata writes each element's headers within a few KiB beside its data: a 1M-cell `.h5ad` reads 1.8 MB before its first view (1.4 s at 80 ms and 2 MB/s).
- **Web app**: A part of an unfiltered HDF5 chunk (anndata's default sparse `.h5ad`) is fetched as the bytes it needs, not the whole chunk: 1,000 cells of a 1M-cell CSR `.h5ad` at 30 ms and 5 MB/s read 4.4 MB (the whole chunks are 121 MB).
- **Web app**: Over a slow link, the first gene of a CSR `.h5ad` keeps the matrix values it reads (its entries spread over the whole matrix), so later genes read nothing (a 1M-cell file at 5 MB/s: the second gene in 6 ms, not 133 MB in 27 s). A cell subset weighs the bytes it would fetch by the link: an analysis of two large groups reads from the gene index in memory (5.5 s) instead of 400 MB of rows (117 s).
- **Web app**: Opening a file reads no category dictionary. Every obs column is described from its metadata, with its category count, and a categorical's categories (and the colours Scanpy stored for them) are read when the column is first shown. A GEO file whose barcodes are stored as a categorical (21,751 categories, 1.1 MB) no longer reads them before its first view.
- **Web app**: Every request into a Zarr store bypasses the HTTP cache (`cache: 'no-store'`), as range reads of a file do: WebKit never finished delivering to the reader worker the bodies of store chunks served from its cache, and a store's metadata is never paired with chunks cached from an earlier version.
- **Web app**: A refusal written for the person (a record with several files, a file with nothing to plot, a folder link with no `datasets.json`, an AnnData too old to read) is shown as it is instead of a generic network sentence. A folder link without `datasets.json` says so for every host: "The linked data could not be opened: https://example.org/atlas/ holds no datasets.json, so it is not a dataset catalog. Link to the folder that holds datasets.json …".
- **Web app**: Analysis computes each gene from its non-zero entries, read as `{length, indices, values}` and never installed in field state: Wilcoxon and Welch statistics treat the zeros as one rank block, bit-identical to the dense kernels for differential expression and Wilcoxon markers (Welch marker p-values within 2.2e-9 relative). Differential expression sends the groups once and each gene as its column; a cell subset is read by whichever of rows or columns costs fewer bytes. The worker pool runs one worker per core less the page and reader (no cap of 8). In Chromium on Suo: marker discovery by `cell_type` (69 groups × 8,192 genes, Wilcoxon) 12.3 → 0.86 s with no long task; differential expression of one 266-cell type against the rest 14.1 → 2.1 s.
- **Web app**: Gene clustering in the Genes panel is exact agglomerative clustering on the analysis workers (Müllner's nearest-neighbour-candidate algorithm, double precision; the same merges as before, trees equal to scipy's to 5e-15): 3,000 genes × 69 groups 5.0 s on the main thread → 0.15 s. Memory is 8 bytes per gene pair (was about 104), so the default budget fits about 11,580 genes instead of 2,270. Gene signatures read each gene once as its sparse column (100 genes 450 → 32 ms), CSV export of page data and signatures is written by the workers, categorical page data is returned as codes, and Correlation computes pages concurrently and refuses a categorical X or Y.
- **Web app**: Box and violin points, and the volcano plot's genes, are one image per trace drawn by a worker at exactly the positions Plotly lays out, instead of one SVG element per point; hover and SVG/PNG export are unchanged. Suo, 544,130 cells in two pages: box with all points 5.0 → 0.13 s, violin 5.2 → 0.23 s, hover no longer blocks for 0.7–0.8 s.
- **Web app**: Figure export computes bounds, draws the point passes, writes vector circles and encodes PNGs in a worker; the page only composes. Where a browser's workers have no WebGL2 -- WebKitGTK and WPE ([WebKit bug 209809](https://bugs.webkit.org/show_bug.cgi?id=209809)), Safari 16.4–16.6 and 17.0–17.1 on macOS Ventura ([WebKit bug 183720](https://bugs.webkit.org/show_bug.cgi?id=183720)) -- the worker says so once and the same code draws the point passes on the page. Output is byte-identical apart from the timestamp, except that Full Vector now honours a depth sort at every point count (it skipped one above 200,000 points). Suo: main-thread long tasks of up to 804 ms → none; 600 dpi PNG about 1.1 → 0.6 s; Full Vector SVG peak memory +239…+581 MB → +15…+114 MB.
- **Web app**: The lasso projects the pane's visible cells once per press and tests only the cells a new vertex can flip, bit for bit the full even-odd test: 0.06–0.11 ms per pointer move on Suo with a 900-vertex lasso (was 0.93–1.3 s). Proximity re-tests only the radius band that changed; KNN expands on the neighbour lists the reader worker built. Selections, highlights, undo history and the Python selection wire share one representation, an ascending `Uint32Array`; selected cell sets are identical before and after over 510 exact comparisons.
- **Web app**: The reader worker builds a neighbour graph's page payload as it checks the matrix: edges arrive already in the renderer's fixed random order as 12-byte records, with one ascending neighbour list per cell, and a transform-feedback pass gathers edge positions on the GPU, so a filter change gathers nothing on the CPU. With edges on (Suo, 6.3M edges, K = 250,000), a legend click takes 6–8 ms instead of 19–57 ms; frames are pixel-identical. Embeddings and vector fields are checked, normalized and bounded in the worker (2D embedding main-thread time 71.6 → 7.3 ms), and no spatial tree is built on the main thread.
- **Web app**: Colouring and filtering touch each cell once, with byte-identical output over 45 interactions. Main-thread time on a 562,000-cell dataset: picking a loaded gene 18 → 5 ms, a legend click 12 → 6 ms, an outlier-slider step 22–24 → 9–12 ms, keeping a view 56 → 40 ms. A view whose filter inputs match an array another view holds reuses it, and GL error reads happen only after new storage is allocated, so restoring a 7-view session runs with no task over 50 ms. Smoke rebuilds only when its inputs change.
- **Web app**: Outlier quantiles are read from the file, one field at a time, on the filter's first use; for any other file or latent space they are computed in the reader worker, every field asked for together in one pass over the latent space, which is then kept within its budget. Quantiles that no longer fit the file's categories are computed instead. The latent space is no longer downloaded when a dataset opens. Sessions never drop an outlier threshold silently: one that cannot be restored is left at 100% and named.
- **Web app**: Memory budgets are disjoint shares of one allowance, half the estimated device memory (`navigator.deviceMemory`, raised past an 8 GiB report by four times the JavaScript heap limit; 8 GiB where neither is reported). Decode workers are sized by cores and that budget; switching expression matrices keeps each built gene index while they fit, so switching back reads nothing; the default analysis budget is a sixteenth of device memory.
- **Web app**: A server or notebook dataset (`cellucid serve`, `serve_anndata`, `show_anndata`) is opened with the browser's AnnData reader, as a linked file is, under the server's dataset id, name and launch defaults; its dataset list is described from metadata and a dataset is opened when chosen. The `remote://` and `jupyter://` schemes and the server-side mapping routes are gone.
- **Web app**: Sessions and community annotations follow the dataset, not the route: a dataset opened from a sample, a link, a server or a local file is the same dataset. A local file is named by `uns['cellucid']['name']`, else `uns['title']`, else its file name.
- **Web app**: The sidebar gene search and the analysis gene selector share one windowed list (`app/ui/components/windowed-listbox.js`) with the same keys: arrows wrap at either end, Home/End jump to the ends, and every option carries `aria-posinset`/`aria-setsize`. Field renames, deletions and restores touch only what changed (Suo: gene delete 20 → 2.3 ms).
- **Datasets**: The datasets are built by `scripts/build_atlas.py` and `scripts/build_pancreas.py` through `scripts/pinned_build.py`, which pins the environment, inputs, analysis and the producer modules `prepare()` executes, reads every staged store back with anndata, and compares its digest with `sources/<id>.json` (schema 4). Starting views were re-captured on each store as format-2 sessions; new `scripts/capture_session_preset.mjs`, `scripts/verify_catalog.mjs` and `scripts/local_site.mjs`. CI opens every store with anndata 0.12.19.
- **Reproducibility**: `cellucid-reproducibility` states that every recorded run measured the export-based Cellucid (web `0ab545aed`, cellucid-python 0.9.1) and what a re-run on the current version needs; its setup checks out the measured commits. Nothing was re-run.
- **Web app**: The Session panel is redesigned. The **Dataset** box describes only the dataset that is open, rewritten on every dataset change: its name, cells and genes, format and host (`Zarr v2 folder · theislab.github.io`, `h5ad file · this computer`), its own description (two lines, **More**), how much has been read in place, and a collapsed **Details** with source and citation links, file, size, obs fields, embeddings, layers, raw, neighbour graph, latent space and vector fields; nothing is cut mid-word. One **Open…** button opens a File… / Folder… menu (a browser dialog picks files or a folder, never both) in place of the **H5AD**, **Zarr ZIP** and **Zarr folder** buttons. One **Link or server** address bar with **Open** replaces **Server** and **Open from link**; what is typed decides (a port, `host:port` or `localhost` is a server, a provider form or a file address is a link, any other web address is a server when one answers `/_cellucid/health`), and a connection shows as one line with **Disconnect**. The **Transfers** block is removed: how much has been read is one line in the Dataset box, and a retried read is a notification. The sample list shows `—` while a file, a single-dataset link or a server is open, and lists a connected catalog's datasets, named "Dataset from <host>", instead of mixing sources. **Data mapping** no longer opens on its own: the best-guess notification's **Show** opens it, its rows are one flat list, and arrays that are off and wider than three columns, or cannot be drawn, are listed under a collapsed **Other arrays**. Without a latent space the outlier filter is one **Choose latent space** button; with one, the slider and at most one short line. A Zarr v2 folder chosen in Safari, whose folder dialog leaves out hidden files, is refused with how to open it instead (drag the folder onto the page, or its `.zarr.zip`).
- **Web app**: While the first dataset is still opening, the sample list already offers every entry of `datasets.json` under the name its path gives it, the one being opened marked *· opening…*; names and cell counts follow once the files are read.
- **Web app**: A dataset that does not open is said once, in one notification of one sentence (*Could not open "Garcia".*), its raw error folded under the notification's **Details** and in the console; the explanation and **Try again** stay next to the selector. Opening and connection failures keep their raw error under **Details** too, never in the sentence.
- **Web app**: The Highlighting tool's gesture line (and *Color by a field first…*) and the Antialiasing line (*Chosen automatically: …*) are the first line of their (i) tips; under the Tool buttons the sidebar shows only what a gesture did.
- **Web app**: The *Signed in to Google* line with **Sign out** goes when the sign-in expires: the token store says when a sign-in is stored, removed or expires.

- **Web app**: A link to a host that ignores range requests opens: its file is downloaded once into memory, with its progress shown, and opened from there. A host that does not let other sites read its files says so and what to do: "`<host>` did not let this site read `<file>`: it does not allow other websites to read its files (CORS), or it could not be reached. Download the file and open it from your computer instead."
- **Web app**: Opening a link starts the AnnData reader worker while the link is resolved, so the worker and its modules are loaded by the time the file is found. A CELLxGENE link probes the file in the Census release the browser read last beside the release list, so it usually opens one round trip sooner; a newer release is followed. A Zenodo link that names its file asks for the file's first bytes beside the record, so it opens one round trip sooner too. Opening a link again in the same page (Disconnect and open, a session reopening its dataset) asks its provider nothing and starts from the first bytes read the first time.
- **Web app**: Analysis windows and floating sidebar panels resize from all four edges and corners, by mouse, pen, finger and keyboard (arrow keys on the focused edge, 16 px or 64 px with Shift, Home/End to the least and largest size; `Alt`+arrows move the window). They are held only to their content's least size and the viewport: the 500 px height cap and the twice-the-opening-width cap are gone. Edges snap to other windows, the sidebar, the middle of the free area, the screen and the content height, with a guide line, and show the size while they move. A double-click on an edge fits it; **Maximize** (or a double-click on the header) fills the free area beside the sidebar. A new window of an analysis opens at the size the last one was given. On a screen 720 px wide or less, analysis windows are bottom sheets at 40%, 75% or 100% of the screen, with a grip, a collapsing header and a tab strip for the others.
- **Web app**: The expanded analysis view moves and resizes by mouse, pen and finger (it was mouse-only) from every edge and corner, its splitters and bottom-right corner answer arrow keys, a double-click on its header maximizes it, its least size is its content's (not 600 × 500), and its size and splits are remembered per analysis. On a screen 768 px wide or less it fills the screen with a touch grip on its split.
- **Web app**: The Analysis panel groups its rows: **Genes** (Marker genes, Differential expression, Dot plot, Gene signature), **Cells** (Composition) and **Explore** (Quick, Detailed, Correlation).
- **Web app, breaking**: The `analysis/windows` session chunk is `{windows, sizes}`: each window records its own box (not the header-high box of a collapsed window, which it used to save), whether it is open and its stacking layer, and `sizes` holds the remembered window and expanded-view sizes per analysis. A chunk without them is refused.
- **Web app, development**: The web contract suite needs Node.js 22.15 or newer (`engines.node`); its Zstandard contracts compress their inputs with `node:zlib`'s Zstandard encoder, which Node.js 20 (end of life) lacks. CI runs Node 22, 24 and 26.
- **Web app**: Marker genes compares every group of a field against the rest or a chosen reference group, for every group or a chosen subset, inside all cells, the visible cells or a page. A group under **Min cells per group** is never an error: it is skipped and listed (default), pooled with the other small groups, or tested anyway and flagged low n, and chips on the result name them. The test follows scanpy's `rank_genes_groups`: Wilcoxon (score = the z of the rank sum, with AUROC), Welch's t-test or the t-test with overestimated variance; p-values exact or continuity-corrected (SciPy, R) or scanpy's normal approximation with or without tie correction; Benjamini–Hochberg or Bonferroni over the genes tested; fold changes as stored or with `log1p` undone (automatic when the data mapping marks the matrix `log1p`). The result keeps every gene's statistics for every group, so direction, filters (adjusted p, |log2FC|, % in group, % in reference), ranking (score, adjusted p, log2FC, AUROC, % difference) and the number shown re-apply without testing again; it shows as a list, a heatmap or a clustered heatmap. Genes can be excluded by family (mitochondrial, ribosomal, hemoglobin) or pattern, or limited to the highly variable genes or a list. It exports the shown markers and every gene as CSV with scanpy's column names plus `auroc`, `small_group` and `tested`, and the markers as GMT; **Score as signature** and **Dot plot** hand them on. The statistics match scanpy and SciPy on the shared fixture.
- **Web app**: Differential expression compares two pages or two sets of categories of a field (disease against healthy), optionally within the visible cells, a page or categories of another field, with **Swap A and B**. A group under **Min cells per group** is tested and flagged low n instead of refused; only overlapping groups and a group under 2 cells disable **Run**, with the reason. It has Marker genes' test options, free adjusted-p, |log2FC| and % expressing filters that redraw the volcano without testing again, and exports scanpy's column names with `pct_nz_group`, `pct_nz_reference`, `auroc`, the group means and counts.
- **Web app**: Gene signature scores one or more named gene sets — typed, pasted with any separators, imported (`.txt`, `.csv`/`.tsv` in long or wide shape, `.gmt`) or taken from Marker genes — by mean, median, mean z-score, scanpy's module score (`tl.score_genes`: control genes drawn from expression bins with a seed) or UCell (ranks, with a max rank). A name is matched exactly, or else to the one gene matching it ignoring case; names the dataset lacks are left out and listed, never an error. Rescaling (z-score, 0–1), the plot and the layout of several sets redraw without scoring again. The CSV is `page,set,cell_index,score`. Module score matches scanpy and UCell a numpy reference on the shared fixture.
- **Web app**: Gene-streaming analyses (Marker genes, Differential expression, Gene signature) size themselves: a run reads ahead as many genes as the workers consume while one gene is fetched, from the measured link and each gene's measured size and test time, within the analysis memory budget, and adapts while it runs (more read-ahead when the workers wait on reads, kept only when it makes the run faster; less on a throughput drop or memory pressure). Analysis workers start one at a time, so building the pool never holds the page. How many requests go out at once is the reader's scheduler's to decide. Nothing it decides changes a statistic.

- **Web app**: One design for every screen of the site: the startup failure screen, the code-load failure screen, the no-JavaScript screen, the 404 page, the privacy policy and the sign-in return page are the app's brand card (wordmark, heading, one sentence, clear actions) in the theme chosen in the app; the standalone pages load a small `page.css` instead of the whole app stylesheet. The design principles are in `assets/css/README.md`.
- **Web app**: The startup failure screen says what kind of problem it was (the data could not be opened, a server did not answer, this browser cannot run Cellucid, Cellucid could not start), folds the exact message (address, host, status) under **Details**, and offers **Try again** and **Open Cellucid without this link**. A module that fails to download shows the same screen, not only a red status line.
- **Web app**: Data a link names that does not open no longer stops Cellucid. A `?data=` link, a `?remote=` server, a `?dataset=` or a catalog's default dataset that cannot be read leaves the app running without it: a panel over the empty viewer says the data could not be opened, with the exact message under **Details** (and in the notification log), **Try again**, and **Open without this link**, which removes those settings from the address. The full-screen card is kept for the app itself failing (no WebGL2, its code did not load, its own server). The panel says the reason itself, as its notification does ("This 10x Cell Ranger .h5 file has no cell embedding …"; "… holds no datasets.json …"; "The owner of x.rds has turned off downloading for it"), and when asking again cannot change the answer (a file or folder that is not there, a refusal, a file Cellucid does not read), **Open without this link** comes first, before **Try again**. Google Drive and OneDrive refusals and a Drive read held back are shown as written wherever they arrive.
- **Docs**: The documentation site takes the app's look: Inter, the app's ink and cyan colours in light and dark, square surfaces, the app's favicon, and the wordmark without a duplicate "Cellucid" beside it.
- **Docs**: Every screenshot is recaptured from the current app by its scenario, with a `captures.json` record naming the clean web commit it shows, so none is left without provenance; captions and alt text describe what each image shows. New figures for the R objects, Loom, Cell Ranger `outs/`, MatrixMarket and `.h5mu` Modality pages, the opening progress and failure panels, the Google Drive and OneDrive sign-in lines, Composition and Dot plot, the R viewer, and the phone layout. The JupyterLab images are replaced by the viewer `show_anndata()` serves and a transcript of a real `viewer.debug_connection()` report, and the Firefox and WebKit captures are gone. The screenshot tool takes a per-scenario `viewport` (with `touch` for a phone) and grows a crop to hold the whole drawn pointer.
- **Docs**: A new home page: the tagline and one sentence, **Open Cellucid** and the 60-second tour as calls to action at a documentation size, the formats it opens, a fresh light and dark capture of the Suo sample in 3D, the three ways in (web app, Python, R), four features with their own captures, compact capability and learning cards, and citation and contribution; it has no sidebars. The outdated close-up screenshot is gone. The images are WebP captures from new `home-*` scenarios of the screenshot tool, which now writes WebP (`format: 'webp'`, a fixed `quality`) as well as PNG, with provenance in `captures.json`.
- **Web app**: One `favicon.ico` (the root one) for every page.
- **Web app**: A published sample view that cannot be applied says so in one sentence; the reason is under the notification's **Details**.

- **Python, R, web app**: Highlights and visibility changes carry their cells in a binary frame, not a JSON list. `viewer.highlight_cells()` and `viewer.set_visibility()` take a list of ints or a one-dimensional NumPy integer array (the fast form), in any order, each once; the cells are checked and sorted once and written straight into the frame. `set_visibility([])` is refused (pass `None` for every cell). Two million cells are a 7.6 MB frame instead of 15 MB of JSON: from NumPy, Python builds it in 3 ms (the JSON took 0.35 s) and the viewer has applied the highlight 0.22 s after the call; ten million cells, 13 ms (1.7 s) and 0.86 s. R builds the frame of ten million cells in 16 ms (the JSON took 4.2 s).
- **Python**: `viewer.send_message()` returns `True` when the command was sent and `False` when it waits for the viewer to connect; it raises only after `viewer.stop()`. `get_session_bundle()` raises at once when the viewer reports it could not capture the bundle, or sent an invalid one. Hook callbacks run on the viewer channel's thread, one event at a time and in order; a hook that raises is logged and the channel stays open. The event is recorded in `viewer.state` and handed to `wait_for_event()` before any hook runs, and the other hooks still run: one buggy `on_message` hook no longer drops the event, so `get_session_bundle()` no longer times out behind it (several failures are logged together, as one `ExceptionGroup`). `debug_connection()` reports `channel_connected`; the `debug_snapshot` event no longer has `parentOrigin`, and no event carries `viewerId`.
- **Docs**: "Architecture: message routing (HTTP vs postMessage)" is "Architecture: the viewer channel" (`e_jupyter_hooks/05_architecture_viewer_channel`).
- **Web app**: In Full quality, the highlight pulse redraws only the layer over the points. While highlights pulse, each pane's base layer (backdrop, grid, every point, the vector-field overlay) is kept and every pulse frame draws the highlights, edges, centroids and orbit anchor over a copy of it, pixel for pixel the frame drawn whole (`rendering/pulse-layer-cache.js`; `tests/browser/highlight-pulse-layer-parity.spec.mjs`, one view and a grid, three engines). GPU time of a pulse frame (Chromium, Metal, antialiasing off): 10M cells with 20,000 highlighted 9.97 → 0.48 ms, with 5M highlighted 38.9 → 24.1 ms; 1.5M cells with 15,000 highlighted 4.91 → 0.69 ms. With antialiasing on, a multisampled pane keeps its own samples instead (WebGL2 cannot copy one multisampled buffer into another): the ring pass marks the samples it covers in the stencil, and a pulse frame resets and draws again only those, every point still submitted but only those under a ring shaded (`rendering/pulse-sample-mask.js`; the same spec, the pulse moved between two phases, three engines). 10M cells, antialiasing on: 20,000 highlighted cells in one region 40.3 → 7.5 ms, 300,000 in one region 70.0 → 15.7 ms, 20,000 scattered over the embedding 37.2 → 25.3 ms, 1M scattered 44.3 → 45.0 ms (their rings cover the embedding).
- **Web app**: Applying a highlight costs the highlighted cells, not the dataset: the page counts no longer copy a page's cells into an array to count them, nor to digest them (5M of 10M cells: 85 ms of the apply); 1M highlighted cells of 10M apply in 24–26 ms (57 ms).
- **Web app**: Startup no longer runs as one ~90 ms task: it yields between its first stages (the interface, the WebGL programs), and the page's first layout (about 25 ms) is done once, where the welcome dialog takes focus, instead of also being forced by checks that measured the page (the dialog's focus check and the notification center's empty card stack). Chromium: no main-thread task over 50 ms while the app starts. A PNG figure's pixels are read in bands of rows with a yield between them, so the read-back of a 300 dpi figure is no longer one long task.
- **Web app**: Publishing a large dataset no longer holds the main thread. Its GPU stores are allocated in a task of their own just before the publication, once the GPU has done the work queued before (a WebGL allocation is checked by a `getError` round trip that waits for the GPU: up to 30 ms of the publishing task), so the publication reads nothing back; its positions and colors go to the GPU with the first frame (the colors were sent twice); and each per-cell array is written once: the first visibility pass adopts the new dataset's all-shown visibility, clearing the active field leaves gray colors gray and republishes no unchanged transparency, the highlight list of a dataset with no highlight is empty instead of being derived from ten million bytes, a fresh mask is proven binary as it is made, and the 0/1 per-cell bytes (validation, the R8 alpha texture of the live and kept views) are read four cells to a word. The first slice of positions for the spatial-tree worker is copied in its own task. Chromium (Metal), 10,000,000 cells in 3D: opening 227–237 ms → no task over 50 ms, switching to it 337–349 ms → none (largest 45 ms).
- **Web app**: Coloring by a category of a large dataset no longer holds the main thread: before the field is shown, its per-category totals and its centroids for the shown representation are computed in slices of about 8 ms (the same additions in the same order, so the same centroids), the codes are proven against their categories once when they load, the colors are written without the transparency the visibility pass writes next, and the visible count of each category is read from the totals when the only filter is the field's own (a category is then all shown or all hidden). Hiding or showing one category reads no mask when every cell is shown. Chromium, 10,000,000 cells: coloring by a 30-category field 187–195 ms → no task over 50 ms; hiding one category in the legend 319–332 → 38 ms; showing it again 232–237 → 42 ms.
- **Web app**: A 10M-cell tree landing no longer holds the main thread: the spatial-tree worker sends the tree as typed arrays (a structured clone of its 38,568 node objects took 35 ms in one task; now 10–11 ms with the install), and the store's adoption of the tree's draw order copies its source words on the GPU, drawn from the store into the source texture, instead of uploading 120 MB of positions and 40 MB of colours from the page (adoption 32–35 → 14–19 ms). Chromium (Metal), 10M cells: no task over 50 ms from opening to the tree in use (was a 77–90 ms frame of two tasks, once one 101 ms task).
- **Web app**: A highlight of one cell list costs that list: a page with one enabled group is painted (runs of consecutive cells as one fill) and published from the group's own list, counted as its length, and a cell list is proven ascending once (`app/utils/cell-set.js`), so the page selectors and analysis no longer walk it on every change; the highlight layer writes each record's color as one word, in a packing loop of its own. Chromium, 10,000,000 cells, highlighting 5,000,000: 114 → 47 ms, next frame 72 → 30 ms; all ten million: 177 → 72 ms, next frame 164 → 30 ms. A highlight of many cells is drawn from the view's point store, already on the GPU in draw order, with a colour word a slot written on the GPU from the highlight and visibility bytes, whenever that writes fewer bytes than records of its own (the CPU gathered ten million positions at random: 72 ms, and a new 160 MB store was cleared by the GPU process before it could be checked: 54 ms); the pulse frame of ten million highlighted cells takes 9 ms of GPU time instead of 54 (`tests/browser/highlight-store-draw-parity.spec.mjs`, three engines).
- **Web app**: Opening an analysis loads Plotly at the page's next idle moment, so its evaluation (one script, about 190 ms of main thread) runs while the analysis is set up instead of inside the first run: a differential-expression run started two seconds after its panel opened no longer has the 234 ms task (Chromium, 1.5M cells).
- **Web app**: A MatrixMarket folder opens in about one and a half times the time plain inflate of the same bytes takes: the entries are parsed in blocks of about 1 MiB on the decode workers into typed arrays and collected a block at a time, a large gzip text or R file is decompressed on a worker of its own (`formats/rds/gunzip-worker.js`; `DecompressionStream` inflates on the calling thread), the transcode space is written 1 MiB at a time, and the cell tables are read while the matrix streams. 100,000 cells (44M counts): Chromium 1.6 s (1.5× plain inflate), Firefox 1.7 s (1.5×), WebKit 1.0 s (1.33×); 500,000 cells in Chromium 7.6 s. `npm run bench:formats` runs the decode workers on `worker_threads` (Node, 100k cells: 1.80×). An entry whose row or column is not a whole number is refused.
- **Web app**: One shared module (`utils/exact-record.js`) holds the record and key checks the app's modules each defined for themselves (about forty copies, with different messages and different ideas of a record): `requirePlainRecord`/`requireKeys`/`requireExactKeys` for data from JSON, a session or a worker, `requireObjectRecord`/`requireObjectKeys` for an owner the app built, `requireClonedRecord` for worker payloads. A key mismatch reads the same everywhere: "… must contain exactly: a, b, c (missing: b; unexpected: d)." The community-annotation wire contract keeps its own path-addressed messages, the ones the annotation repository's validator writes.
- **Web app**: One formatter for sizes and durations (`utils/format-units.js`): the notification center wrote "2 MB" and "850ms" where the dataset card and read progress wrote "2.0 MB" and "850 ms". Every message now writes "512 B", "4.2 KB", "12 MB" and "850 ms", "4.2 s", "3 min 5 s", and no unit rounds into the next one ("1024 KB", "1000 ms", "60.0 s").
- **Web app**: A file whose bytes are not HDF5 names the format it was opened as: "The selected file is not a valid HDF5 file, so it cannot be a Loom (.loom) file. Write it again with loompy …"; a damaged file or Zarr store gives its own format's route to write it again.
- **Web app**: The viewer channel's modules and names are neutral: `data/viewer-host-source.js` (`ViewerHostDataSource`), `app/viewer-host-command-handler.js` and `app/viewer-host-pointer-delivery.js` (they were `jupyter-*`, `JupyterBridgeDataSource`, for a channel that serves Python, Jupyter and R hosts), and their messages say "Viewer host …". The `jupyter` source type of sessions and analytics, and the `jupyter` runtime, are unchanged.
- **Web app**: The "Link or server" ⓘ tip names Figshare, Google Drive and OneDrive, and every other link provider (a contract holds it to the provider registry). The cell-info tip and page state the spatial-index threshold the viewer uses (views of 250,000 cells or more, built in a worker); they said "above two million cells" and "above 256 MB of positions no background build".
- **Web app, development**: `npm run bench:formats` measures every non-AnnData format (Loom, a Cell Ranger `outs/` folder, a MatrixMarket folder, Seurat and SingleCellExperiment `.rds`, a Seurat `.rda`, `.h5mu`, MuData Zarr) against the same data as a Cellucid-ready `.h5ad`: first view, first gene, peak memory, and a stream format's pass against plain inflate of its own bytes (`scripts/benchmark-formats/`).
- **Web app, development**: A monthly CI job (`.github/workflows/r-corpus-drift.yml`, also on changes to the R readers; locally `sh scripts/rds-corpus/drift.sh <Rscript>`) writes the R object corpus with the current CRAN SeuratObject and Matrix and Bioconductor SingleCellExperiment, reads every file with the transcoders against what R read back, and fails on a slot those packages add or remove (`scripts/rds-corpus/known-slots.json`).
- **Python, R, web app**: An h5Seurat file is refused the same way everywhere, in the words of the input format registry (`refusal`, in all three repositories' `input_formats.json`): "This is an h5Seurat file (SeuratDisk), which Cellucid does not open. In R, read it with seu <- SeuratDisk::LoadH5Seurat(path) and save it with saveRDS(seu, "object.rds"), then open the .rds, or show it with cellucid::cellucid_view(seu)." `cellucid serve x.h5Seurat` and `cellucid_serve()` refuse it before serving.
- **Python**: Stopping a server closes the connections it still holds, so a viewer started next on the same port never answers from the previous dataset.
- **Web app**: A gzip or zlib Zarr chunk, and a ZIP entry, are decompressed with their input fed in 32 KiB slices, so the size the format declares bounds the output while the chunk decodes, in every engine (WebKit's `DecompressionStream` inflates everything written to it at once). The session codec and the format readers share this one input feeder (`formats/io/inflate.js`).
- **Web app**: A read's progress counts the bytes of every request as they arrive, so the progress, rate and time left of a long request (an 850 KB gene-name heap over a slow link, an R object or a MatrixMarket folder read once) stay live; the reader worker posts a running read's progress at most every 250 ms and always posts the last change of an interval, and the bytes of an attempt that fails are taken back. A GEO sample file its host sends compressed is downloaded with its bytes counted.
- **Web app**: Picking or dropping a file, or making a connection, is always a new opening, so the file already open opens again (a newer copy with the same dataset id is shown); only choosing a listed dataset again keeps what is open. An `.rds` opened again in the tab reads what was kept the first time.
- **Web app**: A format refusal's advice is what can be done from where the file was opened (the reader worker's `open` is told `fromLink`): a 10x Cell Ranger `.h5` with no embedding opened from a link says to download the whole `outs/` folder and open it from this computer, or to compute an embedding (a link opens one file, and no provider shares a folder Cellucid reads); a linked molecule table says to link `filtered_feature_bc_matrix.h5` instead; a 10x MatrixMarket matrix from a GEO sample, which has no `outs/` folder, says to compute an embedding; and an `outs/` folder whose `analysis/` has no projection says Cell Ranger ran without its secondary analysis.
- **Web app**: An address typed without `http://` on this computer (`127.0.0.1:8000/atlas/`), or a server's address, that answers but is no Cellucid server (a static web host serving a folder: no `/_cellucid/health`, or a web page there) is said to be no Cellucid server, with how to open such a folder: the address of its `datasets.json`.
- **Web app**: The startup phase that opens the dataset is published as `dataset` (`data-startup-phase`, and the `startup_phase` analytics field), not `identity`.
- **Web app, development**: `npm run bench:anndata` benchmarks the AnnData reader, and `scripts/benchmark-real-links.mjs` opens files from this computer through the File… and Folder… buttons (`{local}/path`), serves a link over HTTP/2 (`--link-protocol h2`), and reports each step's peak browser memory and main-thread long tasks.
- **Web app**: Opening a link (`?data=`, the address bar) or a server reads no sample file: the dataset list names the samples from `datasets.json` (one request) and reads their files (title and cell count, about 85 requests) only once the list is reached for. A page opened on a sample, or with nothing connected, lists them as before.
- **Web app**: The address bar names a link the page was opened with in its provider's form, as one opened from the bar is named: a pasted Google Drive link becomes `?data=gdrive:<id>`, a Figshare page `figshare:<id>`, a OneDrive share link `onedrive:<drive>/<item>`, a DOI of another publisher the Figshare item it names, and a Dropbox link its share link (`dropbox:https://www.dropbox.com/…`).

### Removed

- **Python**: The prepared export format and everything that wrote or served it: `prepare()`'s `latent_space`, `out_dir`, `obs_categorical_dtype`, `var_quantization`, `obs_continuous_quantization`, `compression`, `var_gene_name_columns`, `source_citation` (now `citation_text`/`citation_url`), `created_at` and `centroid_outlier_quantile` arguments; `Embedding` (store the whole array and choose columns with `embedding_defaults={"3d": ("X_pca", [0, 1, 2])}`); `generate_datasets_manifest()`; `serve()`, `CellucidServer`, `show()` and `CellucidViewer`; the prepared-export branch of `cellucid serve`; the export lock, transaction journal, manifests, quantization and payload helpers.
- **Python**: `AnnDataAdapter` and the server-side virtual export it computed (manifests, payload and mapping routes, server-side gzip, centroids at startup), the backed read mode, `serve_vector_fields` / `cellucid serve --vector-fields` (the viewer reads a vector field when its overlay is shown), `/_cellucid/health`'s `is_backed`, and `centroid_outlier_quantile` from `serve_anndata`, `show_anndata` and `cellucid serve`.
- **R**: `cellucid_prepare()`'s `latent_space` (name an embedding with `latent_key`), `out_dir`, `obs_categorical_dtype`, `var_quantization`, `obs_continuous_quantization`, `compression`, `var_gene_name_columns`, `source_citation`, `created_at` and `centroid_outlier_quantile` arguments, `cellucid_embedding()`, and the export lock and transaction journal.
- **Web app**: Prepared exports (`dataset_identity.json`, the obs, var and connectivity manifests, quantized and gzip payloads), catalogs of exports, `state-snapshots.json`, the **Prepared** button and the prepared-folder drop; with them gene download-ahead and its payload cache, the quantization worker, the data-mapping "Read from file" section and the analytics `data_version` parameter.
- **Web app**: `defaults['centroidOutlierQuantile']` (a file or server naming it is refused as an unknown default) and the per-centroid outlier array of views and kept-view snapshots; the marker cache's IndexedDB tier (`cellucid_marker_cache`) and its `maxAgeDays`/`cacheVersion` options; the bulk gene-expression cache and its API (`getGeneExpressionSubset`, `fetchBulkGeneExpression`); the data layer's session cache export and import.
- **Datasets**: `exports-v3/` in `cellucid-datasets` (299 MB of prepared exports and their snapshots) and the eight preparation notebooks; `exports/` and the Zarr v3 folder in `cellucid-demo-custom-datasets`.

- **Web app**: The per-file notices "Opened H5AD directly in the browser…" and "Opened a Zarr … directly in the browser…", the loaders' own download notifications and their `showProgress` option, and `ZarrDataSource.loadFromFileMap` (nothing called it): the opening progress shows every file.
- **Web app**: The analysis **Performance settings** (batch size, memory budget, network and compute parallelism, and their presets); Marker genes' mode select (Ranked, Clustered, Custom genes), its **Use cached results** checkbox (the cache is always on, marked **Cached** with **Recompute**) and its ranked and heatmap CSV files; Gene signature's sum; the run-cost line under the run button.

- **Python, web app**: The notebook `postMessage` command path (the script in the notebook output, `AnnDataViewer.viewer_origin`) and the routes `POST /_cellucid/events`, `/_cellucid/selection` and `/_cellucid/session_bundle`: commands, events, selections and session bundles travel on the viewer channel, and the server answers every `POST` with 405. The viewer's health poll every 3 s is gone: the host closing the channel is the stop signal.
- **Web app**: The remote source's unused update socket (`ws_port` in `/_cellucid/info`) and the connection-loss recovery that only it could trigger; a server advertising `ws_port` is refused.

### Fixed

- **Python**: The viewer's web build is checked once per file version: every asset request re-read the inventory and hashed the asset under one lock, 30 ms each (6 s for 200 modules before a first frame); now 0.2-0.75 ms.
- **Python**: An error rebuilt from a `DataMappingError` raised `TypeError` instead of the error.
- **Python**: Stopping a viewer ends every wait still pending on it (`get_session_bundle`, `wait_for_event`, `set_data_mapping`) at once with `RuntimeError` naming the stop; a wait without a timeout used to hang forever.
- **Python**: `get_session_bundle(timeout=None)` waits for the bundle as long as it takes: its request stayed answerable for an hour only, so a bundle that arrived later was refused while the call kept waiting.
- **Python**: `CellucidSessionBundle.save()` writes atomically and saving a bundle onto its own file leaves it as it is.
- **Python**: The viewer channel closes the RFC 6455 way, so the page reads the close code and reason it was sent. The server shut its socket down right after its close frame, with what the page was still sending (the rest of a text frame refused as over 1 MiB, its own close frame) unread, which resets the connection; on Windows the reset discards the close frame before the page reads it, so the page saw an abnormal closure (1006) instead of 1009, 1008 or 4401. Every close (a broken rule, a refused hello, a replaced or stopped viewer) now sends the close frame, shuts the write side, and reads and drops what the page still sends until it closes its side or 5 s pass (`CLOSE_DRAIN_SECONDS`), on the connection's own thread; a close made from another thread reaches a thread waiting to read within 0.25 s. Stopping the server closes an open channel with 1001 ("The server stopped"); it sent no close frame at all. A session bundle cut off by a dropped connection no longer leaves its partial temporary file behind.

- **Web app**: On a portrait screen (a phone held upright) the camera no longer cuts the data off at both sides. The 45° field of view applies to the plot's shorter side (`rendering/viewer/field-of-view.js`), so every camera — the default of 1D, 2D and 3D, a sample's preset, a session's — shows across a portrait screen what it shows down a landscape one; landscape screens draw exactly as before.
- **Web app**: A count and its noun agree everywhere ("1 cell", "1 row", "1 entry"; the legend read "1 cells"): every count with a noun goes through `formatCountOf` or `countNoun` in `utils/format-units.js`, with `formatCount` moved there from the sidebar's formatter so every layer uses one; numbers are grouped in English whatever the browser's locale (a refusal read "40.000.000 points" on a German or Turkish machine). `tests/count-noun-agreement-contract.test.mjs` refuses a `${n} cells` template, a hand-rolled `n === 1 ? '' : 's'` and an unlocalised `toLocaleString()`.
- **Python**: Messages, errors and warnings make a count and its noun agree ("1 row", "1 column", "1 cell"; they read "1 rows"): the noun comes from `cellucid._plural.plural`, and `tests/test_count_noun_agreement.py` refuses an f-string that writes a count before a plural noun and a hand-rolled `"" if n == 1 else "s"`.
- **Web app**: A sign-in line that links the administrator steps (a Google Workspace account its administrator blocks, a Microsoft organisation that lets only an administrator approve Cellucid) shows its reason on a row of its own, then the link, **Use another account** and the (i); beside them the reason was squeezed into a column one letter wide.
- **Web app**: Links in the light theme (dataset details, ⓘ tips, link buttons) are cyan-700, 5.9:1 on white, instead of cyan-600 at 4.1:1, below WCAG AA; the contrast contract now measures the link tokens.
- **Web app, Python, R**: Stored zeros in a sparse neighbour graph are not edges: an `obsp['connectivities']` (or a `prepare()` / `cellucid_prepare()` graph) with explicit zero entries refused the whole dataset; the zeros are dropped. Duplicate coordinates are still refused.
- **Web app, Python**: "Cellucid found nothing to plot" ends with how to fix it: "Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again." A file with no cells says so ("This AnnData has no cells (n_obs is 0), so there is nothing to plot. …") instead of asking to try again, and a catalog dataset that cannot be read names its file and catalog entry.
- **Web app**: Opening another dataset with kept views or a gene shown, or after a session restore, works; it failed with "cannot replace … live points … while snapshot views remain published", and the half-done rollback left the old dataset open. A replaced dataset is released: after a switch from Suo the page's array buffers are back at the pre-open 9 MB.
- **Web app**: A deleted WebGL texture or program is never bound again. A smoke rebuild deleted a texture an owner later restored as its saved binding, raising INVALID_OPERATION (0x502), which a dataset switch then reported as its own failure; every texture and program is deleted through one tracked path, and restores drop deleted objects.
- **Web app**: No request to a web host is left hanging. Each request is made again after a transient failure, with growing waits or as long as `Retry-After` asks: HTTP 408, 429 or 5xx, or a failed or dropped connection. A 503 that GitHub Pages sends without CORS headers arrives as a failed connection. A request whose answer or body stops arriving is aborted and made again. It counts as stalled when neither it nor any other request of its link has received anything for TCP's retransmission timeout (RFC 6298), measured on the link. A failure that outlasts its repeats names the request's URL ("…/obsm/X_umap_3d/3.0: the host could not answer (HTTP 503). Asked 7 times."). Each repeat is logged as a `[cellucid:transfer]` warning naming the URL. A failure that will stay is reported at once: the same error status twice (Zenodo's 500 for a missing record), or no answer twice (a host that does not let this site read it). Asking 7 times took about 16 s. Catalog files and the CELLxGENE, GEO, Hugging Face and Zenodo APIs use the same rules. Cancelling still ends a request, and any wait to repeat it, at once.
- **Web app**: The startup error reads "Startup stalled while opening the dataset (its metadata, field lists and the embedding its first view draws): nothing arrived from the host of dataset "suo" for 60 s, while reading …". It names the reads still open and the last failed request with its URL. Opening a dataset at startup now stalls only when nothing at all arrives for 60 s. Before, any open that took 60 s failed, even one still receiving bytes over a slow link.
- **Web app**: A dataset opens on the embedding its first view draws. The other dimensions named in `uns['cellucid']['defaults']` are read when first shown, not one after another before anything is drawn. The opening embedding is read beside the field lists. A consolidated Zarr v2 store on a web host opens with one request (`.zmetadata`; `zarr.json` is asked for only when it is absent), so no request at open answers 404. The sample list reads the other catalog datasets only after the dataset being opened has opened, and lists that one from its opening, without reading its store twice. Suo from the sample catalog over HTTP/2, first render before → now (median of three, WebKit / Chromium / Firefox): at 30 ms and 5 MB/s, 3.3 / 3.1 / 3.5 s → 1.9 / 1.9 / 2.2 s; at 80 ms and 2 MB/s, 6.9 / 6.8 / 7.2 s → 3.8 / 3.7 / 4.1 s. Before first render: 139 → 27 requests, 10.9 → 5.5 MB, six 404s → none. The 0.9.4 app with the prepared export took 2.1 / 2.0 / 2.4 s and 4.3 / 4.2 / 4.5 s (11 requests, 6.3 MB).
- **Web app**: A gene-column index that a closed or crashed tab spilled to the origin-private file system is removed when the app starts, not only when the next index is written.
- **Web app**: A file opened after another dataset opens on its own view. Its camera is the default framing of the dimension it opens in; it kept the pose the previous dataset was left in, so a 28,000-cell 2D file opened after Suo (whose preset puts the camera inside its 3D cloud) filled the canvas with giant points. The point size follows its cell count, as before; a session or published preset still sets both afterwards.
- **Web app**: A file without a default UMAP opens on the default whose name ranks best (UMAP, then t-SNE, force layouts, PHATE, PaCMAP/TriMap, spatial, diffmap, any other name, then PCA and latent spaces), the higher dimension among equals, instead of always on its highest dimension. A segmentation export with 2-column `X_umap_scale*_latents` and a 3-column t-SNE or spatial array opened in 3D on the 3-column array; it opens on the UMAP in 2D, with the 3D default still a key press away. A file whose defaults are all UMAPs opens in its highest dimension as before.
- **Web app**: In Data mapping, an `obsm` array wider than three columns (`X_pca`, a latent space thousands of columns wide, a 4-column `coordinates`) is offered at its first two columns in 2D and listed after the embeddings; it was offered as a 3D embedding at its first three columns beside the file's real 2D UMAP.
- **Web app**: A dataset chosen while another was still opening stops the replaced opening at once (a catalog dataset's read, or a file's read and its reader worker). Its progress stayed on screen until the stopped read answered, the newer opening's outcome (*User data ready: …*) could be lost, and a sample or file chosen while the startup dataset loaded was announced under the replaced dataset's name.
- **Web app**: The `Full`-quality highlight pulse no longer holds the main thread on a slow renderer. A frame drawn only for the pulse now waits until what the frames before it held the main thread is measured (from the first of them to the first task after them), and then until the page has been free for as long (`rendering/ambient-frame-pacing.js`), so pulse frames take at most half of the main thread; where a frame costs less than a display interval nothing changes. A frame anything else asks for is drawn at once. While the measurement is owed the viewer's loop stays unarmed, so the measuring task always runs, also where the browser runs the next animation frame before queued tasks (Firefox, and WebKit at times). `tests/ambient-frame-pacing-contract.test.mjs` proves the rule on that event loop with 450 ms frames, in both task orders and under a stream of worker replies. Headless Chromium with 2 analysis workers: a 400-gene marker run with highlighted pages went from 7.0 s to 0.5 s (9.7 s to 2.0 s at 4× CPU throttling); on the CI runner the same run had not finished in 30 s.
- **Web app**: A category no cell holds has no centroid. The centroid pass returned one for every category in a field's list, at the origin with `n_points: 0`, and the viewer drew a point and built a label for each: a field listing 70,000 categories of which three hold cells built 70,000 label elements per view (140,000 after restoring a session with a kept view), stacked at (0, 0, 0), and the stats line read "Centroids: 70000". The pass now returns the categories with cells, in category order; its consumers find a centroid's category by identity, never by its place in the list (the visibility pass keyed categories by their text, so a `1` and a `"1"` shared one count). A consensus column rebuilt while it is shown is drawn again from its new codes, colours and centroids, so a category that gains its first cell gets its centroid at once. Colouring by the 70,000-category browser fixture took 0.74 / 1.26 / 1.21 s (Chromium / WebKit / Firefox) and now 0.31 / 0.69 / 0.59 s; restoring its session 2.2 / 3.0 / 3.3 s, now 0.65 / 1.3 / 1.0 s; at 6× CPU throttling in Chromium, 4.3 → 1.8 s and 13.0 → 3.5 s.
- **Web app**: Every view's copy of a field you made is the same field. Deleting or merging one of its categories, or rebuilding a consensus column, replaced the categories, codes and centroid tables of the edited view's copy and of the stored field only. A kept view then switched to an embedding not yet loaded failed with "Normalized user-defined centroids of layout "X_umap_2d" were not published", went on drawing the old categories, and an edit made in a kept view gave the live view the new categories under colours meant for the old ones. Each copy now takes the new definition at once and keeps its own colours and hidden categories, by category; a view that shows the field, or hides some of its categories, is drawn again. A field you made that is deleted, restored or confirmed deleted in one view is so in every view once that view is selected, as a dataset field already was, and Confirm releases its codes from every copy (from the Deleted fields list it only set flags and kept them).

- **Web app**: A data mapping applied while another dataset is opened belongs to the dataset it was asked of. Data mapping's **Apply** reads its neighbour graph, gene list or embedding over the network before it publishes anything; a dataset chosen meanwhile received the old dataset's graph, gene list and mapping label, its refusals and, on failure, the old dataset's state restored into it. The transaction now stops at the first await after another dataset is published (an `AbortError`), publishes and restores nothing, and the panel shows nothing of it.
- **Web app**: An outlier filter being computed stops when another dataset is published. Its read of the old latent space ran on (hundreds of MB at 2M cells), then announced "Outlier filter … ready" for the old dataset, and a latent space it could not use was marked unusable in the new one.
- **Web app**: An analysis releases what it left in the analysis workers however it ends (finished, cancelled, superseded or failed): a marker run's groups no longer stay in every worker or get replayed into replaced workers, a run cancelled during its setup leaves no groups in the workers that took them, and a UCell signature's per-cell ranks (about 18 bytes a cell per set gene) go with it.
- **Web app**: An analysis scoped to **Visible cells** reads the filters: a filter change marks its result **Stale · Run** and reads its groups again; it changed neither. A highlight change also marks stale the result of an analysis mode that is not shown (only the shown mode heard it), which catches up when shown.
- **Web app**: Gene search lists the gene whose name is exactly what was typed first, then names starting with it, then the rest (in the sidebar and in the gene-set editor), so `TNF` + Enter colours by TNF, not by whichever of `C1QTNF1` or `TNFAIP3` came first in the file.
- **Web app**: The notebook's `clear_highlights()` (Python) and `cellucid_clear_highlights()` (R) clear every highlight page, as documented; they cleared only the active one. The Clear button still clears the active page (`clearActivePageHighlights`).
- **Web app**: Main-thread work per interaction at scale. Marker genes' display options (direction, filters, ranking, top N) pick each group's top genes with a bounded heap instead of sorting every gene of every group: 100 groups of 30,000 genes, 716 → 36 ms per change. Differential expression resolves its two groups (a pass over every cell, 8 ms at 2M cells) once per change of the options that choose them, not several times on every option, highlight and page change. The highlight page tabs count a page again only when its enabled cell lists change (each highlight or page event, and every input event of a colour drag, recounted every page); a colour drag publishes once per frame. Mean z-score signatures reuse gene moments a marker run already computed and yield between genes. Marker genes no longer reread its groups on a highlight change when **Cells** is All cells.
- **Web app**: No GPU-pipeline drain on a frame path. In Firefox and WebKit, `isEnabled` and `getParameter` of the active texture unit, viewport, scissor box, colour mask, blend functions and clear colour each wait for every queued GPU command (6–8 ms in Firefox, 16–24 ms in WebKit after one 64 MB upload; bindings and pixel-store parameters answer at once in every engine). The backdrop look asked three of them per pane per frame, the GPU gather one per draw-order change, the smoke density, noise and edge uploads several per build. Every renderer now reads that state from a page-side shadow of the context (`rendering/gl-state-shadow.js`), kept by the context's own setters.
- **Web app**: The viewer redraws less. A smoke volume waiting for its noise textures no longer redraws the scene every frame until they land (the generation wakes the viewer once); reads, listener and handler registrations, `withBorrowedViewData` (every figure-export preview build), `stopInertia` and the navigation preferences no longer mark a frame, which redrew every pane of a grid; a live filter publication that moved no cell keeps the vector-field particles and their trails.
- **Web app**: A background spatial-tree build whose positions were replaced (another dataset, representation or a removed kept view) is dropped with its worker: it ran to the end first, so the new data drew without culling or LOD meanwhile, and kept the old positions alive. A replaced dataset also releases its smoke volume (its 3D textures and the copies of its colours and visibility) at once.
- **Web app**: Keyboard focus survives a rebuild. Toggling or removing a filter, or a highlight group, keeps the focus on that row's control (on the row now in its place after a removal, on the list after the last one); each row's remove button is named for its item ("Remove filter: cell_type: 2 of 3"). A long categorical legend keeps the rows still in view when it scrolls (it rebuilt every row, dropping the focus and an open rename). The phone sheet tabs take arrow keys, Home and End, and keep the focus on the tab chosen.
- **Web app**: The orbit anchor is shaded for the pane's real camera: its eye used another spherical convention, which at the horizon put the camera overhead, so its specular highlights and rim light faced nowhere.
- **Web app**: Smaller fixes: a failed group read in Marker genes or Differential expression says why on Run instead of reading again and posting its error on every redraw; the retired-notification log appends a row instead of rebuilding the list; the drop overlay names no format list that can fall behind the registry; analysis plot names and options are in sentence case.

- **Web app**: Messages name the sidebar's panels and blocks as their headings do ("Turn "Show overlay" off under Vector field", not "under Vector Field Overlay"; "in the Highlighting section", not "Highlighted Cells"): they take the names from one module (`app/ui/place-names.js`), and a contract test holds each name to its heading in `index.html` and fails on a message that names a place in words of its own.
- **Web app**: Analysis plots lay their labels out for the size they are drawn at, and again when resized, so none lands on another at any panel width. The correlation's per-page R² boxes were stacked 15% of the plot apart, less than a box is tall: they stand in a column, else a row, else one compact box beside the legend (which gives way in a narrow preview), and each colour category is listed once across pages, in one colour. The dot plot's colour bar and size legend stand in a column right of the plot, placed in pixels: the size legend's title ran into the colour bar in the sidebar preview and its dots into each other in the expanded view. Its gene names shrink and turn to fit their columns, and a panel too narrow to name every gene scrolls the preview across. A legend above a plot gets the rows it needs, or is left out of a preview too short for it.
- **Web app, Python, R**: An embedding under the R world's conventional UMAP name is the default UMAP, as Scanpy's `X_umap` is: after `X_umap_1d`/`X_umap_2d`/`X_umap_3d`, the first 1- to 3-column one of `X_umap` (also a Seurat object's `umap` reduction), `X_UMAP` (a SingleCellExperiment's `UMAP` reducedDim, as scater's `runUMAP()` names it), `X_wnn.umap` and `X_ref.umap` (a Seurat object's weighted-nearest-neighbour and Azimuth reference UMAPs) fills the slot of its width. The rule is pinned in the shared data-mapping conformance fixture (version 6, SHA `cfc84c8c…`, in all three repositories; web `defaultEmbeddingSlots`/`UMAP_KEYS`, Python `cellucid.data_mapping.default_embedding_slots`/`UMAP_KEYS`, R `.shown_embedding_keys`/`.umap_keys`).
- **Web app**: A sample whose file cannot be read while the list is made (a busy host, a file taken down) no longer empties the dataset picker with "The dataset list could not be loaded": the open dataset stays shown, the notice names the sample that failed with Try again, and a dataset that opens does not clear that notice.
- **Web app, development**: The two design-token validators (`scripts/validate-tokens.js`, `scripts/validate-token-types.js`) run with the contract suite; `types/design-tokens.d.ts` had drifted by sixteen names. Dead stylesheet rules for removed controls, two unused sidebar-bound tokens that disagreed with the real clamp, a duplicate `.gene-chip` rule and raw opacities are gone.

- **Web app**: On a phone the open full-screen menu no longer shows a band of the plot above it: the strip its ✕ sits in is painted as part of the menu (it was transparent). The ☰/✕ toggle no longer stays grey and pressed after a tap on a touch screen; its hover look is kept for a pointer that hovers.

- **Web app**: After a deploy, a browser cache no longer hands the page modules of two builds. A page could run a module its browser had cached from the build before a deploy beside one fetched from the build after it (GitHub Pages caches every file for ten minutes, a static server lets the browser guess a freshness of its own), and the app stopped on "The requested module '../../../data/data-mapping.js' does not provide an export named 'UMAP_KEYS'" with "Starting Cellucid…" shown forever. Every script is now requested at a URL that names its bytes (`<path>?v=<first 12 hex digits of its SHA-256>`): `index.html` (and the sign-in callback page) carries an import map of every module under `assets/`, the `src` of each script tag says the same, and the Content-Security-Policy allows the map by hash. `scripts/generate-web-asset-inventory.mjs` writes all of it, and the JSON-LD hash, beside the inventory; nothing needs a build step to run. A file that did not change keeps its URL and stays cached. A startup module that fails to parse, link or run now shows the "Cellucid could not load" card with **Reload** instead of waiting forever, and when the server already serves a newer build than the page's the card says "Cellucid has been updated". A page left open across a deploy that then loads a part on demand (the benchmark panel, figure export renderers, Plotly) and fails says so in a notification with **Reload**. Workers, which get no import map in any browser, are covered by the next entry.
- **Web app**: A deployed page starts every worker from one bundle named by its bytes, so a worker can no longer pair its own modules from two builds. `scripts/generate-web-asset-inventory.mjs` now builds each worker (AnnData reader, decode, gunzip, spatial tree, point raster, figure export, analysis, benchmark, and the two modules decode workers import by URL) with esbuild 0.28.2, pinned as a dev dependency, into `assets/workers/<name>.<first 12 hex digits of its SHA-256>.js` with a source map to the sources, and writes `assets/workers/scripts.js` naming them; they are committed, listed in the inventory (the Python and R web caches download and verify them) and `npm test` refuses a stale one. On `localhost`, `127.0.0.1` and `[::1]` (development, the browser tests, the Python and R servers) workers still start from their sources, with nothing to build. A worker that fails to start on a page that outlived its build shows **Cellucid has been updated** with **Reload**.

- **Web app**: Opening a file no longer hangs in WebKit. The inflation of a compressed chunk on the reader thread (`formats/io/inflate.js`) released its stream reader's lock, and WebKit (314030@main, also landed on a Safari release branch) runs `ReadableStreamDefaultReader.releaseLock()`'s JavaScript while it holds a lock its garbage collector takes when it visits the reader: a collection that started inside the release never finished, and every later collection in the tab, the page's included, waited behind it. The file stayed "Opening …" for good and closing the page hung (sampled native stacks: the reader thread parked in that lock inside `releaseLock`, the main thread in `JITWorklist::suspendAllThreads`); under CPU load it hung within 1 to 5 opens of the 70,000-category fixture, and the browser suite's WebKit shard failed on it. No browser code releases a stream reader any more (`tests/browser-stream-reader-release-contract.test.mjs`): a read-to-the-end or cancelled reader is dropped, and 100 opens in a row under the same load all opened. Closing the app also stops the dataset change still opening and the file or connection it reads (`retire()` of the dataset-change coordinator) instead of waiting for them to answer, and an opening that makes no progress for 60 s, the startup stall interval, says which step it waits on and the reads still open with what has arrived of each, in its notification and in the console (`[cellucid:opening] Opening a.h5ad stalled while reading the file's structure: nothing has changed for 1 min 0 s, while reading the structure of a.h5ad (0 B of 512 KB arrived).`), without stopping it; while a startup step's own deadline runs, that deadline reports the stall instead, once.
- **Web app**: The genes an analysis computes in its workers (marker genes, differential expression, the gene-moments pass of a module score, the dot plot, UCell) no longer pass through the page. Each gene went from the AnnData reader worker to the page and from the page to an analysis worker, and WebKit counts every array a worker sends the page as the page's own allocation: a run over every gene drove the page's garbage collector, whole-heap collections included, while the analysis ran, each collection holding the page's thread. A run now opens a gene route, one `MessageChannel` per analysis worker, and the reader delivers each task's genes straight to the worker that computes them (`data/gene-delivery.js`, `app/analysis/data/gene-route.js`); the page sends which genes each task takes, hears their bytes and takes in each task's statistics, packed. A task takes several genes (sixteen tasks per worker over the run), one read serves a task on each worker, and a worker the pool replaces during a run has its genes read again. Measured on the 20,000-cell, 200-gene markers run: 13 MB of gene values in and 13 MB out of the page per run became none, the page's allocation per run halved, and 40 reader requests became 13 (8 workers).

## [0.9.4]

Version 0.9.4 is a patch release on 0.9.3. The export format stays 3, so
exports written by 0.9.3 open unchanged; `prepare(obs_categorical_dtype="uint32")`
is the one addition to what an export may hold. The web app, the Python
package and the R package now share one version and are released together.
The main changes are:

- **Six new backgrounds** (Atelier, Nebula, Abyss light and dark, Aurora light
  and dark), prepared in the background and sharpened while the GPU keeps up.
- **No fabricated limits**: any number of kept views, selections of any size,
  direct H5AD and Zarr files of any size, categorical fields with more than
  65,535 categories, and session files bounded only by their format.
- **Analysis on the workers**: summaries, tests, histograms, violins and
  density plans leave the main thread; gene downloads go through one scheduler
  with a payload cache; marker-gene clustering is bounded by memory, not a
  count.
- **Jupyter selections** of any size reach Python as binary indices.
- **Sessions** restore the vector-field controls after the embeddings they
  depend on, and keep every view's embedding.

### Added

- **Web app, version**: the web app declares the release version it shares with the Python and R packages: `version:` in its `CITATION.cff`, a `cellucid-version` meta element and `softwareVersion` in the page's structured data, and the footer prints it before the build id. The three are released together under one version.
- **Web app**: Six new backgrounds, in light and dark pairs: **Atelier (light)** and **Nebula (dark)**, **Abyss (light)** and **Abyss (dark)**, **Aurora (light)** and **Aurora (dark)**, listed under the Light and Dark groups of the **Background** menu; **Studio (light)** stays the default. Each is a painted stage that turns with the camera, with a nearer layer that slides against it for depth: Atelier's paper studio, where the cells cast a soft contact shadow that follows filters (and glows through the paper seen from below); Nebula's deep field, tinted by the colours of the field the cells are coloured by; Abyss's sunlit shallows and deep ocean with glowing marine snow; Aurora's dawn and night skies. A look is prepared in the background in idle time the first time it is chosen (its sky coarse first, then sharpened level by level while the GPU is measured to keep up, so it appears within a fraction of the full work on every machine and reaches 1024² per face on a capable one), costs nothing while the scene is still, keeps its large areas outside the lightness range of the category palettes, and appears without its fade under the system's reduce-motion setting. Only the look on screen keeps GPU memory (at most about 27 MB). Sessions and dataset presets can name the looks, and figure export redraws the look under the points in PNG and SVG (as an image) when the figure's background is Match viewer. (backgrounds)
- **Web app**: Keep any number of views: the four-view ceiling is gone. GPU memory decides instead: when an allocation is refused, the colour and visibility stores of kept views that are not on screen are released first and rebuilt, unchanged, when shown again. The benchmark harness's view-count axis takes any number too.
- **Web app**: Selections of any size: the annotation, lasso, KNN and proximity tools, category highlights and Jupyter commands hold selected cells as one bit a cell instead of a JavaScript `Set`, which Chromium refuses above 16,777,216 entries, so selecting a large cluster in a 30M+ cell dataset failed. Combining selections (union, intersection, difference) works a word at a time.
- **Web app**: Categorical fields with more than 65,535 categories (barcodes, free text) load: they use 32-bit codes (export format `uint32`, missing marker 4,294,967,295) in prepared exports, `cellucid serve`, Jupyter, direct H5AD/Zarr and saved sessions. A direct H5AD/Zarr or served AnnData with such a column used to be refused as a whole. Existing exports are unchanged. (limits)
- **Web app**: Direct H5AD files have no size ceiling: files above 512 MiB used to be refused before they were read. The file is read from disk as each view needs it instead of being copied into memory, in a background worker, so opening it no longer blocks the interface (a 250 MB, 1-million-cell CSR file: about 0.8 s of main-thread loader work and 1.1 GB peak before, about 0 s and 0.8 GB now; it opens in about 0.25 s instead of 0.7 s). A row-major (CSR) matrix is reorganised by gene in the background right after opening, with progress shown, and a gene asked for meanwhile is answered from the same pass; when the reorganised matrix does not fit the memory budget it is kept on disk (the browser's private file storage). The fixed working-set ceilings on sparse matrices, dense matrices and observation columns are removed. (limits)
- **Web app**: Zarr directories and Zarr ZIP archives load through the same background reader, with several chunks decoded at once, and without the fixed per-array and per-archive size ceilings; every chunk is still checked against the exact size its shape and dtype declare. (limits)
- **Web app**: Memory budgets for kept layouts, inactive fields, the AnnData cache, decoded blocks and the gene index scale with the memory the browser reports (Chromium) instead of fixed values, and they decide only what is kept for reuse: no dataset is refused for them any more. Browsers that report nothing keep the previous values. (limits)
- **Web app**: Latent-space outlier quantiles of a direct H5AD/Zarr file work for a latent space of any size: the 512 MiB latent-space budget is removed. The latent space is streamed from the file in row slabs in the reader worker, with the same bit-identical rule, and kept in memory only while it fits the reader's budget. (limits)
- **Web app**: Categorical fields made from highlight pages work at any cell count (they were refused above 16,777,216 cells, a JavaScript `Set` limit), an intersection field combines up to 32 pages (was 12), and the limits of 20 user-defined fields and of 256 characters for a field name or category label are removed. (limits)
- **Web app**: Session files have no 16 MiB manifest or 512 MiB chunk limit; the only limit is the format's 4 GiB − 1 per frame (its length is a u32), and every inflate is still bounded by the chunk's declared size. (limits)
- **Web app**: Prepared exports have no 512 MiB size ceiling per file: a 3D dataset above 44.7 million cells was refused. Every payload is checked against the exact size its manifest declares, and one the browser cannot allocate is reported as too large, naming the bytes. Compressed payloads above 4 GiB load too (their gzip trailer is compared modulo 2^32).
- **Web app**: A sample dataset's published starting view has no size limit (it was refused above 32 KiB, or 64 KiB decompressed): the catalog now declares the state file's exact length in `state_bytes` beside `state_sha256` (`exports-v3/datasets.json`, and `cellucid serve`'s dataset listing), and the viewer reads exactly that many bytes and checks the digest. A catalog entry advertising a state must carry all three fields. (limits)
- **Web app**: Camera Path: no limit of 100 keyframes, 60 seconds per transition or 40 characters per label; a transition still lasts at least 0.1 s. (limits)
- **Web app**: Highlight selection undo reaches every earlier step: the oldest steps go only when every tool's recorded selections together outgrow a share of the device's memory (was 5 steps per tool). (limits)
- **Web app**: The cell info card lists every loaded categorical field (was 4), adding columns when they do not fit the window's height. (limits)
- **Web app**: The gene search list lists every match (built 100 at a time as you scroll or arrow; Home and End jump to the ends) instead of the first 100 with "…and N more". (limits)
- **Web app**: Notifications retired early to keep five cards on screen are kept: a "Show N earlier messages" control under the stack lists every one of them. (limits)
- **Web app**: Community annotation: a Pull has no limit of 10,000 user files or 64 MB in total. Changed files are fetched and cached about 32 MiB at a time, and the merged view is compiled by streaming the cache one file at a time; an interrupted Pull keeps what it cached. CAP lookup results are kept whole (results with more than 100 synonyms, 200 markers or 512-character text were dropped), every fetched result is shown (was 20 of 25), and synonyms, markers and author handles are shown in full. (limits)
- **Web app**: Analysis: the gene list lists every matching gene (built a page at a time while scrolling or arrowing) instead of the first 50; analysis windows and saved plot options have no count limit. (analysis, limits)
- **Web app**: Marker genes: clustering is no longer limited to 500 genes or groups. It is admitted by the memory it needs (about 100 bytes per pair of genes: about 100 MB for 1,000 genes, 1.6 GB for 4,000) against the analysis memory budget. A request past the budget is refused with that size and what to change, instead of "maximum is 500". (analysis, limits)
- **Web app**: Jupyter: selections of any size reach Python (a million cells arrive in about 120–200 ms). A confirmed selection is posted to the notebook's server as binary cell indices (`POST /_cellucid/selection`, little-endian `uint32`, ascending) instead of a JSON event sent with `keepalive`, which the browser refused above 64 KiB (about 8,000 cells). Validating and encoding one million cells takes a few milliseconds instead of 150–300 ms, and no step builds a per-cell `Set`, so selections above 16,777,216 cells are no longer refused. No event is sent with `keepalive` any more. Works with the matching `cellucid` Python version, released next. (DX-17)
- **Web app**: Gene downloads (analysis and the sidebar) go through one download scheduler: a gene the user chose starts at once and no new download of a long run (DE, markers, heatmap) or of a prefetch starts until it has arrived; the browser gets the same order through the fetch `priority` hint, and a prefetched gene the user then chooses is raised, not downloaded twice. Observation-field loads count as foreground too. A finished or aborted run cancels its queued and in-flight downloads, and replacing the dataset ends every download of the old one, so a stalled read can never hold back the next dataset. (analysis, network)
- **Web app**: A dataset load keeps its downloaded gene payloads (compressed, up to 64 MiB) and decodes an evicted gene again from them: a second DE or marker run over the same genes makes no requests (pancreas over 4G: 3,753 revalidation round trips before). Long runs download the payloads of the genes after their decode window ahead, as many at once as the connection needs (the window starts at 8 and doubles while throughput grows, never more than half the cache), so a run is bound by bandwidth instead of round trips. (analysis, network)
- **Web app**: A gene the user waits for shows its download progress in its loading notification. (analysis, network)
- **Web app**: The documentation moved to https://docs.cellucid.com (the old Read the Docs addresses redirect there); the footer, Cite, info-tip and welcome links point at it. Search and sharing metadata: a description that names UMAP, 3D, Scanpy/AnnData and Seurat; `og:image:alt`, `twitter:image:alt` and `og:site_name`; corrected JSON-LD (`WebApplication`, operating system "Web browser", licence, documentation, `sameAs` for GitHub, PyPI and the R package), still allowed by its CSP hash; a crawlable "About Cellucid" text with links to the docs, GitHub and citing; a `sitemap.xml` whose `<lastmod>` is the build date, written by `scripts/generate-web-asset-inventory.mjs`; a custom `404.html`; an `llms.txt`; the manifest name "Cellucid"; and a preload of the interface font.
- **Web app**: Data mapping, Embeddings: one compact row per array -- checkbox, name and dimension; an embedding that is on adds a second line with **Default** and **guessed**, one that is off stays one line. A long name wraps to a second line instead of losing its columns to an ellipsis. The dimension of an AnnData array is a chip with a chevron (`3D ⌄`) that opens the row's editor, closed until opened: **Dimension**, one labelled column menu per axis on one line, and, for an array wider than three columns, **Add columns** (it replaces the "+ Columns of …" pill under every wide array). Nothing in the panel overflows the sidebar on a desktop or a phone. (data mapping)
- **Web app**: Data mapping: every `obsm` array of two or more columns can be shown in any dimension its columns allow, at any of its columns -- a 2-column UMAP in 1D or 2D, a 3-column array in 1D, 2D or 3D, as a wide array already could. The guessed default is unchanged; views, sessions, highlights and figure export take the new selections like any other (`UMAP 1`, `Spatial 1,3`). (data mapping)
- **Python**: Jupyter: `POST /_cellucid/selection` receives selections as binary cell indices (little-endian `uint32`, strictly ascending, with no size limit; a body the server cannot hold is refused with HTTP 413), authenticating before it reads the body; `@viewer.on_selection` still receives `cells` as a `list[int]`, now always ascending. `POST /_cellucid/events` refuses `selection` and `session_bundle` records (HTTP 400): only their own routes produce them. Selections were previously capped by the event endpoint's 1 MiB limit. (DX-17)
- **Python**: `prepare(obs_categorical_dtype="uint32")` writes 32-bit categorical codes (`obs/{index}.codes.u32`, missing marker 4,294,967,295) for fields with more than 65,535 categories, inside export format 3; existing exports are unchanged and an older web app refuses such an export by name. The width error now names the width that fits ("Pass obs_categorical_dtype='uint32'."). `serve()` and `show_anndata()` serve a string or boolean obs column with more than 65,535 distinct values with uint32 codes instead of refusing the whole AnnData, and session bundles carry uint32 user-defined codes (encodings 4 and 5). (limits)
- **Python**: `cellucid serve` lists a published default session with its exact length (`state_bytes`) beside `state_manifest` and `state_sha256`, and reads a bundle of any size (the 32 KiB limit is removed), hashing it as it streams. `debug_connection()` accepts the three fields together. (limits)
- **Python**: Documentation: served from https://docs.cellucid.com, the Read the Docs canonical domain; every link in the package, README, CITATION.cff and docs uses it. Every page declares `<link rel="canonical">` on that host (`html_baseurl`), `sphinx-sitemap` (pinned in the `docs` extra) writes https://docs.cellucid.com/en/latest/sitemap.xml, and a `robots.txt` at the root of the default version advertises it.

### Changed

- **Web app**: When several views' spatial trees land at once (four views of four representations), each view switches to its tree in its own frame instead of all in one long frame.
- **Web app**: Slider, legend and camera-path changes are drawn in the frame they happen in, not the next one: they used to wait for an animation frame that ran after the viewer's, which added one frame of latency to every slider step (about 50–90 ms at 10 million cells in four views).
- **Web app**: Starting a drag or a wheel zoom in another view no longer waits for the sidebar: the camera follows the view at once, and the sidebar switches to it (legend, counts, controls) when the gesture ends (on release, or once the wheel has rested for 150 ms). A wheel turn that cannot zoom further, or a wheel in free-fly mode, no longer redraws the scene.
- **Web app**: The full-window selection overlay (lasso, proximity and neighbour rings, the projectile charge ring) is hidden while nothing is drawn on it, so the browser no longer composites an empty window-sized layer over the plot every frame.
- **Web app**: Category counts in the legend are recounted only while the legend can be seen: with its panel collapsed or the sidebar hidden, a filter change no longer recounts every category; the counts update when the legend is shown.
- **Web app**: A camera drag captures the mouse until it is released, so dragging across the sidebar or an analysis plot no longer makes the page work out hover styles and plot tooltips under the pointer at every move.
- **Web app**: Shader programs are compiled in batches: the point programs compile 9 shaders instead of 14 (shared vertex shaders) and the velocity overlay 8 instead of 14, with every link issued before any status is read, so the driver overlaps the work and the page waits once per batch instead of once per compile and link. Turning velocity on and starting the app pause less.
- **Web app**: Multiview: a view whose camera and contents did not change keeps its pixels instead of being redrawn, so dragging, zooming or flinging the camera of one view of four redraws only that view, and the settle frame after a change in a single view draws nothing.
- **Web app**: Zooming through levels of detail no longer rebuilds the figure-export preview at every step; it rebuilds once the camera settles.
- **Web app**: An enabled velocity overlay with nothing to draw (no field for the layout, nothing visible, particles off when zoomed out, smoke mode) no longer keeps the scene redrawing every frame.
- **Web app**: Filtering no longer blocks on the GPU: writing alpha, colours and draw orders into storage the GPU already holds no longer waits for the previous frame to finish (a filter tick at 10M cells in four views spent about 40 ms waiting).
- **Web app**: Framed (cropped) figure exports are drawn at their full resolution: a frame's points are rasterized at the pixels they land on instead of being enlarged from the whole view, so a framed PNG, or the point layer of a Hybrid SVG, is as sharp as an unframed one.
- **Web app**: When a Remote server connection drops, the viewer says so at once, even while another dataset you chose is still loading; the notice is completed with "Reconnect to continue" once recovery has run.
- **Web app**: Filtering a kept view uploads only the rows of its visibility texture that changed, reuses its own transparency array instead of allocating a new one per update, and skips its highlight, edge and overlay work when nothing moved.
- **Web app**: Recolouring a kept view rewrites its colour store in place instead of allocating a new store and vertex array each time (about 70 MB transient a view at 10M cells).
- **Web app**: A kept view's colour store is RGBA8 on a 4-byte stride, like the live one: Metal (Chrome and Safari on macOS) needs that alignment, and a 3-byte RGB store made the browser keep and re-convert a hidden 4-byte copy.
- **Web app**: Figure export of the Studio background no longer has a ceiling on the dots it considers: a flat embedding's dot surface seen at a grazing angle is enumerated only where it shows, and the exported marks are unchanged.
- **Web app**: An idle page no longer wakes 60–120 times a second: the frame loop stops once the scene has settled and starts again when something asks for a frame (an interaction, a data or setting change, a resize, a held navigation key).
- **Web app**: The KNN selection tool builds its neighbour graph as compact index ranges instead of a JavaScript map per cell, so it works on graphs with more than 16.7 million cells and builds much faster.
- **Web app**: Spatial trees (level of detail, frustum culling, hover) are built in the background worker at every dataset size; above 512 MiB of positions (44.7 million cells in 3D) they were built on the main thread, freezing the page for about 15 s at 50 million cells. Positions go to the worker in 64 MiB slices, so the page never holds a second full copy.
- **Web app**: Analysis keeps the interface responsive on large pages: the Detailed and Gene signature summaries and tests (from 200,000 values), the Quick Insights median and quartiles and the signature-score median (from 500,000 finite values), histogram bins (from 200,000 values) and density-plot plans are computed on the analysis workers, with the same functions and the same numbers as before. Categorical counts are taken once while a page is built and shared by every consumer; gathering a page and aligning two variables for a scatter plot yield to the interface in 12 ms slices. (analysis)
- **Web app**: Gene search, in the sidebar and in analysis: a gene kept under the pointer or highlighted with the arrow keys for 150 ms is loaded in the background, so choosing it usually draws at once; moving on cancels the load and choosing the gene joins it. (analysis)
- **Web app**: Violin plots of large pages no longer freeze the interface: the box statistics Plotly's violin calc needs are computed on the analysis workers with the density, and the vendored Plotly takes that calcdata entry instead of building a point record per value. A ten-million-cell violin held the main thread for 1.5 s and 800 MB; it now holds it for under 50 ms and 90 MB, drawing the identical violin. Violin outliers, like box-plot outliers, draw at most 16 copies of an identical value. (analysis)
- **Web app**: Analysis worker dispatches copy at most 96 MiB per main-thread task; a ten-million-cell density plot used to copy about 1 GB in one 114 ms task. A density over more than 4,096 distinct values sends the page's own values with unit weights instead of a float64 copy and an array of ones (a quarter of the bytes). Quick Insights' category composition yields to the interface like its continuous summary. (analysis)
- **Web app**: Latent-space outlier quantiles are 4 to 6 times faster (10 million cells by 10 columns: 4.6 s to 0.8 s), bit-identical to `cellucid-python` and `cellucid-r`: the latent space is read row by row and large categories are ranked by radix ordering. (analysis)
- **Web app**: Removed the analysis bulk loaders that nothing called (`DataLayer.fetchAnalysisData`, `fetchComprehensiveAnalysisData`, `loadAnalysisBulkData`, `loadLatentEmbeddings`); every analysis gene load goes through the field loader and its cache. (analysis)
- **Web app**: The outlier filter control stays hidden for a field whose every outlier quantile is missing (no category reaches `centroid_min_points`), as for a field without quantiles. (limits)
- **Python**: `prepare()` exports a categorical field where no category reaches `centroid_min_points` with quantized outlier quantiles, publishing every quantile as missing, instead of refusing the export. The mapped latent-space route of `serve()`/`show_anndata()` has no 512 MiB latent-space budget any more and uses `obsm` without caching a copy. Session bundles accept field names and category labels of any length (the 256-character limit is removed). The release build's sdist limit is PyPI's 100 MB upload limit instead of 1 MiB, and the 4,096-member limit is removed. (limits)
- **Python**: `.cellucid-session` reading has no 16 MiB manifest or 512 MiB chunk limit, matching the web app: a frame holds up to 4 GiB − 1 and a length past the end of the file is refused. (limits)
- **Python**: Jupyter: `get_session_bundle()` uploads have no 512 MB limit; the upload streams to disk and is bounded by the free space of the temporary directory (HTTP 413 beyond it). `wait_for_event()` and `set_data_mapping()` are handed their event as it is recorded, so a burst of other events (thousands of hovers) can no longer push it out of the 512-event log they used to scan. (limits)

### Fixed

- **Web app**: Sessions: the vector-field controls (overlay on or off, the field, its colour map and sliders) are restored after the embeddings, the kept views and the active view, because the field select lists the fields bound to the embedding the active view shows. A session saved while that view showed an embedding without vector fields was refused by a page whose view showed one (`Session UI select option "" is not available in current control "velocity-field"`), and the other way round the saved field was overwritten when the embedding came back. The block is marked `data-state-serializer-after-embedding` in `index.html`; any control whose choices follow the embedding can use the same marker. (K3)
- **Web app**: The page preconnects to the sample-catalog origin only when it will fetch from it: a page told another catalog (`?exportsBaseUrl=`, as `cellucid serve`, Jupyter and the browser tests do) no longer opens a connection to `theislab.github.io` it never uses. WebKit reported that unused connection's failed TLS handshake as a console error when the network was poor. (`assets/js/app/ui/core/preconnect-init.js` replaces the static `<link rel="preconnect">`.)
- **Web app**: Benchmark harness: a configuration waits until every presented view draws from its final geometry (its spatial tree, with levels of detail when LOD is on, built and adopted by its store) before a measurement window starts. On a slow host the tree built off the main thread landed inside the window's first frames, so the window measured the landing's uploads and error reads and reported sync stalls that no frame of the configuration makes. The viewer exposes this as `isPresentedGeometrySettled()`.
- **Web app**: Direct H5AD/Zarr: the embedding payload is handed to the viewer instead of being shared with the reader's cache, so evicting a layout now frees its memory and a reload re-reads the file rather than renormalizing an array it already rewrote. Embeddings with a positive range below 1e-8 now get the export format's scale, the same one `cellucid prepare` and `cellucid serve` use, instead of a floor of 1. (DX-14)
- **Web app**: Firefox draws datasets above 30 million points. Every point draw (the point passes, the GPU gather, highlights and figure export) is now split into calls of at most 16,777,216 points; Firefox refused a single draw of more than 30,000,000 and showed nothing. The image is unchanged.
- **Web app**: Per-point textures (visibility, draw-order permutation, the GPU gather's source) are stacked in layers instead of one 2D texture, so no point count is refused for texture size. On GPUs that report a 8,192 texture limit (Firefox on macOS) a 3D dataset above 22.4 million points could not build its draw order, and a dataset above 67.1 million points (268 million elsewhere) was refused at load; the only device limit left is the texture layer count, at least 2^34 points.
- **Web app**: Connectivity edges load and draw for any graph the browser has memory for: a k = 15 neighbour graph was refused from about 12 million edges (1.1 million cells) and aborted the whole dataset. Only the edges a panel shows reach the GPU (the slider's K, 28 bytes each), so GPU memory no longer grows with the graph, and the edge shader reads no texture. The shown edges are the first K visible edges of one fixed random order computed on the fly (a keyed bijection), not a shuffled copy of the graph, so a different, equally random sample of edges is shown than before. Counting all visible edges for the slider no longer blocks a filter change. A graph the browser cannot allocate leaves the dataset loaded without edges.
- **Web app**: Switching back to a representation stays free whenever memory allows: geometry and colour stores kept for switching back are released only when the GPU refuses an allocation (then everything kept is released and the allocation runs again), never ahead of one. A spatial tree whose build failed for memory is built again on the next request once memory has been released, instead of staying off until level of detail or culling was toggled.
- **Web app**: Figure export draws its point passes on WebGL contexts kept for the page instead of creating and losing a context for every pass: exports start faster, and Linux WebKit (WebKitGTK/WPE) no longer crashes the page during an export (it aborted on a GL call made without a current context while a context was torn down).
- **Web app**: Figure export works on synthetic benchmark clouds (it refused them: "normalization for view live is incomplete").
- **Web app**: PNG figure export draws each plot's point raster at exactly the device pixels the plot covers: at 300 dpi a plot's edges fell between device pixels (for example at x = 256.25), so the raster was resampled with a quarter-pixel blur and its bytes depended on the browser's resampling. The raster is now rendered at the snapped size and copied one to one. (figure export)
- **Web app**: The Python documentation names `POST /_cellucid/selection` wherever it describes how the viewer reports back to Python; it had described `/_cellucid/events` alone, which carries every event except selections.
- **Web app**: Analysis loops that hand the event loop back (statistics, scoring, clustering) now yield through the app's one shared macrotask channel. Each yield used to create its own MessageChannel held only by its handler, which WebKit can garbage-collect with the message still pending, leaving the loop waiting forever. (analysis)
- **Web app**: Analysis reads a categorical field at the width it is stored in. A field with more than 65,535 categories (stored as uint32 codes) was copied into 16-bit codes, so every category above 65,535 was counted as another one, in Quick Insights, Detailed counts and marker discovery, without an error; marker discovery also refused more than 32,767 groups. Codes are no longer copied to another width, and a uint8 field is no longer widened for analysis. (analysis, uint32 categories)
- **Web app**: Sessions: restoring a session with several views keeps the live view on the embedding it was turned to. Restoring the views marked the live view as following its dimension's default, so when the session's mapping differed from the one the file opened with (for example a 3-column array saved in 2D), the live view was moved to the default. (data mapping)
- **Web app**: Data mapping: choosing something in the panel (an embedding, a column, Genes, Neighbour graph, Latent space, Vector fields) no longer scrolls the sidebar or drops the focus. The checklist was rebuilt on every change, and Safari's scroll anchoring then scrolled the sidebar by about a screen; rows are now updated in place, and unchanged text, options and lists are not rewritten. Pressing **Default** keeps the keyboard focus on it. (data mapping)
- **Web app**: The antialiasing tooltip and documentation now quote current measurements: in a dense 2D view antialiasing costs about twice the GPU time of a frame (it said about ten times).

## [0.9.3]

Version 0.9.3 is a patch release on 0.9.2. It introduces export format 3, so
exports written by 0.9.2 or earlier must be prepared again: a viewer reads
exactly one format. The main changes are:

- **Export format 3**: any number of embeddings per dimension, each with a
  label and its source columns; several gene-name sets; vector fields bound to
  an embedding id; centroids keyed by embedding.
- **Data mapping**: choose, per dataset, which embeddings are on and which one
  each dimension opens on, the neighbour graph, the gene names and expression
  matrix, the latent space for outlier filters, and the vector-field bindings,
  in the browser, from `cellucid serve` and from a notebook.
- **Several representations at once**: every view shows one of the enabled
  embeddings of any dimension, with its own camera and navigation mode.
- **Rendering**: points are stored in the draw order of their spatial tree,
  built in a worker; levels of detail draw a prefix, culling draws tree-node
  ranges, and colours are gathered on the GPU.

### Added

- **Web app, data loading**: drop one `.h5ad` file, one `.zarr.zip` archive or
  one prepared export folder onto the viewer to open it. An outline and a
  *Drop to open* card cover the viewer while you drag; the drop is read like
  the matching button and replaces any load in progress. Anything else is
  refused with a message saying what can be dropped, and a file dropped
  anywhere on the page, even while Cellucid is starting, never makes the
  browser leave it.
- **Web app, cell info**: **Show cell info** (Coloring & filtering, on by
  default) shows a card for the cell under the pointer: its value on the
  field the view is coloured by, where that sits on the colour scale, and for
  a category how many of its cells are visible in that view (what selecting it
  would pick); **Show recent fields** adds the cell's labels in other loaded
  fields. It works in every view of a grid; press and hold on a touch screen.
  A lookup never loads a field or redraws the plot (1M cells: 1.4 ms mean).
- **Web app, citing and sharing**: a **Cite** dialog (APA and BibTeX, read
  from `CITATION.cff`, with a placeholder note until the preprint has a DOI),
  and a **Share** dialog for datasets that can be opened from a link: the
  link, an "Open in Cellucid" README badge and an embed snippet.
- **Docs**: new pages on cell info, share links, badges and embeds, how to
  visualize single-cell data interactively, and how Cellucid relates to other
  viewers.
- **Citation**: a placeholder citation for the forthcoming Cellucid bioRxiv
  preprint, the same in every repository. It is written once, as the
  `preferred-citation` of `theislab/cellucid`'s `CITATION.cff`, and copied from
  there to this repository's `CITATION.cff`, the README, and the new
  [Citing Cellucid](https://docs.cellucid.com/en/latest/citing.html) page;
  the test suite fails if a copy drifts. Until the preprint has a DOI, every copy
  says it is a placeholder.
- **Web app**: cells can be highlighted on a touch screen. Every selection tool
  started only from an `Alt` press, which a phone or tablet does not have, so touch
  users could not select at all. A **Select mode** switch in the Highlighting panel
  now stands in for a held `Alt`: while it is on, a one-finger tap or drag runs the
  active tool through the same path an `Alt` press takes, and its Intersect / Add /
  Subtract buttons stand in for `Shift` and `Ctrl`. Two fingers always pinch to zoom
  and, newly, drag to pan; a second finger landing mid-selection abandons that
  gesture rather than committing it. `Esc` with nothing in progress leaves the mode.
  Mouse and keyboard behaviour is unchanged while it is off.

- **Web app**: a **selection bar** at the bottom of the plot carries Select mode's
  controls, so a selection can be made, combined, confirmed, cancelled and cleared
  without opening the menu, which covers the whole screen on a phone. It shows the
  same switch and Intersect / Add / Subtract choice as the panel, the active
  tool's cell count with its Undo, Cancel and Confirm while a selection is in
  progress, the active page's highlighted count with Clear, and a one-line note
  when a tap on empty plot changed nothing. Every control on it is the panel's
  control shown again, so the two always agree. It is always shown on touch
  screens (just the switch while the mode is off), shown to a mouse while the
  mode is on, centred in the part of the plot the sidebar leaves uncovered, and
  notifications stack above it.
- **Web app, backgrounds**: two new viewer backgrounds, **Studio (light)** and
  **Studio (dark)** (neutral grey and charcoal): a 3D view stands on a floor
  of constant-width hairlines that fades out before its edge, the same ruled
  floor from every angle above it (below it, the stage is plain); 1D, 2D and
  Planar views sit in front of dotted paper. Which one a view shows follows
  its dimension and navigation mode, never the camera angle, so rotating a 3D
  view never swaps it to dots. The pattern rescales with zoom without popping, and
  figure export reproduces it in SVG and PNG. **Studio (light) is the new
  default** (the sample datasets open on it too), and the Background menu is
  grouped into Light and Dark.
- **Web app, shader quality**: a fourth quality, **Sculpted (soft 3D
  light)**. Every cell is a soft-lit, glossy ball — a wrapped key light, a
  sky-and-ground ambient, a broad gloss instead of a hard white highlight, and
  darkening toward the edge that keeps touching cells apart — computed in
  linear light so the centre of each ball is exactly its legend colour. With
  antialiasing on, its edges are smooth in 3D views too, without the
  background showing through nearer cells; highlighted cells get the same 3D
  ring as `Full`, held still. It costs about what `Full` does (3–5% more in
  3D views on an Apple M1 Pro, the same in 2D). `Full` stays the default.
- `var_gene_name_columns=` on `prepare()`: further `var` columns offered as
  alternative gene names (for example symbols beside an Ensembl index).
- `embedding_keys=` (`--embedding-1d/2d/3d KEY[:c0,c1]`), `expression_source=`
  (`--expression X|raw|layer:NAME`) and `--gene-id-column` for serving AnnData.
  Without `--expression`, a served object reads `X`, or, without `X`, `raw`,
  its only layer, or its best-named layer (`norm`/`log` first, `counts` last),
  and reports that guess -- the rule a file opened in the browser follows.
- `cellucid.data_mapping`: the data-mapping model shared with the web app and
  the R package, with a conformance fixture
  (`tests/fixtures/data_mapping_conformance.json`) and a cross-writer byte
  parity fixture (`tests/fixtures/writer_parity_v3.json`).
- A sixth reconciled surface, **Embedding**: the writer proves `embeddings/`
  holds exactly the files `dataset_identity.json` declares.
- `cellucid serve`, `serve_anndata()` and `show_anndata()` with AnnData publish
  the data mapping they serve at `/_cellucid/mapping/candidates.json` (as
  `defaults`, beside every candidate `obsm`, `obsp`, layer and `var` column,
  read from metadata only), and `AnnDataAdapter` gains `get_data_mapping()` and
  `get_mapping_candidates()`. The web app reads it to record the mapping a
  session was saved with. Every other narrow `obsm` array is turned on beside
  the served defaults and served on request at
  `/_cellucid/mapping/embeddings/{obsm}/{columns}.bin` (see the data-mapping
  routes below), so a served file shows
  any of its embeddings and its checklist works as for a file read in the
  browser.
- **Web app, data mapping**: a **Data mapping** group at the bottom of the
  Session panel lists every embedding the dataset can be shown with as a
  checklist: turn any number on, of any mix of 1D, 2D and 3D, and mark one per
  dimension as its **Default**. For an `.h5ad` or Zarr ZIP opened in the
  browser every numeric `obsm` array is listed, a wider one at any selection
  of 1 to 3 of its columns (`X_pca` columns 2 and 3 as `PCA 3,4`), as many
  selections as you like; a prepared dataset lists the embeddings it was
  exported with, all on. Turning an embedding on loads nothing until a view
  shows it. Choices are staged and applied as one step that loads and checks
  every embedding a view must switch to first; the neighbour graph, genes,
  latent space and vector fields have disclosures of their own (below). The
  data mapping's embeddings are now `{enabled, defaults}`, in
  the web app, `cellucid.data_mapping` and the R package's conformance fixture.
- **Web app, several embeddings at once**: every view shows one of the
  embeddings that are on, of any dimension — two 3D embeddings side by side, a
  3D one beside a 2D one, a t-SNE beside a UMAP — chosen in **View N shows**
  under **Compare views**, grouped by dimension with the default first. Keys
  `1`, `2` and `3` show that dimension's default in the selected view and,
  pressed again, its next embedding ("2D · t-SNE (2 of 3)"). Each view row
  names what its view shows (`UMAP · 3D`), its field, its navigation mode, and
  whether its camera is linked. Views of one embedding share its coordinates
  and GPU buffers; coordinates no view shows are kept in a cache bounded by a
  quarter of the device memory (512 MiB to 2 GiB) and released least recently
  shown first. Up to four views (the live view and three kept views).
- **Web app, per-view cameras**: each view owns its camera and navigation
  mode. **Navigation** is its own accordion and edits the selected view; its
  **Applies to** row names that view and the views linked to it. **Link
  cameras** (on by default) links only views that show the same embedding, and
  is off and disabled while no two views do.
- **Web app, rendering**: points are stored on the GPU in the draw order of
  their spatial tree (by LOD band, then tree leaf or coarse node, then cell),
  in every dimension, whenever the tree exists. A reduced level of detail
  draws a prefix of the stored points, frustum culling draws whole
  tree-node ranges (one multi-draw call of at most 1,024 ranges, 256 one call
  a range without `WEBGL_multi_draw`, bridging the gaps that hide the fewest
  points, with no per-point index upload), and full detail draws faster
  because neighbouring points are neighbours in memory: at 10M points about
  three times faster in 3D and in 2D. LOD no longer keeps a second copy of the
  positions: about 200 MB of GPU memory at 10M points instead of 320. The
  tree, its levels of detail and its draw order are built in a worker, so
  turning LOD or culling on, loading a dataset or showing another
  representation never blocks the page (it took about 2.2 s at 10M points);
  the view draws full detail until the tree lands, then adopts its order in
  one frame, gathered on the GPU. Colours are gathered on the GPU too, so a
  recolour costs one upload of 4 bytes a point and one GPU pass whether LOD is
  on or off: 5-8 ms at 10M points (Chromium), where a recolour with LOD on
  used to take about 68 ms. A view switched back to a representation it
  showed keeps its colours without an upload, and a representation switch
  back to a layout already shown no longer rescans its positions. Translucent
  points in 2D overlap in tree order, the same in every view, every level of
  detail and figure export.
- **Web app**: the in-browser position budget (512 MiB where the browser
  reports no device memory) never refuses the embeddings the views show plus
  the one being shown; it only evicts cached embeddings no view shows.
- **Web app**: a notification can carry one action button; the best-guess
  notification uses it (**Show**) to open the Data mapping section.
- **Web app, usage analytics**: `data_load` reports how the embeddings were
  chosen (`mapping_origin`: `default` or `best_guess`) and how many the
  dataset could show (`embedding_candidates`: `1`, `2-3`, `4-9`, `10+`);
  how many it opened with turned on (`enabled_representations`); `feature_used`
  counts `data_mapping` (`embedding`, `embedding_columns`, `default`) and
  `dimension` when the picker switches a view to another dimension. No key or
  embedding name is ever sent.
- **Further data-mapping routes** beside that catalog: on request, the server
  serves any catalog entry -- an embedding at chosen columns,
  a vector field bound to an embedding, a neighbour graph, the gene list and
  values of `X`, `raw` or a layer named by any `var` column, and outlier
  quantiles in any latent space. Entries are addressed by catalog position,
  never by key. An unknown index is 404; an entry whose values fail their check
  is 422 with `{"error", "key", "reason"}`, so a graph that is not symmetric is
  refused, not repaired. `reason` is the short clause (`NaN at cell 17`) a file
  opened in the browser gets in the same words, and the catalog lists an array
  the best guess already refused as unusable with that reason, as the browser
  does. HEAD answers lengths without reading payloads, and one
  byte-bounded, lock-guarded cache serves the server thread and the notebook
  together. `AnnDataAdapter` gains the `get_mapped_*` methods behind them.
- **Jupyter:** the `setDataMapping` command and the `data_mapping` event, both
  listed by `/_cellucid/protocol`. `viewer.data_mapping` (and
  `viewer.data_mapping_origin`) hold the mapping the viewer last reported,
  `@viewer.on_data_mapping` runs on each report, and
  `viewer.set_data_mapping(**fields)` changes the embeddings, graph, gene names,
  expression, latent space or vector bindings from the notebook, validated before
  anything is sent.
- **Web app, served data mapping**: a file served by `cellucid serve` or shown
  with `show_anndata()` can show any numeric `obsm` array, at any columns, from
  **Session › Data mapping**, as a file opened in the browser can; the server
  reads the array on request, and its centroid labels are computed from its
  positions in the viewer, as for a file read in the browser. A choice
  the server refuses is listed unusable with the server's reason and the
  previous mapping stays.
- **Web app, neighbour graph**: for an AnnData file (in the browser, served or
  in a notebook), **Session › Data mapping › Neighbour graph** chooses the
  graph connectivity edges and the KNN drag tool use: None or any square
  `obsp` matrix with one row per cell. The graph is read and checked when it is
  applied; one that is not symmetric (`obsp['distances']`) is refused in the
  server's words and stays unavailable, never repaired. Shown edges are redrawn
  from the new graph and **Dataset › Connectivity** counts its edges. A session
  restores the graph it was saved with when the dataset has it. In the browser,
  a file keeps its own graph and at most one other in memory, and their edges
  are byte-identical to what `cellucid serve` sends for the same file.
  `feature_used` counts `data_mapping` / `knn`.
- **Web app, gene names and expression**: for an AnnData file (in the browser,
  served or in a notebook), **Session › Data mapping › Genes** chooses the
  expression matrix (`X`, any layer, `raw.X` with its own gene list, or None)
  and the gene names (the var index or any text, categorical or integer var
  column; a floating-point or true/false column is listed with the reason it
  cannot name genes). Empty names take the var-index name and repeated names
  are made unique (`MATR3`, `MATR3-1`) by the rule `prepare()` and the server
  use, and the panel says what that changed and when Ensembl ids have a symbol
  column beside them. Apply reads the new list first; a view coloured by a
  gene keeps it with its colour range, log scale and filter (by position while
  the var axis stays, as between `X` and its layers; by name to or from
  `raw.X`), a gene not on the new
  list loses its colouring and is named, cells and highlights are untouched,
  and analysis windows close. A prepared dataset switches between the
  gene-name sets it was exported with. A session restores the gene names and
  matrix it was saved with and checks its gene count against that matrix.
  The browser reads one matrix at a time (`X`, `layers/<name>` or `raw/X`),
  and its gene lists and values are byte-identical to what `cellucid serve`
  sends for the same file. `feature_used` counts `data_mapping` /
  `gene_names` and `expression`.
- **Web app, latent space**: for an AnnData file, **Session › Data mapping ›
  Latent space** chooses the `obsm` array (or None) a categorical field's
  outlier filter measures distance in. Apply reads and checks the space and
  computes the quantiles of every loaded categorical field before anything
  changes, off the main thread (`data/latent-outliers.js`), bit for bit what
  `prepare()` and `cellucid serve` compute; a space with a bad value is refused
  in the server's words and stays unavailable. Every view keeps its outlier
  threshold; None removes the outlier filters. A space larger than 512 MiB as
  float32 is listed with that reason, as the server lists it.
- **Web app, vector fields**: **Session › Data mapping › Vector fields** binds
  vector arrays to embeddings with one checkbox per array that fits each
  embedding: as wide as its dimension (read whole) or as wide as its array
  (read at the same columns, `velocity_pca` on `X_pca[2,3]`). Arrays named for
  the embedding are marked suggested. A prepared dataset lists the bindings it
  was exported with, any of which can be turned off. Apply rebinds and turns
  the overlay off; vectors are read, checked and scaled only when the overlay
  shows them, through `_cellucid/mapping/vectors/...` on every AnnData path,
  and the overlay holds GPU data only for what it shows. `feature_used` counts
  `data_mapping` / `latent` and `vector_fields`.
- **Web app, Jupyter**: the viewer reads `GET /_cellucid/protocol` when it
  connects. It applies `setDataMapping` as one transaction and answers it with
  a `data_mapping` event naming the request; it reports every other mapping
  change, and reports a refused notebook command as `command_error`, to a
  `cellucid` that lists those events.
- `/_cellucid/info`'s `version` is the documented signal a newer viewer uses to
  tell a user their server predates the viewer's export format ("This Cellucid
  server (0.9.2) is older than the viewer; update it with pip install -U
  cellucid.").

### Changed

- **Latent-space outlier quantiles** follow one stated rule in `prepare()`,
  `cellucid_prepare()` in R (now in C), `cellucid serve` and the web app:
  float32 values widened to float64, centroid sums accumulated in cell order,
  squared distances in column order, and the share at that distance or closer
  as float32 (numpy's pairwise summation and R's long-double sums are no longer
  part of the result). `tests/fixtures/latent_outlier_conformance.json`
  (written by `scripts/generate_latent_outlier_conformance.py`) pins the bytes
  in all three repositories. Quantiles can differ from 0.9.1's in the last rank
  only where two distances tie within rounding.
- **`prepare()` gene names from an integer column**: a `var` column of integers
  names genes by their digits (`7157` is `"7157"`), the rule the viewer,
  `cellucid serve` and `cellucid_prepare()` in R apply to the same column,
  instead of being refused. A floating-point or boolean column is refused with
  the catalog's reason (`floating-point values cannot name genes`).
- **`cellucid serve` / `show_anndata()` vector fields**: without
  `--vector-fields` (`serve_vector_fields=True`) the served data mapping binds
  no vector field, matching what the server declares; with it, the default
  rule binds `<field>_umap_<d>d` as before. The startup line now reads
  `Vector fields: not bound (... to bind them, or bind one under Session > Data
  mapping)`.

- **Export format version 3.** One export can now carry several embeddings per
  dimension and several sets of gene names. This is a clean break: exports in
  the earlier format are not read, and format-3 exports need a web app that
  reads format 3. Re-export with `prepare()` (or `cellucid_prepare()` from the
  matching R release, which writes the same bytes). `dataset_identity.json`
  has `version: 3`, and `embeddings` is now `{default_dimension, defaults, items}`: every exported
  embedding is an item with a `label`, a `dimension`, the exact payload `file`
  `embeddings/{index}_{dimension}d.bin(.gz)` and its `source` (`key` and the
  0-based `columns` it was read from). Items are ordered by dimension, then id;
  `defaults` names the item each dimension shows first. The `points_*d.bin`
  files at the export root are gone.
- **Vector fields bind to embedding ids.** A field's `files` are keyed by
  embedding id, `vectors/{fieldIndex}_{embeddingIndex}.bin(.gz)`, scaled by that
  embedding's normalization factor. `basis`, `available_dimensions` and the
  per-field `default_dimension` are removed; `label` is the field id.
- **Manifests are `compact_v2`.** A categorical field's element `[5]` is
  `centroidsByEmbedding`, keyed by embedding id. `var_manifest.json` replaces
  `var_gene_id_column` with `gene_name_sets`: the default set (from
  `var_gene_id_column` or the var index) and any alternative sets, stored as
  metadata only, so a reader can switch gene names without reading expression
  again.
- **`prepare()` API:** `X_umap_1d`, `X_umap_2d` and `X_umap_3d` are replaced by
  `embeddings={id: array or Embedding}` (required) and `embedding_defaults`
  (required for every dimension with more than one embedding). The new
  `cellucid.Embedding(values, columns=None, label=None, source_key=None)` reads
  chosen 0-based columns of a wide array, such as `X_pca` columns `[0, 1]`.
  `vector_fields` is now `{field: {embedding id: array}}`.
- **Gene names are made unique instead of refused.** Repeated names become
  `MATR3`, `MATR3-1`, ... by the rule anndata's `var_names_make_unique()` uses,
  and an empty name takes its row's var index name. The same rule runs in the
  Python writer, the R writer and the AnnData server, and what changed is
  printed and recorded in `var_manifest.json`. `gene_identifiers` matches the
  unique names.
- **`cellucid serve` / `serve_anndata()` / `show_anndata()` for AnnData** serve
  identity v3. `serve_connectivity` and `--connectivity` are replaced by
  `knn_key` and `--knn KEY`: the named `obsp` matrix is read and checked when it
  is named, and one that is not symmetric (such as `distances`) is refused with
  that reason instead of being symmetrized.
- When an object has no `X_umap_1d`/`X_umap_2d`/`X_umap_3d` and no 1- to
  3-column `X_umap`, the server no longer stops with rename hints. It makes a
  deterministic best guess from the other `obsm` arrays (t-SNE, force-directed
  layouts, PHATE, spatial coordinates and diffusion maps are preferred over
  latent spaces such as `X_pca`) and reports it; when nothing can be plotted it
  says what `obsm` holds and why each entry was not used.
- Sessions: `datasetFingerprint` must carry `cellOrder` (now naming the
  `embedding` its digest was taken on) and the new `dataMapping`. User-defined
  categorical fields no longer store `centroidsByDim` or `normalizedDims`.
  `apply_cellucid_session_to_anndata()` checks the mapping against the AnnData,
  compares `varCount` on the mapped gene axis, and resolves var-sourced fields
  through the mapped gene names and expression matrix.

- **Web app, data**: the viewer reads export format 3 and only format 3 — from a
  prepared folder, the sample datasets, GitHub, a server, and the notebook —
  and opens `.h5ad` and Zarr files as the same format. Every embedding an
  export carries is loaded only when a view shows it; each dimension opens on
  the embedding `defaults` names. An `.h5ad` or Zarr file's vector field now
  reports its id as its label (`velocity_umap`, not `Velocity (UMAP)`), as the
  format requires, and the reader now refuses any export whose vector-field
  label is not its id, from every source.
- **Web app, data**: an `.h5ad` or Zarr ZIP without a UMAP is no longer
  refused. It opens on the same deterministic best guess `cellucid serve`
  makes (t-SNE, force-directed layouts, PHATE, spatial coordinates, …, and last
  the first two columns of a wide array such as `X_pca`), marked *guessed*,
  with a notification naming what is shown. Only a file with nothing drawable
  is refused, with the message `cellucid serve` prints: "Cellucid found
  nothing to plot…", listing what `obsm` holds and why each entry was not used
  (the error code is `NO_PLOTTABLE_EMBEDDING`). Every embedding is checked
  before it is drawn — NaN, infinity, a value outside 32-bit float range, or
  no coordinate variation is refused naming the first cell — so a one-cell
  file can no longer be opened. A one-dimensional `obsm` array is read as one
  column, and a wide array is read only around the columns a view needs.
- **Web app, vector fields**: the overlay is offered by the embedding a view
  shows, not by its dimension: a field is drawn only in views showing the
  embedding it is bound to, and the info line names those embeddings
  ("Vector fields are bound to … Show one of them in this view to enable.").
- **Web app, sessions**: a session records its data mapping (the embeddings
  that are on and each dimension's default), the embedding its cell order was
  taken on, the embedding every view shows, every view's camera and navigation
  mode, and whether cameras are linked. Loading checks the dataset, its data
  mapping and the cell order of the embedding the session names before
  anything changes on screen, loading that embedding if another is shown — so
  a session saved in 2D now opens while 3D is shown, instead of asking you to
  switch back. A session whose embeddings this dataset can show turns them on
  before its views are restored, and one naming an embedding the dataset does
  not have is refused naming it ("This session shows obsm['X_tsne'] columns
  [0, 1], which this dataset does not have."). A
  session saved by an earlier release is refused with a message saying so.
  User-defined categorical fields no longer store their centroids; they are
  recomputed from the codes.

### Removed

- `vector_fields.missing_embedding_message()`, `resolve_embedding_keys()`,
  `embedding_declaration_hints()` and the other helpers of the one-UMAP-per-
  dimension rule, together with `AnnDataAdapter.get_points_binary()`.

### Fixed

- **Web app, analysis performance**: analysis no longer freezes the page on
  large selections, and every result is the same number it was. Medians,
  quartiles, ranks and sorted copies (Quick, Detailed summaries and tests, CSV
  export, Spearman correlation, Kruskal–Wallis, high-cardinality Wilcoxon in
  differential expression and Marker Genes, gene signatures) come from
  linear-time radix ordering instead of comparison sorts; correlation pairs its
  colour field by merging sorted cell indices instead of a per-cell map; box,
  violin and density plots of large pages compute their statistics on the
  analysis workers; a histogram hands the plot a few values per bin instead of
  every value; and a box plot draws at most 16 copies of an identical outlier
  (all on one spot) instead of one SVG element per cell. At one million cells:
  box plot of a zero-inflated gene 1.1 s → 0.12 s (longest main-thread task
  70 ms), histogram 0.17 s → 0.06 s, Detailed summary 0.10 s → 0.03 s,
  three-page Kruskal–Wallis 0.67 s → 0.09 s, Spearman 0.86 s → 0.14 s,
  high-cardinality marker statistics 0.33 s → 0.06 s per gene. The Plotly
  notification watcher no longer runs a subtree query for every element any
  plot inserts.
- **Web app, data loading**: when the connection to a server is lost, Cellucid
  clears only the dataset that came from it, and not while another dataset you
  chose is still loading; it used to clear over that choice. Recovering from a
  lost connection while its dataset was on screen also failed outright
  ("Remote dataset URL requires an active connection"), because the open
  dataset's address was re-read from the closed connection; it is now the
  address the dataset was loaded from. A server with one dataset that answers
  **Connect** after you applied a session no longer replaces that dataset and
  session; its dataset is listed in the picker instead.
- **Web app, data loading**: a dataset chosen while another is still loading
  now replaces it. At startup the data controls stayed disabled until the
  startup dataset had downloaded, and **Load State** did nothing at all in
  that window; afterwards, a file still being read published over a sample or
  **None** chosen after it. Now the controls open once the dataset catalog has
  been read, and every dataset change (the startup dataset, a sample, **None**,
  a file or folder, a server or GitHub repository with one dataset, the
  benchmark) supersedes the load in flight the moment you act: its downloads
  are cancelled, it never appears late, and the dataset picker shows what is
  loading. A connection that answers after a newer choice keeps its dataset in
  the picker instead of replacing that choice, and a session chosen during a
  load waits for that dataset, says so, and is then applied.
- **Web app, rendering**: switching from **Grid (light)** to **Grid (dark)**
  no longer leaves the scene redrawing on every frame (the grid fade only
  moved upwards, so it never reached the lower target). With reduced motion
  requested by the system, background fades and crossfades are instant.
- **Web app, rendering**: connectivity edges declare their per-vertex
  attribute at location 0. Left to the linker, Firefox on macOS could draw
  them without an enabled attribute 0 and emulate one per vertex.
- **Web app, no dataset**: after choosing **None**, Keep view, figure Export,
  Save and Connect repo are unavailable and say why instead of throwing;
  Differential expression, Gene signature and Marker genes show `Requires gene
  expression data.` (as on a dataset without genes) instead of a page error;
  Analyze performance and Copy situation report ask for data; and Load says no
  dataset is open instead of quoting `0 cells and 0 genes`.
- **Web app, analysis**: marker-gene results are cached per gene list, so
  results computed with one set of gene names or one expression matrix are
  never shown for another.
- **Web app, data**: a dataset re-exported at the same address can no longer
  show values from its previous version. Gene and observation payloads were
  reused from the browser cache without asking the server, and payload files
  are numbered, not named: a re-export that added one gene left
  `var/3.values.u8.gz` holding another gene of exactly the same size, so the
  cached file drew the old gene's values under the new gene's name. Every
  dataset file — catalogs, manifests, payloads and saved views — is now
  revalidated with the server before use; an unchanged file costs a short
  "not modified" reply.
- **Web app, interface**: one calmer design. One font (Inter) in sentence case
  with a small type scale; soft borders; compact spacing. Every panel has
  short headings, and every explanation moved from the panel into an (i) tip
  that links to the exact docs page. Navigation is its own panel and shows
  which view it drives when several are open; the benchmark panel is called
  Benchmark; Visualization, Navigation and Compare views start collapsed; the
  Outlier filter sits under Legend; the logo follows the theme. Slider thumbs
  no longer overhang the rows around them, and the legend tick no longer
  touches its swatch.
- **Web app, rendering**: a hover lookup and any click on the page no longer
  redraw the scene for a second.
- **Web app, picking**: hovering and clicking now pick the dot drawn under the
  pointer, at the size it is drawn, instead of the first cell near the view
  ray, which could be a nearer cell beside the pointer. Where dots overlap the
  front-most wins; filtered-out cells and cells the current level of detail
  does not draw are never picked. Cell info, Alt+click selection, proximity
  and KNN all use it, each in the pane under the pointer.
- **Web app, selection**: in a grid of views, Alt+click and Alt+drag in a
  pane other than the focused one now select by that pane's field. They used
  the focused view's field, so a click in a pane colored by `cell_type`
  beside a focused view colored by a continuous field selected a value range
  of the wrong field. The category, the value range, the drag preview and its
  range label now all come from the pressed pane, among the cells it shows;
  view focus still does not move. A pane colored by no field now says so.
- **Web app, selection**: whether an Alt press can start, and the armed
  cursor, now follow the pane under the pointer instead of the focused view,
  so Alt+click works in a colored pane while the focused view has no field,
  and the panel asks for a field only when no pane on screen has one.
- **Web app, cell info**: the card no longer disappears for good. A key whose
  release the page never saw (a Cmd shortcut such as the macOS screenshot
  keys, or a key held while the window lost focus) read as camera motion
  forever, and so did a drag released outside the page; the card waits for
  the camera to settle, so it never came back. Held keys are now only the
  navigation keys, are dropped on blur, and Shift alone no longer counts as
  motion. The card no longer repeats the pane's name in a grid of views.
- **Web app, benchmark**: the presets are now 500K, 1M, 5M, 10M, 20M and 40M.
- **Web app, welcome**: the welcome dialog opens with no choice pre-selected;
  focus rests on the dialog and the first Tab reaches **Choose a dataset**.
  The sidebar footer links GitHub next to Docs and Cite.
- **Web app, onboarding**: closing the welcome dialog, by any path, no longer
  focuses the **Sample dataset** dropdown; **Choose a dataset** opens the
  **Session** panel, and from the keyboard moves focus to its heading.
- **Web app, data loading**: dataset catalogs, identities and manifests are
  revalidated on every load, so for a few minutes after a dataset is updated
  a stale cached copy can no longer clash with the new files (for example a
  sample view failing with a missing control).
- **Web app, pointer**: hover, click and selection read the camera as it is,
  not as the last frame drew it, so a pointer read right after the camera
  moved no longer lands on the previous picture.

- **Web app**: automatic antialiasing now accounts for the view's dimension.
  A 1D or 2D view has no depth test, so every point is shaded and antialiasing
  costs about ten times more there (1M cells: 33.8 ms with it, 3.4 ms without).
  It is now on below 5,000,000 cells in 3D and below 200,000 cells in 1D/2D,
  and the automatic choice is re-decided when the active view switches
  dimension. The status line names both (`Chosen automatically: 3,696 cells in
  3D.`); an explicit choice is never changed.
- **Web app**: less memory per cell, so larger datasets fit. Per-cell visibility
  is one byte (exactly 0 or 1) instead of four, the outlier quantiles are no
  longer copied into the state, categorical codes keep their stored width, and
  a kept view shares the visibility arrays copy-on-write instead of cloning
  them. Loaded fields are released by bytes rather than by count, so browsing
  genes no longer grows memory without bound. A loaded 3D dataset went from
  36 to 22 bytes per cell, and each kept view from 26 to about 6; 150M cells in
  2D and 120M in 3D load and render on a 16 GB machine. Pixels, filters and
  exports are byte-identical. The benchmark's pre-publication GPU estimate now
  follows the dimension (13 B/point in 2D, 17 in 3D) instead of a fixed 28.

- **Web app**: usage analytics rebuilt around one declared schema
  (`assets/js/analytics/schema.js`). It no longer sends dataset names, ids built
  from local file names, exact counts, error text or full page URLs (the page
  view carried `remote=` and `github=` addresses and their tokens); a dataset is
  described by count ranges, dimensions and writer version, with an id only for
  the public samples. Each load is reported once (it was reported twice), timed
  from the click or file pick (local loads excluded the file read), and reported
  when it fails before activation (such failures sent nothing). New events cover
  startup success and failure, feature adoption, and device classes. Every
  payload is validated against GA4's limits, and analytics is never awaited.

- **Web app**: faster, quieter full-detail rendering with identical pixels.
  Hovering the plot or typing into a sidebar field redrew the whole scene for a
  one-second settle window (237 full frames over two seconds of hover at 1M
  cells); only input that changes something now draws. In 2D views — drawn
  without depth testing — the point shaders use a single `discard`, and `Light`
  quality uses a no-depth program that blends a zero-alpha fragment instead of
  discarding: at 1M cells with antialiasing, `Light` 30.7 → 17.0 ms and `Full`
  33.8 → 30.6 ms per frame, byte-identical in 45 compared scenes. Orbit and pan
  inertia, smooth zoom and held-key acceleration now scale with elapsed time, so
  a 120 Hz display no longer moves the camera at double speed, and the render
  loop allocates about half as much per frame.

- **Web app**: with an ad or tracking blocker (uBlock Origin, Brave, AdGuard, iOS
  content blockers — anything shipping EasyPrivacy, which lists
  `/analytics/tracker.js`), Cellucid did not start and showed “Cellucid could not
  load its code…”. The analytics module was a static import, so its refused
  request aborted the whole module graph. Analytics is now loaded with a dynamic
  import and stays off when refused; the dataset-load vocabulary it shared with
  the controls moved to `app/data-load-methods.js`. No other application file
  matches the EasyList, EasyPrivacy, uBlock Origin or AdGuard lists.

- **Web app**: choosing volumetric smoke on a device without float blending
  (Safari on iPhone and iPad, many Android browsers) failed mid-build with
  “Smoke density unavailable: GPU smoke density requires EXT_float_blend for
  exact R32F accumulation.”, and restoring a smoke session there threw. Smoke's
  requirements are now checked up front; such a device gets a plain message
  (“Volumetric smoke isn't available on this device: …”) and stays in Points,
  sessions restore into Points, and other smoke errors no longer show internal
  detail.

- **Web app**: the volumetric smoke cloud is now **beta** rather than alpha. It
  is coloured, bounded and idle-aware after this release's rewrite; the tag beside
  `Render mode:` reads `Beta` in the calmer info tone, and its hover text no
  longer says the mode is still changing.

- **Web app**: analysis plots of large pages draw several times faster, with
  identical results. Plotly's violin density evaluated its kernel once per value
  at every curve point; a local, test-pinned patch to the vendored Plotly
  evaluates it once per run of equal values (expression data is mostly ties) and
  skips terms that are exactly zero, adding the rest in the same order. Box plots
  that draw points hand Plotly exact precomputed statistics; density curves are
  computed once (not twice), reused across cosmetic changes, and split across the
  worker pool when large. At 2,000,000 cells: a gene violin 3.7 → 0.67 s, a
  continuous box 1.15 → 0.26 s, a continuous density 1.39 → 0.56 s. Page values
  keep their own numeric type, the page-data cache is bounded at 256 MiB as well
  as 100 entries, and decodes go to the least-loaded worker.

- **Web app**: on browsers that composite in software (blocklisted GPUs, virtual
  machines, remote desktops), part of the plot vanished beside a floating panel
  whose content scrolled: Chromium culled the canvas under the panel's opaque
  scroller and drew the hole at its mirror image. Each floating window is now its
  own compositing surface.

- **Web app**: multiview and the data path are faster. Kept views drew their
  reduced Level-of-Detail levels through a gather from the full-detail streams;
  they now draw the live view's compact positions with their own compact colors,
  and share the live positions instead of copying them. Measured on an Apple M1
  Pro (Chromium, ANGLE Metal), 10,000,000 cells in 2D with LOD: two views
  15.9 → 6.9 ms a frame, four views 43.6 → 13.5 ms, and each kept view now adds
  38 MB of GPU memory instead of 114 MB. Hiding a category took 225 ms instead
  of 413 ms, a cached recolor 225 ms instead of 321 ms, and loading a prepared
  export 1.30 s instead of 1.43 s (30,000,000 cells: 2.89 s instead of 3.45 s).
  The documentation's claim that a four-view grid costs about as much as one
  view held only for large, overlapping points and is corrected.

- **Web app**: a panel's undock, dock, copy and close buttons no longer sit inside
  its header's open/close toggle. A `<summary>` is itself a button, so each was a
  control nested in a control — announced poorly by screen readers and flagged by
  accessibility checkers. They now stand beside the header, look and sit exactly
  where they did, are named after their panel ("Undock Visualization panel into a
  floating window"), and come just before the header in `Tab` order.

- **Web app**: multiview panes are laid out in the part of the window the sidebar
  leaves uncovered. The grid was cut from the whole canvas, so with the sidebar
  open the left pane was partly hidden under it (at 1440 px with two panes, 430 px
  of the left pane was visible against 720 px of the right one). Every pane now
  has the same size beside the open sidebar, and the layout follows at once when
  the sidebar is closed, opened, resized, or the window or its pixel density
  changes. One layout drives drawing, picking, hover, focus, every selection
  tool, projectiles, the pane title chips, centroid labels, figure export, and
  the selection bar, so what you click is always the pane you see. Wheel zoom
  now focuses the pane under the cursor, free-fly is centred in the uncovered
  area like Orbit and Planar, and the hidden sidebar is `inert`, so `Tab` no
  longer walks into it.

- **Web app**: the sidebar's ✕ works from the first frame. It was wired only once
  the first dataset had loaded, so during a slow first load a phone user, who
  sees the menu full-screen, tapped it and nothing happened.

- **Web app**: smaller layout fixes.
  - On a short screen (a phone on its side), startup notifications covered the
    welcome dialog's buttons; while it is open they now wait beneath it.
  - A long view title in a narrow multiview pane ran into the next pane's
    title; it now ends in an ellipsis inside its own pane, with the full name
    on hover.
  - The benchmark's Frame Time reading wrapped onto three lines. The second
    number now sits on its own line, labelled `CPU` for what it measures: the
    main-thread time spent issuing the frame, not GPU time.
  - Figure export's Background and Font selects now have visible labels, with
    the same trailing colon as every other label.

- **Web app**: in multiview, a finger dragging a pane now focuses that pane first,
  as a mouse press does, instead of moving whichever pane last had focus, and a
  touch drag saves the camera to the focused pane when cameras are unlocked.

- **Web app**: the volumetric smoke cloud is rebuilt and fixed.
  - It drew near-black soot whatever the cells were coloured by. The cloud now
    takes the active field's colours: each region shows the mean display colour
    of the cells there, lit by a key light and sky occlusion pre-computed once
    per build rather than by a shadow ray from every sample.
  - Switching dimension or recolouring (field, colormap, category colour, log
    scale) left the previous volume on screen. Both now rebuild it, and a 2D or
    1D embedding is spread into a soft slab instead of a sheet rays stepped over.
  - `Ray quality` changed almost nothing. It now sets samples per voxel, 0.5 to
    2 (`High`, 1.41, by default), and empty space is skipped at no cost.
  - The cloud shimmered and the GPU never idled: the ray-start jitter moved
    every frame and every frame was redrawn. The jitter is now fixed, and at
    `Animation speed` 0, which reads `Still`, a still view is not redrawn.
  - `Render resolution` reached 2.00x and `Noise detail` 256³, either of which
    could stall one frame long enough to lose the WebGL context. Render
    resolution now runs 0.25x–1.00x (default 0.50x) and a frame marches at most
    2560×1440 rays; noise detail is 32³–128³ (default 128³).
  - `Cloud density`, `Fine detail` and `Light absorption` printed readouts that
    changed with grid and noise resolution, so a recorded value did not
    reproduce a look. Every smoke slider now has one fixed scale.

  The documentation's control tables, recovery order and pipeline description
  are corrected (grid size was never the cause of context loss: the grid is
  capped at 128³, an 8 MiB volume), and the sample datasets' saved views now
  open with the new smoke defaults.

- **Web app, user-defined fields**: a categorical made with
  `Create Categorical`, a category merge or delete, or a field duplicate now
  shows its category labels, keeps them when a view switches dimension (that
  switch threw `Raw 2D positions are unavailable for centroid computation`),
  and gets them back from a saved session exactly where they were; restoring a
  session whose active field is one of these no longer fails with
  `Active categorical field is absent from its rendered selector`.
- **Web app, sessions**: a dataset without genes can be saved to a session;
  capture used to fail with `Dataset fingerprint requires the current var
  field inventory`. Such a dataset now has an empty gene inventory rather than
  none.
- **Web app, camera**: at a low frame rate (a large dataset with level of
  detail off) inertia, zoom and the grid fade now run in real time instead of
  slowing down below 20 fps; only a stall over a quarter second is cut short.
- **Web app, rendering**: switching to another window or tab no longer redraws
  the whole scene for a second; focus changes redraw only when they end a held
  key or a selection in progress.
- **Web app, touch**: on an iPhone or iPad, pressing and holding the plot no
  longer selects the whole page as text or opens the callout menu; the press
  reaches the viewer (Cell info's hold, Select mode, the camera). Safari reads
  only `-webkit-user-select`, so every non-selectable control now declares it
  too.
- **Web app, highlights**: a highlighted selection no longer keeps the whole
  scene redrawing every frame in `Light` and `Ultra-light` shader quality,
  where the ring does not animate; the view goes idle as it does without a
  selection. `Full` still pulses its ring.
- **Web app, fog**: points, highlight rings, connectivity edges, the orbit
  anchor and the velocity overlay share one fog model, evaluated once per
  point instead of per pixel, so every layer of a scene fades at the same rate.
  Fog now thins a point toward whatever is behind it instead of tinting it
  with a fixed fog colour, so fogged points fade correctly over a grid or a
  vignette; over a plain background the picture is unchanged. With fog off a
  dense 10M-cell view draws about 8% faster (Apple M1 Pro).
- **Web app, 2D views**: point edges in 1D and 2D views are smooth (a
  one-pixel antialiased rim) in `Full` and `Light`, instead of a jagged
  staircase that the antialiasing setting could not smooth; `Full` draws 2D
  views as cheaply as `Light`, and figure export draws them the same way.
- **Web app, rendering**: connectivity edges take about a third less GPU time
  (all 6.3M edges of the Suo sample: about 190 to 115–128 ms a frame) and the
  volumetric smoke cloud about 15% less (17.2 to 14.5 ms), on an Apple M1 Pro.
  Every shader now pins a vertex attribute to location 0, which Firefox
  otherwise emulates per vertex. Projectiles no longer allocate new vertex
  arrays every frame.
- **Web app, projectiles**: projectiles bounce off the **Studio** floor under a
  3D view, as they bounce off the grid's walls, and land on it instead of
  falling through. A projectile that had come to rest against a cell no longer
  falls through the grid's floor either.
- **Web app, rendering**: after a change the viewer draws two more frames
  and goes idle, instead of at least twelve frames and a full second; an
  adaptive LOD still stepping toward its level keeps the scene drawing on its
  own account, so nothing is cut short. A fast display no longer draws 60–120
  identical frames after every click, and where a frame takes seconds (a
  software renderer) every change costs seconds less.
- **Web app, smoke**: generating the cloud's noise volumes no longer queues
  GPU work faster than the GPU finishes it. Each batch of slices waits for the
  previous one and shrinks while the GPU lags, so on a slow or software
  renderer the page stays responsive while smoke mode prepares (before, a
  click there could wait half a minute behind the queued noise work).
- **Web app, selection tools**: while a lasso, proximity or KNN selection is
  drawn, a pointer move that leaves the previewed cells unchanged no longer
  re-uploads the highlight or redraws the scene; only a change does.
- **Web app, sessions**: a session saved while a view was coloured by a
  continuous field could not be loaded ("Session UI control … does not
  exist"): it recorded the legend's range sliders, whose ids are made afresh
  on every page. Controls with page-made ids are no longer recorded; the
  colour range and filter they show are saved with the field, as before.
- **Web app, shader quality**: `Full` shading is a little brighter on the
  side of each cell away from the light (ambient 0.5, was 0.4); the lit side
  still shows the colour itself. Highlight rings and centroid markers use the
  same lighting.
- **Web app, Studio backgrounds**: the floor's rules and the dots are lighter
  and finer (major rules no thicker than minor ones), so the stage reads as a
  calm, clean surface behind the data.
- **Web app, rendering**: the plot is drawn at exactly the screen pixels it
  covers. Its drawing buffer used to be rounded down from a fractional layout
  width (a plot beside the sidebar is often half a pixel wider), and the
  browser then rescaled every frame by a fraction of a pixel, which softened
  one-pixel lines and dots.

## [0.9.2] - 2026-10-17

Version 0.9.2 is a patch release on 0.9.1.

### Fixed

- A refused upload or event now reaches the client as the status the server sent.
  Every early refusal of a POST -- an expired or unknown session-bundle request, bad
  viewer credentials, a wrong `Content-Type`, an oversized body, an unknown path, a
  foreign `Host` -- was answered before the request body was read, and the socket
  was then closed with that body unread. A socket closed with unread input is reset
  rather than shut down, so a client still sending or not yet reading lost the
  response: Windows reported `WinError 10053`, macOS and Linux a reset or a broken
  pipe. The server now drains a refused request's declared body, within the
  largest accepted upload and for at most five seconds, before it closes. Both
  servers share one `do_POST` that does this.

## [0.9.1] - 2026-10-17

Version 0.9.1 is the PyPI submission release.

### Added

- Scanpy's own `X_umap` is served as written: when no `X_umap_1d`/`X_umap_2d`/`X_umap_3d`
  key exists, a bare `X_umap` of 1, 2, or 3 columns is read at the dimension its column
  count states, without renaming or writing to the caller's object. A dimensional key
  still decides its own dimension and still takes precedence; any other width is refused
  with its shape named. The same rule governs `add_transition_drift_to_obsm` and the
  browser's h5ad and Zarr readers, so one file opens the same way everywhere.
- A long build reports its phases. The direct-AnnData server runs five reported steps, the
  adapter reports obs classification, embedding resolution and which key resolved it,
  vector-field discovery and the graph read, and the manifest and centroid pass is a step
  of its own instead of silence after a success line. `AnnDataAdapter` takes `quiet`.
- `GET /_cellucid/protocol` answers with the wire capabilities the installed package
  accepts, so a web build can ask instead of assuming.
- A notebook learns when a command it sent was refused, rather than reading `None` as
  success.
- Release-contract validation for package, citation, documentation, and downstream recipe
  metadata; installed wheel and sdist gates on Linux, macOS, and Windows before
  publication; reproducible sdist normalization.
- Antialiasing is a live setting, applied on the next frame instead of at context
  creation, and defaults from the cell count until the checkbox is touched.
- A dataset opens at a point size derived from its cell count; the slider now
  reaches 0.050 and every stored session keeps the size it recorded.
- The volumetric smoke render mode is marked `alpha` in the viewer.

### Changed

- **Web app**: the synthetic benchmark generator's custom point-count ceiling is
  raised from 50,000,000 to 250,000,000. The previous value was a declared product
  limit chosen well below the hardware's; the new one is the hardware's, sitting
  just inside the `MAX_TEXTURE_SIZE²` bound on the renderer's interaction-state
  texture (268,435,456 on a context reporting 16,384). A request that the machine
  cannot allocate now fails as an allocation error rather than being refused in
  advance, because the limit that actually binds is the machine's memory — 16 bytes
  per point, so 4 GB of source arrays at the new ceiling — and the control cannot
  know what that is. The point-count input, both validators and the worker boundary
  all read the one constant, so they cannot drift apart.
- The neighbor graph and the vector-field overlay are read only when asked for:
  `serve_connectivity=True` / `--connectivity` and `serve_vector_fields=True` /
  `--vector-fields` on the direct-AnnData path. Each is otherwise never touched,
  and a dataset serving neither publishes exactly what an object holding neither
  publishes. Both are read and validated in full before the socket binds, because
  the manifests the viewer fetches first declare what they contain — and on a
  large object that is the longest part of startup, for a capability most sessions
  never turn on. A 50-neighbor graph over 18 million cells is roughly 900 million
  stored neighbors. A graph asked for and absent is reported; vector fields asked
  for and undeclared serve none, which is ordinary; `vector_field_default` is
  refused unless the fields it selects among are served. The startup report names
  each state separately. `prepare()` is unchanged in Python and R — both already
  took `connectivities=` and `vector_fields=` as data — and now announces all
  three optional capabilities alike when they are not requested.
- Connectivity validation is measurably cheaper: degrees are counted with `np.bincount`
  rather than the unbuffered `np.add.at` (about twenty times faster on that step), one
  coordinate expansion feeds every check, and a graph already in compressed form is
  compared against its own transpose instead of being rebuilt first. Measured on a
  9M-stored-neighbor graph: 1.2x faster end to end and 36% lower peak memory. A graph
  above 50,000,000 stored neighbors logs what it is about to cost.
- The tested Python range is 3.11 through 3.14, and direct AnnData operation targets
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication uses GitHub OIDC trusted publishing and accepts only a pushed tag
  matching package version metadata exactly.
- Export format: every scientific payload is named by its integer index or a fixed
  neutral name, so no filename carries an observation key, a gene name, or a field id.
- The notebook viewers, the session applier, the data server, and the export writer are
  packages rather than single files. Every public name is imported from where it was, and
  the servers no longer import numpy at startup.
- Adaptive level of detail is bounded by a point budget rather than a fixed ratio.
  A cloud already inside the budget is never reduced, and one above it stops at the
  budget instead of falling to a 44x reduction on any pull-back: at 18.1M cells that
  is 2.4M points rather than 412,000.
- Served gene columns are cached under a byte budget rather than a count of a
  hundred columns, which on an 18.1M-cell object was 7.3 GB of resident memory.

### Fixed

- `cellucid serve --help` no longer fails with `UnicodeEncodeError` on a Windows console
  using a single-byte code page such as CP1252. argparse wrote help, usage, version and
  error text straight to the stream, bypassing the escaping every other console line
  already gets; the CLI's parsers now route that text through the same writer, so any
  character the console cannot show is printed as a backslash escape instead of
  aborting the command. The `serve` help text itself is now plain ASCII.
- A bind address is no longer printed as a URL. `--host 0.0.0.0` used to print
  `http://0.0.0.0:<port>` and try to open it; the banner now prints the loopback origin,
  which is what a browser on the serving machine and the near end of an SSH tunnel both
  use, plus the machine's own name for clients elsewhere. An IPv6 literal bind is
  bracketed as a URL requires.
- An import failure names the condition it is: a distribution absent from the environment
  (`pip install`), one present but too old to hold the name asked of it
  (`pip install --upgrade`), or a standard-library module the running interpreter cannot
  provide -- reported with that interpreter's version, its executable, and the supported
  range, because no package supplies it.
- `cellucid serve` answers an operator mistake with one actionable line instead of a
  traceback, and keeps the traceback for genuine defects.
- Both writers refuse the same numeric input everywhere it can enter. Float32
  underflow (a nonzero value below `2^-149` becoming `0.0`) is now refused for vector
  fields here, and for quantized fields, embeddings, and `latent_space` in the R
  writer, where normalization hid the caller's magnitude from the write-time check.
- A label the browser cannot draw, including an unpaired surrogate, is refused by name.
- `prepare()` reconciles the export root, so it cannot publish coordinates the viewer is
  never told to fetch, and every manifest is checked as written bytes rather than only as
  an in-memory payload.
- Every published payload is little-endian whatever the host is.
- `prepare()` and the web cache both refuse a directory nobody set aside for them,
  including the working directory, the home directory, and paths reached through a link.
- A constant continuous field -- a gene at one level in every cell, or a column a subset
  flattened -- has an encoding instead of failing the export.
- `debug_connection()` reports the dataset probes for a published sample state.
- Publishing a sample's starting view restates its identity. A view is captured
  against a locally served copy of the generation, and the catalog serves the
  same bytes under its own id, so a capture published as-is was refused on every
  open and the sample showed none of its view. The publish step now renames the
  capture -- only after re-deriving the generation's counts and cell-order digest
  from what it ships -- and repoints the catalog digest itself, so the artifact
  and its advertisement cannot part company.
- A refusal that names two identical cell and gene counts as proof of a
  different dataset now names the two dataset identities instead.
- `prepare()` no longer rejects an export because of a gene it was never asked to export.
- A payload that is not entirely finite is refused with counts -- how many NaN,
  infinite, or outside float32, and the first positions -- instead of a bare `500`.
  The server answers `422` with that diagnosis as JSON; the exporter and the R writer
  raise the same report, and the viewer shows it.
- An official sample's published view applies over plain HTTP. Its integrity check used
  `crypto.subtle`, which is secure-context only, so a LAN bind opened the dataset
  uncoloured; the digest is computed the same way with or without WebCrypto.
- `Auto` level of detail follows the camera at every dataset size. Its floor was an
  absolute point budget, unsatisfiable below two million points, so on most datasets
  `Auto` answered full detail at every distance and did nothing; the floor is now the
  smaller of that budget and an 8x reduction cap.
- Level of detail no longer draws square patches. Sorting on a bit-reversed Morton
  code made each level a set of whole lattice cells; sorting on the plain code and
  bit-reversing the rank takes the same fraction from every neighbourhood. At 44x
  reduction: density error 138% to 5%, empty populated regions 42% to none.
- The velocity overlay's `Sync with LOD` read the level scale backwards, running at
  full particle count on the coarsest level and disposing itself at full detail.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts use
  authenticated schemas, bounded payloads, and path confinement.
- Both servers validate the `Host` header ahead of routing, file reads, and event
  delivery. A non-loopback bind is an explicit request for network exposure and has no
  authentication of any kind: every user who can reach the port reads the dataset.

### Documentation

- The Python, R, and web guides are organized around the current executable contracts,
  with real browser screenshots for primary workflows.
- Serving from an HPC cluster is documented end to end, including why a compute node
  commonly refuses SSH, the login-node port forward that works instead, why that requires
  `--host 0.0.0.0`, and why the browser address stays `127.0.0.1`.
- Corrected embedding guidance: an all-identical embedding has no finite normalization
  range and is rejected.
- Added a complete standard scVelo recipe for vector fields.

## [0.0.9] - 2026-01-01

This release graduates Cellucid Python out of alpha (still pre-1.0).
The GitHub release was tagged `v0.9.0`, while its Python distribution metadata
declared `0.0.9`. Release 0.9.1 restores one exact version across both systems.

### Added

- AnnData-first viewing and serving with `.h5ad` and `.zarr` support (lazy/backed loading where possible).
- Unified CLI serving with auto-detection for AnnData files, Zarr stores, and
  pre-exported dataset directories.
- Jupyter notebook integration (`CellucidViewer`, `AnnDataViewer`, `show()`, `show_anndata()`) with event hooks and session export.
- Session bundle support (`.cellucid-session`) including `CellucidSessionBundle` and `apply_cellucid_session_to_anndata()` for round-tripping highlights and user-defined fields back into AnnData.
- Multi-dimensional embedding exports (1D/2D/3D) and vector-field overlays (RNA velocity / drift) via `prepare()` and `vector_fields` helpers.
- Hosted web UI proxy mode with on-disk caching helpers (`get_web_cache_dir()`, `clear_web_cache()`).

### Changed

- Export format now includes explicit dataset identity metadata (`dataset_identity.json`) for reproducible sharing and session compatibility checks.
- Reduced export size and improved load performance with optimized manifests, improved connectivity edge export, and optional quantization + gzip compression knobs.

### Security

- Session bundles are treated as untrusted input with bounds checks and dataset mismatch policies when applying to AnnData.

### Documentation

- Major Read the Docs expansion and restructuring (Python package + web app guides), plus new publishing and contributing documentation.

## [0.0.1a0] - 2025

### Added

- Initial alpha release
- AnnData visualization with UMAP embeddings (1D, 2D, 3D)
- Gene expression overlays with sparse matrix support
- Cell metadata coloring (categorical and continuous)
- Interactive filtering and cell selection
- KNN connectivity visualization
- Multiple deployment modes:
  - Local demo (browser-only)
  - Browser file picker (h5ad and exported formats)
  - Server CLI (`cellucid serve`)
  - Python API (`serve()`, `serve_anndata()`)
  - Jupyter integration (`show()`, `show_anndata()`)
- Export functionality for web deployment

[0.9.6]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.6
[0.9.5]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.5
[0.9.4]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.4
[0.9.3]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.3
[0.9.2]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.2
[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
