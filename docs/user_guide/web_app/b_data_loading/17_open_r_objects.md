# Open R objects and MatrixMarket folders

**Audience:** anyone with a Seurat or SingleCellExperiment object saved from
R, or a count matrix in MatrixMarket form  
**Time:** 5 minutes  
**What you’ll learn:** how to open an `.rds`, `.rda` or `.RData` file and a
MatrixMarket folder directly, what Cellucid reads from each and what it
leaves out, how long a large object takes, and what each refusal asks you to
do  
**Prerequisites:** none

Cellucid opens these files as they are: there is no "convert to AnnData
first" step. Unlike an `.h5ad`, a compressed R file or a gzip text file has
no random access, so it **cannot be read in parts**: it is read whole, once,
from the start into this browser's storage (the origin-private file system),
the expression matrix laid out by gene and the cells in priority order as a
Cellucid-ready AnnData's are (so a subsample opened later reads the first rows
of every array), and
every view afterwards reads what was kept there, as it reads any AnnData
store ({doc}`07_file_format_expectations`). The opening says so before it
reads, with how much it reads. What is kept stays for your next visit when
the file has a version (a file on this computer always has one; a linked
file when its host reports one): the same file opened again, in this tab or
another, today or next week, opens from there in seconds and reads nothing of
the file ({ref}`which-formats-are-read-in-parts`). A file known only by its
size is kept for its tab only.

---

## Seurat and SingleCellExperiment objects (`.rds`, `.rda`, `.RData`)

**Open… › File…**, drop the file on the viewer, or paste a link to it
({doc}`14_open_from_link`; a linked R file is read as it downloads, and
nothing is saved to your Downloads folder). A gzipped copy (`.rds.gz`,
`.rda.gz`, `.RData.gz`, as GEO stores R objects) opens the same way: the
file's content, not its name, says how many times it is compressed, and it
is decompressed as it streams. A file dialog that matches only a file's last
extension may grey out an `.rds.gz`; drop it on the viewer instead.

### Try it with SeuratObject's own sample

In R:

```r
library(SeuratObject)
saveRDS(pbmc_small, "pbmc_small.rds")
```

Open `pbmc_small.rds`. Expected: **80 cells**, **230 genes**, the embeddings
`X_pca` and `X_tsne`, the obs fields of `pbmc_small[[]]` (`orig.ident`,
`nCount_RNA`, `RNA_snn_res.0.8`, `letter.idents`, `groups`, …) plus
`active_ident`, and the `RNA_snn` neighbour graph. The Dataset box says
*Read whole once (it cannot be read in parts), kept in this browser's storage
for next time · 5 parts not read*:
`pbmc_small` carries a `scale.data` matrix, PCA loadings and JackStraw
results, and its `commands`, which Cellucid leaves out
([below](dataset-box-conversion)).

A SingleCellExperiment saved the same way (`saveRDS(sce, "sce.rds")`) opens
the same way.

### What is read

| Seurat (SeuratObject 3.0 and later) | SingleCellExperiment (1.8 and later), SummarizedExperiment | In Cellucid |
|---|---|---|
| `Cells(obj)`, `obj[[]]` (`meta.data`) | `colnames(sce)`, `colData(sce)` | the cells and their obs fields |
| `Idents(obj)` | | the obs field `active_ident` |
| `Embeddings(obj, "<r>")` for every reduction | `reducedDim(sce, "<name>")` | the embedding `X_<name>` (`X_umap`, `X_PCA`: the name's case is kept) |
| the default assay's `data` layer (else `counts`) | `logcounts` (else `counts`) | the expression matrix `X`; the other is a layer |
| `meta.features` / the Assay5 feature table, `VariableFeatures()` | `rowData(sce)`, row names | the gene table; `highly_variable` marks the variable features |
| other assays (`ADT`, `peaks`, `SCT`) | `altExps(sce)` | modalities ({doc}`19_multimodal_mudata`), the default assay first |
| `obj@graphs` (`RNA_snn`, `RNA_nn`) | `colPairs(sce)` (weights: the first metadata column, as `colPair(sce, asSparse = TRUE)` reads them) | neighbour graphs |
| Visium and FOV images | | the cells' positions as the embedding `spatial` |
| `obj@misc$cellucid` | `metadata(sce)$cellucid` | what the dataset says about itself (`uns['cellucid']`) |
| `obj@project.name` | | the dataset's title (unless `SeuratProject`) |

Seurat v5 objects with split layers (`counts.1`, `data.ctrl`, made by
`split()` or `merge()`) are joined as `JoinLayers()` would. Factors keep
their level order and `NA` stays missing; logical and integer columns with
`NA` stay nullable.

**Left out**, and never downloaded twice: `scale.data` (dense, the variable
genes only), `commands`, `tools`, `neighbors`, feature loadings, JackStraw
results and fitted models. The file still has to be read past them once, but
nothing of them is kept, not even a Seurat v5 `scale.data` layer, whose name
R writes only after its values.
Each one that holds something is listed in the Dataset box's Details under
**Not read**, with its size in the file; the empty `scale.data`, `commands`
and `neighbors` every Seurat object carries are not.

Cellucid reads R's own serialization (versions 2 and 3, gzip, bzip2 or xz,
or uncompressed) with its own parser in the browser: R is not needed. An
**uncompressed** file (`saveRDS(obj, file, compress = FALSE)`, or
`readr::write_rds()`, which does not compress) is not copied at all: its
structure is read and every matrix is then read in place from the file, as
from an `.h5ad`. Its cells therefore keep the object's order (a converted
file's are put in priority order), and a subsample reads it in blocks, as any
file read in place.

### A file with several objects

A `.rda` or `.RData` file that holds more than one Seurat object or
SingleCellExperiment (`save(seu_all, seu_tcells, file = "objects.rda")`) is
read once, keeping every one of them, and the opening asks which to open:

> objects.rda holds 2 single-cell objects. Open which one?  
> seu_all (Seurat, 52,113 cells) · seu_tcells (Seurat, 8,201 cells)

```{figure} ../../../_static/screenshots/data_loading/rda-choose-object.png
:alt: The opening panel over the viewer reading "Opening pancreas-objects.rda" and "pancreas-objects.rda holds 2 single-cell objects. Open which one?", with one button per Seurat object, "pancreas (Seurat, 3,696 cells)" and "pancreas_ductal (Seurat, 916 cells)", and Cancel.
:width: 824px

An `.rda` holding two Seurat objects, opened from this computer: the file is
read once, then the opening asks which object to open, with each one's class
and cell count. **Cancel** keeps the dataset that was open.
```

Opening the same file again and choosing the other one reads nothing
again. Other objects in the file (data frames, vectors) are passed
over.

### How long a large object takes

The first opening runs at the speed of the slower of the link and the
decompression; every later one opens what was kept. Measured on an Apple M1
Pro (2026-10-27) with Seurat v5 objects (counts and data, PCA, UMAP, an SNN
graph, meta.data; `scripts/rds-corpus/make_bench.R`):

| Object | First opening | Later visits | A Cellucid-ready AnnData of it |
|---|---|---|---|
| 50,000 cells, 1,500 genes per cell, with `scale.data`; 1.2 GB gzip `.rds` (2.65 GB decompressed) | linked over 25 MB/s and 40 ms: 56 s to the first view, 1.2 GB downloaded, 1.6 GB peak browser memory; from disk (Node): 11 s | 1.2 s to the first view, nothing downloaded; a gene in 0.3 s | 244 MB `.zarr.zip`: 1.3 s, 0.7 MB to the first view; a gene in 0.3 s, 0.2 MB |
| 300,000 cells, 500 genes per cell; 1.05 GB gzip `.rds` | linked over 25 MB/s: 55 s, 1.05 GB, 2.1 GB peak; from disk (Node): 19 s | 1.2 s, nothing downloaded | |
| the 50,000-cell object as xz / bzip2 (0.98 GB / 0.99 GB) | from disk (Node): about 200 s / 115 s (the decompression itself, about 13 / 23 MB/s) | 1.2 s | |

A Seurat object stores its assays first and its reductions near the end
(the 50,000-cell object's embedding starts after 98 % of the file), so no
view can be drawn before the whole file is read. The pass holds a bounded
window in memory, whatever the object's size (under 700 MB in Node for
either object). While it reads, the opening shows how much is read, how
fast and the time left, with **Cancel**; a cancelled pass leaves nothing
behind. To share a large object with others, write it once as a
Cellucid-ready AnnData with `cellucid::cellucid_prepare(path, data = obj)`
(R) or `cellucid.prepare()` (Python): that is read in place, so a viewer
reads only what a view needs.

```{figure} ../../../_static/screenshots/data_loading/opening-progress.png
:alt: The opening panel centred over the empty viewer while a linked Seurat .rds is downloaded: "Opening" and the file's address, the step "Downloading the whole file (9.4 MB): … does not serve parts of files", a progress bar part-way, how much has arrived with the speed and time left, the note saying the host serves the file only whole, how much is downloaded and that it is kept, and a Cancel button.
:width: 856px

A linked Seurat `.rds` opening while no other dataset is on screen: the panel
over the empty viewer names the file, the step under way (here a host that
sends whole files only) with its size, how much has arrived, the speed and the
time left, then, before anything else, why the file is read whole, how much is
read, that it is kept for next time and what is read in parts instead, and
**Cancel**. With a dataset on screen the same lines are a
notification, and that dataset stays until the new one is ready.
```

---

## MatrixMarket folders

**Open… › Folder…** (or drop the folder) on a folder holding a MatrixMarket
count matrix and the two lists that name its rows and columns:

| Written by | Matrix | Features | Cells |
|---|---|---|---|
| Cell Ranger 3 and later | `matrix.mtx.gz` | `features.tsv.gz` | `barcodes.tsv.gz` |
| Cell Ranger 2, STARsolo | `matrix.mtx` | `genes.tsv` or `features.tsv` | `barcodes.tsv` |
| kb-python | `cells_x_genes.mtx` | `cells_x_genes.genes.txt` | `cells_x_genes.barcodes.txt` |
| alevin-fry | `quants_mat.mtx` | `quants_mat_cols.txt` | `quants_mat_rows.txt` |
| a GEO upload | `GSM…_matrix.mtx.gz` | `GSM…_features.tsv.gz` or `_genes` | `GSM…_barcodes.tsv.gz` |

A sample prefix (`GSM5008737_RNA_3P-`, ending in `_`, `-` or `.`) pairs
each matrix with its own lists, so a folder of a GEO series' extracted
`_RAW.tar`, every sample side by side, opens too: when it holds several
samples Cellucid asks which one to open. Lists under a shorter shared prefix
serve every matrix below it (one `GSE1_barcodes.tsv.gz` for
`GSE1_RNA_matrix.mtx.gz` and `GSE1_ADT_matrix.mtx.gz`), but never a matrix
that has lists of its own. A sample whose barcode or feature list is missing,
or that has both a `features` and a `genes` list, is refused naming the
files. A GEO sample can also be opened from its link without downloading it
({doc}`14_open_from_link`).

A list's first line is a header only when it is a column name and nothing
else (`barcode`, `cell_id`, `gene`, `gene_id`, `feature`, `id`, Parse's
`bc_wells`): a first cell named `cell0` or a first gene named `IDH1` is a
cell or a gene.

The matrix may be genes × cells or cells × genes; which it is is read from
its size against the two lists. A filtered matrix is preferred over a raw
one in the same folder. The gene table has `gene_ids` and, for Cell Ranger 3
and later, `feature_types`. Features of several types (`features.tsv`
column 3: Gene Expression with Peaks in a Multiome sample, or with Antibody
Capture) open as one modality per type, Gene Expression first, as a 10x
`.h5` does ({doc}`19_multimodal_mudata`). A matrix written by cell (Cell
Ranger, STARsolo) or by feature (kb-python, alevin-fry) is read as it was
written and one in any other order is sorted in blocks, so the browser's
memory never has to hold the whole matrix; the store keeps it by gene, its
cells in priority order.

A count matrix has **no cell embedding**, and Cellucid draws cells at
embedding coordinates. The folder's own cell tables give it one:

- Cell Ranger's `analysis/` folder: `analysis/umap/*/projection.csv`,
  `tsne`, `pca` (as `X_umap`, `X_tsne`, `X_pca`) and
  `analysis/clustering/*/clusters.csv` (as obs fields named after the
  clustering);
- any `.csv` or `.tsv` table of the cells in the folder (`cell_metadata.csv`,
  a GEO `metadata.tsv`): its cell column is matched to the barcodes, columns
  named like `UMAP_1`/`UMAP_2` or `tSNE-1`/`tSNE-2` become embeddings and the
  others obs fields. A column named for a clustering (`Cluster`,
  `seurat_clusters`, `leiden`, `louvain`, `kmeans_10`) is categorical even
  when its labels are numbers, as Cell Ranger's `clusters.csv` is.

Without either the folder is refused with the lines that compute an
embedding in Python and in R. A Cell Ranger `outs/` folder that also holds
`filtered_feature_bc_matrix.h5` opens through that file instead, read in
place ({doc}`18_open_loom_and_10x`).

```{figure} ../../../_static/screenshots/data_loading/matrix-market-dataset-box.png
:alt: The Dataset box after a MatrixMarket folder with a cell table was opened with Open… › Folder…: the folder's name, its cells and genes, "MatrixMarket folder · this computer", the line saying where the converted data is kept, and Details open.
:width: 488px

A Cell Ranger 3 MatrixMarket folder with a `cell_metadata.csv` beside it,
opened with **Open… › Folder…**: the table's `UMAP_1`/`UMAP_2` columns are the
embedding and its other columns obs fields. Like an R object, the folder is
read whole once and kept in this browser's storage for next time.
```

---

(dataset-box-conversion)=
## What the Dataset box says

Under the line that names the format and where it is read from
(*R object (.rds) · this computer*, *MatrixMarket folder · this computer*),
one more line says where the converted data is kept and how many parts of
the file were not read:

| The line | When |
|---|---|
| Read whole once (it cannot be read in parts), kept in this browser's storage for next time | the first opening of the file |
| Opened from this browser's storage, read whole on an earlier open | the same file opened again (a file with a version: in any tab, on any later visit; one known only by its size: in its own tab): nothing of the file is read |
| Read whole once (it cannot be read in parts), into this browser's storage for this tab | a file known only by its size (a link whose host reports no version): kept until the tab closes. A viewer host's file (`cellucid serve`, a notebook, `cellucid_view()`) is listed with its `ETag`, its version: kept for next time, and read again once it is rewritten |
| Structure read once (arrays read in place), kept in this browser's storage for next time | an uncompressed R file (`saveRDS(compress = FALSE)`) |
| Read whole once (it cannot be read in parts), into memory for this tab | the browser gives the page no storage (a Firefox or Safari private window, a page not served over HTTPS): the converted data is held in memory instead, and is read again on each opening |

followed by *· 5 parts not read* when something was left out. The
**Details** disclosure says the rest:

- **Read**: *Whole, once, from start to end (it cannot be read in parts),
  from this computer; nothing is uploaded*, or *… from* the host of a link;
- **Kept**: where the converted data stays (this browser's storage, the
  origin-private file system, until it is cleared or the browser needs the
  room; or memory, until the tab closes);
- **Read in parts**: what is read in parts instead, a Cellucid-ready AnnData
  written with `cellucid_prepare()` (R) or `cellucid.prepare()` (Python);
- **Not read**: one line per part left out and why, named by its place in the
  R object, for example
  *assays$RNA@scale.data (15 KB): Cellucid does not use it.*,
  *assays$ADT@counts is kept outside the file (MatrixDir at …) and was not
  read.*, or a cell table of a MatrixMarket folder that could not be read.

An `.h5ad`, an AnnData Zarr store and the other formats read in place have
none of these: nothing of them is copied.

**Session › Open data › Stored** lists how many files read whole are kept in
this browser and how much they take, with **Clear**, which removes every one
no tab has open. When a new file does not fit in the storage the browser
grants the site, kept gene indexes are removed first, then the files opened
least recently; a gene index being kept never removes a converted file.

Everything kept is written with a checksum of every 16 KB of it (CRC-32), and
every later read of it is checked against them: converted data, and a file
downloaded whole from a host that serves no parts of files. Damaged bytes are
never shown: the read fails saying the kept copy is damaged, the copy is
removed, and opening the file again reads it again.

```{figure} ../../../_static/screenshots/data_loading/rds-dataset-box.png
:alt: The Dataset box after a Seurat object saved with saveRDS() was opened: its cells and genes, "R object (.rds) · this computer", the line "Read whole once (it cannot be read in parts), kept in this browser's storage for next time", and Details open with Read (whole, once, from start to end, nothing uploaded), Kept (the origin-private file system until it is cleared), Read in parts (a Cellucid-ready AnnData) and what it holds.
:width: 488px

A Seurat object saved with `saveRDS()` and opened from this computer, with
**Details** open: **Read** says it was read whole, once, from start to end and
nothing was uploaded, **Kept** that the converted data stays in this browser's
storage until it is cleared, and **Read in parts** what is read in parts
instead.
This object carries nothing Cellucid leaves out, so no *parts not read* follow.
```

---

## When a file is refused

Each refusal says what the file is and what to do; the sections of
{doc}`08_troubleshooting_data_loading` say more.

| Message says | Why | What to do |
|---|---|---|
| … keeps the RNA assay's counts outside the file (MatrixDir at …) | a BPCells or HDF5Array matrix: the values are not in the `.rds` | in R, `seu[["RNA"]]$counts <- as(seu[["RNA"]]$counts, "dgCMatrix")` (or `assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")`), then `saveRDS()`; or `cellucid::cellucid_view(obj)` |
| … has no layers in the file (SaveSeuratRds() moved them …) | `SaveSeuratRds(move = TRUE)` keeps the layers beside the file | load it with `LoadSeuratRds()` and save the object with its layers as `dgCMatrix` |
| This is a Seurat v2 object | Seurat 2 (2017–2018) classes | `UpdateSeuratObject()`, then `saveRDS()` |
| This is an h5Seurat file | SeuratDisk's own format, not read | `SeuratDisk::LoadH5Seurat()`, then `saveRDS()` |
| … has no cell embedding | no reduction, no `reducedDims` | run PCA and UMAP (`RunPCA()`, `RunUMAP()`; `scater::runUMAP()`) and save again |
| … compressed with zstd | `saveRDS(compress = "zstd")` (R 4.5 and later): browsers have no zstd decompressor | save with the default gzip, or `"xz"` / `"bzip2"` |
| … saved as text / in R's native binary form | `saveRDS(ascii = TRUE)`, `xdr = FALSE` | save with the defaults |
| … holds an R data frame, not a Seurat object or a SingleCellExperiment | the file holds some other R value | save the single-cell object itself |
| This MatrixMarket matrix is symmetric / complex / a dense array | not a count matrix | write the counts as a coordinate matrix |

`.qs` and `.qs2` files (the `qs` packages) are not R serialization: save the
object with `saveRDS(obj, file)`.

## Related pages

- {doc}`07_file_format_expectations`: every format Cellucid reads
- {doc}`18_open_loom_and_10x`: Loom, 10x `.h5` and `outs/` folders, read in place
- {doc}`19_multimodal_mudata`: modalities (Seurat assays, `altExps`)
- {doc}`../../r_package/index`: `cellucid_prepare()` from R
