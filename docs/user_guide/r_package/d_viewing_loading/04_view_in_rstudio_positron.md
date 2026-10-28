# View an Object in RStudio or Positron

**Audience:** anyone with a Seurat or SingleCellExperiment object in R
**Time:** 2 minutes
**What you'll do:** open the object in the Cellucid viewer in the IDE's Viewer
pane, with nothing written to disk and no Python.

## The one line

```r
library(cellucid)
library(SeuratObject)
data("pbmc_small")

viewer <- cellucid_view(pbmc_small)
#> Dataset id 'pbmc_small' (sessions are kept under it); name it with dataset_id = , or in misc$cellucid / metadata(x)$cellucid.
#> Opened pbmc_small in the Viewer pane. If it stays blank (no WebGL in the pane), view it with pane = "browser".
```

`cellucid_view()` takes a Seurat, SingleCellExperiment or SummarizedExperiment
object, or the path of any file the web app opens (an `.h5ad`, a `.zarr` folder
or `.zarr.zip`, a MuData, a Loom or 10x `.h5`, an `.rds`). In RStudio and
Positron it opens in the Viewer pane; in a Jupyter notebook in the cell's
output ({doc}`08_view_in_jupyter_notebooks`); elsewhere in your browser.
Choose with `pane = "viewer"`, `"browser"`, `"notebook"` or `"none"`.

The 80 cells of `pbmc_small` appear on their t-SNE, with the PCA and the
metadata columns (`groups`, `RNA_snn_res.1`, `active_ident`, …) under
**Color by**, and the 230 genes searchable.

## What happens

In RStudio, Positron or an R console (in a Jupyter notebook no server starts:
see {doc}`08_view_in_jupyter_notebooks`):

1. The object is read as `cellucid_prepare(data = pbmc_small)` reads it
   ({doc}`../e_integrations_recipes/01_seurat_recipe`), into the same
   Cellucid-ready AnnData, held in memory. No chunk is made yet.
2. A small web server starts in the background (port 8765, the one local
   address Google sign-in accepts; the next free port when it is taken).
3. The first time, the Cellucid web app is downloaded from www.cellucid.com
   (about 17 MB) into `file.path(tools::R_user_dir("cellucid", "cache"), "web")`, every file
   checked against the build's SHA-256 inventory. Later viewers reuse it, and
   download again only when a new build is published.
4. The viewer reads the store with range requests; each chunk is made when it
   is asked for, uncompressed (the store declares no compressor: the browser
   reads it over loopback, where deflating on R's one thread cost far more than
   the bytes it saved; `cellucid_prepare()` files stay compressed). An object's
   genes x cells matrix is never transposed whole: reading a gene makes the
   block of genes around it, one binary search a cell. A 100,000-cell object
   opens in under a second; on a 1,000,000-cell object (150 values a cell) the
   first gene takes about 0.3 s and 100 more genes about 0.13 s.

The console stays usable while the viewer is open: the server answers the
viewer whenever R is idle. Run `cellucid_stop(viewer)` when you are done, or
`cellucid_stop_all()`.

```{note}
The object is read when the viewer starts. Changing it afterwards does not
change what the viewer shows; view it again.
```

## Choose what is shown

`cellucid_view()` takes the arguments of `cellucid_prepare()` that choose what
is shown:

```r
viewer <- cellucid_view(
  seu,
  dataset_id = "pbmc-v1", dataset_name = "PBMC",
  latent_key = "X_pca",                    # the outlier filter's space
  embedding_defaults = c("2d" = "X_umap"),
  obs_keys = c("seurat_clusters", "orig.ident"),
  modalities = c("RNA", "ADT")             # several assays: see Multimodal assays
)
```

Sessions and community annotations are kept under the dataset id, so give a
dataset you come back to a stable `dataset_id`.

Without one, an object is named by its own record (`misc$cellucid` of a Seurat
object, `metadata(x)$cellucid` of a SingleCellExperiment), else after the
expression you passed (`seu`). A file is named by its own record only when it
is a Zarr store (`uns['cellucid']`); an `.h5ad` or `.h5mu` is named after the
file, because the R package reads no HDF5. Python's `serve_anndata()` and
`show_anndata()` read that record from an `.h5ad` and `.h5mu` too, and refuse a
file that names no dataset, so the same file can be kept under two ids: pass
`dataset_id =` (and `dataset_name =`) to keep its sessions in one place.

## A file given by path

A path is served from disk as it is, and the browser reads the file itself.
The R package has no HDF5 reader (no `rhdf5` or `hdf5r` among its
dependencies) and no reader of Zarr chunks, so it cannot put a file's cells
in **priority order** the way it serves an object, or Python serves a file.
When the file's cells are not in priority order already (a Zarr store
`cellucid_prepare()` wrote is; an `.h5ad` or `.h5mu` may be, which R cannot
tell), `cellucid_view()` and `cellucid_serve()` say so as they start: a
subsample the viewer opens (Session › *Cells to open*) then reads the
storage blocks that hold its cells, with the cells around them, as it does
from a link. For a subsample that reads only its own cells:

- read the file into R and view the object — `Seurat::Read10X_h5()`,
  `SeuratDisk::LoadLoom()`, `zellkonverter::readH5AD()`,
  `MuDataSeurat::ReadH5MU()` — which is served in priority order;
- write it once in priority order with `cellucid_prepare()`, and serve that;
- or serve the file from Python (`cellucid serve <path>`, `show_anndata()`),
  which serves `.h5ad`, `.h5mu`, Loom, Cell Ranger `.h5`/`outs/` and AnnData
  and MuData Zarr stores in priority order, read on the kernel.

An `.rds`, `.rda` or MatrixMarket folder is read whole once by the browser,
into a store in priority order, and needs none of this.

The viewer is told the file's size and its version — a strong `ETag` made of
its size, modification time and status-change time, read when the viewer
asks — and keeps what it reads of the file under that version: its byte
cache, its gene indexes and an `.rds` converted once last across visits, so
opening the same file again reads none of it. A file rewritten while it is
served, in place or replaced, is a new version and is read again.

## Next

- {doc}`06_hooks_selection_and_sessions`: highlight cells from R, react to a selection, save the session.
- {doc}`05_serve_from_hpc_and_ssh`: a viewer for an object on a cluster.
- {doc}`07_multimodal_assays`: RNA and protein (or chromatin) side by side.
