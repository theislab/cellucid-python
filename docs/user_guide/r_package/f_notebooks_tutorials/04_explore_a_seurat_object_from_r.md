# Explore a Seurat Object from R

**Audience:** anyone with a Seurat object in an R session
**Time:** 10 minutes
**Goal:** look at the object in the Cellucid viewer, drive the viewer from the
console, and bring the cells you select back into R

You need the `cellucid` and `SeuratObject` packages, and RStudio, Positron or
any R console with a browser. Nothing is written to disk and no Python is
involved. The object is SeuratObject's own 80-cell example, `pbmc_small`; use
your own in its place.

## Step 1: view the object

```r
library(cellucid)
library(SeuratObject)
data("pbmc_small")

viewer <- cellucid_view(pbmc_small)
```

```text
Dataset id 'pbmc_small' (sessions are kept under it); name it with dataset_id = , or in misc$cellucid / metadata(x)$cellucid.
```

The viewer opens in the RStudio or Positron Viewer pane (elsewhere, in your
browser). The first time, the Cellucid web app is downloaded once (about 17
MB) and kept for later viewers. The 80 cells appear on the t-SNE, the PCA
beside it under **Data mapping**, and every `meta.data` column under
**Coloring & filtering** (`groups`, `RNA_snn_res.1`, …, and `active_ident`).

`viewer` describes itself:

```r
viewer
```

```text
<Cellucid viewer> pbmc_small (pbmc_small), running
  address:  http://127.0.0.1:8765/?viewer=true&viewerId=…&viewerToken=…&anndata=true
  data:     zarr, 80 cells
  channel:  connected
```

## Step 2: drive the viewer from the console

```r
cellucid_color_by(viewer, "groups")
cellucid_highlight(viewer, WhichCells(pbmc_small, idents = "1"), color = "#e41a1c")
```

The cells turn to the colours of `groups`, and the 25 cells of cluster 1 are
highlighted in red. Cells are given by name (as here), by 1-based position, or
as `TRUE`/`FALSE` per cell. A command sent before the page has connected waits
and goes out as soon as it does.

## Step 3: select cells in the viewer, get them in R

Ask R to tell you whenever you select:

```r
cellucid_on(viewer, "selection", function(event) {
  message(length(event$cells), " cells selected with the ", event$source)
})
```

Now, in the viewer's **Highlighting** panel, choose **Lasso**, hold
{kbd}`Alt` ({kbd}`Option` on a Mac), drag around a group of cells and click
**Confirm** ({doc}`../../web_app/f_highlighting_selection/03_highlight_ui`).
The console prints, for a lasso around 20 cells of `g2`:

```text
20 cells selected with the lasso
```

The last selection is always at hand:

```r
selected <- cellucid_selection(viewer)
head(selected$cell_names)
picked <- subset(pbmc_small, cells = selected$cell_names)
table(picked$groups)
```

```text
[1] "ATGCCAGAACGACT" "GAACCTGATGAACC" "TGACTGGATTCTCA" "AGTCAGACTGCACA"
[5] "AGGTCATGAGTGTC" "GGGTAACTCTAGTG"
g2 
20 
```

`selected$cells` holds the same cells as 1-based positions. Keep the
selection as a column of the object:

```r
pbmc_small$picked <- colnames(pbmc_small) %in% selected$cell_names
table(pbmc_small$picked)
```

```text
FALSE  TRUE 
   60    20 
```

The viewer shows the object as it was when it was opened; send it the new
column with `cellucid_push(viewer, pbmc_small$picked, key = "picked")`. Every
highlight, label and score made in the viewer comes back in one call,
`pbmc_small <- cellucid_pull(pbmc_small, viewer)`
({doc}`08_pull_everything_back`).

## Step 4: stop the viewer

```r
cellucid_stop(viewer)
```

This freezes an open page on its last view, as Python's `viewer.stop()` does,
then stops its server and frees its port (in a Jupyter notebook, where the
viewer is the cell's output and no server runs, the output stays, frozen);
`cellucid_stop_all()` stops every viewer of the session.

## Next

- {doc}`06_citeseq_seurat_to_mudata`: an object with RNA and ADT assays
- {doc}`../d_viewing_loading/06_hooks_selection_and_sessions`: every command,
  event and the session file
- {doc}`../d_viewing_loading/05_serve_from_hpc_and_ssh`: an object on a cluster
