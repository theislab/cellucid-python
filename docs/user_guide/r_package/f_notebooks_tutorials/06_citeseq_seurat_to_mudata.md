# CITE-seq: a Seurat Object with RNA and ADT

**Audience:** CITE-seq and multiome users with a Seurat object of several assays
**Time:** 10 minutes
**Goal:** view RNA and surface proteins side by side, switch between them from
R, and write the object as one Cellucid-ready MuData

The object is the small synthetic CITE-seq sample of
[theislab/cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets):
60 cells of three populations, an `RNA` assay of 8 genes and an `ADT` assay of
4 antibodies, a PCA, a UMAP and an SNN graph, saved with `saveRDS()`. Use your
own in its place.

## Step 1: the object

```r
library(cellucid)
library(SeuratObject)

download.file(
  "https://raw.githubusercontent.com/theislab/cellucid-demo-custom-datasets/HEAD/catalog/synthetic-citeseq-seurat.rds",
  "citeseq.rds", mode = "wb"
)
cite <- readRDS("citeseq.rds")
cite
Assays(cite)
```

```text
An object of class Seurat 
12 features across 60 samples within 2 assays 
Active assay: RNA (8 features, 0 variable features)
 2 layers present: counts, data
 1 other assay present: ADT
 2 dimensional reductions calculated: pca, umap
[1] "RNA" "ADT"
```

## Step 2: view both assays

```r
viewer <- cellucid_view(cite)
```

Each assay is a **modality**: the cells, their metadata, embeddings and graph
are shared, and the genes are one modality's at a time. The viewer opens on
the active assay, `RNA`, on the UMAP. Under **Session › Data mapping ›
Genes**, a **Modality** select lists `RNA` (8 features) and `ADT` (4).

Switch to the proteins from R, and colour by cell type:

```r
cellucid_set_data_mapping(viewer, modality = "ADT")
cellucid_color_by(viewer, "cell_type")
```

`cellucid_set_data_mapping()` returns once the viewer has applied the change
(it stops with the viewer's reason when the viewer refuses it). Gene search,
colouring by a feature, marker genes and the dot plot now read the proteins:
search `CD19-ADT` to see the B cells. `cellucid_stop(viewer)` when you are
done.

## Step 3: write one file for others

```r
cellucid_prepare("citeseq.zarr", format = "zarr", data = cite)
```

An object of several assays is written as a **Cellucid-ready MuData**: one
modality per assay, the active assay first, the shared cells, embeddings and
graph at the root. This object names itself (`cite@misc$cellucid`), so no
`dataset_id` is needed. `modalities = c("RNA", "ADT")` chooses and orders the
assays; one name writes a plain AnnData of that assay.

The folder is MuData's own Zarr layout. In Python:

```python
import mudata as md
md.read_zarr("citeseq.zarr")
```

```text
MuData object with n_obs × n_vars = 60 × 12
  obs:	'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'nCount_ADT', 'nFeature_ADT', 'cell_type', 'sample', 'active_ident'
  uns:	'title', 'cellucid'
  obsm:	'X_pca', 'X_umap'
  obsp:	'connectivities'
  2 modalities
    RNA:	60 × 8
    ADT:	60 × 4
```

Publish the folder like any store
({doc}`05_share_a_singlecellexperiment_on_github`); the web app shows the same
**Modality** select. The `.rds` itself opens in the web app too, with no R
({doc}`../../web_app/b_data_loading/17_open_r_objects`): the demo repository
lists it in its catalog.

## Next

- {doc}`../d_viewing_loading/07_multimodal_assays`: which assays are written
  and how
- {ref}`r-cellucid-ready-mudata`: the exact layout
