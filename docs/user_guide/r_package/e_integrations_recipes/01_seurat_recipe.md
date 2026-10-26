# Seurat Recipe

**Audience:** Seurat users (wet lab + computational)
**Time:** 5 minutes
**What you'll do:** view a Seurat object in the Cellucid viewer from R, then
write it as a Cellucid-ready AnnData (or MuData) to share.

```{note}
`cellucid` reads Seurat objects natively: pass the object as `data`.
`SeuratObject` must be installed (it is when Seurat is); `cellucid` itself does
not depend on it.
```

## Prerequisites

- A Seurat object `seu` (Seurat 3 or later; SeuratObject 4 or 5) with at least
  one reduction to plot (UMAP is typical).
- R packages: `cellucid`, `SeuratObject`.

## View it, then write it

```r
library(cellucid)

# seu <- readRDS("my_seurat_object.rds")

viewer <- cellucid_view(seu)          # the RStudio / Positron Viewer pane, else the browser
cellucid_stop(viewer)

store <- cellucid_prepare(
  "seurat_dataset.zarr", "zarr",
  data = seu,
  dataset_id = "seurat-dataset",
  dataset_name = "My Seurat dataset",
  latent_key = "X_pca",
  obs_keys = c("seurat_clusters", "orig.ident", "nCount_RNA")
)
```

Open `seurat_dataset.zarr` on a web host, or write `seurat_dataset.zarr.zip`
with `format = "zarr.zip"` and open it with **Open… › File…**
({doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`).

## What each part becomes

The mapping is the web app's for the same object saved with `saveRDS()`, so a
session saved on the `.rds` in the browser applies to the prepared store and to
`cellucid_view(seu)`.

| Seurat | In the store |
| --- | --- |
| `seu[[]]` (`meta.data`) | `obs`; factors stay categoricals with their level order; `Idents(seu)` as `active_ident` |
| every reduction (`Embeddings(seu, "umap")`, …) | `obsm['X_<name>']` (`X_umap`, `X_pca`, `X_tsne`), every column, float32 |
| the active assay's `data` layer (else `counts`) | `X`, gene-major float32 CSC (one matrix per assay: `counts` beside `data` is not written) |
| `meta.features` (Assay5 `meta.data`), `VariableFeatures()` | `var`, with `highly_variable` |
| `seu@graphs[["<active assay>_snn"]]`, its diagonal removed | `obsp['connectivities']`: the first symmetric graph, `<active assay>_snn`, then any other `_snn` graph, then any other |
| `seu@misc$cellucid` | the dataset's identity and launch defaults (`uns['cellucid']`) |
| `project.name` (unless `"SeuratProject"`) | the dataset name when none is given |
| several assays (RNA, ADT, …) | one modality each: a Cellucid-ready MuData, the active assay first |

Seurat v5 split layers (`counts.1`, `data.ctrl`, after `split()`) are joined on
a copy; your object is not changed. `scale.data` is not read: it is dense, holds
only the variable genes, and is derived from `data`. Layers kept on disk by
BPCells are read into memory block by block.

## Choices and variants

- **Name the dataset once, in the object.**
  `seu@misc$cellucid <- list(id = "pbmc-v1", name = "PBMC", citation = list(text = "…", url = "https://…"), defaults = list(embedding2d = "X_umap", latentKey = "X_pca"))`;
  `cellucid_prepare()` and `cellucid_view()` then need no `dataset_id` or
  `dataset_name`.
- **Replace one part.** Any argument replaces what the object holds:
  `obs = seu[[c("seurat_clusters", "orig.ident")]]` stores only those columns;
  `embeddings = list(X_umap = Embeddings(seu, "umap"))` only that embedding.
- **One assay of several.** `modalities = "RNA"` writes a plain AnnData;
  `modalities = c("SCT", "ADT")` chooses and orders the modalities, and
  `modality = "ADT"` opens the viewer on one. See
  {doc}`../d_viewing_loading/07_multimodal_assays`.
- **Only some genes.** `gene_identifiers = c("CD3E", "MS4A1", …)` (one assay).
- **A 3D UMAP.** `RunUMAP(seu, dims = 1:30, n.components = 3, reduction.name = "umap3d")`
  is stored as `X_umap3d`.

## Troubleshooting

- *"This Seurat holds no embedding (no reductions)"* — run `seu <- Seurat::RunUMAP(seu, dims = 1:30)`
  first, or pass `embeddings = list(X_umap = coordinates)`.
- *"This Seurat object was made by Seurat 2.x"* — `seu <- SeuratObject::UpdateSeuratObject(seu)`.
- *"This Seurat object has Seurat v5 assays … which SeuratObject 4.x cannot read"* —
  `install.packages("SeuratObject")`.
- *"dataset_id and dataset_name are required"* — pass both, or set
  `seu@misc$cellucid`.
- A message *"From the Seurat object: … not read"* names each part left out and
  why (a metadata column of a type AnnData does not store, a reduction that does
  not cover every cell, a graph that is not symmetric).
- More: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
