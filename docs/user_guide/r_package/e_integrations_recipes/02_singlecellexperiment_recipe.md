# SingleCellExperiment Recipe

**Audience:** Bioconductor users
**Time:** 5 minutes
**What you'll do:** view a SingleCellExperiment (or SummarizedExperiment) in the
Cellucid viewer from R, then write it as a Cellucid-ready AnnData (or MuData).

```{note}
`cellucid` reads SingleCellExperiment and SummarizedExperiment objects
natively: pass the object as `data`. `SingleCellExperiment` (with
`SummarizedExperiment`) must be installed; `cellucid` itself does not depend on
them.
```

## Prerequisites

- A SingleCellExperiment `sce` (SingleCellExperiment 1.8 or later) with at least
  one reduced dimension, for example from `scater::runUMAP()`.
- R packages: `cellucid`, `SingleCellExperiment`.

## View it, then write it

```r
library(cellucid)

viewer <- cellucid_view(sce)
cellucid_stop(viewer)

store <- cellucid_prepare(
  "sce_dataset.zarr", "zarr",
  data = sce,
  dataset_id = "sce-dataset",
  dataset_name = "My SCE dataset",
  latent_key = "X_PCA"
)
```

## What each part becomes

| SingleCellExperiment | In the store |
| --- | --- |
| `colnames(sce)`, `colData(sce)` | the `obs` index and columns (`Rle` columns expanded); size factors as `sizeFactor` (or `size_factor`) |
| `reducedDims(sce)` | `obsm['X_<name>']`, the name's case kept (`X_UMAP`, `X_PCA`) |
| `logcounts` (else `counts`, else the first assay) | `X`, gene-major float32 CSC (one matrix per experiment: the other assays are not written) |
| `rownames(sce)`, `rowData(sce)` | the `var` index and columns |
| `altExps(sce)` | one modality each, the main experiment (`mainExpName()`, else `"main"`) first: a Cellucid-ready MuData |
| `colPairs(sce)` | `obsp['connectivities']`: the pair named `connectivities`, then one named like `snn`, then any other: the first that is symmetric, its diagonal removed |
| `metadata(sce)$cellucid` | the dataset's identity and launch defaults |

Assays backed by HDF5Array or another DelayedArray backend are read into memory
block by block, each block a set of whole cells. A **SummarizedExperiment** has no reduced
dimensions: pass `embeddings = list(X_umap = coordinates)`.

## Troubleshooting

- *"This SingleCellExperiment holds no embedding"* — `sce <- scater::runUMAP(sce)`,
  or pass `embeddings`.
- *"This SingleCellExperiment was made by SingleCellExperiment 1.6"* —
  `sce <- updateObject(sce)`.
- *"… has unnamed assays and was not read"* — name them with `assayNames(sce) <- …`.
- More: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
