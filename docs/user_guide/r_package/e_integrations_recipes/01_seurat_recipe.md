# Seurat Recipe

**Audience:** Seurat users (wet lab + computational)
**Time:** 10–20 minutes
**What you'll do:** extract the inputs from a Seurat v5 object → write a
Cellucid-ready AnnData store → open it in the Cellucid web app.

```{note}
`cellucid-r` does not depend on Seurat. This page shows how to *extract* the
inputs from a Seurat object; the code needs `Seurat` (or `SeuratObject`) and
`Matrix` installed.
```

## Prerequisites

- A Seurat object `seu` with at least one reduction to plot (UMAP is typical).
- R packages: `cellucid`, `Seurat`, `Matrix`.

## The whole recipe

```r
library(cellucid)
library(Seurat)

# seu <- readRDS("my_seurat_object.rds")

# Expression: Seurat keeps genes x cells; cellucid_prepare() takes cells x genes.
# A Seurat v5 assay split into layers needs JoinLayers() first.
expr <- Matrix::t(GetAssayData(seu, assay = "RNA", layer = "data"))

# The shared-nearest-neighbour graph links every cell to itself; the stored
# graph must have a zero diagonal.
graph <- seu@graphs[["RNA_snn"]]
Matrix::diag(graph) <- 0

store <- cellucid_prepare(
  out_path = "seurat_dataset.zarr",
  format = "zarr",
  embeddings = list(
    X_umap = Embeddings(seu, reduction = "umap"),
    X_pca = Embeddings(seu, reduction = "pca")
  ),
  dataset_id = "seurat-dataset",
  dataset_name = "My Seurat dataset",
  obs = seu[[]],
  gene_expression = expr,
  connectivities = graph,
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca",
  obs_keys = c("seurat_clusters", "orig.ident", "nCount_RNA")
)
```

Open `seurat_dataset.zarr` on a web host, or write `seurat_dataset.zarr.zip` with
`format = "zarr.zip"` and open it with **Open… › File…** ({doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`).

## What each input becomes

| Seurat | Argument | In the store |
| --- | --- | --- |
| `Embeddings(seu, "umap")` | `embeddings$X_umap` | `obsm['X_umap']`, every column, float32 |
| `Embeddings(seu, "pca")` | `embeddings$X_pca` | `obsm['X_pca']`; `latent_key = "X_pca"` makes it the space outlier quantiles are computed in |
| `seu[[]]` (the metadata) | `obs` | `obs`: factors stay categoricals with their level order, numbers stay numbers |
| `GetAssayData(seu, layer = "data")` | `gene_expression` (transposed) | `X`, gene-major float32 CSC |
| `rownames(seu)` | (the column names of `expr`) | the `var` index (gene names) |
| `seu@graphs[["RNA_snn"]]` | `connectivities` | `obsp['connectivities']`, float32 CSR |

Every input carries the cell names (`colnames(seu)`), and `cellucid_prepare()`
checks that they agree in the same order; they become the `obs` index.

## Choices and variants

- **A 3D UMAP.** Compute it under its own name, for example
  `RunUMAP(seu, dims = 1:30, n.components = 3, reduction.name = "umap3d")`, and
  add `X_umap_3d = Embeddings(seu, "umap3d")` to `embeddings`.
- **Only some genes.** Pass `gene_identifiers = c("CD3E", "MS4A1", …)`; the
  store then holds only those columns of `X`.
- **Gene symbols beside other ids.** Pass a `var` data frame (one row per gene,
  row names `rownames(seu)`) and name the column the viewer shows with
  `var_gene_id_column`.
- **SCTransform data.** Use `GetAssayData(seu, assay = "SCT", layer = "data")`.
- **The kNN graph instead of SNN.** `RNA_nn` is not symmetric, and the stored
  graph must be: symmetrize it first, for example
  `graph <- seu@graphs[["RNA_nn"]]; graph <- (graph + Matrix::t(graph)) / 2;
  Matrix::diag(graph) <- 0`.
- **Dataset description and citation.** Add `dataset_description`,
  `citation_text`, `citation_url`, `source_name`, `source_url`.

## Troubleshooting

- *"gene_expression has … rows, but the embeddings have … cells"* — you passed
  genes x cells; transpose with `Matrix::t()`.
- *"Connectivity diagonal values must all be exactly zero"* — set
  `Matrix::diag(graph) <- 0`.
- *"Connectivity matrix must be exactly symmetric"* — use `RNA_snn`, or
  symmetrize as above.
- *"obs column '…' is a … column, which has no AnnData encoding"* — convert that
  metadata column with `as.character()` or `as.numeric()`, or drop it.
- More: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
