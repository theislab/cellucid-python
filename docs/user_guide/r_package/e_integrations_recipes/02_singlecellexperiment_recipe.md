# SingleCellExperiment Recipe

**Audience:** Bioconductor users
**Time:** 10–20 minutes
**What you'll do:** extract the inputs from a `SingleCellExperiment` → write a
Cellucid-ready AnnData store → open it in the Cellucid web app.

```{note}
`cellucid-r` does not depend on Bioconductor. The code below needs
`SingleCellExperiment`, `SummarizedExperiment`, and `Matrix` installed.
```

## The whole recipe

```r
library(cellucid)
library(SingleCellExperiment)

# sce <- readRDS("my_sce.rds")

# Assays are genes x cells; cellucid_prepare() takes cells x genes.
expr <- Matrix::t(as(logcounts(sce), "CsparseMatrix"))

store <- cellucid_prepare(
  out_path = "sce_dataset.zarr",
  format = "zarr",
  embeddings = list(
    X_umap = reducedDim(sce, "UMAP"),
    X_pca = reducedDim(sce, "PCA")
  ),
  dataset_id = "sce-dataset",
  dataset_name = "My SingleCellExperiment dataset",
  obs = as.data.frame(colData(sce)),
  var = as.data.frame(rowData(sce)),
  gene_expression = expr,
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca"
)
```

`as(logcounts(sce), "CsparseMatrix")` turns a sparse or dense in-memory assay
into a `Matrix` sparse matrix. An on-disk assay (a `DelayedMatrix`, for example
from `HDF5Array`) has to be read into memory this way first: `cellucid_prepare()`
accepts base matrices and `Matrix` matrices.

## What each input becomes

| SingleCellExperiment | Argument | In the store |
| --- | --- | --- |
| `reducedDim(sce, "UMAP")` | `embeddings$X_umap` | `obsm['X_umap']`, every column, float32 |
| `reducedDim(sce, "PCA")` | `embeddings$X_pca` | `obsm['X_pca']`, the latent space with `latent_key = "X_pca"` |
| `colData(sce)` | `obs` | `obs` |
| `rowData(sce)` | `var` | `var`; its row names (`rownames(sce)`) are the gene names |
| `logcounts(sce)` | `gene_expression` (transposed) | `X`, gene-major float32 CSC |

The cell names (`colnames(sce)`) are carried by the reduced dimensions, the
`colData`, and the transposed assay; `cellucid_prepare()` checks that they agree
and stores them as the `obs` index.

## Choices and variants

- **Gene symbols.** When `rownames(sce)` are Ensembl ids and `rowData(sce)` has a
  symbol column, show the symbols with `var_gene_id_column = "Symbol"` (or
  whatever the column is called). A missing symbol takes the row's Ensembl id;
  repeated symbols are made unique (`MATR3`, `MATR3-1`).
- **A neighbour graph.** `bluster`/`scran` graphs are igraph objects; convert
  one to a symmetric sparse matrix with a zero diagonal, for example
  `graph <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE)`,
  then pass `connectivities = graph`.
- **Only some genes.** `gene_identifiers = c("CD3E", "MS4A1", …)`.

## Troubleshooting

- *"obs column '…' is a … column, which has no AnnData encoding"* — `colData` can
  hold list or nested `DataFrame` columns; convert or drop them before passing
  `obs`.
- *"The row names of … and … disagree"* — one input is in a different cell
  order; subset every input with the same `colnames(sce)`.
- *"gene_expression must be a real numeric matrix or a Matrix::Matrix"* — the
  assay is on disk; read it with `as(…, "CsparseMatrix")`.
- More: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
