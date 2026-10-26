# Raw Matrices + Data Frames

**Audience:** anyone whose data is in plain R objects
**Time:** 5 minutes

`cellucid_prepare()` takes base R matrices, numeric data frames, and `Matrix`
matrices, so no single-cell container is needed. This page runs as it is, on
simulated data.

## A complete example

```r
library(cellucid)
set.seed(42)

n_cells <- 500
n_genes <- 200
cells <- sprintf("cell_%03d", seq_len(n_cells))
genes <- sprintf("ENSG%05d", seq_len(n_genes))

# Embeddings: one row per cell. Row names are optional; when present they
# must be the same cells in the same order as every other input.
umap <- matrix(rnorm(n_cells * 2), ncol = 2, dimnames = list(cells, NULL))
pca <- matrix(rnorm(n_cells * 20), ncol = 20, dimnames = list(cells, NULL))

# Cell metadata: one row per cell.
obs <- data.frame(
  cell_type = factor(sample(c("T cell", "B cell", "NK cell"), n_cells, TRUE)),
  sample = sample(c("donor1", "donor2"), n_cells, TRUE),
  n_counts = rpois(n_cells, 2000),
  row.names = cells
)

# Gene expression: cells x genes. A dense matrix works; a sparse one is smaller.
expr <- Matrix::rsparsematrix(n_cells, n_genes, density = 0.05,
                              rand.x = function(k) round(abs(rnorm(k)), 2))
dimnames(expr) <- list(cells, genes)

# Gene metadata: one row per gene, row names = gene ids.
var <- data.frame(
  symbol = sprintf("GENE%d", seq_len(n_genes)),
  row.names = genes
)

store <- cellucid_prepare(
  out_path = file.path(tempdir(), "raw.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap, X_pca = pca),
  dataset_id = "raw-example",
  dataset_name = "Raw matrices example",
  obs = obs,
  var = var,
  gene_expression = expr,
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca",
  var_gene_id_column = "symbol",
  force = TRUE
)
store
```

## Rules worth knowing

- **Orientation.** `gene_expression` is cells x genes. It is never transposed for
  you; a genes x cells matrix is refused because its row count does not match.
- **Names.** Cell names come from the row names the inputs agree on (else `"0"`
  to `"n - 1"`); gene names from `rownames(var)` (else the column names of
  `gene_expression`).
- **Columns.** A factor keeps its level order; a character column with repeats is
  stored as a categorical; numbers stay numbers; `NA` is stored as missing.
- **Values.** Embeddings, expression and vector fields are stored as float32 and
  must be finite.

Each rule in full: {doc}`../c_data_preparation_api/02_input_requirements_global`.
