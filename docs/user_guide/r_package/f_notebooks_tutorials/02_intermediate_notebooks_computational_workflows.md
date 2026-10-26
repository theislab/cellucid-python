# Intermediate Tutorial (Computational Workflows)

**Audience:** computational users
**Time:** 20 minutes
**Goal:** write expression, gene symbols, several embeddings, launch defaults and
citation metadata, and check the result.

## Step 1 — Simulate a dataset

```r
library(cellucid)
set.seed(11)

n_cells <- 2000
n_genes <- 500
cells <- sprintf("cell%04d", seq_len(n_cells))
ensembl <- sprintf("ENSG%011d", seq_len(n_genes))
symbols <- sprintf("GENE%d", seq_len(n_genes))
symbols[c(10, 20)] <- "MATR3"       # a repeated symbol
symbols[30] <- NA                   # a gene without a symbol

expr <- Matrix::rsparsematrix(n_cells, n_genes, density = 0.08,
                              rand.x = function(k) round(log1p(rpois(k, 3) + 1), 3))
dimnames(expr) <- list(cells, ensembl)

var <- data.frame(symbol = symbols, highly_variable = seq_len(n_genes) <= 100,
                  row.names = ensembl)
obs <- data.frame(
  leiden = factor(sample(0:7, n_cells, TRUE)),
  donor = sample(sprintf("D%d", 1:4), n_cells, TRUE),
  pct_mito = runif(n_cells, 0, 0.1),
  row.names = cells
)
pca <- matrix(rnorm(n_cells * 30), ncol = 30, dimnames = list(cells, NULL))
umap <- pca[, 1:2] + matrix(rnorm(n_cells * 2, sd = 0.1), ncol = 2)
umap3 <- pca[, 1:3]
```

## Step 2 — Write the store

```r
store <- cellucid_prepare(
  out_path = file.path(tempdir(), "intermediate.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap, X_umap_3d = umap3, X_pca = pca),
  dataset_id = "intermediate-tutorial",
  dataset_name = "Intermediate tutorial",
  dataset_description = "2,000 simulated cells from four donors.",
  citation_text = "Doe J et al. A simulated atlas. 2026.",
  citation_url = "https://example.org/simulated-atlas",
  source_name = "Simulation",
  obs = obs,
  var = var,
  gene_expression = expr,
  embedding_defaults = list("2d" = "X_umap", "3d" = "X_umap_3d"),
  latent_key = "X_pca",
  obs_keys = c("leiden", "donor", "pct_mito"),
  var_gene_id_column = "symbol",
  centroid_min_points = 20,
  force = TRUE
)
```

What each argument does:

- `embeddings` — all three arrays are stored whole in `obsm`; the viewer plots
  `X_umap` in 2D and `X_umap_3d` in 3D first (`embedding_defaults`) and offers
  `X_pca` too, by any of its columns.
- `latent_key` — the outlier quantiles of `leiden`, `donor` and every other
  categorical `obs` column are computed in `X_pca` and stored in
  `uns['cellucid']['outlierQuantiles']`, so the viewer's outlier filter works
  as soon as the file opens; `centroid_min_points` sets how many cells a
  category needs for them (10 when not set).
- `obs_keys` — the fields the viewer lists, in this order. Every `obs` column is
  stored either way.
- `var_gene_id_column` — the viewer names genes by `symbol`. The gene without a
  symbol is shown under its Ensembl id, and the second `MATR3` as `MATR3-1`.
- The metadata arguments end up in `uns['cellucid']`; only the defaults you set
  are stored.

## Step 3 — Only the genes you need

A gene panel keeps the store small: `gene_identifiers` names the genes by the
names the viewer shows (here, the symbols).

```r
panel <- cellucid_prepare(
  out_path = file.path(tempdir(), "panel.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap),
  dataset_id = "intermediate-panel",
  dataset_name = "Intermediate tutorial, gene panel",
  obs = obs,
  var = var,
  gene_expression = expr,
  var_gene_id_column = "symbol",
  gene_identifiers = c("GENE1", "GENE2", "MATR3", "MATR3-1"),
  force = TRUE
)
dir_size <- function(path) sum(file.size(list.files(path, recursive = TRUE, all.files = TRUE, full.names = TRUE)))
c(all_genes = dir_size(store), panel = dir_size(panel))
```

## Step 4 — Check what was written

```r
cellucid_dir <- file.path(store, "uns", "cellucid")
list.dirs(cellucid_dir, full.names = FALSE)
jsonlite::read_json(file.path(store, "X", ".zattrs"))
```

`X` is a `csc_matrix`: gene-major, so the viewer reads one gene as one
contiguous range. In Python:

```python
import anndata
adata = anndata.read_zarr("intermediate.zarr")
adata.uns["cellucid"]["defaults"]
```

Next: {doc}`03_advanced_notebooks_expert_developer` adds vector fields and a
neighbour graph.
