# Quick Start (3 Levels)

**Audience:** everyone (choose your depth)  
**Time:** 5–30 minutes  
**What you’ll get:** a Cellucid-ready AnnData store that opens in the Cellucid web app

Stop as soon as you have what you need.

- **Level 1**: tiny toy dataset (proves your installation works)
- **Level 2**: real data from Seurat / SingleCellExperiment
- **Level 3**: everything at once — expression, graph, vector field, metadata, launch defaults, one `.zarr.zip` file

```{note}
For a notebook-style, verbose walkthrough, go to {doc}`../f_notebooks_tutorials/index`.
```

## Level 1 — Minimal toy store (5 minutes)

```r
library(cellucid)

obs <- data.frame(
  cluster = factor(c("A", "A", "B")),
  score = c(0.1, 0.2, 0.3)
)
umap <- matrix(c(0, 0,
                 1, 0,
                 0, 1),
               ncol = 2, byrow = TRUE)

store <- cellucid_prepare(
  out_path = file.path(tempdir(), "quickstart_toy.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap),
  dataset_id = "quickstart-toy",
  dataset_name = "Quick Start Toy",
  obs = obs,
  force = TRUE
)

list.files(store, all.files = TRUE, no.. = TRUE)
#>  [1] ".zattrs"    ".zgroup"    ".zmetadata" "layers"     "obs"
#>  [6] "obsm"       "obsp"       "uns"        "var"        "varm"
#> [11] "varp"
```

`cellucid_prepare()` returns the store's absolute path (invisibly).

**Open it:** {doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`.

**What success looks like**
- Cellucid shows 3 points.
- You can colour by `cluster` and `score`.

## Level 2 — From Seurat or SingleCellExperiment (10–20 minutes)

Extract the pieces and pass them; transpose expression to cells × genes.

```r
# Seurat
cellucid_prepare(
  out_path = "pbmc.zarr",
  format = "zarr",
  embeddings = list(
    X_umap = Seurat::Embeddings(seu, "umap"),
    X_pca = Seurat::Embeddings(seu, "pca")
  ),
  dataset_id = "pbmc",
  dataset_name = "PBMC",
  obs = seu[[]],
  gene_expression = Matrix::t(Seurat::GetAssayData(seu, layer = "data")),
  embedding_defaults = c("2d" = "X_umap"),
  latent_key = "X_pca"
)
```

Full recipes, including the neighbour graph and gene symbols:
- {doc}`../e_integrations_recipes/01_seurat_recipe`
- {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`

## Level 3 — Everything at once (15–30 minutes)

```r
library(cellucid)
library(Matrix)

set.seed(1)
n_cells <- 500
n_genes <- 200

expr <- rsparsematrix(n_cells, n_genes, density = 0.05)   # cells x genes
expr@x <- abs(expr@x)
dimnames(expr) <- list(paste0("cell", seq_len(n_cells)), paste0("ENSG", seq_len(n_genes)))
var <- data.frame(symbol = paste0("GENE", seq_len(n_genes)), row.names = colnames(expr))
obs <- data.frame(
  cluster = factor(sample(c("T", "B", "NK"), n_cells, TRUE)),
  row.names = rownames(expr)
)
pca <- matrix(rnorm(n_cells * 20), n_cells)
umap <- matrix(rnorm(n_cells * 2), n_cells)
velocity <- matrix(rnorm(n_cells * 2, sd = 0.1), n_cells)

# A symmetric neighbour graph with a zero diagonal
knn <- rsparsematrix(n_cells, n_cells, density = 0.01)
knn@x <- abs(knn@x)
knn <- knn + t(knn)
diag(knn) <- 0
knn <- drop0(knn)

store <- cellucid_prepare(
  out_path = file.path(tempdir(), "quickstart_advanced.zarr.zip"),
  format = "zarr.zip",
  embeddings = list(X_umap = umap, X_pca = pca),
  dataset_id = "quickstart-advanced",
  dataset_name = "Quick Start Advanced",
  dataset_description = "Simulated cells for the quick start.",
  citation_text = "Doe et al. (2026)",
  citation_url = "https://example.org/doe2026",
  source_name = "Simulation",
  obs = obs,
  var = var,
  gene_expression = expr,
  connectivities = knn,
  vector_fields = list(velocity = list(X_umap = velocity)),
  embedding_defaults = list("2d" = "X_umap", "3d" = list(key = "X_pca", columns = 1:3)),
  latent_key = "X_pca",
  var_gene_id_column = "symbol",
  vector_field_default = "velocity",
  force = TRUE
)
```

What each option does:

- `format = "zarr.zip"` writes one `.zarr.zip` file instead of a `.zarr` folder — for everything but GitHub: single-file hosts, local files, `cellucid serve` and Jupyter.
- `vector_fields` stores `obsm['velocity_umap']` (scVelo's name for velocity on `X_umap`), drawn whenever a view shows `X_umap`
  ({doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`).
- `embedding_defaults`, `latent_key`, `var_gene_id_column`, and `vector_field_default`
  are **launch defaults** stored in `uns['cellucid']['defaults']`; the viewer
  chooses for itself wherever you set none. `latent_key` also stores each
  categorical `obs` column's outlier quantiles in `X_pca`
  (`uns['cellucid']['outlierQuantiles']`), which the outlier filter reads.
- `gene_identifiers = c(...)` would store only those genes.

## Next steps

- What exactly is in the store? {doc}`../b_concepts_mental_models/01_what_is_a_cellucid_ready_anndata`
- Every argument in detail: {doc}`../c_data_preparation_api/index`
- Share it with others: {doc}`../d_viewing_loading/02_host_stores_for_sharing`
