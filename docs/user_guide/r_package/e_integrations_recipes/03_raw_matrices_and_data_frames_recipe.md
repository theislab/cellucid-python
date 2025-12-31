# Raw Matrices + Data Frames (Minimal Dependencies)

**Audience:** everyone  
**Time:** 10–25 minutes  
**Goal:** export data without Seurat/SCE (just matrices + data.frames).

This is the minimal-dependency workflow for `cellucid-r`.

You provide:
- embeddings (UMAP or anything embedding-like)
- `obs` metadata
- a latent space
- optional expression + var
- optional connectivity + vector fields

## Step-by-step: minimal export (no expression)

```r
library(cellucid)

# 1) Embedding (cells x 2)
X_umap_2d <- matrix(rnorm(2000), ncol = 2)

# 2) Latent space (cells x dims)
latent_space <- matrix(rnorm(2000), ncol = 10)

# 3) obs metadata (data.frame with n_cells rows)
obs <- data.frame(
  cluster = factor(sample(letters[1:5], size = nrow(X_umap_2d), replace = TRUE)),
  score = rnorm(nrow(X_umap_2d))
)

out_dir <- file.path(getwd(), "exports", "raw_minimal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  X_umap_2d = X_umap_2d,
  out_dir = out_dir,
  force = TRUE
)
```

## Add gene expression (optional)

Remember: Cellucid expects **cells × genes**.

```r
n_cells <- nrow(X_umap_2d)
n_genes <- 100

gene_expression <- matrix(rexp(n_cells * n_genes, rate = 2), nrow = n_cells, ncol = n_genes)
gene_ids <- paste0("G", seq_len(n_genes))

var <- data.frame(symbol = gene_ids, stringsAsFactors = FALSE)
rownames(var) <- var$symbol

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  var = var,
  gene_expression = gene_expression,
  X_umap_2d = X_umap_2d,
  out_dir = out_dir,
  force = TRUE,
  var_quantization = 8,
  compression = 6
)
```

## Common “raw workflow” pitfalls

- Mismatched row order (if you build pieces from different sources)
- Expression orientation (genes × cells vs cells × genes)
- Metadata columns with unexpected types (character/date becoming categorical)

See:
- global requirements: {doc}`../c_data_preparation_api/02_input_requirements_global`
- obs rules: {doc}`../c_data_preparation_api/04_obs_cell_metadata`
- expression rules: {doc}`../c_data_preparation_api/06_gene_expression_matrix`

## Next steps

- Open in the web app: {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`
- If something fails: {doc}`../i_troubleshooting_index/index`
