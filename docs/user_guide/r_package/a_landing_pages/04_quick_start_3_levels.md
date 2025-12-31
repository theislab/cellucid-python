# Quick Start (3 Levels)

**Audience:** everyone (choose your depth)  
**Time:** 5–30 minutes  
**What you’ll get:** an export folder that loads in the Cellucid web app

This page gives you three “depth levels”. Stop as soon as you have what you need.

- **Level 1**: tiny toy dataset (proves your installation works)
- **Level 2**: “real data” from Seurat / SingleCellExperiment
- **Level 3**: advanced export (compression, quantization, subsets, connectivity, vector fields)

```{note}
If you prefer a notebook-style, highly verbose walkthrough, go to {doc}`../f_notebooks_tutorials/index`.
```

## Level 1 — Minimal toy export (5 minutes)

### Prerequisites

```r
library(cellucid)
```

### Step 1: Create small inputs

```r
latent <- matrix(c(0, 0,
                   1, 1,
                   2, 2),
                 ncol = 2, byrow = TRUE)

obs <- data.frame(
  cluster = factor(c("A", "A", "B")),
  score = c(0.1, 0.2, 0.3)
)

umap2 <- matrix(c(0, 0,
                  1, 0,
                  0, 1),
                ncol = 2, byrow = TRUE)
```

### Step 2: Export

```r
out_dir <- file.path(tempdir(), "cellucid_quickstart_toy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  X_umap_2d = umap2,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)

list.files(out_dir, recursive = TRUE)
```

### Step 3: Open in the web app

Go to {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app` for a click-by-click guide.

**What success looks like**
- Cellucid loads quickly and shows 3 points.
- You can color by `cluster` and `score`.

## Level 2 — Export from Seurat or SingleCellExperiment (10–20 minutes)

You can export from common R single-cell containers by extracting:
- embeddings (UMAP)
- a latent space (often PCA)
- `obs` metadata
- (optional) gene expression and connectivity

Pick one:

- Seurat: {doc}`../e_integrations_recipes/01_seurat_recipe`
- SingleCellExperiment: {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`

## Level 3 — Advanced export (quality + performance + extra features)

This is for large datasets, sharing, or “publication-grade” exports where you want explicit control.

### A) Reduce disk size (recommended)

You have three main levers:

1) **Compression** (gzip)
2) **Quantization** (8-bit / 16-bit) for continuous fields and gene expression
3) **Export fewer genes / fewer metadata columns**

Example:

```r
cellucid_prepare(
  latent_space = latent,
  obs = obs,
  X_umap_2d = umap2,
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  obs_continuous_quantization = 8
)
```

### B) Add gene expression (optional but common)

If you include gene expression, you must also include `var` (gene metadata).

```r
expr <- matrix(c(0, 1,
                 2, 3,
                 4, 5),
               nrow = 3, ncol = 2, byrow = TRUE) # cells x genes

var <- data.frame(symbol = c("G1", "G2"))
rownames(var) <- var$symbol

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  var = var,
  gene_expression = expr,
  X_umap_2d = umap2,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE,
  compression = 6,
  var_quantization = 8
)
```

### C) Add a connectivity graph (optional)

If you have a KNN/SNN graph as an adjacency matrix (dense or sparse), pass it as `connectivities=`.

```r
connectivities <- matrix(0, nrow = 3, ncol = 3)
connectivities[1, 2] <- 1
connectivities[1, 3] <- 1

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  X_umap_2d = umap2,
  connectivities = connectivities,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)
```

### D) Add vector fields (velocity / displacement; optional)

Vector fields are per-cell vectors (1/2/3 components). They are exported as float32 and automatically scaled to match the embedding normalization.

```r
vector_fields <- list(
  velocity_umap_2d = matrix(c(0.2, 0,
                              0.2, 0,
                              0.2, 0),
                            ncol = 2, byrow = TRUE)
)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  X_umap_2d = umap2,
  vector_fields = vector_fields,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)
```

To understand naming rules and edge cases, see {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Quick troubleshooting

If anything fails:
- Export-time issues → {doc}`../i_troubleshooting_index/02_data_preparation_issues`
- Loading in the browser → {doc}`../i_troubleshooting_index/04_web_app_loading_issues`
