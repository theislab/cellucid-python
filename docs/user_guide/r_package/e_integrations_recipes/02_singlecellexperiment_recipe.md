# SingleCellExperiment Recipe

**Audience:** SingleCellExperiment users (wet lab + computational)  
**Time:** 20–50 minutes  
**What you’ll do:** extract reducedDims/colData/assays → export → open in Cellucid

This recipe exports a `SingleCellExperiment` (SCE) object by extracting the arrays that `cellucid_prepare()` needs.

```{note}
`cellucid-r` does not depend on SingleCellExperiment. This page assumes you already have an SCE object.
```

## Prerequisites

- An SCE object `sce` with:
  - a reduced dimension representation you want to visualize (often `"UMAP"`)
  - a latent representation (often `"PCA"`)
  - optional expression assay (e.g., `"logcounts"` or `"counts"`)
- R packages:
  - `cellucid`
  - `SingleCellExperiment`
  - `SummarizedExperiment`
  - `Matrix` (recommended)

## Step 1 — Load packages and the SCE

```r
library(cellucid)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

# Example:
# sce <- readRDS("my_sce.rds")
```

## Step 2 — Choose a stable cell order

```r
cells <- colnames(sce)
stopifnot(length(cells) > 0)
```

## Step 3 — Extract embeddings (reducedDims)

List available reduced dimensions:

```r
reducedDimNames(sce)
```

Extract UMAP (2D example):

```r
umap <- reducedDim(sce, "UMAP")          # rows = cells
umap2 <- umap[, 1:2, drop = FALSE]
umap2 <- umap2[cells, , drop = FALSE]   # align order
```

If you have a 3D UMAP:

```r
# umap3 <- umap[, 1:3, drop = FALSE]
# umap3 <- umap3[cells, , drop = FALSE]
```

## Step 4 — Extract a latent space (PCA is typical)

```r
pca <- reducedDim(sce, "PCA")
pca <- pca[cells, , drop = FALSE]
latent_space <- pca
```

If you don’t have PCA, you can use the embedding as latent space:

```r
# latent_space <- umap2
```

## Step 5 — Build `obs` (cell metadata)

Cell metadata is in `colData(sce)`:

```r
obs_all <- as.data.frame(colData(sce))
obs_all <- obs_all[cells, , drop = FALSE]
```

Optionally select a subset of columns:

```r
obs_keys <- c("cluster", "sample", "nCount", "nFeature")
obs_keys <- intersect(obs_keys, colnames(obs_all))
obs <- obs_all[, obs_keys, drop = FALSE]
```

Make sure important columns are the intended types:

```r
if ("cluster" %in% colnames(obs)) {
  obs$cluster <- factor(obs$cluster)
}
```

## Step 6 — (Optional) Export gene expression + `var`

Declare the optional inputs up front so Step 8 can pass them unconditionally. If
you skip this whole step, they stay `NULL` and no gene expression is exported:

```r
gene_var <- NULL
expr_cxg <- NULL
gene_panel <- NULL
```

### 6.1 Choose which assay to export

List available assays:

```r
assayNames(sce)
```

Pick one (common choices: `"logcounts"`, `"counts"`):

```r
expr_gxc <- assay(sce, "logcounts")   # genes x cells
expr_gxc <- expr_gxc[, cells, drop = FALSE]
```

### 6.2 Convert to a supported matrix type (important)

SCE assays can be dense, sparse, or delayed.

`cellucid-r` supports:
- base `matrix`
- `Matrix` sparse matrices

If `expr_gxc` is already a `dgCMatrix`, you’re good.

If it’s a delayed representation, you may need to realize it:

```r
# For small/medium data only (can be huge memory):
# expr_gxc <- as.matrix(expr_gxc)
```

If you can convert to sparse:

```r
# expr_gxc <- as(expr_gxc, "dgCMatrix")
```

### 6.3 Transpose to cells × genes

```r
expr_cxg <- Matrix::t(expr_gxc)
```

### 6.4 Build `var`

Gene metadata is in `rowData(sce)`; gene IDs are usually `rownames(sce)`.

```r
gene_ids <- rownames(sce)
gene_var <- as.data.frame(rowData(sce))
gene_var$symbol <- gene_ids
rownames(gene_var) <- gene_var$symbol
```

:::{note}
The variable is called `gene_var`, not `var`, because `var()` is a base R
function: a script that only *sometimes* defines `var` will silently pass the
function along when it does not. Name the argument explicitly at the call site
(`var = gene_var`).
:::

### 6.5 Optional: export a gene panel

```r
gene_panel <- c("MS4A1", "CD3D", "NKG7")
gene_panel <- intersect(gene_panel, rownames(gene_var))
```

## Step 7 — (Optional) Connectivities

SCE does not have one canonical place for neighbor graphs; workflows vary.

If you have an adjacency matrix `conn` with shape `(n_cells, n_cells)` (dense or sparse), you can export it:

```r
conn <- NULL
if ("connectivities" %in% names(metadata(sce))) {
  conn <- metadata(sce)$connectivities
  if (!identical(dim(conn), c(length(cells), length(cells)))) {
    stop("metadata(sce)$connectivities must be n_cells x n_cells")
  }
  if (!is.null(rownames(conn)) && !is.null(colnames(conn))) {
    conn <- conn[cells, cells, drop = FALSE]
  }
}
```

If your workflow stores its graph somewhere else, assign that exact matrix to
`conn` and apply the same dimension and cell-order checks before export. Leaving
`conn <- NULL` intentionally exports no connectivity graph.

## Step 8 — Run `cellucid_prepare()`

```r
out_dir <- file.path(getwd(), "exports", "sce_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  X_umap_2d = umap2,
  obs_categorical_dtype = "uint16",
  dataset_name = "SCE Export",
  dataset_id = "sce_export",
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  obs_continuous_quantization = 8,
  var = gene_var,
  gene_expression = expr_cxg,
  gene_identifiers = gene_panel,
  var_quantization = if (is.null(expr_cxg)) NULL else 8,
  connectivities = conn
)
```

## Step 9 — Open in the web app

Follow: {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`

---

## Common gotchas (SCE-specific)

### Assay orientation

Most assays are genes × cells, so you must transpose to cells × genes.

### DelayedArray / HDF5-backed assays

If the assay is delayed, converting to dense can be impossible for large data.

Mitigations:
- export fewer genes (gene panel)
- export fewer cells (subset)
- export only embeddings + obs (no expression) to validate the workflow first

### Missing reducedDims

If `reducedDim(sce, "UMAP")` fails, check `reducedDimNames(sce)` and use the correct name.

---

## Troubleshooting pointers

- Export errors about mismatched shapes → {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
- Browser loading issues → {doc}`../i_troubleshooting_index/04_web_app_loading_issues`
