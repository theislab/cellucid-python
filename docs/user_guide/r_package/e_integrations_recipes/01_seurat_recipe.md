# Seurat Recipe

**Audience:** Seurat users (wet lab + computational)  
**Time:** 20–45 minutes (depends on dataset size)  
**What you’ll do:** extract data from a Seurat object → export → open in the Cellucid web app

This recipe exports a Seurat object using `cellucid::cellucid_prepare()` (alias: `prepare()`).

```{note}
`cellucid-r` does not depend on Seurat. This page shows how to *extract* the needed inputs from Seurat.
```

## Prerequisites

- A Seurat object `seu`
  - must have a UMAP (or other) reduction you want to visualize
  - should have a latent space (PCA is typical)
- R packages:
  - `cellucid`
  - `Seurat`
  - `Matrix` (recommended)

## Step 0 — Decide what you want to export

Before writing code, decide:

1) Which embedding?
   - 2D UMAP (most common)
   - 3D UMAP (if you computed 3 components)
2) Which latent space?
   - PCA is typical
3) Which metadata columns (`obs`) matter?
   - clusters, sample, batch, QC metrics, scores…
4) Do you need gene expression?
   - exporting all genes can be huge; consider a gene panel
5) Do you need a connectivity graph?
   - Seurat stores graphs; exporting them is optional

If you don’t know yet:
- start with embeddings + a handful of obs columns + a small gene panel

## Step 1 — Load packages and the Seurat object

```r
library(cellucid)
library(Seurat)
library(Matrix)

# Example:
# seu <- readRDS("my_seurat_object.rds")
```

## Step 2 — Choose a stable cell order

Cellucid uses row order as cell identity, so we pick a canonical order and reuse it everywhere.

```r
cells <- colnames(seu)
stopifnot(length(cells) > 0)
```

## Step 3 — Extract embeddings (UMAP 2D / 3D)

### 2D UMAP

```r
umap <- Embeddings(seu, reduction = "umap")
umap2 <- umap[, 1:2, drop = FALSE]
umap2 <- umap2[cells, , drop = FALSE]
```

### Optional: 3D UMAP

Only do this if you actually computed a 3D UMAP (e.g., `RunUMAP(..., n.components = 3)`).

```r
# umap3 <- umap[, 1:3, drop = FALSE]
# umap3 <- umap3[cells, , drop = FALSE]
```

## Step 4 — Extract a latent space (PCA is typical)

```r
pca <- Embeddings(seu, reduction = "pca")
pca <- pca[cells, , drop = FALSE]
latent_space <- pca
```

```{note}
If you don’t have PCA, you can use the embedding as latent space (less ideal, but valid):
`latent_space <- umap2`
```

## Step 5 — Build `obs` (cell metadata)

Seurat metadata lives in `seu@meta.data` (rows = cells).

```r
obs_all <- seu@meta.data
obs_all <- obs_all[cells, , drop = FALSE]
```

### Optional: select only a subset of obs columns

Exporting dozens/hundreds of columns is allowed, but can be noisy in the UI.

```r
obs_keys <- c("seurat_clusters", "orig.ident", "nCount_RNA", "nFeature_RNA")
obs_keys <- intersect(obs_keys, colnames(obs_all))
obs <- obs_all[, obs_keys, drop = FALSE]
```

Make sure key columns are the types you want:

```r
obs$seurat_clusters <- factor(obs$seurat_clusters)
```

## Step 6 — (Optional) Export gene expression + `var`

### 6.1 Choose which assay/slot to export

Common Seurat choices:
- `slot="data"`: log-normalized (often best for visualization)
- `slot="counts"`: raw counts (can have huge dynamic range)

```r
assay_name <- DefaultAssay(seu)  # often "RNA"
expr_gxc <- GetAssayData(seu, assay = assay_name, slot = "data")  # genes x cells
```

### 6.2 Align cells and transpose to cells × genes

`GetAssayData` returns **genes × cells**. Cellucid expects **cells × genes**.

```r
expr_gxc <- expr_gxc[, cells, drop = FALSE]
expr_cxg <- Matrix::t(expr_gxc)
```

### 6.3 Build `var` (gene metadata)

The simplest `var` is a data.frame with rownames set to the gene IDs you want users to search.

```r
gene_ids <- rownames(expr_gxc)
var <- data.frame(symbol = gene_ids, stringsAsFactors = FALSE)
rownames(var) <- var$symbol
```

### 6.4 Optional: export a gene panel instead of all genes

Exporting every gene can be huge. A panel often gives a better experience.

```r
gene_panel <- c("MS4A1", "CD3D", "NKG7")
gene_panel <- intersect(gene_panel, rownames(var))
```

## Step 7 — (Optional) Export connectivities (Seurat graphs)

Seurat graphs often live in `seu@graphs`, for example:
- `RNA_snn` (shared nearest neighbor graph)

```r
conn <- NULL
if ("RNA_snn" %in% names(seu@graphs)) {
  conn <- seu@graphs$RNA_snn
  conn <- conn[cells, cells]
}
```

```{note}
Connectivity export symmetrizes and binarizes edges. Weights are not preserved.
```

## Step 8 — Run `cellucid_prepare()`

Pick an output directory (best: a fresh folder per export iteration):

```r
out_dir <- file.path(getwd(), "exports", "seurat_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
```

### Recommended “starter” settings

- `force=TRUE` while iterating
- `obs_continuous_quantization=8` for smaller metadata
- if exporting expression: `var_quantization=8` and/or `gene_identifiers=gene_panel`

```r
cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  X_umap_2d = umap2,
  # X_umap_3d = umap3,  # optional
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  obs_continuous_quantization = 8,
  var = if (exists("var")) var else NULL,
  gene_expression = if (exists("expr_cxg")) expr_cxg else NULL,
  gene_identifiers = if (exists("gene_panel")) gene_panel else NULL,
  var_quantization = if (exists("expr_cxg")) 8 else NULL,
  connectivities = conn
)
```

## Step 9 — Quick validation

```r
list.files(out_dir, recursive = TRUE)
ident <- jsonlite::read_json(file.path(out_dir, "dataset_identity.json"), simplifyVector = TRUE)
ident$stats
```

Confirm:
- `n_cells` matches your Seurat object
- `n_genes` matches the exported gene panel (if you used one)

## Step 10 — Open in the Cellucid web app

Follow: {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`

**What success looks like**
- The web app loads points quickly.
- You can color by your `obs` fields.
- If you exported gene expression, you can search genes and color by them.

---

## Edge cases and common gotchas

### “My metadata columns turned into categories”

If a column is not numeric, it becomes categorical. Explicitly coerce numeric columns.
See: {doc}`../c_data_preparation_api/04_obs_cell_metadata`

### “I exported and got a gigantic folder”

You probably exported too many genes for your cell count.
Use:
- `gene_identifiers` (gene panels)
- `var_quantization=8`
- `compression=6`

### “My UMAP is missing / Embeddings(...) errors”

Check which reductions exist:

```r
Reductions(seu)
```

### “Cell order mismatch”

Always reorder:
- embeddings
- metadata
- expression columns/rows
- connectivity matrix

using the same `cells` vector.

---

## Troubleshooting (symptom → fix)

### Symptom: export fails with shape mismatches

Most common causes:
- expression orientation is wrong (genes × cells vs cells × genes)
- you subset a matrix but not the others

Confirm:
- `dim(expr_cxg)` is `(n_cells, n_genes)`
- `nrow(obs) == nrow(umap2) == nrow(latent_space)`

Fix:
- reorder/transpose explicitly (Steps 2–6)

### Symptom: web app loads but fields look wrong

Likely cause:
- row order mismatch (the export “worked” but you exported misaligned inputs)

Fix:
- rebuild using a single `cells <- colnames(seu)` ordering and reorder everything.

More:
- export-time troubleshooting: {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
- loading-time troubleshooting: {doc}`../i_troubleshooting_index/04_web_app_loading_issues`
