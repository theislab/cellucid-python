# Beginner Notebook (Wet Lab Friendly)

**Audience:** wet lab scientists, beginners, “I just want to see my data”  
**Time:** 15–30 minutes  
**What you’ll learn**
- What “exporting” means in Cellucid
- How to export from an R object (Seurat or SingleCellExperiment)
- How to load the exported folder in the Cellucid web app
- How to recognize success vs common failure modes

```{note}
This is written as a notebook-style walkthrough: copy/paste blocks in order.
If you prefer a shorter page, use the recipes:
- {doc}`../e_integrations_recipes/01_seurat_recipe`
- {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`
```

---

## Part 1 — What you are doing (in plain language)

Cellucid is a **web app** that shows your cells as points in 2D/3D.

To show your dataset, Cellucid needs files on disk (an “export folder”).

This tutorial:

1) takes your R object,
2) writes an export folder using `cellucid-r`,
3) then you open that folder in your browser.

---

## Part 2 — Install once

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

Then:

```r
library(cellucid)
```

If `library(cellucid)` fails, go to {doc}`../a_landing_pages/02_installation`.

---

## Part 3 — Choose your starting point

You probably have one of these:

1) A **Seurat object** (`.rds`)  
2) A **SingleCellExperiment** object (`.rds`)  
3) No object yet (you want to test the tool)

Pick the section that matches you.

---

## Option A — You have a Seurat object

### A1) Load it

```r
library(Seurat)

seu <- readRDS("path/to/my_seurat_object.rds")
```

### A2) Pick what you want to view

In Seurat, the 2D map is usually stored as a “reduction” called `"umap"`.

Check what you have:

```r
Reductions(seu)
```

If you see `"umap"`, great.

### A3) Export (copy/paste block)

This exports:
- UMAP 2D points
- basic metadata (clusters, counts/features if present)
- a small gene panel (optional; but recommended to keep exports smaller)

```r
library(Matrix)

cells <- colnames(seu)

umap <- Embeddings(seu, "umap")
umap2 <- umap[, 1:2, drop = FALSE]
umap2 <- umap2[cells, , drop = FALSE]

# Latent space (used internally for category summaries/outliers)
latent_space <- Embeddings(seu, "pca")[cells, , drop = FALSE]

# Metadata (cell table)
obs_all <- seu@meta.data[cells, , drop = FALSE]
obs_keys <- intersect(c("seurat_clusters", "orig.ident", "nCount_RNA", "nFeature_RNA"), colnames(obs_all))
obs <- obs_all[, obs_keys, drop = FALSE]
if ("seurat_clusters" %in% colnames(obs)) obs$seurat_clusters <- factor(obs$seurat_clusters)

# Gene expression (optional, but useful)
expr_gxc <- GetAssayData(seu, slot = "data")     # genes x cells
expr_gxc <- expr_gxc[, cells, drop = FALSE]
expr_cxg <- Matrix::t(expr_gxc)                  # cells x genes

gene_panel <- c("MS4A1", "CD3D", "NKG7")
gene_panel <- intersect(gene_panel, rownames(expr_gxc))

var <- data.frame(symbol = rownames(expr_gxc), stringsAsFactors = FALSE)
rownames(var) <- var$symbol

out_dir <- file.path(getwd(), "exports", "cellucid_from_seurat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  var = var,
  gene_expression = expr_cxg,
  gene_identifiers = gene_panel,
  X_umap_2d = umap2,
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  var_quantization = 8,
  obs_continuous_quantization = 8,
  centroid_min_points = 10
)
```

### A4) What success looks like

After running, you should see:

```r
list.files(out_dir, recursive = TRUE)
```

You should see files like:
- `dataset_identity.json`
- `points_2d.bin.gz` (or `.bin` if you disabled compression)
- `obs_manifest.json`
- `obs/`
- `var_manifest.json`
- `var/` with a few genes (your panel)

---

## Option B — You have a SingleCellExperiment object

### B1) Load it

```r
library(SingleCellExperiment)

sce <- readRDS("path/to/my_sce.rds")
```

### B2) Check which embeddings exist

```r
reducedDimNames(sce)
```

You want something like `"UMAP"` and ideally `"PCA"`.

### B3) Export (copy/paste block)

```r
library(SummarizedExperiment)
library(Matrix)

cells <- colnames(sce)

umap <- reducedDim(sce, "UMAP")
umap2 <- umap[, 1:2, drop = FALSE]
umap2 <- umap2[cells, , drop = FALSE]

latent_space <- reducedDim(sce, "PCA")[cells, , drop = FALSE]

obs_all <- as.data.frame(colData(sce))[cells, , drop = FALSE]
obs <- obs_all

expr_gxc <- assay(sce, "logcounts")              # genes x cells (typical)
expr_gxc <- expr_gxc[, cells, drop = FALSE]
expr_cxg <- Matrix::t(expr_gxc)                  # cells x genes

var <- data.frame(symbol = rownames(sce), stringsAsFactors = FALSE)
rownames(var) <- var$symbol

gene_panel <- c("MS4A1", "CD3D", "NKG7")
gene_panel <- intersect(gene_panel, rownames(var))

out_dir <- file.path(getwd(), "exports", "cellucid_from_sce")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  var = var,
  gene_expression = expr_cxg,
  gene_identifiers = gene_panel,
  X_umap_2d = umap2,
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  var_quantization = 8,
  obs_continuous_quantization = 8,
  centroid_min_points = 10
)
```

```{warning}
If your SCE assay is delayed/HDF5-backed, `assay(sce, "logcounts")` might not be a `Matrix` object.
If export fails, switch to “export without gene expression first” and then follow the troubleshooting section.
```

---

## Option C — You don’t have an object (test with toy data)

Use the toy quickstart:
- {doc}`../a_landing_pages/04_quick_start_3_levels`

---

## Part 4 — Load the export folder in the Cellucid web app

Follow the click-by-click guide:
- {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`

**What you should try first in the web app**

1) Confirm you see points.
2) Color by a categorical field (e.g., clusters).
3) Color by a continuous field (e.g., nCount_RNA).
4) If you exported gene expression:
   - search for a gene from your panel
   - color by that gene

---

## Massive troubleshooting (beginner-oriented)

### Symptom: “I don’t see any points / the canvas is blank”

**Likely causes**
- The dataset didn’t load (wrong folder selected).
- Your embedding file is missing or corrupt.
- Browser/WebGL issues.

**How to confirm**
- Check you can see `dataset_identity.json` in your export folder.
- Try a different browser (Chrome/Firefox).

**Fix**
- Re-export with `force=TRUE`.
- Ensure you selected the *export folder* (the one containing `dataset_identity.json`).
- See: {doc}`../i_troubleshooting_index/04_web_app_loading_issues`

### Symptom: “Export worked but colors/fields don’t make sense”

**Likely cause**
- Row order mismatch (metadata not aligned to embedding).

**Fix**
- Use the recipe pages and ensure you reorder everything using one `cells <- ...` vector.

### Symptom: “Export folder is huge / took forever”

**Likely cause**
- You exported too many genes for your cell count.

**Fix**
- Export a gene panel (`gene_identifiers=...`) and use `var_quantization=8`.

More:
- {doc}`../c_data_preparation_api/10_performance_tuning_prepare_export`
- {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
