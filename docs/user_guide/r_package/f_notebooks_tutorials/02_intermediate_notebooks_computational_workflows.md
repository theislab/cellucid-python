# Intermediate Notebooks (Computational Workflows)

**Audience:** computational users  
**Time:** 30–90 minutes  
**What you’ll learn**
- How to export *reproducibly* (stable dataset IDs, curated fields, versioning)
- How to keep exports small enough to share (gene panels + quantization + compression)
- How to include connectivities when available
- How to validate an export before sending it to collaborators

This tutorial is “intermediate” because it assumes:
- you can navigate Seurat/SCE structures,
- you care about data curation and reproducibility,
- and you want predictable behavior when sharing exports.

---

## Part 0 — A reproducible export plan (do this before coding)

### Decide:

1) **Dataset identity**
   - `dataset_id`: stable machine-friendly ID (recommended)
   - `dataset_name`: human-friendly name
2) **What cells are included**
   - whole dataset or a subset?
3) **What metadata is exported**
   - keep it focused (the UI gets noisy with 200 columns)
4) **Gene export strategy**
   - all genes (often too big)
   - a gene panel (recommended)
   - or “no expression” (fastest)
5) **Performance knobs**
   - `var_quantization=8`
   - `obs_continuous_quantization=8`
   - `compression=6`

If you have not read it yet, dataset identity matters: {doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility`.

---

## Part 1 — Build a gene panel (recommended)

Start with something practical:

- marker genes you expect to check visually
- HVGs
- pathway gene sets

Example (toy list):

```r
gene_panel <- c("MS4A1", "CD3D", "NKG7", "LYZ", "PPBP")
```

You will intersect this with your object’s gene IDs later.

---

## Part 2 — Export from Seurat (curated, shareable)

This is a compact version of the full Seurat recipe, focused on reproducible exports.

```r
library(cellucid)
library(Seurat)
library(Matrix)

cells <- colnames(seu)

# Embedding + latent
umap2 <- Embeddings(seu, "umap")[cells, 1:2, drop = FALSE]
latent_space <- Embeddings(seu, "pca")[cells, , drop = FALSE]

# Curated obs
obs_all <- seu@meta.data[cells, , drop = FALSE]
obs_keys <- intersect(c("seurat_clusters", "orig.ident", "nCount_RNA", "nFeature_RNA"), colnames(obs_all))
obs <- obs_all[, obs_keys, drop = FALSE]
if ("seurat_clusters" %in% colnames(obs)) obs$seurat_clusters <- factor(obs$seurat_clusters)

# Expression (optional; skip if too big)
expr_gxc <- GetAssayData(seu, slot = "data")[, cells, drop = FALSE]   # genes x cells
expr_cxg <- Matrix::t(expr_gxc)                                       # cells x genes

var <- data.frame(symbol = rownames(expr_gxc), stringsAsFactors = FALSE)
rownames(var) <- var$symbol

gene_panel <- intersect(gene_panel, rownames(var))

# Optional connectivities
conn <- NULL
if ("RNA_snn" %in% names(seu@graphs)) {
  conn <- seu@graphs$RNA_snn
  conn <- conn[cells, cells]
}

out_dir <- file.path(getwd(), "exports", "seurat_shareable_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  var = var,
  gene_expression = expr_cxg,
  gene_identifiers = gene_panel,
  X_umap_2d = umap2,
  connectivities = conn,
  out_dir = out_dir,
  force = TRUE,
  compression = 6,
  var_quantization = 8,
  obs_continuous_quantization = 8,
  dataset_id = "my_project_pbmc_v1",
  dataset_name = "My Project PBMC (v1)",
  dataset_description = "Curated export for collaboration.",
  source_name = "MyLab",
  source_url = "https://example.org",
  source_citation = "Add your paper/preprint citation here."
)
```

### Validate before sharing

```r
ident <- jsonlite::read_json(file.path(out_dir, "dataset_identity.json"), simplifyVector = TRUE)
ident$stats
```

Confirm:
- `n_cells` matches
- `n_genes` equals `length(gene_panel)`
- `has_connectivity` is correct

---

## Part 3 — Export from SingleCellExperiment (curated, shareable)

```r
library(cellucid)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

cells <- colnames(sce)

umap2 <- reducedDim(sce, "UMAP")[cells, 1:2, drop = FALSE]
latent_space <- reducedDim(sce, "PCA")[cells, , drop = FALSE]

obs_all <- as.data.frame(colData(sce))[cells, , drop = FALSE]
obs <- obs_all

expr_gxc <- assay(sce, "logcounts")[, cells, drop = FALSE]  # genes x cells
expr_cxg <- Matrix::t(expr_gxc)

var <- data.frame(symbol = rownames(sce), stringsAsFactors = FALSE)
rownames(var) <- var$symbol

gene_panel <- intersect(gene_panel, rownames(var))

out_dir <- file.path(getwd(), "exports", "sce_shareable_v1")
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
  dataset_id = "my_project_sce_v1",
  dataset_name = "My Project SCE (v1)"
)
```

```{warning}
If your assay is delayed/HDF5-backed, you may need to realize/convert it. For large datasets, prefer exporting a gene panel and/or skipping expression entirely.
```

---

## Part 4 — Hosting and sharing

After export, you can:
- load locally via file picker, or
- host via a static server (better for sharing)

Start here:
- {doc}`../d_viewing_loading/02_host_exports_for_sharing`

---

## Massive troubleshooting (intermediate)

### Symptom: “Export loads, but gene coloring looks washed out”

**Likely cause**
- A few extreme expression values stretched the min/max scaling.

**Confirm**
- Inspect the gene’s value range before export (summary/quantiles).

**Fix**
- Export log-normalized expression rather than raw counts.
- Clip extreme values before export (advanced; do carefully).

### Symptom: “The export is too big to share”

**Fixes (ordered)**
1) Export fewer genes (gene panel)
2) Quantize expression (`var_quantization=8`)
3) Compress (`compression=6`)
4) Export fewer cells (subsample or subset to a cell type)

### Symptom: “My collaborator sees different cells/labels than I do”

**Likely cause**
- You overwrote an export folder or changed cell order without versioning.

**Fix**
- Version your export folders and dataset IDs. See {doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility`.
