# Troubleshooting: Prepare/Export

This page is the deep troubleshooting guide for `cellucid_prepare()` / `prepare()`.

If you want a shorter index by topic, see {doc}`../i_troubleshooting_index/index`.

---

## Debug checklist (do this first)

Before reading long troubleshooting sections, run these checks and write down the results:

### 1) Confirm required inputs exist and have the right shapes

```r
stopifnot(is.data.frame(obs))
stopifnot(is.matrix(latent_space) || is.data.frame(latent_space))

stopifnot(is.matrix(X_umap_2d) && ncol(X_umap_2d) == 2)

n_cells <- nrow(X_umap_2d)
stopifnot(nrow(obs) == n_cells)
stopifnot(nrow(latent_space) == n_cells)

if (!is.null(gene_expression)) {
  stopifnot(nrow(gene_expression) == n_cells)
  stopifnot(!is.null(var))
  stopifnot(nrow(var) == ncol(gene_expression))
}

if (!is.null(connectivities)) {
  stopifnot(nrow(connectivities) == n_cells)
  stopifnot(ncol(connectivities) == n_cells)
}
```

### 2) If you have IDs, confirm row order alignment

```r
if (!is.null(rownames(obs)) && !is.null(rownames(X_umap_2d))) {
  stopifnot(identical(rownames(obs), rownames(X_umap_2d)))
}
if (!is.null(rownames(obs)) && !is.null(rownames(latent_space))) {
  stopifnot(identical(rownames(obs), rownames(latent_space)))
}
```

### 3) Count missing/invalid values (quick sanity)

```r
any_na <- function(x) sum(is.na(x))
any_inf <- function(x) sum(is.infinite(x))

cat("NA in UMAP:", any_na(X_umap_2d), " Inf:", any_inf(X_umap_2d), "\n")
cat("NA in latent:", any_na(latent_space), " Inf:", any_inf(latent_space), "\n")
```

Embeddings should not contain `NA`/`Inf`.

---

## Symptom: “At least one embedding must be provided”

**Likely causes**
- You forgot to pass `X_umap_1d`, `X_umap_2d`, and `X_umap_3d`.
- Your embedding object is `NULL` because extraction failed (common with Seurat/SCE when the reduction name is wrong).

**How to confirm**
- Print `str(umap2)` and ensure it’s a numeric matrix with `ncol == 2`.

**Fix**
- Pass at least one embedding:
  ```r
  cellucid_prepare(..., X_umap_2d = umap2)
  ```

---

## Symptom: “X_umap_2d must have exactly 2 columns”

**Likely causes**
- You passed a vector instead of a matrix.
- Your embedding matrix has the wrong dimension (e.g., 3 columns).

**How to confirm**
```r
dim(X_umap_2d)
```

**Fix**
- Use the right embedding or subset columns:
  ```r
  X_umap_2d <- X_umap[, 1:2, drop = FALSE]
  ```

---

## Symptom: “latent_space is required for outlier quantile calculation”

**Likely cause**
- You passed `latent_space = NULL`.

**Fix**
- Provide a latent space matrix with `n_cells` rows.

Practical choices:
- PCA coordinates (recommended if available)
- the same matrix as your embedding (acceptable if you don’t have PCA)

---

## Symptom: “obs data.frame is required”

**Likely cause**
- You passed `obs = NULL`.

**Fix**
- Provide a cell metadata `data.frame` with `n_cells` rows.

If you truly have no metadata, create a minimal one:

```r
obs <- data.frame(dummy = rep("all", n_cells))
```

---

## Symptom: “obs has X rows, but embeddings have Y cells”

**Likely causes**
- You subset `obs` but not the embeddings (or vice versa).
- You used different cell orderings for different inputs.

**How to confirm**
- If you have cell IDs, compare:
  - `rownames(obs)`
  - `rownames(umap2)`

**Fix**
- Reorder and/or subset explicitly using cell IDs.

---

## Symptom: “gene_expression has X cells, but embeddings have Y cells”

**Likely causes**
- You passed expression in genes × cells orientation.
- You subsetted cells inconsistently.

**How to confirm**
```r
dim(gene_expression)
```

**Fix**
- Ensure the matrix is cells × genes.
  - For Seurat/SCE, transpose with `Matrix::t(...)`.

---

## Symptom: “var has X rows, but gene_expression has Y genes”

**Likely causes**
- `var` and `gene_expression` refer to different gene sets.
- Expression was transposed but var was not adjusted.

**Fix**
- Ensure `nrow(var) == ncol(gene_expression)`.
- Ensure the gene identifiers you want are in `rownames(var)` or your chosen `var_gene_id_column`.

---

## Symptom: “var_gene_id_column '...' not found in var”

**Likely cause**
- The column name is wrong.

**How to confirm**
```r
colnames(var)
```

**Fix**
- Use a real column name, or set `var_gene_id_column = "index"` and populate `rownames(var)`.

---

## Symptom: “Matrix package is required to export sparse gene_expression objects”

**Likely cause**
- You passed a sparse matrix but `Matrix` is not installed.

**Fix**
```r
install.packages("Matrix")
```

---

## Symptom: “Matrix package is required to export connectivity matrices”

**Likely cause**
- You provided `connectivities=...` without `Matrix` installed.

**Fix**
```r
install.packages("Matrix")
```

---

## Symptom: “Connectivity matrix shape (...) does not match number of cells”

**Likely causes**
- Connectivity matrix corresponds to a different subset of cells.
- Row/col names are not aligned to your export cell order.

**How to confirm**
- Compare `dim(connectivities)` to `nrow(X_umap_2d)`.
- If `rownames(connectivities)` exist, confirm they match your cell IDs.

**Fix**
- Reorder/subset the matrix:
  ```r
  connectivities <- connectivities[cell_ids, cell_ids]
  ```

---

## Symptom: “Field 'X' has N categories, but uint8 can only hold 254”

**Likely cause**
- You forced `obs_categorical_dtype = "uint8"` and have >254 categories.

**Fix**
- Use `obs_categorical_dtype = "auto"` (default) or `"uint16"`.

---

## Symptom: “Vector field id '...' contains unsupported characters”

**Likely cause**
- Your `vector_fields` list key contains spaces, slashes, or starts/ends with `_`/`.`.

**Fix**
- Rename the key to a filesystem-safe identifier:
  ```r
  names(vector_fields) <- gsub("[^A-Za-z0-9._-]+", "_", names(vector_fields))
  ```

Better: choose explicit safe names up front (see {doc}`08_vector_fields_velocity_displacement`).

---

## Symptom: “I re-exported but nothing changed”

**Likely cause**
- `force = FALSE` (default) and manifests already existed, so export skipped writing.

**How to confirm**
- Check file timestamps in the export folder.

**Fix**
- Re-run with `force = TRUE`, or export to a fresh `out_dir`.

---

## Symptom: export folder is gigantic / takes forever

**Likely causes**
- You exported too many genes (common).
- You exported float32 instead of quantized values.
- Your filesystem is slow (network drive, cloud-synced folder).

**How to confirm**
- Count files:
  ```r
  length(list.files(file.path(out_dir, "var")))
  ```
- Estimate size:
  - on macOS/Linux: `system(paste("du -sh", shQuote(out_dir)))`

**Fix (ordered)**
1) Export fewer genes (`gene_identifiers=...`)
2) Set `var_quantization=8`
3) Set `compression=6`
4) Export to a local SSD

---

## Symptom: “My gene names are weird numbers”

**Likely cause**
- `rownames(var)` were `NULL`, so gene IDs defaulted to `"0".."n_genes-1"`.

**Fix**
- Set `rownames(var)` explicitly, or set `var_gene_id_column` to the right column.

---

## Symptom: “Some genes/fields disappear or overwrite each other”

**Likely cause**
- Filename sanitization collisions (two keys map to the same `<safe_key>`).

**How to confirm**
- Compare original IDs vs sanitized IDs and look for duplicates.

**Fix**
- Rename genes/fields to be unique after sanitization, or export a smaller subset.

---

## If you still can’t resolve it

Gather these artifacts before asking for help (it speeds everything up):

- The exact R call to `cellucid_prepare(...)` (including args)
- Output of:
  - `sessionInfo()`
  - `dim(...)` for embeddings/latent/obs/expression/connectivities
  - `list.files(out_dir, recursive = TRUE)`
- The first ~50 lines of:
  - `dataset_identity.json`
  - `obs_manifest.json`
