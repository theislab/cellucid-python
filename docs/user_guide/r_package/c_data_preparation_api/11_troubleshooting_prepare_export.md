# Troubleshooting: Prepare/Export

**Audience:** everyone whose export just stopped  
**Time:** find your symptom in under a minute  
**What you'll get:** the cause behind each `cellucid_prepare()` error, and the fix

This page is the deep troubleshooting guide for `cellucid_prepare()`.

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

## Symptom: “cellucid_prepare() is missing 1 required argument: embeddings”

**Likely causes**
- You did not pass `embeddings`.
- Your embedding object is `NULL` because extraction failed (common with Seurat/SCE when the reduction name is wrong).

**How to confirm**
- Print `str(umap2)` and ensure it’s a numeric matrix with `ncol == 2`.

**Fix**
- Pass at least one embedding:
  ```r
  cellucid_prepare(..., embeddings = list(X_umap = umap2))
  ```

---

## Symptom: “Embedding '…' has 10 columns, and an embedding plots 1, 2, or 3 of them”

**Likely causes**
- You passed a wide reduction (a PCA, a 10-component UMAP) as a plain matrix.

**How to confirm**
```r
dim(umap)
```

**Fix**
- Name the 1-based columns to plot:
  ```r
  cellucid_prepare(..., embeddings = list(X_pca = cellucid_embedding(pca, columns = c(1, 2))))
  ```

---

## Symptom: “… so embedding_defaults must name the one shown first”

**Likely cause**
- Two or more embeddings share a dimension.

**Fix**
- Add `embedding_defaults = c("2d" = "X_umap")` for the dimension the message names.

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

## Symptom: “obs must be a data.frame.”

**Likely cause**
- You passed `obs = NULL`, or a matrix / list where a `data.frame` was expected.

**Fix**
- Provide a cell metadata `data.frame` with `n_cells` rows.

If you truly have no metadata, create a minimal one:

```r
obs <- data.frame(cell_group = factor(rep("all_cells", n_cells)))
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
- A repeated gene name is exported as `NAME-1`, `NAME-2`, …; request that unique
  name in `gene_identifiers`.

---

## Symptom: “var_gene_id_column '...' not found in var”

**Likely cause**
- The column name is wrong.

**How to confirm**
```r
colnames(var)
```

**Fix**
- Use a real column name, or set `var_gene_id_column = NULL` and populate
  `rownames(var)`.

---

## Symptom: “Matrix package is required to export sparse gene_expression objects”

**Likely cause**
- You passed a sparse matrix but `Matrix` is not installed.

**Fix**
```r
install.packages("Matrix")
```

---

## Symptom: “Matrix package is required to validate sparse connectivity matrices.”

**Likely cause**
- You provided a sparse `connectivities=...` without `Matrix` installed. A dense
  base R matrix needs no `Matrix` at all.

**Fix**
```r
install.packages("Matrix")
```

---

## Symptom: “Connectivity matrix shape must exactly match the cell axis: expected (…), got (…).”

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

## Symptom: “Field 'X' has N categories, but uint8 supports at most 255. Pass obs_categorical_dtype='uint16'.”

**Likely cause**
- You chose `obs_categorical_dtype = "uint8"` and one field has more than 255
  categories. The top code of each width is reserved for “missing”, so `uint8`
  addresses 255 categories, `uint16` 65,535 and `uint32` 4,294,967,295. The
  message ends with the width that fits.

**How to confirm**
```r
sort(vapply(obs, function(x) length(unique(x)), integer(1)), decreasing = TRUE)[1:5]
```

**Fix**
- Use the width the message names (`"uint16"`, or `"uint32"` above 65,535
  categories), or collapse rare categories. The width
  applies to *every* categorical field in the export, not one field at a time.

---

## Symptom: a `vector_fields` entry is rejected

**Likely causes**
- A field's value is a bare matrix, and the export stops with
  `Vector field 'velocity' must map embedding ids to arrays, for example list(velocity = list(X_umap = vectors)).`
- A field names an embedding you did not export:
  `Vector field 'velocity' is bound to embedding 'X_tsne', which is not one of the exported embeddings: …`
- An array's width is not its embedding's dimension:
  `Vector field 'velocity' on embedding 'X_umap' must have exactly 2 columns, got 3.`

The field id itself carries no character rule beyond the display rule: it names
a field in `dataset_identity.json`, and the payload is written to
`vectors/<field index>_<embedding index>.bin`.

**Fix**
- Bind each array to the embedding it was computed in:
  ```r
  vector_fields <- list(velocity = list(X_umap = velocity_2d))
  ```

See {doc}`08_vector_fields_velocity_displacement`.

---

## Symptom: “Categorical field 'organ' has 2 labels the viewer cannot show as written”

Cause:
- a category label carries a character with no glyph — padding such as
  `"Liver "`, a control character, a zero-width character, or a byte-order mark
  at the front of a column read from a UTF-8 CSV. An empty label (`""`) is
  reported the same way.

Confirm:

```r
labels <- levels(factor(obs$organ))
labels[labels != trimws(labels, whitespace = "[\\h\\v]")]
labels[grepl("[\\x01-\\x1f\\x7f]", labels, perl = TRUE)]
```

Fix:

```r
obs$organ <- factor(trimws(as.character(obs$organ), whitespace = "[\\h\\v]"))
```

Check the result before re-exporting: if the column previously held both
`"Liver"` and `"Liver "`, trimming merges them into one category and moves
cells between them. That may be exactly what you want — but it is a change to
your annotation, which is why `cellucid_prepare()` will not make it for you.

A second message covers two labels that render identically:
`Categorical field 'organ' labels 'T  cell' and 'T cell' are stored as
different categories but are drawn identically`. Rename one of them.

The same rule covers `dataset_name`, `dataset_description`, `source_name`,
`source_url`, and `source_citation`.

---

## Symptom: “I re-exported but nothing changed”

**Likely cause**
- `force = FALSE` (default) and the output generation already existed, so the
  complete candidate was rejected.

**How to confirm**
- Read the exact `out_dir already exists` error; the prior generation remains
  unchanged.

**Fix**
- Re-run with `force = TRUE` for an intentional atomic replacement, or export
  to a fresh `out_dir`.

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

## Symptom: “var has only automatic row names …”

**Likely cause**
- `var_gene_id_column` is `NULL` (the default) and you never set
  `rownames(var)`. R then reports the automatic sequence `"1"` … `"n_genes"`,
  which would name every gene after its own row number, so the export stops
  instead.

**Fix**
- Set `rownames(var)` to the gene identifiers, or pass `var_gene_id_column` to
  name the column that holds them.

---

## Symptom: gene or field identifiers are rejected

**Likely cause**
- an identifier is empty, two observation keys are equal, or an exported
  identifier carries a character with no glyph. There is no filename rule to
  violate: payload files are named by an integer index, so `HLA-DRB1/2`, `CON`,
  and a name with interior spaces all export. Repeated and empty **gene** names
  are not an error: they are made unique, and an empty one takes its row name.

**How to confirm**
- Compare observation keys for duplicates with `anyDuplicated(keys)`.
- Check for invisible characters with `identical(ids, trimws(ids))`.

**Fix**
- Rename the duplicate, or clean the invisible character. The exporter rejects
  the complete candidate before writing any payload rather than trimming, which
  would rewrite an annotation you did not ask to change.

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
