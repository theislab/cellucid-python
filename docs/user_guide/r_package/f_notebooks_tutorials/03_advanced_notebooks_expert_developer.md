# Advanced Notebooks (Expert / Developer)

**Audience:** expert computational users, developers/maintainers  
**Time:** 60–120 minutes  
**What you’ll learn**
- How to reason about the export format at the file level
- How to validate binaries (points, edges, values) and spot corruption/misalignment
- How vector fields are named and scaled
- How exact connectivity validation and index dtype selection work

This tutorial assumes you’ve already successfully exported at least one dataset.

---

## Part 1 — Export with “everything on” (small synthetic example)

We use a tiny synthetic dataset so you can inspect files safely.

```r
library(cellucid)
library(Matrix)

n_cells <- 4

X_umap_2d <- matrix(c(0, 0,
                      1, 0,
                      0, 1,
                      1, 1),
                    ncol = 2, byrow = TRUE)

latent_space <- matrix(rnorm(n_cells * 5), ncol = 5)

obs <- data.frame(
  cluster = factor(c("A", "A", "B", "B")),
  score = c(0.1, 0.2, 0.3, 0.4)
)

gene_expression <- matrix(rexp(n_cells * 3), nrow = n_cells, ncol = 3)
var <- data.frame(symbol = c("G1", "G2", "G3"), stringsAsFactors = FALSE)
rownames(var) <- var$symbol

conn <- Matrix(0, nrow = n_cells, ncol = n_cells, sparse = TRUE)
conn[1, 2] <- 1
conn[2, 1] <- 1
conn[1, 3] <- 1
conn[3, 1] <- 1
conn[4, 3] <- 1
conn[3, 4] <- 1

vector_fields <- list(
  velocity_umap_2d = matrix(c(0.2, 0,
                              0.2, 0,
                              0.2, 0,
                              0.2, 0),
                            ncol = 2, byrow = TRUE)
)

out_dir <- file.path(tempdir(), "cellucid_advanced_debug")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cellucid_prepare(
  latent_space = latent_space,
  obs = obs,
  var = var,
  gene_expression = gene_expression,
  connectivities = conn,
  vector_fields = vector_fields,
  X_umap_2d = X_umap_2d,
  obs_categorical_dtype = "uint8",
  dataset_name = "Cellucid Advanced Debug",
  dataset_id = "cellucid_advanced_debug",
  out_dir = out_dir,
  force = TRUE,
  compression = NULL,
  var_quantization = NULL,
  obs_continuous_quantization = NULL,
  centroid_min_points = 1
)

list.files(out_dir, recursive = TRUE)
```

---

## Part 2 — Validate embeddings (`points_2d.bin`)

`points_2d.bin` is float32 row-major.

```r
path <- file.path(out_dir, "points_2d.bin")
con <- file(path, open = "rb")
on.exit(close(con), add = TRUE)

vals <- readBin(con, what = "numeric", size = 4, endian = "little", n = n_cells * 2)
coords <- matrix(vals, ncol = 2, byrow = TRUE)
coords
```

Sanity expectations:
- values should be finite
- coordinates should be centered around ~0
- range should be roughly within [-1, 1] (after normalization)

Normalization details: {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`

---

## Part 3 — Validate obs encoding (codes + outliers)

Inspect `obs_manifest.json`:

```r
obs_manifest <- jsonlite::read_json(file.path(out_dir, "obs_manifest.json"), simplifyVector = TRUE)
str(obs_manifest, max.level = 3)
```

Look for:
- `_continuousFields` and `_categoricalFields`; every entry begins with its
  payload index, then its key
- file patterns under `_obsSchemas`, which are templates over `{index}`

Read categorical codes. The filename is the field's payload index — take it from
the manifest entry rather than guessing from the key:

```r
field <- obs_manifest$`_categoricalFields`[[1]]
codes_path <- file.path(out_dir, "obs", sprintf("%d.codes.u8", field[[1]]))
con <- file(codes_path, open = "rb")
on.exit(close(con), add = TRUE)
codes <- readBin(con, what = "integer", size = 1, n = n_cells, endian = "little")
codes
```

Remember:
- codes are 0-based
- missing is 255 (uint8) or 65535 (uint16)

---

## Part 4 — Validate gene values

Read one gene file:

```r
g1_path <- file.path(out_dir, "var", "G1.values.f32")
con <- file(g1_path, open = "rb")
on.exit(close(con), add = TRUE)
g1 <- readBin(con, what = "numeric", size = 4, n = n_cells, endian = "little")
g1
```

If you exported with quantization, the file dtype changes and you must read as bytes/uint16 instead.

Gene export details: {doc}`../c_data_preparation_api/06_gene_expression_matrix`

---

## Part 5 — Validate connectivity edge pairs

Connectivity export writes three aligned arrays:
- `edges.src.bin`
- `edges.dst.bin`
- `edges.weights.f64.bin`

The dtype is `uint16` through 65,536 cells and `uint32` through
4,294,967,296 cells. Larger axes are outside the current graph contract.

For small `n_cells` it will be uint16:

```r
src_path <- file.path(out_dir, "connectivity", "edges.src.bin")
dst_path <- file.path(out_dir, "connectivity", "edges.dst.bin")
weights_path <- file.path(out_dir, "connectivity", "edges.weights.f64.bin")

read_u16 <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "integer", size = 2, endian = "little", n = 1000)
}

src <- read_u16(src_path)
dst <- read_u16(dst_path)

read_f64 <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "numeric", size = 8, endian = "little", n = 1000)
}

weights <- read_f64(weights_path)
cbind(src[1:3], dst[1:3], weights[1:3])
```

Remember:
- indices are 0-based
- only unique undirected edges are kept (`src < dst`)
- weights are exact little-endian Float64 values aligned by array index
- the input must already be finite, non-negative, exactly symmetric in
  topology and weight, and zero on the diagonal; invalid input is rejected
  rather than transformed

---

## Part 6 — Vector fields: naming + scaling

Vector binaries live in `vectors/`, named by the field's payload index. The
`vector_fields` section of `dataset_identity.json` records the exact path for
each field and dimension.

```r
vec_path <- file.path(out_dir, "vectors", "0_2d.bin")
con <- file(vec_path, open = "rb")
on.exit(close(con), add = TRUE)
vec_vals <- readBin(con, what = "numeric", size = 4, endian = "little", n = n_cells * 2)
vec <- matrix(vec_vals, ncol = 2, byrow = TRUE)
vec
```

Vectors are automatically scaled by the embedding normalization factor. This often surprises people when they compare exported values to original velocity units.

Details: {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`

---

## Part 7 — “Export format bugs” that are actually user bugs

### 1) Row order mismatches

Symptoms:
- clusters don’t match the point cloud
- gene expression appears on the “wrong” cells

Root cause:
- embeddings/obs/expression were not aligned to the same cell order

Prevention:
- choose a canonical `cells` vector and reorder everything explicitly

### 2) Duplicate identifiers

Symptoms:
- missing genes
- a lookup by gene ID or obs key resolving the wrong row

Root cause:
- two rows of `var`, or two entries of `obs_keys`, carry the same identifier

Prevention:
- check `anyDuplicated(ids) == 0` before export. Visible character content is
  not a problem — identifiers never become filenames, so nothing needs escaping
  — but invisible ones are, so check `identical(ids, trimws(ids))` too.

---

## Next steps

- Full format spec: {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
- Export troubleshooting: {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
