# Advanced Notebooks (Expert / Developer)

**Audience:** expert computational users, developers/maintainers  
**Time:** 60–120 minutes  
**What you’ll learn**
- How to reason about the export format at the file level
- How to validate binaries (points, edges, values) and spot corruption/misalignment
- How vector fields are named and scaled
- How connectivity export works (dtype selection, symmetrization)

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
conn[1, 3] <- 1
conn[4, 3] <- 1

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
- `_continuousFields` and `_categoricalFields`
- file patterns under `_obsSchemas`

Read categorical codes:

```r
codes_path <- file.path(out_dir, "obs", "cluster.codes.u8")
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

Connectivity export writes parallel arrays:
- `edges.src.bin`
- `edges.dst.bin`

The dtype depends on `n_cells` (uint16/uint32/uint64).

For small `n_cells` it will be uint16:

```r
src_path <- file.path(out_dir, "connectivity", "edges.src.bin")
dst_path <- file.path(out_dir, "connectivity", "edges.dst.bin")

read_u16 <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "integer", size = 2, endian = "little", n = 1000)
}

src <- read_u16(src_path)
dst <- read_u16(dst_path)
cbind(src[1:3], dst[1:3])
```

Remember:
- indices are 0-based
- only unique undirected edges are kept (`src < dst`)

---

## Part 6 — Vector fields: naming + scaling

Vector binaries live in `vectors/`.

```r
vec_path <- file.path(out_dir, "vectors", "velocity_umap_2d.bin")
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

### 2) Filename collisions after sanitization

Symptoms:
- missing genes
- overwritten fields

Root cause:
- two IDs map to the same sanitized filename

Prevention:
- enforce uniqueness after sanitization (especially for gene IDs)

---

## Next steps

- Full format spec: {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
- Export troubleshooting: {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
