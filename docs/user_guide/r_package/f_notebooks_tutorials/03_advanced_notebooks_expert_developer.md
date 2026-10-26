# Advanced Tutorial (Expert / Developer)

**Audience:** computational biologists and developers
**Time:** 30 minutes
**Goal:** store vector fields and a neighbour graph, then read the store back
byte for byte to see exactly what `cellucid_prepare()` wrote.

## Step 1 — A trajectory with velocities and a graph

```r
library(cellucid)
set.seed(3)

n <- 1000
t <- sort(runif(n, 0, 2 * pi))
umap <- cbind(cos(t), sin(t)) * 5 + matrix(rnorm(2 * n, sd = 0.2), ncol = 2)
umap3 <- cbind(umap, t)
velocity_2d <- cbind(-sin(t), cos(t))          # along the circle
velocity_3d <- cbind(velocity_2d, 1)

# A symmetric 10-nearest-neighbour graph with a zero diagonal.
d <- as.matrix(dist(umap))
diag(d) <- Inf
nearest <- t(apply(d, 1, order))[, 1:10]
knn <- Matrix::sparseMatrix(i = rep(seq_len(n), 10), j = as.vector(nearest),
                            x = 1, dims = c(n, n))
knn <- (knn + Matrix::t(knn)) / 2

obs <- data.frame(pseudotime = t, phase = cut(t, 4, labels = c("G1", "S", "G2", "M")))
```

## Step 2 — Write the store

```r
store <- cellucid_prepare(
  out_path = file.path(tempdir(), "advanced.zarr"),
  format = "zarr",
  embeddings = list(X_umap = umap, X_umap_3d = umap3),
  dataset_id = "advanced-tutorial",
  dataset_name = "Advanced tutorial",
  obs = obs,
  connectivities = knn,
  vector_fields = list(velocity = list(X_umap = velocity_2d, X_umap_3d = velocity_3d)),
  vector_field_default = "velocity",
  embedding_defaults = c("2d" = "X_umap", "3d" = "X_umap_3d"),
  force = TRUE
)
sort(list.files(file.path(store, "obsm")))
```

The field `velocity` is stored under scVelo's names for it,
`obsm['velocity_umap']` and `obsm['velocity_umap_3d']`. The viewer lists it as
`velocity` and draws it whenever a view shows `X_umap` or `X_umap_3d`. A vector
field can be bound to any embedding named `X_<basis>`, with as many columns as
it; see {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Step 3 — Read an array back by hand

Every array of a Zarr v2 store (the default) has a `.zarray` file describing
it, and one file per chunk, compressed with zlib, which base R decompresses. An
array of values wider than one byte is byte-shuffled first (its `.zarray`
lists the `shuffle` filter): the first byte of every value, then the second,
and so on; undo it before reading the values.

```r
read_f4_chunk <- function(store, key, chunk) {
  meta <- jsonlite::read_json(file.path(store, key, ".zarray"))
  path <- file.path(store, key, chunk)
  bytes <- memDecompress(readBin(path, "raw", file.size(path)), type = "gzip")
  if (length(meta$filters) > 0 && identical(meta$filters[[1]]$id, "shuffle")) {
    bytes <- as.vector(t(matrix(bytes, ncol = 4)))   # undo the byte shuffle
  }
  values <- readBin(bytes, "double", n = prod(unlist(meta$chunks)), size = 4, endian = "little")
  matrix(values, ncol = meta$chunks[[2]], byrow = TRUE)   # C order: row-major
}

meta <- jsonlite::read_json(file.path(store, "obsm", "X_umap", ".zarray"))
str(meta[c("shape", "chunks", "dtype", "compressor", "filters")])

stored <- read_f4_chunk(store, "obsm/X_umap", "0.0")
float32 <- function(x) readBin(writeBin(as.double(x), raw(), size = 4), "double", n = length(x), size = 4)
identical(stored[seq_len(n), ], matrix(float32(umap), ncol = 2))
```

The embedding is stored as float32, at its own coordinates; the viewer scales it
to the view (`normalize_embeddings = FALSE` turns that off).

## Step 4 — The graph

```r
graph_meta <- jsonlite::read_json(file.path(store, "obsp", "connectivities", ".zattrs"))
graph_meta$`encoding-type`
length(list.files(file.path(store, "obsp", "connectivities", "data")))
```

`obsp['connectivities']` is a float32 `csr_matrix`. For a symmetric matrix, the
column-compressed arrays of a `dgCMatrix` are also its row-compressed arrays, so
`cellucid_prepare()` writes them as they are.

## Step 5 — The same store as one file, and read in Python

```r
archive <- cellucid_prepare(
  out_path = file.path(tempdir(), "advanced.zarr.zip"),
  format = "zarr.zip",
  embeddings = list(X_umap = umap, X_umap_3d = umap3),
  dataset_id = "advanced-tutorial",
  dataset_name = "Advanced tutorial",
  obs = obs,
  connectivities = knn,
  vector_fields = list(velocity = list(X_umap = velocity_2d, X_umap_3d = velocity_3d)),
  force = TRUE
)
head(utils::unzip(archive, list = TRUE)$Name)
```

```python
import anndata, zarr
adata = anndata.read_zarr(zarr.storage.ZipStore("advanced.zarr.zip", mode="r"))
adata.obsm["velocity_umap"].shape, adata.obsp["connectivities"].format
```

The full layout of the store is in
{doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.
