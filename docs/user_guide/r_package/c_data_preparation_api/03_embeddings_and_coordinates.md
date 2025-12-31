# Embeddings and Coordinates

**Audience:** everyone exporting (especially if your dataset “looks wrong” in the viewer)  
**Time:** 10 minutes  
**Goal:** understand what embeddings are required and how Cellucid normalizes them.

Cellucid visualizes cells as points. The positions of those points come from one or more **embeddings**.

In `cellucid-r`, the embedding arguments are:
- `X_umap_1d` (shape `(n_cells, 1)`)
- `X_umap_2d` (shape `(n_cells, 2)`)
- `X_umap_3d` (shape `(n_cells, 3)`)

At least one must be provided.

```{warning}
`X_umap_4d` is reserved for future work and must be `NULL`. Passing it raises an error.
```

## What “embedding” means (practical)

An embedding is any numeric coordinate matrix with:
- rows = cells
- columns = 1/2/3 coordinates

Despite the argument names, the values do not have to come from UMAP specifically.
You can export tSNE, PCA, diffusion maps, etc., as long as the shape is correct.

## Embedding validation rules

`cellucid-r` enforces:
- each embedding must be a 2D matrix-like object,
- number of columns must match the dimension (`2D` → exactly 2 columns),
- all provided embeddings must have the same number of rows (`n_cells`).

It does **not** currently sanitize `NA`/`Inf` in embeddings. Avoid missing/invalid coordinates.

## Embedding normalization (important!)

Before writing `points_*d.bin`, `cellucid-r` normalizes each embedding:

1) Compute per-axis min/max.
2) Find the maximum axis range across axes.
3) Center the embedding at the midpoint of the min/max box.
4) Scale so the max range spans ~2 units.

In code terms (matches `cellucid-r` implementation):

- `center = (axis_mins + axis_maxs) / 2`
- `scale_factor = 2 / max(axis_maxs - axis_mins)`
- `coords_normalized = (coords - center) * scale_factor`

If the embedding is (nearly) constant (`max_range < 1e-8`), the exporter uses `max_range = 1` to avoid divide-by-zero.

### Why normalize?

Normalization makes:
- camera defaults behave consistently,
- vector fields scale correctly (they are scaled with the same factor),
- and exports more comparable across datasets.

## File encoding details (for debugging)

`points_2d.bin` contains:
- float32 values
- little-endian
- row-major order

“Row-major” means:
- cell 1’s coordinates are stored first,
- then cell 2’s, etc.

You can sanity-check a tiny export by reading the binary:

```r
path <- file.path(out_dir, "points_2d.bin")
con <- file(path, open = "rb")
on.exit(close(con), add = TRUE)

vals <- readBin(con, what = "numeric", size = 4, endian = "little", n = n_cells * 2)
coords <- matrix(vals, ncol = 2, byrow = TRUE)
head(coords)
```

If you exported with compression, the file is `points_2d.bin.gz` and you should open with `gzfile(...)`.

## Edge cases

### Constant embeddings (all points identical)

If all points are the same (or nearly the same), normalization will produce near-zero coordinates.
In the viewer, this looks like “everything is on top of itself”.

This is not a Cellucid bug; the embedding contains no separation.

### Mixed 2D + 3D embeddings

You can export both `X_umap_2d` and `X_umap_3d`. The dataset identity file records:
- `available_dimensions` (e.g. `[2, 3]`)
- `default_dimension` (prefers 3D if available)

### “My dataset is mirrored / rotated”

Normalization does not rotate embeddings; it only centers and scales.
If a dataset looks rotated, that’s usually just a different embedding convention.

## Troubleshooting pointers

- Export errors about “must have exactly N columns” → check your embedding matrix shape.
- Viewer loads but points appear “squashed” → likely a coordinate issue (range extremes/outliers).
- Vector fields look wrong → see {doc}`08_vector_fields_velocity_displacement` (they are scaled by the normalization factor).
