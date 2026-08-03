# Embeddings and Coordinates

**Audience:** everyone exporting (especially if your dataset “looks wrong” in the viewer)  
**Time:** 10 minutes  
**Goal:** understand what embeddings are required and how Cellucid normalizes them.

Cellucid visualizes cells as points. The positions of those points come from one or more **embeddings**.

In `cellucid-r`, the embedding arguments are:
- `X_umap_1d` (shape `(n_cells, 1)`, or a plain length-`n_cells` numeric vector)
- `X_umap_2d` (shape `(n_cells, 2)`)
- `X_umap_3d` (shape `(n_cells, 3)`)

At least one must be provided. Only `X_umap_1d` accepts a plain vector, as the
single column it declares; `X_umap_2d` and `X_umap_3d` require the full matrix.

## What “embedding” means (practical)

An embedding is any numeric coordinate matrix with:
- rows = cells
- columns = 1/2/3 coordinates

Despite the argument names, the values do not have to come from UMAP specifically.
You can export tSNE, PCA, diffusion maps, etc., as long as the shape is correct.

## Embedding validation rules

`cellucid-r` enforces:
- each embedding must be a 2D matrix-like object — a base R `matrix`, a
  `data.frame` whose every column is a plain numeric vector, or a
  `Matrix::Matrix`,
- number of columns must match the dimension (`2D` → exactly 2 columns),
- all provided embeddings must have the same number of rows (`n_cells`),
- `n_cells` must be between 1 and 2,147,483,647.

It rejects `NA`, `NaN`, infinities, nonnumeric matrices, and embeddings without
a positive finite coordinate range before publication.

:::{note}
A `data.frame` column that carries a class — a `units` column, a `bit64::integer64`
column, or any S3 class of your own — is refused rather than silently unclassed,
because `as.matrix()` would drop the attribute that says what the numbers mean.
Convert deliberately before the call.
:::

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

A zero range has no substitute. If every point shares one coordinate on every
axis, `max(axis_maxs - axis_mins)` is `0`, and rather than invent a scale nobody
chose, `cellucid_prepare()` stops before writing anything:

```text
Embedding coordinates must span a positive finite coordinate range.
```

A range that is merely *small* is still positive, so it is accepted and scaled
up like any other.

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

An embedding in which every point sits at the same coordinate is **rejected**,
not exported: `Embedding coordinates must span a positive finite coordinate
range.` It usually means the extraction returned the wrong object, or the
reduction was never computed.

A *nearly* constant embedding is a different case. Any positive range is scaled
up to span ~2 units, so it exports normally, and what you see in the viewer —
one tight clump — is the separation your embedding actually has.

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
