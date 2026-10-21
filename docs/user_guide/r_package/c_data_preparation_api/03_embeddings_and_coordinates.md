# Embeddings and Coordinates

**Audience:** everyone exporting (especially if your dataset “looks wrong” in the viewer)  
**Time:** 10 minutes  
**Goal:** understand how embeddings are passed and how Cellucid normalizes them.

Cellucid visualizes cells as points. The positions of those points come from one or more **embeddings**.

In `cellucid-r`, every embedding is passed in one argument, `embeddings`: a
named list mapping an embedding id to a coordinate matrix or to a
`cellucid_embedding()`. At least one is required, and **any number may share a
dimension** — a UMAP and a t-SNE in 2D, a 3D UMAP, two principal components.

```r
library(cellucid)

cellucid_prepare(
  latent_space = pca,
  obs = obs,
  embeddings = list(
    X_umap = cellucid_embedding(umap_2d, label = "UMAP", source_key = "umap"),
    X_tsne = cellucid_embedding(tsne_2d, label = "t-SNE", source_key = "tsne"),
    X_pca = cellucid_embedding(pca, columns = c(1, 2), label = "PC 1 / PC 2"),
    X_umap_3d = umap_3d
  ),
  embedding_defaults = c("2d" = "X_umap"),
  out_dir = "exports/my_dataset",
  dataset_name = "My dataset",
  dataset_id = "my-dataset-v1",
  obs_categorical_dtype = "uint16"
)
```

## What “embedding” means (practical)

An embedding is any numeric coordinate matrix with:
- rows = cells
- 1, 2, or 3 plotted columns

The values do not have to come from UMAP specifically. You can export t-SNE,
PCA, diffusion maps, spatial coordinates, and so on. A plain matrix of 1, 2, or
3 columns — or a plain length-`n_cells` numeric vector, for 1D — is plotted
whole, and its column count is its dimension.

## `cellucid_embedding()`: columns, label, provenance

```r
cellucid_embedding(values, columns = NULL, label = NULL, source_key = NULL)
```

- `values`: the coordinate matrix (a base R `matrix`, a `data.frame` whose every
  column is a plain numeric vector, or a `Matrix::Matrix`), or a numeric vector
  for one column.
- `columns`: the **1-based** columns to plot, one per axis, so their count (1,
  2, or 3) is the dimension. Required when `values` has more than three
  columns; otherwise every column is plotted. The export records them 0-based,
  as `cellucid-python` does: `columns = c(4, 2)` is written as `[3, 1]`.
- `label`: the text the viewer shows for this embedding; defaults to the id.
- `source_key`: provenance, such as the reduction name; recorded in
  `dataset_identity.json` and changes nothing else.

Only the selected columns are read, so a non-finite value in an unselected
principal component does not stop the export.

## Which embedding is shown first: `embedding_defaults`

`embedding_defaults` names the embedding each dimension shows first, as a named
character vector or list such as `c("2d" = "X_umap")`. It is required for every
dimension that has more than one embedding; a dimension with exactly one uses it.
The viewer opens in the highest dimension the export carries.

```text
2 2D embeddings were supplied (c("X_tsne", "X_umap")), so embedding_defaults must name the one shown first, for example embedding_defaults = c("2d" = "X_tsne").
```

## Embedding validation rules

`cellucid-r` enforces:
- `embeddings` is a named list with at least one entry, and every id obeys the
  field-identity rule (non-empty, distinct, drawable exactly as stored),
- each value is a 2D matrix-like object or a numeric vector,
- a value wider than 3 columns is passed through `cellucid_embedding(columns = ...)`,
- `columns` are distinct whole numbers within the array,
- all embeddings have the same number of rows (`n_cells`),
- `n_cells` is between 1 and 2,147,483,647.

It rejects `NA`, `NaN`, infinities, nonnumeric matrices, and embeddings without
a positive finite coordinate range before publication.

:::{note}
A `data.frame` column that carries a class — a `units` column, a `bit64::integer64`
column, or any S3 class of your own — is refused rather than silently unclassed,
because `as.matrix()` would drop the attribute that says what the numbers mean.
Convert deliberately before the call.
:::

## Embedding normalization (important!)

Before writing `embeddings/<index>_<dim>d.bin`, `cellucid-r` normalizes **each
embedding on its own**:

1) Compute per-axis min/max of the plotted columns.
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
- vector fields scale correctly (each is scaled with the factor of the embedding it is bound to),
- and exports more comparable across datasets.

## File encoding details (for debugging)

Each embedding is one file, `embeddings/<index>_<dim>d.bin`, where `<index>` is
its position after sorting by dimension and then by id. It contains:
- float32 values
- little-endian
- row-major order

“Row-major” means:
- cell 1’s coordinates are stored first,
- then cell 2’s, etc.

You can sanity-check a tiny export by reading the default 2D embedding:

```r
identity <- jsonlite::read_json(file.path(out_dir, "dataset_identity.json"))
item <- identity$embeddings$items[[identity$embeddings$defaults[["2d"]]]]
path <- file.path(out_dir, item$file)
con <- if (endsWith(path, ".gz")) gzfile(path, open = "rb") else file(path, open = "rb")
on.exit(close(con), add = TRUE)

vals <- readBin(con, what = "numeric", size = 4, endian = "little", n = n_cells * 2)
coords <- matrix(vals, ncol = 2, byrow = TRUE)
head(coords)
```

## Edge cases

### Constant embeddings (all points identical)

An embedding in which every point sits at the same coordinate is **rejected**,
not exported: `Embedding coordinates must span a positive finite coordinate
range.` It usually means the extraction returned the wrong object, or the
reduction was never computed.

A *nearly* constant embedding is a different case. Any positive range is scaled
up to span ~2 units, so it exports normally, and what you see in the viewer —
one tight clump — is the separation your embedding actually has.

### Several embeddings, several dimensions

The dataset identity file records every embedding under `embeddings.items` —
its `label`, `dimension`, exact `file`, and `source` (`key` and 0-based
`columns`) — plus `embeddings.defaults` and `default_dimension` (the highest
dimension exported). Adding an embedding renumbers the files after it, so always
resolve a file through `dataset_identity.json`.

### “My dataset is mirrored / rotated”

Normalization does not rotate embeddings; it only centers and scales.
If a dataset looks rotated, that’s usually just a different embedding convention.

## Troubleshooting pointers

- `Embedding '…' has 10 columns, and an embedding plots 1, 2, or 3 of them` →
  pass `cellucid_embedding(values, columns = c(1, 2))`.
- `… so embedding_defaults must name the one shown first` → add
  `embedding_defaults = c("2d" = "X_umap")`.
- Viewer loads but points appear “squashed” → likely a coordinate issue (range extremes/outliers).
- Vector fields look wrong → see {doc}`08_vector_fields_velocity_displacement` (they are scaled by the normalization factor of their embedding).
- Exports in format 3 need a viewer release that reads format 3.
