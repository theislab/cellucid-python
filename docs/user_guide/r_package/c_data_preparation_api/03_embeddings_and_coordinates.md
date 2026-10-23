# Embeddings and coordinates

**Audience:** everyone
**Time:** 10 minutes
**Goal:** pass embeddings correctly and choose what the viewer shows first.

## What to pass

`embeddings` is a **named list**; each name is an `obsm` key and each value has
one row per cell:

```r
embeddings = list(
  X_umap = umap,        # cells x 2
  X_umap_3d = umap3,    # cells x 3
  X_pca = pca           # cells x 50
)
```

Each value may be a numeric matrix, a numeric `data.frame`, a `Matrix` matrix,
or a plain numeric vector (one column). At least one embedding is required, and
all must have the same number of rows; that number is the cell count of the
whole store.

## How they are stored

Each embedding is stored **whole**, as `obsm[<name>]`:

- every column, not only the first two or three;
- as float32, at its own coordinates (nothing is centred or scaled);
- in chunks of whole rows, every column in each chunk
  ({doc}`09_cellucid_ready_anndata_specification`).

The viewer chooses which columns to plot and, unless told otherwise, scales
them to its view. A 50-column PCA is therefore one array, from which the viewer
can show columns 1-2 in 2D or 1-3 in 3D, and which can also serve as the latent
space for outlier quantiles.

Values must be finite and within the float32 range
({doc}`02_input_requirements_global`). Names follow the path rules on the same
page.

## What the viewer shows first

Without launch defaults the viewer picks an embedding for each view dimension
itself, by name: a UMAP before a t-SNE before other layouts, and PCA, scVI,
Harmony or other latent-space names last. To decide yourself, use `embedding_defaults`:

```r
# The 2D view shows X_umap (its first two columns).
embedding_defaults = c("2d" = "X_umap")

# The 3D view shows columns 1, 2 and 3 of X_pca.
embedding_defaults = list(
  "2d" = "X_umap",
  "3d" = list(key = "X_pca", columns = 1:3)
)
```

- The names are `"1d"`, `"2d"`, and `"3d"`, each at most once.
- A value is an `embeddings` name, plotted by its first columns; it needs at
  least as many columns as the dimension.
- Or `list(key = , columns = )`: exactly as many distinct, whole, 1-based
  column numbers as the dimension. They are stored 0-based, as AnnData indexes
  them (`columns = 1:3` is stored as `[0, 1, 2]`).

These are stored in `uns['cellucid']['defaults']` as `embedding1d`,
`embedding2d`, and `embedding3d`.

## Related launch defaults

- `latent_key = "X_pca"`: the latent space the viewer's outlier filter
  measures in (`defaults['latentKey']`), such as PCA or scVI, never a 2D/3D
  layout. It must be one of the `embeddings`. The outlier quantiles of every
  categorical obs column (each cell's distance to its category's centre, as
  the share of its category at that distance or closer) are computed here and
  stored in `uns['cellucid']['outlierQuantiles']`
  ({doc}`09_cellucid_ready_anndata_specification`), so the viewer filters outliers at
  once instead of reading the whole latent space to compute them.
- `normalize_embeddings = FALSE`: show embeddings at their own coordinates
  instead of scaled to the view (`defaults['normalizeEmbeddings']`).
- `centroid_min_points`: the fewest cells a category needs to get outlier
  quantiles (10 when not set).

Every default is optional and stored only when set. In the viewer, all of them
can be changed under **Session › Data mapping**
({doc}`../../web_app/b_data_loading/13_data_mapping`).

## Common mistakes

- **A matrix with a different number of rows** than the others:
  `All embeddings must have one row per cell: 'X_umap' has 3 rows and 'X_b' has 2.`
- **Non-finite coordinates**: `Embedding 'X_umap' must contain only finite values.`
  followed by counts and the first positions.
- **A default the embedding cannot fill**:
  `embedding_defaults[["3d"]] names 'X_umap', which has 2 columns; a 3d view needs 3.`

## Next

- {doc}`04_obs_cell_metadata`
- {doc}`08_vector_fields_velocity_displacement`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/03_embeddings_and_coordinates`
