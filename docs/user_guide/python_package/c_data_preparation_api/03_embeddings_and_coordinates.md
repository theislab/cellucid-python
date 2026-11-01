# Embeddings and coordinates

**Audience:** everyone
**Time:** 10 minutes
**Goal:** pass embeddings, choose what each dimension opens on, and name the latent space

## `embeddings=`

`embeddings` maps each `obsm` key to a per-cell array. Every array is stored
whole as `obsm[key]` in float32: every column, at its own coordinates.

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={
        "X_umap": adata.obsm["X_umap"],          # 2 columns
        "X_umap_3d": adata.obsm["X_umap_3d"],    # 3 columns
        "X_tsne": adata.obsm["X_tsne"],
        "X_pca": adata.obsm["X_pca"],            # 50 columns
    },
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
)
```

- At least one embedding is required; any number may be passed.
- An array is `(n_cells, k)` with `k ≥ 1`; a plain `(n_cells,)` vector is one
  column. Sparse input is densified.
- Values fit in float32. A row holding NaN (or an infinity) is a cell the
  embedding does not place (one computed on some of the cells): it is written
  as given, the viewer leaves it out of that embedding's views and says how
  many cells it leaves out, and every other view and analysis keeps it. At
  least one cell must have finite coordinates. In the latent space
  (`latent_key`), such a cell gets no outlier quantile.
- Keys follow the path rules of {doc}`02_input_requirements_global`.

Nothing is normalized or cropped. The viewer scales each embedding into its
view (unless `normalize_embeddings=False`) and plots any 1 to 3 columns of any
array: a 50-column PCA can be shown as components 1 and 2 in 2D.

## What each dimension opens on

Without `embedding_defaults`, the viewer's own rule applies: `obsm['X_umap_1d']`,
`['X_umap_2d']`, `['X_umap_3d']` fill their dimensions, or a 1- to 3-column
`obsm['X_umap']` fills the dimension its column count states; with neither, it
makes a deterministic best guess from the other arrays and says so.

`embedding_defaults` names the embedding a dimension opens on:

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_tsne": tsne, "X_pca": pca},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    embedding_defaults={"2d": "X_tsne", "3d": ("X_pca", [0, 1, 2])},
)
```

- A key alone (`"X_tsne"`) is plotted by its first columns; it is stored as
  that key.
- `(key, columns)` names the 0-based columns, one per axis; it is stored as
  `{key, columns}`.
- The key must be one of `embeddings`, wide enough for the dimension, and the
  columns distinct and in range; anything else is refused before writing:

  ```text
  embedding_defaults['2d'] names obsm['X_tsne'], which this store does not have. obsm holds: ['X_pca', 'X_umap'].
  ```

The choice is written to `uns['cellucid']['defaults']['embedding2d']` (and
`embedding1d`/`embedding3d`). The user can switch to any other array in the
viewer; this sets only what it opens on.

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: The Navigation block reading Planar above the Compare views block whose Shows picker reads UMAP · 2D: choosing a 2D embedding for a view switched its navigation to Planar.
:width: 480px

Showing a 2D embedding in a view switches that view's navigation to
**Planar**; a mode you choose for the view yourself stays as you set it.
```

Related UI docs:
- Switching embeddings: {doc}`../../web_app/c_core_interactions/05_dimension_switching_1d_2d_3d`
- Choosing embeddings, the graph and the expression matrix: {doc}`../../web_app/b_data_loading/13_data_mapping`

## The latent space: `latent_key=`

The viewer's outlier filter measures, for each categorical field, how far
each cell lies from its category's centre in a latent space (for example scVI
or PCA), never in the 2D/3D layout. Pass the latent space as one of the
`embeddings` and name it:

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_umap": umap, "X_scvi": scvi_latent},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    latent_key="X_scvi",
)
```

`latent_key` must name one of the `embeddings`. It is written as
`uns['cellucid']['defaults']['latentKey']`, and `prepare()` computes the
outlier quantiles of every categorical `obs` column (categorical, text,
boolean, nullable boolean) in that space and stores them:

```text
uns['cellucid']['outlierQuantiles'] = {
    latentKey: "X_scvi",
    minPoints: 10,
    fields: {<obs column>: float32[n_obs], ...},   # NaN where a cell has none
}
```

A cell's quantile is the share of its category at its distance from the
category's centre or closer. A category with fewer than `centroid_min_points`
cells (10 when None) gets none, so its cells are never filtered; a column none
of whose categories reaches it gets no entry. Each field costs 4 bytes per
cell before compression, and the viewer reads one field's array
on the outlier slider's first move, so the filter works at once.

Without `latent_key`, nothing is stored and no latent space is named: the
viewer's outlier filter is shown disabled until one is chosen in **Data
mapping**, and then computes the quantiles in the browser, reading the whole
latent space once ({doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`).

## `normalize_embeddings=`

`normalize_embeddings=False` writes `normalizeEmbeddings: false`, so the viewer
shows every embedding at its own coordinates instead of scaling it into the
view. Leave it unset unless the coordinates themselves are meaningful (a
spatial layout in micrometres, for example).

## Next

- {doc}`04_obs_cell_metadata`
- {doc}`08_vector_fields_velocity_displacement`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/03_embeddings_and_coordinates`
