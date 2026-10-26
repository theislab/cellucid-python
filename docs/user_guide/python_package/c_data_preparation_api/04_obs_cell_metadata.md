# `obs` (cell metadata)

**Audience:** everyone
**Time:** 10 minutes
**Goal:** know how each `obs` column is stored and shown

`obs` is a pandas DataFrame with one row per cell, in the cell order of the
embeddings. Its index names the cells ({doc}`02_input_requirements_global`).
Every column is stored; the viewer lists them as fields to colour and filter
by.

## How each column is stored

The rules are anndata's own encodings, applied the same way by
`cellucid_prepare()` in R:

| Column | Stored as |
|---|---|
| text, where a value repeats or is missing | `categorical`, categories in code-point order, missing = code -1 |
| text, every value distinct and present (barcodes, ids) | `string-array` |
| pandas `category` | `categorical`, its own categories and order (`ordered` kept) |
| float | its float dtype, NaN as missing |
| integer | its integer dtype |
| boolean | `bool` |
| pandas nullable `Int64` / `boolean` with missing values | `nullable-integer` / `nullable-boolean` |

Refused, with the conversion that makes them storable:

- infinities in a float column (`replace it with NaN`);
- text mixed with other values in one column (`.astype(str)`);
- dates, times, intervals, complex numbers and other types with no AnnData
  encoding Cellucid writes (`.astype(str)` or `.astype(float)`);
- categories that are not text or numbers.

Category labels are stored exactly as given. `'Liver'` and `'Liver '` are two
categories; clean labels before preparing if they should be one.

## Choosing and ordering fields: `obs_keys=`

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_umap": umap},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    obs=adata.obs,
    obs_keys=["cell_type", "leiden", "total_counts"],
)
```

`obs_keys` is written to `uns['cellucid']['defaults']['obsKeys']`: the viewer
lists those columns, in that order. Every `obs` column is stored either way, so
other tools still read them. Each key must be a column of `obs`, named once.
To leave columns out of the file, drop them from `obs` before the call.

## Centroids and outliers

Categorical centroids and their labels are computed by the browser from the
embeddings; nothing about them is stored.

Latent-space outlier quantiles are stored when you name a latent space with
`latent_key`: for every categorical `obs` column (categorical, text, boolean
or nullable boolean), each cell's distance to its category's centre in that
latent space, as the share of its category at that distance or closer, in
`uns['cellucid']['outlierQuantiles']`. A category with fewer than
`centroid_min_points` cells (10 when None) gets none, and a column none of
whose categories reaches it gets no entry. The viewer's outlier filter then
works at once; without them it computes the same values in the browser on
first use ({doc}`03_embeddings_and_coordinates`,
{doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`).

## Next

- {doc}`05_var_gene_metadata`
- {doc}`../../web_app/d_fields_coloring_legends/index`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/04_obs_cell_metadata`
