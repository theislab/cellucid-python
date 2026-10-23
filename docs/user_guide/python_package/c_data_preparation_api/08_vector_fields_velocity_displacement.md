# Vector fields (velocity / displacement)

**Audience:** everyone (exact requirements for computational users)
**Time:** 10–20 minutes
**Goal:** store per-cell displacement vectors the web app draws as an animated overlay

Vector fields are optional per-cell **direction + magnitude** vectors in
embedding space: RNA velocity (scVelo), CellRank drift, or any displacement
you want to see as particle flow over an embedding.

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with its `velocity` field; the field, density, speed, trail, size, opacity, palette, and
LOD controls are visible.
```

UI docs: {doc}`../../web_app/i_vector_field_velocity/index`.

## How a vector field is stored

In an AnnData a vector field is stored under scVelo's name for it: the `obsm`
array **`<field>_<basis>`** is the field `<field>` drawn over the embedding
`obsm['X_<basis>']`, for any basis (`velocity_umap` over `X_umap`,
`velocity_umap_3d` over `X_umap_3d`, `velocity_branches` over `X_branches`,
`velocity_pca` over `X_pca`, `T_fwd_umap` over `X_umap`). It has one row per
cell and exactly as many columns as that embedding. A vector drawn over another
layout would point the wrong way, so that is the only binding the file
expresses; the user can bind another array to another embedding in the viewer
(Session > Data mapping > Vector fields).

The key is read by the vector-field key rule
({doc}`09_cellucid_ready_anndata_specification`): a key that starts with `X_`
is an embedding, never a vector field; any other key is split at one `_` into
a field and a basis whose `X_<basis>` is in the same `obsm`, the longest basis
winning (with `X_umap_3d` and `X_3d` both present, `velocity_umap_3d` is
`velocity` on `X_umap_3d`); a key with no such split is an ordinary array.
When the dataset opens, the viewer binds each field to every enabled view of
its embedding when both have the same number of columns, and reads the vector
at the view's columns (`velocity_pca` over `X_pca` shown at columns `[0, 1]`
reads the vector's columns `[0, 1]`). The keys of one field on every embedding
are one field: `velocity_umap`, `velocity_umap_3d` and `velocity_branches` are
all listed as `velocity`.

## `vector_fields=`

`vector_fields` maps each field name to `{embedding key: array}`, each array
with one row per cell and as many columns as its embedding:

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={
        "X_umap": adata.obsm["X_umap"],          # 2D
        "X_umap_3d": adata.obsm["X_umap_3d"],    # 3D
        "X_pca": adata.obsm["X_pca"],            # 50 components
    },
    dataset_id="pancreas",
    dataset_name="Pancreas",
    vector_fields={
        "velocity": {
            "X_umap": adata.obsm["velocity_umap"],
            "X_umap_3d": adata.obsm["velocity_umap_3d"],
            "X_pca": adata.obsm["velocity_pca"],  # 50 columns
        },
        "T_fwd": {"X_umap": adata.obsm["T_fwd_umap"]},
    },
    vector_field_default="velocity",
)
```

This stores `obsm['velocity_umap']`, `obsm['velocity_umap_3d']`,
`obsm['velocity_pca']` and `obsm['T_fwd_umap']`, and writes
`vectorFieldDefault: "velocity"`: the viewer's field ids are `velocity` and
`T_fwd`. `{"velocity": {"X_branches": v}}` would store `velocity_branches`.

Rules (each refused before writing when broken):

- The embedding key is one of `embeddings` and starts with `X_` (`X_<basis>`);
  any other key is refused ("needs an 'X_<basis>' key").
- An array has exactly as many columns as its embedding (any number: 50 for a
  50-column `X_pca`), one row per cell, finite float32-range values. A plain
  `(n_cells,)` vector is one column.
- No two arrays may be stored under the same key, and no stored key may collide
  with an embedding key.
- Every stored key must read back as its own field and embedding by the
  vector-field key rule: a field named `X_...`, or a key another embedding's
  longer basis claims (`T_fwd` on `X_umap` beside an embedding `X_fwd_umap`
  would read as `T` on `X_fwd_umap`), is refused; rename the field.
- `vector_field_default` names one of the `vector_fields` by its field name
  (`"velocity"`), which is also the viewer's field id.

Vectors are stored as given, in the coordinates of the embedding they are drawn
over; the viewer scales them with that embedding when it normalizes it. Do not
rescale them yourself.

## Computing drift from a transition matrix

{func}`~cellucid.add_transition_drift_to_obsm` computes a CellRank-style drift
from a cell-cell transition matrix in the coordinates of the embedding
`obsm['X_<basis>']`, every column of it, and stores it as
`obsm['<field_prefix>_<basis>']` of an AnnData, ready to pass on:

```python
key = cellucid.add_transition_drift_to_obsm(
    adata, transition_matrix, normalize_rows=True
)  # 'T_fwd_umap', as wide as X_umap
key = cellucid.add_transition_drift_to_obsm(
    adata, transition_matrix, basis="umap_3d", normalize_rows=True
)  # 'T_fwd_umap_3d'
```

It refuses a `basis` whose `X_<basis>` the AnnData does not have, and a
`field_prefix` whose key would not read back as that field on that embedding.
{func}`~cellucid.compute_transition_drift` returns the same drift as an array
without storing it.

The R package computes the same drift with `cellucid_transition_drift()` and
stores it through `cellucid_prepare(vector_fields = )`; both packages are held
to one shared conformance fixture
([Computing drift from a transition matrix in R](../../r_package/c_data_preparation_api/08_vector_fields_velocity_displacement.md#computing-drift-from-a-transition-matrix)).

## Next

- {doc}`09_cellucid_ready_anndata_specification`
- {doc}`../f_notebooks_tutorials/33_vector_fields_and_velocity_overlay_end_to_end`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/08_vector_fields_velocity_displacement`
