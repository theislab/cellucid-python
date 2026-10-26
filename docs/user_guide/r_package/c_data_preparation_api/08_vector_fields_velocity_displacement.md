# Vector fields (velocity / displacement)

**Audience:** computational users with RNA velocity, CellRank drift, or other per-cell vectors
**Time:** 15 minutes
**Goal:** store vectors the viewer draws as an animated overlay on an embedding.

Vector fields are per-cell **direction + magnitude** vectors in embedding
space. The viewer draws them as an animated particle flow over the embedding
({doc}`../../web_app/i_vector_field_velocity/index`).

## What to pass

`vector_fields` maps each field to the embeddings its arrays belong to:

```r
cellucid_prepare(
  ...,
  embeddings = list(X_umap = umap, X_umap_3d = umap3, X_pca = pca),
  vector_fields = list(
    velocity = list(X_umap = velocity_2d, X_umap_3d = velocity_3d, X_pca = velocity_pca),
    drift = list(X_umap = drift_2d)
  ),
  vector_field_default = "velocity"
)
```

Each array has **one row per cell** and **exactly as many columns as its
embedding** (any number: 50 for a 50-column `X_pca`), and must be finite and
within the float32 range. Pass
the vectors in the same coordinates as the embedding itself; the viewer scales
them together with the embedding.

## How they are stored: the `<field>_<basis>` rule

Each array is stored in `obsm` under scVelo's name for it: bound to the
embedding `X_<basis>`, the field `<field>` becomes `obsm['<field>_<basis>']`:

| You pass | Stored as |
|---|---|
| `velocity = list(X_umap = v)` | `obsm['velocity_umap']` |
| `velocity = list(X_umap_3d = v)` | `obsm['velocity_umap_3d']` |
| `velocity = list(X_branches = v)` | `obsm['velocity_branches']` |
| `drift = list(X_umap = v)` | `obsm['drift_umap']` |

The viewer binds such an array by its name: an `obsm` key `<field>_<basis>` is
the field `<field>` drawn over `obsm['X_<basis>']` whenever a view shows that
embedding, provided both have the same number of columns. The arrays of one
field on every embedding are one field in the viewer, listed by its name:
`velocity_umap`, `velocity_umap_3d` and `velocity_branches` are all
**`velocity`**. The exact key rule — `X_` keys are embeddings, the longest
matching basis wins — is in
{doc}`../../python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification`.

That rule decides what can be bound:

- The embedding must be one of the `embeddings` and its name must start with
  `X_` (`X_<basis>`); any other name is refused:

  ```text
  Vector field 'velocity' on embedding 'pca': the viewer binds a vector field to the embedding 'X_<basis>' it is named for ('velocity_umap' to 'X_umap'), so the embedding needs an 'X_<basis>' key. Rename 'pca', for example to 'X_pca'.
  ```

- The array has exactly as many columns as its embedding.
- The stored name must not collide with an embedding name, and two arrays may
  not be stored under the same name.
- The stored name must read back as its own field and embedding: a field named
  `X_...`, or a name another embedding's longer basis claims (`T_fwd` on
  `X_umap` beside an embedding `X_fwd_umap` would read as `T` on
  `X_fwd_umap`), is refused; rename the field.

To draw an array over a layout it is not named for, pass it as one more entry
of `embeddings` (it is stored as `obsm[<name>]` like any other array) and bind
it to that layout in the viewer under **Session › Data mapping › Vector
fields** ({doc}`../../web_app/b_data_loading/13_data_mapping`).

## Computing drift from a transition matrix

`cellucid_transition_drift()` computes a CellRank-style drift from a
cell-cell transition matrix in the coordinates of an embedding, every column
of it: the transition-weighted mean position of each cell's successors minus
its own position, `T %*% E - E`, with the product divided by the row sums of
`T` first when `normalize_rows = TRUE`. Pass the result as a vector field on
that embedding:

```r
drift_2d <- cellucid_transition_drift(transition_matrix, umap, normalize_rows = TRUE)
drift_3d <- cellucid_transition_drift(transition_matrix, umap3, normalize_rows = TRUE)
cellucid_prepare(
  ...,
  embeddings = list(X_umap = umap, X_umap_3d = umap3),
  vector_fields = list(T_fwd = list(X_umap = drift_2d, X_umap_3d = drift_3d))
)  # obsm['T_fwd_umap'] and obsm['T_fwd_umap_3d']
```

`transition_matrix` is cells x cells, finite and non-negative, as a `Matrix`
sparse matrix or a numeric matrix, in the cell order of the embedding;
`embedding` takes any form `embeddings` takes, such as
`Seurat::Embeddings(seu, "umap")` or `SingleCellExperiment::reducedDim(sce, "UMAP")`.
`normalize_rows` has no default: `TRUE` divides every row by its sum and
refuses a row that sums to zero, `FALSE` uses a row-stochastic matrix (such as
CellRank's) as given. The result is a double matrix as wide as the embedding,
rounded to float32 once, when `cellucid_prepare()` writes it.

This is the computation of Python's `cellucid.compute_transition_drift()`.
Python's `add_transition_drift_to_obsm()`, which writes the drift into the
`obsm` of an AnnData held in memory, has no R counterpart: R holds no AnnData
in memory, and `vector_fields` already names the key `<field>_<basis>` and
checks it by the same rule (above). Both packages are held to one shared
conformance fixture computed with exact rational arithmetic: bit for bit where
the product and the row sums are exact in double precision, and otherwise
within the forward-error bound of the computation, since a BLAS, a sparse
kernel and `rowSums()` may add the terms in different orders.

## Which field is selected first

`vector_field_default` names the field selected first, by the name you used in
`vector_fields` (`"velocity"`). It is stored as
`uns['cellucid']['defaults']['vectorFieldDefault']` unchanged, since that name
is also the viewer's field id. Without it the viewer chooses.

## When the overlay does not show

- The view shows another embedding: a field bound to `X_umap` is drawn only
  while the 2D view shows `X_umap`. Choose it with
  `embedding_defaults = c("2d" = "X_umap")` or in the viewer.
- The field has no array for the embedding the view shows: a field bound to
  `X_umap` only is not drawn in a view showing `X_umap_3d`.
- An `obsm[['velocity_umap_2d']]` written by hand beside a bare `X_umap` (no
  `X_umap_2d`) no longer binds: name it `velocity_umap`, scVelo's own name, or
  pass it through `vector_fields = list(velocity = list(X_umap = v))`.
- The vectors are in another coordinate system than the embedding (wrong basis),
  or all zero.

## Next

- {doc}`09_cellucid_ready_anndata_specification`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement`
- Transition drift in Python, `cellucid.add_transition_drift_to_obsm()`: [Computing drift from a transition matrix](../../python_package/c_data_preparation_api/08_vector_fields_velocity_displacement.md#computing-drift-from-a-transition-matrix)
