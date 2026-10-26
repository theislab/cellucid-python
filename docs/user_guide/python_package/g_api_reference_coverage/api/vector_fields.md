# Vector fields (velocity / drift overlays)

```{eval-rst}
.. currentmodule:: cellucid
```

Cellucid “vector fields” are **per-cell displacement vectors in embedding space**.

They enable overlays such as:
- RNA velocity particle flow (scVelo-style, if you already have vectors)
- CellRank drift vectors derived from a transition matrix

This page documents:
- {func}`~cellucid.compute_transition_drift`
- {func}`~cellucid.add_transition_drift_to_obsm`

---

## Mental model (beginner-friendly)

- You have an embedding `X_<basis>` in `obsm`, such as `X_umap` with shape
  `(n_cells, 2)`, or a 50-column `X_pca`.
- You compute (or already have) vectors with the *same* shape.
- Cellucid visualizes those vectors as an animated overlay on top of the embedding.

Vectors are not absolute positions. Each row is a displacement at the
corresponding cell, and the viewer integrates those displacements into animated
particle flow.

---

## Practical path (common workflows)

### 1) From a CellRank transition matrix → drift vectors

```python
from cellucid import compute_transition_drift

drift = compute_transition_drift(T, adata.obsm["X_umap"], normalize_rows=True)
```

Where:
- `T` is `(n_cells, n_cells)` (dense or sparse)
- `adata.obsm["X_umap"]` is `(n_cells, k)`, any number of columns `k`; the
  drift has the same shape

### 2) Store drift in `adata.obsm` using Cellucid naming conventions

```python
from cellucid import add_transition_drift_to_obsm

key = add_transition_drift_to_obsm(
    adata,
    T,
    basis="umap",
    field_prefix="T_fwd",
    normalize_rows=True,
)
print("Wrote vector field to:", key)
```

This computes the drift over `obsm['X_<basis>']`, every column of it, and
writes `obsm['<field_prefix>_<basis>']`:
- `T_fwd_umap` (over `X_umap`, the default `basis="umap"`)
- `T_fwd_umap_3d` with `basis="umap_3d"`, `T_fwd_branches` with
  `basis="branches"`

It refuses a `basis` whose `X_<basis>` is missing, and a `field_prefix` whose
key would not read back as that field on that embedding (for example when
another embedding's longer basis claims it).

### 3) Open it, or store it in a Cellucid-ready AnnData

`show_anndata(adata, dataset_name=..., dataset_id=...)` shows the new field
directly. To ship it in a prepared file:

```python
from cellucid import prepare

path = prepare(
    "my_study.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="my-study-v1",
    dataset_name="My study",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    vector_fields={
        # A field name, then the embedding it is drawn over; stored as obsm["T_fwd_umap"]
        "T_fwd": {"X_umap": adata.obsm[key]},
    },
)
```

---

## Naming conventions (important)

The viewer discovers a vector field by scVelo's naming convention: an `obsm`
key `<field>_<basis>` is the field `<field>` drawn over `obsm['X_<basis>']`, for
any basis — `velocity_umap` over `X_umap`, `T_fwd_umap_3d` over `X_umap_3d`,
`velocity_branches` over `X_branches`. These are the keys
{func}`~cellucid.add_transition_drift_to_obsm` writes, and the keys
{func}`~cellucid.prepare` stores `vector_fields={field: {embedding: array}}`
under. A key that starts with `X_` is an embedding, never a vector field; when
a key splits into several bases the file holds, the longest basis wins. The
field is bound when the vector and its embedding have the same number of
columns, and is listed under its `<field>` part (`velocity` for
`velocity_umap` and `velocity_umap_3d`). The full rule:
{doc}`../../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

---

## API reference

```{eval-rst}
.. autofunction:: compute_transition_drift
```

```{eval-rst}
.. autofunction:: add_transition_drift_to_obsm
```

---

## Edge cases (do not skip)

### Row normalization (transition matrices)
- Choose `normalize_rows=True` to divide each transition row by its sum, or
  `normalize_rows=False` when the matrix already has the exact row semantics you
  intend. The choice is required.
- With `normalize_rows=True`, a row that sums to 0 (a cell with no successor)
  is refused: `Every transition_matrix row must have a nonzero finite row sum
  when normalize_rows=True.` Give the cell a successor (a self-transition), or
  pass an already row-stochastic matrix with `normalize_rows=False`, where an
  empty row is a cell drifting to the origin.

### Shape mismatches
- `transition_matrix` must be `(n_cells, n_cells)`
- `embedding` must be `(n_cells, k)`, any number of columns `k`
- The product `T @ embedding` must produce `(n_cells, k)`

### Missing embeddings for a field
- If a dataset has a 3D UMAP but the field has only an array for the 2D one, the 3D overlay will not be available.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “embedding must be a 2D array, got shape …”
Fix:
- Ensure you pass an array shaped `(n_cells, dim)` (not a flattened vector).

### Symptom: “transition_matrix must have shape (n, n), got …”
Fix:
- Confirm `T` is `(n_cells, n_cells)` and aligned to the same cell order as the embedding.

---

### Symptom: “Vector overlay is not available / not visible in the viewer”
Likely causes:
- The vector field is not in `obsm` under a `<field>_<basis>` key beside its
  `X_<basis>` (a `velocity_umap_2d` beside a bare `X_umap` with no
  `X_umap_2d` no longer binds; scVelo's `velocity_umap` does).
- The vector does not have as many columns as its embedding.
- The view shows another embedding than the one the field is named for.

How to confirm:
- Confirm the vector field exists in `adata.obsm` (or the prepared file's
  `obsm`) with the expected key.

Fix:
- Store the vectors with `add_transition_drift_to_obsm`, or via
  `prepare(..., vector_fields={...})`.
- Name each array after its embedding, such as `velocity_umap` for `X_umap`
  and `velocity_umap_3d` for `X_umap_3d`.

---

### Symptom: “Drift vectors look backwards”
Likely causes:
- You used a backward transition matrix or the wrong convention for `field_prefix`.

Fix:
- Compute/label forward vs backward consistently (e.g., `field_prefix="T_fwd"` vs `"T_bwd"`).
- Validate by checking a few expected transitions in a small subset.

---

## See also

- {doc}`prepare` for storing vector fields via `prepare(...)`
- The R package's `cellucid_transition_drift()`, the same computation, held to
  the same shared conformance fixture:
  {doc}`../../../r_package/c_data_preparation_api/08_vector_fields_velocity_displacement`
