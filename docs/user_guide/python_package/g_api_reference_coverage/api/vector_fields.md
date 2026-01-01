# Vector fields (velocity / drift overlays)

```{eval-rst}
.. currentmodule:: cellucid
```

Cellucid “vector fields” are **per-cell displacement vectors in embedding space**.

They enable overlays such as:
- RNA velocity / flow arrows (scVelo-style, if you already have vectors)
- CellRank drift vectors derived from a transition matrix

This page documents:
- {func}`~cellucid.compute_transition_drift`
- {func}`~cellucid.add_transition_drift_to_obsm`

---

## Mental model (beginner-friendly)

- You have an embedding like UMAP: `X_umap` with shape `(n_cells, 2)` or `(n_cells, 3)`.
- You compute (or already have) vectors with the *same* shape.
- Cellucid visualizes those vectors as an animated overlay on top of the embedding.

Vectors are not “absolute positions”; they are **arrows attached to each cell**.

---

## Practical path (common workflows)

### 1) From a CellRank transition matrix → drift vectors

```python
from cellucid import compute_transition_drift

drift = compute_transition_drift(T, adata.obsm["X_umap"], normalize_rows=True)
```

Where:
- `T` is `(n_cells, n_cells)` (dense or sparse)
- `adata.obsm["X_umap"]` is `(n_cells, dim)`

### 2) Store drift in `adata.obsm` using Cellucid naming conventions

```python
from cellucid import add_transition_drift_to_obsm

key = add_transition_drift_to_obsm(
    adata,
    T,
    basis="umap",
    field_prefix="T_fwd",
)
print("Wrote vector field to:", key)
```

This writes a key like:
- `T_fwd_umap_2d` or `T_fwd_umap_3d`

### 3) Export the vectors so the viewer can load them

```python
from cellucid import prepare

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    X_umap_2d=adata.obsm.get("X_umap_2d", adata.obsm["X_umap"]),
    vector_fields={
        # You can pass the new obsm entry directly
        key: adata.obsm[key],
    },
    out_dir="./my_export",
)
```

---

## Naming conventions (important)

Cellucid’s recommended keys for vector fields follow:
- Explicit: `<field>_<basis>_<dim>d`
  - examples: `velocity_umap_2d`, `T_fwd_umap_3d`
- Implicit: `<field>_<basis>` (only if you’re sure there’s no ambiguity)

Why explicit keys are recommended:
- it avoids clashes when you have both 2D and 3D embeddings,
- it makes export/serve behavior deterministic.

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
- If your transition matrix is not row-stochastic, set `normalize_rows=True` (default).
- Rows with sum 0 are handled safely (division-by-zero is avoided), but the resulting drift can be zero/undefined depending on the matrix.

### Shape mismatches
- `transition_matrix` must be `(n_cells, n_cells)`
- `embedding` must be `(n_cells, dim)`
- The product `T @ embedding` must produce `(n_cells, dim)`

### Dimension mismatch between export and vectors
- If you export `X_umap_3d` but only provide `*_umap_2d` vectors, the 3D overlay will not be available.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “embedding must be 2D”
Fix:
- Ensure you pass an array shaped `(n_cells, dim)` (not a flattened vector).

### Symptom: “T @ embedding produced shape …”
Fix:
- Confirm `T` is `(n_cells, n_cells)` and aligned to the same cell order as the embedding.

---

### Symptom: “Vector overlay is not available / not visible in the viewer”
Likely causes:
- You did not export the vector field (or you are serving in AnnData mode without exposing it).
- The key naming doesn’t match the embedding dimension currently shown (2D vs 3D).

How to confirm:
- In an exported folder, check that `vectors/` exists and contains files.
- In AnnData mode, confirm the vector field exists in `adata.obsm` with the expected key.

Fix:
- Export vectors via `prepare(..., vector_fields={...})`.
- Prefer explicit keys like `velocity_umap_2d` and `velocity_umap_3d` if you use multiple dimensions.

---

### Symptom: “Drift vectors look backwards”
Likely causes:
- You used a backward transition matrix or the wrong convention for `field_prefix`.

Fix:
- Compute/label forward vs backward consistently (e.g., `field_prefix="T_fwd"` vs `"T_bwd"`).
- Validate by checking a few expected transitions in a small subset.

---

## See also

- {doc}`export` for exporting vector fields via `prepare(...)`
