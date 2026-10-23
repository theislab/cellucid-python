# Vector fields and velocity overlay end-to-end

This tutorial explains Cellucid vector fields (velocity/drift overlays) end-to-end:

Data → vectors in embedding space → viewer overlay → debugging when it looks wrong

You will learn:
- what a “vector field” means in Cellucid (per-cell displacement vectors in embedding space)
- how to store vector fields in `AnnData.obsm` with the correct naming conventions
- how to compute drift vectors from a CellRank transition matrix
- how to store vector fields in a Cellucid-ready AnnData via `prepare(...)` (optional)
- the most common failure modes (dimension mismatch, scale, cell-order mismatch)

---

## At a glance

**Audience**
- Computational users (velocity/drift users)
- Developers (format/debugging)

**Time**
- If you already have vectors: ~15–30 minutes
- If you compute drift from a transition matrix: ~30–60 minutes

**Prerequisites**
- An embedding in `adata.obsm` named `X_<basis>` (a 2D or 3D UMAP, a PCA, ...)
- One of:
  - a vector field already computed (e.g. from scVelo) aligned to the embedding
  - a transition matrix (e.g. from CellRank) + an embedding (to compute drift)

---

## What a vector field is in Cellucid

Cellucid visualizes a vector field as:

> for each cell i, a displacement vector Δxᵢ in the same coordinate system as the embedding

So for 2D UMAP:
- embedding: `(n_cells, 2)`
- vector field: `(n_cells, 2)`

For 3D UMAP:
- embedding: `(n_cells, 3)`
- vector field: `(n_cells, 3)`

For a 50-component PCA, the vector field has 50 columns too; the viewer reads
the columns the view shows.

```{important}
The vector field must be aligned to:
- the same cells (same order)
- the same embedding basis (UMAP vs PCA)
- the same number of columns as that embedding
```

---

## Step 1 — Confirm your embedding

Common embedding keys:
- `adata.obsm["X_umap"]` (Scanpy's UMAP, usually 2D)
- `adata.obsm["X_umap_3d"]` (explicit 3D)

```python
print(list(adata.obsm.keys()))
X = adata.obsm["X_umap"]
X.shape
```

---

## Step 2 — Put a vector field into `adata.obsm` (naming matters)

Cellucid detects vector fields using scVelo's naming convention:

**Required:** a key `<field>_<basis>` beside the embedding `X_<basis>` it is
drawn over, with as many columns as that embedding:
- `<field>_umap`  over `X_umap` (shape `(n_cells, 2)` for a 2D UMAP)
- `<field>_umap_3d`  over `X_umap_3d` (shape `(n_cells, 3)`)

Examples:
- `velocity_umap` (the key scVelo writes for the UMAP basis)
- `T_fwd_umap` (forward drift)
- `velocity_branches` over `X_branches`

The viewer lists a field by its `<field>` part, so `velocity_umap` and
`velocity_umap_3d` are the one field `velocity`. A key that names no embedding
in `obsm` — such as `velocity_umap_2d` beside a bare `X_umap` with no
`X_umap_2d` — is not a vector field; name it `velocity_umap`. The exact rule is
in {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

If more than one field is *served*, the direct-viewing APIs require
`vector_field_default` to name the initial field exactly. For example, an object
whose served fields are both example keys above needs
`vector_field_default="velocity"`. Omitting it, or naming an unserved
field, raises `ValueError`.

Serving them at all is opt-in — see the next section.

---

## Option A — You already have velocity vectors in UMAP space

If your pipeline already produced vectors aligned to UMAP (common in scVelo workflows), store them:

```python
adata.obsm["velocity_umap"] = V  # V must have as many columns as X_umap
```

Sanity-check:

```python
import numpy as np

V = adata.obsm["velocity_umap"]
print("V shape:", V.shape)
print("finite fraction:", np.isfinite(V).all(axis=1).mean())
print("median norm:", float(np.median(np.linalg.norm(V, axis=1))))
```

---

## Option B — Compute drift vectors from a transition matrix (CellRank-style)

If you have a transition matrix `T` (shape `(n_cells, n_cells)`), you can compute drift as:

> E[next_embedding | i] − current_embedding

Cellucid provides helpers:
- `cellucid.compute_transition_drift`
- `cellucid.add_transition_drift_to_obsm`

### Compute and store drift

```python
from cellucid import add_transition_drift_to_obsm

# transition_matrix: (n_cells, n_cells), dense or sparse
key = add_transition_drift_to_obsm(
    adata,
    transition_matrix,
    basis="umap",
    field_prefix="T_fwd",      # stored as obsm["T_fwd_umap"], as wide as X_umap
    normalize_rows=True,
    overwrite=False,
)

print("Wrote vector field:", key)
```

---

## Step 3 — View the vector field in Cellucid

### Vector fields on the direct-AnnData path

With `show_anndata()`, `serve_anndata()`, and `cellucid serve`, the viewer reads
the AnnData store itself. It finds every `obsm` key `<field>_<basis>` beside
an `X_<basis>` of the same width from the store's metadata and reads a field's
vectors only when the overlay shows it, so a vector field costs nothing at
startup.

`vector_field_default=` (`--vector-field-default` in the terminal) names the
field the overlay opens on. It must be a field the store declares — naming one
it does not is an error before the server starts. The full rule is in
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.

`prepare()` does not scan `obsm` — Option B below passes the arrays through
`vector_fields=`, bound to the embedding they are drawn over, and stores
each one under scVelo's `<field>_<basis>` key for it.

### Option A: view directly from AnnData (notebook)

```python
from cellucid import show_anndata

# Option A stores obsm["velocity_umap"] (field "velocity"); Option B stores
# obsm["T_fwd_umap"] (field "T_fwd"). Name the field you stored.
default_vector_field = "velocity"

viewer = show_anndata(
    adata,
    height=650,
    dataset_name="My study",
    dataset_id="my-study-v1",
    vector_field_default=default_vector_field,
)
viewer
```

Then in the UI:
- open **Visualization → Vector Field Overlay**
- enable **Show overlay**
- choose the field under **Vector field:** (for example, `velocity` or
  `T_fwd`)
- tune **Particle density:**, **Flow speed:**, and **Trail length:** until
  the animated particle flow is clear without obscuring the cells

Web-app docs for the overlay UI:
- {doc}`../../web_app/i_vector_field_velocity/index`

### Option B: write a Cellucid-ready AnnData via `prepare(...)` (recommended for sharing)

If you are preparing a file anyway:

```python
from cellucid import prepare

# A field name, then the embedding each array is drawn over. The array is
# stored as obsm["T_fwd_umap"], which the viewer reads as field "T_fwd".
vector_fields = {
    "T_fwd": {"X_umap": adata.obsm["T_fwd_umap"]},
}

path = prepare(
    "my_dataset_with_vectors.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="my_dataset_with_vectors",
    dataset_name="My dataset with vectors",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    vector_fields=vector_fields,
    latent_key="X_pca",
)
```

Then open:

```python
from cellucid import show_anndata
viewer = show_anndata(path, height=650)
viewer
```

---

## Interface reference

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with its velocity field; the field, density, speed, trail, size, opacity, palette, and LOD
controls are visible. A `velocity_umap` key is the field `velocity`.
```

---

## Edge cases (the “why is it wrong?” list)

### 1) Width mismatch (most common)

If your embedding is 2D but vectors are 3D (or vice versa), the field is not bound to it by default.

Confirm:
```python
X = adata.obsm["X_umap"]
V = adata.obsm["T_fwd_umap"]
print("X:", X.shape, "V:", V.shape)
```

### 2) Cell-order mismatch (silent but deadly)

If vectors were computed on a different cell ordering than the embedding:
- the overlay may look like random noise

Best practice:
- always compute vectors on the same `adata` object (same row order)
- avoid reindexing between computing `X_umap` and computing `V`

### 3) Vector-magnitude mismatch (particle flow looks too fast or nearly static)

Vectors can be too small or too large relative to the embedding.

Sanity-check magnitudes:
```python
import numpy as np
mag = np.linalg.norm(V, axis=1)
float(np.quantile(mag, 0.5)), float(np.quantile(mag, 0.95))
```

Fix:
- rescale vectors before storing them (multiply by a scalar)
- use **Flow speed:** for presentation only; record any data-space rescaling in
  the analysis that produced the vectors

### 4) Non-finite values (NaN/Inf)

If vectors contain NaN/Inf, many renderers will drop them.

Confirm:
```python
import numpy as np
np.isfinite(V).all()
```

Fix:
- replace non-finite values with 0

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “No vector fields are available in the UI”

Likely causes:
- vector field keys don’t match the naming convention (no `X_<basis>` for the
  key's basis, e.g. `velocity_umap_2d` beside a bare `X_umap` with no
  `X_umap_2d`, which no longer binds)
- vectors are stored in the wrong place (not in `obsm`)
- the view shows an embedding other than the one the field is named for

Fix:
- rename to scVelo's key for the embedding, like `T_fwd_umap` for `X_umap`
- ensure `adata.obsm` (or the prepared file's `obsm`) contains the key before calling `show_anndata`

### Symptom: “Overlay is present but looks like random noise”

Likely causes:
- cell-order mismatch
- wrong basis (vectors in PCA, overlay on UMAP)

Fix:
- recompute vectors on the same embedding you’re visualizing
- store separate fields per basis/dimension explicitly

### Symptom: “Overlay is extremely slow”

Likely causes:
- very large dataset + dense overlay settings

Fix:
- reduce overlay density in the UI
- consider downsampling or focusing on a subset for exploration

---

## Next steps

- {doc}`21_prepare_a_cellucid_ready_anndata` (a shareable file)
- Web app overlay docs: {doc}`../../web_app/i_vector_field_velocity/index`
