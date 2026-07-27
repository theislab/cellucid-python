# Vector fields and velocity overlay end-to-end

This tutorial explains Cellucid vector fields (velocity/drift overlays) end-to-end:

Data → vectors in embedding space → viewer overlay → debugging when it looks wrong

You will learn:
- what a “vector field” means in Cellucid (per-cell displacement vectors in embedding space)
- how to store vector fields in `AnnData.obsm` with the correct naming conventions
- how to compute drift vectors from a CellRank transition matrix
- how to export vector fields via `prepare(...)` (optional)
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
- An embedding (UMAP 2D or 3D) in `adata.obsm`
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

```{important}
The vector field must be aligned to:
- the same cells (same order)
- the same embedding basis (UMAP vs PCA)
- the same dimension (2D vs 3D)
```

---

## Step 1 — Confirm your embedding

Common embedding keys:
- `adata.obsm["X_umap_2d"]` (explicit 2D)
- `adata.obsm["X_umap_3d"]` (explicit 3D)

```python
print(list(adata.obsm.keys()))
X = adata.obsm["X_umap_2d"]
X.shape
```

---

## Step 2 — Put a vector field into `adata.obsm` (naming matters)

Cellucid detects vector fields using a naming convention:

**Required:**
- `<field>_umap_2d`  (shape `(n_cells, 2)`)
- `<field>_umap_3d`  (shape `(n_cells, 3)`)

Unsuffixed keys such as `<field>_umap` are not vector-field declarations.

Examples:
- `velocity_umap_2d`
- `T_fwd_umap_2d` (forward drift)

If the AnnData object contains more than one field ID, the direct-viewing APIs
require `vector_field_default` to name the initial field exactly. For example,
an object containing both example keys above needs
`vector_field_default="velocity_umap"`. Omitting it, or naming an unavailable
field, raises `ValueError`.

---

## Option A — You already have velocity vectors in UMAP space

If your pipeline already produced vectors aligned to UMAP (common in scVelo workflows), store them:

```python
adata.obsm["velocity_umap_2d"] = V  # V must be (n_cells, 2)
```

Sanity-check:

```python
import numpy as np

V = adata.obsm["velocity_umap_2d"]
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
    field_prefix="T_fwd",
    dim=2,                     # set to 3 for 3D, or leave None to infer
    normalize_rows=True,
    overwrite=False,
)

print("Wrote vector field:", key)
```

---

## Step 3 — View the vector field in Cellucid

### Option A: view directly from AnnData (notebook)

```python
from cellucid import show_anndata

# This tutorial's Option A stores "velocity_umap"; Option B stores "T_fwd_umap".
# Name the exact field you stored.
default_vector_field = "velocity_umap"

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
- choose the field under **Vector field:** (for example, `velocity_umap` or
  `T_fwd_umap`)
- tune **Particle density:**, **Flow speed:**, and **Trail length:** until
  the animated particle flow is clear without obscuring the cells

Web-app docs for the overlay UI:
- {doc}`../../web_app/i_vector_field_velocity/index`

### Option B: export the vector field via `prepare(...)` (recommended for sharing)

If you are exporting anyway:

```python
from cellucid import prepare

vector_fields = {
    "T_fwd_umap_2d": adata.obsm["T_fwd_umap_2d"],
}

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    X_umap_2d=adata.obsm["X_umap_2d"],
    vector_fields=vector_fields,
    out_dir="exports/my_dataset_with_vectors",
    dataset_name="My dataset with vectors",
    dataset_id="my_dataset_with_vectors",
    obs_categorical_dtype="uint16",
    compression=6,
    var_quantization=8,
)
```

Then open:

```python
from cellucid import show
viewer = show("exports/my_dataset_with_vectors", height=650)
viewer
```

---

## Interface reference

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity_umap particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with `velocity_umap`; the field, density, speed, trail, size,
opacity, palette, and LOD controls are visible.
```

---

## Edge cases (the “why is it wrong?” list)

### 1) Dimension mismatch (most common)

If your embedding is 2D but vectors are 3D (or vice versa), the overlay will fail or look nonsensical.

Confirm:
```python
X = adata.obsm["X_umap_2d"]
V = adata.obsm["T_fwd_umap_2d"]
print("X:", X.shape, "V:", V.shape)
```

### 2) Cell-order mismatch (silent but deadly)

If vectors were computed on a different cell ordering than the embedding:
- the overlay may look like random noise

Best practice:
- always compute vectors on the same `adata` object (same row order)
- avoid reindexing between computing `X_umap_2d` and computing `V`

### 3) Vector-magnitude mismatch (particle flow looks too fast or nearly static)

Vectors can be too small or too large relative to the embedding.

Sanity-check magnitudes:
```python
import numpy as np
mag = np.linalg.norm(V, axis=1)
float(np.quantile(mag, 0.5)), float(np.quantile(mag, 0.95))
```

Fix:
- rescale vectors before storing/exporting (multiply by a scalar)
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
- vector field keys don’t match the naming convention
- vectors are stored in the wrong place (not in `obsm` for AnnData mode)
- export didn’t include vector fields

Fix:
- rename to an explicit key like `T_fwd_umap_2d`
- verify `dataset_identity.json` advertises `vector_fields` (exported mode)
- in AnnData mode, ensure `adata.obsm` contains the key before calling `show_anndata`

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

- {doc}`21_prepare_exports_with_quantization_and_compression` (shareable exports)
- Web app overlay docs: {doc}`../../web_app/i_vector_field_velocity/index`
