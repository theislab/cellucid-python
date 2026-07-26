# Vector Fields (Velocity / Displacement)

**Audience:** computational users and developers  
**Time:** 15–25 minutes  
**Goal:** export per-cell vector overlays (velocity/displacement) correctly.

Vector fields are optional, but powerful: they let Cellucid render per-cell arrows/flows (for example RNA velocity).

In `cellucid_prepare()`, vector fields are provided as:

- `vector_fields`: a **named list** of arrays

Each entry must be either:
- a numeric vector (interpreted as 1D), or
- a numeric matrix with **1, 2, or 3 columns**

Rows must be cells.

## Naming conventions (important)

`cellucid-r` requires an explicit dimensional suffix.

Use keys like:
- `velocity_umap_2d`
- `velocity_umap_3d`

The exporter groups these by base field id (`velocity_umap`) and writes one file per dimension that matches an available embedding.

## Field IDs must be filesystem-safe

Like every exported identifier, vector-field IDs must already be portable.

Allowed characters are effectively:
- letters, numbers, `.`, `_`, `-`

The ID must start with a letter or digit, must not end with `.`, and must not
be a Windows device name. Nothing is stripped or rewritten.

If an ID is not safe, export fails with an explicit error suggesting a safe alternative.

## Embedding coupling and automatic scaling

Vector fields are dimension-specific:
- 2D vectors require a 2D embedding to be present
- 3D vectors require a 3D embedding to be present

If you provide a 3D vector field but do not export `X_umap_3d`, the complete
candidate is rejected before publication.

### Automatic scaling (critical)

Embeddings are normalized (center + scale). To keep vectors consistent with the normalized coordinates, `cellucid-r` multiplies vectors by the same per-dimension `scale_factor` used for the embedding normalization.

Practical implication:
- provide vectors in the same coordinate system/units as your original embedding coordinates
- the exported vectors will “match” the exported points

## Output files and metadata

Vector files are written under:
- `<out_dir>/vectors/`

Each file is:
- float32, little-endian
- row-major
- shape `(n_cells, dim)`

Example:
- `vectors/velocity_umap_2d.bin` (or `.gz`)

Metadata is written into `dataset_identity.json` under `vector_fields`.

The exporter records one exact `default_field`:
- with one vector field, that field is the default;
- with more than one vector field, pass `vector_field_default` explicitly.

If the field id ends with `_umap`, the identity metadata adds:
- `basis = "umap"`
and the UI label becomes `"<Title> (UMAP)"`.

## Minimal example

```r
vector_fields <- list(
  velocity_umap_2d = matrix(
    c(0.1, 0,
      0.1, 0,
      0.1, 0),
    ncol = 2, byrow = TRUE
  )
)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  X_umap_2d = umap2,
  vector_fields = vector_fields,
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)
```

## Edge cases

### Mismatched shapes

Export fails if:
- a vector field has a different number of rows than your embeddings (`n_cells`)
- a vector field claims to be 2D/3D but the array has the wrong number of columns

### Bad naming

Export fails if:
- the field ID contains spaces or slashes
- the field ID begins/ends with `.` or `_`

### Duplicate definitions

If two input keys resolve to the same `<field>_umap` and dimension, export
rejects the duplicate declaration. No declaration overrides another.

## Troubleshooting pointers

- “Vector field id contains unsupported characters” → rename the list key to a safe identifier.
- “Vector field declared as 2D but has shape mismatch” → check `ncol(...)`.
- “Vectors look too long/short” → remember vectors are scaled by embedding normalization.
