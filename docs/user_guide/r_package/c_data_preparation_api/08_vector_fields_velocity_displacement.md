# Vector Fields (Velocity / Displacement)

**Audience:** computational users and developers  
**Time:** 15–25 minutes  
**Goal:** export per-cell vector overlays (velocity/displacement) correctly.

Vector fields are optional, but powerful: they let Cellucid render per-cell arrows/flows (for example RNA velocity).

In `cellucid_prepare()`, vector fields are provided as:

- `vector_fields`: a **named list** of arrays

Each entry must be a numeric matrix with **exactly the number of columns its key
declares** and one row per cell. A `_1d` key also accepts a plain numeric vector
of length `n_cells` as the single column it declares; `_2d` and `_3d` keys
require the full matrix.

Rows must be cells.

## Naming conventions (important)

`cellucid-r` requires an explicit dimensional suffix. Every key must exactly
match `<field>_umap_<1|2|3>d`.

Use keys like:
- `velocity_umap_2d`
- `velocity_umap_3d`

The exporter groups these by base field id (`velocity_umap`) and writes one file per dimension that matches an available embedding.

An unsuffixed key such as `velocity_umap` is rejected with
`Vector field key 'velocity_umap' must exactly match '<field>_umap_<1|2|3>d'.`
The dimension is always declared, never inferred from an array's width: a field
can carry both a 2D and a 3D version, and a suffix that disagrees with the
array's shape is a mistake worth catching rather than a shape to trust. This is
the same rule the Python package enforces, so the same field is named the same
way from either writer.

## Field IDs are names, not filenames

A field id is recorded in `dataset_identity.json` and is never used as a path.
Its payload is written to `vectors/<index>_<dim>d.bin`, where `<index>` is the
field's position once the ids are sorted by code point. So a field id carries no
filename rule at all — no character restriction, no length limit, no Windows
device-name rule.

What remains is the key grammar above (`<field>_umap_<1|2|3>d`) and the identity
rules every exported identifier obeys: a field id must be a non-empty string,
**distinct**, and text the viewer can draw verbatim — no control or zero-width
characters and no leading or trailing whitespace, because the id becomes a UI
label. Two list names that produce the same field id and dimension are the same
list name, so they are rejected with `vector_fields names must be unique.`
before any array is read.

```{note}
Sorting by code point is what makes `cellucid-r` and the `cellucid` Python
package assign the same index to the same field, so the two writers produce
byte-identical layouts from the same input.
```

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

Example, for the field at payload index `0`:
- `vectors/0_2d.bin` (or `.gz`)

Metadata is written into `dataset_identity.json` under `vector_fields`, and that
is what maps each field id to its files.

The exporter records one exact `default_field`:
- with one vector field, that field is the default;
- with more than one vector field, pass `vector_field_default` explicitly.

Every field id ends with `_umap` by construction, so the identity metadata
always records:
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
  obs_categorical_dtype = "uint8",
  dataset_name = "My dataset",
  dataset_id = "my_dataset",
  out_dir = out_dir,
  centroid_min_points = 1,
  force = TRUE
)
```

## Edge cases

### Mismatched shapes

Export fails if:
- a vector field has a different number of rows than your embeddings (`n_cells`)
- a vector field's array has a different number of columns than its key declares,
  which includes passing a plain vector under a `_2d` or `_3d` key

### Bad naming

Export fails if:
- the key is not exactly `<field>_umap_<1|2|3>d` (no suffix, a suffix other than
  `1`/`2`/`3`, or no `_umap` before it)

The part before `_umap` is otherwise unconstrained: interior spaces, slashes, a
leading `.`, and non-ASCII characters are all accepted, because the id names a
field in the identity metadata and never a file. Only invisible characters are
rejected, and only because the id is drawn as a UI label.

### Duplicate definitions

A key determines its field id and dimension exactly, so two declarations of the
same field and dimension are two identical list names, and `cellucid_prepare()`
rejects them with `vector_fields names must be unique.` before reading any
array. No declaration overrides another.

`Velocity_umap_2d` and `velocity_umap_3d` are two different fields and both
export: they no longer share a payload path, so case alone is not a collision.

## Troubleshooting pointers

- `Vector field key '…' must exactly match '<field>_umap_<1|2|3>d'.` → add or
  correct the dimensional suffix on the list name.
- `vector_fields names must be unique.` → the same list name is declared twice;
  remove or rename one of them.
- `Vector field '…' declares 2D but contains 1 components.` → check `ncol(...)`;
  a plain vector counts as one component.
- `Vector field '…' requires a matching 3D embedding.` → export `X_umap_3d`, or
  drop that vector dimension.
- “Vectors look too long/short” → remember vectors are scaled by embedding normalization.
