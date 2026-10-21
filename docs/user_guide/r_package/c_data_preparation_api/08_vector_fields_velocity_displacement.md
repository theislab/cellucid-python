# Vector Fields (Velocity / Displacement)

**Audience:** computational users and developers  
**Time:** 15–25 minutes  
**Goal:** export per-cell vector overlays (velocity/displacement) correctly.

Vector fields are optional, but powerful: they let Cellucid render per-cell arrows/flows (for example RNA velocity).

In `cellucid_prepare()`, vector fields are provided as:

- `vector_fields`: a **named list of named lists** — field id, then the id of
  the exported embedding each array is drawn over:
  `list(velocity = list(X_umap = v2, X_umap_3d = v3))`

Each array must be a numeric matrix with one row per cell and **exactly the
number of columns of the embedding it is bound to**. A plain numeric vector of
length `n_cells` is accepted for a 1D embedding.

## Binding a field to embeddings

A vector field is drawn over one embedding at a time, so every array names the
embedding it belongs to. One field may be bound to several embeddings —
typically a 2D and a 3D UMAP — and an embedding with no binding simply shows no
overlay: a velocity computed in UMAP space points the wrong way over a t-SNE.

- Every embedding id a field names must be one of the exported `embeddings`.
- A field must bind at least one embedding.
- `vector_field_default` names the field selected first; it is required when
  more than one field is exported.

## Field IDs are names, not filenames

A field id is recorded in `dataset_identity.json` and is never used as a path.
Its payloads are written to `vectors/<i>_<e>.bin`, where `<i>` is the field's
position once the ids are sorted by code point and `<e>` is the payload index of
the bound embedding. So a field id carries no filename rule at all — no
character restriction, no length limit, no Windows device-name rule.

What remains are the identity rules every exported identifier obeys: a field id
must be a non-empty string, **distinct**, and text the viewer can draw verbatim —
no control or zero-width characters and no leading or trailing whitespace,
because the id becomes a UI label. `RNA velocity (steady state)` is a legal id.

```{note}
Sorting by code point is what makes `cellucid-r` and the `cellucid` Python
package assign the same index to the same field, so the two writers produce
byte-identical layouts from the same input.
```

### Automatic scaling (critical)

Each embedding is normalized on its own (center + scale). To keep vectors
consistent with the normalized coordinates, `cellucid-r` multiplies each array
by the `scale_factor` of **the embedding it is bound to**.

Practical implication:
- provide vectors in the same coordinate system/units as the original
  coordinates of their embedding
- for an embedding exported from chosen columns (`cellucid_embedding(columns = ...)`),
  pass the vectors for those same columns, in the same order
- the exported vectors will “match” the exported points

## Output files and metadata

Vector files are written under:
- `<out_dir>/vectors/`

Each file is:
- float32, little-endian
- row-major
- shape `(n_cells, dim)` of the bound embedding

Example, for the field at payload index `0` bound to the embeddings at payload
indices `1` and `3`:
- `vectors/0_1.bin` and `vectors/0_3.bin` (or `.gz`)

Metadata is written into `dataset_identity.json` under `vector_fields`, and that
is what maps each field id to its files. Each field entry records exactly
`label` and `files`:

- `label` is derived by the writer, not chosen by you: it is the field id
  **character for character**. The web app compares it and refuses to load the
  dataset if it differs, so there is no display-name field in this format; to
  show a prettier name, rename the list entry.
- `files` maps each bound embedding id to its exact payload path, in embedding
  order.

The exporter also records one exact `default_field`:
- with one vector field, that field is the default;
- with more than one, `vector_field_default` is required, and passing it without
  `vector_fields` is an error rather than a no-op.

## Minimal example

```r
vector_fields <- list(
  velocity = list(
    X_umap = matrix(
      c(0.1, 0,
        0.1, 0,
        0.1, 0),
      ncol = 2, byrow = TRUE
    )
  )
)

cellucid_prepare(
  latent_space = latent,
  obs = obs,
  embeddings = list(X_umap = umap2),
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
- a vector array has a different number of rows than your embeddings (`n_cells`)
- a vector array has a different number of columns than the dimension of the
  embedding it is bound to

### Bindings

Export fails if:
- a field's value is not a named list of arrays,
- a field binds no embedding,
- a field names an embedding id that is not exported.

## Troubleshooting pointers

- `Vector field 'velocity' on embedding 'X_umap' must have exactly 2 columns, got 3.`
  → bind the 3D vectors to the 3D embedding, or pass the matching 2D vectors.
- `Vector field 'velocity' is bound to embedding 'X_tsne', which is not one of the exported embeddings: …`
  → export that embedding, or bind the field to one you export.
- `Vector field 'velocity' must map embedding ids to arrays, for example list(velocity = list(X_umap = vectors)).`
  → wrap the array in a list named by its embedding id.
- `vector_field_default is required when more than one vector field is declared.`
  → name the field you want selected on load.
- `vector_field_default 'drift' does not match an available field: …` → the value
  must be one of the field ids.
- `vector_field_default was provided but no vector fields exist.` → you named a
  default but passed no vectors.
- “Vectors look too long/short” → remember vectors are scaled by their embedding's normalization.
