# Vector fields (velocity / displacement)

**Audience:** everyone (visual explanation for wet lab/non-technical; exact requirements for computational users)  
**Time:** 30–60 minutes (longer if you need to derive vectors)  
**Goal:** export per-cell displacement vectors that the web app can render as an animated overlay

Vector fields are optional per-cell **direction + magnitude** vectors in embedding space.

In the web app they appear as an animated overlay (often called “velocity overlay”), commonly used for:
- RNA velocity (scVelo-style),
- CellRank drift vectors,
- trajectory/displacement fields,
- any per-cell vector you want to visualize in embedding coordinates.

---

## Fast path (get an overlay to show up)

Minimum requirements:

1) Export the embedding the vectors are drawn over.
2) Pass `vector_fields` as `{field id: {embedding id: array}}`, where each array
   has shape `(n_cells, dim)` of that embedding.

```text
prepare(
    ...,
    embeddings={"X_umap": adata.obsm["X_umap"]},
    vector_fields={
        "velocity": {"X_umap": adata.obsm["velocity_umap"]},
    },
    ...
)
```

Then load the export in the web app and enable the overlay:
- {doc}`../../web_app/i_vector_field_velocity/index`

---

## Practical path (computational users)

### What vector fields enable in the UI

When vector fields are present, the web app can:
- list available vector fields for the embedding a view shows,
- render an animated particle-flow overlay,
- let users tune overlay parameters (density, speed, etc.).

UI docs (highly recommended if you’re exporting vectors for others):
- {doc}`../../web_app/i_vector_field_velocity/index`

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity_umap particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with `velocity_umap`; the field, density, speed, trail, size,
opacity, palette, and LOD controls are visible.
```

### Binding a field to embeddings

A vector field is drawn over one embedding at a time, so each array is bound to
an exported embedding by its id. One field may be bound to several embeddings —
typically a 2D and a 3D UMAP:

```text
prepare(
    ...,
    embeddings={
        "X_umap": adata.obsm["X_umap"],
        "X_umap_3d": adata.obsm["X_umap_3d"],
        "X_tsne": adata.obsm["X_tsne"],
    },
    embedding_defaults={"2d": "X_umap"},
    vector_fields={
        "velocity": {
            "X_umap": adata.obsm["velocity_umap"],
            "X_umap_3d": adata.obsm["velocity_umap_3d"],
        },
        "T_fwd": {"X_umap": adata.obsm["T_fwd_umap"]},
    },
    vector_field_default="velocity",
    ...
)
```

Here `X_tsne` has no vectors bound to it: the velocity overlay is offered for
the two UMAPs and not for the t-SNE, because a vector computed in UMAP space
points the wrong way in any other layout.

- Every embedding id a field names must be one of the exported `embeddings`.
- A field must bind at least one embedding.
- `vector_field_default` names the field selected first. It is required when
  more than one field is exported, and must be one of them.

### Supported shapes and dtypes

Each bound array is a per-cell array:
- shape: exactly `(n_cells, dim)`, where `dim` is the dimension of the embedding
  it is bound to; a plain `(n_cells,)` vector is accepted for a 1D embedding
- dtype: any real numeric dtype; the export scales your values at their own
  precision and rounds the product to `float32` once, at the write
- values: finite, and representable as a nonzero finite `float32`

Vectors must be aligned to the same cell order as embeddings and `obs`. For an
embedding exported from chosen columns of a wider array (`Embedding(...,
columns=[...])`), pass the vectors for those same columns, in the same order.

#### Identity rule: a vector field id is a name, not a path

A vector payload is named by integer positions, never by its field id, so no
filename rule reaches the id. Exactly as for an obs key, an embedding id, and a
gene name, a vector field id must be:

- a non-empty string,
- distinct from every other vector field id, and
- text the viewer can draw exactly as it is stored: no control characters, no
  zero-width characters, and no leading or trailing whitespace.

`cellucid-r` applies the identical rule. Spaces, slashes, and punctuation are
ordinary characters here — `RNA velocity (steady state)` is a valid field id.
The field's label in the viewer is its id.

### Scale and normalization (why overlays sometimes look “too small”)

Remember:
- `prepare()` normalizes each embedding on its own to fit roughly `[-1, 1]`.
- To keep vectors consistent with those normalized coordinates, the exporter
  scales each array by the scale factor of **the embedding it is bound to**.

Practical implications:
- Provide vectors in the same coordinate system as your *original* embedding.
- Do not manually rescale vectors to `[-1,1]` before passing them unless you know exactly what you’re doing.

If vectors look too small/large in the UI, it’s usually because:
- vectors are not in embedding coordinates (wrong basis), or
- vectors have an unexpected magnitude distribution (needs preprocessing), or
- you bound the field to an embedding other than the one the vectors were computed in.

---

## Deep path (formats and metadata)

Vector files are written under:

```text
out_dir/vectors/<field index>_<embedding index>.bin(.gz)
```

- dtype: `float32`, rounded once from the scaled value
- shape: `(n_cells, dim)` of the bound embedding
- `<field index>` is the field's position in the sorted field list, and
  `<embedding index>` the bound embedding's position in the exported embedding
  order; no id ever appears in a path

The presence and locations of vector files are recorded in
`dataset_identity.json["vector_fields"]`, as `files` mapping each bound embedding
id to its path.

Full spec: {doc}`09_output_format_specification_exports_directory`

```{note}
Serving an AnnData object directly discovers vector fields differently: it reads
`obsm` keys of the form `<field>_umap_<n>d` and binds each to the slot of its
dimension when that slot shows a UMAP. See
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.
```

---

## Edge cases and common footguns

- **Wrong basis**: vectors must be in the same coordinate system as the embedding they are bound to.
- **Row-order mismatch**: vectors aligned to a different cell order than points/obs.
- **Mixed availability**: you bound only a 2D embedding but view the dataset in 3D (no overlay is offered there).
- **All-zero vectors**: overlay renders but looks static/empty.
- **NaN/Inf in vectors**: validation fails before export.

---

## Troubleshooting (vector fields)

### Symptom: overlay option missing in the UI

Meaning:
- no `vector_fields` metadata was found in `dataset_identity.json`, or the
  field is not bound to the embedding the view shows.

How to confirm:
- open `<out_dir>/dataset_identity.json` and read `"vector_fields"` → `files`.
- check that `<out_dir>/vectors/` exists.

Fix:
- export with `vector_fields=...` bound to the embedding you view, and re-export with `force=True`.

### Symptom: `Vector field 'velocity' is bound to embedding 'X_tsne', which is not one of the exported embeddings`

Fix:
- export that embedding in `embeddings=`, or bind the field to one you export.

### Symptom: `Vector field 'velocity' on embedding 'X_umap' must have exactly 2 columns, got 3.`

Meaning:
- the array's width is not the bound embedding's dimension.

Fix:
- bind the 3D vectors to the 3D embedding (or pass the matching 2D vectors).

### Symptom: overlay renders but looks wrong (direction/magnitude)

Likely causes:
- vectors are not in embedding space (wrong coordinate system),
- vectors need preprocessing/clipping,
- the field is bound to the wrong embedding.

Fix:
- verify basis and shape, and validate by plotting vectors against embedding in Python before export.

---

## Next steps

- Output spec (including `dataset_identity.json["vector_fields"]`): {doc}`09_output_format_specification_exports_directory`
