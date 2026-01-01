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

1) Export the matching embedding dimension (2D vectors require `X_umap_2d`, 3D vectors require `X_umap_3d`).
2) Provide a `vector_fields` dictionary where the values have shape `(n_cells, dim)`.
3) Use a safe field id and the recommended naming convention.

```python
prepare(
    ...,
    X_umap_2d=adata.obsm["X_umap_2d"],
    vector_fields={
        "velocity_umap_2d": adata.obsm["velocity_umap_2d"],
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
- list available vector fields (per embedding dimension),
- render an animated particle-flow / arrow-like overlay,
- let users tune overlay parameters (density, speed, etc.).

UI docs (highly recommended if you’re exporting vectors for others):
- {doc}`../../web_app/i_vector_field_velocity/index`

<!-- SCREENSHOT PLACEHOLDER
ID: vector-overlay-ui-controls
Where it appears: User Guide → Python Package → Data Preparation API → Vector fields
Capture:
  - Load a dataset export that includes vector fields
  - Open the vector/velocity overlay UI panel
  - Show the field dropdown and at least one obvious parameter slider
Crop:
  - Include: overlay panel + enough of the canvas to show overlay is active
Redact:
  - Remove: sensitive dataset names/IDs if needed
Annotations:
  - Callouts: (1) overlay toggle, (2) field selector dropdown, (3) one “speed/density” control
Alt text:
  - Vector overlay controls panel with a selectable vector field and tuning sliders.
Caption:
  - When vector fields are exported, the viewer can render an animated overlay and let users select and tune the field.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the vector overlay controls.
:width: 100%

When vector fields are exported, the viewer can render an animated overlay and let users select and tune the field.
```

### Supported shapes and dtypes

Each vector field is a per-cell array:
- shape: `(n_cells, dim)` (or `(n_cells,)` for 1D, which is reshaped to `(n_cells, 1)`)
- dtype: converted to `float32` for export

Vectors must be aligned to the same cell order as embeddings and `obs`.

### Naming conventions (recommended)

The exporter supports two styles:

**Preferred (explicit dimensional keys):**
- keys end with `_<dim>d` (e.g., `_2d`, `_3d`)
- the base id becomes the “field id” in `dataset_identity.json`

Examples:
- `velocity_umap_2d`, `velocity_umap_3d` → field id `velocity_umap`
- `T_fwd_umap_2d` → field id `T_fwd_umap`

**Fallback (implicit key, dimension inferred from shape):**
- key has no `_<dim>d` suffix
- exporter infers dimension from the array shape and uses the key as the field id

Recommendation:
- Use the explicit form; it’s clearer and avoids accidental mismatches.

#### Safety rule: vector field ids must already be filesystem-safe

Unlike obs keys and gene IDs, vector field ids are strict:
- only `A–Z`, `a–z`, `0–9`, `.`, `_`, `-` are allowed
- spaces or slashes will raise an error

So prefer ids like:
- `velocity_umap`
- `T_fwd_umap`
- `drift_umap`

### Dimension matching and skipping behavior

Vectors are only exported for dimensions you also exported points for.

Example:
- you provide `velocity_umap_3d`, but did not export `X_umap_3d`
- result: the exporter prints a warning and skips the 3D vector file

To avoid confusion:
- if you want overlay in 3D, export both `X_umap_3d` and the corresponding 3D vectors.

### Scale and normalization (why overlays sometimes look “too small”)

Remember:
- `prepare()` normalizes each embedding dimension to fit roughly `[-1, 1]`.
- To keep vectors consistent with those normalized coordinates, the exporter scales vectors by the same per-dimension scale factor.

Practical implications:
- Provide vectors in the same coordinate system as your *original* embedding.
- Do not manually rescale vectors to `[-1,1]` before passing them unless you know exactly what you’re doing.

If vectors look too small/large in the UI, it’s usually because:
- vectors are not in embedding coordinates (wrong basis), or
- vectors have an unexpected magnitude distribution (needs preprocessing), or
- you are viewing a different embedding dimension than the vectors you exported.

---

## Deep path (formats and metadata)

Vector files are written under:

```text
out_dir/vectors/<field_id>_<dim>d.bin(.gz)
```

- dtype: `float32`
- shape: `(n_cells, dim)`

The presence and locations of vector files are recorded in:
- `dataset_identity.json["vector_fields"]`

Full spec: {doc}`09_output_format_specification_exports_directory`

---

## Edge cases and common footguns

- **Wrong basis**: vectors must be in the same coordinate system as the embedding (UMAP space if you name them `_umap_*`).
- **Row-order mismatch**: vectors aligned to a different cell order than points/obs.
- **Mixed availability**: you export 2D vectors but load the dataset in 3D mode (overlay may appear missing).
- **All-zero vectors**: overlay renders but looks static/empty.
- **NaN/Inf in vectors**: exporter will write them; UI behavior is undefined (avoid).

---

## Troubleshooting (vector fields)

### Symptom: overlay option missing in the UI

Meaning:
- no `vector_fields` metadata was found in `dataset_identity.json`.

How to confirm:
- open `<out_dir>/dataset_identity.json` and search for `"vector_fields"`.
- check that `<out_dir>/vectors/` exists.

Fix:
- export with `vector_fields=...` and re-export with `force=True`.

### Symptom: exporter prints “Skipping vector field … embedding points_<dim>d not provided”

Meaning:
- you provided vectors for a dimension you didn’t export points for.

Fix:
- export the matching `X_umap_<dim>d` embedding, or remove that vector dimension.

### Symptom: overlay renders but looks wrong (direction/magnitude)

Likely causes:
- vectors are not in embedding space (wrong coordinate system),
- vectors need preprocessing/clipping,
- you’re looking at a different dimension than the exported vectors.

Fix:
- verify basis and shape, and validate by plotting vectors against embedding in Python before export.

---

## Next steps

- Output spec (including `dataset_identity.json["vector_fields"]`): {doc}`09_output_format_specification_exports_directory`
