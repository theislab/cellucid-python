# Embeddings and coordinates

**Audience:** everyone  
**Time:** 15–30 minutes  
**Goal:** export embeddings that render correctly and predictably in the web viewer

Embeddings (UMAP coordinates) are the **minimum required geometry** for Cellucid to render anything.

`prepare()` supports exporting the *same dataset* in multiple dimensionalities (1D, 2D, 3D). The web app can then switch between them at runtime.

---

## Fast path (what you must provide)

You must provide **at least one** of:
- `X_umap_1d` with shape `(n_cells, 1)`, or
- `X_umap_2d` with shape `(n_cells, 2)`, or
- `X_umap_3d` with shape `(n_cells, 3)`.

Minimal sanity checks before export:

```python
import numpy as np

assert X_umap_2d.ndim == 2 and X_umap_2d.shape[1] == 2
assert np.isfinite(X_umap_2d).all()  # no NaN/Inf
```

If you only have `adata.obsm["X_umap"]`, map it to the right argument based on its column count:

```python
X_umap = adata.obsm["X_umap"]
X_umap_2d = X_umap if X_umap.shape[1] == 2 else None
X_umap_3d = X_umap if X_umap.shape[1] == 3 else None
```

---

## Practical path (computational users)

### Supported embedding arguments (exact)

| Argument | Shape | Stored as |
|---|---:|---|
| `X_umap_1d` | `(n_cells, 1)` | `points_1d.bin(.gz)` |
| `X_umap_2d` | `(n_cells, 2)` | `points_2d.bin(.gz)` |
| `X_umap_3d` | `(n_cells, 3)` | `points_3d.bin(.gz)` |
| `X_umap_4d` | `(n_cells, 4)` | **Not supported** (raises `NotImplementedError`) |

All embeddings you provide must:
- be 2D arrays (`ndim == 2`),
- have the correct number of columns for their dimensionality,
- and have the same number of rows (`n_cells`) across all provided dimensions.

### Coordinate normalization (what `prepare()` actually does)

Cellucid normalizes each dimensionality **independently** to fit a stable coordinate range for rendering.

For each dimension `dim ∈ {1,2,3}`:

1) Compute per-axis min/max across all cells.
2) Compute the per-axis ranges and take `max_range = max(axis_ranges)`.
3) Compute the bounding-box center `center = (mins + maxs) / 2`.
4) Scale all axes by the same factor `scale_factor = 2 / max_range` (aspect ratio preserved).
5) Write normalized coordinates:

```text
X_normalized = (X - center) * scale_factor
```

Practical implications:
- Distances are preserved up to a single scale factor **within each dimensionality**.
- You should **not** compare distances between 2D and 3D exports (they are normalized independently).
- If your embedding is degenerate (all points identical), `prepare()` avoids division-by-zero and you’ll get a collapsed point cloud.

### Default dimension and UI switching

The export records:
- which dimensions are available (1/2/3),
- and a default dimension (priority: 3D > 2D > 1D).

In the web app, users can switch between available dimensions.

Related UI docs:
- Dimension switching: {doc}`../../web_app/c_core_interactions/05_dimension_switching_1d_2d_3d`

<!-- SCREENSHOT PLACEHOLDER
ID: embeddings-dimension-switch-badge
Where it appears: User Guide → Python Package → Data Preparation API → Embeddings and coordinates
Capture:
  - Load an exported dataset that has both 2D and 3D points
  - Show the dimension badge/control (e.g., “2D/3D”) in the viewer UI
  - Capture a before/after pair if possible (2D view and 3D view)
Crop:
  - Include: the dimension badge + enough of the plot to show the change
Redact:
  - Remove: sensitive dataset names/IDs if needed
Annotations:
  - Callouts: (1) the dimension badge, (2) what changed (camera/plot appearance)
Alt text:
  - Viewer UI showing the dimension switch control for 2D and 3D embeddings.
Caption:
  - If you export multiple embedding dimensionalities, users can switch between them in the viewer.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the dimension switch control.
:width: 100%

If you export multiple embedding dimensionalities, users can switch between them in the viewer.
```

### “UMAP” naming vs what you can actually use

The argument names are UMAP-branded (`X_umap_*d`), but the viewer treats them as **generic coordinates**.

You can export:
- UMAP coordinates (most common),
- t-SNE coordinates (as long as they are `(n_cells, 2)` or `(n_cells, 3)`),
- spatial coordinates (e.g., x/y) as a 2D embedding.

What you cannot do (yet):
- export multiple *different* embeddings of the same dimensionality and choose among them in the UI (e.g., both UMAP and t-SNE 2D).

If you need multiple embeddings, use separate exports for now.

### Spatial coordinates (common pitfalls)

If you use spatial coordinates as `X_umap_2d`:
- Make sure units and axis scaling are meaningful (Cellucid preserves aspect ratio).
- If your dataset is in microns and spans a huge range, normalization will still fit it to `[-1,1]`, but relative geometry is preserved.
- If you want “square pixels”, ensure x and y are in the same units.

---

## Deep path (expert / developer)

### File format: `points_<dim>d.bin(.gz)`

- dtype: `float32`
- shape: `(n_cells, dim)`
- storage: raw row-major float32 bytes (optionally gzip-compressed)

The web app reads these via typed arrays and uses `dataset_identity.json["embeddings"]` to know which files exist.

Full output spec: {doc}`09_output_format_specification_exports_directory`

### Determinism and reproducibility

Embedding normalization is deterministic given:
- the exact embedding arrays,
- and their floating-point values.

If you recompute UMAP (stochastic), you will get a different export even if everything else is unchanged.

---

## Edge cases and common footguns

- **NaN/Inf coordinates**: must be removed before export.
- **Wrong shape**: `(n_cells,)` is not accepted; use `(n_cells, 1)`.
- **Mismatched rows**: embeddings computed on a subset but `obs`/`X` not subset identically.
- **Collapsed embedding**: all points identical → viewer shows a single point; export is “valid” but useless.
- **4D temptation**: `X_umap_4d` exists for future work, but is explicitly not supported today.

---

## Troubleshooting (embeddings)

### Symptom: `prepare()` errors with “must have exactly 2 columns”

Likely causes:
- You passed a 3D embedding to `X_umap_2d` (or vice versa).

Fix:
- Pass the array to the matching argument name (`X_umap_3d` for `(n_cells, 3)`).

### Symptom: viewer loads but points are missing / everything is blank

Likely causes:
- NaN/Inf in embeddings.
- Exported the wrong folder (not the dataset folder root).

How to confirm:
- `np.isfinite(X_umap).all()`
- Check that `out_dir/dataset_identity.json` exists and points file exists.

Fix:
- Filter/recompute embedding; re-export with `force=True` (or to a new directory).

### Symptom: dimension switch control missing

Meaning:
- Only one of 1D/2D/3D was exported.

Fix:
- Export multiple dimensionalities by providing multiple `X_umap_*d` arrays.

---

## Next steps

- Exporting metadata fields (`obs`): {doc}`04_obs_cell_metadata`
