# Embeddings and coordinates

**Audience:** everyone  
**Time:** 15–30 minutes  
**Goal:** export embeddings that render correctly and predictably in the web viewer

Embeddings (UMAP coordinates) are the **minimum required geometry** for Cellucid to render anything.

`prepare()` supports exporting the *same dataset* in multiple dimensionalities (1D, 2D, 3D). The web app can then switch between them at runtime.

---

## Fast path (what you must provide)

You must provide **at least one** of:
- `X_umap_1d` with shape `(n_cells, 1)`, or a plain `(n_cells,)` vector, or
- `X_umap_2d` with shape `(n_cells, 2)`, or
- `X_umap_3d` with shape `(n_cells, 3)`.

Minimal sanity checks before export:

```python
import numpy as np

assert X_umap_2d.ndim == 2 and X_umap_2d.shape[1] == 2
assert np.isfinite(X_umap_2d).all()  # no NaN/Inf
```

`prepare()` takes arrays, not `obsm` keys: pass each set of coordinates to the
argument named for its dimension. Scanpy writes `sc.tl.umap()` output to
`adata.obsm["X_umap"]`, so an export call reads
`X_umap_2d=adata.obsm["X_umap"]` — your object is never renamed or mutated to
make that work.

```{note}
**`prepare()` and serving AnnData directly are two different paths, with two
different rules. Do not carry one over to the other.**

`prepare()` is unchanged. It takes explicit `X_umap_1d=`, `X_umap_2d=`, and
`X_umap_3d=` **array** arguments, and reads no `obsm` keys at all. Which
dimension an array is exported as is decided by the argument you pass it to.

Serving an AnnData object directly ({func}`~cellucid.serve_anndata`,
{func}`~cellucid.show_anndata`, `cellucid serve …`) reads `obsm` keys instead.
There, `X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each name their own dimension
and always decide it, and an object that declares none of them but carries a
bare `X_umap` has it read at the dimension its own column count states — 1, 2,
or 3 columns, resolved by its own width, with no rename and no change to your
object. A bare `X_umap` of any other width is refused with a message naming its
shape, and an explicit dimensional key anywhere in `obsm` means the bare key
never joins the resolved set.

The full serve-path rule is in
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.
```

---

## Practical path (computational users)

### Supported embedding arguments (exact)

| Argument | Shape | Stored as |
|---|---:|---|
| `X_umap_1d` | `(n_cells, 1)` or `(n_cells,)` | `points_1d.bin(.gz)` |
| `X_umap_2d` | `(n_cells, 2)` | `points_2d.bin(.gz)` |
| `X_umap_3d` | `(n_cells, 3)` | `points_3d.bin(.gz)` |

All embeddings you provide must:
- be 2D arrays (`ndim == 2`) — except `X_umap_1d`, which also accepts a plain
  `(n_cells,)` vector and reshapes it to the single column it declares, exactly
  as AnnData does when such an array is put in `obsm`,
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
- An embedding with all points identical is rejected because its largest
  per-axis range is zero, so no finite normalization scale exists.

### Default dimension and UI switching

The export records:
- which dimensions are available (1/2/3),
- and a default dimension (priority: 3D > 2D > 1D).

In the web app, users can switch between available dimensions.

Related UI docs:
- Dimension switching: {doc}`../../web_app/c_core_interactions/05_dimension_switching_1d_2d_3d`

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: Cellucid Compare Views controls showing 2D dimension and Planar navigation selected.
:width: 246px

For a 2D embedding, Cellucid selects 2D with Planar navigation; both settings remain explicit if the user chooses another valid configuration.
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
- **Right array, wrong argument**: the dimension is decided by the argument name,
  not by the array, so a `(n_cells, 3)` array passed to `X_umap_2d` is rejected
  rather than truncated.
- **Mismatched rows**: embeddings computed on a subset but `obs`/`X` not subset identically.
- **All-identical embedding**: rejected before export; recompute or supply an
  embedding with a nonzero finite range on at least one axis.

---

## Troubleshooting (embeddings)

### Symptom: `prepare()` errors with “must have exactly 2 columns”

Likely causes:
- You passed a 3D embedding to `X_umap_2d` (or vice versa).

Fix:
- Pass the array to the matching argument name (`X_umap_3d` for `(n_cells, 3)`).

### Symptom: the object rendered when served directly, but `prepare()` reports no embedding

Meaning:
- The serve path read the object's own bare `X_umap` at the dimension its column
  count states. `prepare()` reads no `obsm` keys, so an array you never passed
  never reaches it.

Fix:
- Pass the array to the argument for its dimension, for example
  `X_umap_2d=adata.obsm["X_umap"]` for a two-column array.

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
