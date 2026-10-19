# Embeddings and coordinates

**Audience:** everyone  
**Time:** 15–30 minutes  
**Goal:** export embeddings that render correctly and predictably in the web viewer

Embeddings are the **minimum required geometry** for Cellucid to render anything.

`prepare()` exports **any number of embeddings, in any of 1D, 2D, and 3D** — a
UMAP and a t-SNE in 2D, a 3D UMAP, two PCA components, spatial coordinates —
each under an id of your choosing. One embedding per dimension is that
dimension's default, and the viewer opens in the highest dimension the export
carries.

---

## Fast path (what you must provide)

Pass at least one embedding in `embeddings=`, keyed by its id:

```python
from cellucid import prepare

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    embeddings={"X_umap": adata.obsm["X_umap"]},
    out_dir="./exports/my_dataset",
    dataset_name="My dataset",
    dataset_id="my-dataset-v1",
    obs_categorical_dtype="uint16",
)
```

A plain array with 1, 2, or 3 columns — or a plain `(n_cells,)` vector, for 1D —
is plotted whole, and its column count is its dimension. `prepare()` takes
arrays, not `obsm` keys: the id is only a name, and your object is never renamed
or mutated.

Minimal sanity checks before export:

```python
import numpy as np

umap = adata.obsm["X_umap"]
assert umap.ndim == 2 and umap.shape[1] in (1, 2, 3)
assert np.isfinite(umap).all()  # no NaN/Inf
```

---

## Practical path (computational users)

### Several embeddings, and which one is shown first

```python
from cellucid import Embedding, prepare

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    embeddings={
        "X_umap": Embedding(adata.obsm["X_umap"], label="UMAP", source_key="X_umap"),
        "X_tsne": Embedding(adata.obsm["X_tsne"], label="t-SNE", source_key="X_tsne"),
        "X_pca": Embedding(adata.obsm["X_pca"], columns=[0, 1], label="PC 1 / PC 2"),
        "X_umap_3d": adata.obsm["X_umap_3d"],
    },
    embedding_defaults={"2d": "X_umap"},
    out_dir="./exports/my_dataset",
    dataset_name="My dataset",
    dataset_id="my-dataset-v1",
    obs_categorical_dtype="uint16",
)
```

- **`embeddings`** maps each id to an array or to an {class}`~cellucid.Embedding`.
  An id obeys the field-identity rule: a non-empty string the viewer can draw
  exactly as stored. Ids are distinct by construction, because they are mapping
  keys.
- **`embedding_defaults`** names the embedding each dimension shows first, as
  `{"1d": ..., "2d": ..., "3d": ...}`. It is required for every dimension that
  has more than one embedding, and optional otherwise: a dimension with exactly
  one embedding uses it. A key for a dimension with no embedding, or a value that
  is not an embedding of that dimension, is rejected.
- The highest dimension present is the one the viewer opens in.

### `Embedding`: columns, label, and provenance

| Field | Meaning |
|---|---|
| `values` | a real `(n_cells, width)` array, or a `(n_cells,)` vector (width 1); sparse input is densified |
| `columns` | the 0-based columns to plot, one per axis, so their count (1, 2, or 3) is the dimension. Required when `width > 3`; defaults to every column otherwise |
| `label` | display text for the embedding, such as `"t-SNE"`; defaults to the id |
| `source_key` | provenance, such as the `obsm` key the array came from; recorded in `dataset_identity.json` and changes nothing else; defaults to `None` |

`columns` may reorder or drop columns of a narrow array too:
`Embedding(adata.obsm["X_umap_3d"], columns=[0, 1])` exports the first two UMAP
axes as a 2D embedding. Columns must be distinct and exist in the array. Only
the selected columns are read, so a non-finite value in an unselected PCA
component does not reject the export.

### Coordinate normalization (what `prepare()` actually does)

Cellucid normalizes **each embedding independently** to fit a stable coordinate
range for rendering:

1) Compute per-axis min/max across all cells of the selected columns.
2) Compute the per-axis ranges and take `max_range = max(axis_ranges)`.
3) Compute the bounding-box center `center = (mins + maxs) / 2`.
4) Scale all axes by the same factor `scale_factor = 2 / max_range` (aspect ratio preserved).
5) Write normalized coordinates:

```text
X_normalized = (X - center) * scale_factor
```

Practical implications:
- Distances are preserved up to a single scale factor **within each embedding**.
- You should **not** compare distances between two embeddings (they are normalized independently).
- An embedding with all points identical is rejected because its largest
  per-axis range is zero, so no finite normalization scale exists.
- A vector field is scaled by the factor of the embedding it is bound to; see
  {doc}`08_vector_fields_velocity_displacement`.

### Default dimension and UI switching

The export records every embedding with its dimension, the default embedding of
each dimension, and the highest dimension as the one shown first. In the web app
every embedding of the export is turned on, so each view can show any of them —
two 2D embeddings side by side, or a 3D one beside a 2D one — and the default of
a dimension is what keys `1`, `2` and `3` show first.

Related UI docs:
- Switching embeddings: {doc}`../../web_app/c_core_interactions/05_dimension_switching_1d_2d_3d`
- Turning embeddings on and off: {doc}`../../web_app/b_data_loading/13_data_mapping`

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: Cellucid Compare Views controls (earlier interface) showing 2D dimension and Planar navigation selected.
:width: 246px

For a 2D embedding, Cellucid selects Planar navigation (earlier interface: the
Dimension dropdown is now the embedding picker); a mode the user chooses for a
view stays explicit.
```

### What you can export

The viewer treats every embedding as **generic coordinates**, so you can export:
- UMAP coordinates (most common),
- t-SNE, force-directed, PHATE, or diffusion-map coordinates,
- two or three principal components of a PCA, through `columns`,
- spatial coordinates (e.g., x/y) as a 2D embedding.

Exports in format 3 carry any number of embeddings per dimension, and need a
viewer release that reads format 3.

### Spatial coordinates (common pitfalls)

If you export spatial coordinates:
- Make sure units and axis scaling are meaningful (Cellucid preserves aspect ratio).
- If your dataset is in microns and spans a huge range, normalization will still fit it to `[-1,1]`, but relative geometry is preserved.
- If you want “square pixels”, ensure x and y are in the same units.

```{note}
**`prepare()` and serving AnnData directly are two different paths.**

`prepare()` reads no `obsm` keys: it exports exactly the arrays you pass in
`embeddings=`. Serving an AnnData object directly
({func}`~cellucid.serve_anndata`, {func}`~cellucid.show_anndata`,
`cellucid serve …`) reads `obsm` keys instead — the slots you name with
`embedding_keys`, then `X_umap_1d`/`X_umap_2d`/`X_umap_3d` or a bare `X_umap`,
then a deterministic best guess. The full serve-path rule is in
{doc}`../d_viewing_apis/08_anndata_mode_show_anndata_and_serve_anndata`.
```

---

## Deep path (expert / developer)

### File format: `embeddings/<index>_<dim>d.bin(.gz)`

- one file per embedding, in item order: dimension ascending, then id in
  code-point order. `{index}` is the position in that order
- dtype: `float32`
- shape: `(n_cells, dim)`
- storage: raw row-major little-endian float32 bytes (optionally gzip-compressed)

`dataset_identity.json["embeddings"]` lists every item — its `label`,
`dimension`, exact `file`, and `source` (`key` and the 0-based `columns`) — with
each dimension's default and the default dimension.

Full output spec: {doc}`09_output_format_specification_exports_directory`

### Determinism and reproducibility

Embedding normalization is deterministic given:
- the exact embedding arrays and selected columns,
- and their floating-point values.

If you recompute UMAP (stochastic), you will get a different export even if everything else is unchanged.

---

## Edge cases and common footguns

- **NaN/Inf coordinates**: must be removed from the selected columns before export.
- **A wide array without `columns`**: a 50-column PCA is rejected rather than
  truncated; name the columns to plot.
- **Two embeddings of one dimension, no default**: rejected until
  `embedding_defaults` names the one shown first.
- **Mismatched rows**: embeddings computed on a subset but `obs`/`X` not subset identically.
- **All-identical embedding**: rejected before export; recompute or supply an
  embedding with a nonzero finite range on at least one axis.
- **Renumbered files**: adding an embedding changes the payload index of every
  item after it; resolve paths through `dataset_identity.json`, never by name.

---

## Troubleshooting (embeddings)

### Symptom: `Embedding 'X_pca' has 50 columns, and an embedding plots 1, 2, or 3 of them`

Fix:
- Pass `Embedding(adata.obsm["X_pca"], columns=[0, 1])`, naming the 0-based
  columns to plot.

### Symptom: `... so embedding_defaults must name the one shown first`

Meaning:
- Two or more embeddings share a dimension.

Fix:
- Add `embedding_defaults={"2d": "X_umap"}` (or the dimension the message names).

### Symptom: `Embedding 'X_pca': Each axis needs a different column.`

Fix:
- List each column once in `columns`.

### Symptom: viewer loads but points are missing / everything is blank

Likely causes:
- NaN/Inf in embeddings.
- Exported the wrong folder (not the dataset folder root).

How to confirm:
- `np.isfinite(values).all()` on the array you passed.
- Check that `out_dir/dataset_identity.json` and `out_dir/embeddings/` exist.

Fix:
- Filter/recompute embedding; re-export with `force=True` (or to a new directory).

### Symptom: the `Shows` picker is greyed out

Meaning:
- Only one embedding was exported (or only one is turned on under
  **Session › Data mapping**); its hover text reads "This dataset shows one
  embedding. Turn more on in Session › Data mapping."

Fix:
- Export more than one embedding, for example a 2D and a 3D UMAP, or a 2D
  t-SNE beside the 2D UMAP, and turn them on in **Session › Data mapping**.

---

## Next steps

- Exporting metadata fields (`obs`): {doc}`04_obs_cell_metadata`
