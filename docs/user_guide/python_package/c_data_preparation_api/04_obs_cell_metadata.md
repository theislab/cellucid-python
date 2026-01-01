# Obs (cell metadata)

**Audience:** everyone (this is the part that makes the viewer useful)  
**Time:** 30–60 minutes (longer if you have messy metadata)  
**Goal:** export cell metadata that behaves correctly for coloring, filtering, legends, and outlier workflows

In Cellucid terminology:
- **`obs` fields** = per-cell metadata (clusters, sample, batch, QC metrics, scores, pseudotime, etc.).
- These are what most users interact with first (color-by, filters, legends).

`prepare()` exports `obs` into:
- `obs_manifest.json` (a compact index of exported fields),
- and one or more binary files per field under `obs/`.

```{important}
Cell identity is row order. `obs` must have the same row order as embeddings and latent space.
If your `obs` is misaligned, the export will load and look “wrong” in ways that are hard to detect later.
```

---

## Fast path (make `obs` work)

If you’re exporting a dataset for the first time, start small.

### A minimal recommended `obs` set

For most single-cell datasets:
- one or two categorical fields:
  - `cluster` (or `leiden`, `louvain`)
  - `sample` / `batch`
- a few continuous fields:
  - `n_counts`, `n_genes`
  - `pct_mito`
  - a score (e.g., `cell_cycle_score`)

Then add more fields only once the export pipeline is stable.

### Quick rules

- If a column is “numeric but categorical” (e.g., cluster labels `0,1,2`), cast it to categorical before export.
- Avoid columns with **one unique value per cell** (cell IDs) unless your dataset is small (see edge cases).
- Avoid mixing types inside one column (e.g., numbers and strings) → make it consistent first.

---

## Practical path (computational users)

### How `prepare()` classifies obs fields (continuous vs categorical)

`prepare()` decides field kind using pandas dtypes:

- **categorical** if:
  - pandas dtype is categorical (`is_categorical_dtype`), or
  - dtype is bool, or
  - dtype is “other” (strings/objects/datetimes/etc.)
- **continuous** if:
  - pandas dtype is numeric

This is critical, because it changes:
- the UI control (category legend vs continuous colormap),
- how values are stored,
- and whether the outlier filter is available (outlier quantiles are computed for categorical fields).

#### Example: cluster labels that look numeric

If `adata.obs["leiden"]` is stored as strings, it will export as categorical (good).
If it is stored as integers, it will export as continuous (usually not what you want).

```python
adata.obs["leiden"] = adata.obs["leiden"].astype("category")
```

### Selecting which columns to export (`obs_keys`)

By default, `prepare()` exports **all** columns in `obs`.

You can (and often should) export only a curated set:

```python
obs_keys = ["leiden", "cell_type", "sample", "n_counts", "pct_mito"]

prepare(
    ...,
    obs=adata.obs,
    obs_keys=obs_keys,
    ...
)
```

If `obs_keys` contains keys that do not exist, `prepare()` raises a clear error.

Why subsetting matters:
- it reduces clutter for non-technical users,
- it reduces export size (especially if you have huge string columns),
- and it prevents accidental leakage of private identifiers (patient IDs, barcodes).

### Continuous columns (numeric)

Continuous obs fields are written as:
- `float32` by default, or
- quantized `uint8/uint16` if `obs_continuous_quantization` is set.

#### Quantization (optional, lossy)

If you set:
- `obs_continuous_quantization=8`, values are stored as `uint8`
- `obs_continuous_quantization=16`, values are stored as `uint16`

Encoding rules (current exporter):
- valid values map to `0..254` (8-bit) or `0..65534` (16-bit)
- invalid values (`NaN`, `Inf`, `-Inf`) map to the reserved missing marker:
  - `255` (8-bit) or `65535` (16-bit)
- per-field `minValue`/`maxValue` are recorded in the manifest and used to dequantize in the browser.

Dequantization in the web app is:

```text
value = minValue + q * (maxValue - minValue) / maxQuant
```

Where `maxQuant = 254` or `65534`.

Practical guidance:
- Use quantization for large exports to reduce size.
- Prefer 16-bit if you care about subtle gradients (scores/pseudotime).
- If you have extreme outliers, quantization uses min/max and can compress the useful range.
  Consider clipping or transforming values before export.

### Categorical columns

Categorical fields are stored as:
- a list of category labels (strings) in the manifest,
- and a per-cell integer code array written to `obs/<field>.codes.u8` or `.u16`.

Missing values:
- pandas missing categories (`NaN`) become a reserved missing code:
  - `255` for `uint8`
  - `65535` for `uint16`

#### `obs_categorical_dtype` (auto vs forcing)

`obs_categorical_dtype` controls how codes are stored:

- `auto` (recommended):
  - if `n_categories ≤ 254` → `uint8` codes
  - else → `uint16` codes
- `uint8`:
  - forces `uint8` and errors if `n_categories > 254`
- `uint16`:
  - forces `uint16` (useful if you want stable dtype across datasets)

Practical guidance:
- Leave this on `auto` unless you are building a pipeline that needs strict consistency.
- Treat fields with tens of thousands of categories as a design smell; the UI will not be usable.

### Outlier quantiles for categorical fields (why `latent_space` is required)

For every categorical field, `prepare()` computes **per-cell outlier quantiles** using `latent_space`.

High-level idea:
- For each category, compute the centroid in `latent_space`.
- Compute each cell’s distance to its category centroid.
- Rank those distances within the category and convert to a `0..1` quantile-like value.
- Export those per-cell values as `<field>.outliers.*`.

These quantiles are used by the web app’s **outlier filter** (when a categorical field is the active field).

UI context (how users experience this):
- {doc}`../../web_app/e_filtering/02_outlier_filtering_per_active_field`

Important details:
- Categories with fewer than `centroid_min_points` cells get `NaN` outlier quantiles (exported as “missing”).
- Outlier quantiles can also be quantized using `obs_continuous_quantization`.

### Centroids (category label positions)

`prepare()` also computes per-category **centroids** in embedding space, for each exported dimension (1D/2D/3D).

Centroid computation has an optional outlier-removal step:
- compute category center,
- remove points beyond `centroid_outlier_quantile` distance (default 0.95),
- recompute center if enough points remain.

Why you care:
- centroids are used by the viewer for category label placement and some UI affordances.

If you set `centroid_outlier_quantile=None`, centroids are exported as empty lists (labels may be missing).

---

## Naming rules (safe filenames)

Obs keys are converted to safe filenames using:

```text
safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
safe = safe.strip("._")
safe = safe or "field"
```

This affects the filenames under `obs/`, and can create collisions.

Recommendation:
- keep your exported keys already-safe and unique (`cluster`, `cell_type`, `sample_id`, `pct_mito`)
- avoid spaces/slashes/punctuation in keys you intend to export long-term

---

## Deep path (formats and schemas)

### `obs_manifest.json` (compact format)

The Python exporter writes a compact manifest:
- `"_format": "compact_v1"`
- schema patterns under `"_obsSchemas"`
- continuous field tuples under `"_continuousFields"`
- categorical field tuples under `"_categoricalFields"`

The web app expands this into a verbose “fields array” at load time.

### Binary payloads under `obs/`

Per continuous field:
- `obs/<safe_key>.values.f32` (or `.values.u8` / `.values.u16`)

Per categorical field:
- `obs/<safe_key>.codes.u8` (or `.codes.u16`)
- `obs/<safe_key>.outliers.f32` (or `.outliers.u8` / `.outliers.u16`)

All files may have a `.gz` suffix if `compression` is enabled.

Full spec: {doc}`09_output_format_specification_exports_directory`

---

## Edge cases and common footguns

### “Cell IDs” in `obs`

If you export a column where every cell has a unique string ID:
- it becomes categorical,
- category count becomes `n_cells`,
- and you hit practical/technical limits quickly:
  - `uint8` max categories: 254
  - `uint16` max categories: 65534
  - large category lists bloat `obs_manifest.json`

Preferred pattern:
- keep cell IDs in your analysis object (AnnData),
- map selection indices back to `adata.obs_names` in Python when needed.

### Datetime columns

Datetimes become categorical unless you convert them. If you want time as continuous:
- convert to numeric (e.g., days since baseline) before export.

### Mixed dtype columns

Object columns that mix numbers and strings can produce confusing categories.
Make them consistent before export.

### Duplicate column names (or collisions after sanitization)

Two different keys can map to the same safe filename key and overwrite files.
Avoid this.

---

## Troubleshooting (obs)

### Symptom: a field is missing in the UI

Likely causes:
- It wasn’t exported because you passed `obs_keys` and forgot it.
- You reused an `out_dir` and `obs_manifest.json` was skipped (`force=False`).
- The key name collides after sanitization and another field overwrote it.

How to confirm:
- Open `<out_dir>/obs_manifest.json` and search for the key.
- Check that a corresponding file exists under `obs/`.

Fix:
- Re-export with `force=True` (or a new output directory).
- Ensure the key exists in `obs` and is included in `obs_keys`.

### Symptom: a cluster label is shown as a continuous gradient

Meaning:
- The column is numeric, so it was exported as continuous.

Fix:
- Convert it to categorical before export:

```python
adata.obs["cluster"] = adata.obs["cluster"].astype("category")
```

### Symptom: outlier filtering does nothing

Likely causes:
- The active field is continuous (outlier quantiles are computed for categorical fields).
- Categories are too small (`centroid_min_points` too high) so outlier quantiles are missing.
- You’re looking at a dataset that was exported without updating manifests (stale files).

Fix:
- Use a categorical active field.
- Consider lowering `centroid_min_points` if your categories are small.
- Re-export with `force=True`.

### Symptom: export succeeds but the UI becomes unusable (huge legend / slow field list)

Likely causes:
- Extremely high category-count field exported (e.g., cell barcodes, per-cell IDs).

Fix:
- Remove that field from `obs_keys` and re-export.
- Replace with a higher-level grouping field.

---

## Next steps

- Gene metadata (`var`) and gene IDs: {doc}`05_var_gene_metadata`
- Gene expression export (size/performance/quantization): {doc}`06_gene_expression_matrix`
