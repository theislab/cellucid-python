# `obs`: Cell Metadata

**Audience:** everyone (especially if you want coloring, legends, filtering, and category summaries)  
**Time:** 15–25 minutes  
**Goal:** export metadata correctly and avoid category/quantization surprises.

`obs` is the cell-level metadata table:

- one row per cell
- one column per metadata field

In the viewer, `obs` powers:
- categorical coloring (clusters, samples, conditions)
- continuous coloring (QC metrics, scores)
- legends and category lists
- some filtering and summarization features

`obs` is required.

## Which columns get exported?

By default, `cellucid_prepare()` exports **all columns** in `obs`.

To export a subset, use:

```r
cellucid_prepare(..., obs_keys = c("cluster", "sample", "score"))
```

## How `cellucid-r` decides “continuous” vs “categorical”

For each `obs[[key]]`, the exporter uses this rule:

- **categorical** if the column is:
  - `factor`
  - `logical`
  - anything else (including `character`, `Date`, `POSIXct`, lists, etc.)
- **continuous** if the column is:
  - `numeric`

### Recommendation (avoid surprises)

Before export, explicitly coerce:
- categories to `factor(...)`
- continuous values to numeric (`as.numeric(...)`)

This is especially important for columns like:
- `"1"`, `"2"`, `"3"` stored as characters (become categorical)
- `Date` columns (become categorical)

## Output files (what gets written)

All `obs` binaries live under:
- `<out_dir>/obs/`

### Continuous fields

For a continuous field `score`:

- float32 (default): `obs/score.values.f32`
- 8-bit quantized: `obs/score.values.u8`
- 16-bit quantized: `obs/score.values.u16`

Quantized exports also record `min_val` and `max_val` in `obs_manifest.json` so the viewer can recover approximate real values.

### Categorical fields

For a categorical field `cluster`:

- codes: `obs/cluster.codes.u8` or `obs/cluster.codes.u16`
- outlier quantiles: `obs/cluster.outliers.f32` (or `.u8/.u16` if quantized)

The **categories list** (levels) is stored in `obs_manifest.json`.

Missing values are encoded as a reserved integer:

| Codes dtype | Missing marker |
|---|---|
| `uint8` | `255` |
| `uint16` | `65535` |

## Quantization (continuous fields and categorical outliers)

### What quantization does

Quantization maps floating values to integer bins to save space:

- 8-bit: `0..254` for valid values, `255` reserved for missing/invalid
- 16-bit: `0..65534` for valid values, `65535` reserved for missing/invalid

Invalid values are:
- `NA`
- `Inf`
- `-Inf`

### Quantizing continuous fields

Use:
- `obs_continuous_quantization = 8` (recommended default for big exports), or
- `obs_continuous_quantization = 16` (higher precision, larger)

```r
cellucid_prepare(..., obs_continuous_quantization = 8)
```

### Quantizing categorical “outlier quantiles”

Categorical outlier quantiles are also continuous values, so they follow the same quantization setting:

- if `obs_continuous_quantization` is set, outlier files are `.u8/.u16`
- otherwise outlier files are float32 (`.f32`)

## Categorical centroids and outlier quantiles (why `latent_space` is required)

For each categorical field:

### A) Centroids (embedding space)

`cellucid-r` computes per-category centroids in embedding space (for each exported dimension).

Behavior:
- categories with fewer than `centroid_min_points` cells are skipped
- optional “inlier-only” centroids using `centroid_outlier_quantile`:
  - compute distances to the centroid
  - keep points up to the distance quantile
  - recompute centroid using inliers (if enough points)

These centroids are stored in the categorical entry in `obs_manifest.json`.

### B) Outlier quantiles (latent space)

`cellucid-r` computes a per-cell “how typical is this cell inside its category?” score:

1) for each category (with at least `centroid_min_points` cells),
2) compute a latent-space centroid,
3) compute each cell’s distance to that centroid,
4) convert distances to quantile ranks within that category.

This is why `latent_space` is required.

Cells in categories smaller than `centroid_min_points` get `NaN` outlier quantiles.

## Categorical dtype selection (`obs_categorical_dtype`)

Categorical codes can be written as `uint8` or `uint16`.

Default:
- `obs_categorical_dtype = "auto"`

Auto behavior:
- ≤ 254 categories → `uint8`
- > 254 categories → `uint16`

If you force `uint8` and a field has too many categories, export fails with an explicit error.

## Naming and filename safety

Obs column names are sanitized to become filenames:
- unsupported characters become underscores
- leading/trailing dots/underscores are removed

This is convenient, but can cause collisions if two column names sanitize to the same string.

Recommendation:
- keep obs column names simple and stable (letters/numbers/underscores)

## Edge cases (common in real datasets)

### Continuous field is all missing / all infinite

If a continuous field has no valid values:
- export still succeeds,
- `min_val`/`max_val` fall back to `0`/`1`,
- all values become the reserved missing marker.

### Continuous field is constant

If all valid values are the same:
- export still succeeds,
- but the recorded range is widened (`max_val = min_val + 1`) to avoid divide-by-zero.

### Massive categorical fields

Fields with thousands of categories are technically exportable, but often unusable in the UI (legends become too large and humans can’t interpret them).

Recommendations:
- collapse rare categories (“Other”)
- export only the fields you need (`obs_keys`)

## Troubleshooting pointers

- “My numeric column became categorical” → check the R type (`str(obs$col)`).
- “Export fails: uint8 can only hold 254” → use `obs_categorical_dtype="auto"` or `"uint16"`.
- “Outlier files are all missing/NaN” → category counts below `centroid_min_points` or latent space issues.
- Full symptom-based troubleshooting: {doc}`11_troubleshooting_prepare_export`
