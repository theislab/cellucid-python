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

Every file is named by the field's **payload index**, not by its key. The index
is the field's position in `obs_keys`, and `obs_manifest.json` is what maps it
back to the column name. Continuous and categorical payloads share one `obs/`
directory, so they share one index space: across `_continuousFields` and
`_categoricalFields` together the indices run `0 … N-1`, each used once.

### Continuous fields

For a continuous field at payload index `0`:

- float32 (default): `obs/0.values.f32`
- 8-bit quantized: `obs/0.values.u8`
- 16-bit quantized: `obs/0.values.u16`

Quantized exports also record `min_val` and `max_val` in `obs_manifest.json` so the viewer can recover approximate real values.

### Categorical fields

For a categorical field at payload index `1`:

- codes: `obs/1.codes.u8` or `obs/1.codes.u16`
- outlier quantiles: `obs/1.outliers.f32` (or `.u8/.u16` if quantized)

The **categories list** (levels) is stored in `obs_manifest.json`.

Missing values are encoded as a reserved integer:

| Codes dtype | Missing marker |
|---|---|
| `uint8` | `255` |
| `uint16` | `65535` |

(r_package-category-label-display-text)=
#### Category labels must read as the value they store

A category label — a factor level, or a distinct value of a character column —
is drawn **verbatim** in the legend, in the field selector, and in every
exported figure. A character with no glyph therefore changes the stored value
without changing the picture, so `cellucid_prepare()` rejects a label that:

- is empty (`""`),
- contains a control character (`U+0001`–`U+001F`, `U+007F`–`U+009F`, which
  includes tab and newline),
- contains a zero-width character (`U+200B` ZERO WIDTH SPACE, `U+2060` WORD
  JOINER, or `U+FEFF`, the byte-order mark a spreadsheet leaves at the front of
  a UTF-8 CSV), or
- starts or ends with whitespace of any kind, including `U+00A0` NO-BREAK SPACE.

It also rejects two labels in one field that a whitespace-collapsing renderer
draws identically — `"T cell"` and `"T  cell"` are two legend rows with one
appearance.

Logical labels are unaffected. Case-only differences are not rejected:
`"Liver"` and `"liver"` look different on screen, so they remain two categories.
An unused factor level is checked like any other, because it is still published
in `obs_manifest.json`.

```{important}
`cellucid_prepare()` never trims a label for you. Trimming rewrites an
annotation you did not ask to change, and where both `"Liver"` and `"Liver "`
exist it would merge two distinct categories into one and move cells between
them.
```

The error names every offending label in the field at once, so one pass fixes
the column:

```r
obs[["organ"]] <- factor(
  trimws(as.character(obs[["organ"]]), whitespace = "[\\h\\v]")
)
```

The same rule applies to `dataset_name`, `dataset_description`, `source_name`,
`source_url`, and `source_citation`, which the viewer also prints verbatim.
`dataset_description` is the one that may be empty: it defaults to `NULL`, and
`NULL`, `""`, and leaving it out all publish the same empty description. The
Python package enforces the identical rule in `prepare()` and accepts `None` and
`""` the same way.

## Quantization (continuous fields and categorical outliers)

### What quantization does

Quantization maps finite floating values to integer bins:

- 8-bit: `0..254`
- 16-bit: `0..65534`

Continuous fields reject `NA`, `NaN`, and infinities. The final integer value
is reserved only for categorical missing codes and generated nullable
categorical outlier quantiles.

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

Choose `obs_categorical_dtype = "uint8"` or `"uint16"` explicitly for every
export. `uint8` stores at most 255 categories; choose `uint16` for larger
categorical fields.

## Naming

Obs column names are recorded exactly in `obs_manifest.json` and are never used
as filenames, so they carry no filename rule: any characters, any length, any
case. `cell type (fine)` and `% mito` both export.

Each exported key must be a non-empty string, **distinct** among the keys you
export so a lookup by key resolves one field, and **text the viewer can draw
verbatim** — no control or zero-width characters and no leading or trailing
whitespace, the same rule category labels obey. A duplicate stops the export
with `Observation field key 'cluster' is duplicated.`, naming the first repeat.

Recommendation:
- keep obs column names stable, because a renamed key is a different field to a
  saved session

## Edge cases (common in real datasets)

### Continuous field contains a non-finite value

Any `NA`, `NaN`, or infinity rejects the complete candidate.

### Continuous field is constant

If all valid values are the same, the field exports normally:
- quantized, it takes the compact contract's constant case —
  `minValue == maxValue` and every code `0` — and the viewer decodes it back to
  the exact constant;
- unquantized, it is a float32 payload like any other.

The field still colours the embedding one flat colour, which is a property of
the data rather than of the export.

### Massive categorical fields

Fields with thousands of categories are technically exportable, but often unusable in the UI (legends become too large and humans can’t interpret them).

Recommendations:
- collapse rare categories (“Other”)
- export only the fields you need (`obs_keys`)

## Troubleshooting pointers

- “My numeric column became categorical” → check the R type (`str(obs$col)`).
- “Export fails because uint8 capacity is exceeded” → use
  `obs_categorical_dtype="uint16"`.
- “Outlier files are all missing/NaN” → category counts below `centroid_min_points` or latent space issues.
- Full symptom-based troubleshooting: {doc}`11_troubleshooting_prepare_export`
