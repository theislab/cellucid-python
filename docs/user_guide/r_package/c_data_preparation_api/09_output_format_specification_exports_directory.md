# Export Directory Format (Specification)

**Audience:** advanced users, maintainers, anyone debugging exports  
**Time:** 20–40 minutes  
**Goal:** know what every file means and how to validate an export folder.

This page is a practical specification of what `cellucid-r` writes to `out_dir`.

```{note}
This describes the export format as produced by `cellucid-r` today.
The format is shared across languages (R/Python), but details can evolve.
Treat JSON manifests as the source of truth.
```

## Top-level layout

At the top level you will typically see:

- `dataset_identity.json` (always)
- `points_1d.bin` / `points_2d.bin` / `points_3d.bin` (at least one; optionally `.gz`)
- `obs_manifest.json` (always)
- `obs/` (always)
- `var_manifest.json` + `var/` (only if gene expression is exported)
- `connectivity_manifest.json` + `connectivity/` (only if connectivities are exported)
- `vectors/` (only if vector fields are exported)

## Binary conventions (applies everywhere)

Unless stated otherwise:

- **little-endian**
- float values are **float32**
- matrices are written in **row-major order**
- if `compression` is enabled, binaries are written as **gzip** (`.gz` suffix)

### Row-major clarification

For a matrix shaped `(n_rows, n_cols)`, the file contains:

```
row1_col1, row1_col2, ..., row1_colN,
row2_col1, row2_col2, ..., row2_colN,
...
```

In R, the exporter achieves row-major by writing `as.vector(t(mat))`.

## Embeddings: `points_*d.bin`

Files:
- `points_1d.bin`
- `points_2d.bin`
- `points_3d.bin`

Encoding:
- float32, little-endian
- row-major

Shape:
- `(n_cells, dim)`

Notes:
- values are normalized (center + scale); see {doc}`03_embeddings_and_coordinates`

## Obs manifest: `obs_manifest.json`

This JSON file tells the viewer:
- how many points exist (`n_points`)
- where obs binary files are (path patterns)
- which fields exist and whether they are categorical/continuous
- how categorical codes and outlier files are encoded

Important keys (high level):
- `_format`: format identifier (currently `"compact_v1"`)
- `n_points`
- `centroid_outlier_quantile`
- `latent_key` (currently `"latent_space"`)
- `compression` (null or integer)
- `_obsSchemas`:
  - schema for continuous and categorical files (path patterns, dtypes, quantization flags)
- `_continuousFields`: list of continuous field entries
- `_categoricalFields`: list of categorical field entries (categories, dtype, missing marker, centroids, outlier ranges if quantized)

## Obs binaries: `obs/`

### Continuous

For each continuous field `key`, one file is written:

- `obs/<safe_key>.values.f32` (float32), or
- `obs/<safe_key>.values.u8` / `.u16` (quantized)

If quantized:
- the manifest entry includes `min_val` and `max_val`
- invalid values (`NA`/`Inf`) are set to the reserved marker

### Categorical codes

For each categorical field `key`, codes are written:

- `obs/<safe_key>.codes.u8` or `.u16`

Codes:
- start at 0
- correspond to `levels(factor(obs[[key]]))`
- missing values are encoded as a reserved marker:
  - `255` (uint8) or `65535` (uint16)

### Categorical outlier quantiles

For each categorical field `key`, outlier quantiles are written:

- `obs/<safe_key>.outliers.f32` (float32), or
- `obs/<safe_key>.outliers.u8` / `.u16` (quantized)

Semantics:
- per-cell value describing the distance rank inside its category in latent space
- categories smaller than `centroid_min_points` may have `NaN` (or missing marker when quantized)

## Var manifest: `var_manifest.json` (only if gene expression exported)

This manifest tells the viewer:
- what the gene IDs are,
- how to find each gene’s values file,
- and (if quantized) the min/max for dequantization.

Important keys (high level):
- `_format` (currently `"compact_v1"`)
- `n_points`
- `var_gene_id_column`
- `compression`
- `quantization` (null, 8, or 16)
- `_varSchema` (path pattern, dtype, quantization flags)
- `fields`: list of gene entries

Each `fields` entry is:
- `[gene_id]` for float32 exports, or
- `[gene_id, min_val, max_val]` for quantized exports

## Var binaries: `var/`

For each exported gene ID, one file is written:

- `var/<safe_gene_id>.values.f32`, or
- `var/<safe_gene_id>.values.u8` / `.u16`

Values:
- length `n_cells` (dense vector)
- float32 or quantized integer

```{warning}
Even sparse input matrices become dense-per-gene output vectors. Plan disk usage accordingly.
```

## Connectivity manifest: `connectivity_manifest.json` (optional)

This is a small JSON file describing an edge-pairs representation.

Important keys:
- `format`: `"edge_pairs"`
- `n_cells`
- `n_edges`
- `max_neighbors`
- `index_dtype` / `index_bytes`
- `sourcesPath`
- `destinationsPath`
- `compression`

## Connectivity binaries: `connectivity/edges.*.bin` (optional)

Files:
- `connectivity/edges.src.bin`
- `connectivity/edges.dst.bin`

Semantics:
- `src[k]` and `dst[k]` form one undirected edge
- only unique edges are kept (`src < dst`)
- indices are **zero-based**
- dtype depends on `n_cells` (`uint16/uint32/uint64`)

## Dataset identity: `dataset_identity.json` (always)

This is the top-level “what’s in this dataset?” file.

Important keys (high level):

- `version` (currently 2)
- `id`, `name`, `description`
- `created_at`
- `cellucid_data_version` (R package version)
- `stats` (counts: cells, exported genes, fields, edges; `n_genes` respects `gene_identifiers`)
- `embeddings`:
  - `available_dimensions`
  - `default_dimension`
  - `files` mapping (e.g. `"2d": "points_2d.bin.gz"`)
- `obs_fields`: summary list of exported obs fields
- `export_settings`: compression/quantization settings
- optional `source` (name/url/citation)
- optional `vector_fields` metadata

## Vector fields: `vectors/` (optional)

Vector binaries are float32 matrices stored row-major:

- `vectors/<field>_2d.bin` etc

The web app discovers them via `dataset_identity.json` (`vector_fields` section).

## Validation checklist (quick)

Before sharing an export folder:

1) Confirm the identity file parses:
   - `jsonlite::read_json(file.path(out_dir, "dataset_identity.json"))`
2) Confirm required files exist:
   - `points_*d.bin[.gz]`
   - `obs_manifest.json`
   - `obs/` directory is non-empty
3) If gene expression is included:
   - `var_manifest.json` exists
   - `var/` contains gene files
4) If connectivities are included:
   - `connectivity_manifest.json` exists
   - `connectivity/edges.src.bin` and `.dst.bin` exist

For deeper validation and binary inspection, see:
- {doc}`../d_viewing_loading/03_validate_exports_and_debug_loading`
