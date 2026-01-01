# Output format specification (exports directory)

**Audience:** everyone (beginners get “what files exist”, developers get exact schemas)  
**Time:** 45–90 minutes  
**Goal:** understand exactly what `prepare()` writes and how the web app interprets it

This page is the **ground truth** spec for the on-disk export format produced by {func}`~cellucid.prepare`.

If you just need a high-level overview, start here:
- {doc}`../../web_app/b_data_loading/07_folder_file_format_expectations_high_level_link_to_spec`

---

## Fast path (directory tour)

### Minimal export folder (required)

```text
my_export/
├── dataset_identity.json
└── points_2d.bin(.gz) or points_3d.bin(.gz)   # at least one points file required
```

### Typical full export folder (recommended)

```text
my_export/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # optional (gene expression)
├── connectivity_manifest.json        # optional (KNN edges)
├── points_1d.bin(.gz)                # optional
├── points_2d.bin(.gz)                # optional
├── points_3d.bin(.gz)                # optional
├── obs/                              # obs field binaries
├── var/                              # gene expression binaries
├── connectivity/                     # graph edge binaries
└── vectors/                          # optional vector fields
```

---

## Key naming and safety rules (applies everywhere)

Many paths in manifests use `{key}` placeholders.

The web app expands `{key}` by applying the same sanitization rule as the Python exporter:

```text
safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
safe = safe.strip("._")
safe = safe or "field"
```

This means:
- manifest keys can be human-readable (`"pct mito"`),
- but files on disk use a safe form (`pct_mito`).

Collisions are possible; avoid them (see {doc}`02_input_requirements_global`).

---

## Required file: `dataset_identity.json`

This file is required for exported datasets.

It provides:
- dataset name/id/description metadata,
- counts (cells/genes/fields),
- embeddings metadata (which points files exist and which dimension is default),
- optional vector field metadata,
- export settings used (compression/quantization),
- and lightweight obs field summaries.

### Schema (current exporter)

```json
{
  "version": 2,
  "id": "pbmc_demo",
  "name": "PBMC demo",
  "description": "",
  "created_at": "2026-01-01T00:00:00Z",
  "cellucid_data_version": "0.0.1a2",
  "stats": {
    "n_cells": 10000,
    "n_genes": 2000,
    "n_obs_fields": 12,
    "n_categorical_fields": 4,
    "n_continuous_fields": 8,
    "has_connectivity": true,
    "n_edges": 123456
  },
  "embeddings": {
    "available_dimensions": [2, 3],
    "default_dimension": 3,
    "files": {
      "2d": "points_2d.bin.gz",
      "3d": "points_3d.bin.gz"
    }
  },
  "obs_fields": [
    { "key": "leiden", "kind": "category", "n_categories": 12 },
    { "key": "pct_mito", "kind": "continuous" }
  ],
  "export_settings": {
    "compression": 6,
    "var_quantization": 8,
    "obs_continuous_quantization": 8,
    "obs_categorical_dtype": "auto"
  },
  "source": {
    "name": "Optional source name",
    "url": "https://...",
    "citation": "..."
  },
  "vector_fields": {
    "default_field": "velocity_umap",
    "fields": {
      "velocity_umap": {
        "label": "Velocity (UMAP)",
        "available_dimensions": [2, 3],
        "default_dimension": 3,
        "files": {
          "2d": "vectors/velocity_umap_2d.bin.gz",
          "3d": "vectors/velocity_umap_3d.bin.gz"
        },
        "basis": "umap"
      }
    }
  }
}
```

Notes:
- `vector_fields` is present only if you exported vectors.
- `source` is present only if you provided `source_*` metadata.
- `n_genes` currently reflects the input matrix gene count (not necessarily the exported gene subset if you used `gene_identifiers`).

---

## Points files: `points_<dim>d.bin(.gz)`

For each exported dimensionality:

- filename: `points_1d.bin`, `points_2d.bin`, `points_3d.bin` (optionally with `.gz`)
- dtype: `float32`
- shape: `(n_cells, dim)`
- coordinate system: normalized independently per dimension to fit approximately `[-1, 1]`

---

## Obs manifest: `obs_manifest.json` (and `obs/` binaries)

`obs_manifest.json` is a compact manifest that describes:
- which obs fields exist,
- how to find their payload files under `obs/`,
- and (for categorical fields) categories, centroids, and outlier quantiles.

### Manifest structure (compact_v1)

Top-level keys:

- `"_format": "compact_v1"`
- `"n_points": <n_cells>`
- `"compression": <gzip level or null>`
- `"_obsSchemas": { ... }`
- `"_continuousFields": [...]`
- `"_categoricalFields": [...]`

#### `_obsSchemas`

Example (quantized continuous + quantized outliers):

```json
{
  "_obsSchemas": {
    "continuous": {
      "pathPattern": "obs/{key}.values.u8.gz",
      "ext": "u8",
      "dtype": "uint8",
      "quantized": true,
      "quantizationBits": 8
    },
    "categorical": {
      "codesPathPattern": "obs/{key}.codes.{ext}.gz",
      "outlierPathPattern": "obs/{key}.outliers.u8.gz",
      "outlierExt": "u8",
      "outlierDtype": "uint8",
      "outlierQuantized": true
    }
  }
}
```

Notes:
- `{key}` is sanitized by the loader before substitution.
- `{ext}` for codes is chosen based on the field’s codes dtype (`uint8 → u8`, `uint16 → u16`).

#### `_continuousFields`

Each entry is either:
- `[key]` (non-quantized float32), or
- `[key, minValue, maxValue]` (quantized)

#### `_categoricalFields`

Each entry is either:
- `[key, categories, codesDtype, codesMissingValue, centroidsByDim]` (float32 outliers), or
- `[key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]` (quantized outliers)

Where:
- `categories` is a list of strings
- `codesDtype` is `"uint8"` or `"uint16"`
- `codesMissingValue` is `255` or `65535`
- `centroidsByDim` is a dict keyed by dimension strings (`"1"`, `"2"`, `"3"`) mapping to centroid lists

### Binary payload files under `obs/`

Per continuous field:
- `obs/<safe_key>.values.f32` or `.values.u8` or `.values.u16` (optionally `.gz`)

Per categorical field:
- `obs/<safe_key>.codes.u8` or `.codes.u16` (optionally `.gz`)
- `obs/<safe_key>.outliers.f32` or `.outliers.u8` or `.outliers.u16` (optionally `.gz`)

Quantized missing markers:
- 8-bit: `255`
- 16-bit: `65535`

---

## Var manifest: `var_manifest.json` (and `var/` binaries)

`var_manifest.json` indexes gene expression vectors stored under `var/`.

### Manifest structure (compact_v1)

Top-level keys:
- `"_format": "compact_v1"`
- `"n_points": <n_cells>`
- `"var_gene_id_column": "index" | "<column>"`
- `"compression": <gzip level or null>`
- `"quantization": 8 | 16 | null`
- `"_varSchema": { ... }`
- `"fields": [...]`

#### `_varSchema`

Non-quantized:

```json
{
  "_varSchema": {
    "kind": "continuous",
    "pathPattern": "var/{key}.values.f32.gz",
    "ext": "f32",
    "dtype": "float32",
    "quantized": false
  }
}
```

Quantized:

```json
{
  "_varSchema": {
    "kind": "continuous",
    "pathPattern": "var/{key}.values.u8.gz",
    "ext": "u8",
    "dtype": "uint8",
    "quantized": true,
    "quantizationBits": 8
  }
}
```

#### `fields`

Each entry is either:
- `[gene_id]` (non-quantized), or
- `[gene_id, minValue, maxValue]` (quantized)

### Binary payload files under `var/`

For each exported gene:
- `var/<safe_gene_id>.values.f32` or `.values.u8` or `.values.u16` (optionally `.gz`)

Each file stores a vector of length `n_cells`.

---

## Connectivity manifest: `connectivity_manifest.json` (and `connectivity/` binaries)

When connectivities are exported, the exporter writes an unweighted undirected edge list:

- `connectivity/edges.src.bin(.gz)`
- `connectivity/edges.dst.bin(.gz)`

`connectivity_manifest.json` contains:

```json
{
  "format": "edge_pairs",
  "n_cells": 10000,
  "n_edges": 123456,
  "max_neighbors": 30,
  "index_bytes": 2,
  "index_dtype": "uint16",
  "sourcesPath": "connectivity/edges.src.bin.gz",
  "destinationsPath": "connectivity/edges.dst.bin.gz",
  "compression": 6
}
```

Notes:
- edge arrays are parallel: `src[i]` connects to `dst[i]`
- indices are 0-based and refer to the exported cell row order

---

## Vector files: `vectors/<field>_<dim>d.bin(.gz)` (optional)

Vector fields are stored as:
- dtype: `float32`
- shape: `(n_cells, dim)`
- files: `vectors/<field_id>_<dim>d.bin(.gz)`

Vector presence and file locations are indexed in `dataset_identity.json["vector_fields"]`.

---

## Compression strategy (gzip)

If `compression` is set to an integer 1–9:
- all binary files are written with a `.gz` suffix,
- the web app decompresses them in the browser.

Practical guidance:
- `compression=6` is a good balance for final exports.
- Disable compression during iteration to reduce CPU time, then enable for publish/share.

Compatibility note:
- modern browsers support native `DecompressionStream` for gzip;
- for older environments, ensure the web app build includes a gzip fallback (or export uncompressed).

---

## Determinism and reproducibility

Exports are deterministic given:
- the exact input arrays and DataFrames,
- and the same version of the exporter.

Common non-deterministic sources are **upstream**:
- stochastic embeddings (UMAP),
- non-deterministic clustering,
- different preprocessing (normalization/logging).

Recommended reproducibility practice:
- set `dataset_id` explicitly and keep it stable across re-exports of the “same dataset”.
- record your pipeline version and the cellucid package version (`cellucid_data_version` is written automatically).

---

## Multi-dataset export roots (`datasets.json`)

For GitHub-hosted multi-dataset roots, the web app expects:

```text
exports_root/
├── datasets.json
├── dataset_a/
│   └── dataset_identity.json
└── dataset_b/
    └── dataset_identity.json
```

The Python package includes a helper to generate `datasets.json`:

```python
from cellucid.prepare_data import generate_datasets_manifest

generate_datasets_manifest("./exports_root", default_dataset="dataset_a")
```

Web app docs:
- {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`

---

## Edge cases and common footguns

- **Partial exports**: interrupted runs can leave missing files; always verify `dataset_identity.json` + at least one `points_*d.bin`.
- **Stale manifests**: reusing an `out_dir` with `force=False` can skip writing manifests while `dataset_identity.json` updates.
- **Mixed compression**: manually mixing compressed and uncompressed files breaks expectations; keep exports consistent.
- **Case sensitivity**: dataset IDs and filenames can behave differently on macOS vs Linux vs Windows; prefer lowercase IDs.
- **Huge category lists**: categorical fields with massive category counts bloat `obs_manifest.json`.

---

## Troubleshooting (format-level)

### Symptom: web app says “not a valid export”

Likely causes:
- missing `dataset_identity.json`,
- loading the wrong folder (one level above/below the actual export root),
- invalid JSON (truncated write).

Fix:
- confirm the folder contains `dataset_identity.json`,
- re-export to a clean directory with `force=True`.

### Symptom: web app loads points but fields/genes are missing

Likely causes:
- `obs_manifest.json` or `var_manifest.json` missing,
- or stale manifests because export was skipped.

Fix:
- re-export with `force=True`,
- confirm manifest files exist and match the exported content.

---

## Next steps

- Performance tuning: {doc}`10_performance_tuning_guide_prepare_export`
- Deep troubleshooting: {doc}`11_troubleshooting_prepare_export`
