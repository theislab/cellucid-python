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
├── dataset_identity.json                      # always required
├── obs_manifest.json                          # always required
├── obs/                                       # payload per declared obs field
└── points_2d.bin(.gz) or points_3d.bin(.gz)   # at least one points file required
```

`obs_manifest.json` is **not optional and never degrades gracefully**. The web
app loads it as a required file for every dataset and fails the whole dataset
when it is absent, exactly as it does for a missing `dataset_identity.json`.
Neither writer can produce an export without it: `obs` is a required argument of
{func}`~cellucid.prepare` and of `cellucid_prepare()`, and both create the
`obs/` directory on every export. `obs/` is empty only when the observation
frame declares no fields at all.

### Typical full export folder (recommended)

```text
my_export/
├── dataset_identity.json             # required
├── obs_manifest.json                 # required
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

Many paths in manifests use `{key}` template tokens. Every exported observation
key and gene ID must already satisfy the one portable identifier contract;
neither the exporter nor loader rewrites them. The exact identifier appears in
the manifest and path. Unsafe identifiers and case-insensitive collisions reject
the complete candidate (see {doc}`02_input_requirements_global`).

The contract covers exactly the identifiers this directory contains. An `obs`
column left out by `obs_keys`, or a `var` gene left out by `gene_identifiers`,
produces no manifest entry and no payload file, so it is not held to it. Gene
IDs are additionally required to be distinct across all of `var`, which is what
makes `gene_identifiers` resolve one row per identifier.

## Exact key sets (applies to every manifest)

Every JSON object in this format is validated against an **exact key set**, not
a minimum one. This is the single rule most often missed when a manifest is
written or patched by hand:

- a **missing** key is fatal, and
- an **unknown extra** key is equally fatal — there is no room for private
  annotations, editor metadata, or comment keys.

The rule governs `dataset_identity.json` (and its nested `stats`, `embeddings`,
`export_settings`, `source`, and `vector_fields` objects), `obs_manifest.json`,
`var_manifest.json`, `connectivity_manifest.json`, the manifest schema objects
(`_obsSchemas.continuous`, `_obsSchemas.categorical`, `_varSchema`), and each
categorical centroid entry. The key lists in this page are therefore complete
specifications, not examples: emit every key shown and nothing else.

Only `dataset_identity.json` has optional keys at all — `created_at`,
`export_settings`, `source`, and `vector_fields` may be absent, and everything
else in it is required.

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

The reader validates this object against an exact key set, so a key that is not
listed here makes the dataset fail to load. In this example
`cellucid_data_version` is `0.9.1`. <!-- CELLUCID_VERSION -->

```json
{
  "version": 2,
  "id": "pbmc_demo",
  "name": "PBMC demo",
  "description": "",
  "created_at": "2026-01-01T00:00:00Z",
  "cellucid_data_version": "0.9.1",
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
    { "key": "pct_mito", "kind": "continuous" },
    { "key": "leiden", "kind": "category", "n_categories": 12 }
  ],
  "export_settings": {
    "compression": 6,
    "var_quantization": 8,
    "obs_continuous_quantization": 8,
    "obs_categorical_dtype": "uint16"
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
- `name`, `description`, and `source.name` / `source.url` / `source.citation`
  are shown to the reader verbatim, so they obey the same display-text rule as
  a string category label: no control characters, no zero-width characters, and
  no leading or trailing whitespace. `description` may be the empty string;
  the others may not be empty. Both writers reject a candidate export that
  breaks the rule, and `generate_datasets_manifest()` rejects a published
  `dataset_identity.json` whose `name` or `description` breaks it.
- `n_genes` is the **exported** gene count, not the input matrix gene count. If
  you narrow the export with `gene_identifiers`, `n_genes` is the size of that
  subset. Both writers derive it from the exported field list, and the reader
  rejects the dataset unless `n_genes` equals `var_manifest.json` `fields`
  length exactly. Use `0` when no gene expression was exported; neither writer
  emits `var_manifest.json` in that case.
- `stats.n_edges` must equal `connectivity_manifest.json` `n_edges`, and
  `stats.has_connectivity` must agree with whether that manifest exists.
- `obs_fields` must list **every continuous field first, then every categorical
  field**, in the same order as `obs_manifest.json`. Both writers always emit
  that order, and the reader compares the two lists position by position for
  `key`, `kind`, and `n_categories` — any other order fails to load. Within each
  group the order is the order the fields were exported in.

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

Top-level keys — all eight are required, and no other key is allowed:

- `"_format": "compact_v1"`
- `"n_points": <n_cells>`
- `"centroid_outlier_quantile": <number in (0, 1) or null>`
- `"latent_key": "latent_space"`
- `"compression": <gzip level or null>`
- `"_obsSchemas": { ... }`
- `"_continuousFields": [...]`
- `"_categoricalFields": [...]`

`centroid_outlier_quantile` records the quantile used to compute the per-category
outlier payloads. Both writers accept `None`/`NULL` or a finite number greater
than `0.5` and less than `1`, and default to `0.95`; the reader accepts `null`
or any finite number strictly between `0` and `1`.

`latent_key` is always the exact string `"latent_space"` in exports written by
either writer. The reader requires the key to be present and to be `null` or a
non-empty string.

A manifest that omits `centroid_outlier_quantile` or `latent_key` is rejected by
both writers' own validators and fails to load in the web app.

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
- `{key}` is the exact previously validated portable identifier.
- `{ext}` for codes is chosen based on the field’s codes dtype (`uint8 → u8`, `uint16 → u16`).
- `_obsSchemas` itself is an exact key set: it holds `continuous` only when
  `_continuousFields` is non-empty, `categorical` only when
  `_categoricalFields` is non-empty, and nothing else.
- `quantizationBits` appears on the continuous schema only when
  `"quantized": true`.

#### `_continuousFields`

Each entry is either:
- `[key]` (non-quantized float32), or
- `[key, minValue, maxValue]` (quantized)

#### `_categoricalFields`

Each entry is either:
- `[key, categories, codesDtype, codesMissingValue, centroidsByDim]` (float32 outliers), or
- `[key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]` (quantized outliers)

Where:
- `categories` is a list of JSON scalars — strings, booleans, or finite numbers
  — that is unique under its exact JSON representation. Every **string** label
  must additionally be text a reader can see: non-empty, free of control
  characters (`U+0000`–`U+001F`, `U+007F`–`U+009F`), free of the zero-width
  characters `U+200B`, `U+2060`, and `U+FEFF`, and without leading or trailing
  whitespace. No two string labels may collapse to the same text when runs of
  whitespace are reduced to one space. Both writers reject a candidate export
  that breaks either rule; neither writer trims a label. See
  [Category labels must read as the value they store](python_package-category-label-display-text).
- `codesDtype` is `"uint8"` or `"uint16"`
- `codesMissingValue` is `255` or `65535`
- `centroidsByDim` is a dict keyed by dimension strings (`"1"`, `"2"`, `"3"`)
  mapping to centroid lists. Each centroid is an exact three-key object
  `{"category": <declared category>, "position": [<dim finite numbers>],
  "n_points": <0..n_points>}`, with at most one entry per category.

The two arrays are read in order — every `_continuousFields` entry first, then
every `_categoricalFields` entry — and that concatenation is the order
`dataset_identity.json` `obs_fields` must reproduce exactly.

### Binary payload files under `obs/`

Per continuous field:
- `obs/<key>.values.f32` or `.values.u8` or `.values.u16` (optionally `.gz`)

Per categorical field:
- `obs/<key>.codes.u8` or `.codes.u16` (optionally `.gz`)
- `obs/<key>.outliers.f32` or `.outliers.u8` or `.outliers.u16` (optionally `.gz`)

Quantized outlier/categorical missing markers:
- 8-bit: `255`
- 16-bit: `65535`

---

## Var manifest: `var_manifest.json` (and `var/` binaries)

`var_manifest.json` indexes gene expression vectors stored under `var/`.

### Manifest structure (compact_v1)

Top-level keys — all seven are required, and no other key is allowed:
- `"_format": "compact_v1"`
- `"n_points": <n_cells>`
- `"var_gene_id_column": null | "<exact-column>"`
- `"compression": <gzip level or null>`
- `"quantization": 8 | 16 | null`
- `"_varSchema": { ... }`
- `"fields": [...]`

`fields` length must equal `dataset_identity.json` `stats.n_genes`. When
`stats.n_genes` is `0` neither writer emits this file, and the reader never
looks for it — `dataset_identity.json` is the index that decides which optional
manifests are fetched at all.

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

When connectivities are exported, the exporter writes a weighted undirected
edge list:

- `connectivity/edges.src.bin(.gz)`
- `connectivity/edges.dst.bin(.gz)`
- `connectivity/edges.weights.f64.bin(.gz)`

`connectivity_manifest.json` contains all twelve keys below and no others:

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
  "weightsPath": "connectivity/edges.weights.f64.bin.gz",
  "weight_dtype": "float64",
  "weight_bytes": 8,
  "compression": 6
}
```

### Cell index width (`index_dtype` / `index_bytes`)

The width is **not fixed at `uint16`**. It is a pure function of `n_cells`, and
both writers and the reader compute the same one:

| `n_cells` | `index_dtype` | `index_bytes` |
| --- | --- | --- |
| `1` … `65536` | `"uint16"` | `2` |
| `65537` … `4294967296` | `"uint32"` | `4` |

Above `4294967296` cells the export is rejected: the browser contract supports
at most the complete unsigned 32-bit index domain.

The reader pins these two fields to the **smallest exact** width for the
declared `n_cells` and rejects anything else, so a `uint32` file cannot be
described as `uint16` and a small dataset cannot be padded up to `uint32`. Of
the datasets shipped in `cellucid-datasets`, `garcia` (219,731 cells), `he`
(71,650 cells), and `suo` (561,947 cells) all use `uint32`/`4`; `pancreas`
(3,696 cells) uses `uint16`/`2`. Any reader or converter written against this
format must handle both widths.

`weight_dtype` and `weight_bytes` are the opposite case: they are always exactly
`"float64"` and `8`.

### Edge array rules

- the three arrays are aligned and equally long (`n_edges` entries each):
  `src[i]` connects to `dst[i]` with Float64 weight `weights[i]`
- all three arrays are little-endian: unsigned integers of `index_bytes` width
  for the endpoints, IEEE-754 binary64 for the weights
- indices are 0-based and refer to the exported cell row order; every index is
  less than `n_cells`
- the graph is undirected and stored once per edge as its upper triangle:
  every edge satisfies **`src[i] < dst[i]`**
- edges are in **strict lexicographic order** by `(src, dst)` — strictly
  increasing, so duplicate pairs are impossible
- exported weights are **finite and strictly positive**. A zero weight is not a
  zero-strength edge, it is an absent edge: dense input drops zeros, sparse
  input must not store an explicit zero coordinate, and the reader rejects any
  weight that is not `> 0`
- the input graph must be finite, free of negative weights, exactly symmetric,
  zero on the diagonal, and free of duplicate sparse coordinates; both writers
  verify all of that before writing anything
- `max_neighbors` is the maximum node degree of that undirected graph, counting
  each edge at both endpoints. It is `0` exactly when `n_edges` is `0`, and it
  must be less than `n_cells`

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

Browser requirement:
- gzip-compressed exports require native `DecompressionStream('gzip')`;
- current Chrome, Edge, Firefox, and Safari releases provide this API.

---

## Determinism and reproducibility

Scientific binary payloads, including their gzip bytes, are deterministic given:
- the exact input arrays and DataFrames,
- the same exporter version and compression level.

Gzip payloads use a fixed timestamp and no filename header, so repeated exports
do not change their compressed bytes because of the clock or output directory.
By default, `dataset_identity.json["created_at"]` records the real UTC
generation time. For a complete byte-identical export, pass one provenance-fixed
timestamp explicitly:

```python
prepare(
    ...,
    dataset_name="Canonical build",
    dataset_id="canonical-build-v1",
    obs_categorical_dtype="uint16",
    created_at="2026-07-26T12:34:56Z",
)
```

The value must be a valid timestamp in exact `YYYY-MM-DDTHH:MM:SSZ` form:
seconds only, UTC `Z`, no offset or fractional seconds. Cellucid preserves an
accepted explicit value exactly. Invalid values reject before the output parent
directory is created.

Common non-deterministic sources are **upstream**:
- stochastic embeddings (UMAP),
- non-deterministic clustering,
- different preprocessing (normalization/logging).

Recommended reproducibility practice:
- set `dataset_id` explicitly and keep it stable across re-exports of the “same dataset”.
- fix `created_at` in the canonical builder when whole-directory byte identity
  is part of the artifact contract.
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

- **Rejected candidates**: generation publication is atomic; a failed or
  interrupted candidate is not adopted.
- **Existing targets**: `force=False` rejects a non-empty target, while
  `force=True` publishes one complete replacement generation.
- **Mixed compression**: manually mixing compressed and uncompressed files breaks expectations; keep exports consistent.
- **Case sensitivity**: dataset IDs and filenames can behave differently on macOS vs Linux vs Windows; prefer lowercase IDs.
- **Huge category lists**: categorical fields with massive category counts bloat `obs_manifest.json`.

---

## Troubleshooting (format-level)

### Symptom: web app says “not a valid export”

Likely causes:
- missing `dataset_identity.json` or missing `obs_manifest.json` — both are
  required, and neither has a degraded mode,
- loading the wrong folder (one level above/below the actual export root),
- invalid JSON (truncated write),
- an extra or missing key in any manifest object (see
  [Exact key sets](#exact-key-sets-applies-to-every-manifest)),
- a hand-edited manifest whose cross-file counts disagree —
  `stats.n_genes` vs `var_manifest.json` `fields` length, `stats.n_edges` vs
  `connectivity_manifest.json` `n_edges`, or `obs_fields` vs the obs manifest
  field order.

Fix:
- confirm the folder contains both `dataset_identity.json` and
  `obs_manifest.json`,
- re-export to a clean directory with `force=True`.

### Symptom: web app rejects a hand-assembled or converted export

A directory that was not produced by the current atomic exporter fails as a
whole rather than loading partially. There is no mode in which the app shows
points but silently skips a missing or malformed `obs_manifest.json`,
`var_manifest.json`, or `connectivity_manifest.json`. `dataset_identity.json`
decides which optional manifests are fetched, so a manifest it declares
(`stats.n_genes > 0`, `stats.has_connectivity`) but that is missing or malformed
is fatal, and a manifest present on disk that it does not declare is never
loaded and its data is invisible.

Fix:
- re-export with `force=True`,
- confirm every manifest exists exactly when `dataset_identity.json` says it
  should, and that its keys and counts match this page exactly.

---

## Next steps

- Performance tuning: {doc}`10_performance_tuning_guide_prepare_export`
- Deep troubleshooting: {doc}`11_troubleshooting_prepare_export`
