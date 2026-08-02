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
  - the path patterns are templates over `{index}` (and `{ext}` for codes), for
    example `obs/{index}.values.u8.gz` and `obs/{index}.codes.{ext}.gz`
- `_continuousFields`: list of continuous field entries
- `_categoricalFields`: list of categorical field entries (categories, dtype, missing marker, centroids, outlier ranges if quantized)

Every field entry begins with its own **payload index**, an integer:

- continuous: `[index, key]`, or `[index, key, min_val, max_val]` when quantized
- categorical: `[index, key, categories, codes_dtype, codes_missing_value,
  centroids_by_dim]`, or the same six followed by
  `outlier_min_value, outlier_max_value` when outliers are quantized

Both arrays write into `obs/`, so their indices are **one shared space**: across
`_continuousFields` and `_categoricalFields` together the indices are exactly
`0 … N-1`, each used once. `cellucid_prepare()` proves that before publishing,
because two fields sharing an index would overwrite one payload with another and
the viewer would then draw one field's values under another field's name.

```{note}
A continuous entry and a gene entry are both four members long when quantized.
Never infer a field's kind from an entry's length — the array it came from is
what says which it is.
```

Every **string** category label published here must be text a reader can
see: non-empty, free of control characters (`U+0001`-`U+001F`,
`U+007F`-`U+009F`), free of the zero-width characters `U+200B`, `U+2060`,
and `U+FEFF`, and without leading or trailing whitespace. No two string
labels in one field may collapse to the same text when runs of whitespace
are reduced to one space. `cellucid_prepare()` rejects the complete
candidate export rather than trimming a label; the Python writer enforces
the identical rule. See
[Category labels must read as the value they store](r_package-category-label-display-text).

## Obs binaries: `obs/`

Payload filenames are the field's **payload index**, never its key. The manifest
entry beside the file is what says which field an index belongs to, so no
exported path depends on a dataset's vocabulary.

### Continuous

For each continuous field with payload index `i`, one file is written:

- `obs/<i>.values.f32` (float32), or
- `obs/<i>.values.u8` / `.u16` (quantized)

If quantized:
- the manifest entry includes `min_val` and `max_val`
- all input values are finite; non-finite values reject publication

### Categorical codes

For each categorical field with payload index `i`, codes are written:

- `obs/<i>.codes.u8` or `.u16`

Codes:
- start at 0
- correspond to `levels(factor(obs[[key]]))`
- missing values are encoded as a reserved marker:
  - `255` (uint8) or `65535` (uint16)

### Categorical outlier quantiles

For each categorical field with payload index `i`, outlier quantiles are written:

- `obs/<i>.outliers.f32` (float32), or
- `obs/<i>.outliers.u8` / `.u16` (quantized)

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
- `_varSchema` (path pattern, dtype, quantization flags); the pattern is a
  template over `{index}`, for example `var/{index}.values.f32.gz`
- `fields`: list of gene entries

Each `fields` entry is:
- `[index, gene_id]` for float32 exports, or
- `[index, gene_id, min_val, max_val]` for quantized exports

The first element is the gene's payload index. Within `var/` the indices are exactly
`0 … N-1`, each used once, and `cellucid_prepare()` proves that before
publishing. There is **one** identity per gene: the entry carries the exact
identifier read from `rownames(var)` or `var_gene_id_column`, and nothing else.
No accession is retained beside it, and no parallel identifier array exists.

## Var binaries: `var/`

For each exported gene with payload index `i`, one file is written:

- `var/<i>.values.f32`, or
- `var/<i>.values.u8` / `.u16`

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
- `weightsPath`
- `weight_dtype` (`"float64"`) / `weight_bytes` (`8`)
- `compression`

## Connectivity binaries: `connectivity/edges.*.bin` (optional)

Files:
- `connectivity/edges.src.bin`
- `connectivity/edges.dst.bin`
- `connectivity/edges.weights.f64.bin`

Semantics:
- `src[k]`, `dst[k]`, and `weight[k]` form one weighted undirected edge
- only unique edges are kept (`src < dst`)
- indices are **zero-based**
- weights are exact little-endian Float64 values
- dtype is `uint16` through 65,536 cells or `uint32` through
  4,294,967,296 cells; larger axes are rejected
- zero-length arrays represent a present, validated zero-edge graph and are
  distinct from absent connectivity

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

`name`, `description`, and every `source` string are shown to the reader
verbatim, so they obey the same display-text rule as a string category
label. `description` may be `""`; the others may not be empty.
- optional `vector_fields` metadata

## Vector fields: `vectors/` (optional)

Vector binaries are float32 matrices stored row-major:

- `vectors/<i>_2d.bin`, `vectors/<i>_3d.bin`, etc

`<i>` is the field's payload index, assigned by sorting the field ids by code
point so that `cellucid-r` and the `cellucid` Python package assign the same
index to the same field. Within `vectors/` the indices are exactly `0 … N-1`,
each used once.

The web app discovers them via `dataset_identity.json` (`vector_fields` section),
which is what maps each field id to its files.

## The writer proves the layout before publishing

Before a candidate generation is published, `cellucid_prepare()` reconciles what
is on disk against what the export declares, on **five** surfaces: the four
payload directories — `Observation`, `Gene`, `Connectivity`, and `Vector field` —
and the `Export` root itself. On each one, every path the export declares must
exist and every file present must be declared. A mismatch stops the export with,
for example:

```text
Gene manifest does not describe the payloads that were written. Declared but
absent: c("var/7.values.u8.gz"). Written but undeclared: c("var/8.values.u8.gz").
```

That is what makes the index space trustworthy: a stale file from an earlier
generation cannot survive beside a current manifest.

The root is checked for a reason of its own. `points_<dim>d.bin(.gz)` is the only
payload this format declares by path *from the export root* — in
`dataset_identity.json`, under `embeddings.files` — so no axis manifest can speak
for it, and its declared name and its written name are two expressions of one
`compression` setting. After a successful export the root holds exactly
`dataset_identity.json`, `obs_manifest.json`, the declared point payloads,
whichever of `var_manifest.json` and `connectivity_manifest.json` the export
wrote, and the payload directories it created — and nothing else. A leftover
scratch file, or a directory the export did not create, is refused rather than
published.

`prepare()` in the `cellucid` Python package performs the same five
reconciliations and refuses the same generations. Only the rendering of the two
lists differs, as `['var/7.values.u8.gz']` rather than `c("var/7.values.u8.gz")`.

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
   - the manifest-declared source, destination, and Float64 weight payloads
     exist

For deeper validation and binary inspection, see:
- {doc}`../d_viewing_loading/03_validate_exports_and_debug_loading`
