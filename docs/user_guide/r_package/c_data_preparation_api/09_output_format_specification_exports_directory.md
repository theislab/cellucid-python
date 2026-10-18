# Export Directory Format (Specification)

**Audience:** advanced users, maintainers, anyone debugging exports  
**Time:** 20–40 minutes  
**Goal:** know what every file means and how to validate an export folder.

This page is a practical specification of what `cellucid-r` writes to `out_dir`.

```{note}
This describes the export format as produced by `cellucid-r` today, in the terms
an R caller works in. The format itself is one contract shared by both writers
and the viewer, and its normative specification —
{doc}`../../python_package/c_data_preparation_api/09_output_format_specification_exports_directory`
— is where a reimplementation, or a question this page does not answer, should
go. Nothing here is R-specific except the examples: an export written by
`cellucid_prepare()` and one written by `cellucid.prepare()` are byte-comparable.
```

## Top-level layout

At the top level you will typically see:

- `dataset_identity.json` (always)
- `embeddings/` (always; one `<index>_<dim>d.bin` per embedding, optionally `.gz`)
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

## Embeddings: `embeddings/<index>_<dim>d.bin`

Files: one per entry of `embeddings`, for example
- `embeddings/0_1d.bin`
- `embeddings/1_2d.bin`
- `embeddings/2_2d.bin`
- `embeddings/3_3d.bin`

`<index>` is the embedding's position after sorting by dimension and then by id
in code-point order (`sort(method = "radix")`, which is what Python's `sorted()`
does), so both writers name the same embedding with the same file.

Encoding:
- float32, little-endian
- row-major

Shape:
- `(n_cells, dim)`, where `dim` is the number of columns plotted

Notes:
- each embedding is normalized on its own (center + one shared scale); see
  {doc}`03_embeddings_and_coordinates`
- only the selected columns are written; `dataset_identity.json` records them
  0-based under `source.columns` (`cellucid_embedding(columns = c(4, 2))` is
  recorded as `[3, 1]`)
- every source coordinate must be finite *and* representable as a nonzero finite
  `float32`. That is checked before normalization, not at the write:
  normalization pulls every coordinate into roughly `[-1, 1]`, so a coordinate
  outside the `float32` range would arrive at the writer already inside it and
  be published as a number you never supplied. `latent_space` is checked the
  same way. The `cellucid` Python package accepts exactly the same inputs

## Obs manifest: `obs_manifest.json`

This JSON file tells the viewer:
- how many points exist (`n_points`)
- where obs binary files are (path patterns)
- which fields exist and whether they are categorical/continuous
- how categorical codes and outlier files are encoded

Important keys (high level):
- `_format`: format identifier (currently `"compact_v2"`)
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

- continuous: `[index, key]`, or `[index, key, minValue, maxValue]` when quantized
- categorical: `[index, key, categories, codes_dtype, codes_missing_value,
  centroidsByEmbedding]`, or the same six followed by
  `outlierMinValue, outlierMaxValue` when outliers are quantized

`centroidsByEmbedding` is an object keyed by embedding id — every exported
embedding, in file order — whose centroid positions have that embedding's
dimension.

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
- the manifest entry includes `minValue` and `maxValue`
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
- what the gene names are, in the default set and in every alternative set,
- how to find each gene’s values file,
- and (if quantized) the min/max for dequantization.

Important keys (all seven are required):
- `_format` (currently `"compact_v2"`)
- `n_points`
- `compression`
- `quantization` (null, 8, or 16)
- `_varSchema` (path pattern, dtype, quantization flags); the pattern is a
  template over `{index}`, for example `var/{index}.values.f32.gz`
- `fields`: list of gene entries
- `gene_name_sets`: the default gene-name set and every alternative set

Each `fields` entry is:
- `[index, gene_id]` for float32 exports, or
- `[index, gene_id, minValue, maxValue]` for quantized exports

The first element is the gene's payload index. Within `var/` the indices are exactly
`0 … N-1`, each used once, and `cellucid_prepare()` proves that before
publishing. The second is the gene's name in the default set: what
`rownames(var)` or `var_gene_id_column` holds, with an empty value (`""` or
`NA`) taking the row name and a repeated name made unique (`MATR3`, `MATR3-1`).

`gene_name_sets[1]` onwards hold the columns named in `var_gene_name_columns`.
Each set is `{source, filled_from_index, renamed, names}`: `source` is the
column (`null` for the row names in the default set), `filled_from_index` and
`renamed` report what the make-unique rule changed over the exported genes, and
`names` is `null` for the default set and otherwise one name per gene, `null`
where it equals the default name. The sets are metadata only — every set names
the same `var/` payloads. The normative definition is the
{doc}`shared specification <../../python_package/c_data_preparation_api/09_output_format_specification_exports_directory>`.

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

- `version` (currently 3)
- `id`, `name`, `description`
- `created_at`
- `cellucid_data_version` (R package version)
- `stats` (counts: cells, exported genes, fields, edges; `n_genes` respects `gene_identifiers`)
- `embeddings`:
  - `default_dimension` (the highest dimension exported)
  - `defaults` (the embedding each dimension shows first, e.g. `"2d": "X_umap"`)
  - `items`: one entry per embedding id, with `label`, `dimension`, `file`
    (e.g. `"embeddings/1_2d.bin.gz"`), and `source` (`key`, and the 0-based
    `columns`)
- `obs_fields`: summary list of exported obs fields
- `export_settings`: compression/quantization settings
- optional `source` (name/url/citation)
- optional `vector_fields` metadata

`name`, `description`, and every `source` string are shown to the reader
verbatim, so they obey the same display-text rule as a string category
label. `description` may be `""`; the others may not be empty.

## Vector fields: `vectors/` (optional)

Vector binaries are float32 matrices stored row-major, one per field and bound
embedding:

- `vectors/<i>_<e>.bin`, for example `vectors/0_1.bin` and `vectors/0_3.bin`

`<i>` is the field's payload index, assigned by sorting the field ids by code
point so that `cellucid-r` and the `cellucid` Python package assign the same
index to the same field, and `<e>` is the payload index of the embedding the
array is drawn over. Each array has that embedding's dimension and is scaled by
that embedding's normalization factor.

The web app discovers them via `dataset_identity.json` (`vector_fields` section),
which maps each field id to `{label, files}`, `files` keyed by embedding id.

Every component must be finite and representable as a nonzero finite `float32`.
The small end matters as much as the large one: a nonzero component below
`2^-149` converts to exactly `0.0`, which *is* finite, so a check that asked
only about finiteness would publish a vector whose component had silently become
zero.

## The writer proves the layout before publishing

Before a candidate generation is published, `cellucid_prepare()` reconciles what
is on disk against what the export declares, on **six** surfaces: the five
payload directories — `Observation`, `Gene`, `Connectivity`, `Vector field`, and
`Embedding` — and the `Export` root itself. On each one, every path the export declares must
exist and every file present must be declared. A mismatch stops the export with,
for example:

```text
Gene manifest does not describe the payloads that were written. Declared but
absent: c("var/7.values.u8.gz"). Written but undeclared: c("var/8.values.u8.gz").
```

That is what makes the index space trustworthy: a stale file from an earlier
generation cannot survive beside a current manifest.

The `Embedding` surface is declared from `dataset_identity.json`: every item's
`file` must be exactly `embeddings/<index>_<dim>d.bin` (plus `.gz` when
compressed) before the directory is compared with it. The root holds no payload
at all: after a successful export it holds exactly `dataset_identity.json`,
`obs_manifest.json`, whichever of `var_manifest.json` and
`connectivity_manifest.json` the export wrote, and the payload directories it
created — and nothing else. A leftover
scratch file, or a directory the export did not create, is refused rather than
published.

`prepare()` in the `cellucid` Python package performs the same six
reconciliations and refuses the same generations. Only the rendering of the two
lists differs, as `['var/7.values.u8.gz']` rather than `c("var/7.values.u8.gz")`.

## Validation checklist (quick)

Before sharing an export folder:

1) Confirm the identity file parses:
   - `jsonlite::read_json(file.path(out_dir, "dataset_identity.json"))`
2) Confirm required files exist:
   - every `embeddings/<index>_<dim>d.bin[.gz]` the identity declares
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
