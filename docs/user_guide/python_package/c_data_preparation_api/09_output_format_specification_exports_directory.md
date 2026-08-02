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

## Binary payload conventions (apply to every payload file)

Every binary file an export contains — `points_<dim>d.bin`, everything under
`obs/`, `var/`, `connectivity/`, and `vectors/` — obeys the same three rules.
None of them is recorded anywhere in the export, so they hold unconditionally:

- **Little-endian.** Every multi-byte value is stored least-significant byte
  first, whatever the byte order of the machine that produced the export.
  `uint8` payloads are single bytes and have no byte order.
- **Row-major (C order).** A payload of shape `(n_cells, dim)` stores all of
  cell 0's components, then all of cell 1's, and so on.
- **Raw and headerless.** No magic number, no length prefix, no padding, no
  alignment gap. The file is exactly `n_elements × itemsize` bytes; the element
  count and the dtype come from the manifest.

Byte order is the rule a reader cannot check, which is why it is stated here
once and normatively rather than left to each payload section. The dtype
strings the manifests publish — `float32`, `float64`, `uint8`, `uint16`,
`uint32` — name a width and a scalar kind and carry **no byte-order
component**, and the web app constructs its typed arrays directly on the bytes
it received, which is host order by definition in JavaScript. A big-endian
payload therefore does not fail to load. It loads as different numbers, and the
coordinates, expression values, and edge weights that result are entirely
plausible. There is no field to check and no error the reader could raise.

Both writers pin the order explicitly rather than inheriting the host's:
{func}`~cellucid.prepare` converts every payload at the write, and
`cellucid_prepare()` names `endian = "little"` at every `writeBin`. An export
produced on a big-endian machine — s390x is the realistic case — is therefore
byte-identical to one produced on x86 or ARM from the same input.

Anything else that writes this format must do the same. Emitting the host's
native order is correct on x86 and ARM by accident, not by contract, and the
mistake is invisible until an export crosses architectures.

---

## Payload naming: a filename never carries dataset content

Every scientific payload file in an export is named by its **integer index**, or
by a fixed neutral name. No filename anywhere in the tree contains an
observation key, a gene name, or a vector field id:

```text
var/0.values.u8.gz            var/1.values.u8.gz            …
obs/0.values.u8.gz            obs/1.codes.u8.gz             obs/1.outliers.u8.gz
vectors/0_2d.bin.gz           vectors/0_3d.bin.gz
connectivity/edges.src.bin.gz  edges.dst.bin.gz  edges.weights.f64.bin.gz
points_2d.bin.gz              points_3d.bin.gz
```

Every manifest path pattern therefore substitutes `{index}`, never a key, and
every field entry declares its own index as element `[0]`.

**Within one axis directory the indices are exactly `0 … N-1`, each used once.**
The index is assigned by position: for `var` it is the position in the exported
gene subset, for `obs` the position in the exported column selection, for
vectors the position in the sorted field list.

`obs` shares **one** index space across `_continuousFields` and
`_categoricalFields`, because both arrays write into `obs/`. That space follows
the `obs` column order, so it interleaves the two arrays — a manifest whose
first categorical entry declares index `1` and whose first continuous entry
declares index `0` is normal and correct. Two entries holding the same index
would write one payload over another and the app would draw one field's values
under another field's name, so both writers assert the rule against the manifest
they have just built, and re-derive every declared path to compare it with the
directory they actually wrote.

**Never infer a field's kind from its entry length.** A `var` entry and an
obs-continuous entry are both length 4 when quantized and both length 2 when
not. The array an entry came from is what says which it is.

(python_package-quantized-continuous-payloads)=
### How a quantized continuous payload decodes

Every quantized continuous payload in the format — an obs-continuous entry, a
`var` entry, and a categorical entry's outlier quantiles — carries a trailing
`minValue`/`maxValue` pair and decodes the same way:

```text
maxQuant = 254   (8-bit)   or  65534  (16-bit)
value    = minValue + code * (maxValue - minValue) / maxQuant
```

`minValue` and `maxValue` are both finite JSON numbers, and `minValue <= maxValue`.
The division is by `maxQuant`, which is a nonzero constant, so nothing in the
decode divides by the bounds.

#### The constant case: `minValue == maxValue`

A continuous field whose every value is the same value is published with
**equal bounds and every code `0`**. This is a named case of the format, not a
degenerate accident:

- A gene detected in no exported cell — routine once an atlas is subset to one
  lineage — is published like any other gene, at `minValue == maxValue == 0`.
- An obs column that a subset flattened, or a source column that genuinely
  carries one value, is published the same way.

Writers **must** take an explicit branch for it and never derive a scale;
readers **must** take the matching branch and return `minValue` directly. The
general formula above is exact for this case as well — `code` is `0` and the
numerator is `0` — so a constant field decodes back to its exact value in
every reader, bit for bit, rather than to an approximation of it.

Both writers implement this: `_quantize_continuous()` in
`cellucid/prepare_data/_quantization.py` and `.quantize_continuous()` in
`R/quantization.R`, each guarded by a single named predicate
(`_is_constant_continuous_range` / `.is_constant_continuous_range`).

A **non**-constant payload is the complementary case, and its two terminal
codes are both used: `minValue` and `maxValue` are the payload's own extremes,
so `0` and `maxQuant` each occur at least once.

The one continuous payload with no encoding at all is a categorical entry's
generated outlier quantiles when *every* quantile is missing — no category held
`centroid_min_points` cells, so there is no value to publish and no bound to
declare. Both writers reject that export and name every affected field.

### What an identifier still has to be

Because an identifier is no longer a path, it carries no filename rule at all:
no ASCII restriction, no length budget, no reserved Windows device name, no
case-insensitive collision rule. `HLA-DRB1/2`, `CON`, `% mito`, `细胞`, and the
two distinct fields `Field` and `field` all export unchanged.

Three rules remain, and each of them is about what the identifier is *for*:

- **non-empty** — it names a field;
- **distinct within its axis** — so a lookup resolves exactly one field;
- **drawable exactly as stored** — no control characters, none of the
  zero-width characters `U+200B`, `U+2060`, `U+FEFF`, and no leading or trailing
  whitespace. This is the same rule every string category label obeys, for the
  same reason: a gene name and a category label are drawn in the same legend.
  See [Category labels must read as the value they store](python_package-category-label-display-text).

The display rule covers exactly the identifiers this directory contains. An
`obs` column left out by `obs_keys`, or a `var` gene left out by
`gene_identifiers`, produces no manifest entry and no payload file, so it is not
held to it. Gene names are additionally required to be distinct across all of
`var`, which is what makes `gene_identifiers` resolve one row per name.

`dataset_id` is the exception, and it is not a field identifier: it names a real
directory in a multi-dataset export root and served URL, so it keeps the full
portable-component rule (see {doc}`02_input_requirements_global`).

### Gene names come from you, not from Cellucid

`var_gene_id_column` selects the column whose values become the exported gene
names, and both writers record whatever it selects, faithfully. Neither writer
performs a symbol lookup, ships a mapping, or knows what HGNC is. If you want
HGNC symbols rather than Ensembl accessions, resolve them in your own
preprocessing, write the result into a `var` column, and pass that column.

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

Optional keys exist in exactly two places, both in `dataset_identity.json`. At
its top level, `created_at`, `export_settings`, `source`, and `vector_fields`
may be absent. Inside the nested `source` object, `name` is required while `url`
and `citation` are each optional, because both writers write them only when the
corresponding `source_*` argument was supplied. Everything else in the format —
every other key of `dataset_identity.json`, every key of every other manifest
object — is required whenever its containing object is present.

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
    "n_obs_fields": 2,
    "n_categorical_fields": 1,
    "n_continuous_fields": 1,
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
    "name": "Source dataset name",
    "url": "https://...",
    "citation": "..."
  },
  "vector_fields": {
    "default_field": "velocity_umap",
    "fields": {
      "velocity_umap": {
        "label": "velocity_umap",
        "available_dimensions": [2, 3],
        "default_dimension": 3,
        "files": {
          "2d": "vectors/0_2d.bin.gz",
          "3d": "vectors/0_3d.bin.gz"
        },
        "basis": "umap"
      }
    }
  }
}
```

Notes:
- `vector_fields` is present only if you exported vectors. Inside a field
  entry, `label` and `basis` are **derived by the writer, not chosen by you**:
  `label` is always the field id character for character, and `basis` is always
  the string `"umap"`. The reader compares both and rejects the dataset when
  either differs, so `"label": "Velocity (UMAP)"` on a field called
  `velocity_umap` fails to load. There is no display-name field in this format.
- `source` is present only if you provided `source_*` metadata, and inside it
  only `name` is required. `url` and `citation` are each written only when you
  supplied that argument, so `"source": {"name": "Cellucid browser CI"}` is a
  complete and valid `source` object. The reader also accepts a `source.filename`
  key, which is **not part of the export format**: the app's in-browser
  H5AD/Zarr adapters set it to record the file a dataset was opened from.
  Neither writer emits it, and an export must not contain it.
- `stats.n_obs_fields` must equal `stats.n_categorical_fields +
  stats.n_continuous_fields`, and `obs_fields` must hold exactly
  `stats.n_obs_fields` elements. The example above is one complete small
  dataset, not an excerpt: two obs fields, one of each kind.
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
- rounding: the extents and the normalization are computed at the input's own
  precision, and the normalized coordinate is rounded to `float32` exactly
  once, at the write. Both writers do this, so one embedding produces one
  `points_<dim>d.bin` whichever writer produced it

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
      "pathPattern": "obs/{index}.values.u8.gz",
      "ext": "u8",
      "dtype": "uint8",
      "quantized": true,
      "quantizationBits": 8
    },
    "categorical": {
      "codesPathPattern": "obs/{index}.codes.{ext}.gz",
      "outlierPathPattern": "obs/{index}.outliers.u8.gz",
      "outlierExt": "u8",
      "outlierDtype": "uint8",
      "outlierQuantized": true
    }
  }
}
```

Notes:
- `{index}` is the entry's own element `[0]`, written in plain decimal with no
  padding: `0`, `1`, … `N-1`.
- `{ext}` for codes is chosen based on the field’s codes dtype (`uint8 → u8`, `uint16 → u16`).
- `_obsSchemas` is always a JSON **object**, including when it is empty — an
  export with no observation fields writes `{}`, never `[]`. The reader refuses
  an array. This is worth stating because the two writers reach it through
  different serialisers, and a language whose empty map and empty list share a
  representation will emit the wrong one unless it is forced.
- `_obsSchemas` is otherwise an exact key set: it holds `continuous` only when
  `_continuousFields` is non-empty, `categorical` only when
  `_categoricalFields` is non-empty, and nothing else.
- `quantizationBits` appears on the continuous schema only when
  `"quantized": true`.

#### `_continuousFields`

Each entry is either:
- `[index, key]` — **2 members**, non-quantized float32, or
- `[index, key, minValue, maxValue]` — **4 members**, quantized, with
  `minValue <= maxValue`; see
  [How a quantized continuous payload decodes](python_package-quantized-continuous-payloads)

#### `_categoricalFields`

Each entry is either:
- `[index, key, categories, codesDtype, codesMissingValue, centroidsByDim]` — **6 members**, float32 outliers, or
- `[index, key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]` — **8 members**, quantized outliers, with
  `outlierMinValue <= outlierMaxValue`; see
  [How a quantized continuous payload decodes](python_package-quantized-continuous-payloads)

Where:
- element `[0]` is this field's payload index in the single `obs/` index space
  shared with `_continuousFields`
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
- `codesMissingValue` is `255` or `65535`. That value is reserved for "missing",
  so it caps how many categories a field may declare: **at most 255 categories
  under `"uint8"` codes and at most 65,535 under `"uint16"`**. Both writers
  reject an export whose field exceeds the cap for the requested
  `obs_categorical_dtype` and name the offending field.
- `centroidsByDim` is a dict keyed by dimension strings (`"1"`, `"2"`, `"3"`)
  mapping to centroid lists. Each centroid is an exact three-key object
  `{"category": <declared category>, "position": [<dim finite numbers>],
  "n_points": <0..n_points>}`, with at most one entry per category.

The two arrays are read in order — every `_continuousFields` entry first, then
every `_categoricalFields` entry — and that concatenation is the order
`dataset_identity.json` `obs_fields` must reproduce exactly. That reading order
is independent of the payload indices: element `[0]` says which file a field owns,
and array position says where the field appears in the field list.

### Binary payload files under `obs/`

Per continuous field:
- `obs/<index>.values.f32` or `.values.u8` or `.values.u16` (optionally `.gz`)

Per categorical field:
- `obs/<index>.codes.u8` or `.codes.u16` (optionally `.gz`)
- `obs/<index>.outliers.f32` or `.outliers.u8` or `.outliers.u16` (optionally `.gz`)

Quantized outlier/categorical missing markers:
- 8-bit: `255`
- 16-bit: `65535`

Worked example — four `obs` columns declared in the order
`leiden` (categorical), `pct_mito`, `total_counts`, `phase` (categorical):

```json
{
  "_continuousFields": [[1, "pct_mito"], [2, "total_counts"]],
  "_categoricalFields": [
    [0, "leiden", ["0", "1"], "uint8", 255, {"2": []}],
    [3, "phase", ["G1", "S"], "uint8", 255, {"2": []}]
  ]
}
```

```text
obs/0.codes.u8      obs/0.outliers.f32     # leiden
obs/1.values.f32                           # pct_mito
obs/2.values.f32                           # total_counts
obs/3.codes.u8      obs/3.outliers.f32     # phase
```

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
    "pathPattern": "var/{index}.values.f32.gz",
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
    "pathPattern": "var/{index}.values.u8.gz",
    "ext": "u8",
    "dtype": "uint8",
    "quantized": true,
    "quantizationBits": 8
  }
}
```

#### `fields`

Each entry is either:
- `[index, name]` — **2 members**, non-quantized, or
- `[index, name, minValue, maxValue]` — **4 members**, quantized, with
  `minValue <= maxValue`; see
  [How a quantized continuous payload decodes](python_package-quantized-continuous-payloads)

```json
{
  "fields": [[0, "TNMD"], [1, "CIITA"], [2, "HLA-DRB1/2"]]
}
```

Element `[0]` is the position in the **exported** gene set, so a narrowed export
is still `0 … n_genes-1` and never inherits the row numbers of `var`. Element
`[1]` is whatever `var_gene_id_column` selected, recorded verbatim.

`var` entries and obs-continuous entries deliberately share a shape. Read the
kind from the array the entry came from, never from its length.

### Binary payload files under `var/`

For each exported gene:
- `var/<index>.values.f32` or `.values.u8` or `.values.u16` (optionally `.gz`)

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
- the endpoints are unsigned integers of `index_bytes` width and the weights are
  IEEE-754 binary64, in the little-endian row-major encoding every payload file
  in the export uses
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

### Browser working-set ceiling (`n_edges`)

Neither writer caps `n_edges`, but the reader does. It budgets a fixed
**536,870,912 bytes (512 MiB)** working set for the edge list and rejects the
manifest before fetching a single edge byte when the declared graph would exceed
it. The budget charges every copy the browser holds at once — the canonical
endpoint and weight arrays, the render-owned copy of both, the GPU topology and
Float32 weight staging buffers, one `uint32` degree counter per cell, and, when
`index_bytes` is not `4`, the raw file buffer as well:

| `index_dtype` | bytes charged |
| --- | --- |
| `"uint32"` | `44 × n_edges + 4 × n_cells` |
| `"uint16"` | `48 × n_edges + 4 × n_cells` |

That is a ceiling of about **12.2 million edges** for a `uint32` export and about
**11.1 million** for a `uint16` one. The largest graph shipped in
`cellucid-datasets` is `suo` — 6,279,148 edges over 561,947 cells, charging
about 266 MiB — so every shipped dataset clears the limit with room to spare. An
export above the line is a well-formed export that the web app still refuses to
open, so thin the neighbour graph before exporting rather than after.

---

## Vector files: `vectors/<index>_<dim>d.bin(.gz)` (optional)

Vector fields are stored as:
- dtype: `float32`
- shape: `(n_cells, dim)`
- files: `vectors/<index>_<dim>d.bin(.gz)`
- rounding: the caller's values are scaled by the embedding's own scale factor
  at the input's precision, and the product is rounded to `float32` exactly
  once, at the write. Rounding the input first and the product second rounds
  twice, which moves a component one ULP away from the correctly rounded value;
  neither writer does that

`<index>` is the field's position in the sorted field list, `0 … N-1`, exactly
as on the other axes; the field id itself never appears in a path. Vector
presence and file locations are indexed in
`dataset_identity.json["vector_fields"]`, whose per-field `files` object states
the complete path.

That path is **not** free-form, and this is the rule a hand-written or converted
export trips on most often. The reader rebuilds each path from the manifest —

```text
vectors/{index}_{dim}d.bin      when export_settings.compression is null
vectors/{index}_{dim}d.bin.gz   otherwise
```

— where `{index}` is the field's position in the key order of the
`vector_fields.fields` object, and it rejects the dataset when the declared
string is anything else. A path that resolves to a real file is not enough: it
must be the exact string the producer rule generates. Two further rules come
with it:

- `available_dimensions` is strictly increasing and is a subset of
  `embeddings.available_dimensions` — a vector field cannot advertise a
  dimension the dataset has no points file for — and `files` holds exactly one
  key per entry in it, named `"1d"`, `"2d"`, or `"3d"`, and no others;
- `default_dimension` must be the **largest** entry of `available_dimensions`,
  not merely one of them. Both writers emit `max(available_dimensions)`. (The
  top-level `embeddings.default_dimension` is the looser case: there, any
  available dimension is accepted.)

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

## The writer proves the layout before publishing

A candidate generation is built in a staging directory and is reconciled against
its own declarations before the transaction publishes it. {func}`~cellucid.prepare`
checks **five** surfaces: the four payload directories — `Observation`, `Gene`,
`Connectivity`, and `Vector field` — and the `Export` root itself. On each one,
every path the export declares must exist and every file present must be
declared. A mismatch aborts the run and rolls the whole generation back, so a
partial or self-contradictory export is never published:

```text
Gene manifest does not describe the payloads that were written. Declared but
absent: ['var/7.values.u8.gz']. Written but undeclared: ['var/8.values.u8.gz'].
```

That is what makes the index space trustworthy: a stale file from an earlier
generation cannot survive beside a current manifest.

The root is checked for a reason of its own. `points_<dim>d.bin(.gz)` is the only
payload this format declares by path *from the export root* — in
`dataset_identity.json`, under `embeddings.files` — so no axis manifest can speak
for it. Its declared name and its written name are two expressions of one
`compression` setting, and a generation whose points disagree would otherwise
publish successfully and then fail in the browser with no coordinates:

```text
Export manifest does not describe the payloads that were written. Declared but
absent: ['points_2d.bin.gz']. Written but undeclared: ['points_2d.bin'].
```

The root check covers the whole directory rather than only the points. After a
successful export the root holds exactly `dataset_identity.json`,
`obs_manifest.json`, the declared point payloads, whichever of
`var_manifest.json` and `connectivity_manifest.json` the export wrote, and the
`obs/`, `var/`, `connectivity/`, and `vectors/` directories it created — and
nothing else. A leftover scratch file or a directory the export did not create is
refused rather than published. A directory standing where a payload file belongs
gets its own message, `Export payload directory holds a non-file entry: …`, which
is also what an axis directory raises for a subdirectory inside it.

`cellucid_prepare()` in the R package performs the same five reconciliations and
refuses the same generations. Only the rendering of the two lists differs, as
`c("var/7.values.u8.gz")` rather than `['var/7.values.u8.gz']`.

---

## Determinism and reproducibility

Within one writer, scientific binary payloads, including their gzip bytes, are
deterministic given:
- the exact input arrays and DataFrames,
- the same exporter version and compression level,
- the same zlib build underneath it.

### What the two writers guarantee about each other

The cross-writer guarantee is **semantic equality of every export, plus byte
identity of the uncompressed binary payloads**. One input written by
`cellucid.prepare()` and by `cellucid_prepare()` produces:

- **byte-identical `points_<dim>d.bin`, `vectors/<index>_<dim>d.bin`, `obs/`,
  `var/`, and connectivity payloads** when `compression=None` — same bytes, same
  length, same order. The single float32 rounding described above is what makes
  that true of the coordinate and vector payloads, and the pinned little-endian
  byte order is what makes it true on any architecture rather than only on the
  little-endian ones;
- **equal** `dataset_identity.json`, `obs_manifest.json`, `var_manifest.json`,
  and `connectivity_manifest.json` — the same keys with the same parsed values,
  which is what a reader consumes.

Byte identity does **not** extend to the JSON files as bytes. It does extend to
the compressed `.gz` payloads, but only as far as the two zlib builds agree.
Every known difference is listed here; none of them changes a value a reader
parses.

The **whole ten-byte gzip member header is identical between the two writers, on
every platform, unconditionally**: the magic `1f 8b`, deflate, no optional header
fields, `MTIME = 0`, the RFC 1952 §2.3.1 extra-flags byte for the level (`0x04`
at level 1, `0x02` at level 9, `0x00` in between), and `OS = 0xff`, the
"unknown" code. Python writes it through `gzip.GzipFile(filename="", mtime=0)`;
`cellucid_prepare()` replaces the header `gzfile()` produced with those same ten
bytes, because `gzfile()` would otherwise stamp the code of the platform zlib was
built for and `XFL = 0x00` at every level. The CRC32 and ISIZE trailer bytes
always match as well.

A **complete compressed member is byte-identical across the writers whenever the
two zlib builds agree** — verified for a 4 KiB float32 payload at all nine
compression levels. When the builds differ, only the deflate stream between that
header and that trailer differs. Measured with R linked against zlib 1.3.2 and
Python against zlib 1.2.12, one 4,495,576-byte payload compressed at level 6 to
4,160,196 bytes under R and 4,163,855 bytes under Python. That version pairing is
an example of build skew rather than a property of the packages: link both
languages against one zlib and the compressed bytes match. The cause is the zlib
build and not either writer, which two measurements pin down — within a single
zlib, chunked deflate equals one-shot deflate, which rules out R's 16 KiB
streaming, and the two writers agree at every level on the 4 KiB payload, which
could not happen if their `deflateInit2` parameters differed.

| # | Where | Python writes | R writes |
| --- | --- | --- | --- |
| 1 | end of every JSON file | no trailing newline (`write_text(json.dumps(...))`) | one trailing `\n` (`writeLines()`) |
| 2 | `dataset_identity.json` layout | `json.dumps(indent=2)` | `jsonlite::prettify(indent = 2)` — a different pretty-printer, with its own array and whitespace layout |
| 3 | separators in the three compact manifests | `, ` and `: ` | `,` and `:` |
| 4 | non-ASCII text in any JSON file | escaped as `\uXXXX` | written as UTF-8; both decode to the same string |
| 5 | categorical centroid `position` | — | may differ from Python in the last float64 digit, because the two languages sum a category's coordinates in a different order |

Difference 5 is the only numeric one; it is a JSON double either way and is ~9
orders of magnitude below the float32 precision the coordinates are drawn at.

To compare two exports across writers, compare parsed JSON rather than JSON
bytes, and compare the payloads either uncompressed or after decompression — that
comparison holds whichever zlib each language was built against.

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
- **Case sensitivity**: `dataset_id` and the directory names in a multi-dataset
  root are real paths and behave differently on macOS vs Linux vs Windows;
  prefer lowercase IDs. Field identifiers are not paths and are unaffected.
- **Renumbering on re-export**: an index is a position, not a stable name.
  Adding one gene or reordering `obs_keys` renumbers the payload files, so
  never bookmark, cache, or hard-code a payload path across exports — resolve it
  through the manifest every time.
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
  `connectivity_manifest.json` `n_edges`, `stats.n_obs_fields` vs the
  `obs_fields` array length or vs
  `n_categorical_fields + n_continuous_fields`, or `obs_fields` vs the obs
  manifest field order.

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
