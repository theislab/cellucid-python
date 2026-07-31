# Export format spec and invariants

This page is the **specification of the Cellucid export folder format** written by `cellucid.prepare(...)`.

It is the contract between:
- `cellucid-python` (writer/server), and
- the Cellucid web app (reader).

If you change this format, you must coordinate with the web app and document compatibility clearly.

Related pages:
- Export pipeline internals: {doc}`07_prepare_export_pipeline_architecture`
- Server behavior: {doc}`09_server_mode_architecture_endpoints_and_security`

---

## What counts as a “valid export folder”

In practice, there are two “levels”:

### Level 1: Minimal viewer-loadable dataset

- `obs_manifest.json` exists
- at least one `points_<dim>d.bin` (or `.bin.gz`) exists

This is what the static server uses to detect datasets.

### Level 2: “Produced by `prepare`” (recommended)

Everything above, plus:

- `dataset_identity.json` exists (used by CLI auto-detection and dataset display)
- optional: `var_manifest.json`, `connectivity_manifest.json`, `vectors/`

Unless you are deliberately building custom exports, treat `dataset_identity.json` as required.

---

## Directory layout (canonical)

```text
my_export/
  dataset_identity.json
  obs_manifest.json
  var_manifest.json                    # optional
  connectivity_manifest.json           # optional
  points_1d.bin(.gz)                   # optional
  points_2d.bin(.gz)                   # optional
  points_3d.bin(.gz)                   # optional
  obs/
    <index>.values.f32(.gz)            # continuous, not quantized
    <index>.values.u8(.gz)             # continuous, quantized 8-bit
    <index>.values.u16(.gz)            # continuous, quantized 16-bit
    <index>.codes.u8(.gz)              # categorical codes
    <index>.codes.u16(.gz)
    <index>.outliers.f32(.gz)          # per-cell outlier quantiles (latent space)
    <index>.outliers.u8(.gz)           # quantized outliers (if enabled)
    <index>.outliers.u16(.gz)
  var/
    <index>.values.f32(.gz)            # gene expression, not quantized
    <index>.values.u8(.gz)             # quantized gene expression
    <index>.values.u16(.gz)
  connectivity/
    edges.src.bin(.gz)                 # integer indices
    edges.dst.bin(.gz)
    edges.weights.f64.bin(.gz)         # little-endian Float64 weights
  vectors/
    <index>_1d.bin(.gz)
    <index>_2d.bin(.gz)
    <index>_3d.bin(.gz)
```

Notes:
- `<index>` is the field's own payload index, a plain unpadded decimal integer
  in `0 … N-1` (see below). No payload filename contains an identifier.
- `.gz` files are **raw gzip files** (the viewer decompresses them; servers typically do not set `Content-Encoding` for these).

---

## Payload indices and identifiers

Payload filenames are integer indices, so an identifier is never a path.

Within one axis directory the indices are exactly `0 … N-1`, each used once,
assigned by position in the exported selection. `obs` shares that one space
across `_continuousFields` and `_categoricalFields`, because both write into
`obs/`. Both writers assert this against the manifest they have just built and
re-derive every declared payload path to compare it against the directory they
actually wrote; a mismatch rejects the complete candidate before publication.

Because it is not a path, an identifier carries no filename rule. `prepare`
still requires each exported observation key, gene name, and vector-field id to
be a non-empty string, distinct within its axis, and drawable exactly as stored
(no control characters, no `U+200B`/`U+2060`/`U+FEFF`, no leading or trailing
whitespace) — the rule string category labels obey, since all of them are drawn
in the same legend. Nothing is ever sanitized: an identifier that breaks the
rule rejects the complete candidate before publication.

`dataset_id` is the exception. It names a real directory in a multi-dataset
export root and in a served URL, so it remains a 1–180-byte portable ASCII
filename component: it starts with a letter or digit, contains only letters,
digits, `.`, `_`, or `-`, does not end with `.`, and is not a Windows device
name.

The contract is scoped to the identifiers that become paths. `obs_keys` and
`gene_identifiers` narrow what is exported and narrow this contract with it, so
an `obs` column or `var` gene that is left out is never checked against it. Gene
IDs carry one additional, unnarrowed rule: every ID in `var` must be a non-empty
string and must be distinct, since `gene_identifiers` addresses `var` rows by
identifier and a repeated ID resolves to no single row.

---

## Compression semantics (do not confuse these)

There are two gzip mechanisms in the ecosystem:

### A) Export-folder gzip (`.gz` suffix)

If `prepare(compression=...)` is enabled:
- files are written as `*.gz` on disk (raw gzip stream)
- manifests point to those `*.gz` file paths
- the viewer is expected to fetch those bytes and decompress in JS

### B) HTTP gzip (`Content-Encoding: gzip`)

In AnnData server mode, if the browser sends `Accept-Encoding: gzip`:
- the server may gzip-compress the response *on the fly*
- and set `Content-Encoding: gzip`
- the browser transparently decompresses before JS sees bytes

Export folders generally do **not** rely on `Content-Encoding`.

---

## Binary file formats

All binaries are raw, headerless arrays. You must know the shape from context + manifests.

Unless stated otherwise:
- numeric endianness is the machine-native numpy endianness (practically little-endian on modern systems)
- the viewer assumes little-endian

### `points_<dim>d.bin`

- dtype: `float32`
- shape: `(n_cells, dim)`
- value range: normalized to approximately `[-1, 1]`

### `vectors/<index>_<dim>d.bin`

- dtype: `float32`
- shape: `(n_cells, dim)`
- scaled into the same normalized coordinate system as the points

### `obs/<index>.values.*`

Continuous field values, either:

- `float32` (`.values.f32`)
- quantized `uint8` (`.values.u8`) or `uint16` (`.values.u16`)

Quantized encoding:
- reserve max value for an explicitly absent quantized value (for example, an
  undefined categorical outlier quantile):
  - u8: `255` is missing; `0..254` are valid quantized values
  - u16: `65535` is missing; `0..65534` are valid quantized values
- manifests store `minValue` and `maxValue` for dequantization, with
  `minValue <= maxValue`
- a field whose values are all the same is the format's constant case:
  `minValue == maxValue` and every code `0`, decoded by returning `minValue`
  directly. See
  {ref}`How a quantized continuous payload decodes <python_package-quantized-continuous-payloads>`.

### `obs/<index>.codes.*`

Categorical codes:
- dtype: `uint8` or `uint16`
- values:
  - `0..(n_categories-1)` are category indices
  - missing marker is stored in the manifest (255 or 65535)

### `obs/<index>.outliers.*`

Per-cell outlier quantiles computed in latent space:
- either `float32` or quantized `uint8/uint16`
- same missing-marker rules for quantized variants

### `var/<index>.values.*`

Gene expression values (one file per gene):
- either `float32` or quantized `uint8/uint16`
- quantized variants use the same missing-marker rules as continuous obs values

### Connectivity edge payloads

Three aligned arrays of equal length:
- dtype: `uint16` for at most 65,536 cells or `uint32` for at most
  4,294,967,296 cells (see manifest); `uint64` is not part of the current
  contract
- semantics: undirected edges stored once as `(src, dst)` with `src < dst`
- arrays are sorted lexicographically by `(src, dst)`
- `edges.weights.f64.bin(.gz)` stores the corresponding exact little-endian
  Float64 weight
- the source graph must already be finite, non-negative, exactly symmetric in
  topology and weight, and exactly zero on the diagonal
- sparse inputs must contain neither stored zeros nor duplicate coordinates;
  producers reject rather than reinterpret them

---

## JSON files

### `dataset_identity.json`

This is the main metadata file written by `prepare`.

Key fields (current):

- `version`: integer (currently `2`)
- `id`, `name`, `description`, `created_at`
- `cellucid_data_version`: version string of the Python package writing the export
- `stats`: counts + connectivity presence
- `embeddings`: available dimensions and file paths
- `obs_fields`: simplified obs field summary for UI display
- `export_settings`: compression + quantization knobs used
- optional `source`: name/url/citation
- optional `vector_fields`: field metadata + per-dimension file paths

### `obs_manifest.json` (format: `compact_v1`)

Fields:

- `_format`: `"compact_v1"`
- `n_points`: number of cells
- `centroid_outlier_quantile`: float or null
- `latent_key`: currently `"latent_space"` (conceptual key)
- `compression`: gzip level or null
- `_obsSchemas`: schema object that defines path patterns and dtypes
- `_continuousFields`: compact list of continuous fields
- `_categoricalFields`: compact list of categorical fields

Compact lists:

- continuous:
  - `[key]` (float32, unquantized)
  - `[key, minValue, maxValue]` (quantized)
- categorical:
  - `[key, categories, codesDtype, codesMissingValue, centroidsByDim]` (outliers float32)
  - `[key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]` (outliers quantized)

Centroids:
- `centroidsByDim` is a dict keyed by dimension (as strings), each value is a list of objects:
  - `{category: <str>, position: <list[float]>, n_points: <int>}`

### `var_manifest.json` (format: `compact_v1`)

Fields:

- `_format`: `"compact_v1"`
- `n_points`: number of cells
- `var_gene_id_column`: `null` for `var.index`, otherwise an exact column name
- `compression`: gzip level or null
- `quantization`: `8`, `16`, or null
- `_varSchema`: schema object with path patterns
- `fields`: list of compact per-gene entries:
  - `[gene_id]` (float32)
  - `[gene_id, minValue, maxValue]` (quantized)

### `connectivity_manifest.json` (format: `edge_pairs`)

Fields:

- `format`: `"edge_pairs"`
- `n_cells`, `n_edges`, `max_neighbors`
- `index_dtype` (`uint16`/`uint32`) and `index_bytes`
- `sourcesPath`, `destinationsPath`
- `compression`: gzip level or null

An explicitly supplied, validated zero-edge graph has `n_edges: 0`,
`max_neighbors: 0`, and two zero-length edge arrays. It is distinct from absent
connectivity.

---

## Reading files in Python (quick examples)

### Read `points_3d.bin(.gz)`

```python
import gzip
import json
import numpy as np
from pathlib import Path

export_dir = Path("my_export")
identity = json.loads((export_dir / "dataset_identity.json").read_text("utf-8"))
n = identity["stats"]["n_cells"]
path = export_dir / identity["embeddings"]["files"]["3d"]
raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
pts = np.frombuffer(raw, dtype=np.float32).reshape(n, 3)
```

### Read a categorical obs codes file

```python
import gzip
import json
import numpy as np
from pathlib import Path

export_dir = Path("my_export")
obs_manifest = json.loads((export_dir / "obs_manifest.json").read_text("utf-8"))

# Find one categorical field entry: [index, key, categories, dtype, ...]
field = obs_manifest["_categoricalFields"][0]
payload_index = field[0]
key = field[1]
codes_dtype = field[3]  # "uint8" or "uint16"
codes_ext = "u8" if codes_dtype == "uint8" else "u16"

path = export_dir / "obs" / f"{payload_index}.codes.{codes_ext}"
raw = gzip.open(str(path) + ".gz", "rb").read() if (path.with_suffix(path.suffix + ".gz")).exists() else path.read_bytes()
codes = np.frombuffer(raw, dtype=np.uint8 if codes_ext == "u8" else np.uint16)
```

---

## Invariants (the “do not break” list)

1) `n_points` must match the number of rows in all per-cell binaries.
2) Category codes must be stable and match the `categories` list in the manifest.
3) Categorical codes and generated nullable categorical outlier quantiles must
   respect their declared missing markers. Gene and continuous-observation
   values are finite-only.
4) If you change any schema/filenames, coordinate with the web app.
5) `.gz` means “raw gzip bytes”, not HTTP `Content-Encoding`.

---

## Troubleshooting

### Symptom: “Viewer loads, but fields are missing”

Likely causes:
- `obs_manifest.json` references a file that isn’t present (partial export),
- identifier/path mismatch,
- compression mismatch (`.gz` expected but missing).

Fix:
- re-export with `force=True`, or delete and re-run `prepare`.

### Symptom: “Continuous field looks blank”

Likely causes:
- the valid finite values occupy a very narrow portion of the declared
  min/max range,
- the active filters hide the cells carrying the visible range.

Confirm:
- inspect min/max in the manifest for that field,
- check the raw values upstream in Python.
