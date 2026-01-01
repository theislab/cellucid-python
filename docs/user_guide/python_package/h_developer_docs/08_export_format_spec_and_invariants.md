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
  points_4d.bin(.gz)                   # optional (viewer may not support yet)
  obs/
    <field>.values.f32(.gz)            # continuous, not quantized
    <field>.values.u8(.gz)             # continuous, quantized 8-bit
    <field>.values.u16(.gz)            # continuous, quantized 16-bit
    <field>.codes.u8(.gz)              # categorical codes
    <field>.codes.u16(.gz)
    <field>.outliers.f32(.gz)          # per-cell outlier quantiles (latent space)
    <field>.outliers.u8(.gz)           # quantized outliers (if enabled)
    <field>.outliers.u16(.gz)
  var/
    <gene>.values.f32(.gz)             # gene expression, not quantized
    <gene>.values.u8(.gz)              # quantized gene expression
    <gene>.values.u16(.gz)
  connectivity/
    edges.src.bin(.gz)                 # integer indices
    edges.dst.bin(.gz)
  vectors/
    <field_id>_1d.bin(.gz)
    <field_id>_2d.bin(.gz)
    <field_id>_3d.bin(.gz)
```

Notes:
- `<field>` and `<gene>` are “safe filename components” derived from keys (see below).
- `.gz` files are **raw gzip files** (the viewer decompresses them; servers typically do not set `Content-Encoding` for these).

---

## Safe filenames: how keys map to paths

`prepare` uses a conservative safe-name function:

- allow: `A–Z`, `a–z`, `0–9`, `.`, `_`, `-`
- replace anything else with `_`
- strip leading/trailing `.` and `_`

Implications:
- the *file name* for an obs field or gene may not match the original key exactly,
- collisions are possible if two different keys sanitize to the same safe name.

Developer recommendation:
- avoid obs field names that differ only by punctuation/whitespace,
- and consider adding collision detection if you extend the format.

Vector field IDs are stricter:
- vector field ids must already be safe; unsafe ids raise.

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

### `vectors/<field_id>_<dim>d.bin`

- dtype: `float32`
- shape: `(n_cells, dim)`
- scaled into the same normalized coordinate system as the points

### `obs/<field>.values.*`

Continuous field values, either:

- `float32` (`.values.f32`)
- quantized `uint8` (`.values.u8`) or `uint16` (`.values.u16`)

Quantized encoding:
- reserve max value for missing (NaN/Inf):
  - u8: `255` is missing; `0..254` are valid quantized values
  - u16: `65535` is missing; `0..65534` are valid quantized values
- manifests store `minValue` and `maxValue` for dequantization

### `obs/<field>.codes.*`

Categorical codes:
- dtype: `uint8` or `uint16`
- values:
  - `0..(n_categories-1)` are category indices
  - missing marker is stored in the manifest (255 or 65535)

### `obs/<field>.outliers.*`

Per-cell outlier quantiles computed in latent space:
- either `float32` or quantized `uint8/uint16`
- same missing-marker rules for quantized variants

### `var/<gene>.values.*`

Gene expression values (one file per gene):
- either `float32` or quantized `uint8/uint16`
- quantized variants use the same missing-marker rules as continuous obs values

### `connectivity/edges.src.bin` and `edges.dst.bin`

Two parallel integer arrays of equal length:
- dtype: `uint16` / `uint32` / `uint64` (see manifest)
- semantics: undirected edges stored once as `(src, dst)` with `src < dst`
- arrays are sorted by `(src, dst)` to improve gzip compression

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
- `var_gene_id_column`: `"index"` or a column name
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
- `index_dtype` (`uint16`/`uint32`/`uint64`) and `index_bytes`
- `sourcesPath`, `destinationsPath`
- `compression`: gzip level or null

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

# Find one categorical field entry
field = obs_manifest["_categoricalFields"][0]
key = field[0]
codes_dtype = field[2]  # "uint8" or "uint16"
codes_ext = "u8" if codes_dtype == "uint8" else "u16"

safe = key  # careful: actual filename uses safe-key transformation
path = export_dir / "obs" / f"{safe}.codes.{codes_ext}"
raw = gzip.open(str(path) + ".gz", "rb").read() if (path.with_suffix(path.suffix + ".gz")).exists() else path.read_bytes()
codes = np.frombuffer(raw, dtype=np.uint8 if codes_ext == "u8" else np.uint16)
```

(In real tooling, you should apply the same safe-key transformation as `prepare`.)

---

## Invariants (the “do not break” list)

1) `n_points` must match the number of rows in all per-cell binaries.
2) Category codes must be stable and match the `categories` list in the manifest.
3) Missing markers must be respected (especially for quantized values).
4) If you change any schema/filenames, coordinate with the web app.
5) `.gz` means “raw gzip bytes”, not HTTP `Content-Encoding`.

---

## Troubleshooting

### Symptom: “Viewer loads, but fields are missing”

Likely causes:
- `obs_manifest.json` references a file that isn’t present (partial export),
- safe-key mismatch (field name sanitization),
- compression mismatch (`.gz` expected but missing).

Fix:
- re-export with `force=True`, or delete and re-run `prepare`.

### Symptom: “Continuous field looks blank”

Likely causes:
- values are mostly NaN/Inf and got mapped to missing marker,
- min/max range is extreme (everything clips).

Confirm:
- inspect min/max in the manifest for that field,
- check the raw values upstream in Python.
