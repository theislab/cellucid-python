# Export format spec and invariants

The export folder format written by `cellucid.prepare(...)` is a contract
between five repositories: two writers (`cellucid-python`, `cellucid-r`), one
reader (the Cellucid web app), and two dataset repositories that ship artifacts
in it.

:::{important}
**The specification is
{doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`.**
It is the single normative description of every file, key, dtype, and byte
order, and it is the page a sixth writer should be implemented against. This
page is the *developer* view of the same format: what the layout looks like at a
glance, what a reader may rely on, and what changing it costs. Where the two
disagree, page 09 is right and this page is a bug.
:::

Related pages:
- Export pipeline internals: {doc}`07_prepare_export_pipeline_architecture`
- Server behavior: {doc}`09_server_mode_architecture_endpoints_and_security`

---

## What counts as a valid export folder

There is one level, not two. `_list_exported_datasets` admits a directory only
when **all** of the following hold, and rejects the whole served root otherwise:

- `dataset_identity.json` is present, readable UTF-8 JSON, an object, and
  carries `"version": 2`;
- its `id` and `name` pass the same contracts the writer applied;
- `obs_manifest.json` is present, readable UTF-8 JSON, and an object;
- at least one `points_<dim>d.bin` or `points_<dim>d.bin.gz` is present and
  non-empty;
- no dimension has *both* a compressed and an uncompressed points file.

`var_manifest.json`, `connectivity_manifest.json`, and `vectors/` are genuinely
optional; `dataset_identity.json` is not. A hand-assembled export that omits it,
or that carries a different `version`, is refused before the server binds — not
loaded with reduced functionality.

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
- every multi-byte value is **little-endian**, pinned by the writer rather than
  inherited from the exporting machine, so an export made on a big-endian host
  is byte-identical to one made on x86 or ARM
- arrays are **row-major** (C order)
- the dtype names below carry no byte-order component, and the viewer builds its
  typed arrays on the received bytes, so a mismatch would decode as wrong values
  rather than as a load failure — see
  {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`

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

The exact key set of each manifest — `dataset_identity.json`,
`obs_manifest.json`, `var_manifest.json`, and `connectivity_manifest.json` — is
specified once, in
{doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`.
Every one of them is validated on both sides: the writer reconciles the manifest
it just built against the directory it just wrote, and the reader refuses a
manifest whose key set is not exactly the current one. Copying the key lists
here is how they drifted before, so this page states only the shape a reader has
to parse, and links for the rest.

**Compact entries are length-typed.** Each manifest entry is a positional array
whose *length* selects the variant, so an entry with an unexpected member count
is rejected rather than read partially:

| Manifest array | Members | Shape |
| --- | ---: | --- |
| `_continuousFields` | 2 | `[index, key]` — float32, unquantized |
| `_continuousFields` | 4 | `[index, key, minValue, maxValue]` — quantized |
| `_categoricalFields` | 6 | `[index, key, categories, codesDtype, codesMissingValue, centroidsByDim]` — outliers float32 |
| `_categoricalFields` | 8 | the six above plus `outlierMinValue, outlierMaxValue` — outliers quantized |
| `fields` (var) | 2 | `[index, name]` — float32, unquantized |
| `fields` (var) | 4 | `[index, name, minValue, maxValue]` — quantized |

The leading member of every entry is the field's payload index, and it is what
the schema `pathPattern` is expanded with. It is not the position of the entry
in the array: `obs` shares one index space across `_continuousFields` and
`_categoricalFields` because both write into `obs/`, so neither list is dense on
its own.

`centroidsByDim` is an object keyed by dimension as a string, whose values are
lists of `{category, position, n_points}` objects.

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
4) `.gz` means “raw gzip bytes”, not HTTP `Content-Encoding`.
5) A schema, key set, filename, dtype, or byte order changed on one side has to
   change on every side in the same pass: both writers (`cellucid-python`,
   `cellucid-r`), the web app's reader, the prepared demo datasets, and this
   documentation. The specification page is
   {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
   — update it first, because it is the page the other implementations are read
   against.

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
