# `prepare()` / export pipeline architecture

This page explains how `cellucid.prepare(...)` works internally: what it writes, why each step exists, and what can break if you change it.

If you only need the **user-facing** “how do I export my AnnData?” guide, start with:
- {doc}`../c_data_preparation_api/index`

:::{admonition} Audience
:class: note

- Beginners/wet lab: you probably don’t need this page unless you hit an export failure; jump to *Troubleshooting*.
- Computational users: focus on “inputs, outputs, and invariants”.
- Developers: read the whole page before changing the export format.
:::

---

## What `prepare(...)` produces (conceptually)

`prepare(...)` turns “Python-side data structures” into a **static export folder** that the web app can load:

- embedding coordinates (`points_2d.bin`, `points_3d.bin`, …)
- obs fields (categorical + continuous)
- gene expression (var fields)
- optional connectivity (KNN edges)
- optional vector fields (velocity/drift overlays)
- metadata (`dataset_identity.json`)

The web app expects this contract (documented fully here):
{doc}`08_export_format_spec_and_invariants`.

---

## High-level pipeline stages (the mental model)

Implementation is in:
- `cellucid-python/src/cellucid/prepare_data.py` (`def prepare(...)`)

The pipeline stages are:

1) **Validate inputs** (shapes, presence of required objects, dtype assumptions)
2) **Normalize embeddings** (each dimension independently to `[-1, 1]`, aspect ratio preserved)
3) **Write embedding binaries** (`points_<dim>d.bin(.gz)`)
4) **Write vector fields** (optional) (`vectors/<field>_<dim>d.bin(.gz)`)
5) **Write obs field binaries + `obs_manifest.json`**
6) **Write gene binaries + `var_manifest.json`**
7) **Write connectivity + `connectivity_manifest.json`** (optional)
8) **Write `dataset_identity.json`** (always)

---

## Stage 1 — Input validation (what’s required and why)

### Required inputs

`prepare(...)` currently requires:

- `latent_space`: used to compute outlier quantiles per category (helps visualization and centroid stability)
- `obs`: the per-cell metadata table
- at least one embedding: `X_umap_1d`/`2d`/`3d`/`4d`

Optional inputs:
- `var` + `gene_expression` for gene overlays
- `connectivities` for KNN edge overlay
- `vector_fields` for velocity/drift overlays

### The key invariant: `n_cells` must match everywhere

Everything is indexed by **cell row order**:

- embeddings have shape `(n_cells, dim)`
- `latent_space` has shape `(n_cells, latent_dim)`
- `obs` has `len(obs) == n_cells`
- `gene_expression` has shape `(n_cells, n_genes)`
- `connectivities` has shape `(n_cells, n_cells)`
- every vector field has shape `(n_cells, dim)`

If these are inconsistent, `prepare` raises early (shape mismatch).

---

## Stage 2 — Embedding normalization

For each dimensional embedding (1D/2D/3D/4D), `prepare` normalizes coordinates to fit within `[-1, 1]`.

Important details:

- **Normalization is per dimension**: 2D and 3D embeddings are normalized independently.
- **Aspect ratio is preserved** within a dimension: a single scale factor is used for all axes.
- The viewer depends on this so “zoom level” is comparable across datasets and embedding dims.

Developer impact:
- if you change normalization, you also change the meaning/scale of vector fields and centroids.

---

## Stage 3 — Writing `points_<dim>d.bin(.gz)`

For each available dimension:
- file name: `points_<dim>d.bin` (or `points_<dim>d.bin.gz` if gzip enabled)
- dtype: `float32`
- layout: row-major, packed `(n_cells, dim)` floats

If you enable gzip (`compression` parameter):
- the file is written as a **raw gzip file** and gets a `.gz` suffix
- the manifest/identity files point to the `.gz` file path

Do not confuse this with HTTP `Content-Encoding: gzip` (that’s a separate mechanism used in AnnData server mode).

---

## Stage 4 — Vector fields (optional)

Vector fields are per-cell displacement vectors (e.g., velocity, CellRank drift) aligned to the embedding space.

Inputs:
- `vector_fields: dict[str, array]`

Accepted shapes:
- `(n_cells,)` or `(n_cells, 1)` for 1D
- `(n_cells, 2)` for 2D
- `(n_cells, 3)` for 3D

Naming conventions:

- Explicit dimension keys: `<field_id>_<dim>d` (e.g. `velocity_umap_2d`, `T_fwd_umap_3d`)
- Implicit keys without suffix are allowed **only if their shape matches a missing dimension**

Export rules:
- vector files are written under `vectors/`
- file name: `vectors/<field_id>_<dim>d.bin(.gz)`
- dtype: `float32`
- scaling: vectors are scaled by the **same scale factor used to normalize the corresponding points**

Important invariant:
- `field_id` must be “filename-safe” (only `[A-Za-z0-9._-]` and underscores). Unsafe IDs raise.

---

## Stage 5 — Obs fields + `obs_manifest.json`

Obs export does three things:

1) decides whether each obs column is **continuous** or **categorical**,
2) writes one or more binaries per field,
3) writes `obs_manifest.json` (a compact schema that tells the viewer what exists).

### 5.1 Field typing rules (current)

`prepare` classifies each obs column as:

- **categorical** if:
  - pandas categorical dtype, or
  - bool dtype, or
  - non-numeric dtype (strings, objects)
- **continuous** if numeric dtype

### 5.2 Continuous obs export

Files:
- `obs/<safe_key>.values.f32(.gz)` when not quantized
- `obs/<safe_key>.values.u8(.gz)` or `.values.u16(.gz)` when quantized

If quantized:
- NaN/Inf values are mapped to a reserved “missing marker” (255 or 65535)
- `obs_manifest.json` stores `minValue` and `maxValue` for dequantization

### 5.3 Categorical obs export

Files:
- `obs/<safe_key>.codes.u8(.gz)` or `obs/<safe_key>.codes.u16(.gz)`
- `obs/<safe_key>.outliers.f32(.gz)` (or `.outliers.u8/.u16` if outlier quantization enabled)

Codes:
- integer codes `0..(n_categories-1)`
- missing values use a reserved marker (255 or 65535)

Centroids:
- computed per category per available embedding dimension
- stored in `obs_manifest.json` as a dict keyed by dimension (`"1"`, `"2"`, `"3"`, …)

Outlier quantiles:
- computed per cell in **latent space** (not embedding space)
- used to support outlier-aware centroid display and UI behavior

### 5.4 The obs manifest “compact format”

`obs_manifest.json` uses a compact schema for size and fast parsing:

- `_continuousFields`: list of `[key]` or `[key, minValue, maxValue]`
- `_categoricalFields`: list containing categories, dtype, missing marker, centroids, and optional outlier min/max
- `_obsSchemas`: path patterns that tell the viewer how to derive filenames

Full spec lives in: {doc}`08_export_format_spec_and_invariants`.

---

## Stage 6 — Gene expression + `var_manifest.json` (optional)

Gene export is the most time-consuming step for large datasets because it writes one file per gene.

Inputs:
- `gene_expression`: dense or sparse matrix, shape `(n_cells, n_genes)`
- `var`: gene metadata table, length `n_genes`

Gene identifiers:
- chosen via `var_gene_id_column` (`"index"` by default)
- optionally subset via `gene_identifiers`

Export:
- one file per gene under `var/`
- file: `var/<safe_gene_id>.values.f32(.gz)` or quantized `.values.u8/.u16(.gz)`
- `var_manifest.json` stores either:
  - `[gene_id]` (not quantized), or
  - `[gene_id, minValue, maxValue]` (quantized)

Important performance note:
- if `gene_expression` is CSR, `prepare` converts to CSC for efficient column access.

---

## Stage 7 — Connectivity export (optional)

Connectivity export turns a KNN matrix into a GPU-friendly edge list:

1) ensure CSR
2) symmetrize: `A + A.T`
3) binarize: set all nonzero entries to 1
4) extract unique undirected edges as `(src, dst)` with `src < dst`
5) sort edges lexicographically for better gzip compression
6) write two binary arrays:
   - `connectivity/edges.src.bin(.gz)`
   - `connectivity/edges.dst.bin(.gz)`

Index dtype is chosen based on `n_cells`:
- `uint16` if `n_cells <= 65535`
- `uint32` if `n_cells <= 4294967295`
- else `uint64`

The manifest stores:
- `n_edges`, `max_neighbors`, `index_dtype`, and file paths

---

## Stage 8 — `dataset_identity.json` (always written)

`dataset_identity.json` is the “front door” metadata:

- dataset id/name/description
- creation timestamp
- stats: cells/genes/fields, connectivity presence
- embeddings: available dims and point file paths
- obs field summary (key + kind + n_categories)
- export settings (compression, quantization choices)
- optional `source` metadata (name/url/citation)
- optional vector field metadata

This is how:
- the CLI detects an export folder (`dataset_identity.json` presence),
- servers list datasets and show names/stats,
- and the UI can display meaningful metadata without scanning all files.

---

## Debugging exports (developer checklist)

When you change `prepare`, verify:

1) `dataset_identity.json` exists and looks sane.
2) `obs_manifest.json` parses and references real files.
3) at least one `points_*d.bin` exists (and has the right shape).
4) optional files (vectors/connectivity/var) are internally consistent with manifests.

Quick Python snippet to validate points shape:

```python
import json
import gzip
import numpy as np
from pathlib import Path

export_dir = Path("my_export")
identity = json.loads((export_dir / "dataset_identity.json").read_text("utf-8"))
n = identity["stats"]["n_cells"]
points_path = export_dir / identity["embeddings"]["files"]["3d"]
raw = gzip.open(points_path, "rb").read() if points_path.suffix == ".gz" else points_path.read_bytes()
pts = np.frombuffer(raw, dtype=np.float32).reshape(n, 3)
print(pts.shape, pts.min(), pts.max())
```

---

## Edge cases and footguns (common)

### Partial exports due to “skip existing files”

By default `prepare` **skips writing** files that already exist unless `force=True`.

This can produce a mixed export folder if you:
- change quantization/compression settings,
- re-run prepare without `force`,
- and end up with manifests referring to files that were not regenerated.

Recommendation:
- during development, either delete the output folder or use `force=True`.

### Too many categories for `uint8`

If you force categorical dtype `uint8` and a field has >254 categories, export will error.

Use:
- `obs_categorical_dtype="auto"` (preferred) or `"uint16"`.

### NaN/Inf handling

Quantized continuous data maps NaN/Inf to a missing marker (255/65535).
If you have many NaNs, the color map may look “empty” even though export succeeded.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “`ValueError: latent_space is required`”

Cause:
- `prepare` uses latent space to compute outlier quantiles.

Fix:
- pass a reasonable latent space (e.g. PCA, scVI):
  - `latent_space=adata.obsm["X_pca"]`

### Symptom: “`ValueError: obs has N rows, but embeddings have M cells`”

Cause:
- you subset/reordered one array but not the others.

Fix:
- ensure all inputs are aligned to the same row order (same cell ordering).

### Symptom: “Export is extremely slow / huge”

Common causes:
- writing one file per gene without quantization or compression,
- dense gene expression matrix,
- exporting too many genes.

Fix options:
- set `var_quantization=8` and `compression=6`,
- export fewer genes (use `gene_identifiers`),
- keep gene expression sparse if possible.

### Symptom: “Vector field export errors about unsupported characters”

Cause:
- vector field ids must be filename-safe (no spaces, slashes, etc.).

Fix:
- rename keys to safe ids (e.g. `T_fwd_umap` instead of `T fwd/umap`).
