# Input requirements (global)

**Audience:** computational users (recommended for anyone exporting real data), plus developers implementing pipelines  
**Time:** 20–30 minutes  
**Goal:** prevent subtle exports that “succeed” but load incorrectly or behave strangely in the viewer

This page documents the **global rules** that apply to *every* input you pass to {func}`~cellucid.prepare`.

---

## Fast path (preflight checklist)

Before you export a real dataset, confirm:

1) You have **at least one embedding** (`X_umap_2d` or `X_umap_3d`).
2) `latent_space` is provided (required by `prepare()`).
3) `obs` has exactly `n_cells` rows and matches the embedding row order.
4) If exporting genes:
   - `gene_expression.shape == (n_cells, n_genes)`
   - `var` has `n_genes` rows and is aligned to gene_expression columns
5) Embeddings contain **no NaN/Inf** (must be finite).
6) Column names you care about are **safe and unique after sanitization** (see Rule 5).
7) You are not accidentally reusing an old `out_dir` with `force=False`.

If you want a copy/paste preflight script, jump to:
- {ref}`python-preflight-ann-data`
- {ref}`python-preflight-raw-arrays`

---

## Rule 1: cell identity is the row order

Cellucid’s export format does not store a separate “cell ID table”.

> Cell `i` is the `i`-th row in every exported array.

This is the single most important invariant.

You must keep a consistent row order across:
- embeddings (`X_umap_1d/2d/3d`)
- `latent_space`
- `obs`
- `gene_expression` rows (if provided)
- `connectivities` rows/cols (if provided)
- every vector field (if provided)

### Practical strategy (recommended)

If you have an AnnData, use it as the canonical source of alignment:
- subset/filter on the AnnData (`adata = adata[mask].copy()`),
- then pull arrays from `adata.obs`, `adata.obsm`, `adata.obsp`, `adata.X`.

If you construct arrays separately (common in pipelines), you must enforce alignment explicitly:
- choose a canonical cell ordering (e.g., a list of cell IDs),
- reorder every array/DataFrame to that ordering **before export**.

---

## Rule 2: required inputs and shapes

`prepare()` requires:
- `latent_space` (for categorical outlier quantiles),
- `obs`,
- and at least one embedding (`X_umap_1d` or `X_umap_2d` or `X_umap_3d`).

### Embeddings (at least one required)

| Argument | Required shape | Notes |
|---|---:|---|
| `X_umap_1d` | `(n_cells, 1)` | Optional |
| `X_umap_2d` | `(n_cells, 2)` | Common |
| `X_umap_3d` | `(n_cells, 3)` | Recommended if you want true 3D navigation |

```{warning}
`prepare()` rejects 1D embeddings shaped `(n_cells,)` and will raise a shape error.
Always pass a 2D array with an explicit second dimension.
```

### Required non-embedding inputs

| Argument | Required shape | Why it exists |
|---|---:|---|
| `latent_space` | `(n_cells, n_latent_dims)` | Used to compute per-cell outlier quantiles for categorical fields |
| `obs` | `n_cells` rows | Cell metadata (coloring/filtering) |

### Optional inputs

| Argument | Required shape | Notes |
|---|---:|---|
| `gene_expression` | `(n_cells, n_genes)` | Dense or sparse; written one file per gene |
| `var` | `n_genes` rows | Required if `gene_expression` is provided |
| `connectivities` | `(n_cells, n_cells)` | Typically `adata.obsp["connectivities"]` |
| `vector_fields[field]` | `(n_cells, dim)` | Optional velocity/displacement vectors aligned to embedding space |

```{warning}
Cellucid export expects **cells × genes** for `gene_expression`.
If you pass `genes × cells`, export will fail (shape mismatch) or silently export nonsense if you happen to also swap `var`.
```

---

## Rule 3: embeddings must be finite (no NaN/Inf)

`prepare()` normalizes each embedding by computing min/max ranges. If your embedding contains NaN/Inf:
- normalization produces NaNs,
- those NaNs get written to `points_*d.bin(.gz)`,
- and the web app may render nothing or crash.

### What to do instead

- Remove cells with invalid coordinates (most common).
- Recompute the embedding.
- Avoid “placeholder” values like `1e20` for missing coords (they break normalization).

---

## Rule 4: `obs` column types decide how fields behave

`prepare()` classifies each `obs` column as:

- **categorical** if:
  - pandas dtype is categorical, or
  - dtype is bool, or
  - dtype is “other” (strings, objects, datetimes, etc.)
- **continuous** if:
  - pandas dtype is numeric

This means a `datetime64` column becomes categorical unless you convert it.

Practical advice:
- Make “category” fields categorical explicitly (`astype("category")`) if you care about category order.
- Convert continuous-looking strings to numeric before export if you want a continuous color scale.

Details: {doc}`04_obs_cell_metadata`

---

## Rule 5: naming rules and filename collisions (do not ignore)

On disk, field keys are converted to filesystem-safe names using the same logic as the web app:

```text
safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
safe = safe.strip("._")
safe = safe or "field"
```

This is convenient (you can have spaces), but it can cause collisions:

| Original key | Safe filename key |
|---|---|
| `sample id` | `sample_id` |
| `sample/id` | `sample_id` |
| `sample_id` | `sample_id` |

If two different keys map to the same safe key, **files can overwrite each other**.

Recommendations:
- Use simple keys for exported fields (letters/numbers/underscore/dash).
- Ensure uniqueness *after* sanitization.
- Avoid leading/trailing dots/underscores (they are stripped).

This applies to:
- `obs` column names,
- gene identifiers (from `var_gene_id_column`),
- vector field ids (vector fields are stricter: unsafe ids raise an error).

---

## Rule 6: output directory semantics (`out_dir`, `force`)

`out_dir` defaults to `./exports/` under your current working directory, but you should **always set it explicitly**.

Important behavior:
- Many files are **skipped** if they already exist and `force=False`.
- `dataset_identity.json` is written at the end of every run (so it can become inconsistent with older manifests if you reused an `out_dir`).

Practical rules:
- During iteration: use `force=True` or a fresh `out_dir` per run.
- For publishable exports: use a clean new `out_dir` and keep it immutable.

---

## Optional: controlling what gets exported (`obs_keys`, `gene_identifiers`)

Subsetting is both a **performance tool** and a **reproducibility tool**:

- `obs_keys`: export only the metadata columns you want users to see.
- `gene_identifiers`: export only a curated gene set (markers/HVGs).

If you export everything by default, it’s easy to create:
- massive folders (especially genes),
- slow browser workflows (huge field lists),
- and confusing “where did this field come from?” situations.

---

## Deep path (what `prepare()` validates and how it fails)

`prepare()` validates (and raises errors for):
- missing required inputs (`latent_space`, `obs`, at least one embedding),
- wrong embedding dimensionality (must be exactly 1/2/3 columns),
- `n_cells` mismatch across inputs,
- missing `var` when `gene_expression` is provided,
- wrong connectivity shape (must be square),
- vector field shape mismatch and unsupported ids.

It does **not** currently validate:
- collisions after safe-key sanitization for obs/gene ids,
- pathological category counts (e.g., > 65534 categories),
- NaNs/Inf in embeddings (it will proceed and write broken outputs).

Treat the preflight checklist as mandatory.

---

## Preflight checks (copy/paste)

(python-preflight-ann-data)=
### AnnData-based preflight

```python
import numpy as np

X_umap = adata.obsm.get("X_umap_3d", adata.obsm.get("X_umap_2d", adata.obsm.get("X_umap")))
if X_umap is None:
    raise ValueError("Missing embedding: provide adata.obsm['X_umap_2d'], ['X_umap_3d'], or ['X_umap'].")

if X_umap.ndim != 2 or X_umap.shape[1] not in (2, 3):
    raise ValueError(f"Unexpected embedding shape: {X_umap.shape} (expected (n_cells, 2) or (n_cells, 3))")

n_cells = X_umap.shape[0]
if adata.n_obs != n_cells:
    raise ValueError("AnnData n_obs does not match embedding row count (alignment bug).")

if not np.isfinite(X_umap).all():
    bad = np.where(~np.isfinite(X_umap).any(axis=1))[0][:10]
    raise ValueError(f"Embedding contains NaN/Inf; fix before export. Example bad rows: {bad.tolist()}")

latent = adata.obsm.get("X_pca", X_umap)
if latent is None or latent.shape[0] != n_cells:
    raise ValueError("latent_space required and must match n_cells. Use adata.obsm['X_pca'] or a fallback latent.")

obs = adata.obs
if obs.shape[0] != n_cells:
    raise ValueError("obs must have n_cells rows.")

# Optional genes
if adata.X is not None:
    if adata.X.shape[0] != n_cells:
        raise ValueError("gene_expression must have n_cells rows.")
    if adata.var is None or adata.var.shape[0] != adata.X.shape[1]:
        raise ValueError("var must have n_genes rows matching gene_expression columns.")
```

(python-preflight-raw-arrays)=
### Raw arrays + DataFrames preflight

```python
import numpy as np
import pandas as pd
from scipy import sparse

def _is_array(x):
    return isinstance(x, (np.ndarray,)) or sparse.issparse(x)

assert _is_array(latent_space)
assert isinstance(obs, pd.DataFrame)
assert _is_array(X_umap_2d) or _is_array(X_umap_3d)

X = X_umap_3d if X_umap_3d is not None else X_umap_2d
X = np.asarray(X)
assert X.ndim == 2 and X.shape[1] in (2, 3)

n_cells = X.shape[0]
assert np.isfinite(X).all()
assert np.asarray(latent_space).shape[0] == n_cells
assert len(obs) == n_cells

if gene_expression is not None:
    assert _is_array(gene_expression)
    assert gene_expression.shape[0] == n_cells
    assert var is not None and len(var) == gene_expression.shape[1]

if connectivities is not None:
    assert connectivities.shape == (n_cells, n_cells)
```

---

## Troubleshooting pointers

- Shape mismatch errors / missing required inputs → {doc}`11_troubleshooting_prepare_export`
- “Viewer loads but everything looks wrong” (often row-order mismatch) → {doc}`11_troubleshooting_prepare_export`

## Next steps

- Embedding details (normalization, multi-dim exports): {doc}`03_embeddings_and_coordinates`
