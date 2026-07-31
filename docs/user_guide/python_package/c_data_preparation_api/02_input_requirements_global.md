# Input requirements (global)

**Audience:** computational users (recommended for anyone exporting real data), plus developers implementing pipelines  
**Time:** 20–30 minutes  
**Goal:** prevent subtle exports that “succeed” but load incorrectly or behave strangely in the viewer

This page documents the **global rules** that apply to *every* input you pass to {func}`~cellucid.prepare`.

---

## Fast path (preflight checklist)

Before you export a real dataset, confirm:

1) You have **at least one embedding** (`X_umap_1d`, `X_umap_2d`, or
   `X_umap_3d`).
2) `latent_space` is provided (required by `prepare()`).
3) `obs` has exactly `n_cells` rows and matches the embedding row order.
4) If exporting genes:
   - `gene_expression.shape == (n_cells, n_genes)`
   - `var` has `n_genes` rows and is aligned to gene_expression columns
5) Embeddings contain **no NaN/Inf** (must be finite).
6) Exported identifiers are non-empty, distinct within their axis, and text a
   reader can see verbatim (see Rule 5).
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

`prepare()` validates embeddings before publication. `NaN`, infinities,
complex values, constant coordinate domains, and values outside finite
`float32` range reject the complete candidate.

### What to do instead

- Remove cells with invalid coordinates (most common).
- Recompute the embedding.
- Avoid sentinel values like `1e20` for missing coordinates; they break
  normalization.

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

## Rule 5: identifiers are names, not filenames

Every payload file an export writes is named by an integer index, so an
identifier is never a path. Use your real names: `HLA-DRB1/2`, `% mito`,
`细胞`, and the two distinct columns `Field` and `field` all export unchanged.
There is no ASCII restriction, no length limit, no reserved-name list, and no
case-insensitive collision rule.

`prepare()` never rewrites an identifier, so three rules remain. Observation
keys, gene names, and vector-field IDs must each:

- be a non-empty string,
- be distinct within their axis, and
- read on screen as the value they store — no control characters, none of the
  zero-width characters `U+200B`, `U+2060`, `U+FEFF`, and no leading or
  trailing whitespace of any kind, `U+00A0` included.

The third rule is the same one string category labels obey, for the same
reason: a gene name, an obs key, and a category label are all drawn verbatim in
the legend and the field selector, so `'Liver '` is a different value from
`'Liver'` while looking identical on screen. Nothing is trimmed for you —
trimming would rewrite an annotation you did not ask to change.

The display rule applies to the identifiers you **export**. `obs_keys` and
`gene_identifiers` narrow the export, and they narrow this rule with it: a
column or gene you leave out reaches no manifest and is never drawn. Gene names
carry one extra rule that is *not* narrowed — every name in `var` must be a
non-empty string and must be distinct, because `gene_identifiers` addresses
`var` rows by name.

`dataset_id` is the one exception, because it is not a field identifier: it
names a real directory in a multi-dataset export root and in a served URL, so
it must still be 1–180 ASCII bytes, begin with an ASCII letter or digit,
otherwise use only ASCII letters, digits, `.`, `_`, or `-`, not end in `.`, and
not be a Windows device name.

---

## Rule 6: output directory semantics (`out_dir`, `force`)

`out_dir` defaults to `./exports` and is resolved against your current working
directory **at the moment you call `prepare()`**, so a later `os.chdir()` changes
where the next call writes. You should **always set it explicitly**.

Important behavior:
- A leading `~` is expanded to your home directory, so
  `out_dir="~/exports/my_dataset"` writes under `$HOME`.
- A non-empty target raises `FileExistsError` when `force=False`.
- `force=True` writes one complete sibling stage, validates it, and atomically
  replaces the previous generation.
- A rejected or interrupted candidate leaves the previous generation unchanged.

Practical rules:
- During iteration: use `force=True` for an intentional complete replacement or
  a fresh `out_dir` per run.
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

It also validates exact identifier portability and collisions, category
capacity, finite scientific arrays, vector domains, and atomic output
publication. The preflight remains useful because it catches data-preparation
mistakes earlier and lets you attach study-specific context to the error.

---

## Preflight checks (copy/paste)

(python-preflight-ann-data)=
### AnnData-based preflight

```python
import numpy as np

embedding_key = next(
    (key for key in ("X_umap_3d", "X_umap_2d", "X_umap_1d") if key in adata.obsm),
    None,
)
if embedding_key is None:
    raise ValueError("Provide X_umap_1d, X_umap_2d, or X_umap_3d in adata.obsm.")

X_umap = np.asarray(adata.obsm[embedding_key])
expected_dimension = int(embedding_key[-2])
if X_umap.ndim != 2 or X_umap.shape[1] != expected_dimension:
    raise ValueError(
        f"{embedding_key} must have shape (n_cells, {expected_dimension}); "
        f"got {X_umap.shape}."
    )

n_cells = X_umap.shape[0]
if adata.n_obs != n_cells:
    raise ValueError("AnnData n_obs does not match embedding row count (alignment bug).")

if not np.isfinite(X_umap).all():
    bad = np.where(~np.isfinite(X_umap).all(axis=1))[0][:10]
    raise ValueError(f"Embedding contains NaN/Inf; fix before export. Example bad rows: {bad.tolist()}")

latent = adata.obsm.get("X_pca", X_umap)
if latent is None or latent.shape[0] != n_cells:
    raise ValueError("latent_space required and must match n_cells. Use adata.obsm['X_pca'] or another aligned latent representation.")

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
