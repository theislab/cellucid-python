# Input requirements (global)

**Audience:** computational users (recommended for anyone exporting real data), plus developers implementing pipelines  
**Time:** 20–30 minutes  
**Goal:** prevent subtle exports that “succeed” but load incorrectly or behave strangely in the viewer

This page documents the **global rules** that apply to *every* input you pass to {func}`~cellucid.prepare`.

---

## Fast path (preflight checklist)

Before you export a real dataset, confirm:

1) You pass **at least one embedding** in `embeddings={...}`, and name the
   default with `embedding_defaults` for every dimension that has more than one.
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
- every embedding in `embeddings`
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
- and at least one embedding in `embeddings`.

### Embeddings (at least one required)

`embeddings` maps an embedding id to an array or to an
{class}`~cellucid.Embedding`:

| Value | Required shape | Dimension |
|---|---:|---|
| a plain array | `(n_cells, 1)`, `(n_cells, 2)`, or `(n_cells, 3)` | its column count |
| a plain vector | `(n_cells,)` | 1D |
| `Embedding(values, columns=[...])` | `(n_cells, width)`, any width | the number of `columns` (1–3) |

Any number of embeddings may share a dimension; `embedding_defaults` then names
the one shown first, for example `{"2d": "X_umap"}`. `prepare()` reads no `obsm`
key: the id is a name you choose, and `adata.obsm["X_umap"]` is simply the array
you pass under it. Details: {doc}`03_embeddings_and_coordinates`.

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
| `vector_fields[field][embedding_id]` | `(n_cells, dim)` of that embedding | Optional velocity/displacement vectors, bound to an exported embedding |
| `var_gene_name_columns` | list of `var` columns | Optional alternative gene-name sets |

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

`prepare()` never rewrites an observation key, embedding id, or vector-field
ID, so three rules remain. Each must:

- be a non-empty string,
- be distinct within their axis, and
- read on screen as the value they store — no control characters, none of the
  zero-width characters `U+200B`, `U+2060`, `U+FEFF`, no unpaired surrogate,
  and no leading or trailing whitespace of any kind, `U+00A0` included.

The third rule is the same one string category labels obey, for the same
reason: a gene name, an obs key, and a category label are all drawn verbatim in
the legend and the field selector, so `'Liver '` is a different value from
`'Liver'` while looking identical on screen. Nothing is trimmed for you —
trimming would rewrite an annotation you did not ask to change.

The display rule applies to the identifiers you **export**. `obs_keys` and
`gene_identifiers` narrow the export, and they narrow this rule with it: a
column or gene you leave out reaches no manifest and is never drawn.

Gene names are the one identifier Cellucid completes for you, by one rule shared
with the R writer, the AnnData server, and the viewer: an empty value (`""`,
`None`, NaN) takes the row's `var.index` name, and a repeated name is made
unique the way anndata's `var_names_make_unique()` does it (`MATR3`, `MATR3-1`).
The rule runs over the whole `var`, so `gene_identifiers` addresses rows by the
unique names, and the export prints and records what it changed. The resulting
names must still be drawable exactly as stored. See
{ref}`python_package-make-unique-gene-names`.

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

`out_dir` must name a directory **dedicated to this one export**, because
publishing replaces the whole directory and `force=True` removes everything the
previous one held. These are refused with a `ValueError` before any path is
created, locked, written, or removed:

- the filesystem root (`/`, or a drive root on Windows),
- your current working directory (`.`, `./`, or its absolute path),
- your home directory (`~` or its absolute path),
- the directory that holds every home (`..`, `/Users`, `/home`).

A symbolic link is resolved first, so an alias to one of those directories is
refused too. Name a child directory of your own instead:
`prepare(..., out_dir="./exports/my_dataset")` is correct, while
`prepare(..., out_dir=".", force=True)` raises before anything is touched.

The R package's `cellucid_prepare()` refuses the same set, so both writers agree
on which directories are never an export target.

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

embedding_key = "X_umap"  # the obsm array you will pass in embeddings={...}
if embedding_key not in adata.obsm:
    raise ValueError(f"adata.obsm has no {embedding_key!r}; compute an embedding first.")

X_umap = np.asarray(adata.obsm[embedding_key])
if X_umap.ndim != 2 or X_umap.shape[1] not in (1, 2, 3):
    raise ValueError(
        f"{embedding_key} has shape {X_umap.shape}; pass Embedding(values, columns=[...]) "
        "to plot chosen columns of a wider array."
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
assert isinstance(embeddings, dict) and embeddings

X = np.asarray(next(iter(embeddings.values())))
assert X.ndim == 2 and X.shape[1] in (1, 2, 3)

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

- Embedding details (normalization, several embeddings per dimension): {doc}`03_embeddings_and_coordinates`
