# Var / gene metadata

**Audience:** everyone exporting gene expression (computational users will care most)  
**Time:** 15–30 minutes  
**Goal:** ensure gene identifiers are stable, unique, and match your expression matrix

`var` is the gene/feature metadata table associated with `gene_expression`.

In the current Python exporter, `var` is primarily used to determine the **gene identifiers** that:
- appear in the UI (gene search / gene overlay names),
- and become the “keys” in `var_manifest.json`.

```{important}
The exporter assumes `var` row order matches `gene_expression` column order.
If you reorder one without the other, the viewer will show the wrong gene values under the wrong names.
```

---

## Fast path (minimum viable gene metadata)

1) Choose which identifier you want the UI to use:
   - **Gene symbols** (human-readable, but can be ambiguous/duplicated), or
   - **Stable IDs** (e.g., Ensembl IDs; better for reproducibility).
2) Ensure it is:
   - **present for every gene**, and
   - **unique** (no duplicates).
3) Export with:
   - `var_gene_id_column=None` if `var.index` is the identifier you want, or
   - `var_gene_id_column="<column_name>"` if the identifier lives in a column.

---

## Practical path (computational users)

### Required alignment with gene expression

If `gene_expression.shape == (n_cells, n_genes)`, then:
- `len(var)` must equal `n_genes`,
- and `var.iloc[j]` must describe `gene_expression[:, j]`.

AnnData makes this easy because `adata.var` is aligned to `adata.X` by construction,
but alignment bugs often happen after manual filtering/reindexing.

### Choosing gene identifiers (`var_gene_id_column`)

`prepare()` chooses gene IDs as follows:

- If `var_gene_id_column is None` (default), gene IDs come from `var.index`.
- Every string is an exact column selector, including the literal `"index"`.
- Identifiers must already be native non-empty strings.

Recommendation:
- For reproducible exports intended for sharing, prefer stable identifiers.
- For wet lab-facing demos, gene symbols may be friendlier (if you can guarantee uniqueness).

#### Example: use gene symbols from a column

```text
prepare(
    ...,
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="gene_symbol",
    ...
)
```

### Uniqueness (do not skip this)

`prepare()` requires unique gene IDs, and requires the **exported** ones to be
exact portable filename components. Duplicates, unsafe IDs, and case-insensitive
filesystem collisions fail before the export is published.

The two rules have different scopes, because they answer different questions:

| Rule | Scope | Why |
| --- | --- | --- |
| IDs are non-empty strings and **distinct** | every row of `var` | `gene_identifiers` addresses rows by ID, so a repeated ID names no single row |
| IDs are **portable filename components**, unique case-insensitively | only the genes actually exported | the rule is about the file the gene is written to, and a deselected gene is written to none |

So a `var` that carries `HLA-DRB1/2` still exports cleanly as long as
`gene_identifiers` leaves that gene out. This is the same scoping `obs_keys`
already has (see {doc}`04_obs_cell_metadata`): narrowing the export narrows the
filename contract with it. `cellucid_prepare()` in the R package applies the
identical scopes.

Preflight check:

```python
import pandas as pd

ids = adata.var.index.astype(str)  # or adata.var["gene_symbol"].astype(str)
dupes = pd.Index(ids).duplicated(keep=False)
if dupes.any():
    raise ValueError(f"Duplicate gene IDs detected. Examples: {sorted(set(ids[dupes]))[:10]}")
```

### Subsetting genes (`gene_identifiers`)

Exporting all genes can be huge (see {doc}`06_gene_expression_matrix`).

Use `gene_identifiers` to export a curated list:

```text
marker_genes = ["MS4A1", "CD3D", "LYZ", "NKG7"]

prepare(
    ...,
    gene_expression=adata.X,
    var=adata.var,
    var_gene_id_column=None,
    gene_identifiers=marker_genes,
    ...
)
```

If any requested gene is absent, the exporter raises `KeyError` and publishes
nothing. A repeated request raises `ValueError` naming the repeated identifier.

`gene_identifiers` is also the supported way past a gene the exporter cannot
write: an ID that is not a portable filename component, or that collides
case-insensitively with another ID, only has to be resolved when the gene is
among the exported ones.

Reproducibility tip:
- store the exact gene list used for export in your pipeline (and ideally version-control it).

### Large `var` tables

The current exporter does **not** export arbitrary `var` columns into the viewer.

So, large additional columns in `adata.var` do not directly affect export size unless you use them as identifiers.
However:
- large object columns can still slow down your own pre-processing,
- and they can be a privacy risk if you accidentally use them as IDs.

---

## Naming rules (exact portable IDs)

Exported gene IDs are not rewritten. They must be 1–180-byte ASCII components
that begin with a letter or digit, contain only letters, digits, `.`, `_`, or
`-`, do not end with `.`, are not dot segments or Windows device names, and are
unique under case-insensitive comparison among the exported IDs. IDs containing
spaces, slashes, or other unsafe characters are rejected before publication.
Collisions are possible (rare for Ensembl IDs, more common for messy symbols).

An ID that stays in `var` but is left out by `gene_identifiers` is never checked
against these rules, because it is never written to a path.

---

## Troubleshooting (var / gene ids)

### Symptom: gene search returns nothing

Likely causes:
- You did not export `gene_expression` at all (no `var_manifest.json`).
- You exported a small `gene_identifiers` list and the queried gene isn’t included.
- Gene IDs in the UI are not what you expected (wrong `var_gene_id_column`).

How to confirm:
- Does `<out_dir>/var_manifest.json` exist?
- Open it and see which gene IDs are present in `fields`.

Fix:
- Export gene expression, or export the genes you need, and re-export with `force=True`.

### Symptom: the “wrong gene” appears (values don’t match expectations)

Likely cause:
- `var` row order does not match `gene_expression` column order.

Fix:
- Rebuild `var` and `gene_expression` together (AnnData usually prevents this bug if used correctly).

---

## Next steps

- Expression export details (quantization, size, sparse vs dense): {doc}`06_gene_expression_matrix`
