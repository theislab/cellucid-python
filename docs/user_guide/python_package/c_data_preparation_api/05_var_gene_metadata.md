# Var / gene metadata

**Audience:** everyone exporting gene expression (computational users will care most)  
**Time:** 15–30 minutes  
**Goal:** ensure gene identifiers are stable, unique, and match your expression matrix

`var` is the gene/feature metadata table associated with `gene_expression`.

In the current Python exporter, `var` is primarily used to determine the **gene
names** that:
- appear in the UI (gene search / gene overlay names),
- and are recorded in `var_manifest.json`.

A gene name is never a filename. Payload files under `var/` are named by integer
index, so `var_manifest.json` is the only place a gene's identity lives.

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

`prepare()` requires distinct gene names, and requires the **exported** ones to
read on screen as the value they store. Duplicates and names carrying
characters with no glyph fail before the export is published.

The two rules have different scopes, because they answer different questions:

| Rule | Scope | Why |
| --- | --- | --- |
| names are non-empty strings and **distinct** | every row of `var` | `gene_identifiers` addresses rows by name, so a repeated name names no single row |
| names are **drawable exactly as stored** | only the genes actually exported | the rule is about text the viewer shows, and a deselected gene is never shown |

`HLA-DRB1/2` is a real HGNC symbol and exports fine, exported or not: it is
never a path. `'MS4A1 '` with a trailing space does not, because the legend
would draw it identically to `'MS4A1'`. This is the same scoping `obs_keys`
already has (see {doc}`04_obs_cell_metadata`). `cellucid_prepare()` in the R
package applies the identical scopes.

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
draw: a name carrying invisible characters only has to be cleaned when the gene
is among the exported ones.

Reproducibility tip:
- store the exact gene list used for export in your pipeline (and ideally version-control it).

### Large `var` tables

The current exporter does **not** export arbitrary `var` columns into the viewer.

So, large additional columns in `adata.var` do not directly affect export size unless you use them as identifiers.
However:
- large object columns can still slow down your own pre-processing,
- and they can be a privacy risk if you accidentally use them as IDs.

---

## Naming rules

Exported gene names are never rewritten and never become paths, so there is no
character set, length limit, reserved name, or case-collision rule. `Gene` and
`gene` are two distinct genes and export as two payloads.

An exported name must be text the viewer can draw exactly as it is stored: no
control characters, none of the zero-width characters `U+200B`, `U+2060`,
`U+FEFF`, and no leading or trailing whitespace of any kind. A name that breaks
this rejects the complete candidate before publication and the message names it.

A name that stays in `var` but is left out by `gene_identifiers` is never
checked against the display rule, because it is never drawn. Distinctness across
all of `var` is still required either way.

---

## Troubleshooting (var / gene ids)

### Symptom: gene search returns nothing

Likely causes:
- You did not export `gene_expression` at all (no `var_manifest.json`).
- You exported a small `gene_identifiers` list and the queried gene isn’t included.
- Gene IDs in the UI are not what you expected (wrong `var_gene_id_column`).

How to confirm:
- Does `<out_dir>/var_manifest.json` exist?
- Open it and read `fields`: each entry is `[index, name]`, or
  `[index, name, minValue, maxValue]` when quantized.

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
