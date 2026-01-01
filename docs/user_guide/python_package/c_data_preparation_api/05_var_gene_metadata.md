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
   - `var_gene_id_column="index"` if `var.index` is the identifier you want, or
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

- If `var_gene_id_column == "index"` (default), gene IDs come from `var.index`.
- Otherwise, gene IDs come from `var[var_gene_id_column]` (cast to string).

Recommendation:
- For reproducible exports intended for sharing, prefer stable identifiers.
- For wet lab-facing demos, gene symbols may be friendlier (if you can guarantee uniqueness).

#### Example: use gene symbols from a column

```python
prepare(
    ...,
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="gene_symbol",
    ...
)
```

### Uniqueness (do not skip this)

`prepare()` does not currently enforce uniqueness of gene IDs.
Duplicates can cause:
- overwritten files on disk (after safe filename sanitization),
- confusing behavior where multiple manifest entries map to the same gene payload,
- “wrong gene” values in the UI.

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

```python
marker_genes = ["MS4A1", "CD3D", "LYZ", "NKG7"]

prepare(
    ...,
    gene_expression=adata.X,
    var=adata.var,
    var_gene_id_column="index",
    gene_identifiers=marker_genes,
    ...
)
```

If you request genes that are not found, the exporter prints a warning and skips them.

Reproducibility tip:
- store the exact gene list used for export in your pipeline (and ideally version-control it).

### Large `var` tables

The current exporter does **not** export arbitrary `var` columns into the viewer.

So, large additional columns in `adata.var` do not directly affect export size unless you use them as identifiers.
However:
- large object columns can still slow down your own pre-processing,
- and they can be a privacy risk if you accidentally use them as IDs.

---

## Naming rules (safe filenames)

Gene IDs are converted to safe filenames using the same sanitization rules as obs keys:

```text
safe = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_id)
safe = safe.strip("._")
safe = safe or "field"
```

This means gene IDs containing slashes/spaces/etc. will be rewritten on disk.
Collisions are possible (rare for Ensembl IDs, more common for messy symbols).

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
