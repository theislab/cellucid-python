# Var / gene metadata

**Audience:** everyone exporting gene expression (computational users will care most)  
**Time:** 15–30 minutes  
**Goal:** ensure gene names are stable, readable, and match your expression matrix

`var` is the gene/feature metadata table associated with `gene_expression`.

In the current Python exporter, `var` is used to determine the **gene names**
that:
- appear in the UI (gene search / gene overlay names),
- are recorded in `var_manifest.json`,
- and, through `var_gene_name_columns`, the alternative name sets the viewer
  can switch to (for example Ensembl IDs beside gene symbols).

A gene name is never a filename. Payload files under `var/` are named by integer
index, so `var_manifest.json` is the only place a gene's identity lives.

```{important}
The exporter assumes `var` row order matches `gene_expression` column order.
If you reorder one without the other, the viewer will show the wrong gene values under the wrong names.
```

---

## Fast path (minimum viable gene metadata)

1) Choose which names the UI shows first:
   - `var_gene_id_column=None` if `var.index` holds them, or
   - `var_gene_id_column="<column_name>"` if they live in a column.
2) Optionally offer more name sets with `var_gene_name_columns=[...]`, for
   example `["gene_ids"]` beside a symbol default.
3) Export. Empty and repeated names are completed for you by one shared rule,
   and the export prints what it changed.

---

## Practical path (computational users)

### Required alignment with gene expression

If `gene_expression.shape == (n_cells, n_genes)`, then:
- `len(var)` must equal `n_genes`,
- and `var.iloc[j]` must describe `gene_expression[:, j]`.

AnnData makes this easy because `adata.var` is aligned to `adata.X` by construction,
but alignment bugs often happen after manual filtering/reindexing.

### Choosing the default names (`var_gene_id_column`)

`prepare()` chooses the default gene names as follows:

- If `var_gene_id_column is None` (default), names come from `var.index`.
- Every string is an exact column selector, including the literal `"index"`.
- A column of text or categories names genes as it is, and a column of
  integers by their digits (`7157` is `"7157"`) — the rule the viewer,
  `cellucid serve` and `cellucid_prepare()` in R apply to the same column. A
  floating-point or boolean column cannot name genes and is refused with the
  reason (`var column 'score': floating-point values cannot name genes.`). An
  empty value (`""`, `None`, or NaN) takes the var index name of its row.

Recommendation:
- For reproducible exports intended for sharing, prefer stable identifiers as
  the default, and offer symbols as an alternative name set (or the reverse for
  wet lab-facing demos).

#### Example: symbols by default, Ensembl IDs as an alternative

```text
prepare(
    ...,
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="gene_symbols",
    var_gene_name_columns=["gene_ids"],
    ...
)
```

To offer the var index itself as an alternative while a column is the default,
copy it into a column first: `adata.var["ensembl_id"] = adata.var.index`.

### Alternative name sets (`var_gene_name_columns`)

Each column listed in `var_gene_name_columns` becomes one more set of names for
the same genes. The sets are metadata only: every set points at the same
expression payload of each gene, so adding a set duplicates no expression data
and the viewer can switch sets without reloading expression. Each column is
resolved with the same rule as the default (below). A column must exist in
`var`, may be listed once, and may not be the `var_gene_id_column` itself.

### Empty and repeated names: the make-unique rule

Gene symbols repeat (`MATR3` has two loci in some annotations) and symbol
columns have gaps, so every name set is completed by one rule, shared by
`cellucid.prepare()`, `cellucid_prepare()` in R, the Python AnnData server, and
the viewer — the rule anndata's `var_names_make_unique()` applies:

1. An empty value takes the var index name of the same row.
2. The first occurrence of a name keeps it.
3. Each later occurrence becomes `name-1`, `name-2`, …, skipping any value that
   already exists, so `A, A, A-1, A` becomes `A, A-2, A-1, A-3`.

The rule runs over the **whole** `var`, so a gene keeps its name however you
narrow the export, and `gene_identifiers` is matched against the unique names.
A set that was already unique is unchanged. The export prints what changed for
each set:

```text
Gene names from var column 'gene_symbols': 312 genes have no 'gene_symbols' value and use their var index name. 41 duplicate names were made unique, e.g. 'MATR3' → 'MATR3-1'.
```

and records it in `var_manifest.json` under `gene_name_sets` (see
{doc}`09_output_format_specification_exports_directory`).

### Drawability (the rule that still refuses)

The **exported** names must read on screen as the value they store. A name
carrying a control character, a zero-width character, or leading or trailing
whitespace fails before the export is published, and the message names it:
`'MS4A1 '` with a trailing space would be drawn identically to `'MS4A1'`.
`HLA-DRB1/2` is a real HGNC symbol and exports fine: it is never a path. This is
the same scoping `obs_keys` already has (see {doc}`04_obs_cell_metadata`), and
`cellucid_prepare()` in the R package applies the identical rule.

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

Identifiers are matched against the default names after they were made unique,
so the second `MATR3` is requested as `MATR3-1`. If any requested gene is
absent, the exporter raises `KeyError` and publishes nothing. A repeated request
raises `ValueError` naming the repeated identifier.

`gene_identifiers` is also the supported way past a gene the exporter cannot
draw: a name carrying invisible characters only has to be cleaned when the gene
is among the exported ones.

Reproducibility tip:
- store the exact gene list used for export in your pipeline (and ideally version-control it).

### Large `var` tables

The exporter exports no `var` column except the ones you name as gene-name
sources.

So, large additional columns in `adata.var` do not affect export size unless you
use them as names. However:
- large object columns can still slow down your own pre-processing,
- and they can be a privacy risk if you accidentally use them as names.

---

## Naming rules

Exported gene names never become paths, so there is no character set, length
limit, reserved name, or case-collision rule. `Gene` and `gene` are two distinct
genes and export as two payloads.

Names are completed only by the make-unique rule above — an empty value takes
its index name, and a repeated one gains a `-k` suffix. Nothing else is
rewritten: an exported name must be text the viewer can draw exactly as it is
stored, and a name that breaks this rejects the complete candidate before
publication.

---

## Troubleshooting (var / gene ids)

### Symptom: gene search returns nothing

Likely causes:
- You did not export `gene_expression` at all (no `var_manifest.json`).
- You exported a small `gene_identifiers` list and the queried gene isn’t included.
- Gene names in the UI are not what you expected (wrong `var_gene_id_column`).
- The name was repeated in `var` and is exported as `NAME-1`.

How to confirm:
- Does `<out_dir>/var_manifest.json` exist?
- Open it and read `fields`: each entry is `[index, name]`, or
  `[index, name, minValue, maxValue]` when quantized. `gene_name_sets` says
  which column each name set came from and what was renamed.

Fix:
- Export gene expression, or export the genes you need, and re-export with `force=True`.

### Symptom: `var column '…': Name at position 0 must be a string, got int.`

Meaning:
- the selected column (or `var.index`) holds non-string values.

Fix:
- convert it with `.astype(str)` if those values really are the names to show.

### Symptom: the “wrong gene” appears (values don’t match expectations)

Likely cause:
- `var` row order does not match `gene_expression` column order.

Fix:
- Rebuild `var` and `gene_expression` together (AnnData usually prevents this bug if used correctly).

---

## Next steps

- Expression export details (quantization, size, sparse vs dense): {doc}`06_gene_expression_matrix`
