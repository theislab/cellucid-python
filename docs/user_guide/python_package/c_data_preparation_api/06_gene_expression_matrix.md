# Gene expression matrix

**Audience:** everyone exporting genes (computational users will care most)  
**Time:** 30–90 minutes depending on dataset size and number of genes exported  
**Goal:** export gene expression in a way that is accurate for visualization and manageable in size/time

Exporting gene expression is optional — but it’s what enables:
- gene search,
- gene overlays (color-by gene),
- and gene-driven analysis features in the web app.

It is also the easiest way to create a gigantic, slow export if you do it naively.

---

## Fast path (a practical, safe default)

For small/medium datasets, a good starting point is:

- export a curated gene list (markers/HVGs),
- enable 8-bit quantization,
- enable gzip compression.

```python
from cellucid import prepare

marker_genes = ["MS4A1", "CD3D", "LYZ", "NKG7"]

prepare(
    ...,
    gene_expression=adata.X,
    var=adata.var,
    var_gene_id_column="index",
    gene_identifiers=marker_genes,
    var_quantization=8,
    compression=6,
    force=True,
)
```

If you try to export all genes for a large dataset and things get slow or huge, jump to:
- {doc}`10_performance_tuning_guide_prepare_export`

---

## Practical path (computational users)

### Supported input types and required shape

`gene_expression` may be:
- a dense `numpy.ndarray`, or
- a `scipy.sparse` matrix.

Required shape:

```text
(n_cells, n_genes)
```

AnnData’s `adata.X` is typically already `(n_cells, n_genes)`.

If you are not using AnnData, be careful:
- many pipelines store expression as `(n_genes, n_cells)` and require a transpose before export.

### Alignment with `var`

If you provide `gene_expression`, you must also provide `var`, and:

- `len(var) == gene_expression.shape[1]`
- `var.iloc[j]` corresponds to `gene_expression[:, j]`

This is the most common cause of “wrong gene values”.

### What gets written on disk (important for performance)

Cellucid exports expression as **dense per-gene vectors**:

- for each gene, write a vector of length `n_cells` under `var/`
- write an index (`var_manifest.json`) so the web app can fetch genes on demand

This means:
- the number of files under `var/` is approximately the number of exported genes,
- exporting “all genes” can create tens of thousands of files,
- even if your input matrix is sparse, the on-disk representation is per-gene dense vectors (gzip helps for many zeros).

### Quantization (`var_quantization`)

`var_quantization` controls how values are stored:
- `None` → store `float32` (lossless for float32)
- `8` → quantize each gene to `uint8` (lossy; ~4× smaller than float32)
- `16` → quantize each gene to `uint16` (lossy; ~2× smaller than float32)

Quantization rules (current exporter):
- quantization is **per gene** (each gene gets its own min/max),
- NaN/Inf are mapped to a reserved missing marker:
  - `255` for 8-bit, `65535` for 16-bit,
- `minValue`/`maxValue` are stored in `var_manifest.json` and used for dequantization.

Dequantization in the web app is:

```text
value = minValue + q * (maxValue - minValue) / maxQuant
```

Where `maxQuant = 254` (8-bit) or `65534` (16-bit).

Practical guidance:
- For visualization and interactive exploration, 8-bit is usually fine.
- If users rely on subtle gradients (scores, near-zero expression differences), prefer 16-bit or float32.

### Missing/invalid values

- NaN/Inf are allowed and will be treated as missing in the UI.
- Negative values are allowed (quantization uses min/max and will encode them).

If you see many NaNs/Inf:
- confirm your preprocessing (e.g., division by zero, log of negative values),
- consider exporting a different representation.

### Choosing what expression to export (counts vs normalized)

Cellucid does not decide what “expression” means.
You choose the matrix you export.

Common choices:
- log1p-normalized expression (good for visualization)
- scaled/z-scored expression (good for certain contrasts; includes negative values)
- raw counts (often dominated by library size; not ideal for color-by without normalization)

Make the choice explicit in your pipeline and in dataset metadata.

---

## Size and performance expectations (rules of thumb)

The raw data volume is roughly:

- float32: `4 * n_cells * n_genes` bytes
- 8-bit: `1 * n_cells * n_genes` bytes
- 16-bit: `2 * n_cells * n_genes` bytes

But real-world exports differ because:
- the exporter writes **one file per gene** (filesystem overhead),
- gzip compression can greatly reduce size for sparse-ish data,
- but compression increases CPU time.

Practical implication:
- Large datasets (hundreds of thousands of cells) should not export “all genes” unless you have a very specific reason.

If you need full gene access on large datasets, prefer server mode:
- {doc}`../d_viewing_apis/09_server_mode_advanced`

---

## Edge cases and common footguns

- **Wrong orientation** (`genes × cells`): fix by transposing to `(cells × genes)`.
- **Mismatch between `var` and matrix columns**: leads to wrong gene names/values.
- **Duplicate gene IDs**: can overwrite files and produce confusing UI results (see {doc}`05_var_gene_metadata`).
- **All-zero genes**: export is valid but gene overlays will be flat.
- **NaNs introduced by preprocessing**: common after invalid log transforms or normalization artifacts.
- **Huge file counts**: tens of thousands of gene files can be slow on some filesystems (especially networked).

---

## Troubleshooting (gene expression)

### Symptom: gene search is missing / disabled

Meaning:
- `var_manifest.json` is missing (you didn’t export genes or export was skipped).

Fix:
- export with `gene_expression` + `var`, and re-export with `force=True`.

### Symptom: export folder exploded in size

Likely causes:
- exported too many genes,
- quantization disabled (float32),
- compression disabled.

Fix:
- export a curated gene list,
- enable `var_quantization=8` and `compression=6`,
- or use server mode instead of static exports for full gene access.

### Symptom: export is extremely slow

Likely causes:
- large `n_genes` (file count),
- compression enabled on a slow CPU,
- writing to a slow disk/network filesystem.

Fix:
- export fewer genes,
- disable compression during iteration,
- write to a local SSD and move the export later.

---

## Next steps

- Optional KNN graph export: {doc}`07_connectivities_knn_graph`
- Performance and scaling guidance: {doc}`10_performance_tuning_guide_prepare_export`
