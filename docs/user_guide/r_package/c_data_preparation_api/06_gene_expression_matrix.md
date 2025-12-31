# Gene Expression Matrix

**Audience:** computational users  
**Time:** 15–30 minutes  
**Goal:** export gene expression correctly without accidentally producing terabytes of files.

Gene expression export is optional, but it is what enables:
- gene search
- “color by gene” overlays

In `cellucid_prepare()`:
- `gene_expression` is the values matrix
- `var` provides the gene IDs / metadata

## Required shape and orientation (critical)

`gene_expression` must be shaped:

> `(n_cells, n_genes)` = **cells × genes**

That means:
- each **row** is one cell
- each **column** is one gene

```{warning}
Many R containers store expression as genes × cells. You almost always need to transpose.

- Seurat: `GetAssayData(...)` is typically genes × cells → use `Matrix::t(...)`
- SingleCellExperiment: `assay(...)` is genes × cells → use `Matrix::t(...)`
```

If you provide `gene_expression`, export enforces:
- `nrow(gene_expression) == n_cells`
- `ncol(gene_expression) == nrow(var)`

## Supported types

Supported inputs:
- base R `matrix`
- `Matrix` sparse matrices (recommended)

Special case:
- `dgCMatrix` is handled efficiently (direct column access without densifying the whole matrix).

## How data is written (what files are produced)

Gene expression values are written as **one file per gene** under:

- `<out_dir>/var/`

For a gene ID like `MS4A1`, you will get one of:

- float32: `var/MS4A1.values.f32`
- 8-bit quantized: `var/MS4A1.values.u8`
- 16-bit quantized: `var/MS4A1.values.u16`

The manifest:
- `<out_dir>/var_manifest.json`

maps gene IDs to file patterns and (if quantized) stores the min/max needed to dequantize.

## Quantization (`var_quantization`)

Quantization is the main disk-size lever for gene expression.

Set:
- `var_quantization = 8` (smallest; fastest I/O; lower precision)
- `var_quantization = 16` (more precision; still much smaller than float32)
- `var_quantization = NULL` (float32; largest)

```r
cellucid_prepare(..., var_quantization = 8)
```

Missing/invalid values (`NA`, `Inf`, `-Inf`) are mapped to a reserved marker:

| Quantization | Valid range | Missing marker |
|---|---|---|
| 8-bit | `0..254` | `255` |
| 16-bit | `0..65534` | `65535` |

## The “sparse matrix” misconception

Even if your input matrix is sparse, the exported per-gene file is a **dense vector of length `n_cells`**.

That means the total exported size scales like:

> `O(n_cells * n_genes)`

Approximate per-gene sizes (before gzip compression):

| dtype | bytes per cell | per gene @ 100k cells |
|---|---:|---:|
| `uint8` | 1 | ~0.10 MB |
| `uint16` | 2 | ~0.20 MB |
| `float32` | 4 | ~0.40 MB |

Multiply by the number of exported genes:
- 20k genes × 0.10 MB ≈ **2 GB** (100k cells, 8-bit)
- 20k genes × 0.40 MB ≈ **8 GB** (100k cells, float32)

For 1M cells, multiply those numbers by 10.

### Practical mitigation strategies

1) Export fewer genes (recommended)
   - pass `gene_identifiers = ...`
   - export HVGs, marker genes, or a curated panel
2) Use 8-bit quantization (`var_quantization = 8`)
3) Use gzip compression (`compression = 6`)

See the full performance guide: {doc}`10_performance_tuning_prepare_export`

## Choosing which expression values to export

From Seurat/SCE you often have choices like:
- raw counts
- log-normalized values
- scaled values

Cellucid can visualize any numeric values, but interpretation differs.

Recommendation:
- export log-normalized expression for visualization
- keep raw counts for analysis pipelines (not necessarily for the viewer)

## Edge cases

### Negative values

Negative expression values export fine. Quantization uses min/max scaling so negatives are representable.

### Extremely large outliers

A single extreme value can stretch min/max and reduce effective contrast for most cells.
If you see “gene coloring is washed out”, consider:
- clipping values before export, or
- exporting a transformed expression (e.g., log1p)

### Invalid values (`NA`, `Inf`)

If you export float32 (`var_quantization = NULL`), invalid values are written as-is (NaN/Inf).
The viewer may handle them, but it’s safer to clean them ahead of time.

## Troubleshooting pointers

- “var has X rows but gene_expression has Y genes” → orientation mismatch.
- “Export is huge / takes forever” → you exported too many genes; use `gene_identifiers` + quantization.
- “Gene IDs look weird” → your `rownames(var)` are missing; set them explicitly.
- Full troubleshooting: {doc}`11_troubleshooting_prepare_export`
