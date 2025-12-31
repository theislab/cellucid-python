# Performance Tuning (Prepare/Export)

**Audience:** anyone exporting medium/large datasets  
**Time:** 15–25 minutes  
**Goal:** avoid exports that are unreasonably slow, huge, or brittle.

Export performance is mostly about:

1) **How much data you export** (cells × genes × fields)
2) **How many files you create** (especially one file per gene)
3) **Your filesystem speed** (SSD vs network drive)

This page gives practical guidance for `cellucid-r`.

## What scales with what?

### Embeddings

- Work scales as `O(n_cells * dim)` (cheap)
- Disk size is tiny compared to gene expression

### Obs fields

- Work scales as `O(n_cells * n_obs_fields)`
- Categorical fields also compute centroids/outlier quantiles (still manageable for typical sizes)

### Gene expression (the big one)

Gene expression export scales like:

- **Time:** roughly `O(n_cells * n_exported_genes)`
- **Disk:** roughly `n_cells * n_exported_genes * bytes_per_value`
- **File count:** `n_exported_genes` files

If you export 20k genes, you create ~20k files in `var/`.

### Connectivities

Cost depends on graph density.
Sparse graphs export fine; dense graphs explode.

## The three biggest levers (do these first)

### 1) Export fewer genes (best lever)

Use `gene_identifiers` to export a curated gene panel:

```r
panel <- c("MS4A1", "CD3D", "NKG7")
cellucid_prepare(..., gene_identifiers = panel)
```

Recommended panels:
- marker genes from your analysis
- HVGs
- curated pathway panels

### 2) Use quantization

For gene expression:
- `var_quantization = 8` is usually a good default for visualization

For continuous metadata:
- `obs_continuous_quantization = 8`

```r
cellucid_prepare(
  ...,
  var_quantization = 8,
  obs_continuous_quantization = 8
)
```

### 3) Use gzip compression (if disk is the bottleneck)

Compression reduces size but can increase CPU time.

Rule of thumb:
- `compression = 1` for speed
- `compression = 6` for a reasonable tradeoff (often a good default)
- `compression = 9` for smallest output (slowest)

```r
cellucid_prepare(..., compression = 6)
```

## Recommended settings by dataset scale (rule-of-thumb)

| Dataset scale | Suggested export strategy |
|---|---|
| ≤ 50k cells | exporting many genes may be OK; still prefer `var_quantization=8` |
| 50k–200k cells | export a gene panel; use quantization; consider gzip |
| > 200k cells | strongly prefer gene panels; full gene export can become impractical |

```{note}
If you truly need “all genes at very large cell counts”, consider the Python server workflow where gene expression can be loaded lazily from an AnnData/Zarr backend.
```

## Filesystem and OS tips

- Export to a **local SSD** (not a network drive) whenever possible.
- Avoid paths synced by cloud clients (Dropbox/OneDrive) during export.
- If you hit “too many files” limits on some filesystems, export fewer genes or split datasets.

## Avoiding “stale export” confusion

If you re-export to the same `out_dir` with `force=FALSE` (default), the exporter may skip writing files.

Recommendations:
- set `force=TRUE` while iterating, or
- export to a new `out_dir` each time (best for reproducibility)

## Troubleshooting pointers

- Export is huge → you exported too many genes; use `gene_identifiers` + quantization.
- Export is slow → same root cause, plus filesystem limitations.
- See: {doc}`11_troubleshooting_prepare_export`
