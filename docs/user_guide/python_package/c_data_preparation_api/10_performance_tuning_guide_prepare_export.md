# Performance tuning guide (prepare/export)

**Audience:** computational users and power users  
**Time:** 30–60 minutes (plus time to run exports)  
**Goal:** make exports fast enough to iterate on and small enough to share

This guide focuses on export-time performance for {func}`~cellucid.prepare`.

If you’re new: defaults are okay for small datasets, but large datasets require intentional choices.

---

## Fast path (safe defaults)

For many small-to-medium datasets, a good starting point is:

- `compression=6`
- `var_quantization=8`
- `obs_continuous_quantization=8`
- export a curated gene list via `gene_identifiers`
- export a curated field list via `obs_keys`

And always:
- use `force=True` while iterating, or export to a new `out_dir` each run.

---

## Mental model: what gets big and why

Roughly:

- **Embeddings (points)**: small (`n_cells × 2/3` float32)
- **Obs fields**: usually manageable (one vector per field)
- **Gene expression**: potentially enormous (**one vector per gene**, each length `n_cells`)
- **Connectivities**: moderate to huge depending on graph density (`n_edges`)
- **Vector fields**: moderate (one vector per field per dimension)

The two most common failure modes:
1) exporting too many genes
2) re-exporting repeatedly without controlling `out_dir` / `force`

---

## Compression (`compression`)

`compression` is gzip level 1–9.

Tradeoffs:
- Higher compression → smaller exports but slower export time (CPU-bound).
- Lower compression → faster export time but larger exports (I/O-bound).

Practical guidance:
- While iterating: disable compression (`compression=None`) to save time.
- For final shared exports: `compression=6` is a good balance.

Compatibility note:
- exporting with compression requires browser-side gzip decompression support.
  Modern browsers support this; if you need maximum compatibility, export uncompressed.

---

## Quantization (`var_quantization`, `obs_continuous_quantization`)

Quantization is the biggest “size knob” that does not change feature availability.

### `var_quantization` (gene expression)

- `None` → float32 (largest, most faithful)
- `8` → uint8 per gene (smallest)
- `16` → uint16 per gene (middle ground)

When to prefer 16-bit:
- subtle gradients are important (e.g., low-expression genes),
- users will compare similar genes visually and you want less banding.

### `obs_continuous_quantization` (continuous obs + categorical outlier quantiles)

Same idea, applied to:
- continuous obs fields
- outlier quantiles (for categorical fields)

Practical tip:
- if you enable `obs_continuous_quantization`, you also reduce the size of outlier quantile files.

```{important}
Quantization is lossy by design. It is appropriate for visualization, not for downstream analysis.
Keep analysis in your source object (AnnData) and treat exports as view artifacts.
```

---

## Subsetting (the most important knob for large datasets)

### Subset obs fields (`obs_keys`)

Export only fields you intend users to interact with.

This improves:
- load time (smaller manifest + fewer files),
- UI usability (shorter field list),
- privacy (avoid accidental leakage of identifiers).

### Subset genes (`gene_identifiers`)

This is usually required for large datasets.

Exporting all genes:
- creates one file per gene,
- can create tens of thousands of files,
- and can generate massive disk usage.

Recommended strategies:
- export marker genes + HVGs (e.g., 1k–5k genes),
- export a panel tailored to the expected audience (wet lab-friendly marker set),
- keep a “full gene” workflow in server mode instead of static exports.

---

## Scaling guidance by dataset size (rule of thumb table)

These are *workflow* recommendations, not strict limits.

| Dataset scale | Recommended export strategy |
|---|---|
| `< 50k` cells | Export genes if you want; `var_quantization=8`, `compression=6` often fine |
| `50k–200k` cells | Export curated genes only; avoid exporting all genes; quantize + compress |
| `200k–1M` cells | Strongly prefer curated genes or server mode; compression may become CPU-heavy |
| `> 1M` cells | Prefer server mode / AnnData-backed loading; static export should be “metadata-first” and very selective |

If you need full gene access on a large dataset, use server mode:
- {doc}`../d_viewing_apis/09_server_mode_advanced`

---

## Iteration workflow (recommended)

A practical export workflow that saves time:

1) Export **embeddings + a few obs fields** to confirm geometry and alignment.
2) Add more obs fields and verify they behave as intended (categorical vs continuous).
3) Add vector fields (if needed) and validate in the viewer.
4) Only then add gene expression, starting with a tiny gene list.
5) Expand gene list gradually; watch export size/time and UI responsiveness.

---

## Deep path (expert notes)

### File count matters as much as byte size

Even if quantized + compressed, exporting 20,000 genes means 20,000 files.
On some systems (especially networked filesystems), many small files are slower than fewer large ones.

If you hit file-count pain:
- export fewer genes,
- or use server mode (lazy loading).

### Connectivity export can dominate runtime

The exporter currently extracts edges in Python loops.
For very large graphs, this can be slow.

If connectivities are not essential, skip them in static exports.

---

## Edge cases and common footguns

- Exporting to a slow filesystem (network mount) makes the “many files” problem much worse.
- Enabling high gzip compression while exporting many genes can be painfully CPU-bound.
- Reusing `out_dir` without `force=True` can produce stale manifests and confusing mismatches.

---

## Troubleshooting (performance)

### Symptom: export is extremely slow

Likely causes:
- exporting many genes (file count),
- high compression level,
- exporting connectivities on a huge graph,
- writing to slow storage.

Fixes (in order):
- export fewer genes (`gene_identifiers`),
- disable compression during iteration,
- skip connectivities,
- export to local SSD.

### Symptom: export folder is huge

Likely causes:
- exporting too many genes,
- quantization disabled,
- compression disabled.

Fix:
- subset genes, enable quantization and compression.

---

## Next steps

- Deep troubleshooting index: {doc}`11_troubleshooting_prepare_export`
