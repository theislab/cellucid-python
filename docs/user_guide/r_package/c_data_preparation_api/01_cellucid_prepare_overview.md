# `cellucid_prepare()` / `prepare()` Overview

**Audience:** everyone  
**Time:** 10–20 minutes  
**Goal:** understand what the exporter does and how to call it correctly

`cellucid::cellucid_prepare()` is the main entry point of `cellucid-r`. It writes an **export folder** containing binaries + JSON manifests that the Cellucid web app can load.

`cellucid::prepare()` is an alias for the same function.

## Minimal “mental model”

1) You provide arrays/data.frames.
2) Cellucid normalizes/encodes them for the viewer.
3) It writes files to `out_dir`.
4) You open `out_dir` in the browser.

The function returns `NULL` (invisibly). It is called for side effects (writing files).

## A minimal call (embeddings + obs only)

```r
cellucid::cellucid_prepare(
  latent_space = latent,  # (n_cells, n_latent_dims)
  obs = obs,              # data.frame with n_cells rows
  X_umap_2d = umap2,      # (n_cells, 2)
  out_dir = "exports/my_dataset",
  force = TRUE
)
```

## What is required vs optional?

### Required

- `latent_space` (matrix-like, `(n_cells, n_dims)`)
- `obs` (data.frame with `n_cells` rows)
- at least one embedding (`X_umap_1d` or `X_umap_2d` or `X_umap_3d`)

### Optional (common)

- `gene_expression` + `var` (enables gene overlays/search)
- `connectivities` (enables graph-based features)
- `vector_fields` (enables velocity/displacement overlays)

### Optional (performance/size knobs)

- `compression` (gzip level 1–9)
- `var_quantization` (8/16-bit for gene expression)
- `obs_continuous_quantization` (8/16-bit for continuous obs + outlier quantiles)
- `gene_identifiers` (export a subset of genes)
- `obs_keys` (export a subset of metadata columns)

## Key behaviors you should know (before you export real data)

### 1) Embeddings are normalized

Each embedding is centered and scaled so it fits a stable range for rendering.
Vector fields are scaled consistently with this normalization.

Details: {doc}`03_embeddings_and_coordinates`

### 2) `latent_space` is used for categorical “outlier quantiles”

For each categorical field, `cellucid-r` computes a per-cell outlier quantile (distance rank inside its category in latent space). This requires `latent_space`.

Details: {doc}`04_obs_cell_metadata`

### 3) Export can “silently skip” work unless `force=TRUE`

By default `force=FALSE`. If manifests already exist, export may skip:
- `obs_manifest.json`
- `var_manifest.json`
- `connectivity_manifest.json`
- and embedding/vector files that already exist

This is useful for incremental work, but it can be confusing when you expect outputs to change.

Rule of thumb:
- use a new `out_dir` for each export iteration, or
- set `force=TRUE` while you are iterating.

### 4) Continuous quantization uses a reserved missing marker

If you quantize continuous values:
- valid values map to `0..254` (8-bit) or `0..65534` (16-bit)
- invalid values (`NA`, `Inf`, `-Inf`) map to `255` / `65535`

Details: {doc}`04_obs_cell_metadata` and {doc}`06_gene_expression_matrix`

### 5) Gene expression export is dense-per-gene

Even if your input expression matrix is sparse, each exported gene file is a **dense vector of length `n_cells`**. This is the main reason large exports can be huge.

Details and mitigation strategies: {doc}`06_gene_expression_matrix` and {doc}`10_performance_tuning_prepare_export`

## What files are written?

At minimum (obs + embeddings):
- `points_*d.bin[.gz]`
- `obs_manifest.json`
- `obs/*`
- `dataset_identity.json`

If you include gene expression:
- `var_manifest.json`
- `var/*`

If you include connectivities:
- `connectivity_manifest.json`
- `connectivity/edges.src.bin[.gz]`
- `connectivity/edges.dst.bin[.gz]`

If you include vector fields:
- `vectors/*.bin[.gz]`
- metadata inside `dataset_identity.json`

Full format spec: {doc}`09_output_format_specification_exports_directory`

## Next steps

- Want the global alignment rules first? {doc}`02_input_requirements_global`
- Exporting from Seurat/SCE? {doc}`../e_integrations_recipes/index`
- Something failed? {doc}`11_troubleshooting_prepare_export`
