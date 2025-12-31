# What is cellucid-r?

**Audience:** wet lab scientists, computational biologists, engineers  
**Time:** 5–10 minutes  
**Goal:** understand what the R package does, and what workflow it supports

`cellucid-r` (package name: `cellucid`) is a **minimal R package** whose job is to export your single-cell data to the **on-disk format** consumed by the **Cellucid web app**.

If you remember one sentence:

> `cellucid-r` turns *R objects* into an *export folder* that you open in your browser (Cellucid).

```{note}
Cellucid itself is a **web app** (the viewer UI). `cellucid-python` and `cellucid-r` are **helper packages** that prepare data and (in Python) provide notebook/server integrations. `cellucid-annotation` is for community annotation workflows.
```

## What `cellucid-r` does (today)

### Fast path (wet lab / non-technical)
- You (or a collaborator) run one R function: `cellucid::cellucid_prepare(...)`.
- It writes an **export folder** to disk (a directory full of `.json` + binary files).
- You open the Cellucid web app and load that folder.
- You interactively explore: embeddings, metadata coloring, gene expression overlays, etc.

### Practical path (computational users)
`cellucid_prepare()` writes:
- **Embeddings**: `points_1d.bin`, `points_2d.bin`, `points_3d.bin` (float32, row-major; optionally gzipped)
- **Cell metadata (“obs”)**: `obs_manifest.json` + `obs/` binary files
- **Gene expression (“var”)** (optional): `var_manifest.json` + `var/` binary files (one file per gene)
- **Connectivity graph** (optional): `connectivity_manifest.json` + `connectivity/edges.*.bin`
- **Dataset identity**: `dataset_identity.json` (what the viewer uses to understand “what’s in this folder”)
- **Vector fields** (optional): `vectors/*.bin` + metadata in `dataset_identity.json`

It is intentionally dependency-light:
- Hard dependency: `jsonlite`
- Optional (recommended) dependency: `Matrix` for sparse matrices and connectivity export

### Deep path (expert / maintainer)
`cellucid_prepare()` also:
- **Normalizes embeddings** to a stable coordinate scale (so different datasets render consistently).
- Computes **categorical centroids** (in embedding space) and **outlier quantiles** (in latent space) to support UI features that need category summaries.
- Supports **quantization** (8-bit or 16-bit) for continuous values to reduce disk size and speed up loading.

## What `cellucid-r` does NOT do (yet)

Be explicit about current boundaries:

- It does **not** embed Cellucid inside RStudio / Shiny (no widget yet).
- It does **not** run a dedicated R server for lazy gene loading (Python has this).
- It does **not** depend on Seurat or SingleCellExperiment; you must **extract matrices/data.frames** yourself (we provide recipes).

If you want a “viewer object” you can control from a notebook, see {doc}`../../python_package/index`.

## The core mental model

Cellucid uses a simple, robust mapping:

1) **Row order defines cell identity.**  
   Cell `i` in your embeddings must be the same cell `i` in your `obs`, `latent_space`, gene expression matrix, connectivity matrix, and vector fields.
2) **The export folder is the dataset.**  
   The browser doesn’t “load an R object” — it loads files.
3) **The web app is the UI.**  
   R is for preparation. The browser is for interaction.

## Common misconceptions (quick clarifications)

- “Can I export from Seurat/SCE directly?”  
  Not as a single call (by design), but you can export by extracting the data. Start with {doc}`../e_integrations_recipes/index`.

- “Do I need gene expression to use Cellucid?”  
  No. You can export embeddings + metadata only. Gene expression enables gene overlays/search.

- “Why is `latent_space` required?”  
  `cellucid-r` uses `latent_space` to compute per-cell **outlier quantiles** for categorical fields. This supports UI features that want a robust “how typical is this cell for its category?” measure.

## Next steps

- If you want a copy/paste minimal example: {doc}`04_quick_start_3_levels`
- If you want Seurat or SingleCellExperiment recipes: {doc}`../e_integrations_recipes/index`
- If you want the detailed API + file format documentation: {doc}`../c_data_preparation_api/index`
