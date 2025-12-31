# What Is an Export Folder?

**Audience:** everyone  
**Time:** 5–10 minutes

An **export folder** is the on-disk representation of a dataset that the Cellucid web app can load.

Think of it as a “dataset bundle”:
- It contains **binary arrays** for fast loading (points, values, edges).
- It contains **small JSON manifests** that tell the viewer what those binaries mean.

## What you typically see in an export folder

Minimal export (no gene expression, no connectivity):

```
<out_dir>/
  points_2d.bin
  obs_manifest.json
  dataset_identity.json
  obs/
    <field>.values.f32
    <field>.codes.u8
    <field>.outliers.f32
    ...
```

If you include gene expression and connectivities, you will also see:

```
<out_dir>/
  var_manifest.json
  var/
    <gene>.values.f32
    ...
  connectivity_manifest.json
  connectivity/
    edges.src.bin
    edges.dst.bin
```

If you include vector fields (velocity/displacement):

```
<out_dir>/
  vectors/
    <field>_2d.bin
    <field>_3d.bin
    ...
```

```{note}
If you set `compression=...`, the binary files are written as `*.gz` (for example `points_2d.bin.gz`).
The JSON manifests record whether compression is enabled.
```

## Why so many files?

This export format optimizes for:

- **Interactive viewing** (load only what you need)
- **Simple hosting** (static file hosting works)
- **Compatibility** across languages (R/Python produce the same format)

The biggest design consequence is:

> Gene expression is stored as **one file per gene** (`var/<gene>.values.*`).

That makes gene lookup simple in the browser, but it can create tens of thousands of files for typical datasets.

## “Cell identity” is the row order

In the export format:
- Cell `0` is the first row of your embedding and metadata.
- Cell `1` is the second row.
- …and so on.

There is no separate “cell IDs” file; the browser just assumes every exported array shares the same row order.

This is why alignment matters so much in Seurat/SCE exports (see {doc}`../a_landing_pages/03_supported_inputs_and_workflows`).

## How to inspect an export folder in R

### Quick inventory

```r
list.files(out_dir, recursive = TRUE)
```

### Read the “identity” file

```r
ident <- jsonlite::read_json(file.path(out_dir, "dataset_identity.json"), simplifyVector = TRUE)
str(ident, max.level = 2)
```

What to look for:
- `stats$n_cells`, `stats$n_genes`
- `embeddings$available_dimensions`
- `obs_fields` (summary of categorical vs continuous fields)
- `export_settings` (compression/quantization)

## Common edge cases

- **No gene expression**: `var_manifest.json` and `var/` are absent; gene search/overlays are unavailable.
- **Only 3D UMAP**: you will have `points_3d.bin` but not `points_2d.bin`.
- **Compression enabled**: binaries are `.gz`; make sure you host them with correct MIME types (most servers do).

## Troubleshooting pointers

- “The folder exists but the web app can’t load it” → {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`
- “I’m missing expected files” → {doc}`../i_troubleshooting_index/03_export_format_and_validation_issues`
