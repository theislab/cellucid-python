# Validate Exports and Debug Loading

**Audience:** computational users and anyone sharing datasets  
**Time:** 15–30 minutes  
**Goal:** catch export problems *before* you debug “web app issues”.

When the viewer fails to load a dataset, the root cause is often:
- a missing file,
- a broken manifest,
- or a cell-order/shape mismatch that only becomes obvious in the UI.

This page gives a practical validation checklist.

## 1) Validate the folder structure

At minimum, your dataset folder should contain:
- `dataset_identity.json`
- `obs_manifest.json`
- `obs/`
- at least one `points_*d.bin` (or `.gz`)

Quick inventory:

```r
list.files(out_dir)
list.files(file.path(out_dir, "obs"))
```

## 2) Validate `dataset_identity.json`

```r
ident <- jsonlite::read_json(file.path(out_dir, "dataset_identity.json"), simplifyVector = TRUE)
str(ident, max.level = 3)
ident$stats
```

Confirm:
- `stats$n_cells` is correct
- `embeddings$available_dimensions` matches what you exported (2D/3D)
- `stats$n_genes` matches what you expect (0 if no expression)

## 3) Validate `obs_manifest.json`

```r
obs_manifest <- jsonlite::read_json(file.path(out_dir, "obs_manifest.json"), simplifyVector = TRUE)
names(obs_manifest)
```

Confirm:
- `_continuousFields` and/or `_categoricalFields` exist (depending on your obs)
- `_obsSchemas` contains file patterns for `obs/`

## 4) Validate embedding binary length (quick check)

If you exported 2D points without compression:

```r
n_cells <- ident$stats$n_cells
path <- file.path(out_dir, "points_2d.bin")
stopifnot(file.exists(path))

con <- file(path, open = "rb")
on.exit(close(con), add = TRUE)

vals <- readBin(con, what = "numeric", size = 4, endian = "little", n = n_cells * 2)
stopifnot(length(vals) == n_cells * 2)
```

If you used compression, the path is `points_2d.bin.gz` and you should open with `gzfile(...)`.

## 5) Validate gene expression export (if present)

If `stats$n_genes > 0`, confirm:
- `var_manifest.json` exists
- `var/` contains files

```r
if (ident$stats$n_genes > 0) {
  stopifnot(file.exists(file.path(out_dir, "var_manifest.json")))
  gene_files <- list.files(file.path(out_dir, "var"))
  cat("gene files:", length(gene_files), "\n")
}
```

```{tip}
If you expected 20,000 genes and you see only a handful, you probably used `gene_identifiers` (gene panel) — which is often a good thing.
```

## 6) Validate connectivities (if present)

```r
if (isTRUE(ident$stats$has_connectivity)) {
  stopifnot(file.exists(file.path(out_dir, "connectivity_manifest.json")))
  stopifnot(file.exists(file.path(out_dir, "connectivity", "edges.src.bin")) ||
              file.exists(file.path(out_dir, "connectivity", "edges.src.bin.gz")))
}
```

## 7) Debug loading failures in the browser (practical)

### A) Start by choosing the right loading method

- Local single dataset folder → use file picker ({doc}`01_open_exports_in_cellucid_web_app`)
- GitHub sharing → exports root + `datasets.json` ({doc}`02_host_exports_for_sharing`)
- Private/large dataset → server mode ({doc}`../../web_app/b_data_loading/04_server_tutorial`)

### B) Use the browser developer tools

Open DevTools:
- Console: look for CORS errors, JSON parse errors
- Network: look for 404s for `dataset_identity.json`, `obs_manifest.json`, or binaries

Common signatures:

- **404 not found**:
  - wrong path / wrong folder selected / missing file
- **CORS blocked**:
  - you are trying to load from a custom host that doesn’t send CORS headers
- **JSON parse error**:
  - corrupted JSON or truncated download

### C) “It loads but looks wrong”

If points load but:
- clusters don’t match the point cloud,
- gene expression is “on the wrong cells”,

the most common cause is a **row order mismatch** in the export inputs.

Fix:
- re-export after explicitly aligning all inputs to a single `cells` order (see {doc}`../c_data_preparation_api/02_input_requirements_global`).

## Troubleshooting pointers

- Export-time problems: {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`
- Web app data-loading troubleshooting: {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`
