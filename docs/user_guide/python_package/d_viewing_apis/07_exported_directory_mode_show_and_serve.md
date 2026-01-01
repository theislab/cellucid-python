# Exported directory mode (`show()` / `serve()`)

This page explains the “exported directory” workflow: you run `cellucid.prepare(...)` once, and then you can view that folder quickly and reproducibly in Cellucid.

If you do not have an export directory yet, start with:
- export overview: {doc}`../c_data_preparation_api/01_prepare_export_overview`
- input requirements + edge cases: {doc}`../c_data_preparation_api/02_input_requirements_global`

## At a glance

**Audience**
- Wet lab / non-technical: exports are a “dataset package” you can open without worrying about AnnData internals.
- Computational: exports are the best default for performance + sharing + reproducibility.

**Prerequisites**
- An export directory produced by `cellucid.prepare(...)`

## Why exports are recommended

Exports are:
- **fast** in the browser (binary files optimized for the viewer),
- **stable** (the same folder loads the same way every time),
- **shareable** (copy the folder, host it, archive it with a paper),
- and the easiest way to avoid RAM-heavy notebook workflows.

AnnData direct mode is great for exploration, but exports are what you want for production, collaboration, and very large datasets.

## Quick start: view an export directory

Pick the environment you’re in.

### Terminal (recommended default)

```bash
cellucid serve /path/to/export_dir
```

### Python (script)

```python
from cellucid import serve

serve("/path/to/export_dir")  # blocks until Ctrl+C
```

### Jupyter (embed)

```python
from cellucid import show

viewer = show("/path/to/export_dir", height=600)
```

## What’s inside an export directory (high-level)

An export directory is a folder containing:

- `dataset_identity.json` (dataset metadata + embeddings + stats)
- `obs_manifest.json` (cell metadata fields list + schemas)
- `var_manifest.json` (gene IDs list + gene expression schema)
- one or more embedding binaries:
  - `points_1d.bin` / `points_2d.bin` / `points_3d.bin` (optionally `.gz`)
- optional subfolders for larger artifacts:
  - `obs/` (per-field values/codes)
  - `var/` (per-gene values)
  - `connectivity/` (edge lists)
  - `vectors/` (vector fields like velocity/drift)

```{note}
Compression: `prepare(..., compression=N)` writes `.gz` files on disk. The viewer reads the exact file paths from manifests/identity, so you don’t “turn on gzip” manually at view time.
```

## Validate an export folder (before debugging the viewer)

If you want a fast “is this export even complete?” check:

### 1) Check required files exist

At minimum, you should see:

- `dataset_identity.json`
- `obs_manifest.json`
- at least one of:
  - `points_1d.bin` / `points_1d.bin.gz`
  - `points_2d.bin` / `points_2d.bin.gz`
  - `points_3d.bin` / `points_3d.bin.gz`

### 2) Confirm the server can serve the identity file

After starting the server, open:

```text
http://127.0.0.1:<port>/dataset_identity.json
```

If `dataset_identity.json` fails to load, the viewer cannot load the dataset.

### 3) Sanity-check stats and embeddings (optional)

Open `dataset_identity.json` and confirm:

- `stats.n_cells` matches your expectation
- `embeddings.available_dimensions` matches what you exported
- `embeddings.files` contains `points_2d...` and/or `points_3d...` paths that actually exist on disk

## Serving multiple exported datasets (one server, many folders)

You can point the server at a directory that contains multiple export subfolders:

```text
exports/
  pbmc_demo/
    dataset_identity.json
    obs_manifest.json
    points_3d.bin.gz
    ...
  pancreas_demo/
    dataset_identity.json
    obs_manifest.json
    points_2d.bin.gz
    ...
```

Then:

```bash
cellucid serve /path/to/exports
```

Debug endpoint:

```text
http://127.0.0.1:<port>/_cellucid/datasets
```

```{note}
The server detects export directories by checking for key files (manifests + points). If your folder is missing those, it may not show up in the dataset list even if the directory exists.
```

## Sharing exports (without giving someone your Python environment)

Two common sharing modes:

1) **Share privately** (zip, shared drive, internal server) and have collaborators open it via:
   - `cellucid serve <export_dir>` (fast), or
   - web app folder picker (easy): {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`

2) **Host publicly** (GitHub Pages / static hosting):
   - end-to-end tutorial: {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`

## Edge cases (exports)

- **Missing gene expression**: you can export embeddings + obs without exporting `gene_expression`. The viewer can load, but gene search/expression overlays won’t work.
- **Huge categorical fields**: very large numbers of categories can be unusable in the UI; consider coarsening or filtering before export.
- **Dimension mismatch**: exporting only 2D points means 3D view and 3D vector fields cannot appear.

## Troubleshooting

- Viewer loads but some data is missing → validate `dataset_identity.json` and manifests.
- Viewer fails to load at all → start with {doc}`15_troubleshooting_viewing`.

## Next steps

- AnnData direct mode: {doc}`08_anndata_mode_show_anndata_and_serve_anndata`
- Performance guide: {doc}`14_performance_scaling_and_lazy_loading`
