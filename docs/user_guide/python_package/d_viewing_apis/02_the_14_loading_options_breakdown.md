# The “14 loading options” breakdown

This page lists the **canonical 14 ways** Cellucid can load data.

Why this matters:
- It helps you choose the safest workflow for your situation (local, remote, notebook, sharing).
- It clarifies when you get **true lazy loading** (critical for large datasets).
- It reduces “it loads on my laptop but not on the cluster” surprises.

If you want a decision tree instead of a matrix, go to {doc}`03_choose_your_workflow_decision_tree`.

## Definitions (so the table is readable)

- **Exported**: a folder produced by `cellucid.prepare(...)` (fast, reproducible, shareable).
- **AnnData**: `.h5ad`, `.zarr`, or an in-memory `AnnData` object.
- **Lazy genes**: gene expression is fetched **on demand** (gene-by-gene) instead of “load the whole matrix”.
  - Lazy genes is the difference between “works for 50k cells” vs “works for 2M cells”.

## The 14 loading options (complete matrix)

```{note}
This matrix is referenced throughout the docs. Most people only use 2–3 of these in practice.
```

| # | Where you run things | How you point Cellucid to the data | Data format | Lazy genes | Best for |
|---:|---|---|---|---|---|
| 1 | Cellucid web app | Built-in demo dataset picker | Exported | ✅ | learning the UI with known-good data |
| 2 | Cellucid web app | Public GitHub export (`?github=...`) | Exported | ✅ | sharing publicly without running a server |
| 3 | Cellucid web app | Browser **folder picker** | Exported | ✅ | quick local viewing of prepared exports |
| 4 | Cellucid web app | Browser **.h5ad picker** | `.h5ad` | ❌* | quick preview of small `.h5ad` |
| 5 | Cellucid web app | Browser **Zarr ZIP picker** | `.zarr.zip` / `.zip` containing one Zarr v2 store | ✅† | portable Zarr viewing without Python |
| 6 | Terminal (CLI) | `cellucid serve <export_dir>` | Exported | ✅ | reliable, fast local viewing |
| 7 | Terminal (CLI) | `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset` | `.h5ad` | ✅ | large `.h5ad` with read-only backed access |
| 8 | Terminal (CLI) | `cellucid serve data.zarr --dataset-name "My dataset" --dataset-id my-dataset` | `.zarr` | ✅ | Zarr stores that fit an eager Python load |
| 9 | Python | `cellucid.serve(<export_dir>)` | Exported | ✅ | scripted server startup |
| 10 | Python | `cellucid.serve_anndata(<data.h5ad>, dataset_name="My dataset", dataset_id="my-dataset")` | `.h5ad` | ✅ | scripted server startup with read-only-backed access |
| 11 | Python | `cellucid.serve_anndata(<data.zarr>, dataset_name="My dataset", dataset_id="my-dataset")` | `.zarr` | ✅ | scripted server startup with eager loading |
| 12 | Jupyter | `cellucid.show(<export_dir>)` | Exported | ✅ | notebook exploration of exports |
| 13 | Jupyter | `cellucid.show_anndata(<data.h5ad>, dataset_name="My dataset", dataset_id="my-dataset")` | `.h5ad` | ✅ | notebook exploration of `.h5ad` |
| 14 | Jupyter | `cellucid.show_anndata(<data.zarr or AnnData>, dataset_name="My dataset", dataset_id="my-dataset")` | `.zarr` / in-memory | ✅ | notebook exploration of eagerly materialized `.zarr` or in-memory data |

\* Browser `.h5ad` loading is typically **not truly lazy** (the browser ends up holding the file in memory).

† Browser Zarr ZIP loading validates and indexes archive metadata before
adoption, then reads gene-expression chunks on demand within the browser's
archive and decoded-chunk memory limits.

## What you should actually use (recommended defaults)

If you don’t have a strong reason otherwise:

- **Local machine, terminal:** #6–#8 (`cellucid serve …`)
- **Notebook-based analysis:** #12–#14 (`show` / `show_anndata`)
- **Sharing with collaborators:** export once (`prepare`) then #2 (public) or #6/#12 (private/local)

## How the CLI maps to the Python API

`cellucid serve <path>` auto-detects format:

- a directory that looks like an export (or contains exported subfolders) → treated as **exported** → served by `serve(...)`
- a file with the exact `.h5ad` suffix → treated as **AnnData** and opened
  read-only-backed by `serve_anndata(...)`
- a complete Zarr v2 directory (both valid `.zgroup` and `.zattrs`, without
  `zarr.json`) → treated as **AnnData** and materialized eagerly by
  `serve_anndata(...)`; the directory name does not need a `.zarr` suffix

Details: {doc}`04_cli_cellucid_serve_quickstart`.

## Where to read next (by option family)

- Options #1–#5 (web app loading): {doc}`../../web_app/b_data_loading/index`
- Options #6–#8 (CLI): {doc}`04_cli_cellucid_serve_quickstart`
- Options #9–#11 (Python servers): {doc}`05_python_serve_and_serve_anndata_quickstart`
- Options #12–#14 (Jupyter): {doc}`06_jupyter_show_and_show_anndata_quickstart`

## Edge cases (high-signal)

- You can *serve a parent directory* containing multiple exported datasets, but your folder must be structured correctly: {doc}`07_exported_directory_mode_show_and_serve`.
- Remote notebooks (JupyterHub, cloud) often require a proxy URL (HTTPS) instead of direct `http://127.0.0.1:<port>`: {doc}`10_notebook_widget_mode_advanced`.

## Troubleshooting

If anything in this table “should work” but doesn’t, start with:
- {doc}`15_troubleshooting_viewing` (symptom → diagnosis → fix)
