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
| 5 | Cellucid web app | Browser **.zarr folder picker** | `.zarr` | ✅† | quick preview of `.zarr` without Python |
| 6 | Terminal (CLI) | `cellucid serve <export_dir>` | Exported | ✅ | reliable, fast local viewing |
| 7 | Terminal (CLI) | `cellucid serve <data.h5ad>` | `.h5ad` | ✅ | large `.h5ad` with backed mode |
| 8 | Terminal (CLI) | `cellucid serve <data.zarr>` | `.zarr` | ✅ | large `.zarr` with chunked access |
| 9 | Python | `cellucid.serve(<export_dir>)` | Exported | ✅ | scripted server startup |
| 10 | Python | `cellucid.serve_anndata(<data.h5ad>)` | `.h5ad` | ✅ | scripted server startup with backed mode |
| 11 | Python | `cellucid.serve_anndata(<data.zarr>)` | `.zarr` | ✅ | scripted server startup with chunked access |
| 12 | Jupyter | `cellucid.show(<export_dir>)` | Exported | ✅ | notebook exploration of exports |
| 13 | Jupyter | `cellucid.show_anndata(<data.h5ad>)` | `.h5ad` | ✅ | notebook exploration of `.h5ad` |
| 14 | Jupyter | `cellucid.show_anndata(<data.zarr or AnnData>)` | `.zarr` / in-memory | ✅ | notebook exploration of `.zarr` or in-memory |

\* Browser `.h5ad` loading is typically **not truly lazy** (the browser ends up holding the file in memory).

† Browser `.zarr` loading can be “effectively lazy” for gene expression, but still reads metadata up front and can be limited by browser file APIs.

## What you should actually use (recommended defaults)

If you don’t have a strong reason otherwise:

- **Local machine, terminal:** #6–#8 (`cellucid serve …`)
- **Notebook-based analysis:** #12–#14 (`show` / `show_anndata`)
- **Sharing with collaborators:** export once (`prepare`) then #2 (public) or #6/#12 (private/local)

## How the CLI maps to the Python API

`cellucid serve <path>` auto-detects format:

- a directory that looks like an export (or contains exported subfolders) → treated as **exported** → served by `serve(...)`
- a `.h5ad` file → treated as **AnnData** → served by `serve_anndata(..., backed=True)`
- a `.zarr` directory (or directory with `.zattrs`/`.zgroup`) → treated as **AnnData** → served by `serve_anndata(...)`

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
