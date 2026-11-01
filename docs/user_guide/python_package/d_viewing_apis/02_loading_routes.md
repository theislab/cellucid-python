# Every loading route

This page lists every way data reaches the Cellucid viewer, grouped by **who
reads the data**: the browser itself, or a Cellucid server that Python starts.
Every route reads with the same reader in the browser, and every route takes
the formats of the input format registry ({doc}`../../web_app/b_data_loading/07_file_format_expectations`):
AnnData (`.h5ad`, a Zarr store folder or `.zip`), MuData (`.h5mu`, MuData
Zarr), Loom, 10x Cell Ranger `.h5` files and `outs/` folders, MatrixMarket
folders, and R objects (`.rds`, `.rda`/`.RData`, gzipped or not); through
Python also an in-memory `AnnData` or `MuData`.

If you want a decision tree instead of a table, go to {doc}`03_choose_your_workflow_decision_tree`.

## Definitions

- **AnnData**: `.h5ad`, a Zarr store (folder or `.zip`), or an in-memory `AnnData` object.
- **MuData**: `.h5mu`, a MuData Zarr store, or an in-memory `MuData` object
  (several modalities over the same cells; needs `pip install "cellucid[mudata]"`
  only to write one).
- **Cellucid-ready AnnData**: an AnnData written by `cellucid.prepare(...)`
  (or R's `cellucid_prepare()`): gene-major, chunked for range reads, and
  naming itself in `uns['cellucid']`
  ({doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`).
- **Lazy genes**: gene expression is fetched **on demand**, gene by gene,
  instead of loading the whole matrix. Every route below reads lazily.

## The browser reads it (no Python)

| Route | How you point Cellucid at the data | Formats | Best for |
|---|---|---|---|
| **Sample dataset** | the **Sample dataset** list, Cellucid's own catalog | Cellucid-ready AnnData | learning the UI with known-good data |
| **Link** | a link or id in the **Link or server** bar, or `?data=…` ({doc}`../../web_app/b_data_loading/14_open_from_link`) | every registry format (GEO MatrixMarket samples and `.rds.gz` included), a dataset catalog folder | opening data where it is already published; sharing it as one link |
| **File or folder** | **Open…**, then **File…** or **Folder…**, or a file or folder dropped on the page | every registry format | viewing a file on your computer without Python |

## A Cellucid server serves it

| Route | Command | Formats | Best for |
|---|---|---|---|
| **Terminal** | `cellucid serve <file or folder> --dataset-name "My dataset" --dataset-id my-dataset` | every registry format | reliable viewing of large data; a remote machine over an SSH tunnel |
| **Python** | `cellucid.serve_anndata(<a path, an AnnData or a MuData>, dataset_name="My dataset", dataset_id="my-dataset")` | every registry format, in-memory AnnData and MuData | scripted server startup |
| **Jupyter** | `cellucid.show_anndata(<a path, an AnnData or a MuData>, dataset_name="My dataset", dataset_id="my-dataset")` | every registry format, in-memory AnnData and MuData | notebook exploration |

A Cellucid-ready AnnData or MuData declares its own name and id. Any other
dataset is named after its file (its id the file's name without the
extension), and the id made so is printed with how to set it; the identity
arguments name it otherwise. Sessions are kept under the id.

The web app reaches a server through its **Link or server** bar (type
`localhost:8765`, or just the port), or the server opens it: `cellucid serve`
opens your browser. In Jupyter, `show_anndata()` starts no server: a widget in
the cell's output frames the app and hands it the bytes through the notebook's
own connection. The viewer reads the data itself exactly as it reads a linked
file — in place with range requests where the format allows it, once from the
start for an R object or a MatrixMarket folder; Python serves the bytes (an in-memory AnnData,
and an AnnData file not in priority order, as a store in priority order made as the viewer reads
it; a MuData object the same way, as a MuData store).

## What you should actually use (recommended defaults)

- **Local machine, terminal:** `cellucid serve …`
- **Notebook-based analysis:** `show_anndata`
- **Sharing with collaborators:** write a Cellucid-ready AnnData once
  (`prepare`), then publish the file anywhere a link can be read from, or
  serve it

## How the CLI maps to the Python API

`cellucid serve <path>` detects the format by the input format registry the
web app reads too, and serves it with `serve_anndata(...)`:

- a file whose name a format names (`.h5ad`, `.zip`/`.zarr.zip`, `.h5mu`,
  `.loom`, `.h5`, `.rds[.gz]`, `.rda[.gz]`, `.RData[.gz]`) → served as one
  file; its first bytes and HDF5 root decide what it is (an HDF5 root marked
  `mudata` is a MuData, whatever its name)
- a Zarr store folder (`zarr.json` for Zarr v3, `.zgroup` for Zarr v2, at its
  root) → served key by key, as AnnData or (by its root's `encoding-type`)
  MuData; the directory name does not need a `.zarr` suffix
- a folder of files: a Cell Ranger `outs/` (its
  `filtered_feature_bc_matrix.h5`), or a MatrixMarket folder (a `.mtx[.gz]`
  file, whatever its name, such as GEO's `GSM4567_matrix.mtx.gz`, with its
  `barcodes`/`features` or `genes` lists) → served file by file, the files listed in
  `/_cellucid/datasets`

Anything else is refused: `Unable to detect format for: <path>`.

Details: {doc}`04_cli_cellucid_serve_quickstart`.

## Where to read next (by route)

- The browser routes: {doc}`../../web_app/b_data_loading/index`
- Terminal: {doc}`04_cli_cellucid_serve_quickstart`
- Python servers: {doc}`05_python_serve_anndata_quickstart`
- Jupyter: {doc}`06_jupyter_show_anndata_quickstart`

## Edge cases (high-signal)

- Remote notebooks (JupyterHub, VS Code Remote, Colab) need nothing extra with
  `show_anndata()`: {doc}`10_notebook_widget_mode_advanced`. A browser-tab
  server on a remote machine needs an SSH tunnel:
  {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

## Troubleshooting

If anything in this table “should work” but doesn’t, start with:
- {doc}`15_troubleshooting_viewing` (symptom → diagnosis → fix)
