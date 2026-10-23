# Data Loading Overview (Every Route)

This page is the **map** of how Cellucid loads data.

If you are new, start here. If you already know what you want, jump to the page that matches your workflow:

- {doc}`02_local_demo_tutorial`: publish a dataset catalog on GitHub or any static host, with no server
- {doc}`03_browser_file_picker_tutorial`: open a file or folder from your computer
- {doc}`04_server_tutorial`: serve a file from a local or remote Python process
- {doc}`05_jupyter_tutorial`: embed Cellucid in a notebook and interact programmatically
- {doc}`11_custom_dataset_repository`: publish and validate a public multi-dataset
  repository, starting from a complete reference implementation
- {doc}`12_sample_dataset_provenance`: what each built-in sample is, what to cite,
  and how far its build is pinned

## At A Glance

**Audience**
- Wet lab / non-technical: choose a safe workflow and avoid common traps.
- Computational: understand formats, lazy reading, and performance tradeoffs.
- Power users: understand what happens under the hood so you can debug failures.

**Time**
- Fast path decision: ~5 minutes
- Full read (with edge cases): ~20 minutes

**What you’ll learn**
- What Cellucid reads (AnnData, in three containers) and what a dataset catalog is
- What “lazy reading” is and why it matters for large datasets
- How to choose between **Files**, **Open from link**, **Server Mode** and **Jupyter**
- The top failure modes, and how to diagnose them quickly

## What Cellucid Reads

Cellucid reads **AnnData**, and nothing else. The browser has its own AnnData
reader, which reads an AnnData dataset in any of three containers:

1) **An `.h5ad` file** (HDF5).
2) **An AnnData Zarr store folder** (`.zarr`, Zarr v2 or v3).
3) **A zipped Zarr store** (`.zarr.zip`, or `.zip` holding one store).

Any AnnData file written by anndata, Scanpy, scVelo, CELLxGENE, SeuratDisk or
R's anndata opens as it is. `cellucid.prepare()` in Python and
`cellucid_prepare()` in R ({doc}`../../python_package/c_data_preparation_api/index`,
{doc}`../../r_package/index`) write a **Cellucid-ready AnnData**: a standard
AnnData file (anndata reads it) laid out so the browser reads it fast over a
network (the expression matrix gene-major, arrays chunked for range requests)
and carrying its name, description, citation and launch defaults in
`uns['cellucid']`. It is not a separate format: every route below opens it,
and every other tool reads it as AnnData.

A **dataset catalog** lists several AnnData files for one link or for the
**Sample dataset** list: a folder holding `datasets.json` and the files it
names ({doc}`07_file_format_expectations`).

For deeper context:
- why stable ids matter: {doc}`06_dataset_identity_why_it_matters`
- what files and keys are expected: {doc}`07_file_format_expectations`
- when loading fails: {doc}`08_troubleshooting_data_loading`

## What “Lazy Reading” Means (In Plain Language)

When you open a dataset, Cellucid needs several kinds of data:

- **Embeddings** (UMAP coordinates): n_cells × 1 to 3 values, read when a view shows them
- **Obs fields** (metadata for colouring and filtering): read when a field is first used
- **Gene expression** (n_cells × n_genes): huge, read **one gene at a time, on demand**

**Lazy reading** means Cellucid reads only what the current view needs, from
wherever the file is: from your disk, or from a web host with HTTP range
requests. Nothing is copied into memory up front, so there is no file size
limit; the only bound is the memory the arrays actually shown need.

How fast a gene reads depends on how the matrix is stored:
- A **gene-major (CSC)** or dense matrix reads one gene directly. A
  Cellucid-ready AnnData stores its matrix this way.
- A **cell-major (CSR)** matrix is first indexed by gene, once. When the file
  arrives as fast as a disk does (a local file, or a server on the same
  machine or network, as the browser measures while opening), the index is
  built in the background as soon as the file opens; over a slower link, when
  the first gene is shown, so opening downloads only what the first view
  draws.

## Vector fields (velocity/drift overlay): what you need for it to appear

Vector fields are optional per-cell vectors (e.g. RNA velocity) that Cellucid
draws as an animated overlay on an embedding. The browser discovers them in
`obsm` by scVelo's naming convention, `<field>_<basis>` (for example
`velocity_umap`), binds each to the embedding `X_<basis>` it is named for when
both have the same number of columns, and reads a field's vectors only when the
overlay shows it. Nothing has to be switched on.

- Vector fields are **embedding-specific**: one component per column of the embedding they are drawn over.
- They must have the **same row order** as the cells.
- The overlay only shows fields bound to the embedding the view **shows**.

If you expect vector fields but don’t see the overlay toggle or dropdown, see
{doc}`07_file_format_expectations` and
{doc}`../i_vector_field_velocity/index`.

## Fast Path: Choose A Workflow (Decision Tree)

Pick the first row that matches you.

| Your situation | Recommended workflow | Why |
|---|---|---|
| “I just want to look at my data, no Python.” | **Files** (H5AD, Zarr folder or Zarr ZIP) | Nothing installed or uploaded; any file size |
| “The data is already on CELLxGENE, GEO, Zenodo, Hugging Face, Dropbox or any web host.” | **Open from link** | Read in place where it is published ({doc}`14_open_from_link`) |
| “I want to share datasets publicly without running a server.” | A dataset catalog on GitHub or a static host | One shareable link; no running server ({doc}`02_local_demo_tutorial`) |
| “The data is on a remote machine.” | Server Mode over an SSH tunnel | The server hands the file to the browser, which reads only what each view needs |
| “I’m already working in a notebook.” | Jupyter (`show_anndata`) | Tight analysis loop; programmatic control |
| “I want the fastest possible web experience.” | Write a Cellucid-ready AnnData with `prepare()` | Gene-major matrix and chunking made for range reads |

## Every Route Into Cellucid

Every dataset reaches the viewer in one of two ways: **the browser reads it
itself**, from your computer or from a web address, or **a Cellucid server
hands it over** (one you start in a terminal, from Python, or that a notebook
starts for you). Either way, the same browser AnnData reader reads it.

### The browser reads it (nothing to install)

| Route | How you point Cellucid at the data | What it opens |
|---|---|---|
| **Sample dataset** | Pick from the **Sample dataset** list: Cellucid's own catalog on www.cellucid.com | The catalog's AnnData files |
| **Open from link** | Paste a link or id into **Open from link**, or put it in `?data=...` ({doc}`14_open_from_link`) | An `.h5ad`, a Zarr ZIP, a Zarr store folder, or a dataset catalog folder |
| **Files** | The **H5AD**, **Zarr ZIP** and **Zarr folder** buttons under **Files**, or a file or folder dropped on the viewer | An `.h5ad`, a `.zarr.zip`, or a `.zarr` folder |

### A Cellucid server hands it over

The server runs where the data is; the web app connects to it, either because
the server opened it (`cellucid serve` opens your browser, a notebook embeds
it) or through the **Server** field under **Open data**. The server computes
nothing: it serves the file's bytes as they are on disk, and the browser reads
them with range requests, exactly as it reads a linked file (an in-memory
AnnData is first written once to a temporary `.h5ad`).

| Route | Command | Best for |
|---|---|---|
| **Terminal** | `cellucid serve data.h5ad` (or `data.zarr`, `data.zarr.zip`) | A dataset on a remote machine (over an SSH tunnel) |
| **Python** | `cellucid.serve_anndata("data.h5ad")` | Scripting server startup |
| **Jupyter** | `cellucid.show_anndata(adata_or_path)` | Exploring while you analyze |

The dataset's id and name come from the file's `uns['cellucid']`; a file that
does not state them needs `--dataset-id` and `--dataset-name`
(`dataset_id=`, `dataset_name=` in Python).

:::{note}
The **Sample dataset** list and a linked catalog folder read the same thing:
a folder holding `datasets.json` and the AnnData files it lists. The sample
catalog is the one the page is deployed with (Cellucid's own on
www.cellucid.com); a linked one is any catalog you point at, on GitHub, a
static host, a CDN or an intranet server, and it behaves identically: it opens
on the dataset its catalog names as default, and applies the default session
its catalog names for that dataset. Worked example:
{doc}`11_custom_dataset_repository`.
:::

:::{tip}
A saved session (**Save** under **Session file**) records which of these routes its dataset
came from when another person could take it too: a sample, a link or a server.
Loading the session then opens that dataset first. See
{doc}`../l_sessions_sharing/03_save_restore_ux`.
:::

## Minimal Commands (Copy/Paste)

### Install

```bash
pip install cellucid
```

### 1) Write a Cellucid-ready AnnData (optional, fastest)

```python
import cellucid

cellucid.prepare(
    "pbmc.zarr",
    format="zarr",                     # "zarr" (folder), "zarr.zip" or "h5ad"
    embeddings={"X_umap": adata.obsm["X_umap"]},
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    dataset_id="pbmc-demo-v1",
    dataset_name="PBMC demo",
)
```

Open `pbmc.zarr` with the **Zarr folder** button, serve it, or publish it in a
catalog. See {doc}`../../python_package/c_data_preparation_api/index` for every
argument.

### 2) Serve a file

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My study" --dataset-id my-study-v1
```

A store written by `prepare()` states its own id and name in
`uns['cellucid']`, so `cellucid serve /path/to/pbmc.zarr` needs neither flag.

### 3) Notebook quick start

```python
from cellucid import show_anndata

viewer = show_anndata("data.h5ad", dataset_name="My study", dataset_id="my-study-v1")
```

## Interface reference

For the complete verified loading gallery, see {doc}`09_screenshots`.

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Cellucid Session panel with a dataset loaded, showing the Dataset summary table; the Open data group with the Sample dataset dropdown, the file buttons under Files, the Server and Open from link fields with their Connect buttons, and the Transfers list; and the Session file group with its Save and Load buttons.
:width: 524px

The browser routes in the table above are the controls in this panel; a
server started in a terminal or from Python is reached through the
**Server** field, and a notebook needs no control at all because it supplies the
address.
```

## Edge Cases (Read This If Something Feels “Weird”)

### Data edge cases
- **No UMAP**: an AnnData without `X_umap_1d`, `X_umap_2d`, `X_umap_3d` or
  `X_umap` opens on a best guess from its other `obsm` arrays
  ({doc}`13_data_mapping`); only a file with nothing drawable is refused.
- **NaN/Inf in embeddings**: an embedding holding one is refused, naming the
  cell, before it can replace what is shown.
- **Duplicate variable names**: the reader makes repeated gene names unique as
  AnnData's `var_names_make_unique()` does (`MATR3`, `MATR3-1`) and says so.
- **Huge categorical fields** (e.g., 100k unique categories): they load (32-bit
  codes) and the legend scrolls, but their colours are not distinguishable.
- **A malformed `uns['cellucid']`**: an unknown key or a wrong value refuses the
  file with a message naming it ({doc}`07_file_format_expectations`).

### Scale edge cases
- **Very large `.h5ad` with a cell-major (CSR) matrix**: opens at any size; the
  first genes wait for the background index, which reports its progress. A
  gene-major (CSC) matrix needs no index.
- **Millions of cells + dense gene matrix**: even with lazy reading, gene queries can be slow or memory-heavy.

### Environment edge cases
- **Zarr folder picker**: it requires the browser's folder-input capability,
  which every supported browser has. H5AD and Zarr ZIP use ordinary
  single-file inputs.
- **Corporate environments**: content blockers can break `raw.githubusercontent.com` (GitHub-hosted catalogs).

## Troubleshooting (Quick Index)

### Symptom: “Nothing loads / it spins forever”
**Likely causes**
- CORS/network blocked (a link, server mode or GitHub)
- The first gene of a large cell-major matrix is waiting for its index
- Browser GPU/WebGL context lost due to memory pressure

**How to confirm**
- Open Developer Tools → Console → look for network errors / memory warnings
- Open the **Transfers** list under **Open data**: it shows what is being read
- Try a built-in sample: if samples load, your environment is probably fine

**Fix**
- Write the file again with `prepare()` (gene-major matrix), or serve it with
  `cellucid serve`

---

### Symptom: “It says Cellucid found nothing to plot”
**Likely causes**
- AnnData has no drawable `obsm` array: only tables, sparse arrays, or arrays
  holding missing or infinite values (the message lists each entry and why).

**Fix**
- Compute an embedding (for example `sc.tl.umap(adata)`) before viewing, or
  fix the array the message names. See {doc}`13_data_mapping`.

---

### Symptom: “Open from link says the folder holds no datasets.json”
**Likely causes**
- The link points at the wrong folder, or at a page instead of the folder
- The folder holds AnnData files but no `datasets.json`

**Fix**
- Link to one file or Zarr store folder directly, or add a `datasets.json`
  listing the files ({doc}`02_local_demo_tutorial`).

## Next Steps

- Want to share datasets publicly (no server)? → {doc}`02_local_demo_tutorial`
- Want a complete public repository to inspect and adapt? → {doc}`11_custom_dataset_repository`
- Want to open a file from your own machine right now? → {doc}`03_browser_file_picker_tutorial`
- Working with data on another machine? → {doc}`04_server_tutorial`
- Already in notebooks? → {doc}`05_jupyter_tutorial`
