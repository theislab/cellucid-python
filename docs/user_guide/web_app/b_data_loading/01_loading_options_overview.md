# Data Loading Overview (Every Route)

This page is the **map** of how Cellucid loads data.

If you are new, start here. If you already know what you want, jump to the page that matches your workflow:

- {doc}`02_local_demo_tutorial`: publish a dataset catalog on GitHub or any static host, with no server
- {doc}`03_browser_file_picker_tutorial`: open a file or folder from your computer
- {doc}`04_server_tutorial`: serve a file from a local or remote Python process
- {doc}`05_jupyter_tutorial`: show Cellucid in a notebook cell and interact programmatically
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
- What Cellucid reads (AnnData, and the formats it reads as AnnData) and what a dataset catalog is
- What “lazy reading” is and why it matters for large datasets
- How to choose between a file or folder (**Open…**), a link or a server (the
  **Link or server** bar), and **Jupyter**
- The top failure modes, and how to diagnose them quickly

## What Cellucid Reads

Cellucid reads **AnnData**. The browser has its own AnnData reader, which
reads an AnnData dataset in any of three containers:

1) **An `.h5ad` file** (HDF5).
2) **An AnnData Zarr store folder** (`.zarr`, Zarr v2 or v3).
3) **A zipped Zarr store** (`.zarr.zip`, or `.zip` holding one store).

It reads the other common formats as AnnData, as they are: MuData (`.h5mu`
or MuData Zarr), Loom, a 10x Cell Ranger `.h5` or `outs/` folder, a Seurat or
SingleCellExperiment `.rds`/`.rda`, and a MatrixMarket folder. The full list,
and how each one is read, is in
[What Cellucid reads](07_file_format_expectations.md#what-cellucid-reads).

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
limit; the only bound is the memory the arrays actually shown need. An R
object (`.rds`, `.rda`) and a MatrixMarket folder are the exception: they
cannot be read in parts, so they are read whole, once, into this browser's
storage, kept there (for your next visit, when the file has a version), and read from there as views need
them ([below](which-formats-are-read-in-parts)).

How fast a gene reads depends on how the matrix is stored:
- A **gene-major (CSC)** or dense matrix reads one gene directly. A
  Cellucid-ready AnnData stores its matrix this way.
- A **cell-major (CSR)** matrix is first indexed by gene, once. When the file
  arrives as fast as a disk does (a local file, or a server on the same
  machine or network, as the browser measures while opening), the index is
  built in the background as soon as the file opens; over a slower link, when
  the first gene is shown, so opening downloads only what the first view
  draws.

(which-formats-are-read-in-parts)=
## Which formats are read in parts

| Format | Read in parts? | How it opens |
|---|---|---|
| AnnData `.h5ad`, AnnData Zarr (folder or `.zarr.zip`), MuData (`.h5mu`, Zarr), Loom, 10x `.h5` and `outs/` | Yes | Only what each view needs, from your disk or with HTTP range requests from a host; a matrix stored cell by cell has its gene index kept in this browser once built, so the same cells opened again read genes at once ({doc}`21_subsampled_loading`). Served from Python by path, each is served with its cells in priority order, so a subsample sends only its cells |
| A Cellucid-ready AnnData (`prepare()` in Python, `cellucid_prepare()` in R) | Yes, fastest | As above, with the matrix stored by gene and chunked for range reads, so one gene is one small read |
| Seurat or SingleCellExperiment `.rds`, `.rda`, `.RData` (gzip, bzip2, xz; GEO's `.rds.gz`) | No | Read whole, once, from start to end, and kept in this browser's storage; the next opening of the same file (even days later) opens from there in seconds, reading nothing of the file |
| An uncompressed R file (`saveRDS(compress = FALSE)`) | Its arrays, yes | Its structure is read once and kept; its arrays are read in place |
| MatrixMarket folder (`matrix.mtx[.gz]` with its lists) | No | As an R object: read whole once, kept |
| Any file on a host that serves no parts of files (no range requests, or `Content-Encoding` on every answer) | No | Downloaded whole, once, into this browser's storage as it arrives (memory where the page cannot write its storage), and kept: the next opening downloads nothing (for a file with a version; below) |
| `.h5ad.gz`, zstd-compressed R files (`compress = "zstd"`), h5Seurat | Not opened | The message says how to write a file Cellucid opens |

Before a file is read whole, the opening says so: that the format cannot be
read in parts, how much will be read, that what is read is kept for next
time, and that a Cellucid-ready AnnData is read in parts instead. A file is
known again by what it is and its version: a file on this computer by its
name, size and modification time; a linked file by the provider's id for it
(for a plain web address, the address without its signature or token
parameters) and the version its host reports: an ETag the host lets other
sites read (`Access-Control-Expose-Headers: ETag`), or the provider's
revision or checksum (CELLxGENE, Zenodo, Figshare, Hugging Face, GitHub,
Google Drive, OneDrive, when they state one). A changed file is read again.
A file a viewer host serves (`cellucid serve`, `serve_anndata()`, a
notebook widget, R's `cellucid_view()` / `cellucid_serve()`) has a version
too: the host lists it with the `ETag` its bytes are answered with (made from
the file's size, modification and change times), so a file rewritten in place
is a new version. Only a file with a version is kept for later visits. One
known only by its size is kept for its tab only, and only that tab opens it
again: a file on a host that reports no readable version, a GEO file whose
size comes from its series' file list.

What is kept stays until you clear it, or until the browser needs the room:
when a new file does not fit in the storage the browser grants this site,
kept gene indexes are removed first (they are the cheapest to make again),
then the files opened least recently, never one a tab has open.
**Session › Open data › Stored** says how many files are kept and how much
they take, and of how many files stored cell by cell (most atlases, every
CELLxGENE file) the gene index is kept, so the same cells opened again read
their genes at once ({doc}`21_subsampled_loading`), with **Clear**, which
removes whatever no open tab uses. In a private window nothing outlives the window:
Firefox and Safari give the page no storage (WebKit's private browsing has
neither the origin-private file system nor IndexedDB, measured in Playwright's
WebKit, whose test windows are private ones), so the file is held in memory
for the tab; Chrome gives it a small storage held in memory, so a large R object there
stops with *This browser's storage for this site is full*: open it in a
normal window.

To share a large R object or MatrixMarket matrix, write it once as a
Cellucid-ready AnnData: in R, `cellucid::cellucid_prepare(path, data = obj)`;
in Python, `cellucid.prepare()`. Every viewer then reads only what each view
needs, with nothing read whole.

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
| “I just want to look at my data, no Python.” | **Open…** under **File or folder** (an `.h5ad`, a `.zarr` folder or a `.zarr.zip`) | Nothing installed or uploaded; any file size |
| “The data is already on CELLxGENE, GEO, Zenodo, Figshare, Hugging Face, Dropbox, a public Google Drive file or any web host.” | The **Link or server** address bar | Read in place where it is published ({doc}`14_open_from_link`) |
| “My collaborator shared a file with me on Google Drive or OneDrive.” | The **Link or server** address bar, then one sign-in | Read in place from Google or Microsoft with your own sign-in ({doc}`16_private_drive_onedrive`) |
| “I want to share datasets publicly without running a server.” | A dataset catalog on GitHub or a static host | One shareable link; no running server ({doc}`02_local_demo_tutorial`) |
| “The data is on a remote machine.” | Server Mode over an SSH tunnel | The server hands the file to the browser, which reads only what each view needs |
| “I’m already working in a notebook.” | Jupyter (`show_anndata`) | Tight analysis loop; programmatic control |
| “I want the fastest possible web experience.” | Write a Cellucid-ready AnnData with `prepare()` | Gene-major matrix and chunking made for range reads |

## Every Route Into Cellucid

Every dataset reaches the viewer in one of three ways: **the browser reads it
itself**, from your computer or from a web address; **a Cellucid server hands
it over** (one you start in a terminal, from Python or R); or **a notebook
hands it over** through the widget it shows the viewer in. Every way, the same
browser AnnData reader reads it.

### The browser reads it (nothing to install)

| Route | How you point Cellucid at the data | What it opens |
|---|---|---|
| **Sample dataset** | Pick from the **Sample dataset** list: Cellucid's own catalog on www.cellucid.com | The catalog's AnnData files |
| **Link** | Paste a link or id into the **Link or server** bar and click **Open**, or put it in `?data=...` ({doc}`14_open_from_link`) | An `.h5ad`, a Zarr ZIP, a Zarr store folder, any other file Cellucid reads (an `.rds`, a `.loom`, an `.h5mu`, a GEO MatrixMarket sample), or a dataset catalog folder |
| **File or folder** | **Open…** under **File or folder**, then **File…** (an `.h5ad`, `.zarr.zip`, `.rds`, `.loom`, 10x `.h5` or `.h5mu`) or **Folder…** (a `.zarr` store, a Cell Ranger `outs/`, a MatrixMarket folder); or a file or folder dropped anywhere on the page | Every format of {doc}`07_file_format_expectations` |

### A Cellucid server hands it over

The server runs where the data is; the web app connects to it, either because
the server opened it (`cellucid serve` opens your browser) or through the
**Link or server** bar under **Open data** (type its
address, for example `localhost:8765`, or just the port). The browser reads
the data with range requests, exactly as it reads a linked file: a file
`prepare()` wrote is served as its bytes on disk; an in-memory object (Python
or R), and an AnnData file in any other order (Python), are served as a Zarr
store in priority order built as the browser asks for it, so a subsample
transfers only its cells.

| Route | Command | Best for |
|---|---|---|
| **Terminal** | `cellucid serve data.h5ad` (or `data.zarr`, `data.zarr.zip`, any other format Cellucid reads) | A dataset on a remote machine (over an SSH tunnel) |
| **Python** | `cellucid.serve_anndata("data.h5ad")` | Scripting server startup |
| **R** | `cellucid::cellucid_view(obj)` (a Seurat, SingleCellExperiment or SummarizedExperiment object, or a path) in RStudio, Positron or an R console | Exploring from an R session ({doc}`../../r_package/index`) |

### A notebook hands it over

In a Jupyter notebook (JupyterLab, Notebook 7, VS Code, Colab, marimo,
Databricks), the viewer is a widget in the cell's output: it frames the app
from www.cellucid.com and hands it the data from the kernel through the
notebook's own connection. No server, no port, no address to type.

| Route | Command | Best for |
|---|---|---|
| **Python** | `cellucid.show_anndata(adata_or_path)` | Exploring while you analyze ({doc}`05_jupyter_tutorial`) |
| **R** (IRkernel) | `cellucid::cellucid_view(obj)` | The same from R ({doc}`../../r_package/d_viewing_loading/08_view_in_jupyter_notebooks`) |

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

Open `pbmc.zarr` with **Open…** › **Folder…**, serve it, or publish it in a
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
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

The browser routes in the table above are the controls in this panel; a
server started in a terminal or from Python is reached through the **Link or
server** bar, and a notebook needs no control at all because its widget
hands the data over itself.
```

## The Dataset box

The **Dataset** box at the top of the **Session** panel describes the dataset
that is open, and only that one; opening another rewrites all of it.

- **The name**, then one line with its cells and genes; with a share of the
  cells open, the cells read "228K of 2.3M" ({doc}`21_subsampled_loading`).
- **Where it is read from:** the format (`h5ad file`, `Zarr v2 folder`,
  `Zarr v3 ZIP`, `Loom (.loom)`, `R object (.rds)`, …) and the host
  (`theislab.github.io`, `localhost:8765`), `this computer` for a file you
  opened, or the program that opened the viewer (`Python`, `Jupyter`, `R`).
- **For an R object or a MatrixMarket folder**, one more line: *Read whole
  once (it cannot be read in parts), kept in this browser's storage for next
  time*, or *Opened from this browser's storage, read whole on an earlier
  open* (or *into memory for this tab*, when the browser gives the page no
  storage, as in a private window), and how many parts of the file were not
  read.
- **The description** the dataset states (`uns['cellucid']['description']`),
  two lines with **More**.
- **How much has been read**, for data read in place from a host:
  "5.5 MB read in 43 requests", or "Reading… 5.5 MB read in 43 requests" while
  bytes arrive. What was read stays in this tab's memory (within its memory
  budget), so opening another share of the same file fetches nothing twice,
  and is released with the tab. Only what is listed under
  {ref}`which-formats-are-read-in-parts` is kept on disk: a file read whole,
  and the gene index of a matrix stored cell by cell.
- **Details**, collapsed: the source and citation as links, the file, its
  format and size, how it is read (and, for a format read whole, where it is
  kept, what is read in parts instead and what was not read), cells (with a
  subsample, "228,245 of 2,282,447 (uniform, seed 0)", and ", fitted to this
  device" when the device chose the share), genes, obs fields, the embeddings shown
  (of how many `obsm` arrays), the modality of a MuData dataset, layers, raw,
  the neighbour graph (with a subsample, once it is read, how much of the
  open cells' graph it keeps), the latent space and the vector fields. A change in
  **Data mapping** updates them.

While a dataset opens, the box shows a spinner and what is opening. With no
dataset open it says so. Progress of long reads and anything that went wrong
are notifications, not text in the panel.

Every opening (a link, a file, a sample, a catalog dataset, a server or a
notebook) is shown from its first moment: over the empty viewer while no
dataset is shown, or as a notification with **Cancel** while another dataset
stays on screen. It names what is opened, the step under way, how much has
arrived, how fast and how long is left when that is known; a request made
again shows as *Retrying…* with why. **Cancel** stops the opening where it is
and keeps the dataset on screen. How it ended goes to the notification center.

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
- **Folder picker**: **Open… › Folder…** uses the browser's folder-input
  capability, which every supported browser has (Safari's leaves out the
  hidden files of a Zarr v2 store: drop the folder instead, see
  {doc}`03_browser_file_picker_tutorial`). **Open… › File…** is an ordinary
  single-file input for every file format Cellucid reads.
- **Corporate environments**: content blockers can break `raw.githubusercontent.com` (GitHub-hosted catalogs).

## Troubleshooting (Quick Index)

### Symptom: “Nothing loads / it spins forever”
**Likely causes**
- CORS/network blocked (a link, server mode or GitHub)
- The first gene of a large cell-major matrix is waiting for its index
- Browser GPU/WebGL context lost due to memory pressure

**How to confirm**
- Open Developer Tools → Console → look for network errors / memory warnings
- Look at the **Dataset** box: for linked data it shows how much has been
  read ("Reading…" while bytes arrive), and a retried request is a notification
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

### Symptom: “A link says the folder holds no datasets.json”
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
