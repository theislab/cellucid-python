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
- How to open a dataset again with one click from **My datasets**
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
- **Obs fields** (metadata for coloring and filtering): read when a field is first used
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
| AnnData `.h5ad`, AnnData Zarr (folder or `.zarr.zip`), MuData (`.h5mu`, Zarr), Loom, 10x `.h5` and `outs/` | Yes | Only what each view needs, from your disk or with HTTP range requests from a host; what a host sent is kept in this browser, so opening the same version again downloads nothing read before (below); a matrix stored cell by cell has its gene index kept in this browser once built, so the same cells opened again read genes at once ({doc}`21_subsampled_loading`). Served from Python by path, each is served with its cells in priority order, so a subsample sends only its cells |
| A Cellucid-ready AnnData (`prepare()` in Python, `cellucid_prepare()` in R) | Yes, fastest | As above, with the matrix stored by gene and chunked for range reads, so one gene is one small read |
| Seurat or SingleCellExperiment `.rds`, `.rda`, `.RData` (gzip, bzip2, xz, zstd; GEO's `.rds.gz`) | No | Read whole, once, from start to end, and kept in this browser's storage; the next opening of the same file (even days later) opens from there in seconds, reading nothing of the file |
| An uncompressed R file (`saveRDS(compress = FALSE)`) | Its arrays, yes | Its structure is read once and kept; its arrays are read in place |
| MatrixMarket folder (`matrix.mtx[.gz]` with its lists) | No | Its lists and cell tables first: the cells are drawn from them. Then the matrix, read once from its start, keeping only the rows of the cells genes are read for (as many as this device holds a gene index of); both kept ({doc}`17_open_r_objects`). A matrix no larger than its tables, or of several modalities, is read with them |
| Any file on a host that serves no parts of files (no range requests, or `Content-Encoding` on every answer) | No | Downloaded whole, once, into this browser's storage as it arrives (memory where the page cannot write its storage), and kept: the next opening downloads nothing (for a file with a version; below) |
| `.h5ad.gz`, h5Seurat | Not opened | The message says how to write a file Cellucid opens |

Before a file is read whole, the opening says so: that the format cannot be
read in parts, how much will be read, that what is read is kept for next
time, and that a Cellucid-ready AnnData is read in parts instead. A file is
known again by what it is and its version: a file on this computer by its
name, size and modification time; a linked file by the provider's id for it
(for a plain web address, the address without its signature or token
parameters) and its version, asked of the source at every opening, never
remembered from an earlier one: the provider's name for its content where it
has one (a CELLxGENE Census build, a Zenodo version record and checksum, a
Figshare version and checksum, a Hugging Face file's own content hash, the GitHub commit that
last wrote it or its Git LFS object's hash, a Google Drive checksum or
revision, a OneDrive content tag), else its exact size with the
`Last-Modified` date and the strong ETag its host sends (an ETag only when
the host lets other sites read it, `Access-Control-Expose-Headers: ETag`; a
weak `W/` one says nothing of the bytes). A changed file is read again.
A file a viewer host serves (`cellucid serve`, `serve_anndata()`, a
notebook widget, R's `cellucid_view()` / `cellucid_serve()`) has a version
too: the host lists it with the `ETag` its bytes are answered with (made from
the file's size, modification and change times), so a file rewritten in place
is a new version. Only a file with a version is kept for later visits. Of
one known only by its size nothing is used by a later opening, in this tab
or another: a file on a host that sends no `Last-Modified` date and no
readable strong ETag, or a GitHub file whose commit GitHub's API would not
name.

A file read in parts from a host on the network (a link, CELLxGENE, a
provider) keeps what was downloaded of it — its
structure, embeddings, obs columns, the expression of the genes read — in
this browser, under its version, as it arrives: a visit closed the moment
it showed the picture keeps all it read (a browser that crashes loses at
most the piece being written). Opening the same version again, another day
or in another tab, downloads none of it: only its version is asked of the
host (Google Drive and OneDrive read one byte of it beside the request that
finds the file; Figshare reads its first 512 KiB again there), and only what no earlier opening
read is downloaded. A store folder on a web host keeps its files the same
way only when its files are named by their content (Hugging Face; GitHub,
read at the commit that last changed the folder); a folder on a plain
host is read anew at each opening. A file on this machine or served from a
notebook (`cellucid serve`, a widget) keeps no downloaded parts: it costs
no download. Before anything kept is used, every opening
asks the source for the file's version; what is kept of another version is
removed then, and a file that changes while it is open is refused with
"open it again to read the new version". What cannot be noticed, and how
each source names its version:
{doc}`../p_developer_docs/09_data_loading_pipeline_and_caching`, *When a
dataset changes*.

(stored-in-this-browser)=
### Stored in this browser

What is kept stays until you remove it, or until the browser needs the room:
when a new file does not fit in the storage the browser grants this site,
kept gene indexes and downloaded parts are removed first (they are the
cheapest to make again),
then the files opened least recently, never one a tab has open.
**Session › Stored**, a section of its own right below Open data, sums it up
on one line beside **Manage…**, never wrapped: the datasets and what they take
(for example *2 datasets · 540 MB kept*), or with no dataset kept, the one kind
kept (*1 connection speed*) or *3 items*. Hovering the line names every kind,
saved states and remembered connection speeds included (*2 datasets · 540 MB ·
2 saved states · 1 connection speed*), and **Manage…** lists them; its
(i) says what is kept and why. The line does not count **My datasets**,
which has its own control under **Open data** ({ref}`my-datasets`).
**Manage…** opens the *Stored in this browser* dialog, which lists each
dataset something is kept of, by its file's name and host (*this computer*
for a local file; never its link) and the size it all takes, with a row per
kind kept of it: the parts downloaded of a file read in parts (*Downloaded
parts*, with what keeping them saves the next opening), a converted copy (an
R object, a MatrixMarket folder), a whole download (a file its host serves
only whole), or the gene index of a
matrix stored cell by cell (most atlases, every CELLxGENE file) for one
choice of cells, so the same cells opened again read their genes at once
({doc}`21_subsampled_loading`). Each row shows its size, when it was last
used, and whether an open tab uses it. Tick rows and **Remove selected**,
**Remove** one dataset's whole, or **Remove all**; a row an open tab uses
cannot be ticked and stays until that tab closes.

Under the datasets come the other things this browser keeps:

- **Saved states** (*histories of this browser*, with their size): the
  history of each dataset this browser kept as you worked, and the states you
  saved under a name ({doc}`../m_state_history/03_saved_states_and_autosave`).
  A row each, named by its dataset, with its size, what it is (*History kept
  in this browser* or *Saved as "…"*), its steps and when it was last saved;
  ticked and removed like the others, or all at once with the group's
  **Remove**.
- **My datasets**, as one row: how many datasets it lists (*7 datasets opened
  in this browser*) and, in Chrome and Edge, how many of them are files or
  folders of this computer it keeps the handle of (*2 files to reopen*).
  **Show…** opens the My datasets dialog, the one place each dataset is
  listed ({ref}`my-datasets`); **Remove** forgets the whole list with the
  handles of its files (never the files themselves, and never the saved
  states, which are their own rows).
- *Connection speeds remembered (2 hosts)*: the hosts data was read from whose
  speed is kept (in `localStorage`, to size the next opening from each;
  {doc}`21_subsampled_loading`; its tooltip gives each host's round trip and
  rate), with its own **Remove**.
- *Decoder choices measured (4 chunk sizes)*: which decoder (WebAssembly or
  JavaScript, the browser's own or fflate) decodes each compression fastest in
  this browser (in its IndexedDB, so a later visit skips measuring them; no
  dataset is kept), with its own **Remove**.

**Remove all** removes everything no open tab uses: the kept files, every
saved state, My datasets with its handles, the connection speeds and the
decoder choices.

The dialog's foot, **Where it's kept**, answers where all of it lives, in
any browser on any system:

| Row | What it says |
|---|---|
| Place | This browser's private storage for this site: the *origin private file system*, IndexedDB and `localStorage` of the address you opened Cellucid at, for example `www.cellucid.com`. It lives inside your browser profile and is managed by the browser — Chrome and Edge, Firefox and Safari each keep a site's storage in their own profile data — so it is not a folder you can open from your file manager. |
| Holds | *files, saved states, My datasets, connection speeds*: files in its private file system; the saved states, My datasets (in Chrome and Edge with the browser's handles of files to reopen) and the decoder choices in its IndexedDB; the connection speeds in its `localStorage`. All of them are this site's storage in this browser. |
| Removed by | This dialog, or the browser's *Clear site data* (also clearing cookies and site data for the site, or removing the browser profile). |
| Space | What the whole site uses and what its quota still leaves free, marked *browser's figures*: the browser's own estimate (`navigator.storage.estimate()`). Each browser decides how much one site may use, from the disk's size and free space; Cellucid sets no limit of its own. |
| Kept when space is low | *yes* when the browser keeps the site's storage when the disk runs low; *no, the browser may clear it* otherwise (`navigator.storage.persisted()`). **Keep it when space is low** asks the browser (`navigator.storage.persist()`), only when you press it: Firefox asks you, other browsers decide from how the site is used. What a browser clears is read again on the next opening. |

A private or incognito window gives the site no storage that outlives the
window: the Stored row says *Kept for this tab only*, and **Where
it's kept** says the same. Cellucid finds this by asking the browser for its
storage, not by its name. My datasets and the saved states go with the window
too ({ref}`my-datasets`).

In a private window nothing outlives the window:
Firefox and Safari give the page no file storage (WebKit's private browsing has
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
| “I opened it before, in this browser.” | **My datasets** under **Open data** | One click; nothing to find again ({ref}`my-datasets`) |
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
A shared session file (**Share…** under **State › Session**) records which of these routes its dataset
came from when another person could take it too: a sample, a link or a server.
Loading the session then opens that dataset first. See
{doc}`../l_sessions_sharing/03_save_restore_ux`.
:::

(my-datasets)=
## My datasets

**My datasets**, under **Open data** in the Session panel, lists the
datasets you opened in this browser, most recently opened first, so you can
open one again without finding it again. In another browser you can add them
from a library file ({ref}`my-datasets-library-file`).

What is listed:

- a sample you chose under **Sample dataset**;
- a link or a repository id you opened, in **Link or server** or in the
  address that opened the page;
- a Cellucid server (`cellucid serve`, or an address typed in **Link or
  server**);
- a file or folder of this computer, chosen with **Open…** or dropped on the
  page.

What is not listed:

- the sample Cellucid opens by itself when it starts with nothing else to
  open, because you did not choose it (choose it under **Sample dataset** and
  it is listed);
- a dataset a notebook or a Python or R viewer hands to its viewer
  (`show_anndata()`, `cellucid_view()`): it is served to that viewer only and
  cannot be opened again from here.

A dataset is listed once. Opening it again, or opening other cells of it
(**Cells to open**), updates its entry.

### Open a dataset again

The button names the dataset open now when it is listed, and reads *Choose a
dataset…* otherwise. Its menu offers the five most recently opened datasets.
Under each name it says the kind (*Sample*, *Link*, *Server*, *File* or
*Folder*), where it comes from (a link's provider, such as Zenodo, or its
host; a server's host) and when it was last opened. The dataset open now reads
*open now* and is not offered again. **All…**, at the bottom, lists every
dataset. While nothing is listed yet, the menu offers **All…** alone, for its
**Import…**.

Choosing one opens it, in place of the dataset open now:

| Kind | What one click does |
|---|---|
| Sample | Opens it. |
| Link | Opens it. What the link carried to grant access (a key, a signature, a token) is never kept, so a link that needs one may be refused: paste the whole link again in **Link or server**. |
| Server | Opens it while the server runs at that address. Start the server again first (`cellucid serve`). |
| File or folder, in Chrome or Edge | Reads it where it is, through the handle the browser kept of it. The browser may ask you first (*Allow this time*, *Allow on every visit*). A file that moved or changed since is not opened: choose it again ({ref}`reopen-local-file`). |
| File or folder, in Firefox or Safari | These browsers keep no handle, so the menu says *choose it again*. The click opens **Open… › File…** (or **Folder…**), and a notification names the file to choose. |

A dataset opens with the cells you chose for it: the same subsample method,
size and seed. A dataset last open with every cell, or with as many cells as
this device fits, opens as this tab opens datasets now (**Cells to open**,
{doc}`21_subsampled_loading`).

### All…: every dataset

**All…** opens a dialog of every dataset, most recently opened first. Each
row shows the name, the kind and when it was last opened, and one line under
it: format, cells, genes and saved states. Type in **Filter by name…** to
narrow the list. Click a name to see the facts behind it: where it comes from
and how it opens again, its format and size, its cells and genes and the
cells it opens with, when it was first and last opened (and how many times),
its saved states, and what this browser stores of it to open it faster.

Each row has three buttons:

- **Open** opens it again, as the menu does.
- **Continue** opens it, then its latest saved state, so its history is back
  at once and nothing is asked
  ({doc}`../m_state_history/03_saved_states_and_autosave`). It is grayed out
  when this browser keeps no saved state of the dataset, and when a file must
  be chosen again first: choose it, and the offer to continue where you left
  off appears as it opens.
- **Forget** (the bin) asks in the row first, saying what goes and what
  stays. The dataset leaves My datasets and, for a file or folder, this
  browser forgets its handle; the file itself is never touched. Its saved
  states and what is stored of it stay, unless you tick *Also delete its 2
  saved states*. **Keep** leaves everything as it was.

(my-datasets-library-file)=
### Take them to another browser: Export… and Import…

My datasets and the saved states are kept in one browser. To take them to
another browser or another computer, or to keep a copy, use the two buttons
beside **Filter by name…** in **All…**.

**Export…** shows, in place of the rows, every dataset with its saved states
under it, all ticked except saved states kept by another version of Cellucid.
Saved states of datasets no longer listed follow under their own heading,
unticked. Untick what should not go. The line under the list says how many
datasets and saved states are ticked and how much room the saved states take
here. **Export** downloads one `.cellucid-library` file
(*cellucid-library-20261104_103000.cellucid-library*).

Each saved state goes whole: every tracker and branch, every step with its
values and analysis results, notes, **Recently deleted** and the
conversations with agents, with their exact times. Conversations can be
private, so keep the file to yourself. Changes made while recording was paused
and not yet kept as a step cannot go into the file: the export leaves them
out, the history stays paused, and the export says so. Resume recording in
that history first to keep them as a step.

**Import…** asks for a `.cellucid-library` file and adds what it holds
beside what this browser keeps. Nothing kept is replaced:

- a dataset this browser does not list is added; one it lists stays as it is;
- each saved state is kept as a named saved state, *Figure 2 (imported)*, or
  *Autosave of Pancreas (imported)* for an autosave. Importing the same file
  twice keeps two copies.

The line above the list then says what was added and what was already there:
*Added 2 datasets and 3 saved states; 1 dataset was already in My datasets.*
**Continue** and **State › Session › Load…** open an imported
saved state like any other.

A library file dropped on the page, or chosen in **Open… › File…**, is
imported the same way: **All…** opens and its line says what was added. A
session file dropped or chosen there loads as **State › Session › Load… ›
Open a session file…** loads it: its dataset opens first when another is
shown, then its history is added beside yours.

The whole file is read and checked before anything is added, and an import
adds everything or nothing: the saved states are kept in one go, and the
datasets are added only once they are; a damaged or truncated file, a browser
with no room left, or a list that cannot be written adds nothing, and the line
says why. A session file chosen here is refused with *This
file is a Cellucid session file (the history of one dataset), not a library
file. Open it with State › Session › Load… › Open a session file…*, and a library file chosen in
**Load… › Open a session file…** with *This file is a Cellucid library file (datasets and saved states
of My datasets), not a session file. Import it with Session › Open data › My
datasets › All… › Import….*

The file is the session file's container under its own signature
(`CELLUCID_LIBRARY/1`): a manifest of the dataset entries and the saved
states, then each saved state as one complete `.cellucid-session` file
(format 4). What a dataset entry holds is what My datasets keeps of it (below):
a link in its public form, never what it carried to grant access.

### Kept in this browser only

My datasets is kept in this browser's own storage (IndexedDB), never
uploaded: there is no server and no account. It leaves the browser only in a
library file you export. For each dataset it keeps its
name and how it opens again (a link in its public form, never the key,
signature or token it carried; a server's address; a file's or folder's
name, never its path or content), its format, size, cells and genes, the
cells it was last open with, and when it was opened. A file's handle is the
browser's own reference to it, usable only by this site in this browser. The
full list is in {doc}`../o_accessibility_privacy_security/02_privacy_model`.

To clear it, **Forget** datasets one by one in **All…**, or press **Remove**
on the My datasets row under **Session › Stored › Manage…**
({ref}`stored-in-this-browser`). Clearing this site's data in the browser
removes it too.

A private or incognito window keeps My datasets only until the window closes,
like everything else of the site. The line under the menu says when the list
will not outlast the window or the page: *Kept only until this window
closes.*, or *Kept for this page only: this window keeps nothing for this
site.* when the browser gives the page no IndexedDB at all (the list then
goes with a reload).

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
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, My datasets naming the open sample, the Open… button and the Link or server bar; the Stored line; and the collapsed Subsampling and Data mapping sections.
:width: 532px

The Session panel with a sample open: the **Dataset** box, **Open data**,
and the collapsed **Subsampling** and **Data mapping** sections.
```

The browser routes in the tables above are the controls of the **Session**
panel in the sidebar. From the top:

- **Dataset**: the dataset open now ([The Dataset box](#the-dataset-box),
  below).
- **Open data**, four ways to open a dataset. Opening one replaces the
  dataset open now.

  | Control | What it opens |
  |---|---|
  | **Sample dataset** | A public dataset of the page's sample catalog, read in place as each view needs it. While a catalog of several datasets is connected, it lists that catalog's datasets instead. |
  | **My datasets** | A dataset you opened in this browser, again with one click ({ref}`my-datasets`); **All…** also exports and imports them with their saved states ({ref}`my-datasets-library-file`). |
  | **File or folder** | **Open…**, then **File…** or **Folder…**: a file or folder of this computer; nothing is uploaded ({doc}`03_browser_file_picker_tutorial`). Dropping one on the page does the same. |
  | **Link or server** | A link or an id of published data, or the address of a Cellucid server (`localhost:8765`, or just the port), then **Open** ({doc}`14_open_from_link`). |

- **Stored**, once something is kept: one line of what this browser keeps,
  and **Manage…** ({ref}`stored-in-this-browser`).
- **Subsampling** and **Data mapping**, each collapsed to one line that says
  what it holds ({doc}`21_subsampled_loading`, {doc}`13_data_mapping`).

A server started in a terminal or from Python is reached through the **Link
or server** bar, and a notebook needs no control at all because its widget
hands the data over itself. Sharing and loading a session file is in the
sidebar's **State** section, under **Session** (**Share…** and **Load…**;
{doc}`../m_state_history/index`).

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
  bytes arrive. This tab holds what was read (within its memory budget), so
  opening another share of the same file fetches nothing twice. For the next
  visit, this browser also keeps, in this site's storage in your browser
  profile, the downloaded parts of a file read in parts, the gene index of a
  matrix stored cell by cell, and the copy of a file read whole
  ({ref}`which-formats-are-read-in-parts`): only when the source names the
  file's version, and never the parts of data served from this computer or a
  notebook. {ref}`stored-in-this-browser` lists them and removes them, and so
  does clearing this site's data; a private window keeps none of it.
- **Details**, collapsed: the source and citation as links, the file, its
  format and size, how it is read (and, for a format read whole, where it is
  kept, what is read in parts instead and what was not read), cells (with a
  subsample, "228,245 of 2,282,447 (uniform, seed 0)", and ", fitted to this
  device" when the device chose the share), genes, obs fields, the embeddings shown
  (of how many `obsm` arrays), the modality of a MuData dataset, layers, raw,
  the neighbor graph (with a subsample, once it is read, how much of the
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
again shows as *Retrying…* with why, and a request to a host that answers
with errors says so on its line ("… is answering with errors (HTTP 502); trying
again in 1.4 s (attempt 2)…"). **Cancel** stops the opening where it is
and keeps the dataset on screen. How it ended goes to the notification center.

With no dataset open and nothing opening (**None** chosen, a link opened
without, an opening canceled), the same panel over the viewer says **No
dataset open**, where to open one (**Session › Open data**, or a file dropped
on the page) and offers **Show Open data**, which opens the sidebar there.

## Edge Cases (Read This If Something Feels “Weird”)

### Data edge cases
- **No UMAP**: an AnnData without `X_umap_1d`, `X_umap_2d`, `X_umap_3d` or
  `X_umap` opens on a best guess from its other `obsm` arrays
  ({doc}`13_data_mapping`); only a file with nothing drawable is refused.
- **NaN/Inf in embeddings**: a cell whose row holds NaN or an infinity has
  no position in that embedding (one computed on some of the cells, a MuData
  modality that did not measure it). The embedding opens; the cell is not
  drawn, picked or selected in its views, a short notice counts such cells and
  names the first, and every other view, filter and analysis keeps them. In a
  latent space such a cell gets no outlier quantile; in a vector field it has
  no arrow. Only an embedding with no finite row, or a value past float32, is
  refused, naming the cell.
- **Duplicate variable names**: the reader makes repeated gene names unique as
  AnnData's `var_names_make_unique()` does (`MATR3`, `MATR3-1`) and says so.
- **Huge categorical fields** (e.g., 100k unique categories): they load (32-bit
  codes) and the legend scrolls, but their colors are not distinguishable.
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
- AnnData has no drawable `obsm` array: each dense one has another row count,
  places no cell (every row missing or infinite), or holds a value past
  float32 (the message lists each entry and why). An `obsm` with no dense
  numeric array at all (empty, or only tables and sparse arrays) is refused
  sooner as having no cell embedding
  ({ref}`no-cell-embedding`).

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
