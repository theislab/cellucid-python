# Data Loading in the Web App (All Paths)

These pages cover every current route into Cellucid: public samples,
files and folders from your computer, links, local or remote Python and R
servers, Jupyter, and published dataset catalogs. Every route opens the same
way: AnnData, or a file Cellucid presents as AnnData (a Seurat or
SingleCellExperiment `.rds`, Loom, 10x, MatrixMarket, MuData), read in the
browser as each view needs it.

They are written for **mixed audiences**:
- **Wet lab / non-technical**: click-by-click, “what success looks like”, and safe defaults.
- **Computational users**: formats, performance tradeoffs, parameter choices, and reproducibility.
- **Power users / developers**: edge cases, failure modes, and how to debug loading issues.

:::{important}
If you are not sure which method to use, start with {doc}`01_loading_options_overview`.
:::

:::{note}
Vector fields (velocity/drift overlays): if your dataset includes per-cell
vectors, Cellucid can draw them as an animated overlay after loading. Nothing
has to be switched on: however the file is opened — a browser picker, a drop,
a link, a catalog, `cellucid serve`, `serve_anndata()` or `show_anndata()` —
the browser discovers `<field>_<basis>` keys in `obsm` (`velocity_umap` for
`X_umap`, scVelo's naming convention) and reads a field's
vectors only when the overlay shows it.

See {doc}`../i_vector_field_velocity/index` and {doc}`07_file_format_expectations`.
:::

---

## Fast Path (Choose Your Workflow)

Use this as a decision tree. You can always switch later.

| You have… | Best first choice | Why | Next page |
|---|---|---|---|
| A Zarr store folder (`.zarr`), for example from `prepare()` | **Open… › Folder…**, or drop the folder | Read key by key from disk, nothing zipped or uploaded | {doc}`03_browser_file_picker_tutorial` |
| A `.loom`, a 10x Cell Ranger `outs/` folder, or a MuData `.h5mu` | **Open… › File…** (Loom, `.h5mu`) or **Folder…** (`outs/`) | Read in place as it is, no conversion | {doc}`18_open_loom_and_10x` or {doc}`19_multimodal_mudata` |
| A Seurat or SingleCellExperiment object (`.rds`, `.rda`, `.RData`), or a MatrixMarket folder | **Open… › File…** (R objects) or **Folder…** (MatrixMarket) | Cannot be read in parts: read whole once into this browser's storage and kept for next time, no R needed | {doc}`17_open_r_objects` |
| A portable `.zarr.zip` or `.zip` | **Open… › File…** | One validated file selection in every supported browser | {doc}`03_browser_file_picker_tutorial` |
| A `.h5ad` of any size, and no Python | **Open… › File…** | One file selection, nothing installed, nothing uploaded; read from disk as each view needs it | {doc}`03_browser_file_picker_tutorial` |
| A `.h5ad` on a remote machine, or one you drive from a notebook | Server mode or Jupyter | Python hands the file to the browser, which reads it in place as each view needs it | {doc}`04_server_tutorial` or {doc}`05_jupyter_tutorial` |
| A Zarr store folder on a remote machine | Server mode or Jupyter | Served key by key as it is on disk, nothing zipped or copied | {doc}`04_server_tutorial` or {doc}`05_jupyter_tutorial` |
| An in-memory `AnnData` in a notebook | Jupyter | Fastest way to iterate while analyzing | {doc}`05_jupyter_tutorial` |
| One dataset you want to share publicly, of any size | A Zarr store folder on GitHub | Free hosting, no file over GitHub's 100 MB limit, no Git LFS quota; one link opens it | {doc}`15_host_on_github` |
| A dataset collection you want to share publicly | A dataset catalog on GitHub | Exact shareable URLs, no running server | {doc}`11_custom_dataset_repository` |
| A static host, CDN, or intranet server you control | A dataset catalog on that host | `?data=https://host/datasets/` opens your catalog from one shareable link | {doc}`02_local_demo_tutorial` |

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

Every loading path has its own label under **Open data** in the Session panel,
and **Save** and **Load** sit under **Session file**, below it. The **Dataset**
box at the top always describes the dataset currently open, whichever path
opened it ([The Dataset box](01_loading_options_overview.md#the-dataset-box)).
```

:::{note}
The Session panel keeps every loading action and connection state visible.
Its compact `i` buttons only open explanatory help: click one, or focus it with
{kbd}`Tab` and press {kbd}`Enter` or {kbd}`Space`. Press {kbd}`Escape`, click
elsewhere, or move focus outside to close it. See
{doc}`../a_orientation/04_ui_glossary_terminology` and
{doc}`../o_accessibility_privacy_security/01_accessibility`.
:::

(changing-dataset-while-loading)=
### Changing your mind while a dataset loads

Every dataset control is live while a dataset loads, including the dataset
Cellucid opens when it starts. Whatever you choose last is what you get:

- Choosing another sample, **None**, a file or folder with **Open…**, an
  address in the **Link or server** bar, or running the benchmark, stops the
  load in progress at once (its downloads are cancelled) and loads your new
  choice instead. The **Dataset** box shows a spinner and what is opening, for
  example *Opening data.h5ad…*.
- Connecting to a server or GitHub repository that holds exactly one dataset
  opens that dataset, unless you chose another dataset, or applied a session
  with **Load** under **Session file**, after clicking **Open**. Then what you
  did later stays and the connection is kept: its line under the bar still
  names it, with **Disconnect**.
- **Load** chosen while a dataset loads waits for it, says so in a
  notification, and then applies the session to the dataset that loaded. A
  session saved from a different dataset is refused with an explanation, as
  always.

A load you replaced never appears afterwards and reports nothing; only the
outcome of your latest choice is reported. If the dataset Cellucid opens at
startup fails and you had not chosen anything else, Cellucid starts without
it: a panel over the empty viewer says the data could not be opened and why,
with the raw error under **Details**, **Try again** and **Open without this
link**. When asking again cannot change the answer (a file without an
embedding, a folder without `datasets.json`, a file that is not shared),
**Open without this link** comes first.

(drop-to-open)=
### Dropping a file or folder onto the viewer

Instead of a button, you can drag one file (`.h5ad`, `.zarr.zip`, `.rds`,
`.rda`/`.RData`, `.loom`, `.h5`, `.h5mu`), one `.zarr` store folder, one Cell
Ranger `outs/` folder or one MatrixMarket folder from your file manager and
drop it anywhere on the viewer. While you drag, an outline and a *Drop to open* card cover the viewer.
A drop is read exactly as **Open… › File…** or **Open… › Folder…** would
read it, and it replaces any load in progress like any other choice. A
dropped folder's files are asked for several at a time, as many as the
browser answers faster with, and in Chrome and Edge through the folder's
handle, which those browsers list faster than its entries: measured on a
real drop in Chrome, a Zarr folder of 22,023 files is listed in 2.7 s
rather than 4.4 s; in Safari's engine the listing takes about a quarter
less time than one file after another. In Chrome and Edge a dropped file or
folder is kept to reopen with one click, as a chosen one is
({ref}`reopen-local-file`).

Anything else is refused with a notification that says what can be dropped:
several items at once (drop the `.zarr` folder itself, not the files inside
it), any other file, or a folder that holds no Zarr metadata (`.zgroup`,
`.zarray` or `zarr.json`), no Cell Ranger matrix and no MatrixMarket matrix. A file
dropped onto the page never makes the browser leave Cellucid, including while
it is still starting; until the data controls are ready, the pointer shows the
drop as not allowed. The buttons remain the way to open data from the keyboard
or with a screen reader.

---

## Getting Started

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Loading Options Overview
:link: 01_loading_options_overview
:link-type: doc

Every route into Cellucid — the browser itself (a sample, a link, a file) or a Cellucid server — and how to choose between them.
:::

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Publish a Dataset Catalog
:link: 02_local_demo_tutorial
:link-type: doc

Publish a folder of AnnData files with a `datasets.json` catalog on GitHub Pages or any static host, and share it with one link — no server.
:::

:::{grid-item-card} {octicon}`repo;1.5em;sd-mr-1` Custom Dataset Repository
:link: 11_custom_dataset_repository
:link-type: doc

Follow a complete five-dataset reference repository, validate its exact
catalog contract, and publish stable dataset-specific links.
:::

:::{grid-item-card} {octicon}`lock;1.5em;sd-mr-1` Google Drive, OneDrive and Private Sharing
:link: 16_private_drive_onedrive
:link-type: doc

Open a file from Google Drive or OneDrive in place, and share unpublished data
with a collaborator: one sign-in, nothing stored on a server.
:::

:::{grid-item-card} {octicon}`mark-github;1.5em;sd-mr-1` Host a Dataset on GitHub
:link: 15_host_on_github
:link-type: doc

Publish a dataset as a Zarr folder on GitHub (Pages or raw), with the measured
speed of each host, GitHub's limits, and Git LFS's bandwidth trade-off.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Standard Pancreas Sample
:link: 10_standard_pancreas_dataset
:link-type: doc

Load the scVelo pancreas sample and verify its 1D/2D/3D embeddings,
dimension-matched velocity, metadata, genes, connectivity, and provenance.
:::

:::{grid-item-card} {octicon}`verified;1.5em;sd-mr-1` Sample Provenance
:link: 12_sample_dataset_provenance
:link-type: doc

What to cite for each built-in sample, which one can be rebuilt and checked
from its recipe, and what is recorded for the four that cannot.
:::

::::

---

## Viewing Methods

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`upload;1.5em;sd-mr-1` Browser File Picker
:link: 03_browser_file_picker_tutorial
:link-type: doc

Load a file (`.h5ad`, `.zarr.zip`, `.rds`, `.loom`, 10x `.h5`, `.h5mu`) or a
folder (Zarr, Cell Ranger `outs/`, MatrixMarket) with **Open…**, or by
dropping it on the viewer.
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Server Mode
:link: 04_server_tutorial
:link-type: doc

Run a local HTTP server that hands an `.h5ad` or a Zarr store to the viewer, on your machine or over an SSH tunnel.
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` Jupyter Integration
:link: 05_jupyter_tutorial
:link-type: doc

Visualize AnnData objects directly within Jupyter notebooks with embedded interactive widgets.
:::

::::

---

## Concepts & Troubleshooting

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`tag;1.5em;sd-mr-1` Dataset Identity
:link: 06_dataset_identity_why_it_matters
:link-type: doc

What makes a dataset “the same dataset” in Cellucid, and why it matters for sessions, sharing, and community annotation.
:::

:::{grid-item-card} {octicon}`git-compare;1.5em;sd-mr-1` Data Mapping
:link: 13_data_mapping
:link-type: doc

Choose which `obsm` array each dimension and each view
shows, and what happens when a file has no UMAP.
:::

:::{grid-item-card} {octicon}`file-code;1.5em;sd-mr-1` R Objects and MatrixMarket
:link: 17_open_r_objects
:link-type: doc

Open a Seurat or SingleCellExperiment `.rds`, an `.rda` of several objects, or a MatrixMarket count matrix folder.
:::

:::{grid-item-card} {octicon}`file-binary;1.5em;sd-mr-1` Loom and 10x Cell Ranger
:link: 18_open_loom_and_10x
:link-type: doc

Open a `.loom`, a Cell Ranger `.h5` or its `outs/` folder with Cell Ranger's own UMAP and clusters.
:::

:::{grid-item-card} {octicon}`stack;1.5em;sd-mr-1` Multimodal Data (MuData)
:link: 19_multimodal_mudata
:link-type: doc

Open an `.h5mu` or MuData Zarr and switch modality (RNA, protein) under Data mapping.
:::

:::{grid-item-card} {octicon}`mortar-board;1.5em;sd-mr-1` Tutorials: Published Data As It Is
:link: 20_tutorials_open_published_data
:link-type: doc

A GEO sample's MatrixMarket files and a Seurat object on Zenodo by link, and a
Cell Ranger `outs/` folder from disk.
:::

:::{grid-item-card} {octicon}`file;1.5em;sd-mr-1` Format Expectations
:link: 07_file_format_expectations
:link-type: doc

Every format Cellucid reads and how, what `uns['cellucid']` may say, and the `datasets.json` catalog contract.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 08_troubleshooting_data_loading
:link-type: doc

Symptom → diagnosis → fix for file picker, server mode, links and catalogs, and common data issues.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Verified captures
:link: 09_screenshots
:link-type: doc

Current loading controls, a successful direct H5AD load, and browser-engine
acceptance captures.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Vector Fields Overlay
:link: ../i_vector_field_velocity/index
:link-type: doc

If your dataset includes vector fields (e.g. RNA velocity), enable the overlay after loading and pick a field bound to the embedding the view shows.
:::

::::

---

## Quick Comparison

| Method | Best For | Separate server command? |
|--------|----------|--------------------------|
| Public sample | Learning with known-good data | No |
| Browser file picker | Local viewing without Python | No |
| GitHub catalog | Public, dataset-specific share links | No |
| Self-hosted catalog | A catalog on your own static host or CDN | No — any HTTP server will do |
| Server mode | Data on a remote machine, or a Zarr store folder | Yes — `cellucid serve` |
| Jupyter | Interactive analysis workflows | No — the viewer is a notebook widget; no server, no port |

:::{note}
Server mode binds to `127.0.0.1` by default; use the documented SSH-tunnel
workflow for a remote machine. Jupyter uses no server at all: the widget hands
the data to the viewer through the notebook's own connection, so a remote
notebook needs nothing extra.
:::

---

```{toctree}
:maxdepth: 1
:hidden:

01_loading_options_overview
03_browser_file_picker_tutorial
04_server_tutorial
05_jupyter_tutorial
14_open_from_link
16_private_drive_onedrive
13_data_mapping
21_subsampled_loading
17_open_r_objects
18_open_loom_and_10x
19_multimodal_mudata
20_tutorials_open_published_data
10_standard_pancreas_dataset
12_sample_dataset_provenance
15_host_on_github
11_custom_dataset_repository
02_local_demo_tutorial
06_dataset_identity_why_it_matters
07_file_format_expectations
08_troubleshooting_data_loading
09_screenshots
```
