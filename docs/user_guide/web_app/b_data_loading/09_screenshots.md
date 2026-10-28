# Verified data-loading captures

Every capture on this page came out of the running web application, at a fixed
1440×1000 viewport, and each crop is taken from the bounding rectangle of the
control it shows rather than from hand-chosen pixels. The scenarios that
reproduce them are in `docs/_tooling/screenshots/scenarios.mjs`.

The files opened here are the public Pancreas sample (scVelo's E15.5
endocrinogenesis data) written in each format, and every link points at a
server on the same computer, so the captures show what the app does with each
format rather than with any one host.

## Loading and session controls

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

The Session panel with the Pancreas sample open. The **Dataset** box at the
top describes the open dataset; **Open data** holds every way to open
another: the **Sample dataset** list, **Open…** for a file or folder, and the
**Link or server** bar. **Save** and **Load** sit under **Session file**, and
**Data mapping** below them.
```

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel-dark.png
:alt: The same Session panel in the dark theme: the Dataset box, Open data with the Sample dataset dropdown, Open… and the Link or server bar, Session file with Save and Load, and the collapsed Data mapping group.
:width: 532px

The same panel in the dark theme.
```

```{figure} ../../../_static/screenshots/data_loading/choose-local-data-control.png
:alt: Close-up of File or folder: the full-width Open… button with its menu open, showing File… (AnnData, Seurat/SCE .rds/.rda, Loom, 10x .h5, .h5mu) and Folder… (Zarr, 10x outs, MatrixMarket), with the pointer on File….
:width: 496px

**Open…** under **File or folder**, with its menu open. **File…** opens a
single-file chooser for every file format Cellucid reads; **Folder…** opens a
folder chooser (a Zarr store, a Cell Ranger `outs/`, a MatrixMarket folder). A
browser's file dialog picks files or a folder, never both, which is why the
button asks first.
```

```{figure} ../../../_static/screenshots/data_loading/connect-remote-server.png
:alt: The Link or server bar holding the address of a cellucid serve server on this computer, with one line under it reading Server, the server's address and a compact Disconnect button.
:width: 496px

The **Link or server** bar connected to a running `cellucid serve`: the bar
holds its address (the port alone, `127.0.0.1:8765`, or the **Local URL**
`http://127.0.0.1:8765` all reach it), and the line under it reads **Server**
with **Disconnect**. The **Viewer URL** does not go here: its `?…` query makes
the bar open it as a link to data. The Viewer URL belongs in the browser
address bar. See {doc}`04_server_tutorial`.
```

```{figure} ../../../_static/screenshots/data_loading/connect-github-catalog.png
:alt: The Link or server bar with the GitHub path "theislab/cellucid-demo-custom-datasets/catalog" typed into it and the pointer on its Open button.
:width: 496px

The **Link or server** bar holding a GitHub catalog folder path; **Open**
reads its `datasets.json` ({doc}`15_host_on_github`).
```

## Opening each format

```{figure} ../../../_static/screenshots/data_loading/h5ad-current-loaded.png
:alt: The whole Cellucid window after the Pancreas sample was opened as an .h5ad file with Open… › File…: the Session panel's Dataset box names the file, its cells and genes and "h5ad file · this computer", and the canvas shows its UMAP.
:width: 1440px

An `.h5ad` opened with **Open… › File…**. It is read in place, so only what a
view needs is read; the Dataset box names the format and *this computer*.
```

```{figure} ../../../_static/screenshots/data_loading/zarr-zip-loaded.png
:alt: The Dataset box after a Zarr ZIP archive was opened from this computer: the dataset's name, its cells and genes, "Zarr v2 ZIP · this computer", its description, and Details open with its source and citation, the file, its format and size, how it is read, and what it holds.
:width: 488px

A Zarr ZIP (`.zarr.zip`) opened from this computer, with **Details** open.
```

```{figure} ../../../_static/screenshots/data_loading/rds-dataset-box.png
:alt: The Dataset box after a Seurat object saved with saveRDS() was opened: its cells and genes, "R object (.rds) · this computer", the line "Read whole once (it cannot be read in parts), kept in this browser's storage for next time", and Details open with Read (whole, once, from start to end, nothing uploaded), Kept (the origin-private file system until it is cleared), Read in parts (a Cellucid-ready AnnData) and what it holds.
:width: 488px

A Seurat object saved with `saveRDS()`. An R file has no random access, so it
is read whole once into this browser's storage and kept for next time, and the
Dataset box says so on its own line ({doc}`17_open_r_objects`).
```

```{figure} ../../../_static/screenshots/data_loading/rda-choose-object.png
:alt: The opening panel over the viewer reading "Opening pancreas-objects.rda" and "pancreas-objects.rda holds 2 single-cell objects. Open which one?", with one button per Seurat object, "pancreas (Seurat, 3,696 cells)" and "pancreas_ductal (Seurat, 916 cells)", and Cancel.
:width: 824px

An `.rda` holding two Seurat objects: the opening asks which one to open.
```

```{figure} ../../../_static/screenshots/data_loading/loom-dataset-box.png
:alt: The Dataset box after a Loom file was opened: the title stored in the file, its cells and genes, "Loom (.loom) · this computer", and Details open with the file, its size, how it is read, and what it holds: four obs fields and a UMAP.
:width: 488px

A Loom file, read in place; its stored title names the dataset
({doc}`18_open_loom_and_10x`).
```

```{figure} ../../../_static/screenshots/data_loading/tenx-outs-window.png
:alt: The whole Cellucid window after a Cell Ranger outs/ folder was opened with Open… › Folder…: the Dataset box reads "10x Cell Ranger outs folder · this computer", and the canvas shows Cell Ranger's UMAP coloured by its clustering, gene_expression_clusters.
:width: 1440px

A Cell Ranger `outs/` folder opened with **Open… › Folder…**: the matrix is
read in place, Cell Ranger's UMAP is the embedding, and its clustering is a
field.
```

```{figure} ../../../_static/screenshots/data_loading/matrix-market-dataset-box.png
:alt: The Dataset box after a MatrixMarket folder with a cell table was opened with Open… › Folder…: the folder's name, its cells and genes, "MatrixMarket folder · this computer", the line saying where the converted data is kept, and Details open.
:width: 488px

A MatrixMarket folder whose cell table holds a UMAP
({doc}`17_open_r_objects`).
```

```{figure} ../../../_static/screenshots/data_loading/data-mapping-modality.png
:alt: The Data mapping group of a CITE-seq .h5mu with Genes open: a Modality select reading "rna · 28 features" above Expression, Gene names, Values and Highly variable genes, and Revert and Apply below.
:width: 488px

A CITE-seq `.h5mu`: **Data mapping › Genes** shows a **Modality** select
({doc}`19_multimodal_mudata`).
```

## The opening itself

```{figure} ../../../_static/screenshots/data_loading/opening-progress.png
:alt: The opening panel centred over the empty viewer while a linked Seurat .rds is downloaded: "Opening" and the file's address, the step "Downloading the whole file (9.4 MB): … does not serve parts of files", a progress bar part-way, how much has arrived with the speed and time left, the note saying the host serves the file only whole, how much is downloaded and that it is kept, and a Cancel button.
:width: 856px

While nothing else is on screen, an opening is shown over the empty viewer:
what is opened, the step under way, how much has arrived, and **Cancel**.
```

```{figure} ../../../_static/screenshots/data_loading/opening-failed.png
:alt: The panel over the empty viewer after a linked file could not be opened: "The data could not be opened", one sentence saying why (the file was not found at its host, HTTP 404) and what to do, a folded Details, and the buttons Open without this link, focused, and Try again.
:width: 856px

Data the page's address names that cannot be opened stays as one panel over
the empty viewer, with the reason, the error under **Details**, **Open
without this link** and **Try again**. The rest of the app works.
```

## What a refusal looks like

Each capture below is the notification Cellucid actually raised for that exact
input, taken from the running application.

```{figure} ../../../_static/screenshots/data_loading/fail-h5ad-wrong-file.png
:alt: A Cellucid notification reading "The selected file is not a valid HDF5 file, so it cannot be an AnnData (.h5ad) file. Save it again with anndata…", with a dismiss button at its right.
:width: 760px

A CSV renamed to `.h5ad`. The leading bytes are checked, not the extension.
```

```{figure} ../../../_static/screenshots/data_loading/fail-remote-server-unreachable.png
:alt: A Cellucid notification saying the Cellucid server could not be reached at a loopback address, and what to check.
:width: 760px

A server address with nothing listening behind it.
```

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: A Cellucid notification saying the linked data could not be opened because the GitHub folder holds no datasets.json, after a repository path rather than its catalog folder was opened.
:width: 760px

A GitHub path pointing at a repository rather than at its catalog folder: the
folder holds no `datasets.json`.
```

```{figure} ../../../_static/screenshots/data_loading/keep-dataset-after-failed-connect.png
:alt: The whole Cellucid window after a server address with nothing behind it was opened: the Pancreas sample is still drawn and still named in the Dataset box, and one error notification sits in the lower right.
:width: 1440px

None of these discard the dataset you already had. After a refusal the previous
dataset is still open and still interactive.
```

## Server and notebook

The server workflow is captured against real data in {doc}`04_server_tutorial`,
with its Viewer URL, the complete `--help` output, and the common startup
failures. The viewer a notebook embeds, and a real `viewer.debug_connection()`
report, are in {doc}`05_jupyter_tutorial`.
