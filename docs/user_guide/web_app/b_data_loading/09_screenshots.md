# Verified data-loading captures

Every capture on this page came out of the running web application, at a fixed
1440×1000 viewport, and each crop is taken from the bounding rectangle of the
control it shows rather than from hand-chosen pixels.

The browser-engine acceptance captures used a deterministic synthetic AnnData
fixture. That acceptance path replaced the prepared dataset with its Zarr ZIP
representation, replaced that with the H5AD representation, and then restored
the prepared dataset through the visible controls. Each replacement asserted the
120-cell, 6-gene dataset shape, 2D Planar navigation, hidden camera transport,
settled notifications, and the absence of browser-console, page, and HTTP
errors.

The control and refusal captures used the published Pancreas sample as the
starting dataset, so the panel state around them is a real loaded dataset rather
than an empty application.

## Loading and session controls

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreas sample loaded: a dataset summary table, a Sample datasets dropdown, H5AD, Zarr ZIP and Prepared buttons, Remote server and GitHub data fields each with a Connect button, and Save State and Load State buttons.
:width: 524px

The Session panel in its settled state, with the standard Pancreas sample
loaded. Every loading path has its own labelled group, and Save State and Load
State stay beside the dataset controls. Both connection buttons read
**Connect** once the app is ready.
```

```{figure} ../../../_static/screenshots/data_loading/choose-local-data-control.png
:alt: Close-up of the Local data controls with a mouse pointer on the full-width Prepared button below the H5AD and Zarr ZIP buttons.
:width: 496px

The three local-data controls. **Prepared** opens a directory chooser; **H5AD**
and **Zarr ZIP** open single-file choosers.
```

```{figure} ../../../_static/screenshots/data_loading/connect-remote-server.png
:alt: Close-up of the Remote server control with a loopback URL typed into the field and a mouse pointer on the Connect button.
:width: 496px

The **Remote server** control, holding the address a running `cellucid serve`
printed.
```

```{figure} ../../../_static/screenshots/data_loading/connect-github-catalog.png
:alt: Close-up of the GitHub data control with an owner/repo/exports path typed into the field and a mouse pointer on the Connect button.
:width: 496px

The **GitHub data** control, holding an exports-root path.
```

## Direct H5AD load

```{figure} ../../../_static/screenshots/data_loading/h5ad-current-loaded.png
:alt: A current-schema H5AD file loaded in Cellucid with dataset metadata visible.
:width: 1280px

A current-schema H5AD file loaded through the browser picker, with its dataset
metadata and available dimensions reported in the Session panel.
```

```{figure} ../../../_static/screenshots/data_loading/h5ad-current-visualization.png
:alt: Cellucid visualizing a current-schema H5AD dataset in two dimensions.
:width: 1280px

The same H5AD dataset rendered as a 2D interactive view with Planar navigation.
```

## Direct Zarr ZIP load

```{figure} ../../../_static/screenshots/data_loading/zarr-zip-loaded.png
:alt: A current-schema AnnData Zarr ZIP archive loaded in Cellucid through the browser picker.
:width: 1280px

The current-schema Zarr ZIP representation loaded through the browser picker.
The Session panel identifies the source as a Zarr ZIP archive, and the
visualization uses the archive's coordinates and categorical metadata.
```

## Browser-engine checks

```{figure} ../../../_static/screenshots/data_loading/h5ad-firefox.png
:alt: Firefox showing Cellucid with the 120-cell synthetic H5AD fixture loaded through the browser picker, its Session panel reporting the dataset summary, and the coloured points rendered in the canvas.
:width: 1280px

The same direct H5AD replacement path exercised in Firefox on macOS.
```

```{figure} ../../../_static/screenshots/data_loading/h5ad-webkit.png
:alt: WebKit showing Cellucid with the same 120-cell synthetic H5AD fixture loaded through the browser picker, its Session panel reporting the dataset summary, and the coloured points rendered in the canvas.
:width: 1280px

The same direct H5AD replacement path exercised in Playwright WebKit on macOS.
WebKit is Safari's browser engine; the repository's platform matrix separately
runs the supported browser projects on their available CI operating systems.
```

## What a refusal looks like

Each capture below is the notification Cellucid actually raised for that exact
input, taken from the running application.

```{figure} ../../../_static/screenshots/data_loading/fail-h5ad-over-size-limit.png
:alt: A notification stating that H5AD direct browser files must have a positive safe size no larger than 512 MiB and suggesting the Cellucid server or the prepared format.
:width: 776px

An `.h5ad` above the 512 MiB browser ceiling, refused before any byte is read.
```

```{figure} ../../../_static/screenshots/data_loading/fail-h5ad-wrong-file.png
:alt: A notification stating that the selected file is not a valid HDF5 or H5AD file and asking for an AnnData .h5ad file.
:width: 760px

A CSV renamed to `.h5ad`. The leading bytes are checked, not the extension.
```

```{figure} ../../../_static/screenshots/data_loading/fail-missing-umap-embedding.png
:alt: A notification listing the expected X_umap_1d, X_umap_2d and X_umap_3d keys and then naming X_pca as the only obsm key present in the selected file.
:width: 760px

A valid `.h5ad` with no supported embedding. The message ends by listing the
`obsm` keys the file does contain.
```

```{figure} ../../../_static/screenshots/data_loading/fail-remote-server-unreachable.png
:alt: A notification stating that the remote server could not be reached because nothing readable came back, naming both a failed connection and a wrong address as possibilities.
:width: 776px

A remote-server address with nothing listening behind it.
```

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: Two stacked notifications, the upper naming the raw.githubusercontent.com URL that returned not found and the lower stating that the GitHub repository could not be reached.
:width: 776px

A GitHub path pointing at a repository rather than at its exports root. The
first notification names the exact URL that 404'd.
```

```{figure} ../../../_static/screenshots/data_loading/keep-dataset-after-failed-connect.png
:alt: The whole Cellucid window after a failed connection: the previously loaded dataset is still rendered and still summarised in the Session panel, with two error notifications in the lower right corner.
:width: 1440px

None of these discard the dataset you already had. After a refusal the previous
dataset is still loaded and still interactive.
```

## Real server evidence

The picker captures above isolate browser-file behavior. The server workflow is
captured against real data in {doc}`04_server_tutorial`, which includes the
prepared-export startup, the direct-H5AD startup and its `?anndata=true` Viewer
URL, the complete `--help` output, and both common startup failures.

Jupyter embedding is captured in
{doc}`../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery`,
which includes the running notebook embed and its exact connection report.
