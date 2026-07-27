# Verified data-loading captures

These captures were produced through the running web application with a
deterministic synthetic AnnData fixture. The acceptance path replaced the
prepared dataset with its Zarr ZIP representation, replaced that with the H5AD
representation, and then restored the prepared dataset through the visible
controls. Each replacement asserted the 120-cell, 6-gene dataset shape, 2D
Planar navigation, hidden camera transport, settled notifications, and the
absence of browser-console, page, and HTTP errors.

## Loading and session controls

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 246px

The Session panel presents each loading path separately and keeps Save State and
Load State beside the dataset controls.
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
:alt: A current-schema H5AD dataset loaded in Firefox.
:width: 1280px

The same direct H5AD replacement path exercised in Firefox on macOS.
```

```{figure} ../../../_static/screenshots/data_loading/h5ad-webkit.png
:alt: A current-schema H5AD dataset loaded in WebKit.
:width: 1280px

The same direct H5AD replacement path exercised in Playwright WebKit on macOS.
WebKit is Safari's browser engine; the repository's platform matrix separately
runs the supported browser projects on their available CI operating systems.
```
