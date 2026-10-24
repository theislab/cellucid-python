# Data Flow: R → Store → Browser

**Audience:** everyone  
**Time:** 5 minutes  
**Goal:** choose the right way to get your store into Cellucid.

R writes the store. The **browser** reads it, chunk by chunk, wherever it lives.

## Three ways to open a store

### Option 1 — Open it from your computer (simplest)

1) Write a store with `cellucid_prepare(out_path = "my_dataset.zarr", format = "zarr", ...)`.
2) Open the Cellucid web app.
3) Open the `.zarr` folder (or the `.zarr.zip` file) from the data panel.

Nothing leaves your machine, and no server runs. A browser folder picker may
leave out the hidden metadata files (`.zgroup`, `.zmetadata`, …) on some
systems; a `.zarr.zip` (`format = "zarr.zip"`) avoids that.

Start here: {doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`

### Option 2 — Open it from a link (for sharing)

Put the `.zarr` folder on GitHub, or the `.zarr.zip` file on any host that
serves files over HTTPS (a lab server, a bucket, Zenodo, Figshare, Hugging
Face), and paste its address into the **Link or server** bar. The viewer reads it in place with range
requests, so only what a view needs is downloaded.

Start here: {doc}`../d_viewing_loading/02_host_stores_for_sharing` and
{doc}`../../web_app/b_data_loading/14_open_from_link`

### Option 3 — Continue in Python

A store written by R is a standard AnnData store: `anndata.read_zarr()` opens
it, and the Cellucid Python package's AnnData workflows (a local server,
Jupyter embedding) take the resulting object like any other:
{doc}`../../python_package/index`.

## Where `cellucid-annotation` fits

Community annotation (`cellucid-annotation`) keys a dataset's annotations on
its `dataset_id`, so keep that id stable (see
{doc}`03_dataset_identity_and_reproducibility`).

## Privacy and security (do not skip)

- **Opening from your computer** keeps the data on your machine.
- **Hosting at a public address** is publishing; treat it like making the
  dataset public unless the host enforces access control.

## Next steps

- Write a store: {doc}`../c_data_preparation_api/index`
- Open it: {doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`
