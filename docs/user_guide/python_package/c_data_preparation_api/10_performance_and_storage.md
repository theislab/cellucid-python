# Performance and storage

**Audience:** anyone hosting or sharing large datasets
**Time:** 10 minutes
**Goal:** choose a container and know what a prepared file costs to write and to read

## Which container

`prepare()` asks for it: `format` is required and has no default.

| `format` | Writes | Best for |
|---|---|---|
| `"zarr"` | a `.zarr` folder of small files | GitHub: no file over 100 MB, no LFS, and GitHub answers many small files at once over HTTP/2. |
| `"h5ad"` | one `.h5ad` file | Everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve` and Jupyter; read by every AnnData tool. |
| `"zarr.zip"` | one `.zarr.zip` file | Wherever an `.h5ad` goes, as Zarr. |
| `"h5mu"` | one `.h5mu` file, for a MuData (`data=mdata`) | A MuData wherever an `.h5ad` goes; `"zarr"` and `"zarr.zip"` write a MuData too. |

A single file is read in ranges that the browser merges into few requests. A
folder is one file per chunk: on a host that answers over HTTP/1.1 (six
requests at a time), reading many genes or a whole matrix from it costs many
round trips; over HTTP/2 (GitHub) or from disk it costs little.
{ref}`which-format` has the measurements.

## Measured

The Suo atlas (561,947 cells × 8,192 genes of raw counts), written by
`prepare()`, read by the web app's own reader through every step a session
takes (open, every `obs` field, the embeddings, genes, 1,000 cells, the graph,
every gene in two cell types). Full method and the per-step numbers:
{ref}`which-format` and {doc}`09_cellucid_ready_anndata_specification`.

| | input `.h5ad` (CSR float64) | `.h5ad` | `.zarr` (v2) | `.zarr.zip` |
|---|---|---|---|---|
| Size | 337 MB | 227.7 MB | 198.1 MB | 198.3 MB |
| First gene, HTTP/1.1, 80 ms, 50 MB/s | 4,110 ms (117.5 MB) | 86 ms | 87 ms | 85 ms |
| Every step, HTTP/2, 80 ms, 50 MB/s | — | 8.3 s | 8.0 s | 8.2 s |
| Every step, HTTP/1.1, 80 ms, 50 MB/s | — | 9.0 s | 46.8 s | 9.0 s |
| Every step, local file in Chromium | — | 2.5 s | 2.6 s | 2.7 s |

## Writing large datasets

- `prepare()` holds the converted arrays in memory next to its inputs: the
  expression matrix as float32 CSC, the embeddings as float32. Plan for a few
  times the size of the expression matrix's stored values.
- A sparse `X` that is not in memory is never read whole: a backed AnnData's
  (`read_h5ad(path, backed="r")`) is read from its file in blocks, and a
  lazily read one (a dask array of sparse blocks, `anndata.experimental.read_lazy`)
  is decoded block of rows by block of rows into a scratch folder in the
  temporary folder (removed afterwards), then laid out from there; the same
  holds for an AnnData `show_anndata()` serves (its X, layers and raw). A
  dense dask array is computed.
- Compressing the chunks is most of the time `prepare()` takes. An `.h5ad`'s
  chunks are compressed on every core and stored with HDF5's direct chunk
  write, in the layout HDF5 itself writes; a Zarr store's are compressed
  concurrently by zarr-python. A `.zarr.zip` adds the time to pack the folder.
- `gene_identifiers=[...]` keeps chosen genes; subsetting cells before the call
  (`adata[mask]`) shrinks everything else.

## Reading in the browser

- One gene is one range of `X/data` and `X/indices` (chunks of 16,384 values),
  so colouring by a gene costs a few requests whatever the dataset size; in
  priority order a subsample reads only the first values of each gene.
- Embeddings and `obs` columns are read by rows, in chunks of 65,536 values
  (the first rows only, for a subsample); the graph is read whole, in chunks
  of 262,144 values.
- Centroids are computed by the browser. Outlier quantiles are stored with
  `latent_key` (4 bytes per cell per categorical column, before compression), and the browser
  reads one column's on the outlier slider's first move; without them, the
  browser computes them on that first move by reading the whole latent space
  once.

## Next

- {doc}`11_troubleshooting_prepare`
- {doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/10_performance_and_storage`
