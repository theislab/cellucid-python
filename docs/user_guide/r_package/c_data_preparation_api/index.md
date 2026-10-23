# Data Preparation API

This section documents the R package's writer in full:

- `cellucid::cellucid_prepare(...)`

and, with the vector fields it stores, `cellucid::cellucid_transition_drift(...)`,
which computes one from a cell-cell transition matrix
({doc}`08_vector_fields_velocity_displacement`).

`cellucid_prepare()` writes a **Cellucid-ready AnnData**: a standard AnnData
store in the Zarr format (a `.zarr` folder, or one `.zarr.zip` file) that the
Cellucid web app opens directly and Python's `anndata` reads like any other
AnnData file. These pages cover:

- what each input must be (shapes, types, ordering, names),
- where each input ends up in the store,
- the exact store layout,
- size and speed, and
- what each error message means.

## Which format?

`format` has no default, because where the store goes decides it:

| The store goes to | `format` | Why |
|---|---|---|
| a GitHub repository, or its GitHub Pages site | `"zarr"`: a `.zarr` folder | GitHub refuses files over 100 MB; a folder of small files publishes any size without Git LFS. GitHub answers over HTTP/2, 100 requests at a time, so one request per chunk costs a folder nothing. |
| everything else: a single-file host (Zenodo, Figshare, Hugging Face, S3, a lab web server), your own disk (the web app's file picker), `cellucid serve`, Jupyter | `"zarr.zip"`: one `.zarr.zip` file | One file, read in ranges that the reader merges into few requests. Most single-file hosts answer over HTTP/1.1, six requests at a time, where a folder's one request per chunk makes many-gene and whole-matrix reads several times slower (below). |

This package writes Zarr only, so it needs no HDF5 library; Python's
`cellucid.prepare(format = "h5ad")` writes the `.h5ad` that single-file hosts
and AnnData tools also take, and reads as fast as a `.zarr.zip`. In a browser,
a local Zarr folder reads as fast as a single file.

Measured on the Suo atlas (561,947 cells × 8,192 genes; the store Python's
`prepare()` writes from the same data, element for element the store this
function writes) with the web app's own reader, every step of opening it and
reading `obs`, the embeddings, genes, 1,000 cells, the graph and every gene in
two cell types:

| Every step, in total | `.zarr` folder | `.zarr.zip` |
|---|---|---|
| Size | 198.1 MB, 1,961 files | 198.3 MB |
| Local file, in Chromium | 2.6 s | 2.7 s |
| HTTP/2, 80 ms round trip, 50 MB/s (GitHub) | 8.0 s | 8.2 s |
| HTTP/1.1, 80 ms, 50 MB/s (a typical single-file host) | **46.8 s** | 9.0 s |
| HTTP/1.1, 150 ms, 5 MB/s | **103.5 s** | 48.7 s |

The Python package's guide has the full table, with `.h5ad`:
{ref}`which-format`.

```{note}
For a quick first store, start with {doc}`../a_landing_pages/04_quick_start_3_levels`
and come back here for details.
```

**Recommended reading order**

1) {doc}`01_cellucid_prepare_overview`
2) {doc}`02_input_requirements_global`
3) {doc}`03_embeddings_and_coordinates`
4) {doc}`04_obs_cell_metadata`
5) {doc}`06_gene_expression_matrix` (+ {doc}`05_var_gene_metadata`)
6) {doc}`07_connectivities_knn_graph` (optional)
7) {doc}`08_vector_fields_velocity_displacement` (optional)
8) {doc}`09_cellucid_ready_anndata_specification`
9) {doc}`10_performance_and_storage`
10) {doc}`11_troubleshooting_prepare`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
