# Preparing an AnnData for Cellucid

The Cellucid web app reads AnnData itself: an `.h5ad` file, a Zarr store
folder (`.zarr`), or a zipped Zarr store (`.zarr.zip`). Any AnnData opens.
{func}`~cellucid.prepare` writes a **Cellucid-ready AnnData**: the same
standard format, laid out so the browser reads it fast (the expression matrix
gene-major, chunks sized for HTTP range reads, codecs browsers decode
natively) and carrying the dataset's name, citation and launch defaults in
`uns['cellucid']`. Given a MuData it writes a Cellucid-ready MuData
(`.h5mu`, or MuData Zarr) the same way
({doc}`12_prepare_from_anndata_or_mudata`). The web app opens other formats
too — Loom, 10x Cell Ranger files and `outs/` folders, MatrixMarket folders,
R objects — without preparation ({doc}`../d_viewing_apis/index`).

```python
import cellucid

cellucid.prepare(
    "pbmc.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    latent_key="X_pca",
)
```

The result is an ordinary AnnData (`anndata.read_h5ad("pbmc.h5ad")` reads it)
and opens with `cellucid serve pbmc.h5ad`, {func}`~cellucid.show_anndata`, the
web app's file picker, or a link to wherever you host it. The R package's
`cellucid_prepare()` writes the same format from Seurat or
SingleCellExperiment objects; see {doc}`../../r_package/index`.

(which-format)=
## Which format?

`format` has no default, because where the data goes decides it:

| The data goes to | `format` | Why |
|---|---|---|
| a GitHub repository, or its GitHub Pages site | `"zarr"`: a `.zarr` folder | GitHub refuses files over 100 MB; a folder of small files publishes any size without Git LFS. GitHub answers over HTTP/2, 100 requests at a time, so one request per chunk costs a folder nothing. |
| a single-file host: Zenodo, Figshare, Hugging Face, S3, a lab web server | `"h5ad"` (or `"zarr.zip"`) | One file, read in ranges that the reader merges into few requests. Most such hosts answer over HTTP/1.1, six requests at a time, where a folder's one request per chunk makes many-gene and whole-matrix reads several times slower (below). |
| your own disk (the web app's file picker), `cellucid serve`, Jupyter | `"h5ad"` (or `"zarr.zip"`) | Every AnnData tool opens an `.h5ad`. In a browser, a local Zarr folder reads as fast as a single file. |
| from R | `"zarr"` for GitHub, `"zarr.zip"` for everything else | `cellucid_prepare()` writes Zarr only, so it needs no HDF5 library. |

A MuData (`data=mdata`) follows the same table with `"h5mu"` in place of
`"h5ad"`; `"zarr"` and `"zarr.zip"` write it as MuData Zarr.

Measured on the Suo atlas (561,947 cells × 8,192 genes, written by
`prepare()`) with the web app's own reader (`scripts/benchmark-anndata.mjs` in
the web repository): open the file, read every `obs` field and the embeddings,
a first and a next gene, 99 genes at once, 1,000 cells, the neighbour graph,
and every gene in the cells of the two largest cell types (what a differential
expression test reads). Remote links are simulated: HTTP/2 with 100 requests
at a time as GitHub answers, HTTP/1.1 with six as most web servers do. Local
files were read in Chromium with the reader's decode workers.

| Every step, in total | `.h5ad` | `.zarr` folder | `.zarr.zip` |
|---|---|---|---|
| Size | 227.7 MB | 198.1 MB, 1,961 files | 198.3 MB |
| Local file, in Chromium | 2.5 s | 2.6 s | 2.7 s |
| HTTP/2, 80 ms round trip, 50 MB/s (GitHub) | 8.3 s | 8.0 s | 8.2 s |
| HTTP/2, 150 ms, 5 MB/s | 46.9 s | 46.1 s | 46.7 s |
| HTTP/1.1, 80 ms, 50 MB/s (a typical single-file host) | 9.0 s | **46.8 s** | 9.0 s |
| HTTP/1.1, 150 ms, 5 MB/s | 49.3 s | **103.5 s** | 48.7 s |
| Opening, at 80 ms | 0.25 s | 0.19 s (HTTP/2), 0.27 s (HTTP/1.1) | 0.13 s |

Over HTTP/1.1 a folder is slow wherever a view reads many chunks: the two
cell types' genes took 20.5 s from the folder and 2.2–2.4 s from one file,
1,000 cells 21.3 s and 2.9–3.0 s. A gene, an embedding or the graph cost the
same from every format. `.h5ad` and `.zarr.zip` read equally fast; the
`.zarr.zip` is 13% smaller (HDF5 stores text — cell names, categories —
uncompressed) and opens with the one read of its end, while every AnnData tool
opens an `.h5ad` directly. In the whole app, over the HTTP/1.1 link above, Suo
opened in 1.8 s from the `.h5ad` (2.6 s before this layout), 1.7 s from the
`.zarr.zip` and 1.4 s from the folder, and showed a gene in 0.3 s from each.

```{note}
`prepare()` is optional. To look at an AnnData you already have, open it
directly: {doc}`../d_viewing_apis/index`. Prepare a file when you will share or
host it, when it is large (a gene-major matrix makes every gene one small
read), or when it should open with a name, a citation and chosen defaults.
```

## Reading order

1) {doc}`01_prepare_overview` — what `prepare()` writes, and the shortest path to a file
2) {doc}`02_input_requirements_global` — shapes, cell order, names, validation
3) {doc}`03_embeddings_and_coordinates` — `embeddings=` and `embedding_defaults=`
4) {doc}`04_obs_cell_metadata` — `obs=` and `obs_keys=`
5) {doc}`05_var_gene_metadata` and {doc}`06_gene_expression_matrix` — `var=`, `gene_expression=`, `gene_identifiers=`
6) {doc}`07_connectivities_knn_graph` — `connectivities=` (optional)
7) {doc}`08_vector_fields_velocity_displacement` — `vector_fields=` (optional)
8) {doc}`09_cellucid_ready_anndata_specification` — the written contract every Cellucid repository shares
9) {doc}`10_performance_and_storage` — sizes, writing large datasets, `.h5ad` vs `.zarr` vs `.zarr.zip`
10) {doc}`11_troubleshooting_prepare` — every refusal and what to do about it
11) {doc}`12_prepare_from_anndata_or_mudata` — `data=`: an AnnData or a MuData whole, a Cellucid-ready MuData

## API reference

- {doc}`../g_api_reference_coverage/api/prepare` ({func}`~cellucid.prepare`, with its full docstring)

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
