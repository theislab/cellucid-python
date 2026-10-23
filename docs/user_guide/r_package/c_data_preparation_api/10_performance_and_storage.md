# Performance and storage

**Audience:** anyone preparing large datasets
**Time:** 10 minutes
**Goal:** know how big a store gets, how long writing takes, and what makes it smaller.

## A measured example

Measured on 2026-10-23 on an Apple M1 Pro (macOS 26.7.1, R 4.5.3, `Matrix`
from the conda-forge build), with simulated data:

- 100,000 cells x 2,000 genes, sparse expression at 5% density
  (10,000,000 stored values, exponential values rounded to two decimals);
- `obs` with 4 columns (a 20-level factor, a 3-value character column, an
  integer, a double);
- `X_umap` (2 columns) and `X_pca` (50 columns), standard normal values;
- a symmetric 15-neighbour graph (2,999,642 stored values).

| Output | Size | Write time |
|---|---|---|
| `.zarr` folder (426 files) | 62.3 MB | 3.5 s |
| `.zarr.zip` | 62.3 MB | 2.8 s |
| `.zarr` folder, Zarr format 3 (`zarr_format = 3`) | 62.3 MB | 4.0 s |
| `.zarr` with expression and `X_umap` only | 28.5 MB | 2.1 s |

Where the bytes of the folder went:

| Part | Size | Uncompressed |
|---|---|---|
| `X` (data + indices + indptr) | 27.6 MB | 80 MB |
| `obsm` (`X_umap`, `X_pca`) | 17.6 MB | 20.8 MB |
| `obsp/connectivities` | 16.0 MB | 24 MB |
| `obs` | 1.1 MB | |

The ZIP adds only its entry headers to the same files. Real data compresses
differently from simulated data -- random normal embeddings barely compress
at all -- so treat these numbers as a scale, not a promise.

## What decides the size

- **Stored expression values.** Each stored value costs 4 bytes of float32
  plus 4 bytes of int32 index before compression. The cell count is fixed;
  the gene count is yours to choose (`gene_identifiers`, for example highly
  variable genes plus marker genes). Store the values you colour by --
  normalized, log-transformed -- rather than several layers.
- **Wide embeddings.** Every column is stored. A 50-column PCA is
  50 x 4 bytes per cell; store the columns you will plot or use as a latent
  space.
- **The graph.** Each edge is stored twice (both directions) as 4-byte weight
  plus 4-byte index. A 15-neighbour graph of 100,000 cells is about 19 MB.
  Leave out `connectivities` if you will not draw the graph.
- **obs columns** are small: categorical columns store one small integer per
  cell plus their labels once.

## Folder or ZIP

The two hold the same files; the choice is about where the store goes
({doc}`index` has the measurements):

- **`.zarr` folder:** for GitHub, which refuses files over 100 MB and answers
  many small files at once over HTTP/2. On a host that answers over HTTP/1.1
  (six requests at a time), a folder's one request per chunk makes reading
  many genes or a whole matrix much slower than one file.
- **`.zarr.zip`:** one file, for everything else: single-file hosts (Zenodo,
  Figshare, Hugging Face, a bucket, a lab server), local files,
  `cellucid serve` and Jupyter. The entries are stored uncompressed, so the
  web app still reads only the chunks a view needs, with range requests, and
  the metadata comes last, so the archive opens with one read of its end.

See {doc}`../d_viewing_loading/02_host_stores_for_sharing`.

## How reading is fast

The store is laid out for the viewer's reads:

- **gene-major expression:** one gene's values are one contiguous range, in
  chunks of 65,536 values, so showing a gene reads that gene and not the
  matrix;
- **whole-row embedding chunks:** an embedding is read in a few large chunks;
- **byte shuffle and deflate,** which browsers decompress natively (the
  shuffle puts the similar high bytes of neighbouring values together; the
  expression values themselves compress best unshuffled);
- **consolidated metadata,** so the structure of the store is one request.

Chunk sizes are listed in {doc}`09_cellucid_ready_anndata_specification`.

## Memory while writing

`cellucid_prepare()` encodes and compresses one chunk at a time, but it checks
and converts each input whole (finite values, the float32 range, a sparse
matrix of another class converted to `dgCMatrix`, a `gene_identifiers`
subset), so expect a few temporary vectors the size of the largest input on
top of the inputs themselves. What costs far more is preparing the inputs
badly: converting a sparse matrix to dense, or `t()` of a large dense matrix. Transpose a sparse matrix with `Matrix::t()`, which stays sparse, and
never pass a dense copy of a large sparse matrix.

## Next

- {doc}`11_troubleshooting_prepare`
- {doc}`../i_troubleshooting_index/05_performance_and_disk_usage_issues`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/10_performance_and_storage`
