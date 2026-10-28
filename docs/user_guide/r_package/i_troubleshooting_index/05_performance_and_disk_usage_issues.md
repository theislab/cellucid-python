# Performance and Disk Usage Issues

**Audience:** anyone writing medium or large datasets  
**Time:** 5–10 minutes  
**What you'll get:** why writing is slow or a store is large, and the levers that fix it

## How big and how fast to expect

The store holds what you pass: `X` as float32 (the data and int32 indices of a
sparse matrix; float64 for integer counts above 2^24), embeddings and vector fields as float32, `obs` and `var`
columns as written, every chunk zlib-compressed. On a simulated 50,000 cells x
2,000 genes sparse matrix (5 million non-zero values) with a 2-column UMAP, a
50-column PCA and two `obs` columns, `cellucid_prepare()` wrote a 21.7 MB store
in about 1 second, as 225 files or as one `.zarr.zip` of the same size (macOS,
Apple silicon, R 4.5.3). Time and size grow with the number of stored values.

## Symptom: writing takes “forever”

**Likely causes**
- A dense `gene_expression` with many zeros: every value is stored and
  compressed.
- A slow or synced target folder (network drive, cloud-synced folder).

**Fix**
- Pass expression as a sparse `Matrix` (for example
  `Matrix::Matrix(expr, sparse = TRUE)`); it is stored as a gene-major CSC
  matrix with only the non-zero values.
- Store only the genes you need with `gene_identifiers`.
- Write to a local disk and copy the finished store afterwards.

## Symptom: the store is too large

**Fix**
- Store a gene panel (`gene_identifiers`) or a subset of cells.
- Drop `obs` / `var` columns you do not need before calling; every column is
  stored (`obs_keys` only chooses what the viewer shows first).
- Embeddings are stored with every column: pass only the PCA components you
  want to keep.

## Symptom: too many files for my host or filesystem

A `.zarr` folder holds one file per chunk plus the metadata files. Write one
file instead with `format = "zarr.zip"` and an `out_path` ending in `.zarr.zip`; the
entries are stored without further compression, so the web app still reads any
chunk with one range request. See
{doc}`../c_data_preparation_api/10_performance_and_storage`.

## Symptom: the web app is slow to show a gene

Expression is stored gene-major, so one gene is one small range of the store,
read in chunks of 16,384 values. If a gene still loads slowly from a remote
host, the host is usually the bottleneck; see
{doc}`../d_viewing_loading/02_host_stores_for_sharing`.
