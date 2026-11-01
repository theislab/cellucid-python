# What is a Cellucid-ready AnnData?

**Audience:** everyone (layered: quick intuition first, details later)  
**Time:** 10 minutes  
**Goal:** know what `cellucid_prepare()` writes to disk and why it is shaped that way.

## Fast path (wet lab / non-technical)

`cellucid_prepare()` writes **one dataset** as **one store**:

- a **`.zarr` folder** (`format = "zarr"`), or
- a **`.zarr.zip` file** (`format = "zarr.zip"`), the same files packed into one ZIP.

`format` has no default: the folder suits GitHub, the single file everything
else ({doc}`../c_data_preparation_api/index`).

You open that folder or file in the Cellucid web app. You can also hand it to a
Python user: it is a standard **AnnData** store, the format Scanpy and most of
the Python single-cell ecosystem use, and `anndata.read_zarr()` opens it. An
object with several assays is written as a **Cellucid-ready MuData** instead,
each assay a Cellucid-ready AnnData under `mod/<name>`, which
`mudata.read_zarr()` opens
({ref}`its layout <r-cellucid-ready-mudata>`).

## Practical path (computational users)

### What is inside

```text
my_dataset.zarr/
  .zgroup .zattrs .zmetadata     root group, AnnData encoding, consolidated metadata
  X/                             expression: gene-major float32, float64 for counts above 2^24 (csc_matrix or dense)
  obs/                           cell metadata data frame (_index = cell names)
  var/                           gene metadata data frame (_index = gene names)
  obsm/                          every embedding, every vector field (float32)
  obsp/connectivities/           neighbour graph (float32 csr_matrix), when supplied
  varm/ varp/ layers/            empty, as anndata writes them
  uns/title                      the dataset name
  uns/cellucid/                  schema, id, name, description, citation, source, defaults
```

`X` and `obsp/connectivities` are present only when you pass
`gene_expression` and `connectivities`.

### Why AnnData Zarr

- **One format everywhere.** The web app, Python's `anndata`, and the Cellucid
  Python package all read it; there is no Cellucid-only conversion step.
- **Gene-major expression.** A sparse matrix is stored as CSC (genes are
  columns), so one gene's values are one contiguous slice: the viewer reads a
  gene with one small range read instead of scanning the matrix.
- **Chunked for range reads.** Every array is split into compressed chunks; the
  viewer fetches only the chunks a view needs, from a local folder or from a web host.
- **Byte shuffle and deflate.** Arrays of values wider than one byte are
  byte-shuffled (except the expression values), then deflate-compressed, which
  browsers decompress natively.
- **Consolidated metadata.** `.zmetadata` holds every `.zgroup`, `.zattrs`, and
  `.zarray` document, so a store on a web host opens with one request.

### What the viewer reads from `uns['cellucid']`

`uns['cellucid']` is what the store says about itself: its `id` (sessions and
community annotations are keyed on it), `name`, `description`, `citation`,
`source`, the **launch defaults** you set (which embedding each view shows
first, the latent space for outlier quantiles, the gene-name column, …), and,
with `latent_key`, the **outlier quantiles** of every categorical `obs` column
in that latent space (`outlierQuantiles`). A key you did not set is absent, and
the viewer chooses for itself. See
{doc}`03_dataset_identity_and_reproducibility`.

## Deep path (maintainers)

- Every element carries the `encoding-type` / `encoding-version` attributes of
  anndata 0.8 and later (`anndata` 0.1.0 root, `dataframe` 0.2.0,
  `categorical` 0.2.0, `csc_matrix` / `csr_matrix` 0.1.0, `string-array` 0.2.0,
  `dict` 0.1.0, …), so anndata reads the store without a warning.
- In Zarr format 2, chunks are compressed with the compressor
  `{"id": "zlib", "level": 6}`, after numcodecs' `shuffle` filter on every
  array of values wider than one byte except the expression values; format 3
  uses the specification's `blosc` codec (shuffle, zlib level 6) for those
  arrays and `gzip` for the others. Arrays with one row per cell hold up to
  65,536 values per chunk, the other arrays read whole 262,144, the data and
  indices of the expression matrix 16,384.
- A `.zarr.zip` stores every entry uncompressed (the chunks are compressed
  already), with the store at the top level of the archive.
- The same inputs give a byte-identical store, folder or ZIP.

The full layout, dtype by dtype: {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

## Common gotchas

- **The suffix matters.** `out_path` must end in `.zarr` (folder) or, with
  `format = "zarr.zip"`, `.zarr.zip`.
- **A store is replaced only on purpose.** An existing store at `out_path` is
  an error unless `force = TRUE`, and a folder that is not a Zarr store is never
  replaced.
- **Interrupted calls leave nothing behind.** The store is written beside
  `out_path` under a temporary name and moved into place only once complete.
- **Do not edit files inside the store by hand.** Write it again from R.

## Next steps

- How the store reaches the browser: {doc}`02_data_flow_r_to_browser`
- Identity and reproducibility: {doc}`03_dataset_identity_and_reproducibility`
