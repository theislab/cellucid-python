# Cellucid-ready AnnData specification (R writer)

**Audience:** developers, people writing their own readers, and anyone checking a store by hand
**Time:** 20 minutes
**Goal:** the exact layout of the store `cellucid_prepare()` writes.

A store written by `cellucid_prepare()` is a **Cellucid-ready AnnData**: an
ordinary AnnData store in Zarr format 2 (the default) or Zarr format 3
(`zarr_format = 3`), following the AnnData on-disk specification that
`anndata` 0.8 and later write and read (format 3: `anndata` 0.11 and later).
Nothing in it is specific to Cellucid except the contents of
`uns['cellucid']`. Python reads it with `anndata.read_zarr("name.zarr")`
without a warning.

Python's `cellucid.prepare(format = "zarr")` writes the same store, element
for element and value for value, from the same input; the shared contract is
{doc}`../../python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification`.
What only Python writes is what R has no type or library for: `obs` and
`var` columns of NumPy types R lacks (`int64`, `int8`, unsigned, `float32`)
keep those types there, and the `.h5ad` output needs Python (this package
writes no HDF5, so it depends on no HDF5 library).

## Layout

```text
name.zarr/
  .zgroup  .zattrs  .zmetadata          root group, encoding "anndata" 0.1.0
  X/                                    gene expression (only with gene_expression)
    data/  indices/  indptr/            sparse: csc_matrix
  obs/                                  dataframe: _index + one element per column
  var/                                  dataframe: _index + one element per column
  obsm/                                 dict: embeddings and vector fields
  obsp/                                 dict
    connectivities/                     csr_matrix (only with connectivities)
  varm/  varp/  layers/                 empty dicts
  uns/                                  dict
    title                               string
    cellucid/                           dict: the dataset's description
```

Every group has a `.zgroup` (`{"zarr_format": 2}`); every AnnData element
carries `encoding-type` and `encoding-version` in its `.zattrs`. In Zarr
format 3 every group and array has one `zarr.json` instead, its attributes in
the `attributes` member, and the root `zarr.json` holds the consolidated
metadata.

## Elements

| Element | `encoding-type` / version | Contents |
|---|---|---|
| root | `anndata` / `0.1.0` | |
| `X` (sparse input) | `csc_matrix` / `0.1.0`, `shape: [n_obs, n_vars]` | `data` `<f4` (`<f8` for integer counts above 2^24), `indices` `<i4`, `indptr` `<i4` (length `n_vars + 1`); no stored zeros |
| `X` (dense input) | `array` / `0.2.0` | `<f4` (`<f8` for integer counts above 2^24), shape `[n_obs, n_vars]` |
| `obs`, `var` | `dataframe` / `0.2.0`, `_index: "_index"`, `column-order: [...]` | `_index` (`string-array`) and one element per column |
| factor or repeating/missing character column | `categorical` / `0.2.0`, `ordered` | `codes` (`array`, the type pandas gives them: `\|i1` for fewer than 127 categories, `<i2` for fewer than 32,767, else `<i4`; -1 is missing), `categories` (`string-array`) |
| unique character column | `string-array` / `0.2.0` | `\|O` with the `vlen-utf8` filter |
| logical column | `array` / `0.2.0` | `\|b1` |
| logical column with `NA` | `nullable-boolean` / `0.1.0` | `values` and `mask` (`array`, `\|b1`) |
| integer column | `array` / `0.2.0` | `<i4` |
| integer column with `NA` | `nullable-integer` / `0.1.0` | `values` (`<i4`) and `mask` (`\|b1`) |
| double column | `array` / `0.2.0` | `<f8`, missing values as NaN |
| `obsm/<name>` | `array` / `0.2.0` | `<f4`, shape `[n_obs, n_columns]`, every column of the input |
| `obsm/<field>_<basis>` | `array` / `0.2.0` | `<f4`, shape `[n_obs, n_columns]`, as many columns as `obsm/X_<basis>` |
| `obsp/connectivities` | `csr_matrix` / `0.1.0`, `shape: [n_obs, n_obs]` | `data` `<f4`, `indices` `<i4`, `indptr` `<i4` |
| `obsm`, `obsp`, `varm`, `varp`, `layers`, `uns`, `uns/cellucid` and its sub-dicts | `dict` / `0.1.0` | |
| a text value in `uns` | `string` / `0.2.0` | 0-dimensional `\|O`, `vlen-utf8` |
| a number or flag in `uns` | `numeric-scalar` / `0.2.0` | 0-dimensional `<i8` or `\|b1` |
| a list in `uns` | `string-array` or `array` / `0.2.0` | `\|O` or `<i8` |
| `uns/cellucid/outlierQuantiles/fields/<column>` | `array` / `0.2.0` | `<f4`, shape `[n_obs]` |

The arrays inside a sparse matrix (`data`, `indices`, `indptr`) carry no
attributes (an empty `.zattrs`), as anndata writes them. `indices` and `indptr` would be `<i8` for a
matrix with 2^31 or more stored values.

## `uns['cellucid']`

| Key | Type | Present |
|---|---|---|
| `schema` | integer, `1` | always |
| `id` | text (`dataset_id`) | always |
| `name` | text (`dataset_name`) | always |
| `description` | text (`dataset_description`) | when supplied |
| `citation` | dict: `text`, and `url` when supplied | when `citation_text` is supplied |
| `source` | dict: `name`, and `url` when supplied | when `source_name` is supplied |
| `defaults` | dict, below | when any default is set |
| `outlierQuantiles` | dict, below | when `latent_key` is supplied |

`uns['title']` repeats the name, so tools that read a file's title (CELLxGENE
files carry one) show it too.

`defaults` holds only the keys whose argument was set, in this order:

| Key | Type | Argument |
|---|---|---|
| `embedding1d`, `embedding2d`, `embedding3d` | text (an obsm key), or dict `{key: text, columns: int64 array, 0-based}` | `embedding_defaults` |
| `geneIdColumn` | text | `var_gene_id_column` |
| `latentKey` | text | `latent_key` |
| `vectorFieldDefault` | text (the field name, `"velocity"`) | `vector_field_default` |
| `obsKeys` | string array | `obs_keys` |
| `normalizeEmbeddings` | boolean | `normalize_embeddings` |
| `centroidMinPoints` | int64 | `centroid_min_points` |

In Python:

```python
>>> adata.uns["cellucid"]
{'schema': 1, 'id': 'pbmc-3k', 'name': 'PBMC 3k',
 'citation': {'text': 'Doe et al. 2026', 'url': 'https://example.org/paper'},
 'defaults': {'embedding2d': 'X_umap',
              'embedding3d': {'key': 'X_pca', 'columns': array([0, 1, 2])},
              'latentKey': 'X_pca'}}
```

`outlierQuantiles` holds the outlier quantiles of every categorical `obs`
column in the latent space `latent_key` names, which the viewer's outlier
filter reads instead of computing them from the whole latent space:

| Key | Type | Value |
|---|---|---|
| `latentKey` | text | `latent_key` |
| `minPoints` | int64 | `centroid_min_points`, or 10 |
| `fields` | dict of float32 arrays, one value per cell | one entry per `obs` column the viewer shows as categorical (a factor, a character column, a logical column with or without `NA`), in `obs` column order, for a column with a category of at least `minPoints` cells |

A cell's value is the share of its category at its distance from the
category's centre or closer, as float32; `NaN` for `NA` and for a category of
fewer than `minPoints` cells. `?cellucid_prepare` (section *Outlier
quantiles*) states the rule, which Python's `cellucid.prepare()` and the web
app follow bit for bit.

## Arrays and chunks

Every array's `.zarray` (Zarr format 2) declares, for `X/indices` for example:

```json
{"chunks": [65536], "compressor": {"id": "zlib", "level": 6},
 "dimension_separator": ".", "dtype": "<i4", "fill_value": 0,
 "filters": [{"elementsize": 4, "id": "shuffle"}], "order": "C",
 "shape": [450000], "zarr_format": 2}
```

and its `zarr.json` (Zarr format 3):

```json
{"shape": [450000], "data_type": "int32",
 "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [65536]}},
 "chunk_key_encoding": {"name": "default", "configuration": {"separator": "/"}},
 "fill_value": 0,
 "codecs": [{"name": "bytes", "configuration": {"endian": "little"}},
            {"name": "numcodecs.shuffle", "configuration": {"elementsize": 4}},
            {"name": "gzip", "configuration": {"level": 6}}],
 "attributes": {}, "zarr_format": 3, "node_type": "array", "storage_transformers": []}
```

- **Shuffle:** every chunk of an array of values wider than one byte is
  byte-shuffled first -- the first byte of every value, then every second
  byte, and so on -- numcodecs' `Shuffle` filter (format 3: the
  `numcodecs.shuffle` codec), as HDF5's shuffle filter does in an `.h5ad`.
  The similar high bytes of neighbouring values then compress together: the
  gene indices of a sparse `X` shrink by about a third, embeddings and graph
  weights by 7-9%. The expression values themselves -- `X/data`, or a dense
  `X` -- are not shuffled: real expression values repeat whole, which deflate
  finds only unshuffled (shuffled, they grew 14-45% on three atlases, for raw
  counts, `normalize_total` + `log1p`, and proportions alike). Strings and
  one-byte values (`|b1`, `|i1`) are not shuffled either; the `.zarray` of
  such an array declares `"filters": null`.
- **Compressor:** every chunk is then DEFLATE level 6, made by base R's
  `memCompress(type = "gzip")`: a zlib stream (RFC 1950) in format 2, exactly
  what numcodecs' `Zlib` codec reads and writes; the same deflate data in gzip
  framing (RFC 1952) in format 3. Browsers inflate both natively.
- **Strings:** dtype `|O`, `filters: [{"id": "vlen-utf8"}]`, `fill_value: ""`
  (format 3: `data_type` `"string"` and the `vlen-utf8` codec).
- **Fill values:** `0` for numbers, `false` for `|b1`, `""` for strings. Every
  chunk is stored at the full chunk shape; the last chunk along an axis is
  padded with the fill value, as both Zarr formats require.
- **Chunk keys:** format 2: `0`, `1`, ... for one dimension, `<row>.<column>`
  for two (`0.0`, `1.0`, ...), and `0` for a 0-dimensional value; format 3:
  `c/0`, `c/1`, ..., `c/<row>/<column>`, and `c` for a 0-dimensional value.

A chunk is the unit of every read, so its size is chosen by how the viewer
reads each array:

| Array | Chunk |
|---|---|
| `X/data` and `X/indices` (read one gene's range at a time) | 65,536 values (2^16) |
| `X/indptr`, every array of `obsp/connectivities`, `obs`/`var` columns and index, `uns` lists (read whole) | 262,144 values (2^18) |
| `obsm` arrays (read a block of rows, every column) | `[rows, n_columns]`, rows = 262,144 / n_columns (rounded down) |
| dense `X` (read one gene at a time) | `[rows, genes]`, rows = up to 262,144, genes = as many as keep the chunk at 262,144 values (at least 1) |

An array shorter than its chunk size is one chunk of exactly its length; an
empty array has chunks of 1. For
example, a 300,000-cell store holds `obsm/X_umap` in chunks of
`[131072, 2]`, `obsm/X_pca` (50 columns) in chunks of `[5242, 50]`, and a
`string-array` index of 300,000 names in two chunks of 262,144.

## Consolidated metadata

`.zmetadata` holds every `.zgroup`, `.zattrs`, and `.zarray` of the store in
one document, as `zarr.consolidate_metadata()` writes it:

```text
{"metadata": {".zgroup": {"zarr_format": 2}, ".zattrs": {...}, "X/.zattrs": {...}, ...},
 "zarr_consolidated_format": 1}
```

In format 3 the root `zarr.json` carries every other node's document in its
`consolidated_metadata` member (`{"kind": "inline", "must_understand": false,
"metadata": {"X": {...}, "X/data": {...}, ...}}`), as
`zarr.consolidate_metadata()` writes it. A reader on a web host fetches this
one file instead of one request per element.

## `.zarr.zip`

With `format = "zarr.zip"` the same files, byte for byte, are entries of one ZIP archive,
with the store at the top level (`X/data/0`, ..., `.zgroup`, `.zmetadata`; no
enclosing folder):

- every chunk comes first, then the metadata documents, the consolidated
  metadata last: the web app reads the end of an archive when it opens it (the
  central directory), and with it the documents that describe the store, so
  the store opens with that one read;

- every entry is **stored** (method 0): the chunks are compressed already, and
  a stored entry can be read with one HTTP range request;
- entry names are UTF-8 (general-purpose flag bit 11);
- every timestamp is 1980-01-01 00:00 and every entry a 0644 regular file, so
  the same data gives the same bytes;
- ZIP64 records are written only where a size, an offset, or the entry count
  needs them.

Python opens it with `zarr.storage.ZipStore(path, mode="r")` and
`anndata.read_zarr()`; the web app opens it as a file or a link
({doc}`../d_viewing_loading/01_open_stores_in_cellucid_web_app`).

## Writing

The store is written under a temporary hidden name beside `out_path`
(`.<name>.<process id>.partial`) and renamed to `out_path` only when complete;
see {doc}`01_cellucid_prepare_overview`.

## Next

- {doc}`10_performance_and_storage`
- {doc}`../b_concepts_mental_models/01_what_is_a_cellucid_ready_anndata`
