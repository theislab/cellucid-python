# Cellucid-ready AnnData specification

**Audience:** developers, people writing their own writers or readers, and anyone checking a file by hand
**Time:** 20–30 minutes
**Goal:** the exact contract a Cellucid-ready AnnData satisfies, and the layout the two writers choose

This page is **the written contract** shared by every Cellucid repository:
{func}`~cellucid.prepare` (Python) and `cellucid_prepare()` (R) write it, the
web app reads it, and the dataset repositories ship it. The R package documents
its own writer in {doc}`../../r_package/c_data_preparation_api/09_cellucid_ready_anndata_specification`;
both writers produce what this page specifies.

A **Cellucid-ready AnnData** is an ordinary AnnData: an `.h5ad` file, a Zarr
store folder (`.zarr`, Zarr format 2 or 3), or that folder zipped
(`.zarr.zip`). It follows the AnnData on-disk specification that `anndata` 0.8
and later write and read, with no custom encodings, so

```python
import anndata as ad

adata = ad.read_h5ad("pbmc.h5ad")   # or ad.read_zarr("pbmc.zarr")
```

reads it without a warning. Nothing in it is specific to Cellucid except the
contents of `uns['cellucid']` and the way its arrays are laid out for the
browser. Any AnnData opens in Cellucid; a Cellucid-ready one reads faster (§3) and
names itself.

---

## 1. Elements

| Element | Contract |
|---|---|
| `X` | The gene expression, cells × genes, **gene-major**: a `csc_matrix` when sparse, a dense array otherwise. `float32`; `float64` only for integer counts above 2^24, which `float32` cannot hold exactly (below). Every stored value finite; a sparse `X` has no stored zeros. Absent when the dataset has no expression. |
| `obs` | One row per cell. The index (`_index`) holds unique cell names. Text columns whose values repeat or are missing are `categorical` (categories in code-point order, missing = code -1); other text columns are `string-array`; categoricals stay categoricals; numbers stay numbers in their own dtype (no infinities; NaN is missing); booleans stay booleans; integer and boolean columns with missing values are `nullable-integer` / `nullable-boolean`. |
| `var` | One row per column of `X`. The index holds unique gene names (usually gene ids); a gene-symbol column, when present, is named by `uns['cellucid']['defaults']['geneIdColumn']`. |
| `obsm/<key>` | Embeddings: dense `float32`, two-dimensional, one row per cell, every column of the input at its own coordinates (the viewer selects columns and normalizes). The latent space for outlier measurement is one of them. |
| `obsm/<field>_<basis>` | A vector field, dense `float32`, one row per cell and exactly as many columns as the embedding `obsm['X_<basis>']` it is named for (scVelo's naming: `velocity_umap` over `X_umap`, `velocity_branches` over `X_branches`). The viewer reads it as the field `<field>` drawn over that embedding: the vector-field key rule below. |
| `obsp/connectivities` | The neighbour graph, when present: `csr_matrix`, `float32`, cells × cells, exactly symmetric, finite non-negative weights, zero diagonal, no stored zeros. |
| `varm`, `varp`, `layers` | Present and empty. |
| `uns['title']` | The dataset name, so tools that read a file's title (every CELLxGENE file carries one) show it too. |
| `uns['cellucid']` | The dataset's description and launch defaults (§2). |

Element names that become paths in a store — obs and var column names, obsm
keys — are portable path components on every platform: no `/ \ < > : " | ? *`
or control character, not `.`/`..` or a Zarr metadata name (`.zarray`,
`.zattrs`, `.zgroup`, `.zmetadata`, `zarr.json`), no trailing space or dot, no
Windows device name (`CON`, `PRN`, `AUX`, `NUL`, `COM1`–`9`, `LPT1`–`9`), and
names in one group differ by more than letter case. A column may not be named
`_index`.

### Vector-field keys

Which `obsm` arrays are vector fields, and where they are drawn, is decided by
their keys alone, following scVelo's naming convention. This rule is
normative: the web app (`parseVectorFieldKey`), the Python package
(`cellucid.data_mapping.parse_vector_field_key`) and the R package apply
it identically.

1. **`X_` names an embedding.** An `obsm` key that starts with `X_` is never a
   vector field.
2. **Split.** Any other key is split at one `_` into a non-empty `<field>` and
   a non-empty `<basis>`. A split qualifies when `X_<basis>` is a key of the
   same object's `obsm`. Matching is exact and case-sensitive.
3. **Longest basis.** When several splits qualify, the one with the longest
   basis wins. With both `X_umap_3d` and `X_3d` present, `velocity_umap_3d` is
   the field `velocity` on `X_umap_3d`; with `X_fwd_umap` present, `T_fwd_umap`
   is the field `T` on `X_fwd_umap`.
4. **No split.** A key with no qualifying split is an ordinary `obsm` array,
   which can be plotted as an embedding.

| `obsm` keys | Vector field | Drawn over |
|---|---|---|
| `velocity_umap`, `X_umap` | `velocity` | `X_umap` |
| `velocity_umap_3d`, `X_umap_3d` | `velocity` | `X_umap_3d` |
| `T_fwd_umap`, `X_umap` | `T_fwd` | `X_umap` |
| `velocity_branches`, `X_branches` | `velocity` | `X_branches` |
| `velocity_pca`, `X_pca` | `velocity` | `X_pca` |
| `velocity_umap_2d`, `X_umap` (no `X_umap_2d`) | none: an ordinary array | — |

**Column count.** A vector field is bound when the dataset opens only when it
and its embedding are usable dense arrays with the **same number of columns**,
any number of them. The viewer binds it to every enabled representation of
`X_<basis>` and reads the vector at that representation's columns:
`velocity_pca` (50 columns) over `X_pca` (50 columns) shown at columns
`[0, 1]` reads the vector's columns `[0, 1]`. A key named for an embedding of
another width is still a vector field (never offered as an embedding) but
stays unbound until it is bound by hand under **Session › Data mapping ›
Vector fields**. Several fields of one basis all bind (`velocity_umap` and
`T_fwd_umap` both on `X_umap`).

**Field id.** A key bound to the embedding it is named for is listed under its
`<field>` part, so the keys of one field on every basis are one field:
`velocity_umap`, `velocity_umap_3d` and `velocity_branches` are all the field
`velocity`. A key bound by hand to another embedding is listed under its own
full key. `vectorFieldDefault` (§2) names a field by this id.

`prepare(vector_fields={field: {embedding_key: array}})` writes each array as
`obsm['<field>_<basis>']` for the embedding `X_<basis>`, and refuses an
embedding key that does not start with `X_`, an array whose column count
differs from its embedding's, and any key that would not read back as its own
field and embedding under this rule
({doc}`08_vector_fields_velocity_displacement`).

Nothing else is stored in `obs`, `obsm` or `obsp`: no centroids, no quantized
or pre-normalized copies. The browser computes centroids itself. The one
derived result the writers store is the latent outlier quantiles of the
categorical fields, in `uns` ({ref}`§2 <outlier-quantiles-spec>`), because
computing them needs the whole latent space.

### AnnData encodings

| Element | `encoding-type` / `encoding-version` | Arrays |
|---|---|---|
| root | `anndata` / `0.1.0` | |
| `X` (sparse) | `csc_matrix` / `0.1.0`, `shape: [n_obs, n_vars]` | `data` (`<f4`, or `<f8` per §1), `indices` `<i4`, `indptr` `<i4` (length `n_vars + 1`) |
| `X` (dense) | `array` / `0.2.0` | `<f4` (or `<f8`), shape `[n_obs, n_vars]` |
| `obs`, `var` | `dataframe` / `0.2.0`, `_index: "_index"`, `column-order` | `_index` (`string-array`) and one element per column |
| categorical column | `categorical` / `0.2.0`, `ordered` | `codes` (`array`, the type pandas gives them: `\|i1` for fewer than 127 categories, `<i2` for fewer than 32,767, `<i4` for fewer than 2^31 − 1, else `<i8`), `categories` (`string-array`) |
| text column | `string-array` / `0.2.0` | variable-length UTF-8 |
| number or boolean column | `array` / `0.2.0` | its own dtype (`<f8`, `<i8`, `\|b1`, ...) |
| integer or boolean column with missing values | `nullable-integer` / `nullable-boolean`, `0.1.0` | `values` and `mask` |
| `obsm/<key>` | `array` / `0.2.0` | `<f4`, shape `[n_obs, n_columns]` |
| `obsp/connectivities` | `csr_matrix` / `0.1.0`, `shape: [n_obs, n_obs]` | `data` `<f4`, `indices` `<i4`, `indptr` `<i4` |
| `obsm`, `obsp`, `varm`, `varp`, `layers`, `uns` and every mapping in it | `dict` / `0.1.0` | |
| a text value in `uns` | `string` / `0.2.0` | 0-dimensional |
| a number or flag in `uns` | `numeric-scalar` / `0.2.0` | 0-dimensional `<i8`, `<f8` or `\|b1` |
| a list in `uns` | `string-array` or `array` / `0.2.0` | one-dimensional |

The arrays inside a sparse matrix carry no attributes (an empty `.zattrs`), as
anndata writes them. `indptr` is `<i8` for a matrix with 2^31 or more stored
values; `indices` stays `<i4` unless the axis it indexes (the cells of a CSC
`X`) has 2^31 entries or more, so a matrix of billions of values keeps 4-byte
positions. Readers take each array's own type.

**The type of `X`.** `X` is `float64` when every value it stores is a whole
number and one exceeds 2^24 in magnitude (an integer counts matrix `float32`
would round); otherwise `float32`. The values considered are the stored ones (a
sparse matrix's non-zeros, a dense matrix's every value) of the genes that are
stored, after `gene_identifiers` has chosen them, so a large count in a gene
left out does not widen the genes kept. Integer and boolean inputs are whole
numbers by type; a `float32` input (Python only) stays `float32`.

**The type of an `obs` or `var` column** is the column's own: each writer stores
the dtype it was given. The two writers produce the same element for the same
column wherever the languages share the type — an R `integer` and a NumPy
`int32` column are both `<i4` (`nullable-integer` with `<i4` values when they
hold a missing value), an R `double` and a `float64` column `<f8` whether or
not their values are whole numbers, an R `logical` and a `bool` column `\|b1`,
an R `factor` and a pandas categorical with text categories the same
`categorical`. R has no other number types, so the R writer never writes the
types only Python has: a pandas `int64`, `int8`, unsigned or `float32` column
keeps that dtype in the Python writer (`<i8`, `\|i1`, `<u4`, `<f4`, ...), and a
pandas categorical may have number categories. Under a missing value (`mask`
set), the `values` of a nullable column are unspecified to a reader; both
writers store 0 or `false` there.

---

## 2. `uns['cellucid']`

What the file says about itself. Every key is optional to a reader, and a file
without `uns['cellucid']` opens as any AnnData does; a key that is present must
be exactly right, and a key not listed here makes the web app refuse the file,
so a misspelt default never passes unnoticed. The writers always write
`schema`, `id` and `name`.

| Key | Type | Meaning |
|---|---|---|
| `schema` | integer `1` | The version of this mapping. The web app reads schema 1 and asks for a newer Cellucid for a higher one. |
| `id` | text | The dataset's stable id: 1–180 ASCII letters, digits, `.`, `_` or `-`, starting with a letter or digit and not ending with `.`. Sessions and community annotations are keyed by it. |
| `name` | text | The name shown for the dataset. |
| `description` | text (may be empty) | One paragraph about the dataset. |
| `citation` | mapping `{text, url}` | How to cite the dataset; `url` optional, an absolute `http(s)` address. |
| `source` | mapping `{name, url}` | Where the data comes from; `url` optional, an absolute `http(s)` address. |
| `defaults` | mapping, below | The launch defaults the dataset opens with. |
| `outlierQuantiles` | mapping, below | Every categorical obs field's latent outlier quantiles in one latent space. |
| `session` | mapping | The record {func}`~cellucid.apply_cellucid_session_to_anndata` keeps of a session applied to the file; the viewer never reads it. |

Text is stored as one `string` scalar (the web app also reads a one-element
list of text as that text) and must be drawable exactly as stored: no control or zero-width
character and no leading or trailing white space.

### `defaults`

Only the keys that were set are written, in this order; an absent key keeps
the viewer's own choice. The names are the web app's
(`requireAnnDataLaunchDefaults`), which checks every value again.

| Key | Type | Meaning | `prepare()` argument |
|---|---|---|---|
| `embedding1d`, `embedding2d`, `embedding3d` | text (an `obsm` key, plotted by its first columns), or mapping `{key, columns}` (columns a 0-based integer list, one per axis) | The embedding the dimension opens on. | `embedding_defaults` |
| `knn` | text (an `obsp` key) | The neighbour graph drawn. (Not written by `prepare()`; the viewer's own rule finds `obsp['connectivities']`.) | — |
| `modality` | text (a modality: a `mod/<name>` of a Cellucid-ready MuData, §5) | The modality the viewer opens on; `geneIdColumn` then names a column of its `var`. Only a MuData has one. Absent: the first gene-expression modality (`rna`, `gex`, `gene expression`, `gene_expression`, `sct`, any case), else the first in `mod-order`. | `modality` |
| `expression` | `'X'`, `'raw/X'` or `'layers/<name>'` | The matrix genes are read from. (Not written by `prepare()`, which stores only `X`.) | — |
| `geneIdColumn` | text (a `var` column) | The column whose values name the genes. | `var_gene_id_column` |
| `latentKey` | text (an `obsm` key) | The latent space the outlier filter measures in (never a 2D/3D layout). | `latent_key` |
| `vectorFieldDefault` | text (a field id, `<field>` of an `obsm['<field>_<basis>']` key: `velocity`) | The vector field the overlay opens on. | `vector_field_default` |
| `obsKeys` | list of text | The obs columns the viewer lists, in this order. | `obs_keys` |
| `normalizeEmbeddings` | boolean | `false` shows embeddings at their own coordinates. | `normalize_embeddings` |
| `centroidMinPoints` | positive integer | The fewest cells a category needs to get outlier quantiles. | `centroid_min_points` |

A server may pass the same keys for one session (`cellucid serve --embedding-2d`,
`serve_anndata(embedding_keys=...)`); they override the file's own.

There is no `defaults` key for how an expression matrix is transformed or which
genes are highly variable. Those are Scanpy's own records, which the viewer
reads from any AnnData: `uns['log1p']` (written by `sc.pp.log1p`; its `base`,
none for e) makes `X`, every layer and a non-integer `raw.X` open as `log1p`,
and a true/false `var['highly_variable']` (written by
`sc.pp.highly_variable_genes`) is the highly variable genes. A second record
of the same facts here could disagree with them. `prepare()` and
`cellucid_prepare()` write neither: `prepare()` keeps the `var` columns it is
given (a `highly_variable` column among them) but writes no `uns['log1p']`,
so a prepared file opens with its values as stored, which the Data mapping
panel can change.

(outlier-quantiles-spec)=
### `outlierQuantiles`

Written by both writers whenever `latent_key` is given, so the web app's
outlier filter reads a few MB instead of the whole latent space. For each
categorical obs field, each cell's quantile is the share of its category that
lies at its distance from the category's centre in the latent space, or closer
(the rule, bit for bit: `cellucid.prepare_data._outliers`; the web app's
`assets/js/data/latent-outliers.js`; the numpy reference
`scripts/generate-latent-outlier-conformance.py` in the web repository):

1. over the latent values as stored (`float32`, widened exactly to `float64`),
   a category's centre is its cells' column sums, accumulated cell by cell in
   cell order, divided by its cell count;
2. a cell's distance is the square root of its squared differences from the
   centre, accumulated column by column in column order;
3. its quantile is the number of the category's cells at that distance or
   closer divided by the cell count, both as `float32`.

Only categories with at least `minPoints` cells are measured; every other cell
(a missing value, a smaller category) is `NaN`.

| Key | Type | Meaning |
|---|---|---|
| `latentKey` | text (an `obsm` key) | The latent space they were computed in: `defaults.latentKey`. |
| `minPoints` | integer (`<i8`) | The fewest cells a category needed: `centroid_min_points`, else 10. |
| `fields` | mapping | `{<obs column>: <f4 array of n_obs>}`, one array per field, in obs column order. |

Each field is an ordinary `array` element (`<f4`, chunks of 2^18 values,
byte-shuffled and deflated like every array read whole, §3), so a reader
reads one field without the others. A field is every obs column the viewer
shows as categorical — a categorical, a text column, a boolean column or a
nullable boolean column, partitioned by its values — that has at least one
category of `minPoints` cells; any other column has no entry. On the
561,947-cell Suo atlas each field is 1.8 MB stored (14 fields, 25 MB, against
the 54 MB of its 24-column latent space read whole for every computation).

The web app uses them when the latent space in use is `latentKey`, its own
fewest-cells rule is `minPoints`, the array has one value per cell, and the
values fit the field's categories as the file holds them now: in every scored
category of `n` cells each value is `r / n` for a whole `r` from 1 to `n`,
exactly `r` cells have a rank of at most `r`, and every other cell is `NaN`.
Otherwise (another latent space, a file whose obs was changed after it was
written) the web app computes them from the latent space, with the same bytes.
Other AnnData tools ignore the mapping, as they ignore all of `uns['cellucid']`.

In Python:

```pycon
>>> adata.uns["cellucid"]
{'schema': 1, 'id': 'pbmc-3k', 'name': 'PBMC 3k',
 'citation': {'text': 'Doe et al. 2026', 'url': 'https://example.org/paper'},
 'defaults': {'embedding2d': 'X_umap',
              'embedding3d': {'key': 'X_pca', 'columns': array([0, 1, 2])},
              'latentKey': 'X_pca'},
 'outlierQuantiles': {'latentKey': 'X_pca', 'minPoints': 10,
                      'fields': {'leiden': array([0.62, 0.15, ..., 0.98], dtype=float32)}}}
```

---

## 3. Layout for the browser

The web app reads the file itself with HTTP range requests (from disk, from
`cellucid serve`, or from any web host): only the bytes a view needs cross the
wire, and every chunk is one request and one decode. These are the writers'
choices, measured with the web app's own reader; a file laid out differently is
still read, only more slowly. Both writers make the same choices for the Zarr
formats (2 and 3); the `.h5ad` rows are the Python writer's alone, since the R
package writes no HDF5.

### Chunks

| Array | Chunk |
|---|---|
| `X/data`, `X/indices` (read one gene's range at a time) | 65,536 values (2^16) |
| `X/indptr`, obs and var columns and indexes, every array of `obsp/connectivities`, `uns` lists (read whole) | 262,144 values (2^18) |
| `obsm` arrays (read whole, by rows) | `[rows, n_columns]`, rows = 262,144 / n_columns (at least 1) |
| dense `X` (read one gene at a time) | `[rows, genes]`: rows up to 262,144, genes as many as keep the chunk at 262,144 values (at least 1) |

An array shorter than its chunk size is one chunk of exactly its length; an
empty array has chunks of 1.

### Codecs

| Container | Codec |
|---|---|
| `.h5ad` | HDF5 deflate (gzip) level 6 on every chunked array, after HDF5's shuffle filter on every array of values wider than one byte except the expression values; scalars uncompressed, as HDF5 requires (Python writer). |
| Zarr format 2 | `{"id": "zlib", "level": 6}` on every array, `uns` scalars and empty arrays included, after the filter `{"id": "shuffle", "elementsize": <bytes>}` on every array of values wider than one byte except the expression values; fill value 0, `false` or `""` |
| Zarr format 3 | `{"name": "blosc", "configuration": {"typesize": <bytes>, "cname": "zlib", "clevel": 6, "shuffle": "shuffle", "blocksize": 0}}` on every array Zarr format 2 shuffles, `{"name": "gzip", "configuration": {"level": 6}}` on every other array; only codecs of the Zarr v3 specification (and `vlen-utf8` for strings); unsharded |

**Shuffle.** A byte shuffle stores the first byte of every value, then every
second byte, and so on, so the similar high bytes of neighbouring values
compress together: on the Suo atlas the gene indices of `X` shrank by a third
and the embeddings and graph weights by 7–9%. The expression values (a sparse
`X`'s `data`, a dense `X`) are the exception: real expression values repeat
whole, which deflate finds only unshuffled — shuffled, the gene-major values
of three atlases grew 14–45% as raw counts, after `normalize_total` + `log1p`,
and as proportions alike. Strings and one-byte values (`|b1`, `|i1`) are not
shuffled either.

Zarr format 2 shuffles with numcodecs' `shuffle` filter (a format 2 filter is
any numcodecs codec). Zarr format 3 has no shuffle codec of its own:
`numcodecs.shuffle` is outside its specification, so another Zarr
implementation need not read it. A format 3 array is shuffled by the
specification's `blosc` codec instead: Blosc byte-shuffles each block of a
chunk (the whole chunk, or blocks of 1 MiB in a larger one) and stores each
byte plane as a zlib stream at level 6. On the He atlas (71,650 cells) the
format 3 store is 29.00 MB so, against 29.11 MB with `numcodecs.shuffle` and
gzip and 37.62 MB with gzip unshuffled.

Browsers inflate deflate, zlib and gzip natively (`DecompressionStream`), the
zlib streams of a Blosc frame included (the web app reads the frame and
un-shuffles in JavaScript); the web app also decodes Blosc's other
compressors, LZ4 and Zstandard in JavaScript, so a file written with those is
read too.

### Consolidated metadata

A Zarr store carries consolidated metadata — `.zmetadata` in format 2, the
`consolidated_metadata` member of the root `zarr.json` in format 3 — so a
reader on a web host learns every array and group with one request; a browser
cannot list a folder. (`cellucid serve` builds it in memory for a store without
it.)

### Where things sit in the file

A single file (`.h5ad`, `.zarr.zip`) is read in ranges, and ranges that lie
close together are fetched as one request. The writers place what is read
together side by side:

- **`.h5ad`: what opening reads first.** HDF5 gathers the file's metadata —
  object headers, attributes, chunk indexes — in blocks of 256 KiB
  (`meta_block_size`), the first at the start of the file (the rebuilt Suo's
  metadata — 5,103 genes, stored outlier quantiles — is 55 KB). Every group and array is declared, with its attributes, before any
  value is written, so the strings of the attributes (every `encoding-type`),
  which HDF5 keeps in global heap collections placed among the values written
  so far, follow the metadata block. Then come the values opening reads —
  `uns` but its stored outlier quantiles, `var`'s index (the gene names),
  `X/indptr` and the labels (`categories`) of every categorical `obs` column —
  then every other value in element order, each array's chunks together, the
  bulk of `X` last. Suo's opening region is the first 431 KB: the web app's
  first request of a file is its first 512 KiB, sent with the probe of its
  size, so the file opens one round trip after the link is followed. File
  format 1.10 (`libver`), which HDF5 1.10 and later — h5py, anndata, rhdf5,
  hdf5r, h5wasm — read, keeps a one-chunk array's chunk address in its header
  and indexes a fixed-size array's chunks in one fixed array, so no chunk
  waits for an index read. The file reads the same as one anndata writes
  (anndata 0.8 through 0.12, scanpy, h5py: every element equal, attribute for
  attribute). An `.h5ad` written with HDF5's defaults scatters its metadata in
  2 KiB blocks between the arrays, and opening it waits one round trip per
  level of its tree; one written by anndata keeps the strings of most
  attributes among `obs`.
- **`.zarr.zip`: the chunks first, the metadata last.** Each array's chunks
  follow one another in grid order (`X/data/9` before `X/data/10`), then come
  the metadata documents, the consolidated metadata the very last entry. A
  reader opening an archive reads its end — the central directory — and with
  it the documents that describe the store, so the store opens with that one
  read.

### `.zarr.zip`

The store folder, with the store at the top level of the archive (`X/data/0`,
…, `.zgroup`, `.zmetadata`; no enclosing folder). Every entry is **stored**
(no ZIP compression): the chunks are compressed already, and a stored entry is
read with one range request at its offset.

### Measured choices

Measured on 23 October 2026 on the Suo atlas (561,947 cells × 8,192 genes,
raw counts) and the He atlas (71,650 cells), with
`scripts/benchmark-anndata.mjs` of the web repository — the web app's reader
opening the file, reading every `obs` field, the embeddings, a first and a next
gene, 99 genes at once, 1,000 cells, the neighbour graph and every gene in the
cells of the two largest cell types — from local files and over simulated
links (0, 20, 80 and 150 ms round trips; HTTP/1.1 with six requests at a time,
HTTP/2 with GitHub's 100), and with the reader in a Chromium worker and its
decode pool for local files. Each step was compared with its bound: one round
trip plus the fewest bytes any layout fetched for it over the link's
bandwidth. Choices were made on the worst step against its bound, not on
totals. Suo, before → after:

| | `.h5ad` | `.zarr` folder | `.zarr.zip` |
|---|---|---|---|
| Size | 230.3 → 227.7 MB | 243.4 → 198.1 MB | 243.6 → 198.3 MB |
| Open, 80 ms link | 542 → 246 ms (16 → 5 requests) | 192 → 194 ms | 206 → 127 ms (3 → 2 requests) |
| First / next gene, 80 ms link | 169 / 165 → 85 / 83 ms | 87 / 86 → 86 / 83 ms | 85 / 85 → 86 / 83 ms |
| Every step, HTTP/2, 80 ms, 50 MB/s | 8.9 → 8.3 s; worst step 6.7× → 3.0× its bound | 8.3 → 8.0 s; 2.4× → 2.4× | 9.2 → 8.2 s; 2.5× → 2.1× |
| Every step, HTTP/2, 150 ms, 5 MB/s | 48.5 → 46.9 s; 11.1× → 5.0× | 59.6 → 46.1 s; 2.2× → 2.1× | 59.9 → 46.7 s; 2.5× → 1.5× |
| Every step, HTTP/1.1, 80 ms, 50 MB/s | 9.9 → 9.0 s; 6.7× → 4.4× | 47.3 → 46.8 s; 15.2× → 15.1× | 12.1 → 9.0 s; 4.5× → 4.7× |
| Every step, local file in Chromium | 2.6 → 2.5 s | 2.6 → 2.6 s | 2.7 → 2.7 s |

What each change bought, and what was measured and not chosen:

- **Byte shuffle** (every array but the expression values): the Zarr store
  18.6% smaller; every read of those arrays fetches fewer bytes (the graph
  76.7 → 66.5 MB; on the 5 MB/s link the Zarr total fell 59.6 → 46.1 s).
  Shuffling the expression values too made them 14–45% larger on real data.
- **`.h5ad` metadata at the start** (a 256 KiB metadata block, file format
  1.10, opening's values written first): opening 16 → 5 requests, 542 →
  246 ms; a gene one round trip instead of two. HDF5's paged file space
  (pages of 64 KiB to 4 MiB) gathered the metadata too but spread each array's
  chunks over pages shared with the others: a whole-matrix read over HTTP/1.1
  took 224 requests instead of 33 (3.3× its bound instead of 1.6×), and the
  file grew 1–21%.
- **`.h5ad` attributes and opening's values in its first 512 KiB** (every
  element declared before any value is written; the `obs` labels with what
  opening reads), read with the probe of the file's size: the rebuilt Suo's
  (5,103 genes) probe and opening over 80 ms at 50 MB/s 327 → 98 ms (7 requests in three dependent
  waves → the probe and one read, in one round trip), 3.9× → 1.2× its bound
  (one round trip and the 208 KB opening needs); over 20 ms at 100 MB/s 4.4×
  → 1.5×, over 150 ms at 5 MB/s 5.5× → 1.5× (the rest of 512 KiB is read
  with what opening needs). In the app in Chromium the dataset's first
  embedding request comes 0.81 → 0.22 s after the link is followed, and the
  dataset is ready 1.9 → 1.3 s after it.
- **`.zarr.zip` metadata last, chunks in grid order**: opening 3 → 2 requests;
  a whole-matrix read over HTTP/1.1 4 requests instead of 56 (entries had been
  in name order, `0, 1, 10, 100, …`).
- **Alternating the chunks of `X/data` and `X/indices`** (a gene's two ranges
  side by side): 99 genes up to a fifth faster over HTTP/1.1, but reading one
  array alone then fetched the other too — 1,000 cells over HTTP/1.1 read
  113 MB instead of 74 and took 6.3 s instead of 4.9. Not chosen.
- **Chunks of `X`**: 2^15 values doubled the requests of a whole-matrix read
  (over HTTP/1.1, 41 s instead of 21 s for the two cell types); 2^17 and 2^18
  cut a folder's whole-matrix reads over HTTP/1.1 (8.5× and 4.1× their bound
  instead of 17×) but fetched 1.7× and 3× the bytes of 99 genes, which on a
  slow link costs time (2.3× and 3.9× their bound instead of 1.4×), and made
  each gene fetch 0.2–0.4 MB. 2^16 is the
  smallest worst step over HTTP/2 and in every single file; a folder over
  HTTP/1.1 is not a layout to choose (§4).
- **Codecs**: on shuffled data deflate level 6 was the smallest of everything
  tried — Zstandard (level 3) 3% larger, Blosc-Zstandard 2%, Blosc-LZ4 30% —
  and the others are decoded in JavaScript rather than natively; no read was
  faster. Deflate level 9 made the expression values 6% (everything else
  under 1%) smaller and wrote them twenty times slower.
- **Zarr v3**: read exactly as fast as v2 (8.0 against 7.8 s over HTTP/2);
  v2 stays the default because anndata 0.8–0.10 and zarr-python 2 read only
  v2. A sharded v3 store (64 chunks per shard) cut a folder's whole-matrix
  requests over HTTP/1.1 (7.8× instead of 15×) but every first gene of a
  shard waited for the shard's index: 351 ms instead of 87 over HTTP/2.
- **A cell-major copy of `X`** (CSR beside the CSC, for reads by cell): the
  web app reads every cell subset by gene columns, each byte once; a copy
  would add 36% to the file (71 MB on Suo) and spare at most about half of a
  read by cell (a median cell type's rows still touch 43% of the copy's
  chunks at 2^16 values, the largest cell type's 57%, 1,000 random cells
  68%). Not written.
- Kept as they were, measured: consolidated metadata (a folder opens in one
  request for its structure), `float32` embeddings, and the category code
  widths pandas gives (`|i1` for fewer than 127 categories).

---

## 4. Writing

Both writers take the output format as a required argument, with no default,
because each suits a different place and only the caller knows where the file
goes; the output path must end in that format's suffix, and a mismatch is
refused with the path to write instead.

| `format` | Writes | For | Writer |
|---|---|---|---|
| `"zarr"` | a `.zarr` folder of small files | GitHub, a repository or its Pages site (no file over 100 MB, no LFS; HTTP/2) | Python and R |
| `"h5ad"` | one `.h5ad` file | everything else: single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local files, `cellucid serve`, Jupyter | Python only: the R package writes no HDF5, so it depends on no HDF5 library |
| `"zarr.zip"` | the same Zarr store as one `.zarr.zip` file | wherever an `.h5ad` goes; the R package's single file | Python and R |
| `"h5mu"` | one `.h5mu` file: a Cellucid-ready MuData (§5) | a MuData wherever an `.h5ad` goes | Python only (R writes a MuData as Zarr) |

Both writers take `zarr_format` (2 or 3; 2 by default) for the Zarr formats.
Why each format goes where, with the measurements:
{ref}`which-format`.

In Python, `format=` is a keyword argument (`prepare(path, format="zarr", ...)`);
in R it follows `out_path` (`cellucid_prepare(path, format = "zarr", ...)`) and
is matched as `match.arg()` matches.

The output is written beside its target under a temporary hidden name
(`.<name>.<process id>.partial`) and moved into place only once complete, so an
interrupted or failed write never leaves a half-written file where a viewer
would open it. An existing output is replaced only with `force=True` (R:
`force = TRUE`), and a folder that is not a Zarr store is never replaced.

(cellucid-ready-mudata)=
## 5. Cellucid-ready MuData

Several feature spaces over the same cells -- RNA and surface proteins
(CITE-seq), RNA and chromatin (multiome) -- are one **Cellucid-ready
MuData**: a standard MuData (the [MuData on-disk
specification](https://github.com/scverse/mudata/blob/main/docs/io/spec.md)) that
`mudata.read_h5mu` / `mudata.read_zarr` read, whose modalities are each a
Cellucid-ready AnnData. {func}`~cellucid.prepare` writes one from a MuData
(`data=mdata`) as an `.h5mu` file (`format="h5mu"`), a MuData Zarr folder
(`"zarr"`) or that folder zipped (`"zarr.zip"`); `cellucid_prepare()` in R
writes the Zarr forms from a multi-assay object
({ref}`its layout <r-cellucid-ready-mudata>`). The web app reads all three.

| Element | Contract |
|---|---|
| root attributes | `encoding-type` `MuData`, `encoding-version` `0.1.0`, `axis` `0`: the modalities share the cells. A MuData whose modalities share features (`axis` 1 or -1) is not a Cellucid-ready MuData. |
| `.h5mu` user block | 512 bytes before the HDF5 signature, starting `MuData (format-version=0.1.0;creator=cellucid;...)` as mudata writes it. |
| `obs` | The cells, one row each, unique names: the global `obs`. A column mudata copies from a modality (`<name>:<column>`) is not stored again; the viewer shows each modality's own columns under that name, missing where the modality did not measure the cell. |
| `obsm` | The shared embeddings and vector fields, laid out as an AnnData's (§1); and mudata's per-modality masks: `obsm/<name>`, boolean `(n_cells, 1)`, which cells the modality measured (the viewer shows them as `<name>:measured`). |
| `obsp/connectivities` | The shared neighbour graph, as in an AnnData (§1). |
| `uns['title']`, `uns['cellucid']` | At the root only, as in an AnnData (§2); `defaults` may name the `modality` the viewer opens on. |
| `var` | Every modality's feature names, in `mod-order`: mudata's global `var` (index only). |
| `varm`, `varp` | Present and empty. |
| `obsmap/<name>`, `varmap/<name>` | For every modality: `uint32`, one entry per global cell (feature), the 1-based row of the modality that holds it, `0` where the modality has none. |
| `mod` | The attribute `mod-order` lists the modalities; each `mod/<name>` is a Cellucid-ready AnnData (§1, §3) without `uns['cellucid']` and `uns['title']`: gene-major `X` (required), its `var` and `obs`, its own embeddings in `obsm` (shown as `<name>:<key>`) and `obsp/connectivities`. |

Modality names follow the element-name rule of §1. Cell names are unique
across the MuData and within every modality: MuData places a modality's cells
among the global ones by name. A modality has at least one cell and one
feature. The launch defaults are checked against what the viewer shows: the
cells' `obs` (with every modality's columns), `obsm` (with every modality's
arrays, `<name>:<key>`), and the `var` of the modality it opens on;
`latentKey` names a root embedding or that of a modality that measured every
cell.

`prepare()` writes the modalities one after the other and converts each
modality's matrix only while it is written, so its memory beyond the MuData
itself is one modality's converted matrix.

## Next

- {doc}`12_prepare_from_anndata_or_mudata`
- {doc}`01_prepare_overview`
- {doc}`10_performance_and_storage`
