# `prepare()` writer architecture

This page explains how `cellucid.prepare(...)` works internally: what it
writes, why each step exists, and what breaks if you change it.

User-facing guide: {doc}`../c_data_preparation_api/index`. The contract the
output follows, shared with cellucid-r and the web app:
{doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

:::{admonition} Audience
:class: note

- Computational users: focus on “inputs, outputs, and invariants”.
- Developers: read the whole page before changing the writer or the contract.
:::

---

## What `prepare(...)` produces

One standard AnnData — `.h5ad`, a `.zarr` store folder (format 2 or 3), or a
`.zarr.zip` — that `anndata` reads without a warning, laid out for the web
app's own HDF5/Zarr reader:

- `X` gene-major (CSC) float32, chunked so one gene is one small range read;
- `obs`, `var` with anndata's encodings (repeated text as categoricals);
- `obsm` embeddings and `<field>_<basis>` vector fields (scVelo's name for the
  field on `X_<basis>`), float32;
- `obsp['connectivities']`, symmetric CSR float32;
- `uns['cellucid']` (schema, id, name, description, citation, source, launch
  defaults) and `uns['title']`;
- with `latent_key`, `uns['cellucid']['outlierQuantiles']`: every categorical
  `obs` column's latent-space outlier quantiles (`{latentKey, minPoints,
  fields: {<column>: float32[n_obs]}}`, NaN where a cell has none), so the
  outlier filter reads one field's array instead of the whole latent space.

Nothing else is precomputed: the browser normalizes embeddings and computes
the centroids it draws itself, and computes the outlier quantiles of a file
that stores none.

---

## Where the code is

`src/cellucid/prepare_data/`:

| Module | Responsibility |
|---|---|
| `_prepare.py` | {func}`~cellucid.prepare`: the required `format` (a call without it is refused with the choices), the output target, argument checks, the launch defaults, `uns['cellucid']`, then one call to the writer |
| `_inputs.py` | every input check and conversion: element keys, embeddings, vector fields, `obs`/`var` columns, cell and gene names, the expression matrix and its dtype, the graph, linked text |
| `_outliers.py` | the latent-space outlier quantiles stored with `latent_key`: every categorical, text, boolean and nullable boolean `obs` column, by the rule the web app and cellucid-r compute to the same bytes |
| `_store.py` | the element writer (`h5py` or `zarr`), chunk sizes, codecs, consolidated metadata, the `.zarr.zip` packing, and the staged, atomic placement |

Shared modules it uses: `_launch_defaults.py` (the same launch-default check
`serve_anndata` runs, producing the web app's key names), `_contracts.py`
(dataset id and display-text rules), `_non_finite.py` (the counted NaN/Inf
diagnosis), `data_mapping.py` (making names unique, the shared rule).

---

## Pipeline stages

1) **Resolve the target** (`_resolve_target`): `format` → container; refuse a
   value outside `"h5ad"`, `"zarr"`, `"zarr.zip"`, an `out_path` whose suffix
   disagrees with it (naming the path to write instead), an existing path
   without `force=True`, a non-Zarr folder even with it, and a `zarr_format`
   with `format="h5ad"`.
2) **Check identity and text**: `dataset_id` (portable id), `dataset_name`,
   description, citation and source (display text; URLs absolute `http(s)`).
3) **Resolve inputs** (`_inputs.py`): embeddings → float32 `obsm`; vector
   fields → `obsm['<field>_<basis>']` for an embedding `X_<basis>` of the same
   width, each key checked to read back as its field and embedding
   (`resolve_vector_fields`, `parse_vector_field_key`); `obs` → stored columns + unique cell names; expression →
   CSC (or dense) float32, float64 for integer counts above 2^24; `var` + gene
   names (+ `gene_identifiers` subset); graph → symmetric CSR float32.
4) **Launch defaults** (`launch_defaults` in `_launch_defaults.py`), checked
   against an `AnnDataContents` built from the resolved elements, so `prepare`
   and the server refuse the same mistakes with the same messages.
5) **Write** (`write_store`): into a hidden sibling
   (`.<name>.<pid>.partial`), then `os.replace` into place; an existing output
   is moved aside first and removed only after the new one is in place.

Every input is checked before stage 5 creates the first path, so a refusal
leaves nothing on disk.

---

## The writer (`_store.py`)

`_Writer` creates every group and array itself, with anndata's encoding
attributes, so it controls each array's chunks and codec:

- numeric arrays: `create_dataset` (HDF5: gzip level 6 after the shuffle
  filter) or `create_array` (Zarr: numcodecs' `Shuffle` filter for a value of
  more than one byte, then v2 `Zlib(6)` / v3 `GzipCodec(6)`, the shuffle as
  the `numcodecs.shuffle` codec), chunked by `SLICED_CHUNK_ELEMENTS` (2^16,
  `X/data` and `X/indices`) or `WHOLE_CHUNK_ELEMENTS` (2^18, everything read
  whole); 2-D `obsm` arrays in row blocks across all columns; dense `X` in
  tall gene strips;
- nullable integer and boolean columns: a group of `values` and `mask`
  arrays, each written as the numeric arrays above;
- string arrays: in a Zarr store `anndata.io.write_elem` with the same chunk
  and codec options; in an `.h5ad` as anndata writes them (variable-length
  UTF-8 strings, `string-array` 0.2.0), declared first like every other
  array (below); strings are not shuffled;
- the arrays inside a sparse matrix carry no attributes, as anndata writes them;
- an `.h5ad` declares every group and array, with its attributes, before it
  writes any value (`_Writer`): HDF5 keeps a string attribute's value (every
  `encoding-type`) in a global heap collection placed among the values
  written so far, so declared first they sit right after the metadata. The
  values opening reads are written next -- `uns` but its stored outlier
  quantiles, the gene names (the `var` index), the pointers of a sparse `X`
  and the labels of every categorical `obs` column -- then the rest in element
  order (`obs`, `obsm`, `obsp`, the outlier quantiles, the bulk of `X`
  last), each array flushed as it is written (HDF5 places a chunk smaller than
  its chunk cache only when its dataset is flushed);
- an `.h5ad` is written in file format 1.10 with HDF5's metadata gathered
  in blocks of `H5AD_META_BLOCK_BYTES` (256 KiB): every object header and
  chunk index sits at the start of the file, so with the attributes' strings
  and opening's values the web app reads all opening needs with its first
  request of the file, its first 512 KiB (Suo's opening region is 431 KB);
  anndata 0.8 through 0.12, scanpy and h5py read the file exactly as one
  anndata writes;
- Zarr writes run with `anndata.settings.override(auto_shard_zarr_v3=False)`:
  one chunk is one file, read with one request;
- a Zarr store is consolidated (`zarr.consolidate_metadata`) so a browser
  learns the whole store from one document; a `.zarr.zip` is the finished
  folder packed with stored (uncompressed) entries at the archive root, every
  chunk first (each array's chunks in grid order), then the metadata
  documents, the consolidated metadata last, so the archive's end -- read
  when it opens -- describes the store.

The chunk and codec constants were chosen by measurement with the web app's
reader; the numbers are in the specification page. Change them only with new
measurements, and keep cellucid-r's writer in step.

---

## Invariants (do not break)

- **Standard AnnData.** No custom encodings: `anndata.read_h5ad` /
  `read_zarr` must read every output without a warning
  (`tests/test_prepare_anndata_contract.py` reads them with warnings as errors).
- **Browser-decodable.** Every array uses a codec in the server's
  `BROWSER_HDF5_FILTERS` / `BROWSER_ZARR_*` lists; the contract test asserts
  `browser_undecodable(...) == []` for all four containers.
- **Exact counts.** Integer counts read back exactly (float32 up to 2^24,
  float64 above).
- **`uns['cellucid']` keys** are exactly the web app's
  (`anndata-cellucid-metadata.js`, `requireAnnDataLaunchDefaults`), and only
  keys that were given are written.
- **Same as R.** cellucid-r's `cellucid_prepare()` writes the same contract;
  a change to element layout, names, or `uns['cellucid']` is a change to both
  writers, the web reader, and the specification page together.

---

## Tests

- `tests/test_prepare_anndata_contract.py`: round trips on `.h5ad`, Zarr v2,
  Zarr v3 and `.zarr.zip` (values, dtypes, encodings, chunks, codecs,
  consolidated metadata, `uns['cellucid']`), serving a prepared file without
  identity arguments, and every refusal leaving nothing on disk.
- `tests/test_docs_current_api_contract.py`: the documented signature and the
  specification's key tables against the code.

## Troubleshooting

Every refusal and its fix: {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
