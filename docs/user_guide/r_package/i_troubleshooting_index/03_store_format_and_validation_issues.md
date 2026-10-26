# Store Format and Validation Issues

**Audience:** everyone with a written store that does not look right  
**Time:** 10 minutes  
**What you'll get:** how to check a store, and what the common surprises mean

`cellucid_prepare()` writes a standard AnnData Zarr store, so the quickest
validation is to read it with anndata. The layout is specified in
{doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.

## Check a store with Python's anndata

```python
import warnings, anndata
warnings.simplefilter("error")          # any warning is an error
adata = anndata.read_zarr("my_dataset.zarr")
print(adata)                             # cells x genes, obs, var, obsm, obsp
print(adata.uns["cellucid"])             # id, name, citation, source, defaults
```

For a `.zarr.zip`:

```python
import zarr
adata = anndata.read_zarr(zarr.storage.ZipStore("my_dataset.zarr.zip", mode="r"))
```

## Check a store from R

```r
store <- "my_dataset.zarr"
list.files(store, all.files = TRUE, no.. = TRUE)   # .zgroup .zattrs .zmetadata X obs var obsm ...
jsonlite::read_json(file.path(store, "uns/cellucid/name/.zarray"))
utils::unzip("my_dataset.zarr.zip", list = TRUE)   # the entries of a .zarr.zip
```

## Symptom: “the folder looks empty”

The metadata files of a Zarr store (`.zgroup`, `.zattrs`, `.zarray`,
`.zmetadata`) start with a dot, so a file browser or `ls` hides them. Use
`ls -a` or `list.files(all.files = TRUE)`.

## Symptom: “there is no `X` / `var` is empty”

You did not pass `gene_expression`; the store then has no `X` and zero genes,
which is a valid AnnData. Pass `gene_expression` (cells x genes) to store
expression.

## Symptom: “gene search finds the Ensembl ids, not the symbols”

Genes are named by the `var` row names. To show another column, pass
`var_gene_id_column = "symbol"`; it is stored as
`uns['cellucid']['defaults']['geneIdColumn']`. See
{doc}`../c_data_preparation_api/05_var_gene_metadata`.

## Symptom: “a default I expected is not in `uns['cellucid']['defaults']`”

Only the launch defaults you passed are written. Without `embedding_defaults`,
`latent_key`, `obs_keys` and the rest, the viewer makes its own choice. See
{doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility`.

## Symptom: “a column is categorical (or a string array) when I did not expect it”

The encodings follow anndata: a factor is a categorical; a character column is
a categorical when a value repeats or is missing, and a string array when every
value is unique; logical and integer columns with `NA` become nullable. Convert
the column before calling. See
{doc}`../c_data_preparation_api/04_obs_cell_metadata`.

## Symptom: “there is a hidden `.<name>.<pid>.partial` next to my store”

The store is written there first and renamed into place when complete. A
leftover one means R was killed during the write; it is not a store, so delete
it.

## Symptom: “the `.zarr.zip` does not open”

- `cellucid_prepare(format = "zarr.zip")` writes the store at the top level of the
  archive, which is what anndata's `zarr.storage.ZipStore` expects. A ZIP you
  make yourself around the `.zarr` folder puts the store one level down, which
  anndata cannot open; write it with `format = "zarr.zip"` instead.
- Its entries are stored, not compressed again; that is expected (each chunk is
  already zlib-compressed).

Next: {doc}`../d_viewing_loading/03_validate_stores_and_debug_loading`.
