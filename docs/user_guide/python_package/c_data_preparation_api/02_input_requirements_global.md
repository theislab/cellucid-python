# Input requirements (global)

**Audience:** computational users preparing their own arrays
**Time:** 10 minutes
**Goal:** know the rules every input obeys before `prepare()` writes anything

## One cell order

Every per-cell input — each embedding, `obs`, `gene_expression`,
`connectivities`, each vector field — has **one row per cell, in the same
order**. The first embedding sets the cell count; any other input with a
different count is refused with both counts named:

```text
obs has 40 rows, but the embeddings have 41.
```

`prepare()` cannot check that the rows *mean* the same cells, only that they
count the same. Pass every array from one AnnData (or subset them together)
so the order is the same by construction.

## Everything is checked first

Every input is validated before the first byte is written: a defect anywhere
fails the call and leaves nothing on disk, not even a partial file. The output
is then written beside its target under a temporary hidden name
(`.<name>.<process id>.partial`) and moved into place only when complete.

## Values

- Numbers must be finite and fit in float32. An expression matrix or a graph
  holding NaN, an infinity, or a magnitude above float32's range is refused
  with a count of what was found and the first positions affected.
- An embedding or a vector field may hold rows of NaN (or an infinity): a cell
  the array has no value for, written as given and not drawn in that
  embedding (or drawn without an arrow). Every other value fits float32, and
  at least one cell has finite values.
- `obs` and `var` numbers may hold NaN (a missing value) but not infinities.
- Integer counts in `gene_expression` are stored exactly: as float32 up to
  2^24, and as float64 when any count is larger.

## Names

- **Cells** are named by `obs.index`. Without `obs`, or with pandas' default
  `0 … n-1` index, they are named `"0"` … `"n - 1"`, as anndata names them.
- **Genes** are named by `var.index`; when `var` has the default index,
  `var_gene_id_column` names them.
- Repeated cell or gene names are made unique with anndata's rule
  (`MATR3`, `MATR3-1`) and the change is printed:

  ```text
  Gene names: 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.
  ```

- Index values must be text; convert other types first
  (`obs.index = obs.index.astype(str)`).

## Keys that become paths

`obs` and `var` column names and `embeddings` keys become paths inside the
store (folders of a Zarr store), so on every platform they must:

- not contain `/ \ < > : " | ? *` or a control character,
- not be `.`, `..`, or a Zarr metadata name (`.zarray`, `.zattrs`, `.zgroup`,
  `.zmetadata`, `zarr.json`),
- not end with a space or a dot,
- not be a Windows device name (`CON`, `PRN`, `AUX`, `NUL`, `COM1`–`COM9`,
  `LPT1`–`LPT9`),
- differ from their neighbours by more than letter case (`Sex` and `sex` are
  one folder on macOS and Windows).

A column named `_index` is refused: AnnData stores the row names under it.

## Dataset identity and text

- `dataset_id`: 1–180 ASCII letters, digits, `.`, `_` or `-`, starting with a
  letter or digit and not ending with `.`.
- `dataset_name`, `dataset_description`, `citation_text`, `source_name`: text
  the viewer draws exactly as stored — no control or zero-width characters, no
  leading or trailing white space. `dataset_description` may be empty.
- `citation_url`, `source_url`: absolute `http(s)` addresses, each only with
  its text (`citation_text`, `source_name`).

## Subsetting

`prepare()` writes what it is given. Subset cells before the call
(`adata[mask]`), pass only the embeddings and `obs` columns you want, and use
`gene_identifiers=[...]` to keep chosen genes in a chosen order
({doc}`05_var_gene_metadata`). `obs_keys` does not drop columns: every column
is stored, and `obs_keys` sets which ones the viewer lists, in order.

## Next

- {doc}`03_embeddings_and_coordinates`
- {doc}`11_troubleshooting_prepare`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/02_input_requirements_global`
