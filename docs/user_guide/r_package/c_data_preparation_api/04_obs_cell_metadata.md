# `obs`: cell metadata

**Audience:** everyone
**Time:** 10 minutes
**Goal:** know how each kind of `obs` column is stored and what the viewer makes of it.

## What to pass

`obs` is an optional `data.frame` with **one row per cell**, in the cell order
of the embeddings. Every column is stored. Without `obs` the store holds an
`obs` data frame with the cell names and no columns.

```r
obs <- data.frame(
  cell_type = factor(cell_type),
  sample = sample_id,          # character, repeats
  n_counts = n_counts,         # integer
  pct_mito = pct_mito,         # double
  is_doublet = is_doublet,     # logical
  row.names = cell_barcodes
)
```

Row names you set become the cell names (they must agree with the row names of
`gene_expression` and the embeddings; see
{doc}`02_input_requirements_global`). Column names become paths in the store
and follow the path rules on that page; `_index` is reserved.

## How each column is stored

The rules are anndata's own, so Python reads the same column types you had in R:

| R column | Stored as (anndata encoding) | In Python |
|---|---|---|
| factor | `categorical`: integer codes plus its levels, **in their order** | `pandas.Categorical` |
| ordered factor | `categorical` with `ordered: true` | ordered `Categorical` |
| character with a repeated or missing value | `categorical`, levels in code-point order | `Categorical` |
| character with every value distinct and none missing | `string-array` | object strings |
| logical without `NA` | `array`, bool | `bool` |
| logical with `NA` | `nullable-boolean` (values + mask) | `boolean` |
| integer without `NA` | `array`, int32 | `int32` |
| integer with `NA` | `nullable-integer` (values + mask) | `Int32` |
| double | `array`, float64; `NA` and `NaN` stored as NaN | `float64` |

Every column keeps its own type: an integer column is int32 because R's
integer is int32, and a double column is float64 even when every value is a
whole number. Python's `cellucid.prepare()` follows the same rule, so the same
column written from R or from Python is the same element.

The character rule is what anndata's `strings_to_categoricals()` does on write:
a column that repeats (a sample, a batch) becomes a categorical, and a column of
unique strings (a barcode copy) stays a string array.

A missing value in a factor or a character column is stored as code `-1` (no
category). Category codes take the integer type pandas gives them: int8 for
fewer than 127 categories, int16 for fewer than 32,767, int32 beyond. There is no limit on the
number of categories.

## What is refused

- **A classed column with no AnnData encoding**, such as `Date`, `POSIXct`,
  `difftime`, or a list column. Convert it first:

  ```text
  obs column 'day' is a Date column, which has no AnnData encoding Cellucid writes. Convert it first, for example with as.character() or as.numeric().
  ```

- **An infinite value** in a double column (a colour scale has no position for
  it); use `NA` for a missing value instead:

  ```text
  obs column 'score' holds 1 infinite value (first at row 2). A colour scale has no position for an infinity; replace it with NA, which is written as a missing value.
  ```

- **A different number of rows** than the embeddings:
  `obs has 2 rows, but the embeddings have 3 cells.`

Both column checks happen while the store is being written, so the call stops
with nothing at `out_path` (see {doc}`01_cellucid_prepare_overview`).

## Choosing what the viewer shows

`obs_keys` names the obs columns the viewer lists, in this order; every column
is stored either way:

```r
obs_keys = c("cell_type", "sample", "pct_mito")
```

It is stored as `uns['cellucid']['defaults']['obsKeys']`. Each name must be a
column of `obs`, at most once.

## Tips

- Turn a numeric cluster id into a factor (`factor(seurat_clusters)`) if you
  want it coloured by category rather than on a continuous scale.
- Set factor levels in the order you want the legend to show.
- Category labels are stored exactly as given: `"T cell"` and `"T cell "` are
  two categories.

## Next

- {doc}`05_var_gene_metadata`
- {doc}`06_gene_expression_matrix`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/04_obs_cell_metadata`
