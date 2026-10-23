# Input requirements (global rules)

**Audience:** everyone preparing data
**Time:** 10 minutes
**Goal:** the rules every input of `cellucid_prepare()` shares, so a call succeeds the first time.

## One cell axis

The number of cells is the number of rows of the embeddings. Every per-cell
input must have exactly that many rows, in the same cell order:

| Input | Shape |
|---|---|
| each `embeddings` matrix | cells x columns |
| `obs` | cells x fields |
| `gene_expression` | cells x genes |
| `connectivities` | cells x cells |
| each `vector_fields` array | cells x 1, 2, or 3 |

`gene_expression` is **always cells x genes**. Seurat and SingleCellExperiment
store genes x cells, so transpose with `Matrix::t()` first. The function checks
the row count but never guesses the orientation: a square matrix in the wrong
orientation cannot be detected from its shape.

## Cell names

The obs index (the cell names) comes from the row names your inputs carry:

- `rownames(obs)`, when you set them (a data frame's automatic `"1"`, `"2"`, ...
  row names do not count),
- `rownames(gene_expression)`,
- the row names of each embedding matrix.

Where several are set they must be **the same names in the same order**,
because a different order means the rows describe different cells:

```text
The row names of gene_expression and obs disagree at row 2 ('q' and 'y'). Every input must list the same cells in the same order.
```

With no names at all the cells are named `"0"` to `"n - 1"`, as anndata names
them. Repeated cell names are made unique the way anndata's
`obs_names_make_unique()` does it (the second `"AAAC"` becomes `"AAAC-1"`),
with a message saying what changed.

## Gene names

The var index comes from `rownames(var)` (when set), else
`colnames(gene_expression)`; if both are set they must agree. See
{doc}`05_var_gene_metadata`.

## Numbers

Every matrix -- embeddings, expression, vector fields, the graph -- is stored as
**float32**. Values must be finite (no `NA`, `NaN`, or `Inf`) and within the
float32 range (magnitude up to about 3.4e38); a larger value would become an
infinity and is refused. `obs` and `var` columns keep their own types
({doc}`04_obs_cell_metadata`).

A data frame is accepted in place of a matrix when every column is a plain
numeric vector. A classed numeric column (for example `units::units` or
`bit64::integer64`) is refused, because converting it to a matrix drops the
attribute that says what its numbers mean.

## Names that become paths

Every `obs` and `var` column name and every `obsm` key (embedding and vector
field) becomes a folder or file name in the store, and a `.zarr.zip` may be
unpacked anywhere. Each must therefore be a portable file name on every
platform, whichever platform writes it:

- no `/`, `\`, `<`, `>`, `:`, `"`, `|`, `?`, `*`, or control character;
- no trailing space or dot (Windows drops them);
- not a Windows device name (`CON`, `PRN`, `AUX`, `NUL`, `COM1`-`COM9`,
  `LPT1`-`LPT9`, also with an extension);
- not `.`, `..`, `.zarray`, `.zattrs`, `.zgroup`, or `.zmetadata`;
- for a data frame column, not `_index` (AnnData keeps the row names there);
- different from its neighbours by more than letter case: `"Sex"` and `"sex"`
  are one folder on macOS and Windows.

A name that breaks a rule is refused with the rule and a way to rename it.

## Dataset text

`dataset_name`, `dataset_description`, `citation_text`, `citation_url`,
`source_name`, and `source_url` are shown to readers verbatim. Each must be one
string free of control characters (U+0001-U+001F, U+007F-U+009F), of the
zero-width characters U+200B, U+2060, and U+FEFF, and of leading or trailing
whitespace. Nothing is trimmed for you. `dataset_name` must not be empty;
`dataset_description` may be `""`. A URL must be an absolute `http://` or
`https://` address, and is accepted only beside its text: `citation_url` needs
`citation_text`, `source_url` needs `source_name`.

`dataset_id` must be a portable identifier: 1-180 ASCII letters, digits, `.`,
`_`, or `-`, starting with a letter or digit, not ending with `.`, and not a
Windows device name. Sessions and community annotations are keyed on it, so
keep it stable across versions of the same dataset.

## Output path

- `format = "zarr"`: `out_path` ends in `.zarr` and names a folder.
- `format = "zarr.zip"`: `out_path` ends in `.zarr.zip` and names a file.
- `format` has no default, and a path that does not end in its suffix is
  refused with the path to write instead.
- If `out_path` exists, the call stops unless `force = TRUE`; even then an
  existing folder is replaced only if it is a Zarr store (it holds `.zgroup`).

## Next

- {doc}`03_embeddings_and_coordinates`
- {doc}`11_troubleshooting_prepare`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/02_input_requirements_global`
