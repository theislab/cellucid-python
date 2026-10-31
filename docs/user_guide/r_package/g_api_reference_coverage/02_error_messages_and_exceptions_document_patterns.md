# Error Messages

**Audience:** everyone who hit an error, and maintainers  
**Time:** 10 minutes  
**What you'll get:** every message `cellucid_prepare()` and
`cellucid_transition_drift()` raise, grouped by the input it is about, with
what to change

## How errors are reported

- Every problem stops the call with an R error (`stop(..., call. = FALSE)`), so
  the message is the whole report; there are no warnings. (Of the viewer
  functions, only a failing `cellucid_on()` handler warns, and the other
  handlers still run.)
- Arguments are checked before anything is written. A problem found while
  writing (a column with no AnnData encoding, for example) happens inside the
  temporary store beside `out_path`, which is removed: a failed call never
  leaves a half-written store, and with `force = TRUE` the store it was going
  to replace is untouched.
- A message that lists values prints them as R code, for example
  `c("X_umap")`, so you can paste them back.
- Making repeated names unique is not an error; it is a `message()`, see the
  end of this page. So is each note naming a part of an object (`data =`)
  that is not written.

The messages below are copied from the package. Parts in angle brackets (`<n>`,
`<path>`) are filled in with your values.

## Required arguments

```text
cellucid_prepare() is missing 5 required arguments:
  - out_path: the path to write: a folder ending in ".zarr" for format = "zarr", or a file ending in ".zarr.zip" for format = "zarr.zip".
  - format: the kind of output to write, "zarr" (a .zarr folder of small files: for GitHub, no file over 100 MB, no LFS) or "zarr.zip" (the same store in one .zarr.zip file: for single-file hosts such as Zenodo, Figshare, Hugging Face, S3 or a lab server, local files, cellucid serve and Jupyter).
  - embeddings: a named list mapping each obsm key to a matrix with one row per cell, for example list(X_umap = umap).
  - dataset_id: one portable identifier of 1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a letter or number and not ending with '.'.
  - dataset_name: one non-empty name for the dataset, shown to the reader verbatim.
```

Only the arguments you left out are listed.

## `data`, `modalities`, `modality`

| Message | Fix |
|---|---|
| `data must be a Seurat, SingleCellExperiment or SummarizedExperiment object; got <a class> object. Pass matrices and data frames as the other arguments instead (embeddings, obs, gene_expression, ...).` | Pass an object of one of those classes, or the matrices themselves. |
| `data is a <class> object from the <package> package, which is not installed, so it cannot be read. Install it with <the install line>.` | Install the package it names. |
| `data is a Seurat v2 object (class 'seurat'), which Cellucid does not read. Update it first: data <- SeuratObject::UpdateSeuratObject(data).` / `This Seurat object was made by Seurat <version>, which Cellucid does not read. ...` | Update the object as the message says. |
| `This Seurat object has Seurat v5 assays (<assays>), which SeuratObject <version> cannot read. Install SeuratObject 5 or later: install.packages("SeuratObject").` | Install SeuratObject 5. |
| `This SingleCellExperiment was made by SingleCellExperiment <version>, which Cellucid does not read. Update it first: data <- updateObject(data).` | `updateObject()`. |
| `This <class> holds no embedding (no reductions) for Cellucid to place its cells at. Compute one first, for example <the line>, or pass embeddings = list(X_umap = coordinates).` | Compute a UMAP, or pass `embeddings`. |
| `data holds no expression matrix Cellucid reads.` | The object has no assay with `data` or `counts` (Seurat), or no named assay (SingleCellExperiment). |
| `dataset_id and dataset_name are required: ... An object can declare its own in misc$cellucid (Seurat) or metadata(x)$cellucid (SingleCellExperiment): list(id = , name = ).` | Pass them, or set the object's record. |
| `<path> must be a list such as list(id = "my-study", name = "My study").` | Fix `misc$cellucid` or `metadata(x)$cellucid`. |
| `modalities chooses among the modalities of data, and no data was supplied.` | Pass `data`, or leave out `modalities`. |
| `modalities must name at least one modality of data, each once.` | Fix the list. |
| `modalities names c(...), which data does not hold. Its modalities: c(...).` | Use the names it lists. |
| `c("gene_expression", ...) describe the genes of one feature space, and data holds several (c(...)). Change the assay in the object first, or write one modality with modalities = "<name>".` | Write one modality, or change the object's assay. |
| `modality must name one of the modalities: c(...).` | Use one of those names. |
| `modality names a modality, and the data has one feature space, which is written as an AnnData with no modalities.` | Leave out `modality`. |
| `The obsm key '<key>' is also a modality's name, under which MuData keeps the modality's mask of measured cells. Rename the embedding.` | Rename the embedding. |

## `out_path`, `format`, `force`

| Message | Fix |
|---|---|
| `format must be "zarr" or "zarr.zip":` followed by what each writes, `(.h5ad is written by Python's cellucid.prepare(format = "h5ad"); this package writes no HDF5.)` and `Got <value>.` | Pass `"zarr"` or `"zarr.zip"`; an unambiguous prefix is matched, as `match.arg()` matches. |
| `format = "zarr" writes a .zarr folder, and out_path '<path>' does not end in ".zarr": write '<path>.zarr' instead.` | Write the path it names. |
| `format = "zarr" writes a .zarr folder, and out_path '<name>.zarr.zip' ends in ".zarr.zip", which format = "zarr.zip" writes: write '<name>.zarr' instead.` | Write the path it names, or pass `format = "zarr.zip"`. |
| `format = "zarr.zip" writes one .zarr.zip file, and out_path '<name>.zarr' ends in ".zarr", which format = "zarr" writes: write '<name>.zarr.zip' instead.` | Write the path it names, or pass `format = "zarr"`. |
| `format = "zarr" writes a .zarr folder, and out_path '<name>.h5ad' ends in ".h5ad", which this package does not write (Python's cellucid.prepare(format = "h5ad") does): write '<name>.zarr' instead.` | Write a Zarr store; convert it with Python's anndata if you need an `.h5ad`. |
| `out_path must be exactly one string.` | Pass one path. |
| `zarr_format must be 2 or 3, got <value>.` | Pass `2` (the default), `3`, or leave it `NULL`. |
| `<path> already exists; pass force = TRUE to replace it.` | Choose another path, or pass `force = TRUE`. |
| `<path> exists and is not a Zarr store, so force = TRUE does not replace it; remove it yourself if it should go.` | The folder has neither a `.zgroup` nor a `zarr.json`; `force` never deletes it. |
| `<path> is a folder, and format = "zarr.zip" writes a file; remove it yourself if it should go.` | Remove the folder, or choose another path. |
| `force must be exactly TRUE or FALSE.` | Pass `TRUE` or `FALSE`. |

## Dataset metadata

| Message | Fix |
|---|---|
| `dataset_id '<id>' is not a portable identifier. Use 1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a letter or number and not ending with '.'.` | For example `pbmc-3k`. |
| `dataset_id '<id>' is reserved on Windows.` | Avoid `CON`, `PRN`, `AUX`, `NUL`, `COM1`-`COM9`, `LPT1`-`LPT9`. |
| `dataset_name must be a non-empty string.` | Give the dataset a name. |
| `<argument> must be exactly one string.` | One string, not a vector or `NULL` (for required text). |
| `<argument> is displayed verbatim, so it must not carry characters that have no glyph: '<text>' ends with U+0020 SPACE. Cellucid does not remove them for you, because that would change text you did not ask to change. Pass the exact text you want shown.` | Remove leading/trailing whitespace, control characters, and U+200B, U+2060, U+FEFF from `dataset_name`, `dataset_description`, `citation_text`, `source_name`, or the URLs. |
| `citation_text is required whenever citation_url is supplied.` | Add the citation text. |
| `source_name is required whenever source_url is supplied.` | Add the source name. |
| `<argument> must be an absolute http:// or https:// address; got '<url>'.` | For example `https://doi.org/...`. |

## `embeddings`

| Message | Fix |
|---|---|
| `embeddings must be a named list mapping each obsm key to a matrix, for example list(X_umap = umap).` | Wrap the matrix in a named `list()`. |
| `embeddings must contain at least one embedding.` | The viewer needs something to plot. |
| `embeddings must be a named list: every embedding needs an obsm key.` | Name every element. |
| `Embedding '<key>' must be a real numeric matrix or vector.` | Pass numbers. |
| `Embedding '<key>' has no columns.` | Pass at least one column. |
| `All embeddings must have one row per cell: '<a>' has <n> rows and '<b>' has <m>.` | Align the cells. |
| `Embedding '<key>': every cell's values are missing (NaN or infinite). A cell whose row holds NaN is a cell the array has no value for; at least one cell needs finite values.` | Pass the embedding computed over the cells. A row of `NA`/`NaN` alone is not refused: that cell is written as given and not drawn in that embedding. |
| `Embedding '<key>' holds values beyond the float32 range (magnitude above 3.4e38), which would be stored as infinities.` | Rescale. |

## Names that become paths

Every `obsm` key and every `obs` / `var` column name is a folder in the store.
The same messages appear with `Embedding`, `obs column`, or `var column` as the
subject:

| Message | Fix |
|---|---|
| `Embedding '<key>' contains a character that cannot appear in a file name on every platform (/ \ < > : " \| ? * or a control character). It becomes a path in the Zarr store, so rename it, for example with names(x)[names(x) == '<key>'] <- 'new_name'.` | Rename. |
| `Embedding '<key>' ends with a space or a dot, which Windows removes from file names. It becomes a path in the Zarr store, ...` | Rename. |
| `Embedding '<key>' is reserved by the Zarr store layout. ...` | `.`, `..`, `.zarray`, `.zattrs`, `.zgroup`, `.zmetadata` cannot be names. |
| `Embedding '<key>' is a reserved device name on Windows. ...` | Rename. |
| `Embedding names '<a>' and '<b>' differ only in letter case. They become folders of one Zarr group, and on macOS and Windows those are the same folder, so rename one of them.` | Rename one. |
| `obs column name '<name>' is used twice.` | Make the column names unique. |
| `obs column '_index' is reserved: AnnData stores the row names under that name. Rename the column.` | Rename it. |
| `<what> names must be non-empty strings.` | Name every column. |

## `obs` and `var` columns

| Message | Fix |
|---|---|
| `obs must be a data.frame with one row per cell.` | Pass a `data.frame` (for example `as.data.frame(colData(sce))`). |
| `obs has <n> rows, but the embeddings have <m> cells.` | Align the cells. |
| `The row names of <a> and <b> disagree at row <i> ('<x>' and '<y>'). Every input must list the same cells in the same order.` | Reorder so `obs`, `gene_expression` and the embeddings list the same cells in the same order. |
| `Cell names must not be missing (<source>).` | Fix the row names. |
| `obs column '<name>' is a Date column, which has no AnnData encoding Cellucid writes. Convert it first, for example with as.character() or as.numeric().` | Convert dates, times, and other classed columns. |
| `obs column '<name>' is of type complex, which has no AnnData encoding Cellucid writes. ...` | Convert it. |
| `obs column '<name>' holds 1 infinite value (first at row 2). A colour scale has no position for an infinity; replace it with NA, which is written as a missing value.` | Replace `Inf` with `NA`. |

## `gene_expression`, `var`, `gene_identifiers`

| Message | Fix |
|---|---|
| `var describes the genes of gene_expression, and no gene_expression was supplied.` (also for `gene_identifiers`, `var_gene_id_column`) | Pass `gene_expression`, or drop the argument. |
| `gene_expression must be a real numeric matrix, a Matrix::Matrix, or a matrix held outside R's memory (BPCells, DelayedArray).` | Pass a numeric matrix, a `Matrix` sparse matrix, or a BPCells / DelayedArray matrix. |
| `gene_expression must hold numbers, not a pattern of non-zeros.` | Convert a pattern matrix (`ngCMatrix`) to numbers. |
| `gene_expression has <n> rows, but the embeddings have <m> cells. gene_expression is cells x genes; transpose a genes x cells matrix with t().` | Usually `Matrix::t(expr)`. |
| `gene_expression must contain only finite values. gene_expression cannot be stored: of <n> stored non-zero values, 1 +Inf. ...` | Remove `NA`, `NaN`, `Inf`. |
| `gene_expression holds values beyond the float32 range (magnitude above 3.4e38), which would be stored as infinities.` | Rescale. |
| `Gene <i> holds more values than were counted: the matrix changed while it was read.` / `... fewer values ...` | A matrix held outside R's memory (BPCells, DelayedArray) changed on disk between its passes; keep it unchanged while `cellucid_prepare()` or `cellucid_view()` reads it. |
| `var must be a data.frame with one row per gene.` | Pass a `data.frame`. |
| `var has <n> rows, but gene_expression has <m> genes (columns).` | One `var` row per column of `gene_expression`. |
| `rownames(var) and colnames(gene_expression) disagree at gene <i> ('<x>' and '<y>'). var must describe the columns of gene_expression in their order.` | Reorder `var`. |
| `The genes have no names: set rownames(var) or colnames(gene_expression) to the gene identifiers, or pass var_gene_id_column.` | Name the genes. |
| `Gene names must not be missing or empty.` | Fix the names. |
| `var_gene_id_column '<name>' not found in var. Available columns: c(...).` | Use an existing column. |
| `var column '<name>': floating-point values cannot name genes.` / `var column '<name>': true/false values cannot name genes.` | Use a text, factor, or integer column. |
| `gene_identifiers must not contain duplicate identifiers.` | De-duplicate. |
| `gene_identifiers contains identifiers not found in var: c("<g>"). Identifiers are matched against the gene names the viewer shows (var_gene_id_column, or the var row names), after repeated names were made unique.` | Check the spelling, and which column names the genes. |

## `connectivities`

| Message | Fix |
|---|---|
| `connectivities must be a real numeric or logical matrix or a Matrix::Matrix sparse matrix.` | Pass a matrix. |
| `connectivities must be <n> x <n> (one row and one column per cell), got <a> x <b>.` | Pass the cells x cells graph. |
| `Connectivity weights must all be finite.` / `Connectivity weights must all be non-negative.` | Fix the weights. |
| `Connectivity matrix must be exactly symmetric: a directed neighbour graph ...` | A directed graph (Seurat's `RNA_nn`): pass a symmetric one (`RNA_snn`), or symmetrize, for example `(graph + Matrix::t(graph)) / 2` (`pmax()` makes a large sparse matrix dense). A diagonal is left out, never refused. |
| `The Matrix package is required to store connectivities.` | `install.packages("Matrix")`. |

## `vector_fields`

| Message | Fix |
|---|---|
| `vector_fields must map each field to a named list of embeddings and arrays, for example list(velocity = list(X_umap = velocity_umap)).` | Use the nested list. |
| `vector_fields must name every field.` / `vector_fields names one field more than once.` | Name each field once. |
| `Vector field '<field>' must map embeddings to arrays, for example list(<field> = list(X_umap = vectors)).` | Bind each array to an embedding. |
| `Vector field '<field>' on embedding '<key>': '<key>' is not one of the embeddings: c(...).` | Bind to an embedding you pass. |
| `Vector field '<field>' on embedding '<key>': the viewer binds a vector field to the embedding 'X_<basis>' it is named for ('velocity_umap' to 'X_umap'), so the embedding needs an 'X_<basis>' key. Rename '<key>', for example to 'X_<key>'.` | Name the embedding `X_<basis>`; see {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`. |
| `Vector field '<field>' on embedding '<key>' has <n> columns, and '<key>' has <m>: a vector field has one component per embedding column.` | Match the embedding's column count. |
| `Vector field '<field>' on embedding '<key>' must have one row per cell (<n>), got <m>.` | Align the cells. |
| `Vector field '<field>' on embedding '<key>' and vector field '<other>' on '<other key>' would both be stored as obsm[['<field>_<basis>']]. Rename one field.` | Rename one field. |
| `Vector field '<field>' on embedding '<key>' would be stored as obsm[['<field>_<basis>']], which the viewer reads as <the field '<f>' on 'X_<b>', or an embedding>. Rename the field.` | Another embedding's longer basis claims the key, or the field starts with `X_`; rename the field. |
| `A vector field would be stored as obsm[['<key>']], which is also the key of an embedding. Rename the field or the embedding.` | Rename. |

## `cellucid_transition_drift()`

| Message | Fix |
|---|---|
| `transition_matrix is required: ...` / `embedding is required: ...` / `normalize_rows is required: TRUE divides every row of transition_matrix by its sum, FALSE uses it as given.` | Pass all three; `normalize_rows` has no default. |
| `normalize_rows must be exactly TRUE or FALSE.` | Pass a flag. |
| `transition_matrix must be a real numeric matrix or a Matrix::Matrix matrix.` / `transition_matrix must contain real numeric values.` | Pass numbers, not a logical or pattern matrix or a data frame. |
| `transition_matrix must be <n> x <n> (one row and one column per row of embedding), got <a> x <b>.` | Use the cells x cells matrix of the embedding's cells, in the same order. |
| `transition_matrix must contain only finite values. ...` / `transition_matrix values must be non-negative.` | Fix the matrix. |
| `Every transition_matrix row must have a nonzero finite row sum when normalize_rows = TRUE.` | A cell has no successor: give it one (a self-transition), or pass a matrix that is already row-stochastic with `normalize_rows = FALSE`. |
| `embedding must be a real numeric matrix or vector.` / `embedding must contain at least one row.` / `embedding must contain only finite values. ...` / `embedding holds values beyond the float32 range ...` | As for `embeddings` above. |
| `<embedding or The computed transition drift> holds <n> nonzero value(s) below the smallest float32 (magnitude at most 2^-150), which would be stored as zero.` | Rescale the embedding. |
| `The computed transition drift holds values beyond the float32 range ...` | Rescale the embedding or the matrix. |

## Launch defaults

| Message | Fix |
|---|---|
| `embedding_defaults must map "1d", "2d", or "3d" to an obsm key, for example c("2d" = "X_umap").` | Use a named vector or list. |
| `embedding_defaults key '<k>' must be "1d", "2d", or "3d".` | Use one of those names. |
| `embedding_defaults names one dimension more than once.` | One value per dimension. |
| `embedding_defaults[["2d"]] must name one of the embeddings: c(...).` | Use an `embeddings` key. |
| `embedding_defaults[["3d"]] names '<key>', which has 2 columns; a 3d view needs 3.` | Choose a wider embedding. |
| `embedding_defaults[["2d"]] must be an obsm key or list(key = , columns = ).` | Use exactly those two names. |
| `embedding_defaults[["2d"]]$columns must name 2 different 1-based columns of '<key>', which has columns 1 to <n>.` | Fix the columns. |
| `latent_key must name one of the embeddings: c(...).` | Pass the latent space as an embedding too. |
| `obs_keys contain columns not in obs: c(...). Available columns: c(...).` | Use existing columns. |
| `obs_keys must name at least one obs column.` / `obs_keys names '<name>' twice.` | Fix the list. |
| `vector_field_default names a field, and no vector_fields were supplied.` | Pass `vector_fields`, or drop `vector_field_default`. |
| `vector_field_default must name one of vector_fields: c(...).` | Use a field name as written in `vector_fields` (`"velocity"`). |
| `normalize_embeddings must be exactly TRUE or FALSE.` | Pass a flag. |
| `centroid_min_points must be one positive integer.` | For example `10`. |

## Messages that are not errors

When names repeat, they are made unique the way anndata's
`obs_names_make_unique()` / `var_names_make_unique()` do it, and the call says
so:

```text
Gene names: 1 duplicate name was made unique, e.g. 'A' → 'A-1'.
Cell names: 1 duplicate name was made unique, e.g. 'c' → 'c-1'.
```

The `var_gene_id_column` column is stored unchanged; the viewer fills its empty
values from the var index and makes repeated values unique when it reads it,
and the call says ahead of time what that changes:

```text
Gene names from var column 'symbol': 1 gene has no 'symbol' value and uses its var index name. 1 duplicate name was made unique, e.g. 'CD3E' → 'CD3E-1'.
```

Silence them with `suppressMessages()` once you have checked them.

## Debugging an error

1. Read the subject of the message: it names the argument or column.
2. Check shapes: `dim()` of every matrix, `nrow(obs)`, `nrow(var)`.
3. Check orientation: `gene_expression` is cells x genes.
4. Check names: `rownames()` must agree across inputs.

More symptoms: {doc}`../c_data_preparation_api/11_troubleshooting_prepare` and
{doc}`../i_troubleshooting_index/index`.
