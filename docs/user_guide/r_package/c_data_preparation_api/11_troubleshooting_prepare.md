# Troubleshooting `cellucid_prepare()`

**Audience:** everyone whose call stopped with an error
**Time:** look up your message
**Goal:** what each error means and how to fix it.

Every check runs before the store is moved into place, so after an error there
is nothing at `out_path` (and with `force = TRUE`, the old store is still
there). Fix the input and run the call again. The messages below are copied
from the package; paths and values in them are examples.

## Arguments

**`cellucid_prepare() is missing 1 required argument:`** followed by the
argument and what it takes. `out_path`, `format`, `embeddings`, `dataset_id`,
and `dataset_name` have no default; pass each one.

**`format must be "zarr" or "zarr.zip": ...`**
Choose `"zarr"` (a `.zarr` folder of small files, for GitHub) or `"zarr.zip"`
(one `.zarr.zip` file, for everything else: single-file hosts, local files,
`cellucid serve` and Jupyter; {doc}`index`). `"h5ad"` is refused: this package writes no HDF5; Python's
`cellucid.prepare(format = "h5ad")` writes `.h5ad`.

**`format = "zarr" writes a .zarr folder, and out_path 'x.zarr.zip' ends in ".zarr.zip", which format = "zarr.zip" writes: write 'x.zarr' instead.`**
`out_path` must end in the suffix of `format`: write the path the message
names, or choose the format that matches your path.

**`... already exists; pass force = TRUE to replace it.`**
Choose another path, or pass `force = TRUE` to replace the store.

**`... exists and is not a Zarr store, so force = TRUE does not replace it; remove it yourself if it should go.`**
`out_path` is a folder that does not hold a Zarr store (no `.zgroup`).
Cellucid never deletes such a folder; pick another path or remove it yourself.

**`dataset_id 'a b' is not a portable identifier. Use 1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a letter or number and not ending with '.'.`**
Use an id such as `pbmc-3k` or `liver_atlas.v2`.

**`dataset_name is displayed verbatim, so it must not carry characters that have no glyph: 'A ' ends with U+0020 SPACE. ...`**
Remove the leading or trailing whitespace, control character, or zero-width
character it names. The same rule applies to `dataset_description`,
`citation_text`, `citation_url`, `source_name`, and `source_url`.

**`citation_text is required whenever citation_url is supplied.`** (and
**`source_name is required whenever source_url is supplied.`**)
A link is stored only beside the text it belongs to.

**`source_url must be an absolute http:// or https:// address; got 'x.org'.`**
Write the full address, `https://x.org`.

**`centroid_min_points must be one positive integer.`**
For example `10`, the number used when it is not set.

**`force must be exactly TRUE or FALSE.`** (and the same for
`normalize_embeddings`)

## Embeddings

**`embeddings must contain at least one embedding.`** /
**`embeddings must be a named list: every embedding needs an obsm key.`**
Pass `embeddings = list(X_umap = umap)`.

**`All embeddings must have one row per cell: 'X_umap' has 3 rows and 'X_b' has 2.`**
Every embedding must describe the same cells in the same order.

**`Embedding 'X_umap' must contain only finite values. ... of 6 values, 1 NA. First affected position: [2, 1]. ...`**
Remove or impute the cells with missing coordinates in every input.

**`Embedding 'X/umap' contains a character that cannot appear in a file name on every platform ...`**
Rename it; the message shows how. See the path rules in
{doc}`02_input_requirements_global`.

**`Embedding names 'X_umap' and 'x_UMAP' differ only in letter case. ...`**
On macOS and Windows those are one folder; rename one.

**`embedding_defaults[["3d"]] names 'X_umap', which has 2 columns; a 3d view needs 3.`**
Pick an embedding with enough columns, or a lower dimension.

**`embedding_defaults[["3d"]]$columns must name 3 different 1-based columns of 'X_pca', which has columns 1 to 50.`**

**`latent_key must name one of the embeddings: c("X_umap").`**
Pass the latent space (for example the PCA) in `embeddings` too.

## obs

**`obs has 2 rows, but the embeddings have 3 cells.`**

**`The row names of gene_expression and obs disagree at row 2 ('q' and 'y'). Every input must list the same cells in the same order.`**
Reorder one input to match the other (for example
`obs <- obs[rownames(expr), , drop = FALSE]`).

**`obs column 'day' is a Date column, which has no AnnData encoding Cellucid writes. Convert it first, for example with as.character() or as.numeric().`**

**`obs column 'score' holds 1 infinite value (first at row 2). A colour scale has no position for an infinity; replace it with NA, which is written as a missing value.`**

**`obs column '_index' is reserved: AnnData stores the row names under that name. Rename the column.`**

**`obs_keys contain columns not in obs: c("b"). Available columns: c("a").`**

## Gene expression and var

**`gene_expression has 2000 rows, but the embeddings have 5000 cells. gene_expression is cells x genes; transpose a genes x cells matrix with t().`**
Use `Matrix::t(expr)`.

**`gene_expression must contain only finite values. ...`**
Replace `NA`, `NaN`, and `Inf` values first.

**`var describes the genes of gene_expression, and no gene_expression was supplied.`**
Pass `gene_expression`, or leave out `var`, `gene_identifiers`, and
`var_gene_id_column`.

**`var has 3 rows, but gene_expression has 4 genes (columns).`**

**`The genes have no names: set rownames(var) or colnames(gene_expression) to the gene identifiers, or pass var_gene_id_column.`**

**`rownames(var) and colnames(gene_expression) disagree at gene 1 ('g1' and 'ENSG1'). var must describe the columns of gene_expression in their order.`**

**`var_gene_id_column 'nope' not found in var. Available columns: c("symbol").`**

**`var column 'score': floating-point values cannot name genes.`**
Use a text, factor, or integer column.

**`gene_identifiers contains identifiers not found in var: c("NOPE"). Identifiers are matched against the gene names the viewer shows (var_gene_id_column, or the var row names), after repeated names were made unique.`**
A repeated gene is addressed by its unique name (`"CD3E-1"`); see
{doc}`05_var_gene_metadata`.

## Connectivities

**`connectivities must be 3 x 3 (one row and one column per cell), got 2 x 2.`**

**`Connectivity matrix must be exactly symmetric.`**
Symmetrize it, for example `(knn + Matrix::t(knn)) / 2`.

**`Connectivity diagonal values must all be exactly zero.`**
Remove self-loops: `Matrix::diag(graph) <- 0`.

**`Connectivity weights must all be non-negative.`** /
**`Connectivity weights must all be finite.`**

## Vector fields

**`Vector field 'velocity' on embedding 'pca': the viewer binds a vector field to the embedding 'X_<basis>' it is named for ('velocity_umap' to 'X_umap'), so the embedding needs an 'X_<basis>' key. ...`**
Name the embedding `X_<basis>` (`X_pca`); see
{doc}`08_vector_fields_velocity_displacement`.

**`Vector field 'velocity' on embedding 'X_umap' has 3 columns, and 'X_umap' has 2: a vector field has one component per embedding column.`**
Pass an array exactly as wide as its embedding.

**`Vector field 'T_fwd' on embedding 'X_umap' would be stored as obsm[['T_fwd_umap']], which the viewer reads as the field 'T' on 'X_fwd_umap'. Rename the field.`**
Another embedding's longer basis (here `X_fwd_umap`) claims the stored name;
rename the field.

**`Vector field 'velocity' on embedding 'X_tsne': 'X_tsne' is not one of the embeddings: c("X_umap").`**

**`vector_field_default must name one of vector_fields: c("velocity").`**

## Messages that are not errors

**`Gene names: 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.`**
and **`Cell names: ...`**: repeated names were made unique the way anndata does
it. The store was written.

## Still stuck

- {doc}`../i_troubleshooting_index/02_data_preparation_issues`
- {doc}`../g_api_reference_coverage/02_error_messages_and_exceptions_document_patterns`
- The same topic for Python, `cellucid.prepare()`: {doc}`../../python_package/c_data_preparation_api/11_troubleshooting_prepare`
