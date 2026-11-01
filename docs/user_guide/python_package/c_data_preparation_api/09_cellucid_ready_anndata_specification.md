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
| `obs` | One row per cell: in {ref}`priority order <priority-order>`, `obs['cellucid_original_position']` (`<u4`) holding each row's original position, or in the order given with `priority_seed=None`. The index (`_index`) holds unique cell names. Text columns whose values repeat or are missing are `categorical` (categories in code-point order, missing = code -1); other text columns are `string-array`; categoricals stay categoricals; numbers stay numbers in their own dtype (no infinities; NaN is missing); booleans stay booleans; integer and boolean columns with missing values are `nullable-integer` / `nullable-boolean`. |
| `var` | One row per column of `X`. The index holds unique gene names (usually gene ids); a gene-symbol column, when present, is named by `uns['cellucid']['defaults']['geneIdColumn']`. |
| `obsm/<key>` | Embeddings: dense `float32`, two-dimensional, one row per cell, every column of the input at its own coordinates (the viewer selects columns and normalizes). The latent space for outlier measurement is one of them. A row holding `NaN` (or an infinity) is a cell the embedding does not place (one computed on some of the cells, a MuData modality that did not measure the cell): it is written as given; the viewer draws it in no view of that embedding and says how many such cells there are, naming the first; other embeddings, filters and analyses keep it. Every other value is finite and within the `float32` range, and at least one cell is placed. |
| `obsm/<field>_<basis>` | A vector field, dense `float32`, one row per cell and exactly as many columns as the embedding `obsm['X_<basis>']` it is named for (scVelo's naming: `velocity_umap` over `X_umap`, `velocity_branches` over `X_branches`). The viewer reads it as the field `<field>` drawn over that embedding: the vector-field key rule below. A `NaN` row is a cell without a vector: drawn without an arrow. |
| `obsp/connectivities` | The neighbour graph, when present: `csr_matrix`, `float32`, cells × cells, exactly symmetric, finite non-negative weights, zero diagonal, no stored zeros. The writers leave a given graph's diagonal out (a self-loop is no edge: Seurat's `RNA_snn` has 1 there) and refuse a directed one (`RNA_nn`, `distances`); from an object without `connectivities` they write the graph the viewer's default rule opens with: Seurat's `<assay>_snn` of the default assay, else the first other `*_snn`. A reader ignores a diagonal the same way. |
| `varm` | Present; empty, except `varm['cellucid_prefix_counts']` beside a sparse `X` of a file in priority order ({ref}`priority order <priority-order>`). |
| `varp`, `layers` | Present and empty. |
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

0. **One namespace at a time.** A MuData modality's arrays are named
   `<modality>:<key>` (the viewer's MuData view, and §7), and each modality's
   names are its own: the rules below apply to the part after the first `:`,
   and the embedding a key names is in the same namespace. `rna:velocity_umap`
   is the field `rna:velocity` over `rna:X_umap`, never over the dataset's
   `X_umap`; `rna:X_umap` is an embedding. A key without `:` is the dataset's
   own.
1. **`X_` names an embedding.** An `obsm` key (its part after the namespace)
   that starts with `X_` is never a vector field.
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
| `rna:velocity_umap`, `rna:X_umap`, `X_umap` | `rna:velocity` | `rna:X_umap` |

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
| `priority` | mapping `{seed, specVersion, nObs, levels}` | The file stores its cells in priority order ({ref}`priority order <priority-order>`): `seed` (integer, 0 to 2^32 − 1) the seed of the ranks, `specVersion` (integer `1`) the version of the {doc}`subsampling specification <../../web_app/p_developer_docs/19_subsampling_specification>` the ranks follow, `nObs` (integer) the cell count the file was written with, `levels` (integer list, ascending, each at least 1: those of `nObs`) the prefix sizes `varm['cellucid_prefix_counts']` counts at. Unlike every other key, a `priority` that does not describe the file's cells is never refused: the file is read as one Cellucid did not write ({ref}`priority order <priority-order>`, "Verified before use"). |

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
| `geneIdColumn` | text (a `var` column) | The column whose values name the genes. Absent: the `var` index. | `var_gene_id_column` |
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

A cell whose latent row is not finite (`NaN` in any column: a cell the latent
space does not place) is left out as a cell with a missing value is: it adds
nothing to its category's centre or count. Only categories with at least
`minPoints` such placed cells are measured; every other cell (a missing value,
an unplaced cell, a smaller category) is `NaN`.

| Key | Type | Meaning |
|---|---|---|
| `latentKey` | text (an `obsm` key) | The latent space they were computed in: `defaults.latentKey`. |
| `minPoints` | integer (`<i8`) | The fewest cells a category needed: `centroid_min_points`, else 10. |
| `fields` | mapping | `{<obs column>: <f4 array of n_obs>}`, one array per field, in obs column order. |

Each field is an ordinary `array` element (`<f4`, chunks of 2^16 values, as obs columns,
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
category, whose `n` cells with a value are its placed cells (at least
`minPoints`), each value is `r / n` for a whole `r` from 1 to `n`, exactly `r`
cells have a rank of at most `r`, and every other cell is `NaN`.
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

(priority-order)=
### Priority order

Every cell has a **priority**: its `uniform` rank in the
{doc}`subsampling specification <../../web_app/p_developer_docs/19_subsampling_specification>`,
`rank(seed, p)` of its original position `p` under the file's seed. For one
seed the rank is a bijection of the uint32 positions, so no two cells share
one. A file in priority order stores its cells by ascending rank: row `r`
holds the cell of the `r`-th lowest rank. The `uniform` subsample of `k` cells
under the file's seed is then exactly the first `k` rows of every array with
one row per cell, and the first values of every gene of `X`, so opening it
reads one range of each array and never a byte of a cell it left out, at any
chunk size and compression (at most one chunk is read in part).

Both writers store cells this way by default (`priority_seed=0`;
`priority_seed=None` stores them as given, with no `priority`, no position
column and no prefix counts). The contract, which the web app checks as it
opens the file (`assets/js/data/anndata/priority-layout.js`):

1. **The order.** `uns['cellucid']['priority']['seed']` is `s` and
   `uns['cellucid']['priority']['nObs']` is `n_obs`, the cell count the file
   was written with. `obs['cellucid_original_position']` (`<u4` as written;
   a reader also takes `<i8`) holds each row's original position: a
   permutation of `0 … n_obs − 1` whose ranks `rank(s, position)` strictly
   ascend down the rows. So the first `k` rows hold exactly the cells
   `subsample_cells(n_obs, k, seed=s)` keeps, and the column is a function of
   `s` and `n_obs` alone: the positions `0 … n_obs − 1` sorted by rank. The
   web app computes it from the seed (an `n_obs`-element sort, 26 ms for a
   million cells, instead of a 3 MB read), once per opened store however
   often its cells are opened again, and checks it against the stored
   column's first chunk of rows, which it reads with the record (in an
   `.h5ad`, with the file's first request: the writers place the column with
   what opening reads).
2. **Every per-cell element in that order.** `obs` (every column, categorical
   codes, nullable values and masks, the index), `obsm`, the stored outlier
   quantiles, the rows of a dense `X`, and **both axes** of `obsp` graphs
   (row `r`, column `c` is the edge between the cells of rows `r` and `c`; the
   induced subgraph of a prefix is rows `[0, k)` with columns below `k`).
   Derived values are computed on the cells as stored: the outlier quantiles'
   sums run in file row order.
3. **Gene-major `X`.** Within each gene the row indices ascend (as in every
   Cellucid-ready `X`), so a gene's values in the first `k` rows are the first
   ones of its range `X/indptr[g] … X/indptr[g+1]`.
4. **Prefix counts.** `uns['cellucid']['priority']['levels']` lists the prefix
   sizes, ascending: from `n_obs`, each level is `floor(previous × 3 / 4)`,
   down to 1; `n_obs` itself is not a level, and a level the rounding repeats
   is kept once (2,500,000 cells: 48 levels). Beside a sparse `X`,
   `varm['cellucid_prefix_counts']` is `n_vars × len(levels)` `<u4`, column
   `j` each gene's values in rows below `levels[j]`; every column is
   non-decreasing across levels per gene and at most the gene's value count.
   It is chunked one column at a time (`[min(n_vars, 2^18), 1]`), so a reader
   opening a prefix of `k` cells reads the one column of the smallest level
   `≥ k` once, and each gene's read is the first `counts[g]` values of its
   range (filtered to rows below `k`); with no level `≥ k` the gene is read
   whole. A reader checks each count against the gene's range and the rows it
   read: every counted row below the level, the next value (when there is one)
   at or past it.
5. **MuData.** The root `obs` is in priority order. Each modality's rows are
   ordered by the file rows of the cells they measured (`obsmap` renumbered),
   so a global prefix of `k` cells is a prefix of every modality: its rows
   whose `obsmap` entry lies in the first `k`. A modality's
   `obs['cellucid_original_position']` holds each row's original row in the
   modality; its `varm['cellucid_prefix_counts']` counts its values in its
   rows of cells below each global level.

**Verified before use.** The record is a claim about the rows, and anndata
keeps `uns` as it was when a file is subset, combined or reordered and saved
again. Every reader therefore verifies it before relying on it, by one rule
(web app: `readPriorityLayout`; Python: `stored_order`; R: `.stored_order`):

- `seed` is a whole number from 0 to 2^32 − 1 and `specVersion` the version
  the reader follows;
- `nObs` is the file's cell count (a cell removed or added changes it);
- `levels`, when present, are the levels of `nObs` (point 4);
- the position column, when present, is the seed's order: the web app checks
  its first chunk of rows (each a position of the file, ranks ascending — a
  reordered file fails at its first moved row), Python and R, which hold the
  whole column in memory, check all of it.

A record that verifies is relied on, column or not (a file whose position
column was dropped still opens in its priority order: the positions are
computed from the seed). A record that does not verify — or a position
column with no record beside it — is never refused and never misread: the
file is read as one Cellucid did not write, each row its own original
position, with a foreign file's default sampling (`blocks`); the position
column is then an ordinary field, and the opening says why (a note in
**Cells to open**). Python and R take such an object the same way, with a
warning: {func}`~cellucid.restore_original_order` returns it as it is, a
session applies to its rows as they are (as the web app read the file the
session was saved on), and `prepare()` writes its rows as given without the
stale column. The one case the web app's check cannot see is a file whose
first chunk of rows (65,536 cells) and cell count were kept while later rows
were reordered; restore the original order before reordering a file's cells.

**Identity.** Everything Cellucid hands back — viewer events and commands,
sessions, exports, pulls — names cells by original position (and obs name),
never by file row. The web app maps an open row to its file row (a prefix:
the same row) and that to its original position (point 1). A saved session
holds every per-cell chunk by ascending original position and its cell-order
digest over the cells in that order, so it does not depend on the order a file
stores its cells in ({doc}`../../web_app/p_developer_docs/10_sessions_persistence_and_serialization`). The viewer does
not list the position column as a field.

**Reading it elsewhere.** `anndata.read_h5ad` and `anndata.read_zarr` (and
mudata) give the cells in the stored order, which is as valid an AnnData as
any. {func}`~cellucid.restore_original_order` (R:
`cellucid_restore_original_order()`) returns them in the order given, without
the position column, `priority` and the prefix counts; a backed AnnData
(`read_h5ad(path, backed="r")`) is read into memory in that order:

```python
adata = cellucid.restore_original_order(anndata.read_zarr("pbmc.zarr"))
```

In R the record is where the importers that keep `uns` put it (anndataR,
zellkonverter): `metadata(sce)$cellucid$priority` of a SingleCellExperiment,
`misc$cellucid$priority` of a Seurat object, the `"cellucid"` attribute of a
data frame. A Seurat object is reordered in every part — each assay's layers
(or matrices), reductions, graphs (both axes), neighbour indices, cell
metadata and identities.

Preparing an object in a verified stored order again (`prepare(data=adata)`,
`cellucid_prepare(data = object)`) puts it back in its original order first,
so every cell keeps its original position whatever `priority_seed` is: with a
seed the file is in that seed's order, with `priority_seed=None` in the
original order. A position column with no record that describes the cells is
not written (with a warning): it names no order.

**Measured.** On the Suo atlas (561,947 cells × 5,103 genes) the priority
order costs 6.4% of the stored size (186.0 → 197.9 MB of the arrays that
change): the gene indices of `X` 44.4 → 50.6 MB, the graph's indices
26.0 → 27.8 MB, the embeddings 0–5%, and the obs columns that were sorted by
sample (about 1 KB each) about 0.3 MB each, 1.3 → 3.5 MB in all. The level
factor 3/4 was chosen against halving and 8/7 on Suo and a 1,000,000-cell
synthetic `X` (2^16-value chunks): a gene's prefix at 1, 5, 10 and 30% of the
cells read up to the next level fetched at most 1.01× the bytes of its exact
prefix (halving: up to 1.09×; 8/7: 1.00×), for a table of 0.09–0.27% of `X`
(halving 0.04–0.11%, 8/7 0.19–0.58%); the one column a subsample reads is
1.6–4.1 KB stored. Before this layout a subsample's genes were read whole and
filtered (each gene's whole range, at every subsample size), and a
cell-block offset table was measured and not chosen: compressed chunks of
2^16 values hold several genes of a typical atlas, so a block's range still
fetched its whole chunk (Suo at 1%: 79 of 103 KB a gene), and making it exact
needed either chunks 16 times smaller (21,334 chunk files for Suo's `X`) or no
compression (`X` 59 → 350 MB), and a round trip for each gene's table row.

Read with the web app's reader over a simulated link (`scripts/benchmark-anndata`
steps, HTTP/2; synthetic `.h5ad` files of 2,000 genes written by `prepare()`,
cells grouped by sample as atlases are): the writer before the priority order
(cells as given, chunks of 2^16 and 2^18; a subsample by `blocks`) against
the priority order with its chunks (a `uniform` subsample under the file's
seed). Every step's bytes summed, and the steps' times:

| Cells, link, share | Before: open / 99 genes / all steps | Priority order: open / 99 genes / all steps |
|---|---|---|
| 300,000; 80 ms, 50 MB/s; 1% | 138 ms / 0.36 s / 13.2 MB | 150 ms / 0.22 s / 3.7 MB |
| 300,000; 150 ms, 5 MB/s; 1% | 345 ms / 2.78 s / 15.0 MB | 338 ms / 0.72 s / 3.4 MB |
| 1,000,000; 80 ms, 50 MB/s; 1% | 141 ms / 0.45 s / 16.4 MB | 156 ms / 0.24 s / 3.4 MB |
| 1,000,000; 80 ms, 50 MB/s; 10% | 144 ms / 0.44 s / 20.3 MB | 170 ms / 0.25 s / 4.8 MB |
| 1,000,000; 150 ms, 5 MB/s; 1% | 299 ms / 3.02 s / 16.4 MB | 359 ms / 0.71 s / 3.4 MB |
| 1,000,000; 150 ms, 5 MB/s; 10% | 303 ms / 3.02 s / 20.0 MB | 348 ms / 0.78 s / 4.9 MB |
| 1,000,000; 80 ms, 50 MB/s; every cell | 139 ms / 0.46 s / 27.4 MB | 172 ms / 0.43 s / 26.8 MB |
| 1,000,000; 150 ms, 5 MB/s; every cell | 332 ms / 3.02 s / 23.6 MB | 370 ms / 1.99 s / 18.9 MB |

Opening reads the same 0.52 MB either way: the original positions are
computed from the seed, and only the stored column's first chunk is read, with
the file's first request. What opening costs more (15–60 ms) is the reader's
own work, ranking every cell. The obs columns that were sorted by sample
compress to nearly nothing as given and to 0.05–0.7 MB each in priority order.

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
| `X/data`, `X/indices` (read one gene's range, or its first values, at a time) | 16,384 values (2^14) |
| every array with one row per cell: obs columns and index, categorical codes, nullable values and masks, the stored outlier quantiles, MuData's masks and cell maps | 65,536 values (2^16) |
| `obsm` arrays (read by rows) | `[rows, n_columns]`, rows = 65,536 / n_columns (at least 1) |
| `X/indptr`, var columns and index, categories, every array of `obsp/connectivities`, `uns` lists (read whole) | 262,144 values (2^18) |
| dense `X` (read one gene at a time) | `[rows, genes]`: rows up to 262,144, genes as many as keep the chunk at 262,144 values (at least 1) |

An array shorter than its chunk size is one chunk of exactly its length; an
empty array has chunks of 1.

A store a server presents for data it holds (`serve_anndata`, `show_anndata`,
`cellucid_serve`, `cellucid_view`) is not a file: it holds the same arrays,
uncompressed, chunked alike except an expression matrix's `data` and `indices`,
which it serves in chunks of 262,144 values (2^18): over loopback a request
costs more than a chunk of 2^14 values moves, and a viewer reads an
uncompressed chunk in part (a `Range`), so a gene or a subsample's first
entries of one read no more than with smaller chunks.

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
- **Chunks of `X`**: 2^14 values since the priority order (27 October 2026).
  Measured on a 1,000,000-cell `.h5ad` in priority order over 150 ms and
  5 MB/s: against 2^16, a 1% subsample's 99 genes fetched 2.4 instead of
  8.9 MB (0.68 instead of 1.94 s) and every cell's 99 genes 8.7 instead of
  14.5 MB (1.90 instead of 3.06 s), each gene dragging fewer of its
  neighbours' values; 2^12 fetched less still (0.8 MB) but its chunk index,
  16 times larger, cost the first gene 0.64 MB and a round trip. The cost is
  requests: a whole-matrix read of a `.zarr` folder fetches four times as many
  files as with 2^16 (before, on Suo: 2^15 doubled them over HTTP/1.1, 41 s
  instead of 21 s for the two cell types); a single file (`.h5ad`,
  `.zarr.zip`) reads neighbouring chunks as one range.
- **Chunks of the per-cell arrays**: 2^16 values (an embedding's rows: 2^16
  over its width), so a subsample — the first rows in priority order — reads
  little beyond its cells: on the file above a 1% subsample's 2-D embedding
  fetched 0.21 MB against 0.84 MB with 2^18, its obs fields 0.05 against
  0.18 MB, its opening 0.70 against 1.20 MB; reading every cell cost the same.
  2^14 halved those bytes again (0.10 MB) for four times the files.
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
