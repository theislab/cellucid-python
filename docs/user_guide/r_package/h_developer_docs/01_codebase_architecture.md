# Codebase Architecture

**Audience:** contributors  
**Time:** 10 minutes  
**What you'll get:** what each source file does, and how one call flows through them

`cellucid-r` is a small package: two exported functions, a hand-written Zarr
(format 2 and 3) and ZIP writer in base R, and one C routine. Its only hard dependency is
`jsonlite`; `Matrix` is suggested and used only when you pass a `Matrix` object
(or a base-matrix neighbour graph).

## Files

| File | Responsibility |
|---|---|
| `R/cellucid_prepare.R` | The exported `cellucid_prepare()`: reports missing required arguments, validates every input, resolves genes (`.resolve_genes()`) and launch defaults (`.resolve_defaults()`), then writes the store. |
| `R/inputs.R` | Embeddings (`.resolve_embeddings()`), `embedding_defaults`, vector fields and their `obsm['<field>_<basis>']` keys (`.resolve_vector_fields()`), the neighbour graph (`.resolve_connectivities()`), and the cell names every input must agree on. |
| `R/matrix-input.R` | What counts as a matrix (base matrix, numeric data frame, `Matrix`); the expression matrix as a `dgCMatrix` or a dense matrix; the finite / float32-range check. |
| `R/var.R` | Which `var` column names the genes (`var_gene_id_column`), and the names it gives them. |
| `R/make-unique.R` | anndata's `make_index_unique()` rule for repeated names, and the one-line note that reports what changed (pinned by the shared conformance fixture). |
| `R/element-keys.R` | Names that become paths in the store: `obsm` keys and `obs` / `var` column names must be portable file names and differ by more than letter case. Also the vector-field key rule (`.parse_vector_field_key()`), shared with the web app and cellucid-python. |
| `R/dataset-identity.R`, `R/identifiers.R`, `R/validate-display-text.R` | `uns['cellucid']` metadata: the portable `dataset_id`, text shown verbatim, and `http(s)` links. |
| `R/validate-arguments.R` | Flags, strings, numbers, and the missing-required-arguments message. |
| `R/error-messages.R` | How values are quoted back (`c("a", "b")`), and the count of non-finite values in a refused matrix. |
| `R/anndata-elements.R` | AnnData's encodings (`encoding-type` / `encoding-version`): dataframe 0.2.0, categorical 0.2.0, string-array 0.2.0, nullable-boolean / nullable-integer 0.1.0, csc/csr_matrix 0.1.0, dict 0.1.0, string and numeric-scalar 0.2.0; the column-encoding rules; dense chunk shapes. |
| `R/zarr-store.R` | The Zarr writer, format 2 (`.zgroup` / `.zattrs` / `.zarray` JSON, `.zmetadata`) and 3 (`zarr.json`, `c/…` chunk keys, consolidated metadata in the root `zarr.json`): little-endian dtypes, `vlen-utf8` strings, the byte shuffle (every array wider than one byte except the expression values), deflate via `memCompress()` (zlib framing for format 2, gzip framing for 3), edge chunks padded to full size, chunk sizes (2^18 values for arrays read whole, 2^16 for the data and indices of `X`); in a `.zarr.zip` the metadata documents are written after every chunk, the consolidated metadata last. |
| `R/zip.R` | The `.zarr.zip` container: stored (method 0) entries with UTF-8 names, fixed 1980-01-01 timestamps, ZIP64 records only where needed. |
| `R/output.R` | `out_path` rules (`.zarr` / `.zarr.zip`, `force`), and writing the store to a temporary sibling that is moved into place once complete. |
| `R/outlier-quantiles.R`, `src/latent_outliers.c` | With `latent_key`, `uns['cellucid']['outlierQuantiles']`: which `obs` columns are categorical and how they partition the cells (`.outlier_partition()`), and the per-category quantile rule in C, bit for bit the one of cellucid-python's `_outliers.py` and the web app's `data/latent-outliers.js` (pinned by the web app's numpy conformance fixture). |
| `R/native.R`, `src/crc32.c`, `src/init.c` | The two `.Call()`s, registered in `src/init.c`: the outlier quantiles of one column, and the CRC-32 of a raw vector, which each ZIP entry and each gzip-framed chunk needs (base R has none, and an R loop is too slow for expression chunks). |
| `R/transition-drift.R` | The exported `cellucid_transition_drift()`: `T %*% E - E` from a cell-cell transition matrix, operation for operation cellucid-python's `compute_transition_drift()` (pinned by the conformance fixture the two packages share). |
| `R/zzz.R` | Unloads the shared library. |
| `man/cellucid_prepare.Rd`, `man/cellucid_transition_drift.Rd`, `NAMESPACE` | Hand-written and authoritative; no roxygen. `NAMESPACE` holds `useDynLib(cellucid, .registration = TRUE, .fixes = "C_")`. |

## One call, in order

1. **Required arguments.** Missing `out_path`, `format`, `embeddings`, `dataset_id`,
   `dataset_name` are reported together.
2. **Scalars.** `format`, `force`, `out_path` (suffix, existing target),
   `dataset_id`, the displayed text, citation and source.
3. **obsm.** Embeddings (shape, finite, float32 range, names), then vector
   fields bound to `X_<basis>` embeddings of the same width, each key checked
   to read back as its field and embedding (`.parse_vector_field_key()`).
4. **obs.** A `data.frame` with one row per cell; column names checked.
5. **X and var.** Orientation, gene names (row names of `var`, else column
   names of `gene_expression`), made unique, `gene_identifiers` subset,
   `var_gene_id_column`.
6. **Cell names.** Row names from `obs`, `gene_expression` and the embeddings
   must agree; made unique; else `"0"` to `"n - 1"`.
7. **obsp.** The neighbour graph: square, symmetric, zero diagonal,
   non-negative.
8. **uns['cellucid'].** `schema`, `id`, `name`, optional `description`,
   `citation`, `source`, only the defaults that were set, and with
   `latent_key` the outlier quantiles of every categorical `obs` column.
9. **Write.** Root group, `X`, `obs`, `var`, `obsm`, `obsp`, empty `varm`,
   `varp`, `layers`, then `uns` — into `.<name>.<pid>.partial` beside the
   target; `.zmetadata` last. Then the existing target (with `force`) is moved
   aside, the new store renamed into place, and the old one deleted. On any
   error the partial store is removed.

Column-type problems (a `Date` column, an infinite value) are found while
writing, so they exercise the cleanup in step 9.

## Design rules

- The store must read without a warning in Python's anndata (0.8+ encodings)
  and open in the web app; `tests/testthat/test-anndata-interop.R` checks the
  first.
- No other dependency: JSON with `jsonlite`, compression with base R.
- A change to what the store holds is coordinated with `cellucid-python` and the
  web app; see {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.
