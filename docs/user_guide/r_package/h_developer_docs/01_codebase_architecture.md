# Codebase Architecture

**Audience:** contributors  
**Time:** 10 minutes  
**What you'll get:** what each source file does, and how one call flows through them

`cellucid-r` has three parts that share one core: the writer
(`cellucid_prepare()`), the viewer (`cellucid_view()`, `cellucid_serve()` and
the functions that talk to an open viewer), and the session reader
(`cellucid_apply_session()`). The core is the **dataset plan**: every
component checked once and a function that puts the whole store into a
backend, a folder, a ZIP file or memory. One plan for both writing and serving
is what makes the served store the same bytes as the written one.

Everything is base R plus a hand-written Zarr (format 2 and 3) and ZIP writer
and a few C routines. The imports are `jsonlite` (JSON), and `httpuv`, `later`
and `curl` for the viewer; `Matrix`, `SeuratObject`, `SingleCellExperiment`,
`SummarizedExperiment`, `S4Vectors` and `DelayedArray` are suggested, each
loaded only when an object of its class arrives.

## Files

### The writer and the plan

| File | Responsibility |
|---|---|
| `R/cellucid_prepare.R` | The exported `cellucid_prepare()`: reports every missing required argument at once, then checks `format`, `zarr_format` and `force`, resolves `out_path`, builds the dataset components, plans the dataset, and writes it. |
| `R/dataset-from.R` | What a Seurat, SingleCellExperiment or SummarizedExperiment object holds, as the plan's components: `obs`, `obsm['X_<name>']`, one modality per assay or experiment, the neighbour graph, and the object's own record (`misc$cellucid`, `metadata(x)$cellucid`). The mapping is the web app's for the same object saved with `saveRDS()`. Seurat v5 split layers are joined on a copy; matrices held outside R's memory (BPCells, DelayedArray backends) are read block by block. |
| `R/dataset-plan.R` | `.dataset_components()` (an object's components, replaced by the arguments given), `.plan_dataset()` (every check, the one or several feature spaces, the default modality, `uns['cellucid']`), `.resolve_genes()` and `.resolve_defaults()`; the plan's `write()` puts the AnnData, or the MuData of several modalities, into a store backend. |
| `R/expression-input.R` | One expression matrix with its orientation: a caller's cells x genes, or an object's genes x cells, never transposed whole: its gene-major indices and values are made block by block (`.gene_major_arrays()`, `src/chunks.c`), one binary search a cell; a served store keeps the blocks in its cache. |
| `R/mudata.R` | The Cellucid-ready MuData: MuData's own Zarr layout, `mod/<name>` AnnData groups, `obsmap`/`varmap`, the per-modality masks in `obsm`. |
| `R/inputs.R` | Embeddings, `embedding_defaults`, vector fields and their `obsm['<field>_<basis>']` keys (`.resolve_vector_fields()`), the neighbour graph (`.resolve_connectivities()`), and the cell names every input must agree on. |
| `R/matrix-input.R` | What counts as a matrix (base matrix, numeric data frame, `Matrix`); the finite / float32-range check. |
| `R/var.R` | Which `var` column names the genes (`var_gene_id_column`), and the names it gives them. |
| `R/make-unique.R` | anndata's `make_index_unique()` rule for repeated names, and the one-line note that reports what changed (pinned by the shared conformance fixture). |
| `R/element-keys.R` | Names that become paths in the store: `obsm` keys, `obs` / `var` column names and modality names must be portable file names and differ by more than letter case. Also the vector-field key rule (`.parse_vector_field_key()`), shared with the web app and cellucid-python. |
| `R/data-mapping.R` | A session's data mapping (embeddings, graph, gene names, expression matrix, latent space, vector fields, modality), checked by the web app's and cellucid-python's rules; also the default modality's shared rule. |
| `R/dataset-identity.R`, `R/identifiers.R`, `R/validate-display-text.R` | `uns['cellucid']` metadata: the portable `dataset_id`, text shown verbatim, and `http(s)` links. |
| `R/validate-arguments.R` | Flags, strings, numbers, and the missing-required-arguments message. |
| `R/error-messages.R` | How values are quoted back (`c("a", "b")`), and the count of non-finite values in a refused matrix. |
| `R/encoding.R` | Text that is the same in every R session, whatever its native encoding (Windows R before 4.2 runs in a code page): every message, warning and error is raised as a condition object (`.inform()`, `.warn()`, `.abort()`), so no step re-encodes its text; bytes read as UTF-8 text (`.utf8_text()`); code points and code-point order whatever a string's encoding; environment keys made from a string's UTF-8 bytes (`.symbol_keys()`). |
| `R/anndata-elements.R` | AnnData's encodings (`encoding-type` / `encoding-version`): dataframe 0.2.0, categorical 0.2.0, string-array 0.2.0, nullable-boolean / nullable-integer 0.1.0, csc/csr_matrix 0.1.0, dict 0.1.0, string and numeric-scalar 0.2.0; the column-encoding rules; dense chunk shapes. |
| `R/outlier-quantiles.R` | With `latent_key`, `uns['cellucid']['outlierQuantiles']`: which `obs` columns are categorical and how they partition the cells (`.outlier_partition()`); the per-category rule is in C. |
| `R/zarr-store.R` | The Zarr writer, format 2 (`.zgroup` / `.zattrs` / `.zarray`, `.zmetadata`) and 3 (`zarr.json`, `c/…` chunk keys, consolidated metadata in the root `zarr.json`): little-endian dtypes (each chunk made in one native pass, `src/chunks.c`), `vlen-utf8` strings, the byte shuffle (every array wider than one byte except the expression values; format 3 through the `blosc` codec), deflate via `memCompress()`, edge chunks padded to full size, chunk sizes (2^18 values for arrays read whole, 2^16 for the data and indices of `X`). Three backends: a folder, a `.zarr.zip`, and memory -- uncompressed, no shuffle (it is served over loopback) -- where a chunk is kept as the function that makes it and run when the viewer asks for it; only chunks made by a transform (text, a matrix block, a gene-major block) are kept, in a least-recently-used cache (`cache_bytes`). |
| `R/zip.R` | The `.zarr.zip` container: stored (method 0) entries with UTF-8 names, fixed 1980-01-01 timestamps, ZIP64 records only where needed; the metadata documents after every chunk, the consolidated metadata last. |
| `R/output.R` | `out_path` rules (`.zarr` / `.zarr.zip`, `force`), and writing the store to a temporary sibling that is moved into place once complete. |
| `R/transition-drift.R` | The exported `cellucid_transition_drift()`: `T %*% E - E` from a cell-cell transition matrix, operation for operation cellucid-python's `compute_transition_drift()` (pinned by the conformance fixture the two packages share). |

### The viewer

| File | Responsibility |
|---|---|
| `R/viewer-api.R` | The exported viewer functions: `cellucid_view()`, `cellucid_serve()`, `cellucid_stop()`, `cellucid_stop_all()`, `cellucid_viewers()`, the commands (`cellucid_highlight()`, …, `cellucid_set_data_mapping()`), the events (`cellucid_on()`, `cellucid_off()`, `cellucid_selection()`, `cellucid_wait_ready()`, `cellucid_session()`) and `cellucid_debug_connection()`. An object is planned once and held as a memory store; a path is classified and served from disk as it is. |
| `R/viewer-server.R` | One `httpuv` server per viewer, with the routes of cellucid-python's AnnData server (`/_cellucid/health`, `info`, `datasets`, `anndata/…`, the web build): exact `Range` answers, the memory store's chunks as they are (no gzip), the `Host` (421) and `Origin` (403) checks, and the channel at `/_cellucid/ws` (hello and welcome, a ping every 20 seconds, events delivered to hooks). |
| `R/viewer-protocol.R` | The channel's commands and events, checked exactly (key sets and types), against the shared fixture `tests/testthat/fixtures/viewer_protocol.json`; a highlight and a visibility change go out as binary frames (their cells as little-endian `uint32`, checked and written in one native pass), and selections and session bundles arrive as binary frames. |
| `R/web-cache.R` | The web app the viewer serves: one published build of www.cellucid.com, downloaded at first use into `tools::R_user_dir("cellucid", "cache")/web` and checked file by file against its SHA-256 inventory, with the rules of cellucid-python's web cache; the exported `cellucid_clear_web_cache()`. |
| `R/input-formats.R` | The input-format registry, `inst/input_formats.json`, byte for byte the web app's and cellucid-python's (SHA-256 pinned in all three): every format the viewer opens, and how `cellucid_view(path)` classifies a file or folder. |

### Sessions

| File | Responsibility |
|---|---|
| `R/session-apply.R` | The exported `cellucid_apply_session()`: a saved session's highlight groups and user-defined fields written onto the object as cell-metadata columns, after every check passed; cellucid-python's `apply_cellucid_session_to_anndata()` rule for rule. |
| `R/session-bundle.R` | Reading a `.cellucid-session` file (session format 2): the signature, the manifest, and every chunk. |
| `R/session-chunks.R` | The chunks a session holds, checked before any is applied. |
| `R/session-view.R` | The object as the viewer showed it, through the same extraction `cellucid_view()` serves, so a session is checked against what it was saved on. |
| `R/session-highlights.R`, `R/session-fields.R` | The columns a session's highlight groups and user-defined fields become. |

### Native code and registration

| File | Responsibility |
|---|---|
| `R/native.R`, `src/init.c` | Every `.Call()` of the package in one place, registered in `src/init.c`. |
| `src/crc32.c` | The CRC-32 of a raw vector, for each ZIP entry and each gzip-framed format 3 chunk. |
| `src/sha256.c` | SHA-256, which checks every file of the downloaded web build against its inventory. |
| `src/random.c` | Bytes from the operating system's cryptographic generator: the viewer's id and token, never R's own generator. |
| `src/extent.c` | One-pass finite / float32-range / whole-number checks of large matrices, and the `X/indptr` of an object's genes x cells matrix without transposing it. |
| `src/chunks.c` | A chunk's bytes from an R vector in one pass (canonical NaN, zero padding), one block of the gene-major arrays of a genes x cells `dgCMatrix` (one binary search a cell), a chunk cut from a block, and a binary command frame of the viewer channel. |
| `src/latent_outliers.c` | The outlier quantiles of one column, bit for bit the rule of cellucid-python's `_outliers.py` and the web app's `data/latent-outliers.js` (pinned by the web app's numpy conformance fixture). |
| `src/session_codecs.c` | Inflating a session's chunks (linked against zlib: `src/Makevars` holds `PKG_LIBS = -lz`), decoding its cell sets and user-defined codes, and the web app's cell-set digest. |
| `R/zzz.R` | Stops every viewer's server, then unloads the shared library. |
| `man/*.Rd`, `NAMESPACE` | Hand-written and authoritative; no roxygen. `NAMESPACE` holds `useDynLib(cellucid, .registration = TRUE, .fixes = "C_")`. |

## One `cellucid_prepare()` call, in order

1. **Required arguments.** Missing `out_path` and `format`, and without `data`
   missing `embeddings`, `dataset_id` and `dataset_name`, are reported together.
2. **Output.** `format`, `zarr_format`, `force`, `out_path` (suffix, existing
   target).
3. **Components.** With `data`, the object's components
   (`.dataset_from()`), each replaced by the argument given; `modalities`
   chooses and orders the feature spaces.
4. **Plan.** `dataset_id` and the displayed text, citation and source;
   `obsm` (embeddings, then vector fields bound to `X_<basis>` embeddings of
   the same width); `obs`; each feature space's `X`, gene names (made unique),
   `gene_identifiers` subset and `var_gene_id_column`; the cell names, which
   must agree across inputs; the neighbour graph; the launch defaults and the
   default modality; with `latent_key` the outlier quantiles.
5. **Write.** The plan's `write()` puts the root group, `X`, `obs`, `var`,
   `obsm`, `obsp`, the empty groups and `uns` (or the MuData groups) into
   `.<name>.<pid>.partial` beside the target, the consolidated metadata last.
   Then the existing target (with `force`) is moved aside, the new store
   renamed into place, and the old one deleted. On any error the partial store
   is removed.

Column-type problems (a `Date` column, an infinite value) are found while
writing, so they exercise the cleanup in step 5. `cellucid_view()` runs steps
3 and 4 the same way and writes into a memory store instead.

## Design rules

- The store must read without a warning in Python's anndata (0.8+ encodings)
  and mudata, and open in the web app; `tests/testthat/test-anndata-interop.R`
  and `test-mudata-interop.R` check the first two.
- The store served for an object is the store written for it, key for key and
  byte for byte (`tests/testthat/test-served-store.R`).
- A change to what the store holds is coordinated with `cellucid-python` and the
  web app; see {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`.
