# Testing and CI

**Audience:** contributors  
**Time:** 10 minutes  
**What you'll get:** how to run the suite locally, what it covers, and what CI checks

## Running the tests

From the `cellucid-r/` directory:

```r
install.packages(c("testthat", "jsonlite", "httpuv", "later", "curl", "Matrix", "SeuratObject"))
testthat::test_local()
```

The full check, as CRAN runs it (the package compiles a few C files and links
zlib, so a C toolchain is needed):

```sh
R CMD build .
R CMD check --as-cran cellucid_<version>.tar.gz
```

`Matrix`, `SeuratObject` and the Bioconductor packages (`SingleCellExperiment`,
`SummarizedExperiment`, `S4Vectors`, `DelayedArray`) are `Suggests`
dependencies; the tests that need one skip without it, so install them to run
those tests.

### The anndata round trip

`tests/testthat/test-anndata-interop.R` writes a small dataset as a `.zarr`
folder and as a `.zarr.zip`, reads each back with Python's `anndata`
(`tests/testthat/python/read_with_anndata.py`, every warning an error), and
compares every value with what R put in: cell and gene names, every `obs`
column, `X`, every `obsm` array, the graph, and `uns`. It uses the Python named
by the `CELLUCID_ANNDATA_PYTHON` environment variable, else `python3` on the
`PATH`, provided it can `import anndata, zarr`; otherwise, and on CRAN, it is
skipped.

```sh
CELLUCID_ANNDATA_PYTHON=/path/to/venv/bin/python \
  Rscript -e 'testthat::test_local(filter = "anndata")'
```

It also writes integer counts above 2^24 and reads them back exactly, as
float64.

### The cross-writer check

`tests/testthat/test-cross-writer.R` writes several inputs (every kind of
column, stored zeros, counts above 2^24, dense and sparse expression, only
embeddings) with `cellucid_prepare()`, describes the same R objects in JSON,
and has `tests/testthat/python/cross_writer.py` write them with Python's
`cellucid.prepare()` and compare the two stores document for document and
value for value (floating values byte for byte, so a NaN must be the same NaN). It needs a Python that imports `cellucid`, `anndata` and
`zarr`, named by `CELLUCID_PYTHON` (else `python3`):

```sh
CELLUCID_PYTHON=/path/to/venv/bin/python \
  Rscript -e 'testthat::test_local(filter = "cross-writer")'
```

Run both whenever you change what the store holds, in either writer.

## What the suite covers

| File | What it checks |
|---|---|
| `test-prepare.R` | The written store, read back with an independent base-R reader (`helper-store.R`): root and element encodings; `obs` column encodings (categorical, string array, float64 with NaN, bool, nullable boolean and integer); `X` as float32 CSC with int32 indices, stored zeros dropped, float64 for integer counts above 2^24; category-code types at pandas' thresholds; chunk sizes of `X` and the graph; `format` and the `out_path` it requires; `obsm`, vector fields, `obsp`; `uns['cellucid']` with only the defaults that were set, and the outlier quantiles a `latent_key` stores; consolidated metadata; the `.zarr.zip` holding the same bytes with every entry stored; chunk padding and `vlen-utf8` strings; format 3's `blosc` codec on shuffled arrays, each frame (header, block size, blocks, one stream per byte plane, whole copies) the one c-blosc writes; dense gene-major `X`; unique names; `gene_identifiers`; every validation error; the write-then-move and `force` behaviour |
| `test-anndata-interop.R` | Python's anndata reads both outputs without a warning, value for value |
| `test-cross-writer.R` | The same input written by `cellucid_prepare()` and by Python's `cellucid.prepare()` (through `tests/testthat/python/cross_writer.py`) gives identical stores: every `.zgroup`, `.zattrs` and `.zarray` document (encodings, dtypes, shapes, chunks, codecs, fill values) and every value. Runs with the Python named by `CELLUCID_PYTHON` (else `python3`) when it imports `cellucid`, `anndata` and `zarr`; skipped otherwise and on CRAN |
| `test-latent-outliers-conformance.R` | The outlier-quantile kernel against the web app's numpy conformance fixture (`fixtures/latent_outlier_conformance.json`, the same SHA-256 pinned in both repositories), byte for byte; how each kind of `obs` column is partitioned |
| `test-transition-drift.R` | `cellucid_transition_drift()` against the conformance fixture shared, byte for byte, with `cellucid-python` (`fixtures/transition_drift_conformance.json`, written by cellucid-python's `scripts/generate_transition_drift_conformance.py` with exact rational arithmetic; the same SHA-256 pinned in both repositories): exact cases bit for bit, general cases within the forward-error bound the fixture states, each with the matrix dense and sparse; every input form; the drift stored by `cellucid_prepare(vector_fields = )`; every refusal |
| `test-data-mapping-conformance.R` | The make-unique rule and its note, the vector-field key rule and the default modality against the conformance fixture shared, byte for byte, with `cellucid-python` and the web app |
| `test-current-contract.R` | Required files in a checkout and in the built tarball, no retired wording, README links, workflow branches, `DESCRIPTION` agreeing with every in-repo version site, pkgdown navigation, and each hand-written `man/` page's `\usage` matching the real formals of every function it documents |
| `test-dataset-from.R` | `cellucid_prepare(data = )` with Seurat (SeuratObject 4 and 5; v5 split layers joined on a copy; several assays as a MuData), SingleCellExperiment (reducedDims, altExps, colPairs, size factors, a DelayedArray assay read block by block) and SummarizedExperiment objects; the object's own record; arguments replacing components; every refusal with its fix |
| `test-plan-robustness.R` | The plan where an object departs from the simple case: an assay holding its cells in another order or only some of them, a factor with an `NA` level, a dense genes x cells assay, a matrix held outside R's memory read in many blocks |
| `test-apply-session.R` | `cellucid_apply_session()` against the session file the web app's own encoders wrote (`fixtures/session_format2/`, a byte-identical copy of cellucid-python's), so both packages are held to one session and one expected result; every refusal |
| `test-served-store.R` | The store `cellucid_view()` serves from memory is the prepared store, key for key and byte for byte, for explicit inputs, Seurat, MuData and SingleCellExperiment; an object's genes x cells matrix stored as the same cells x genes matrix; an object changed after serving leaves the served bytes unchanged |
| `test-mudata-interop.R` | Python's mudata reads R's Cellucid-ready MuData, value for value (`CELLUCID_MUDATA_PYTHON`) |
| `test-viewer-protocol.R` | The viewer channel's declarations against the shared fixture `fixtures/viewer_protocol.json` (the web app's `data/viewer-protocol.js`, SHA-256 pinned in all three repositories): every example event accepted, every example command produced, refusals, a two-million-cell selection decoded at once |
| `test-viewer-server.R` | The server end to end on free loopback ports: routes and documents, the web build, every store key against the prepared store, exact ranges (206, 416, suffix, `If-Range`, 304), the chunk cache budget, the `Host` (421) and `Origin` (403) checks, `allowed_hosts`, the channel's hello/welcome, 4401, 4409, 1008 and 1001, commands queued until the viewer connects, events to hooks, a million-cell selection, session bundles, stop and restart on the same port, several viewers, a taken port, files and folders by path, a file's listed size and `ETag` (the one its bytes are answered with) and a new one once it is rewritten in place (its modification time put back), replaced or gone |
| `test-viewer-robustness.R` | The server and channel off the happy path: `HEAD` answers, long ranges, stale answers, cell indices out of range, sockets that never say hello or keep sending after a violation, the package unloading, the random credentials |
| `test-web-cache.R` | SHA-256 against FIPS 180-4 vectors; the inventory schema; the build id; a cache checked file by file; a build downloaded from a local source, kept, replaced, refused when a file is wrong, and used offline |
| `test-input-formats.R` | The packaged registry `inst/input_formats.json` against the SHA-256 pinned in the three repositories; paths classified as the web app classifies them |
| `test-citation-contract.R` | The generated citation regions match `CITATION.cff` |
| `test-native-encoding.R` | An R session whose native encoding is not UTF-8, as Windows R before 4.2 is (`LC_CTYPE` switched to `C` for a block): the same store bytes and notes as a UTF-8 session, exact condition text, every condition raised as an object (checked in the parsed sources), repeated names, JSON and browser frames read as UTF-8, a memory store's keys, an object's column named in Greek served; in any session, file names handed to httpuv in UTF-8, code-point order of categories, and the sizes the C code refuses before they could wrap a 32-bit `R_xlen_t` |
| `test-size-limit-allowlist.R` | Every size constant in `R/` and `src/` is a real format or platform limit, with its reason |
| `test-minimum-r-version.R` | Every base R function the package, its tests, `inst/CITATION` and the vignettes call existed in the oldest R `DESCRIPTION` admits (`fixtures/r-4.1-api.txt`, written by R 4.1), and no syntax is newer (the pipe placeholder `_` is R 4.2) |

## CI

`.github/workflows/R-CMD-check.yaml` runs `R CMD check --as-cran` for every
push and pull request to `main`: on Ubuntu with R 4.1, 4.2, 4.3, oldrel-1,
release and devel, and on macOS and Windows with R 4.1 and release. Every job
installs all suggested packages (SeuratObject, SingleCellExperiment and the
other Bioconductor ones) except three: one Ubuntu job with none of them, one
with SeuratObject 4.1.4, and the macOS and Windows R 4.1 jobs, which check with
the hard dependencies plus knitr, rmarkdown and testthat (the build and check
tools, installed explicitly rather than taken from the runner image). Windows R
4.1 checks both its sub-architectures, 32-bit i386 (where `R_xlen_t` is an
`int`) and x64, in a CP1252 session. CRAN no longer serves R 4.1 Windows binaries, so
the Bioconductor stack would build RCurl from source without curl headers, and
Matrix 1.6-5 has no R 4.1 arm64 binary and does not link from source there; the
suggested-package paths are covered at R 4.1 on Ubuntu. The R 4.1 and 4.2
Ubuntu jobs pin Matrix 1.6-5, the last release for R older than 4.4, which
SeuratObject 5 needs. The Ubuntu release job installs cellucid-python, anndata
and mudata, so the cross-writer and interoperability tests run there; elsewhere
they skip, naming what is missing.
Bioconductor packages come from bioconductor.org (`R_BIOC_MIRROR`): pak's
default mirror, bioconductor.posit.co, keeps only Bioconductor 3.19 and later,
and R 4.1, 4.2 and 4.3 need 3.14, 3.16 and 3.18.

`release.yaml` builds the source tarball and runs `R CMD check --as-cran` on it
for a published release; see {doc}`03_release_process` and
`cellucid-r/publishing.md`.
