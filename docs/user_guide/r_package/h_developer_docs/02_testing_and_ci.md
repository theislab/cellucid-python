# Testing and CI

**Audience:** contributors  
**Time:** 10 minutes  
**What you'll get:** how to run the suite locally, what it covers, and what CI checks

## Running the tests

From the `cellucid-r/` directory:

```r
install.packages(c("testthat", "jsonlite", "Matrix"))
testthat::test_local()
```

The full check, as CRAN runs it (the package compiles one C file, so a C
toolchain is needed):

```sh
R CMD build .
R CMD check --as-cran cellucid_<version>.tar.gz
```

`Matrix` is a `Suggests` dependency; the tests that need it call
`skip_if_not_installed("Matrix")`, so install it to run them.

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
| `test-prepare.R` | The written store, read back with an independent base-R reader (`helper-store.R`): root and element encodings; `obs` column encodings (categorical, string array, float64 with NaN, bool, nullable boolean and integer); `X` as float32 CSC with int32 indices, stored zeros dropped, float64 for integer counts above 2^24; category-code types at pandas' thresholds; chunk sizes of `X` and the graph; `format` and the `out_path` it requires; `obsm`, vector fields, `obsp`; `uns['cellucid']` with only the defaults that were set, and the outlier quantiles a `latent_key` stores; consolidated metadata; the `.zarr.zip` holding the same bytes with every entry stored; chunk padding and `vlen-utf8` strings; dense gene-major `X`; unique names; `gene_identifiers`; every validation error; the write-then-move and `force` behaviour |
| `test-anndata-interop.R` | Python's anndata reads both outputs without a warning, value for value |
| `test-cross-writer.R` | The same input written by `cellucid_prepare()` and by Python's `cellucid.prepare()` (through `tests/testthat/python/cross_writer.py`) gives identical stores: every `.zgroup`, `.zattrs` and `.zarray` document (encodings, dtypes, shapes, chunks, codecs, fill values) and every value. Runs with the Python named by `CELLUCID_PYTHON` (else `python3`) when it imports `cellucid`, `anndata` and `zarr`; skipped otherwise and on CRAN |
| `test-latent-outliers-conformance.R` | The outlier-quantile kernel against the web app's numpy conformance fixture (`fixtures/latent_outlier_conformance.json`, the same SHA-256 pinned in both repositories), byte for byte; how each kind of `obs` column is partitioned |
| `test-transition-drift.R` | `cellucid_transition_drift()` against the conformance fixture shared, byte for byte, with `cellucid-python` (`fixtures/transition_drift_conformance.json`, written by cellucid-python's `scripts/generate_transition_drift_conformance.py` with exact rational arithmetic; the same SHA-256 pinned in both repositories): exact cases bit for bit, general cases within the forward-error bound the fixture states, each with the matrix dense and sparse; every input form; the drift stored by `cellucid_prepare(vector_fields = )`; every refusal |
| `test-data-mapping-conformance.R` | The make-unique rule and its note against the conformance fixture shared, byte for byte, with `cellucid-python` |
| `test-current-contract.R` | Required files in a checkout and in the built tarball, no retired wording, README links, workflow branches, `DESCRIPTION` agreeing with every in-repo version site, pkgdown navigation, and each hand-written `man/` page matching the real formals of its function (`cellucid_prepare()`, `cellucid_transition_drift()`) |
| `test-citation-contract.R` | The generated citation regions match `CITATION.cff` |
| `test-size-limit-allowlist.R` | Every size constant in `R/` and `src/` is a real format or platform limit, with its reason |

## CI

`.github/workflows/R-CMD-check.yaml` runs `R CMD check` on macOS (release),
Windows (release), and Ubuntu (devel, release, oldrel-1) for every push and
pull request to `main`. The runners' `python3` normally has no `anndata`, so
the round trip is skipped there; run it locally.

`release.yaml` builds the source tarball and runs `R CMD check --as-cran` on it
for a published release; see {doc}`03_release_process` and
`cellucid-r/publishing.md`.
