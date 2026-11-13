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

The core run tests each part of the package once (eight files, listed in
`tools/test-core.R`: preparing from the objects people hold, the viewer
channel, the served store, sessions, the data-mapping and subsampling
fixtures, the generated actions, and the comment-history check). Run it after
every change; `testthat::test_local()` is the full run:

```sh
Rscript tools/test-core.R
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

The MCP tests wait for the server's start and its end on those events alone,
never on a clock: they are R's start, its package load and its teardown, as
long as the computer makes them, and an antivirus starved by CPU load once
stretched a start past a minute on a developer's Mac (by 2026-11-06). The
other waits are bounded by what the test sends, and their deadline only
turns a hang into a failure.

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
`test-served-cross-language.R` does the same for what `cellucid_view()` and
Python's `show_anndata()` serve, with the same `CELLUCID_PYTHON`.

## What the suite covers

| File | What it checks |
|---|---|
| `test-prepare.R` | The written store, read back with an independent base-R reader (`helper-store.R`): root and element encodings; `obs` column encodings (categorical, string array, float64 with NaN, bool, nullable boolean and integer); `X` as float32 CSC with int32 indices, stored zeros dropped, float64 for integer counts above 2^24; category-code types at pandas' thresholds; chunk sizes of `X` and the graph; `format` and the `out_path` it requires; `obsm`, vector fields, `obsp`; `uns['cellucid']` with only the defaults that were set, and the outlier quantiles a `latent_key` stores; consolidated metadata; the `.zarr.zip` holding the same bytes with every entry stored; chunk padding and `vlen-utf8` strings; format 3's `blosc` codec on shuffled arrays, each frame (header, block size, blocks, one stream per byte plane, whole copies) the one c-blosc writes; dense gene-major `X`; unique names; `gene_identifiers`; every validation error; the write-then-move and `force` behavior |
| `test-anndata-interop.R` | Python's anndata reads both outputs without a warning, value for value |
| `test-cross-writer.R` | The same input written by `cellucid_prepare()` and by Python's `cellucid.prepare()` (through `tests/testthat/python/cross_writer.py`) gives identical stores: every `.zgroup`, `.zattrs` and `.zarray` document (encodings, dtypes, shapes, chunks, codecs, fill values) and every value. Runs with the Python named by `CELLUCID_PYTHON` (else `python3`) when it imports `cellucid`, `anndata` and `zarr`; skipped otherwise and on CRAN |
| `test-latent-outliers-conformance.R` | The outlier-quantile kernel against the web app's numpy conformance fixture (`fixtures/latent_outlier_conformance.json`, compared with the web repository's file), byte for byte; how each kind of `obs` column is partitioned |
| `test-transition-drift.R` | `cellucid_transition_drift()` against the conformance fixture shared, byte for byte, with `cellucid-python` (`fixtures/transition_drift_conformance.json`, written by cellucid-python's `scripts/generate_transition_drift_conformance.py` with exact rational arithmetic; the copy compared with cellucid-python's file): exact cases bit for bit, general cases within the forward-error bound the fixture states, each with the matrix dense and sparse; every input form; the drift stored by `cellucid_prepare(vector_fields = )`; every refusal |
| `test-data-mapping-conformance.R` | The make-unique rule and its note, the vector-field key rule and the default modality against the conformance fixture shared, byte for byte, with `cellucid-python` and the web app |
| `test-current-contract.R` | Required files in a checkout and in the built tarball, workflow branches, `DESCRIPTION` agreeing with every in-repo version site, the sibling repositories CI fetches, the pkgdown site contract (`.github/scripts/validate-pkgdown.R`: every vignette and help page indexed once), each hand-written `man/` page's `\usage` matching the real formals of every function it documents, and R code being ASCII outside comments |
| `test-dataset-from.R` | `cellucid_prepare(data = )` with Seurat (SeuratObject 4 and 5; v5 split layers joined on a copy; several assays as a MuData), SingleCellExperiment (reducedDims, altExps, colPairs, size factors, a DelayedArray assay read block by block) and SummarizedExperiment objects; the object's own record; arguments replacing components; every refusal with its fix |
| `test-backed-matrix.R` | A matrix held outside R's memory (BPCells by cells and by genes, DelayedArray over a sparse and over a dense backend) written and served byte for byte as its `dgCMatrix` (folder and `.zarr.zip`, priority order or not, kept genes), with blocks of a few cells and passes of a few values; the passes counted (count, then one or two a share, each over its genes' rows; never the whole matrix); float64 counts decided on the genes kept; `NaN`, `NA`, infinities and out-of-range values refused at counting, on every Bioconductor release; a cells x genes backed `gene_expression` (a DelayedArray read under its transposition, each block transposed in memory, never through DelayedArray's transposed view, which Bioconductor 3.14 to 3.16 cannot read with matrixStats 1.x); a Seurat v5 BPCells layer; a matrix changed between passes refused (skipped without BPCells or DelayedArray) |
| `test-plan-robustness.R` | The plan where an object departs from the simple case: an assay holding its cells in another order or only some of them, a factor with an `NA` level, a dense genes x cells assay, a matrix held outside R's memory read in many blocks |
| `test-category-colors.R` | Cellucid's own category colors (`R/category-colors.R`), which an apply writes as a user field's or a recolored dataset column's colors when nothing else gives a category one, against the shared fixture `fixtures/category_default_colors.json` (the web app's `scripts/generate-category-default-colors.mjs`; the copy compared with the web repository's file), byte for byte |
| `test-sessions.R` | Session format 4 against the shared golden history (`fixtures/history_golden/`: `history_golden.cellucid-session` and `history_golden.expected.json`, written by the web app's `scripts/generate-history-golden.mjs` with the history model's own constructors, byte-identical in the three repositories: the copies here, `inst/extdata/`'s included, are compared with the web repository's files), case for case cellucid-python's `tests/test_sessions.py`: the header, every tracker, branch, step, note, turn and journal entry, the state at a step, diffs, the methods paragraph and pseudocode, records found by name or id; the shared apply plans of a sequence of applies, applying again changing nothing, `keep_session`, a column the person made never replaced, a plan refused once its object changed, the cells and the dataset checked, only what the chosen step uses required, labels, user fields (written with the colors they are shown in, the view in focus's), the id an object was shown under (`dataset_id = `, a viewer's own request), highlight columns named by page and label, a subsample's cells, `overwrite`, gene sets, presentation reported and never applied, a Seurat object and a data frame; and, on files `helper-session-files.R` makes, a steps-only share (read, never applied), damaged files, a newer Cellucid's file (newer, never damaged) and another format (the web app's sentence) |
| `test-session-layout.R` | The chunk layout the web app's writer makes (`app/session/history-file.js`), on files `helper-session-files.R` lays out the same way, case for case cellucid-python's `tests/test_session_layout.py`: parts stored once and joined once, a part that names itself through others, a value that would expand past 2^29 − 24 characters (refused before it is expanded), a part nothing refers to, a part renamed (read) or changed (refused by the value's hash), a shown value naming a lazy part, malformed part references; steps, parts and values in many pieces; chunks stored as they are (`codec` `none`); parameters naming a lazy array; a broken layout; the constants against the web app's source |
| `test-history-model.R` | `inst/history_model.schema.json` byte for byte the web repository's `tests/fixtures/history-model.schema.json`; `R/history-model-generated.R` equal to what `tools/generate-history-model.R` writes now; every record of the schema and every field exposed by the reading model; names compared as the web app compares them; the constants; a record that breaks the model refused naming the path and the rule |
| `test-history-words.R` | The methods text and the pseudocode for what the golden branches do not reach: a squashed step of each cause (recording paused, steps merged, the history cleared, a state shared), the sentence that a squashed step holds only its result, everyone who acted while paused, and the Before/After group headings; the same strings as cellucid-python's `tests/test_history_words.py` and the web app's exports |
| `test-viewer-history.R` | The history live, against a viewer played by a WebSocket client in this process: the actions inside `cellucid_turn()` carrying its id and label (and turns that do not nest), a session file pushed with `cellucid_load_session()` as one binary frame answered by `session_loaded`, every `history_changed` kept (`cellucid_last_history_change()`) and handed to the hooks, and the history's actions as the generated `cellucid_history_*()` functions |
| `test-actions.R` | The action catalog (`inst/action_catalog.json`) byte for byte the web repository's `tests/fixtures/action-catalog.json`; `R/actions-generated.R`, the generated help pages, the NAMESPACE export block and `_pkgdown.yml`'s group of the actions' pages equal what `tools/generate-actions.R` writes now (run it after the web app's catalog changes); no generated name clashes with a hand-written export; parameters checked and converted (cells, values, omitted and `NULL`), the request answered by `action_result` or `command_error`; a change answered without a step keeps why (`request$unrecorded`) and warns it when collected (`cellucid_action_not_recorded`), a step or a read warns nothing |
| `test-served-store.R` | The store `cellucid_view()` serves from memory is the prepared store, key for key and byte for byte, for explicit inputs, Seurat, MuData and SingleCellExperiment; an object's genes x cells matrix stored as the same cells x genes matrix; an object changed after serving leaves the served bytes unchanged |
| `test-mudata-interop.R` | Python's mudata reads R's Cellucid-ready MuData, value for value (`CELLUCID_MUDATA_PYTHON`) |
| `test-viewer-protocol.R` | The viewer channel's declarations against the shared fixture `fixtures/viewer_protocol.json` (the web app's `data/viewer-protocol.js`; the copy compared with the web repository's file): every example event accepted (`action_result`, `history_changed` and `session_loaded` included), every example command produced (each binary action and `load_session` byte for byte, an action's `turn` and `turnLabel` together or neither), refusals, a two-million-cell selection decoded at once |
| `test-viewer-server.R` | The server end to end on free loopback ports: routes and documents, the web build, every store key against the prepared store, exact ranges (206, 416, suffix, `If-Range`, 304), the chunk cache budget, the layouts' scratch files let go when the viewer stops, the `Host` (421) and `Origin` (403) checks with each refused socket closed (1008), `allowed_hosts`, the channel's hello/welcome, 4401, 4409, 1008 and 1001, commands queued until the viewer connects, events to hooks, a million-cell selection, session bundles, stop and restart on the same port, several viewers, a taken port, files and folders by path, a file's listed size and `ETag` (the one its bytes are answered with) and a new one once it is rewritten in place (its modification time put back), replaced or gone |
| `test-viewer-robustness.R` | The server and channel off the happy path: `HEAD` answers, long ranges, stale answers, cell indices out of range, sockets that never say hello or keep sending after a violation, the package unloading, the random credentials |
| `test-web-cache.R` | SHA-256 against FIPS 180-4 vectors; the inventory schema; the build id; a cache checked file by file; a build downloaded from a local source, kept, replaced, refused when a file is wrong, and used offline |
| `test-input-formats.R` | The packaged registry `inst/input_formats.json` byte for byte the web repository's `tests/fixtures/input-formats.json`; paths classified as the web app classifies them |
| `test-citation-contract.R` | The generated citation regions match `CITATION.cff` |
| `test-native-encoding.R` | An R session whose native encoding is not UTF-8, as Windows R before 4.2 is (`LC_CTYPE` switched to `C` for a block): the same store bytes and notes as a UTF-8 session, exact condition text, every condition raised as an object (checked in the parsed sources), repeated names, JSON and browser frames read as UTF-8, a memory store's keys, an object's column named in Greek served; in any session, file names handed to httpuv in UTF-8, code-point order of categories, and the sizes the C code refuses before they could wrap a 32-bit `R_xlen_t` |
| `test-size-limit-allowlist.R` | Every size constant in `R/` and `src/` is a real format or platform limit, with its reason |
| `test-minimum-r-version.R` | Every base R function the package, its tests, `inst/CITATION` and the vignettes call existed in the oldest R `DESCRIPTION` admits (`fixtures/r-4.1-api.txt`, written by R 4.1), and no syntax is newer (the pipe placeholder `_` is R 4.2) |
| `test-priority-order.R` | The priority order (`R/priority.R`): cells stored by ascending rank under the file's seed with every per-cell element, prefix levels stepping down by 3/4 as cellucid-python's, a genes x cells assay counted by cell, `priority_seed = NULL`, the record, a record that does not describe the cells, a file prepared again keeping its original positions, and `cellucid_restore_original_order()` |
| `test-subsample-conformance.R` | `cellucid_subsample_cells()` against the shared fixture `fixtures/subsample_conformance.json` (the web app's; the copy compared with the web repository's file): ranks, sizes, a fraction's shortest decimal (ties to the even digit, subnormals, powers of two), every case's cells and digest, the session encoding of a cell set, nesting, code-point category order, and an object in a file's stored order refused until restored |
| `test-round-trip.R` | `cellucid_pull()`, `cellucid_push()` and subsets against a viewer played by a WebSocket client in this process: every highlight, label and score onto a Seurat object, pruning, a new page's results, another dataset's refused, an object whose cells are not the viewer's in order refused, a subsampled viewer's unopened cells left `NA`, and the web app's exports pulled onto a data frame |
| `test-round-trip-widget.R` | The same round trip in a Jupyter notebook against stand-in comms: nothing waits for the viewer, and a session asked in one cell is collected in a later one |
| `test-viewer-widget.R` | The viewer in a notebook (`viewer-widget.R`) against stand-in comms: the anywidget model with no server and no credentials, the document and store reads by range, the channel over the comm, a cell that cannot wait, a second output replacing the first, and the widget's front end byte for byte cellucid-python's |
| `test-viewer-rerun.R` | A notebook cell run again stops the viewer it showed before (its matrices and chunk cache released), unless another cell shows it; without a cell id, showing a dataset again stops its earlier viewer and says so; a served store keeps a 64th of the machine's memory by default, as cellucid-python's |
| `test-served-priority.R` | The store `cellucid_view()` and `cellucid_serve()` serve holds its cells in priority order, the arrays `cellucid_prepare()` writes; each cell keeps its original position; the expression matrix and the neighbor graph are served in priority order without a copy |
| `test-served-layout-progress.R` | A served expression matrix laid out gene-major in the background on every core: a chunk asked for before the layout ends is made from its progress (values already laid out copied, only the rest passed over), counted, never timed. |
| `test-gene-sample.R` | The gene sample as R reads it: a session's `genes` record and the `color_by` event's `cells_read`, checked as the web app and cellucid-python check them. |
| `test-mcp.R` | `cellucid_mcp()`: `inst/mcp_tools.json` byte for byte the web repository's `tests/fixtures/mcp-tools.json` (the web connector's list for this package's catalog), the shared fixture `local_agents.json` (byte for byte cellucid-python's), the agent commands and events, a viewer's agent port, its publication and endpoint, and the server over stdio against a fake page (initialize, lists, calls, errors, cancellation, the person asked in the app, a tool list that never changes, a viewer it opens itself); and cellucid-python's `cellucid mcp` driving an R viewer. On Windows CI too (the stdin code). |
| `test-served-cross-language.R` | The same object served by R and by Python holds the same bytes in every array both serve, metadata and every chunk, with float64 values, stored zeros and duplicate coordinates stored alike (`CELLUCID_PYTHON`) |
| `test-count-noun-agreement.R` | A count and its noun agree in every message: `.plural()`, and the package's parsed sources refused where a plural noun follows a count |
| `test-total-counts.R` | `obs['total_counts']` and `uns['cellucid']['totalCounts']` (`R/total-counts.R`): the sums over the stored genes in every matrix form (sparse, dense, genes x cells, a DelayedArray), after `gene_identifiers`; a matrix that holds no counts recording no totals; the caller's column (`total_counts =`) recorded over all genes or, with `total_counts_over =`, over a source panel, and each refusal (below the sum of the genes stored, a missing value, not a number, not a column of obs, a matrix without counts); the data's own `total_counts` kept; one record per counts modality of a multimodal object, named by `total_counts` and `total_counts_over` as a vector; an object's own record kept. cellucid-python's `tests/test_prepare_total_counts_contract.py` pins the same and `test-cross-writer.R` the same bytes |
| `test-viewer-runtime.R` | What opened a viewer, as each entry point names it: the shared protocol declares both R runtimes (`r`, `r-mcp`), `cellucid_serve()` and `cellucid_view()` name `r` and `cellucid_mcp(open = )` names `r-mcp` in `_cellucid/info`, in the agents record and in the agent port's description |
| `test-viewer-token-in-output.R` | The viewer's token never reaches what is logged or pasted: `.public_viewer_url()` drops it wherever it stands, as Python's does; every line of `R/` that names a viewer's address is classified, so a new one is seen; a printed viewer and its errors say the address without it; `cellucid_view(pane = "none")` says the address with its token once, to the person |
| `test-comment-history.R` | The comments of `R/`, `tests/` and `tools/` and of `src/`'s C, and the descriptions of the tests, hold no history: no CI run address, finding or plan id, bug-hunt or round story, or pointer into the private plans (the rule is "Comments say what the code does" in `CONTRIBUTING.md`); a first test checks that the checker finds each shape planted in R, C and `test_that()` text and ignores a string |


Two files in `tests/testthat/` run before the tests and hold no test.
`setup-agents.R` points `CELLUCID_AGENTS_DIR` at a folder of this test run's
own, so every viewer a test starts is published to agents there and never to
the person's computer; the child processes (`cellucid_mcp()` over stdio,
Python's `cellucid mcp`) inherit it. `helper-viewer.R` is how a test talks to a
viewer's server from the test process, whose event loop is the server's: HTTP
through curl's multi interface and the WebSocket through a non-blocking socket
(`http_request()`, `http_json()`, `ws_open()` and the frame helpers), both
turning the event loop between polls, so a client never blocks the server it
asks; `viewer_fixture_web_cache()` is the small web build in `fixtures/web`
(nothing is downloaded) and `test_viewer()` starts a viewer on a free port from
it, stopped when the calling test ends.

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
SeuratObject 5 needs. The Ubuntu release job installs anndata, mudata and cellucid-python's
`main` (`CELLUCID_PYTHON_BRANCH` in the workflow's job `env`), so the
cross-writer, served cross-language and interoperability tests run there;
elsewhere they skip, naming what is missing. The two packages write and serve
the same bytes at their current sources, so the job never installs a PyPI
release: when a change touches what both write or serve, push cellucid-python
first (until its `main` has the change, those tests fail here and say what
differs).
Every job also fetches the siblings the cross-repository tests compare with
(`tests/testthat/helper-repository.R`), each at `main`: the web app (its
modules, generators and shared fixtures, fetched after the dependencies are
installed so a day's pushes of the three repositories have landed) and
cellucid-python (its package sources and fixtures), outside the package
directory, named by `CELLUCID_WEB_ROOT` and `CELLUCID_PYTHON_ROOT` (under
`R CMD check` nothing is beside the sources), with Node installed and
`CELLUCID_REQUIRE_SIBLINGS=1`: there a missing sibling or Node fails those
tests instead of skipping them, as it skips them only where no sibling is
beside the package (a local checkout alone, an installed package, a CRAN
check).

The web repository's workflow ("Web app tests") tests this package from the
other side. Its Node job checks out the files of this repository its contracts
read (the shared files and the widget copy) and, under the same
`CELLUCID_REQUIRE_SIBLINGS=1`, fails where one is missing. Its `viewers` job
installs this package from the checkout (`R CMD INSTALL`, with its hard
dependencies and SingleCellExperiment) and runs the web specs that drive the
viewer through it (`r-viewer-channel`, `r-viewer-history`) in Chromium, with
`CELLUCID_RSCRIPT` naming the `Rscript` and `CELLUCID_REQUIRE_VIEWERS=1`
making a missing interpreter a failure, not a skip. The variables are described
in {ref}`the web app's developer docs <sibling-repositories-and-environment-variables>`.

Every file this package shares is compared byte for byte with its home
through one helper, `expect_shared_file()` (`helper-repository.R`): the web
repository's file for what the three repositories share (the action catalog,
`mcp_tools.json`, the history model, the input formats, the viewer protocol,
the conformance fixtures, the category colors, the golden
history, the round-trip fixtures), cellucid-python's for what only the two
packages share (the notebook widget's front end, `local_agents.json`, the
transition-drift fixture, the MuData cross-writer harness). No test holds a
digest of one: a difference names the one command that makes every copy
agree, `node scripts/sync-shared-files.mjs` in the web repository
({ref}`the shared files <shared-files>`).
Bioconductor packages come from bioconductor.org (`R_BIOC_MIRROR`): pak's
default mirror, bioconductor.posit.co, keeps only Bioconductor 3.19 and later,
and R 4.1, 4.2 and 4.3 need 3.14, 3.16 and 3.18.

### What CI sees that a local run does not

The CI run of 2026-10-31 failed on 8 of
12 jobs while `testthat` and `R CMD check --as-cran` passed on a current
macOS. Each failure came from a difference between the machines, and each
is now covered where a local run can reach it:

- **Old Bioconductor.** Up to Bioconductor 3.18 (R 4.1 to 4.3), DelayedArray's
  coercion of a block to a `dgCMatrix` went through a dense block and
  `which(block != 0)`, which drops `NA` and `NaN`: a missing value in an
  HDF5Array or other DelayedArray assay was written as a zero, and the test
  that expects it refused failed on those jobs only. A block is now read with
  `DelayedArray::read_block()` (sparse when its backend is) and converted
  without that step. To run the suite on those releases, make an Intel
  conda environment (Apple silicon runs it under Rosetta), for example
  `CONDA_SUBDIR=osx-64 mamba create -p ./r41 --platform osx-64 -c conda-forge -c bioconda "r-base=4.1" "bioconductor-delayedarray=0.20.0" bioconductor-singlecellexperiment r-testthat r-jsonlite r-httpuv r-later r-curl r-matrix`
  (`r-base=4.3` with `bioconductor-delayedarray=0.28.0` for Bioconductor
  3.18), install testthat 3.2.1 or later from CRAN into it, and install the
  package with `SDKROOT` set to `xcrun --show-sdk-path` and `MAKE` running
  `/usr/bin/make` as arm64 (`arch -arm64`), since the Intel `make` launcher
  cannot load Xcode's arm64-only `libxcrun`.
- **Windows' C library.** Its `printf` breaks a decimal tie as it likes: the
  subsample size of a fraction came from its shortest decimal, found with
  `snprintf()`, and 2^-25 was written `…313` instead of the even `…312`, so
  2^-25 of 50,331,648 cells came out 2. The shortest decimal is now found in
  exact integer arithmetic (`src/subsample.c`), checked against Python's
  `repr()` on about 50,000 doubles (every power of two in (0, 1] and its
  neighbors, subnormals, short decimals, random bit patterns).
- **Windows' files.** An open file keeps its name there, so the layouts'
  scratch files (`FILE_FLAG_DELETE_ON_CLOSE`) are listed in the session folder
  until the viewer stops; the test asserted the Linux and macOS behavior
  everywhere.
- **R 4.1 on Windows.** Rtools 4.0's toolchains declare no Windows API newer
  than XP by default, so `GetTickCount64()` did not link in the 32-bit build;
  the layouts' clock is now the performance counter. Every C file can be
  compiled for 32- and 64-bit Windows at that API level from any machine with
  `zig cc -target x86-windows-gnu -D_WIN32_WINNT=0x0502` (and
  `x86_64-windows-gnu`), with the R and zlib headers on the include path.
- **httpuv.** After `onHeaders` answers a WebSocket handshake (403, 421),
  httpuv writes its own 101 and opens the socket all the same; the server now
  closes such a socket as it opens (`.on_server_ws_open()`). A test that
  closed its client before that open ran made httpuv skip the open but still
  deliver the close, whose handler warns "object '…' not found" in whichever
  test turns the event loop next (seen on Windows); the tests now wait for
  the server's close.
- **httpuv's writes.** httpuv sends a frame as a header, a payload and an
  empty last buffer, and libuv leaves a write whose last buffer is empty
  unfinished until the next turn of its loop; a frame sent right after
  another waits for that turn, and stopping the server first drops it.
  `cellucid_stop()` sent the freeze and the close together and stopped at
  once: the page got the freeze, then the connection's end with no close
  (macOS release; every time under Rosetta, where it reproduces locally).
  It now waits for the page to answer the close (`.viewer_close_seconds`,
  one second, bounds a page that never answers) before stopping.
- **cellucid-python's `main` behind.** The job installs cellucid-python's
  `main`, which on 2026-10-31 did not yet hold the served store's
  `SERVED_ENTRY_CHUNK_ELEMENTS` that cellucid-r's push was written against.
  A change to what both packages write or serve is pushed to cellucid-python
  first (above).

Later runs showed two more (both by 2026-11-12):

- **Old Bioconductor, the transposed view.** Up to Bioconductor 3.16
  (DelayedArray 0.20 to 0.24, R 4.1 and 4.2), a sparse block read through
  DelayedArray's transposed view subsets a `SparseArraySeed` with
  `MatrixGenerics::rowAnyNAs(useNames = NA)`, which matrixStats 1.x, the one
  CRAN serves those R versions, refuses ("Argument 'useNames' must be either
  TRUE or FALSE"). `cellucid_prepare()` of a cells × genes DelayedArray failed
  on those jobs only. Blocks are now read from the matrix under the
  transposition and transposed in memory (`.delayed_untransposed()`), never
  through the transposed view.
- **The site's contract.** `.github/scripts/validate-pkgdown.R` (every
  vignette and help page indexed once in `_pkgdown.yml`) gated only the
  site's build (`pkgdown.yaml`): the hand-kept reference list of the actions'
  pages kept two removed pages and missed a new one, which failed the site's
  deploy and nothing before it. The suite now runs the same script
  (`test-current-contract.R`, "the pkgdown site contract holds"), so it fails
  `testthat` and `R CMD check` too.

Two checks once passed without checking anything (both by 2026-08-02):

- **Tests that stopped running under `R CMD check`.** There the tests'
  repository root resolved to the `.Rcheck` directory, where none of
  `DESCRIPTION`, `NEWS.md`, `README.md`, `NAMESPACE`, `man/` or `R/` exists.
  Every assertion in `test-current-contract.R` sat inside an
  `if (file.exists(...))` or `if (dir.exists(...))`, so forty-two assertions
  about the release identity stopped running while the suite stayed green;
  deleting a page, a workflow or a config did the same. The two layouts (a
  checkout, a tarball under check) are now recognized by which one matches,
  an unrecognized layout is an error, and the required paths are declared
  (`cellucid_required_paths()`) and asserted in both directions.
- **The release workflow's ordering check.** It searched for the bare phrase
  `R CMD check --as-cran`; the file opens with a comment naming it and
  `regexpr()` reports the first match, so the ordering held even with the
  check step renamed out of existence, and an absent step (position −1)
  satisfied it too. The check is anchored on the step headers, and both
  positions must be found.

Two warnings in those runs are not the package's: loading
SummarizedExperiment warns "replacing previous import
'S4Arrays::makeNindexFromArrayViewport' by
'DelayedArray::makeNindexFromArrayViewport'" with the Bioconductor binaries
those jobs install (R 4.3 with S4Arrays 1.2.1 and DelayedArray 0.28.0, and on
oldrel-1, release and devel), in the first test that loads it; DelayedArray
re-exports that function from S4Arrays, and S4Arrays 1.2.0 with DelayedArray
0.28.0 (R 4.3), or 1.12.1 with 0.38.2 (R 4.5), load without the warning. The
server's log line
"a request for /assets/app.js failed: Web asset 'assets/app.js' byte length
is 9" is the test that tampers with the web build on purpose.

`release.yaml` builds the source tarball and runs `R CMD check --as-cran` on it
for a published release, then attaches it to the release and submits it to
CRAN; see {doc}`03_release_process` and
`cellucid-r/publishing.md`.
