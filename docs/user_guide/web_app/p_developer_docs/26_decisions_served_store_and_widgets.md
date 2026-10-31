# Decision record: served stores and notebook widgets

This page answers: **"How do Python's and R's viewers serve a dataset to the
web app, why is it built that way, what was measured, and what was tried and
replaced?"** It covers `cellucid serve`, `serve_anndata()`, `show_anndata()`
and `cellucid.Widget` in Python, `cellucid_view()` and `cellucid_serve()` in
R, and the web app's side of both. It records decisions and their evidence.
The walkthroughs live on the pages linked at the end. Update this page in the
same change as the code it describes.

## At a glance

**Audience**
- Developers changing the Python or R served store, the
  gene-major layout, either notebook widget, or the web app's served source.

**Time**
- 15 minutes.

**Prerequisites**
- {doc}`20_read_path_architecture` (how the browser reads any store) and
  {doc}`19_subsampling_specification` (priority order, prefixes).

---

## What it does now

1. **The server computes nothing per request.** It names the dataset, serves
   bytes by `Range`, states the launch defaults in `/_cellucid/datasets`, and
   carries the viewer channel over one WebSocket (`/_cellucid/ws`). The
   browser's own AnnData reader reads the store in place, as it reads a
   linked file (`cellucid/assets/js/data/served-anndata.js`). Routes are the
   same in both languages (`cellucid-python/src/cellucid/anndata_server.py`,
   `cellucid-r/R/viewer-server.R`).
2. **What is served as it is, and what goes through a priority store.**
   - *Python.* A file already in priority order (what `prepare()` writes) is
     served as it is. So are R objects and MatrixMarket folders, which the
     browser reads once from the start. An AnnData or MuData object goes
     through `PriorityStore` (`cellucid-python/src/cellucid/_served_store.py`),
     and so do an AnnData file not in priority order, Loom, 10x Cell Ranger
     and MuData files, and any file holding a codec the browser cannot decode.
     The result is a Zarr v2 store in memory, cells in priority order under
     seed 0, with a gene-major `X`. Every chunk is made when the browser asks
     for it.
   - *R.* An R object (Seurat, SingleCellExperiment, SummarizedExperiment,
     matrices) becomes the same store, held in memory
     (`cellucid-r/R/zarr-store.R`, `.serve_object` in
     `cellucid-r/R/viewer-api.R`). A file or Zarr folder given by path is
     always served from disk as it is.
3. **One gene-major layout pass.** `prepare()`, a backed `prepare()` and the
   Python served store share one layout
   (`cellucid-python/src/cellucid/prepare_data/_gene_major.py`,
   `GeneMajorLayout`). It runs a counting pass in the source's stored order,
   then a placing pass that writes each gene's entries in ascending served
   row with no global sort. R does the same in C
   (`cellucid-r/src/chunks.c`, shared kernels in
   `cellucid-r/src/gene_major.h`). The layout goes into memory-mapped
   scratch files (`cellucid-python/src/cellucid/_scratch.py`; R maps files
   it unlinks at once), so the matrix is never held twice in memory.
4. **The layout runs in the background, `X` first.** The store, its metadata
   and every other element are ready at once, so the viewer opens and draws
   without waiting. Python runs the layout on the store's worker thread. R
   runs it on native threads, one segment per thread: by cells for a
   genes × cells matrix, by columns of about equal values for a renumbered
   one. Each thread publishes its progress every 256 cells
   (`LAYOUT_STEP`).
5. **A chunk asked for before its matrix is laid out costs only the part not
   laid out yet.** The values already placed are copied. The rest come from
   one scan of the source's stored values in their own order: sequential
   memory, on every core in R (`gene_major_scan`), and in `_scanned` in
   Python. The result has the same bytes the layout will write. R's request
   holds the layout threads while it copies (`cellucid_layout_block`).
   Python's pass is not held, because the GIL interleaves the two anyway.
6. **The ski-rental wait.** A request waits for the layout instead of making
   its chunk when, at the layout's own pace, the seconds left are no more
   than the last request took, or than all requests together took
   (`sooner()` in Python, `layout_sooner` in C). A burst of genes therefore
   costs at most twice what waiting for the layout would.
7. **Python file matrices before counting.** Where each gene's values go is
   known only after the counting pass. Until that pass ends, keys of a file's
   matrix are answered `503` with `Retry-After` (`NotReady`). For an object,
   counting runs before the store opens, so a value that cannot be stored is
   refused at once.
8. **An object is served from a snapshot.** In Python, while the worker lays
   out or copies an array of the object, numpy marks that array read-only
   (`_ReadOnly`): an in-place change is refused instead of mixing into what
   is served. In R, copy-on-change gives the same result. A change made
   after viewing does not show; view the object again.
9. **Chunk sizes for loopback.** Served expression entries come in chunks of
   2^18 values (`SERVED_ENTRY_CHUNK_ELEMENTS`, `.chunk_elements_served`) and
   are uncompressed, so a `Range` of a chunk makes only those bytes. Files
   that `prepare()` writes keep 2^14-value compressed chunks.
10. **`cache_bytes`.** This is what a served store keeps of what it makes,
    least recently used first out. The default is 1/64 of physical memory
    (256 MiB on a 16 GB machine): at most the share the web app gives its
    decoded chunks (a 32nd of what its working set leaves of half the
    device). `0` keeps nothing. Both languages use the same default
    and the same check (`default_cache_bytes`, `.default_cache_bytes`).
11. **Dataset identity.** A given `dataset_id`/`dataset_name` wins, then
    `uns['cellucid']` (in R, read only from a Zarr store; R has no HDF5
    reader), then the file's name. The id is the file name without its
    extension, made portable by one shared rule (`_label_id`, `.label_id`).
    When an id is made rather than given, it is printed in one line together
    with how to set it, because sessions are kept under it. Objects are named
    after a label: in R the call's expression (`deparse(substitute(data))`),
    in Python `AnnData (<n> cells)`.
12. **The browser keeps nothing of a served store across visits.**
    `servedResource` marks a served file `keptParts: false`, and
    `formats/io/kept-parts.js` refuses loopback and the widget's host
    (`sameMachine`). The server already reads the data from beside it, so a
    copy saves no download, and a store made from an object gets a new tag on
    every run.
13. **Notebook widgets.** The page opens with `?viewer=embedded` inside the
    output. Every request and the viewer channel go through the widget's
    comm (`cellucid/assets/js/data/viewer-embed.js`). URLs use the reserved
    origin `https://host.cellucid.invalid`
    (`cellucid/assets/js/formats/io/host-fetch.js`, `routedFetch`), so
    `Range`, `206` and `Content-Range` behave as they do over HTTP.
    - *Python* (`cellucid-python/src/cellucid/jupyter/_widget.py`,
      `_Answers`): each request is answered on a thread pool, so a slow one
      holds back no other. Answers that are ready together share one message
      of at most 8 MiB, with at least 2 ms between messages
      (`_ANSWER_MESSAGE_INTERVAL`, half of Jupyter Server's limit of 1,000
      messages a second).
    - *R* (`cellucid-r/R/viewer-widget.R`): bodies and binary frames go as
      binary buffers of the Jupyter message, not base64. The page keeps at
      most two fetch messages outstanding (`_fetch_messages` = 2); requests
      asked meanwhile join the next message. The answers' JSON is written in
      one pass.
    - Both packages ship the same front end, byte for byte
      (`cellucid-python/src/cellucid/jupyter/static/widget.js`,
      `cellucid-r/inst/widget/widget.js`), and both test suites pin its
      digest.

## Decisions and why

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-28 | One gene-major layout for `prepare()`, backed `prepare()` and the served store, in both languages | Same bytes everywhere (3,628 files identical before and after); one place to optimise | A separate writer and server layout; `lexsort` of the whole matrix |
| 2026-10-28 | Layout in a background worker, `X` first; `503` + `Retry-After` only for a file matrix not yet counted | The viewer opens without waiting for a pass over every value | Blocking the first view on the layout |
| 2026-10-28 | Object arrays read-only while laid out or copied | An in-place change must not mix into what is served | Copying the object before serving it |
| 2026-10-28 | The browser keeps no parts of same-machine or served stores; after a hidden tab released memory, a reopen reads loopback again (measured: 15 requests, 4.2 MB, no network) | A copy duplicates data already on this disk; object stores are new every run | A browser cache of served bytes |
| 2026-10-28 | No per-run secret on served data routes; the docs state that other users on a shared machine can read served data | A secret in a URL leaks into logs and sessions; the limit is stated, not hidden | A token in the data route |
| 2026-10-28 | `show_anndata()` names the dataset and makes its id like R does, and prints the id it made | Sessions are keyed by the id; a made id must be visible | Silent random or fixed ids |
| 2026-10-30 | A chunk asked for during the layout copies what the layout placed and makes only the rest, from one scan in stored order | Cost is the share left, and nothing near the end | A second full pass per chunk |
| 2026-10-30 | Ski-rental wait (last request, or all requests together) | A burst costs at most twice the layout's remaining time | Unbounded repeated passes |
| 2026-10-30 | R lays out on every core, one segment per thread | 10M × 10: 7.56 s on 1 thread, 1.15 s on 8 | A single-threaded layout |
| 2026-10-30 | R widget: binary buffers, at most two fetch messages outstanding, one-pass JSON | IRkernel spends about 1.7 ms of its own on each message; base64 cost 60 MB/s on R's thread | One message per request; base64 bodies |
| 2026-10-30 | Python widget: answers on a thread pool, batched, rate-bounded | One 2.74 s answer blocked the opening's 560 requests; Jupyter drops messages past 1,000 a second | Answering in order on the kernel thread |
| 2026-10-30 | One pass size for every gene-major pass: `PASS_ENTRIES` in `prepare_data/_gene_major.py`, sized so a step's temporaries (about 40 B a value) take a 64th of physical memory (6.7M values on 16 GB), imported by `_served_store.py`; one `machine_memory_bytes()` (`_scratch.py`), and in R one `.machine_memory_bytes()` (`R/native.R`) for the default cache and the backed passes | The served store used 2^23 and the layout 2^22 with no reason; layout time is flat from 2^20 to 2^24 at 3M × 10 (0.66–0.69 s) and the bytes are identical | Two constants for one concept |
| 2026-10-30 | R backed matrices (BPCells, DelayedArray) stay on disk: one counting pass, then placing passes of at most ¼ of physical memory, 64-bit positions | The old path joined every block into one double `dgCMatrix` (12–24 B a value, int32 pointers) | The 2^31 − 1 value cap; a dgCMatrix of the whole matrix |
| 2026-10-30 | Speed and less download outrank memory, as long as memory is counted and bounded | Project decision | Trading speed for memory that is already explained |
| 2026-10-30 | App timings are taken on a real GPU only | SwiftShader holds the main thread on synchronous GL calls | Timing on a software rasteriser |

## Measurements that justify it

All on an Apple M1 Pro (16 GB), synthetic 10M cells × 2,000 genes, 10 values
per cell, 2026-10-30, unless stated. R timings taken with
`pkgload::load_all()` run C compiled at `-O0`. An installed package is
`-O2`; one chunk from a parked 4M × 10 layout took 33–38 ms at `-O0` and
9 ms at `-O2`.

**A gene asked during the layout, measured in process (no browser):**

| | Before | After |
|---|---|---|
| R layout alone | 6.8 s (one thread) | 7.56 s (1 thread) → 1.15 s (8) / 1.19 s (10) |
| R first gene during layout | 2.66–2.98 s | 0.64 s (layout at 0 %); same chunk next 0 ms; others 0.22–0.30 s |
| R 100 genes during layout | 4.9 s | 1.91 s (to the layout's end, 3.10 s) |
| R written block (`cellucid_prepare()`) | 4.73 s (1 thread) | 1.05 s (10 threads), same bytes |
| Python first gene during layout | 2.02–2.20 s | 0.13 s; later ones 0.13–0.26 s |
| Python 100 genes during layout | (layout already over) | 2.06 s (to the layout's end, 2.59 s) |

At 2M × 10 (light load), the R widget through a real IRkernel took 450 small
requests from 1.12 s (round-trip median 629 ms) to 0.26 s (262 ms). Moving
1M cells' embedding and obs (18.9 MB) went from 0.41 to 0.17 s.

**In Chromium on Metal, gene 1 asked right after the first picture (one run
each):**

| Path | First picture | Gene 1 / next | 100 genes |
|---|---|---|---|
| Python `serve_anndata`, earlier run | 2.58 s | 120 / 96–131 ms | 0.38 s |
| Python `serve_anndata`, this design | 2.60 s | 0.87 s (during layout) / 67–69 ms | 0.22 s |
| R `cellucid_serve`, earlier run | 2.30 s | 1.49 s / 78–101 ms | 1.46 s |
| R `cellucid_serve`, this design | 1.84 s (983,040 cells) | 1.03 s / 113–483 ms | 0.48 s |
| Python widget, earlier run | 6.09 s (3.18M) | 80 / 53 ms | 0.65 s |
| Python widget, this design | 4.91 s (2.29M); every cell 8.81 s | 97 / 72 ms | 1.90 s |
| R widget, earlier run | 11.96 s (2.46M); every cell 26.2 s | 72 / 46 ms | 0.59 s |
| R widget, this design | 12.87 s (1.64M); every cell 21.45 s | 79 / 92 ms | 0.53 s |

The R widget's opening round-trip median went from 974 to 708 ms (Python:
176 ms). Earlier widget figures (100 genes 5.8 s in Python and 7.7 s in R,
22.5 s and 42 s to every cell, a 30× gap to HTTP) were measured on
SwiftShader and do not describe the app. On Metal, before the widget fixes,
100 genes took 0.62 s in the Python widget and 2.13 s in the R widget.

**R backed matrices:**

| Case | Result |
|---|---|
| 1M cells × 2.2·10⁹ values (BPCells, 1.59 GB), `cellucid_prepare()` | 413 s, peak RSS 7.26 GB, 1.22 GB written, verified gene by gene |
| 1M × 4·10⁸ object, served | laid out in 17.8 s, peak 5.4 GB (3.2 GB of it mapped scratch) |
| 100k × 2.2·10⁸, prepared | 64 s / 3.8 GB (through a dgCMatrix: 67 s / 6.1 GB) |
| Real 10x 1.3M neurons, BPCells from its `.h5`, 345 genes (4·10⁸ values), prepared | 324 s, peak 6.9 GB, 423 MB written |
| The same file, served | started in 0.13 s, first keys and a range in 0.36 s |
| The same object, 93 genes (1.2·10⁸ values), served | store ready after 250 s (laid out on R's thread), peak 3.7 GB |
| Parallel deflate, 200k × 200 prepare | 24.2 → 7.7 s; deflate 39 → 295 MB/s on 1 → 10 threads |

**Loopback request cost**, which sets the 2^18-value chunk: Python's
`http.server` answers at most about 4,700 requests a second, and R's handlers
about 1,600. A browser's `fetch` from a worker moves 46–79 MB/s in 64 KiB
requests and 300–510 MB/s in 1 MiB requests (Chromium, Firefox, WebKit).

## Tried and reverted

- **A second pass per chunk during the layout (R, before 2026-10-30).**
  Each asked chunk was made by a single-threaded pass over every cell beside
  the layout, 1.3–2.98 s per gene. It was replaced by progress sharing.
- **Waiting by one comparison (R, 2026-10-30).** An interim rule compared
  one pass with the layout's rest, so a burst paid one pass per gene
  (100 genes 4.9 s). It was replaced by the two-condition rule.
- **Making the rest in output order (R, round 1 of the shared design).** A
  pass over each segment's rest in output order touched random columns, with
  two cache misses per cell: 0.2–0.3 s per chunk at 10M on ten threads. The
  scan in stored order is about 10× cheaper and replaced it.
- **Python `entries_now` ignoring progress.** It laid the chunk's genes out
  over every row block, so a request took as long as the layout itself
  (0.42 s per gene at 2M; 2.0–2.2 s at 10M). It was replaced by copying what
  was placed plus `_scanned`.
- **Base64 bodies in the R widget.** These took 17.6 ms per MiB on R's
  thread; binary buffers take 0.3 ms. One message per request was replaced
  by the two-message window. Earlier versions of both widgets also put a
  two-batch gate in the front end, which was removed.
- **Joining backed blocks into a dgCMatrix (R).** `.in_memory_matrix()` was
  replaced by `.readable_matrix` and `cellucid-r/R/backed-matrix.R`.
  Python's dask sparse `X` is now decoded by row blocks into scratch instead
  of `.compute()`d whole.
- **To try progress sharing differently:** byte identity must hold. The
  randomized stress test (262 blocks at random moments, 6 thread counts,
  both kinds) compares each block with the transpose.

## Known limits and open items

- **A served R backed matrix is laid out on R's thread before the viewer
  opens.** Only R's thread may read BPCells or DelayedArray. For the real
  1.3M-neuron object above, that took 250 s. A served object above 2^31
  values has not been run (it needs about 17.6 GB of scratch).
- **R serves files as they are.** A file not in priority order is not
  re-laid out by R the way Python's `PriorityStore` does it, and R cannot
  read an `.h5ad`'s `uns['cellucid']` identity.
- **Object labels differ:** R uses the call's expression, Python
  `AnnData (<n> cells)`. Both print the id they made.
- **Dense dask arrays** are still computed whole in Python.
- **Notebook kernels run comms between cells.** While a cell runs, the
  viewer waits; a call that would wait on the viewer says so at once.
  Jupyter Server drops kernel messages beyond 1,000 a second over 3 s (comm
  messages count): a real platform limit, which is why answers are batched.
- **The R widget's opening is slower than Python's** (first picture 12.9
  against 4.9 s at 10M; round-trip median 708 against 176 ms). The Python
  widget's 100 genes are bound by the notebook page's relay (0.62–0.67 s
  against 0.18 s over HTTP). In the latest single run the Python widget's
  100 genes and DE read (1.90 s, 0.86 s) were slower than in the earlier run
  (0.65 s, 0.33 s). This is not yet explained and needs repeated runs.
- The 10M browser rows are single runs, and the 10M rows were taken before
  the request's scan in stored order replaced the pass in priority order
  (that change measured only at 2M and 4M: one chunk from a parked 4M × 10
  layout 0.20–0.30 s before, 33–38 ms at `-O0`, 9 ms at `-O2`). The 10M run
  of the final design (R direct, served and widget at `-O0` and `-O2`;
  Python direct and served) has not been taken yet: re-measure before citing
  the 10M rows as the final design's.
- The thread start and join code in `cellucid-r/src/chunks.c` (`run_tasks`)
  is a fourth copy beside two in `src/backed.c` and one in `src/deflate.c`;
  one helper in `src/gene_major.h` for all of them is open work.

## Known issues (pre-commit)

Defects, regressions and unverified behaviour known when this version was
committed (2026-10-30), each with where it is recorded above. An item is
removed in the change that fixes or verifies it.

- **10M rows:** taken before the request's scan in stored order; the final design's 10M run (R at `-O0` and `-O2`, Python) is not taken.
- **Python widget:** 100 genes and the DE read were slower in the latest single run (1.90 s, 0.86 s) than in the earlier one (0.65 s, 0.33 s), not yet explained.
- **R served backed matrix:** laid out on R's thread before the viewer opens (250 s for the 1.3M-neuron object).

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| Python server routes, identity, `cache_bytes` | `cellucid-python/src/cellucid/anndata_server.py` | `cellucid-python/tests/test_anndata_server_store_contract.py`, `cellucid-python/tests/test_serve_defaults_and_reporting_contract.py` |
| Python priority store, worker, `NotReady` | `cellucid-python/src/cellucid/_served_store.py`, `cellucid-python/src/cellucid/_served_formats.py` | `cellucid-python/tests/test_served_priority_store.py`, `cellucid-python/tests/test_served_foreign_formats.py` |
| Gene-major layout, progress sharing, `sooner()` | `cellucid-python/src/cellucid/prepare_data/_gene_major.py` | `cellucid-python/tests/test_served_layout_progress.py`, `cellucid-python/tests/test_prepare_anndata_contract.py` |
| Source arrays, scratch, file versions | `cellucid-python/src/cellucid/_source_arrays.py`, `cellucid-python/src/cellucid/_scratch.py`, `cellucid-python/src/cellucid/_file_versions.py` | `cellucid-python/tests/test_served_priority_store.py` |
| Python widget | `cellucid-python/src/cellucid/jupyter/_widget.py`, `cellucid-python/src/cellucid/jupyter/static/widget.js` | `cellucid-python/tests/test_notebook_widget_contract.py`, `cellucid-python/tests/test_round_trip_widget_contract.py` |
| R viewer, server, identity | `cellucid-r/R/viewer-api.R`, `cellucid-r/R/viewer-server.R` | `cellucid-r/tests/testthat/test-viewer-server.R`, `cellucid-r/tests/testthat/test-served-priority.R` |
| R served store, layout jobs | `cellucid-r/R/zarr-store.R`, `cellucid-r/R/expression-input.R`, `cellucid-r/R/native.R` | `cellucid-r/tests/testthat/test-served-store.R`, `cellucid-r/tests/testthat/test-served-layout-progress.R` |
| R layout threads, scan, wait rule | `cellucid-r/src/chunks.c`, `cellucid-r/src/gene_major.h`, `cellucid-r/src/init.c` | `cellucid-r/tests/testthat/test-served-layout-progress.R` |
| R backed matrices | `cellucid-r/R/backed-matrix.R`, `cellucid-r/src/backed.c` | `cellucid-r/tests/testthat/test-backed-matrix.R` |
| R widget | `cellucid-r/R/viewer-widget.R`, `cellucid-r/inst/widget/widget.js` | `cellucid-r/tests/testthat/test-viewer-widget.R`, `cellucid-r/tests/testthat/test-round-trip-widget.R` |
| Web served source, viewer host, channel | `cellucid/assets/js/data/served-anndata.js`, `cellucid/assets/js/data/viewer-host-source.js`, `cellucid/assets/js/data/viewer-channel.js` | `cellucid/tests/served-anndata-contract.test.mjs`, `cellucid/tests/viewer-host-source-contract.test.mjs` |
| Widget transport (web) | `cellucid/assets/js/formats/io/host-fetch.js`, `cellucid/assets/js/data/viewer-embed.js` | `cellucid/tests/notebook-widget-requests-contract.test.mjs`, `cellucid/tests/viewer-embed-contract.test.mjs`, `cellucid/tests/browser/notebook-widget-jupyterlab.spec.mjs` |
| No kept parts for served stores | `cellucid/assets/js/formats/io/kept-parts.js` (`sameMachine`, `keptPartsOf`) | `cellucid/tests/kept-parts-contract.test.mjs` |
| Real-GPU launch for timings | `cellucid/scripts/lib/browser-launch.mjs` | `cellucid/tests/browser-launch-contract.test.mjs` |

## Related pages

- {doc}`20_read_path_architecture` and {doc}`09_data_loading_pipeline_and_caching`
  (the `remote` and `jupyter` sources, caching layers)
- {doc}`19_subsampling_specification` (priority order, seed 0, prefixes)
- {doc}`23_decisions_partial_reading` (share costs, kept data and revisits)
- {doc}`27_decisions_memory_model` (how the page counts what it holds)
- {doc}`21_decisions_decoders_and_codecs` (what the browser decodes, and why
  a file with an undecodable codec goes through the priority store)
- {doc}`01_codebase_map_and_entry_points`
- {doc}`../../python_package/h_developer_docs/10_jupyter_embedding_architecture`,
  {doc}`../../python_package/h_developer_docs/16_performance_profiling_and_scaling`,
  {doc}`../../python_package/e_jupyter_hooks/05_architecture_viewer_channel`,
  {doc}`../../python_package/d_viewing_apis/14_performance_scaling_and_lazy_loading`
- {doc}`../../r_package/h_developer_docs/01_codebase_architecture`
