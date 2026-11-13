# Testing, CI, and release process

This page is the reference for **running and writing tests** of the Cellucid
web app: the commands, what CI runs, the rules a contract test and a browser
spec follow, the device every spec runs on, and the checklists before a
change is shared. Why each rule exists (which run showed what, what was
tried, the numbers) is in {doc}`38_decisions_testing_and_hosted_runners`.

Cellucid is a high-performance interactive app; many regressions are not caught by “does it run?” checks.
The goal here is to make validation repeatable.

## At a glance

**Audience**
- Contributors: run the suites, write a test that holds on every runner, and follow the manual checklist before opening PRs or sharing builds.
- Maintainers: use the release checklist to avoid “it worked locally” deployments.

**Time**
- 20 minutes to read once; then a page to come back to.

**Prerequisites**
- {doc}`02_local_development_setup`.

---

## Validation contract

- The web app is shipped as static files (no required bundler build).
- `npm run check` type-checks every module of `assets/js` from its JSDoc
  (TypeScript's compiler as a development tool: `checkJs`, nothing emitted,
  about 11 s). It fails on the kinds of diagnostic that mean a defect
  whatever the typing style (an unknown name, a wrong argument count, a
  comparison that can never be true, a `@private` member reached from
  outside, an unused name), on any diagnostic in a folder declared fully
  checked, and on a directive comment that silences the checker. The codes
  and the fully checked folders are listed in `scripts/check-types.mjs`;
  `tests/type-check-contract.test.mjs` runs it inside `npm test`, so CI runs
  it with no step of its own.
- `npm test` runs the complete Node contract suite.
- `npm run test:browser` runs the Playwright browser suite.
- `npm run test:core` runs the core contract tests: for each layer of the
  application, the files that run that layer's modules and say the most for
  their time (about 15 seconds). It is what to run after a change, while
  working.
- `npm run test:browser:smoke` runs the smoke specs, one walk from start to
  export (start, open each format, color, filter, select, analyze, undo,
  save, export, an agent's turn), in every engine or the one named
  (`-- --project=chromium`).
- The two short sets only answer first: `npm test` and the full browser
  suite remain the gate. Both are lists in `scripts/test-suites.mjs`, held
  to the tree by `tests/test-suites-contract.test.mjs`; the table "What to
  run when" in the web repository's `CONTRIBUTING.md` (Testing & validation)
  says which to run at each point.
- `npm run test:browser:repeat` runs the browser specs a change touches many
  times in every engine; see
  [Browser specs that wait for state, not time](#browser-spec-stability).
- `npm run test:worker-bundle` checks that the community-annotation worker
  (a Cloudflare Worker, not a web worker) bundles (a `wrangler deploy --dry-run`).
- `npm test` also refuses web worker bundles (`assets/workers/`) that are not
  exactly what `node scripts/generate-web-asset-inventory.mjs` builds from the
  sources. The browser suite runs the app from `127.0.0.1`, where workers
  start from their sources; `tests/browser/deploy-skew.spec.mjs` opens it as a
  deployed site, and
  `CELLUCID_BROWSER_TEST_APP_HOST=deployed.localhost npm run test:browser -- <specs>`
  runs specs that way, every worker started from its bundle (a spec that
  holds a module a worker imports, such as `dataset-load-supersede`, needs
  the sources, and runs from `127.0.0.1`).
- `npm test` refuses in the same way a code fingerprints module
  (`assets/js/data/anndata/transcode/kept-code-fingerprints.js`, which the
  generator writes before the bundles) that is not what the code gives
  (`tests/kept-code-fingerprints-contract.test.mjs`;
  {doc}`09_data_loading_pipeline_and_caching`).
- UI and state changes also require the manual smoke checks below.

`npm test` runs `scripts/run-web-tests.mjs`, which discovers every
`tests/*.test.mjs` file, so a new contract file is picked up by adding it.
Each file runs in a process group of its own, and the runner stops every file
it started when it is interrupted (Ctrl-C, a canceled CI job). A file is
bounded by its progress, not its duration: it is stopped only after finishing
no test for 5 minutes (`CELLUCID_TEST_STALL_MS` sets another bound), and then
fails naming the last test it finished.

In CI (the workflow "Web app tests", `cellucid/.github/workflows/test.yml`),
`npm test` runs on Ubuntu, macOS, and Windows across Node 22, 24, and 26
(Node.js 22.15 or newer is required). Every push and pull request runs the
browser suite in Chromium, Firefox, and WebKit on Ubuntu, plus WebKit on
macOS, in eight shards each; every engine on every operating system runs
nightly and on manual dispatch, except WebKit on Windows (Playwright's
Windows WebKit is the WinCairo port, not Safari). A failed browser job
uploads its Playwright screenshots and error context as a workflow artifact.
A further job, `viewers`, runs the Python and R viewer specs against the
packages installed from their checkouts ({ref}`sibling-repositories-and-environment-variables`).

Each browser job runs one shard through `scripts/run-browser-shard-bounded.mjs`,
which refuses a shard whose files need more browser processes than its
runner allows. A test belongs to the spec file whose top-level suite holds
it, not to the file its `test()` call is written in (the control walk's four
files declare theirs through `helpers/control-walk.mjs`), and before each
batch runs the runner lists exactly what it will pass to Playwright and
refuses a selection other than the plan.
`tests/browser-shard-capacity-contract.test.mjs` plans every shard of the
matrix for every engine and runner from one Playwright listing:
`tests/helpers/shard-plan-reporter.mjs` replicates Playwright 1.65.0's
`createTestGroups` and `filterForShard` (a group weighed by its tests not
expected to be skipped), refuses another version, and is compared test for
test with real `--shard` listings of every other shard and the last, which
hold every boundary.

(contract-test-rules)=
### Contract tests that hold on every Node version and operating system

The contract suite runs on three operating systems, three Node versions and
runners of very different speed. A contract test follows these rules; the
failure behind each is in
{ref}`the record <testing-decisions-contract-suite>`.

- **A file path from a module URL** is
  `fileURLToPath(new URL('…', import.meta.url))`, never `.pathname`: on
  Windows that is `/D:/a/…`, which Node opens as `D:\D:\a\…`.
- **Every `.js` of the repository is an ES module** (`package.json` declares
  `"type": "module"`), the token validators in `scripts/` included.
- **An exact power of ten** is parsed from its literal, ``Number(`1e${n}`)``:
  `Math.pow(10, n)` is not exact on every Node version.
- **A test waits for the event it needs**, never a fixed pause: timers and
  message delivery are slower and coarser on Windows than on Linux and
  macOS. Every wait in the suite is one of four kinds:
  - *For delivery or settling* (a message crossing a `MessageChannel` or
    `BroadcastChannel`, a worker's answer, a file written, work the code
    starts on its own): wait for what proves it — an event, a promise the
    code returns, or the state it changes, asked again until it holds
    (`until(condition, what)` in `tests/helpers/until.mjs`). A port delivers in order, so a question
    asked after a message is answered after it: a test that must read the
    receiving side's state once everything sent before has been handled
    sends one and awaits the answer (`pageHandled` in
    `formats-io-contract`, `quiet()` in `notebook-widget-requests-contract`).
  - *A turn after work that is only promise continuations* (fakes, no I/O,
    timers or channels), or after a timer the code set first with an equal
    or shorter delay: Node runs these in a fixed order, so
    `setTimeout(resolve, 0)` is exact there.
  - *A passage of time the test is about* (a debounce, an interval, an idle
    timeout, a token's expiry, a simulated round trip): kept, with its
    assertion bounded only as the product's contract is.
  - *"Nothing happens"*: no event proves an absence, so the test first waits
    for what ends the work that could do it (the answer to a later message,
    the read that settles it, past the product's own delay), then asserts.
- **A test waits for the work, not for a number of turns of the event loop.**
  Crypto, zlib and file work run on other threads, so no count of `setImmediate`
  turns is enough for them. A fake schedule hands the write's own promise back
  (`manualSchedule` in `tests/helpers/memory-history-database.mjs`); a test with
  mocked timers ticks them once the code has set the timer it means to end.
  Virtual time (`tests/helpers/virtual-time.mjs`) counts zlib's thread-pool
  work, so its clock does not move while a chunk is still inflating.
- **A link the test gives a reader is priced as the test means it**, round
  trip included: the page cache is given one round trip of the link to
  answer, so a link of 1 ms lets the runner's own speed decide what is read
  again. Price bytes dear on a slow link instead.
- **A test whose assertions follow from a plan holds this machine's speed.**
  An opening's plan weighs the speed of the machine besides the link: how
  fast chunks decode (`progressive-open.js` `decodingRate`, averaged over
  everything the process decoded), how long making every cell's order takes
  (`subsample-plan.js` `geneSampleOf`), how fast a kept gene index was read
  and what a gene index pass's work cost. Measured, each follows how busy
  the machine is. Every such rate is read through one module,
  `assets/js/utils/machine.js` (`measureOf`), and held there by a test
  (`holdMachine`). Such a test file imports
  `tests/helpers/machine-at-rest.mjs` for its effect (`MACHINE_AT_REST`):
  decoding held at `DECODE_BYTES_PER_MS_AT_REST` (1 GB/s, fast enough never
  to bind these small stores' plans), the order made in no time worth
  pricing, no kept read and no index pass measured yet. A test of a
  measurement itself leaves it measured
  (`holdMachineAtRest({decodeBytesPerMs: undefined})`) and puts the rest back
  when it ends. The link is not held by that import: those files model
  their links themselves (a clock of their own, a rate set on the meter,
  answers delayed). A test whose host is a function in the test's own
  process, with no time of its own, measures only how busy the machine is,
  and holds the link too (`holdMachineAtRest({link: LINK_AT_REST})`). A byte
  reader whose requests a test holds at a gate takes a clock of its own
  (`createByteReader(source, {clock: virtualTime().now})`): timed, a request
  held for as long as a busy machine takes to open the gate is a slow link's.
- **A test undoes what it set up in reverse order.** `node:test` runs a
  test's `after` hooks first registered first, so a test that fakes
  `fetch` (or `Worker`, `window`, `document`) and then makes a source with
  it would put the real global back before the source is torn down. Such a
  test registers its undoings with `undoLater(t, undo)`
  (`tests/helpers/teardown.mjs`), which runs them last registered first,
  each awaited.
- **No contract test reaches the network.** `npm test` preloads
  `scripts/network-guard-preload.mjs` into every test file's process: a
  request to any host but this machine's loopback (`fetch`, and every TCP
  connection beneath it) is refused before a DNS lookup, written to stderr
  with the stack that asked, and the file that made it fails when its
  process exits, naming each URL -- also when the code caught the refusal
  and retried, or asked after its test had ended. A test fakes every host it
  reads. The one opt-in is for a test whose subject is a real host (none is
  today): `allowRealHost(t, host)` (`tests/helpers/real-network.mjs`) lets
  that host through for that test only. A file run on its own takes the
  guard with `node --import ./scripts/network-guard-preload.mjs --test
  tests/<file>.test.mjs`; with the slow-delivery preload, `NODE_OPTIONS`
  adds to it (`network-guard-contract` checks the guard itself).
- **A contract about cost counts, it does not time.** A threshold on a clock
  is flaky by construction on a shared runner. `blocksRun`
  (`tests/helpers/block-counts.mjs`) counts the blocks of JavaScript a piece of
  work runs, with V8's precise coverage: the same on every machine and every
  run, and a loop per item adds a count per item. A contract says how that
  count grows with the input (a thousand times the bytes adds a few blocks per
  chunk; four times the steps, about four times the blocks). A native pass runs
  no block, so what matters of one is counted too: the input lengths native
  functions were given (`community-annotation-worker-cpu-contract`), ranges
  shifted (`kept-parts-contract`), questions asked and the length of the timer
  each set (`formats-io-contract`). A test's own `timeout` is a fence for a hang,
  never a budget.
- **Starting a process is measured, never assumed**: a test of a time bound
  times a short run first and scales from it
  (`tests/web-test-runner-contract.test.mjs`). **A rate is judged only on a
  watchdog tick that came on time** (`scripts/install-linux-engine-libraries.mjs`).
  **A process started in a group of its own** (`detached`) does not hear the
  terminal's Ctrl-C or a canceled job's SIGTERM: whoever started it forwards
  them.

Two runs find a test that breaks these rules before CI does:

- **Slow delivery and coarse timers.**
  `NODE_OPTIONS="--import ./scripts/slow-delivery-preload.mjs" npm test`:
  every port and channel message arrives up to 25 ms late, in order, and
  every timer runs on a 16 ms tick, in order.
- **A slow clock.** Run the suite with every measured interval made longer:
  preload a module that replaces `performance.now` by
  `base + (now - base) * N` (`node --import ./slow-clock.mjs --test …`, N =
  40 and 400). A test that passes at 1 and fails at 40 asserts something the
  machine's speed decides. (A deadline a test measures on the same clock
  expires N times sooner: read the failure before believing it.)

### Running two browser suites at once

The browser suite binds two loopback ports: `4173` for the application and
`4174` for the same tree with permissive CORS. A second run on the same machine
collides, reports the taken address, and stops rather than waiting on it. Move
it instead:

```bash
CELLUCID_BROWSER_TEST_PORT=4183 npm run test:browser
```

`CELLUCID_BROWSER_TEST_SAMPLE_PORT` defaults to `CELLUCID_BROWSER_TEST_PORT + 1`
and can be set separately. The servers and the specs resolve the address from
the same module, so a moved run also fetches its fixtures from the port it was
moved to, and failure artifacts land in `test-results/<port>/` where concurrent
runs cannot overwrite each other's evidence.

(sibling-repositories-and-environment-variables)=
### Sibling repositories and environment variables

Several contract tests and specs read the other repositories of the workspace:
a shared fixture, a mirrored catalog, the documentation tree, a starting
view. Each finds its sibling beside this checkout
(`../cellucid-python`, `../cellucid-r`, `../cellucid-datasets`,
`../cellucid-annotation`) unless an
environment variable names another place (a worktree, a CI checkout
elsewhere):

| Variable | Default | Read by | When the sibling is absent |
|---|---|---|---|
| `CELLUCID_PYTHON_ROOT` | `../cellucid-python` | `scripts/sync-shared-files.mjs`, the generators and `tests/shared-files-contract.test.mjs` (cellucid-python's copy of every shared file, {ref}`shared-files`), `tests/notebook-widget-requests-contract.test.mjs` (the Python widget copy), the Python viewer specs (`viewer-channel`, `viewer-history`, `viewer-round-trip`, `mcp-local`, `notebook-widget-jupyterlab`: the source they run, named by `PYTHONPATH`), `tests/anndata-copy-cross-language.test.mjs` (the package that reads the web app's copies) | the copies are named and skipped (`t.diagnostic`); a failure under `CELLUCID_REQUIRE_SIBLINGS=1` |
| `CELLUCID_R_ROOT` | `../cellucid-r` | `scripts/sync-shared-files.mjs`, the generators and `tests/shared-files-contract.test.mjs` (cellucid-r's copy of every shared file), `tests/notebook-widget-requests-contract.test.mjs` (the R widget copy), `tests/browser/notebook-widget-jupyterlab.spec.mjs` (the R package source), `scripts/efficiency-floor/served.R` | the copies are named and skipped; a failure under `CELLUCID_REQUIRE_SIBLINGS=1` |
| `CELLUCID_ANNOTATION_ROOT` | `../cellucid-annotation` | `scripts/shared-files.mjs` and `tests/shared-files-contract.test.mjs`: the community-annotation contract corpus, whose home is that repository (`tests/fixtures/contract-corpus.json` there; held here as `tests/fixtures/community-annotation-contract-corpus.json`) | the copy is named and skipped; a failure under `CELLUCID_REQUIRE_SIBLINGS=1` |
| `CELLUCID_REQUIRE_SIBLINGS` | unset | `scripts/shared-files.mjs` (`siblingsRequired()`; `absentSibling(present, what)`, which fails where siblings are required and otherwise says the contract may go without; `skipAbsentSibling(t, present, what)`, which then skips the running test; `pythonDocsRoot()`, the documentation tree), and the packages' test helpers: `1` makes a sibling repository (or what a contract reads of it) that is not there a failure where it would be a skip | — (each repository's CI sets it; in this one the Node job) |
| `CELLUCID_DATASETS_ROOT` | `../cellucid-datasets` | `tests/session-preset-control-inventory-contract.test.mjs`: the starting views the sample catalog advertises (`anndata/datasets.json`), read the way the viewer reads them, held to the control inventory, and `tests/state-record-added-keys.test.mjs` | skipped when the sibling is not checked out, a failure under `CELLUCID_REQUIRE_SIBLINGS=1`; a directory named by the variable that is not a `cellucid-datasets` checkout fails |
| `CELLUCID_PYTHON_DOCS` | `../cellucid-python/docs` | `pythonDocsRoot()` (`scripts/shared-files.mjs`), which `tests/info-tips-contract.test.mjs` (every info tip's link resolves to a heading of these docs), `tests/citation-contract.test.mjs`, `tests/analysis-options-schema-contract.test.mjs`, `tests/action-reference-contract.test.mjs`, `tests/keyboard-shortcut-table-contract.test.mjs`, `tests/category-palette-color-vision-contract.test.mjs`, `tests/seo-head-contract.test.mjs`, `tests/shared-files-contract.test.mjs` and the generators of the action and shortcut reference pages call | the documentation checks skip, saying so; a failure under `CELLUCID_REQUIRE_SIBLINGS=1` |
| `CELLUCID_PYTHON` | none | the Python viewer specs (`viewer-channel`, `viewer-history`, `viewer-round-trip`, `mcp-local`, `notebook-widget-jupyterlab`) through `viewerInterpreter()` (`tests/browser/helpers/viewer-interpreters.mjs`): an interpreter with cellucid's dependencies (JupyterLab and anywidget too for `notebook-widget-jupyterlab`); the package itself comes from `CELLUCID_PYTHON_ROOT`, never an installed copy. `tests/anndata-copy-cross-language.test.mjs` uses it, or the workspace's environment, or `python3`, if it has anndata, mudata and zarr | the spec skips; a failure under `CELLUCID_REQUIRE_VIEWERS=1` |
| `CELLUCID_RSCRIPT` | none | `tests/browser/r-viewer-channel.spec.mjs` and `tests/browser/r-viewer-history.spec.mjs`, through `viewerInterpreter()`: an `Rscript` whose library has the package (installed from the checkout in CI) and its dependencies | the spec skips; a failure under `CELLUCID_REQUIRE_VIEWERS=1` |
| `CELLUCID_REQUIRE_VIEWERS` | unset | `tests/browser/helpers/viewer-interpreters.mjs` and `tests/helpers/anndata-copy-python.mjs`: `1` makes an interpreter that is not named, or a Python that cannot run the copy checks, a failure where it would be a skip | — (only the CI job `viewers` sets it) |
| `CELLUCID_R` | none | `tests/browser/notebook-widget-jupyterlab.spec.mjs`, R part: an R with IRkernel, cellucid-r's dependencies and pkgload | that test skips; a failure under `CELLUCID_REQUIRE_VIEWERS=1`, through `viewerInterpreter()` (the `viewers` job names the `R` it installed) |

In CI the Node job names no place: it checks out, beside this one and at their
default branch, only the files of each sibling its contracts read
(`.github/workflows/test.yml`, sparse: the two widget copies; both packages'
copies of every shared file; this site's Markdown pages; the sample
catalog `anndata/datasets.json` with its starting views; and the annotation
repository's contract corpus), and sets
`CELLUCID_REQUIRE_SIBLINGS=1`.
`tests/ci-sibling-checkouts-contract.test.mjs` holds those lists to the list of
shared files (`scripts/shared-files.mjs`), so a new shared file cannot pass on
a whole checkout and fail on every CI job. What a contract reads of a sibling
and does not find is a failure there (`absentSibling()`), never a skip, so a
path left out of a list is a red job. Only
`tests/anndata-copy-cross-language.test.mjs` reads more (a whole package and a
Python with its dependencies); it runs in the `viewers` job, never in the Node
job.

The `viewers` job (`.github/workflows/test.yml`; Ubuntu, Chromium) runs what
needs the real packages. The viewer specs of Python (`viewer-channel`,
`viewer-history`, `viewer-round-trip`, `mcp-local`,
`notebook-widget-jupyterlab`) and of R (`r-viewer-channel`,
`r-viewer-history`) skip on a machine that names no interpreter, so the browser
shards never run them. The job installs cellucid-python (with the `mcp` and
`mudata` extras) and JupyterLab with its kernel from the checkout, installs
cellucid-r from its checkout with its dependencies, SingleCellExperiment,
IRkernel and pkgload, names the interpreters (`CELLUCID_PYTHON`,
`CELLUCID_RSCRIPT`, and `CELLUCID_R` for the R kernel of the notebook widget
spec, whose kernel the spec registers for its own JupyterLab), points
`CELLUCID_PYTHON_ROOT` and `CELLUCID_R_ROOT` at the checkouts inside its
workspace, and sets `CELLUCID_REQUIRE_VIEWERS=1`: a missing interpreter is a
failure there, not a skip. It then runs the seven specs, the R half of
`notebook-widget-jupyterlab` included, and
`tests/anndata-copy-cross-language.test.mjs`.

(shared-files)=
### Shared files: one home, copies compared byte for byte

What the web app, the Python package and the R package must agree on is one
file with one home and a copy in each repository that reads it. Nothing writes
down a digest of such a file, or a commit of a sibling repository: the copies
are compared, byte for byte, and nobody updates a number.

- **A file the three repositories share has its home in the web repository.**
  `scripts/shared-files.mjs` lists each one, the script that makes it when one
  does, and each copy's place: the action catalog, `mcp-tools.json`, the
  history model, the input formats, the viewer protocol, the category colors,
  the data-mapping, subsampling and latent-outlier conformance fixtures, the
  `gene-sets.gmt`, `obs-columns.csv` and the three golden
  history files (all in `tests/fixtures/`, hyphenated there and underscored in
  the packages; cellucid-r also ships the golden session as its help pages'
  example, `inst/extdata/`), and the nine AnnData files in older layouts
  among `tests/browser/fixtures/` (cellucid-python's
  `tests/fixtures/older_anndata/`).
- **A file only the two packages share has its home in cellucid-python**, and
  cellucid-r holds the copy: the notebook widget's front end
  (`src/cellucid/jupyter/static/widget.js` and `widget.css`; `inst/widget/` in
  cellucid-r), `tests/fixtures/local_agents.json`,
  `tests/fixtures/transition_drift_conformance.json` (written by
  cellucid-python's `scripts/generate_transition_drift_conformance.py`) and
  the MuData cross-writer harness `tests/cross_writer_mudata.py`
  (`tests/testthat/python/` in cellucid-r).
- **A file held here whose home is another repository** is listed beside
  them (`HELD_COPIES` in `scripts/shared-files.mjs`): the community-annotation
  contract corpus has its home in cellucid-annotation
  (`tests/fixtures/contract-corpus.json`, written there by
  `python scripts/generate_contract_corpus.py --copy <the web repository's file>`),
  and the web repository holds the copy
  (`tests/fixtures/community-annotation-contract-corpus.json`).
- **Every side compares.** The web suite compares what a generator makes
  with the home it writes, and every copy with its home
  (`tests/shared-files-contract.test.mjs`). Each package's suite compares its
  copy with the home's file through one helper: `assert_shared_file()` in
  cellucid-python's `tests/_web_root.py`, `expect_shared_file()` in
  cellucid-r's `tests/testthat/helper-repository.R`. cellucid-python, as the
  home of the files the two packages share, compares cellucid-r's copies too.
- **One command makes them agree.** `node scripts/sync-shared-files.mjs`, run
  in the web repository, makes every home that has a generator and writes
  every copy (`--check` writes nothing and names what is stale). The packages
  then write what they build from the copies:
  `python scripts/generate_actions.py` and
  `python scripts/generate_history_model.py` in cellucid-python,
  `Rscript tools/generate-actions.R` and
  `Rscript tools/generate-history-model.R` in cellucid-r. A comparison that
  fails names the command.
- **Where a sibling is.** `CELLUCID_WEB_ROOT`, `CELLUCID_PYTHON_ROOT`,
  `CELLUCID_R_ROOT` or `CELLUCID_ANNOTATION_ROOT`, else the checkout beside
  this one. Where it is not there
  (an installed package, an sdist, a CRAN check) the comparison is skipped and
  says why. Each repository's CI checks its siblings out at `main` and sets
  `CELLUCID_REQUIRE_SIBLINGS=1`, which makes a missing sibling a failure.

A reference fixture an outside tool made (scanpy, edgeR and limma) is the file
in git, with no digest beside it. What does keep a digest does a job with it:
the SHA-256 test vectors, the downloaded compiler and codec sources and the
module they build, the web asset list the packages cache the app by, the
sample datasets' build records, the content hashes a session file is made of,
the conda recipe's source, and the image digest in a documentation
screenshot's record.

(browser-spec-stability)=
### Browser specs that wait for state, not time

A spec that pauses a fixed time for work it cannot see (a redraw, a read, a
notification, a toggle event) passes where that work is quick and fails where
it is not: a hosted runner draws a frame of software GL in a second, macOS
Firefox in seconds. The browser suite keeps `retries: 0`, because a retry
would hide exactly these failures. A spec follows these rules; the failure
behind each is in {ref}`the record <testing-decisions-waits>`.

- **Wait for the product's own state.** A locator's state
  (`toBeVisible`, `toHaveText`, `toHaveAttribute`), `expect.poll` on what the
  page publishes (`getFrameScheduleState()`, `data-startup-phase`, a count the
  state holds), or an event or promise the product exposes. A test that checks
  that something does *not* happen first waits for what ends the work that
  could do it (the opening done, the read settled, a page timer of the
  product's delay armed after the product's own, `outlastPageTimers` in
  `ui-keyboard-lifecycle.spec.mjs`), then asserts. A held request is held on
  a deferred the test releases (`holdOpening` or `holdRequests`, and `releaseAndDrain`, in
  `tests/browser/helpers/dataset-loading.mjs`), not for a set time.
- **No fixed sleeps.** `tests/browser-spec-fixed-waits-contract.test.mjs`
  (`npm test`) finds every `page.waitForTimeout(…)` and every promise a timer
  resolves (`new Promise(resolve => setTimeout(resolve, ms))`, and a timer
  callback that resolves its promise) in `tests/browser/`, page code in
  template literals included. Each must be listed in the contract by file and
  source line. `ALLOWED_FIXED_WAITS` holds only waits whose length is not
  standing in for a state, each with its kind and reason: `time` when time is
  the subject (an observation window in which nothing may be drawn, an
  auto-hide countdown, a press held past a threshold, a simulated slow link),
  `poll` for the interval of a loop that asks the product's state again and
  fails at its deadline, and `deadline` for a bound on waiting for an event
  whose expiry fails the test. `PENDING_RACE_FIXES` is empty and only ever
  removed from: a sleep is never added to it, and its `todo` test reports it
  whenever it is not empty. A timer of no delay (a turn of
  the event loop) and a timer that rejects are not fixed waits.
- **A wait is as long as the machine, and the test's timeout is its fence.**
  A loop that asks the product's state once an animation frame does not count
  its frames, and sets itself no clock.
  `tests/browser-spec-counted-waits-contract.test.mjs`
  (`npm test`) finds every `for` loop that counts to a number, awaits a frame
  or a timer and leaves early, and every loop bounded by a deadline of its
  own; each must be listed with why its count is the subject. The same test
  holds the time a spec gives an assertion or a wait: the suite's own is a
  minute for an assertion and the test's timeout for everything else, and a
  spec that gives one less says how fast the machine must be.
- **No spec reaches the network.** `tests/browser/helpers/outside-requests.mjs`
  watches every context of every spec: a response that names a server address
  outside this machine fails the test, naming it (a route's answer names none).
  It only watches, because a catch-all route would turn off the HTTP cache for
  every spec. The test server serves `index.html` with the suite's own sample
  catalog in place of the published one (`scripts/serve-browser-tests.mjs`), so
  opening the app never depends on the published site, whatever the spec's
  address. A spec plays the hosts it reads from with `page.route`
  (`playAnthropic` for the AI provider).
- **A page is at rest before a spec acts on it.** One helper waits for all the
  work a page starts on its own after an opening, `waitForPageAtRest`
  (`tests/browser/helpers/dataset-ready.mjs`): startup `ready`, the gene list
  read, the frame loop settled (the opening's frame-cost measurement
  included; a scene drawing for a highlight's pulse alone is at rest too) and
  every view on its final geometry. A spec that sets the camera or the GPU
  waits for it first. What that work is:
  - The cell info card asks for a large pane's lookup tree once the page is
    idle, so "no tree before the first hover" holds only while the card is
    off (`cell-info.spec.mjs` turns it off while 300,000 cells open, then
    switches it on and hovers in one task).
  - The gene list is read after the first view, and listing it refreshes the
    interface, which publishes the active view again: a spec that counts
    publications starts once `data-startup-phase` is `ready` and
    `state.whenGenesListed()` has settled.
  - The opening's frame-cost measurement asks for its next sample frame once
    the last one's fence passes: `getFrameScheduleState().settled` is false
    while a sample's fence is awaited or a measurement is armed
    (`frameCostPhase` names the measurement's phase).
- **A change made through the interface is done when the dispatcher says
  so.** A control's change is an action, applied when those before it have
  finished and with what it reads on the way. A spec that goes on through
  the page's internals (a domain restored directly, a pick) waits for
  `getActions().idle()` first (`whenActionsDone`, `dataset-ready.mjs`).
- **A selection gesture is made once, when the page can take it.** A press on
  the plot is read by the tool in use against the frame on screen; the tool's
  mode is switched by an action, and a pick is of the frame that was drawn.
  `tests/browser/helpers/selection-gesture.mjs`
  waits for every action asked so far (`getActions().idle()`) and for the
  page to be at rest, makes the gesture once, and waits for the tool's step
  controls. A gesture is never repeated until its tool arms.
- **A click lands on a page that has stopped moving.** Playwright checks what
  a press lands on and releases a moment later at the same point; if a layout
  change moves the control in between, the release lands on something else,
  the browser sends the `click` to their common ancestor, and the control
  pressed hears nothing, with no error.
  `tests/browser/helpers/split-clicks.mjs` watches every context of every test
  for a trusted mouse press and release at one point whose pressed element is
  no longer under it, and fails the test naming both elements
  (`split-clicks.spec.mjs`). The cure is the spec's: click once the page has
  stopped moving (`waitForPageAtRest`). `openSessionSection` (`panels.mjs`),
  which opens a header under the Session panel's Dataset box whatever a spec
  is about, waits for the interface to be built (`window._cellucidUi`),
  clicks, and clicks again only for a split click reported on its own header.
  A drag has the same hazard one step earlier: wait for what moves the page
  (the steps recorded), move the pointer onto the target itself and check
  what the press landed on before dragging. The product keeps its own
  controls still while it changes a line by itself
  ({doc}`35_decisions_interface`).
- **A box the product takes back is clicked.** A box the page unticks itself
  (a failed velocity allocation) is clicked, not `check()`ed: `check()` fails
  when the page's answer lands before it looks again.
- **A page is closed at a defined point.** A change of the history is handed
  to the browser's storage some tens of milliseconds after it is made, and is
  then kept whole whatever becomes of the page
  ({ref}`all of a write or none <history-autosave>`). A
  spec about what a closed or reloaded page kept waits for that first
  (`window._cellucidHistories.flush()`, or the saved state listed with what
  it holds); a spec about a page that goes at once asserts what the product
  promises then: all of the last write or none of it, never a part
  (`saved-state-close.spec.mjs`, which closes a tab at several moments of one
  write in every engine).
- **A page is parsed only once it is styled.** `index.html` carries a style
  gate (`CELLUCID_STYLE_GATE`): a script straight after the stylesheet link,
  which every engine runs only once the stylesheets before it are in, so
  nothing of `<body>` is parsed and DOMContentLoaded does not come before
  them (WebKit does not hold DOMContentLoaded for a stylesheet; Chromium and
  Firefox do); a stylesheet that did not load is reported on the code-load
  failure screen instead of starting the code without it.
  `style-gate.spec.mjs` holds the stylesheet, and one it imports, for 2.5 s
  in every engine; `tests/startup-code-graph-contract.test.mjs` holds the
  gate's place in the page. `dismissWelcome` waits for the state the product
  keeps (`#welcome-modal[open]`: the welcome is a `<dialog>`, never displayed
  until opened, styled or not).
- **Offline only through `goOffline`.** A page still loading something of its
  own when the connection goes reports that failure too (the analysis tools
  it imports as it starts, a read of the opening still under way).
  `tests/browser/helpers/offline.mjs`
  waits on every page of the context for the page to be at rest and for
  `window._comparisonModule`, which the page sets in the continuation of that
  import, then takes the context offline;
  `tests/browser-spec-offline-contract.test.mjs` allows `setOffline(true)`
  nowhere else.
- **A navigation's leftovers are set aside by the fixture.** Every page of
  every test context, a spec's own `browser.newContext()` included, is
  prepared by `tests/browser/helpers/test.mjs` (`setAsideLeavingReports`,
  `tests/browser/helpers/leaving-document.mjs`): while its `goto`, `reload`,
  `goBack` or `goForward` from a document runs, the wordings WebKit reports
  for the loads that document still had going, and only those, reach none of
  its listeners; no spec handles them itself
  (`tests/browser-leaving-document-contract.test.mjs`). What WebKit reports
  and why: {doc}`38_decisions_testing_and_hosted_runners`.
- **Teardown ends work in flight.** A spec may close a page whose opening is
  still awaited: application teardown retires the live dataset publication,
  whose signal ends that wait and every read of it, and viewer disposal drops
  the deferred WebGL error checks still waiting for their fence.
- **Pixels are compared only after their inputs.** A spec holds, exactly,
  every input to a draw (cameras, order, the bytes of the stores, cells
  drawn) and asserts both frames' draw order before comparing pixels:
  translucent points overlap in draw order, so two frames are comparable only
  when drawn in the same one. A difference attaches both frames, their
  difference and what each was drawn from (camera, level of detail, cells
  drawn, draw order, GL renderer) to the report.
- **What two frames submitted is asserted before their pixels.**
  `tests/browser/helpers/gl-submissions.mjs` records every draw call of a
  frame: program, mode, vertex count, viewport, whether it rasterizes, and
  every uniform of its program at the call; calls that follow one another
  under the same state count as one, so a store drawn one call a segment
  submits what a store drawn in one call does. `segmented-point-stores` and
  `lod-rank-order-resize` assert the two sides equal on every cell of the
  matrix. Specs that count or classify draws tell a point-pass draw by its
  program (one of the renderer's `programs`), not by its count: the GPU
  gather of a tree's draw order (`high-perf/gpu-permute.js`) draws
  `pointCount` points too, with blending off, into its own target; they skip
  transform-feedback passes (`RASTERIZER_DISCARD`), and they count per
  animation frame (`document.timeline.currentTime`), not per pointer move: a
  slow host draws several frames per move.
- **Under a highlight a scene never settles; its frames repeat.** A
  highlight's pulse draws a frame whenever one is due, so `settled` never
  holds. A spec that only acts next waits with `waitForPageAtRest`, for which
  a scene drawn for the pulse alone is at rest. A spec that reads pixels
  waits for a steady pulse frame (`steadyPulseFrame` in `gl-submissions.mjs`):
  nothing invalidated, and the last two frames drawn since the change
  submitted the same, with the pulse's phase pinned
  (`setHighlightAnimationTime`). A spec that waits for the pulse alone
  asks `animationCause` with nothing invalidated, not the cause of the last
  frame drawn (why: {doc}`38_decisions_testing_and_hosted_runners`).
- **A one-shot GL failure** injected into an allocation (`createFramebuffer`,
  `createTransformFeedback`) is aimed at the caller under test by its stack
  (`velocity-multiview-lifecycle.spec.mjs`): the spatial tree's GPU gather
  allocates both at a time the GPU's speed decides.
- **A slow machine on a developer machine.**
  `CELLUCID_BROWSER_TEST_CPU_THROTTLE=6` slows the main thread of every page
  six times (Chromium only, the DevTools protocol's CPU throttling;
  `tests/browser/helpers/test.mjs`), so work is ordered as on a hosted
  runner: what a spec read right after an action is then not there yet. It
  slows a page's own thread, not its workers or the GPU.
- **Repeat before pushing.** `npm run test:browser:repeat` runs specs many
  times in each engine, Chromium, Firefox and WebKit one after another, with
  one worker and no retries as CI runs them, in browser processes no
  longer-lived than CI's (six files and forty tests, the bounded runner's
  limits), and fails on any single failed
  run, printing each engine's passed and failed runs and the tests that
  failed (`scripts/run-browser-specs-repeat.mjs`):

  ```bash
  npm run test:browser:repeat                     # changed specs, 5 times, every engine
  npm run test:browser:repeat -- tests/browser/cell-info.spec.mjs --repeat=10
  npm run test:browser:repeat -- --project=webkit --grep="hover"
  ```

  Without paths it takes the specs changed relative to `origin/main`
  (committed, uncommitted and untracked) and every spec that imports a
  changed helper of `tests/browser/`, directly or through another helper; a
  helper path given explicitly stands for the specs that import it. Each
  engine takes a free port pair from the bounded runner's planner, so it runs
  beside another browser-test run. A machine that is not a hosted runner does
  not show every runner's race, so a green repeat is evidence, not proof: a
  spec that failed in CI is fixed by finding what it raced, then repeated.

(fixed-device)=
### The device every browser spec runs on

The page measures its device and acts on what it measures: level of detail
and frustum culling turn on when the GPU's share of a frame of every cell is
measured over the interactive budget, and stay on only where their frames are
measured faster (`rendering/full-detail-frame-cost.js`, `render-controls.js`);
an opening is progressive when showing its cells is predicted slow from what
showing a cell cost earlier openings (`data/opening-cost.js`); worker pools
follow the cores the browser reports, memory budgets the memory it reports.
A spec about none of those choices must not see the page change state under
it at a time the runner decides, so **every spec runs on a fixed device**,
the `device` option of the `test` fixture
(`tests/browser/helpers/device.mjs`, installed on every context the test
opens, `browser.newContext` included). Why, and what it replaced, is in
{ref}`the record <testing-decisions-device>`.

- `gpu`: the GPU's time for a frame, given exactly through the GPU timer
  the page reads where a browser offers one
  (`EXT_disjoint_timer_query_webgl2`, simulated in every engine): 1 ms by
  default, so the page never turns level of detail or culling on by
  itself. A spec about that choice sets its own (`{frameMs,
  perThousandCellsMs, culledMs}`: flat, following the cells the live view
  drew, or cheaper while the viewer culls), at start
  (`test.use({ device: { gpu: { frameMs: 500 } } })`) or later
  (`setSimulatedGpu(page, …)`). Times far above any frame's own script are
  exact on any runner, since a frame lasts the longer of its GPU's and its
  page's share; the timer's answers are ready at the first read, so a slow
  simulated GPU costs a spec nothing. A spec whose subject is the real
  GPU's time says `gpu: 'real'` (`snapshot-split-vbo-rendering`).
  `full-detail-frame-cost.spec.mjs` asserts each verdict exactly, number for
  number.
- `hardwareConcurrency` (4) and, where an engine reports it (Chromium),
  `deviceMemory` (8), so the pools and budgets are sized alike everywhere:
  the page's `navigator` reports them, and the workers, whose own
  `navigator` an init script cannot reach, are handed the cores (the
  machine's `cores`, next) and take the page's memory budgets
  (`data/memory-budget.js` `connectMemoryWorker`).
- The machine's measures (`assets/js/utils/machine.js`): **the one place
  a rate measured of the machine or the link enters a plan**. The fixture
  holds them in the page's module (`holdMachine`), and the page hands what
  it holds to every worker it starts, in the fragment of the worker's
  address (`workers/worker-script.js` `handToWorkers`: the same script is
  asked for, a fragment being no part of a request, and the worker keeps
  its name); the reader worker hands them on to its decode workers. A worker's module begins with them,
  before any of its code plans. Each is a value, or `'measured'` in a spec
  whose subject is that measurement:
  - `link`: the link the reader plans on, `{latencyMs, bytesPerMs}`: 2 ms
    and 32 MB/s by default (a host on a fast local network; below the rate
    a web host is read ahead at, so fixtures are read on demand, as a link
    opens). A host on this machine is priced as this machine at rest on
    top of it, as in the app. The requests themselves are still timed and
    shaped by the source's own meter: how long a stalled request is waited
    for must follow the real link, on a slow runner above all. A spec
    about a slow link gives its link here (`{latencyMs: 200, bytesPerMs:
    100}`) instead of delaying answers.
  - `decodeBytesPerMs` (1 GB/s), `orderMadeMs` (0),
    `keptReadBytesPerMs` and `indexPassMs` (null: none measured).
  - `openingMsPerCell`: what showing a cell costs (the machine's
    `shownCellMs`): nothing measured by default, so every opening is
    priced by its reads alone, and the page's own records are not made.
    `browser.newContext({ device: { openingMsPerCell: 2 } })` makes an
    opening progressive; `'measured'` lets the page record (the record's
    own spec in `open-cells-in-place`).

`held-machine.spec.mjs` proves the holds in every engine; what it asserts and
the runs that made the link and the decoding rate held:
{doc}`38_decisions_testing_and_hosted_runners`.

A spec that sets level of detail or culling through the viewer
(`setAdaptiveLOD`, `setFrustumCulling`) still makes that opening's choice,
and the render controls leave the rest of its reports alone.

WebKit only: a spec of `persistentStorage` runs with one worker. WebKit
keeps a persistent context's origin storage outside its profile folder
(`helpers/test.mjs` `clearOriginStorage`), and two workers running such
specs at once lose what they keep (CI runs one worker).

The holds reach a page that runs its sources (this machine). A page
opened as a deployed site runs its bundle, whose module instances an
init script cannot name, and measures everything: a spec opened there
asserts nothing that follows from a plan.

(measured-rates)=
**Every rate the app measures, and what a test does with it.** A new
measured rate that feeds a plan gets a name in `utils/machine.js`
(`tests/machine-seam-contract.test.mjs` holds the list to the sources);
one that only paces work is added to this table with why a slow runner is
safe with it.

| Measured | Decides | In a test | Safe on a slow runner because |
| --- | --- | --- | --- |
| The link's round trip, rate, round trip under load, cost of a request (`formats/io/link.js`), as the reader plans on it | The gene sample, the block rung of the cell order, a progressive opening, values read on demand or ahead, a kept index against building it again | Held: `link` (browser: the device; Node: `LINK_AT_REST` where the host has no time of its own, else the file's own clock) | held |
| What earlier visits and the page's own downloads measured of the host (`data/host-link-memory.js`) | The prices above until the opening measures more | Not weighed while `link` is held (`plannedHostLink`) | held |
| The same link as requests are shaped by it (merge gaps, request size, the concurrency search, the protocol guess) | How many requests carry a read, and how long each is | Measured | the bytes read are the same: no spec asserts a count of requests, and a Node test of requests gives its byte reader a clock |
| Response times (`ResponseClock`) | When a stalled request is given up and asked again | Measured | it must follow the real link: its patience grows with what the runner measures |
| A host's pace (`formats/io/hosts.js`) | The spacing of requests to a host that refused some | Measured | learned only from refusals (429, 403), which a fixture gives only where that is its subject |
| How fast chunks decode (`formats/io/decode-pool.js`) | The gene sample, a progressive first view | Held: `decodeBytesPerMs` | held |
| Which decoder is fastest (`formats/io/decoder-choice.js`) | The decoder of a codec family | Measured | every decoder gives the same bytes |
| How long making every cell's order took | File order or a `blocks` order for every cell | Held: `orderMadeMs` | held |
| How fast kept gene indexes were read | A kept index read, or built again | Held: `keptReadBytesPerMs` | held |
| What a gene index pass's work cost an entry (`PassRates`) | Which pass indexes a share spread through the store | Held: `indexPassMs` | held |
| What showing a cell cost (`data/opening-cost.js`) | A progressive opening | Held: `openingMsPerCell` | held |
| How long the first gene and each top-up read took (`field/gene-sample.js`) | How many more cells genes are read for after the first gene | Measured | it runs only on a gene sample fitted to the gene time, which the held link and decoding never make; a spec that gives a slow link asserts what holds for any top-up |
| What an analysis kind cost (`analysis/shared/progressive-cells.js`) | Whether a later run of it shows a sample pass first | Measured | the final result is the same cells' and is what specs wait for |
| Analysis read and task times (`run-tuner.js`) | How far ahead genes are read, how many tasks are in flight | Measured | scheduling only |
| The GPU's time for a frame (`rendering/full-detail-frame-cost.js`) | Level of detail and culling by itself | Held: `gpu` (a simulated timer) | held |
| What frames held the main thread, how late the page's frames ran behind the GPU, and the display interval (`rendering/ambient-frame-pacing.js`) | Which ambient frames (pulse, smoke, velocity particles, backdrop motion) are held back | Measured | specs wait for what was drawn (an overlay draw, the page at rest), never for a frame count |
| The GPU's time for a level of a backdrop's sky (`rendering/backdrops/backdrop-look-renderer.js`) | How far the sky is refined (its faces' resolution) | Measured | not verified for every spec (a slow GPU refines the sky less; the record's Known limits say what a spec may assert) |
| What a neighbor row cost from the file, and a pass over the held edges (`app/ui/modules/connectivity-controls.js` `knnRowReader`) | Which of the two gives a cell's neighbors | Measured | both give the same rows |
| The cores | Pool sizes | Held: `hardwareConcurrency` (page `navigator`, the machine's `cores` in workers) | held |
| The memory the browser reports | Budgets, the cells an auto-fit opens | Held where reported (`deviceMemory`, Chromium); Firefox and WebKit report none and take 8 GiB | determined |
| An allocation the engine refused; a tab that crashed while opening | The memory allowance halved | Measured | only where the runner really cannot allocate, which the page must then follow |
| Free storage (`storage.estimate()`) | A gene index on disk, a copy kept or refused | Measured | the fixtures are megabytes; a runner without room fails by a named storage error |
| Elapsed time in loops that yield (index passes, selection queries, noise generation, inflate) | When a loop gives the thread back | Measured | scheduling only: the slices are bounded by time, so a slow runner does less a slice, not longer slices |

(long-suites)=
### Long, sharded and opt-in suites

- **Every browser spec file declares tests.** Playwright says nothing about a
  spec file that holds none; `tests/browser-specs-hold-tests-contract.test.mjs`
  (`npm test`) fails on a spec file under `tests/browser/` that calls no
  `test(…)` / `test.describe(…)` (or hands `test` to a helper that declares
  them, `declareWalkTests`) or does not end with a newline.
- **The control walk** (`tests/browser/helpers/control-walk.mjs`) walks every
  control the page shows (discovered from its interactive roles, never from a
  list; what a control opens, makes or enables is walked too; every setting
  `settings.list` names that is shown must be met) on a fixture dataset in
  five scenes. It fails on a control that changes state with no
  step, an outside change, more than one step after the step shown, a step
  without a label, or an undo that does not put back exactly (the state
  fingerprint; {ref}`drift guards <drift-guards>`); new steps on
  another branch none of which follows the step shown (Clear history keeping
  the state, Branch without this step replaying the later steps) are judged
  as an edit of the history.
  It is one test per sidebar section, one for everything outside them and the
  graph-and-vectors test, spread over `control-walk-1-of-4.spec.mjs` …
  `4-of-4`, as many tests in each (Playwright's `--shard` splits by test
  count) and filled by weight, so the longest tests are in different files (a
  test weighs the controls its walk operated when last recorded:
  `tests/browser/helpers/control-walk-weights.json`, written by a run of the
  walk with `CELLUCID_CONTROL_WALK_WEIGHTS=record`, never by hand; a part the
  record does not name weighs as the average one, and every part is in
  exactly one file whatever the record says,
  `tests/browser-workflow-sharding-contract.test.mjs`), and runs on Chromium only (headless Chromium has no volume
  renderer, so its smoke settings are never shown). In CI it runs in a job
  of its own, a file a job (`control-walk` in `.github/workflows/test.yml`:
  `CELLUCID_CONTROL_WALK=only`, `--shard=k/4`), and the eight shards of the
  browser job leave it out (`CELLUCID_CONTROL_WALK=apart`): Playwright splits
  shards by test count, so four walk files of minutes a test would all fall
  into one shard. Run one part with
  `npx playwright test tests/browser/control-walk- --project=chromium --grep "<section id>"`,
  or the four files side by side with `--shard=<i>/4`; each part prints
  where its time went, by section, its slowest controls and phase (about 15
  minutes for the whole walk on a GPU). A control is judged by what it does,
  never by how long it takes. What bounds the walk is its progress: a wait
  of the walk that has not ended after ten times the slowest wait the part
  has seen end (and at least two minutes) is a hang, which names the control
  and ends the part; the part's test timeout moves with each wait that ends;
  a Playwright action waits 20 of the page's frames as the last settle
  measured them (at least 5 s). After Escape it waits until the page is
  still, never a fixed time.
- **History stress** (`tests/stress/`, `npm run stress:history`, which runs
  `scripts/history-stress.mjs`: `--steps 1000,10000,50000,100000`,
  `--cells`, `--quick`, `--readers` for the Python and R session readers):
  `tests/history-stress-harness.test.mjs` is the short check in `npm test`;
  the long runs and the opt-in browser spec
  (`CELLUCID_HISTORY_STRESS=1 npx playwright test tests/browser/history-stress.spec.mjs --project=chromium`,
  a 1M-cell page, run alone) measure against the budgets and results in
  `tests/stress/RESULTS.md`.
- **Live agents** (`npm run test:live-agents`): the specs in
  `tests/browser/live-agents/` against the real Anthropic and OpenAI APIs,
  Chromium, one worker. They need the developer's keys (the environment's
  `ANTHROPIC_API_KEY` / `OPENAI_API_KEY`, else the file `CELLUCID_AGENT_ENV`
  names, else the checkout's `.env`, which is git-ignored and never printed)
  and skip with the reason when a key is missing; they never run in CI, since
  they spend real credit ({doc}`32_agents_in_the_tab`).

(hosted-runners)=
### What the hosted CI runners differ in

A spec that passes on a developer machine can fail only in CI, for reasons no
change of the spec's subject explains. What a hosted runner is, as far as a
spec can tell; the runs that showed each, with their numbers, are in
{ref}`the record <testing-decisions-runners>`.

- **Software GL.** Hosted Linux runners draw on a software rasterizer:
  Chromium on SwiftShader, Firefox on llvmpipe, WebKit behind the masked name
  "Apple GPU" at llvmpipe's speed. The shard
  running `webgl2-runtime.spec.mjs` logs the renderer of each cell. macOS
  WebKit draws on the virtual machine's GPU; macOS Firefox's GL is about 2,500
  times slower than a developer Mac's, and its cost is per pixel, not per
  cell.
- **One cell of the matrix draws less, for time.** Firefox on GitHub's
  macOS runner takes a minute for a frame of full size, so a spec there draws
  less where it says so (`tests/browser/helpers/hosted-gl.mjs`:
  `segmented-point-stores` draws 12,000 points, still three segments), and a
  test whose budget follows the GPU measures it (`smoke-dataset-switch`,
  `renderer-smoke`, `segmented-point-stores`). What it draws, it compares
  exactly like every other cell: `tests/browser/helpers/pixel-parity.mjs`
  allows no tolerance, and the spec holds every input to a draw and what the
  frames submitted to the GPU equal first. That cell with an exact
  comparison is not established by a run on any machine here.
- **Core counts.** The analysis worker pool's size follows the cores, so a
  hosted runner may have one worker;
  `analysis-density-workers-identity.spec.mjs` requires the slice count the
  pool's size plans, and the bit-for-bit identity at any count.
- **Browser-reported metrics.** A browser records no largest contentful
  paint if the first input (the Escape that closes the welcome dialog) came
  before a text or image paint, so `analytics-production.spec.mjs` requires
  LCP only when the browser recorded one. A routed response the page let go
  of by navigating away is disposed, and fulfilling it fails: the route
  handler lets that one error pass.
- **Sibling checkouts.** The Node job checks out `cellucid-python`,
  `cellucid-r`, `cellucid-datasets` and `cellucid-annotation` at their default
  branch when it starts.
  A cross-repository change whose sibling push lands later fails until the next
  run; `notebook-widget-requests-contract.test.mjs` says which commit each
  checkout is at and runs its behavior contracts on both copies, so the
  failing copy is the one behind. Push the siblings first.
- **No WebKit on Windows.** Playwright's Windows WebKit is the WinCairo port,
  not Safari: it has no origin-private file system and no `OffscreenCanvas`.
  Safari is covered by WebKit on macOS (every push) and Linux.

**Software GL on a developer machine.** Playwright's headless Chromium draws
with SwiftShader on a Mac too (`ANGLE (Google, Vulkan 1.3.0 (SwiftShader
Device …))`) unless launched with `--use-angle=metal`, and it composites the
page in software: the configuration Chromium chooses itself when SwiftShader
is its only GL. The plain `chromium` project is therefore the local
reproduction of software GL. It is slower than a GPU but faster than a hosted
runner, so a race that needs a runner's speed may still not show.

Do not force SwiftShader onto the whole GPU process instead (`--use-gl=angle
--use-angle=swiftshader --enable-unsafe-swiftshader`): that is a machine no
runner and no reader has, and a spec that fails only under it measures that
configuration, not the app (what it does is in the record, "Tried and
reverted"). To run the automatic level of detail and culling in every engine
as a runner does, run the specs on a copy of the tree with
`INTERACTIVE_FRAME_BUDGET_MS` set to 0 in
`rendering/full-detail-frame-cost.js` (every frame of every cell is then too
slow), never on the working tree itself.

---

## Manual smoke test checklist (recommended)

Run these steps on at least one small dataset (fast) and optionally one larger dataset (stress).

### A) Boot + basic UI

1) Open the app (local or hosted) and confirm:
   - no console errors at startup
   - WebGL renders a point cloud after dataset load

2) Toggle sidebar open/closed; resize sidebar.

3) Switch theme (light/dark) and confirm plots and UI update.

### B) Data + fields

1) Load a dataset via at least one relevant path (local-demo / local-user / remote / link), in the format your change touches.
2) Switch between:
   - a categorical field
   - a continuous field
   - a gene expression field (if present)
3) Confirm legend updates and color buffers update (no stale colors).

### C) Filtering + visibility

1) Apply at least one filter (category toggle or continuous range).
2) Confirm:
   - visible point count changes
   - points actually hide/show (alpha mask updates)
3) Toggle outlier filtering (if available) and confirm it affects visibility.

### D) Highlights + pages

1) Create a highlight group; add/remove cells.
2) Create a second highlight page; switch between pages.
3) Confirm highlight overlay and counts update correctly.

### E) Multiview + embeddings

1) Click “Keep view” to create a kept view.
2) Switch view layout (grid/single) if available.
3) In `View 2 shows`, choose another embedding; press 1, 2 and 3 on each view
   (the default, then the next embedding of that dimension).
4) Confirm **Link cameras** links only views that show the same embedding.

### F) History and sessions

1) Make a few changes through the controls (a field, a filter, a highlight,
   the camera); confirm each is one step in the State panel's branch graph,
   and that undo (Ctrl/⌘+Z) and redo put each back.
2) Undo two steps and make another change; confirm a branch is started and
   both branches can be shown (double-click a row; a single click only
   selects it).
3) Write a session file (State › Session › **Share…**), reload the page and
   load it (**Load…** › **Open a session file…**); confirm its tracker is
   added beside the one there was.
4) Confirm:
   - the step shown is the one saved, and every step and tracker is back
   - camera/layout, active fields, filters and highlight pages/groups are back
   - a jump to an earlier step works right after the load
5) Close and reopen the tab on the same dataset; confirm the Continue offer
   puts the history back.

### G) Analysis

1) Open the **Analysis** panel and run at least one plot.
2) Switch highlight pages and confirm analysis UI reacts.
3) (Optional) Open a copied analysis window and close it; confirm resources clean up.

### H) Figure export

1) Export a PNG and a small SVG.
2) For larger datasets, confirm large-export prompts work and hybrid/optimized exports complete.

### I) Community annotation (only if touched)

1) Sign in and confirm identity is shown.
2) Connect to a repo where the GitHub App is installed.
3) Pull and confirm cached files update.
4) Publish a small change (or verify the PR flow when lacking push access).

---

## Automated checks

### Protected smoke-density contracts

The browser matrix runs the density transaction in Chromium, Firefox, and
WebKit. Its retained checks prove:

- byte equality of the density channel (red of the RGBA8 volume, within one
  8-bit step) against the historical R8 GPU pipeline, including adversarial
  Float32 inputs;
- zero product `readPixels` calls;
- a maximum 262,144-cell staging upload (one position buffer and one RGBA8
  color buffer) and bounded draw submissions;
- exact GL-state restoration, including samplers and separate blend equations;
- complete resource cleanup after an injected mid-publication failure;
- atomic renderer publication, prior-volume rollback, and no owned texture in
  the public return value; and
- context-loss invalidation plus the viewer's reload-required contract.

`smoke-volume-reconstruction.spec.mjs` checks what the ray march does with that
volume: the smoke is drawn at the world positions the splat wrote, one density
column reconstructs as a single monotone tent, the silhouette stays inside the
voxels written, and empty-space skipping through the occupancy volume loses no
smoke. In the Node suite, `smoke-control-scales.test.mjs` pins every smoke
slider's shipped readout to its fixed scale, and `smoke-transaction.test.mjs`
pins the per-frame ray budget, the 0.25×–1.00× resolution range, and that a
still smoke view asks for no frames once its noise exists.

Keep these checks cross-engine. Do not replace the historical GPU authority
with a JavaScript voxel model: binary64 arithmetic does not reproduce GLSL
Float32 cutoff behavior.

(drift-guards)=
### Drift guards for actions, state domains and the history

The guards that fail when the layers of
{doc}`29_action_layer_domains_and_history` drift apart:

| Guard | Fails when |
|---|---|
| `defineAction`, `assertActionSchema` (at import) | A definition has an unknown field, a schema outside the subset, an open object, an array without `items`, an unknown effect, domain or `confirm` kind, a read or history move that touches a domain, a change that touches none, or no `describe` or `code` |
| `tests/action-catalog-contract.test.mjs` | An action's own example is invalid, its label or pseudocode is not one trimmed line, it names an unknown domain, an action that records no step (`read`, `control`, `history`, `dataset`, `outside`) touches a domain or has no `result`, a change touches none; or `history.undo`, `history.redo`, `history.jump` stop being actions of effect `history` |
| `tests/action-catalog-document-contract.test.mjs` | `tests/fixtures/action-catalog.json` is not what the definitions make; the document is not canonical; a tool name or parameter is not a valid Python method, R function or keyword |
| `tests/shared-files-contract.test.mjs` | A package's copy of a shared file differs from its home in the web repository (compared when the sibling is beside the checkout; a missing sibling fails under `CELLUCID_REQUIRE_SIBLINGS=1`) |
| cellucid-python `tests/test_action_catalog.py` | The packaged catalog is not byte for byte the web repository's `tests/fixtures/action-catalog.json`, or not what the definitions make now; `_actions_generated.py` differs from what the generator writes; a method does not take the action's parameters, or hides another attribute |
| cellucid-r `tests/testthat/test-actions.R` | The catalog is not byte for byte the web repository's, or not what the definitions make now; a generated function, help page or NAMESPACE block differs from what the generator writes, or a page is left over; a generated name takes a hand-written export's or a reserved word |
| `tests/relay-bundle-contract.test.mjs` | The relay bundle is not what the web sources make now (an action added, a parameter changed, an instruction edited) |
| `tests/agent-bridge.test.mjs` | The MCP tools are not the catalog's, or the agent instructions name a tool that does not exist |
| `tests/viewer-protocol-contract.test.mjs` | The protocol's copy of the action type pattern is not the catalog's |
| `tests/action-history-core.test.mjs` | The dispatcher and history semantics on fake domains: one step per change with its origin, no step for a no-op or read, coalescing, rollback of a failed handler, refusals, agent confirmation, branches, jumps restoring only what differs, delete a step with a flagged dependant, outside changes as steps, observed changes, re-capture after `cells`, history moves as actions, turns, ordering, trackers saved and loaded; and the value hash of NaN, ±Infinity and -0 |
| `tests/state-domains-contract.test.mjs` | A domain of `domains.js` has no owner module in `OWNERS` (or two); a value shape, a cell by original position, a field by name, equality and the pool hash, an exact restore, a describe line or the fingerprint regresses; `controlDomainOf` names a control that does not exist or splits a control between two domains |
| `tests/browser/state-domains-roundtrip.spec.mjs` | A domain changed through the interface does not come back exactly (value, visible state, fingerprint), for the live and a kept view, after a data-mapping change, and on a subsample |
| `applyUndoRedo` in `tests/browser/actions-fields.spec.mjs`, `actions-camera-path.spec.mjs`, `actions-highlights-views.spec.mjs`, `actions-settings-layout-analysis.spec.mjs` | An action changes a domain it does not list in `touches`, or undo or redo does not give the state fingerprint back |
| `tests/browser/actions-*.spec.mjs` ("every control … dispatches the action it stands for") | A control of the area stops dispatching its action with exact parameters |
| `tests/history-autosave.test.mjs`, `tests/browser/history-autosave.spec.mjs` | The autosave writes more or less than what changed, loses a step on a closed page, or stops pausing on a full storage |
| `tests/session-history-file.test.mjs`, `tests/history-golden-fixture.test.mjs` (with cellucid-python `tests/test_sessions.py` and cellucid-r `tests/testthat/test-sessions.R`) | The session file stops round-tripping a history, or the shared golden history no longer reads as `history-golden.expected.json` says, is no longer what the generator makes, or a package's copy differs from the web repository's ({doc}`10_sessions_persistence_and_serialization`) |
| `tests/history-exports.test.mjs` | The pseudocode, methods paragraph or report drops a step, an outside change, a squashed step (any cause) or a note, loses a group heading, a bookmark or the abstract, or cannot export a branch that is not shown |
| `tests/action-method-classification-contract.test.mjs` | A public `DataState` method (its mixins in `app/state/managers/`) or a key of the viewer's public API is not classified in `app/actions/classification.js` as `reads()`, `internal(reason)` or `mutates(domains, actions)`; a classification names a method that is gone, or an owning action the catalog lacks or that does not touch the method's domains |
| `tests/action-mutation-owners-contract.test.mjs` (the static mutation scan) | A `mutates` method is called outside the modules `MUTATION_OWNERS` names (an action's implementation, a domain owner, an observation recorder, a dataset's publication), naming the file, the line and the action to dispatch instead; an owner entry no longer holds; an underscore member of the state is called outside `app/state/` (`NAME_COLLISIONS` lists same-named methods of other objects) |
| The control walk (`tests/browser/helpers/control-walk.mjs`; `tests/browser/control-walk-1-of-4.spec.mjs` … `4-of-4`) | A control the page shows changes state with no step, keeps a "Changed outside the history" step, makes more than one step after the step shown or a step without a label, or its step's undo does not put back exactly (the state fingerprint). What it walks and how it runs: {ref}`long-suites` |
| `tests/browser/history-keyboard-chrome.spec.mjs` | A global key shortcut, a movement key, a key no shortcut uses or a chord changes state with more or fewer than one step, or keeps an outside change; a keyboard camera gesture is not one step; the phone sheet changes the history |
| `tests/browser/action-catalog-coverage.spec.mjs` | An action of the catalog has no implementation on a page with a dataset open (`unprovided()`), a domain has no owner (`missing()`), or a method of the live state or viewer is not classified |
| `tests/action-reads-contract.test.mjs`, `tests/browser/action-reads.spec.mjs` | A `reads` or `context` declaration is unknown, duplicated, a touched domain, both, non-empty for a stepless action, or leaves a parameter's name unresolved; an action applied in two states that differ only in domains it neither touches, reads nor resolves names in gives different after-values |
| `tests/action-versions-contract.test.mjs` | An action's parameter shape changed without a version bump and an upgrader (rewording does not count), an older example (`tests/fixtures/action-versions/`) no longer upgrades to valid parameters, or a version of a recorded action has no example file (a file is written once: the generator never writes over one that exists) |
| `tests/history-catalog-contract.test.mjs`, `tests/history-actions-contract.test.mjs` | A `history.*` action and its provider disagree (a result outside its schema, an edit not journaled, a refusal without its sentence) |
| `tests/history-model-contract.test.mjs`, `tests/history-model-schema.test.mjs` | A record's validator accepts what the model's rules refuse, or `history_model.schema.json` is not what `historyModelDocument()` makes, here or in a sibling |
| `tests/history-v2-core.test.mjs`, `tests/history-tree-properties-*.test.mjs` | The core's semantics regress; 2,100 seeded random operation sequences (three shards) break an invariant after some operation: older restores, interleaved turns, two writers at once, loads beside, each step's before-values its parent's, no loose value in the pool |
| `tests/action-reference-contract.test.mjs` | {doc}`30_action_reference` differs from what `scripts/generate-action-reference.mjs` writes from the catalog |
| `tests/agent-data-guard.test.mjs` | An action that answers with data is not classified in `app/agent/bridge/data-guard.js`, or the guard lets a cell value, a picture or a single-cell statistic through under "names and summaries only" ({doc}`32_agents_in_the_tab`) |

No control calls a mutating state method itself: the static mutation scan
holds it at the source, and the control walk and the area specs drive every
control and compare the step it makes.

### CSS design system validation (recommended when touching CSS)

Docs:
- `cellucid/assets/css/README.md`

Scripts:

```bash
node cellucid/scripts/validate-token-types.js
node cellucid/scripts/validate-tokens.js
```

These catch:
- token/type drift
- accidental raw hex usage outside token files
- theme mapping contract issues

---

## Documentation build sanity check

If you changed docs pages, build the `cellucid-python` docs locally, as CI
and Read the Docs do, with warnings as errors:

```bash
sphinx-build -W --keep-going -b html docs docs/_build/html
```

This ensures every `{doc}` and `{ref}` role and every toctree is valid.

---

## Release checklist (web app)

`www.cellucid.com` serves `main`; a deployment is a merge to `main`. A GitHub
Release `vX.Y.Z` in `theislab/cellucid` marks the version and publishes the
npm name-reservation package (`npm-publish.yml`); the order with the Python
and R packages is in {doc}`/contributing/releasing`.

Because Cellucid is a static web app, releases often fail due to “small” deployment details.
Before deploying:

1) Run the manual smoke tests on a clean browser profile.
2) Confirm no console errors on load.
3) Confirm dataset loading works in the target deployment environment (CORS, mixed content, CSP).
4) Confirm the dataset host answers range requests with `206` and does not re-encode AnnData files (see {doc}`03_build_run_and_deployment`).
5) Confirm session save/load and export still work.
6) If community annotation is enabled, confirm worker origin + CORS allowlist are correct (see `cellucid/assets/js/app/community-annotations/REPO_SETUP.md`).

---

Next: {doc}`15_extension_points_overview`.
