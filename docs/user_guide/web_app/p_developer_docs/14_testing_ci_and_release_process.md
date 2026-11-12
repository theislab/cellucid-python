# Testing, CI, and release process

This page defines how to **validate changes** to the Cellucid web app before you share them with others.

Cellucid is a high-performance interactive app; many regressions are not caught by “does it run?” checks.
The goal here is to make validation repeatable.

## At a glance

**Audience**
- Contributors: follow the manual checklist before opening PRs or sharing builds.
- Maintainers: use the release checklist to avoid “it worked locally” deployments.

---

## Validation contract

- The web app is shipped as static files (no required bundler build).
- `npm test` runs the complete Node contract suite.
- `npm run test:browser` runs the Playwright browser suite.
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
- UI and state changes also require the manual smoke checks below.

In CI, `npm test` runs on Ubuntu, macOS, and Windows across Node 22, 24, and
26 (Node.js 22.15 or newer is required). Every push and pull request runs the browser suite in Chromium,
Firefox, and WebKit on Ubuntu, plus WebKit on macOS; every engine on every
operating system runs nightly and on manual dispatch. A failed browser job
uploads its Playwright screenshots and error context as a workflow artifact.
Each job runs one shard through `scripts/run-browser-shard-bounded.mjs`,
which refuses a shard whose files need more browser processes than its
runner allows. A test belongs to the spec file whose top-level suite holds
it, not to the file its `test()` call is written in (the control walk's four
files declare theirs through `helpers/control-walk.mjs`; grouped by that
helper, the runner once ran none of their 15 tests and reported success),
and before each batch runs the runner lists exactly what it will pass to
Playwright and refuses a selection other than the plan.
`tests/browser-shard-capacity-contract.test.mjs` plans every
shard of the matrix for every engine and runner from one Playwright listing
(`tests/helpers/shard-plan-reporter.mjs` replicates Playwright 1.65.0's
`createTestGroups` and `filterForShard`, a group weighed by its tests not
expected to be skipped, refuses another version, and is
compared test for test with real `--shard` listings of every other shard
and the last, which hold every boundary). It listed each of the 21 shards
in a process of its own before: 26 s idle and 166 s under load, past the
suite's former bound of 120 s a file (the contract suite now bounds progress
instead: a file stops only after finishing no test for 5 minutes,
`scripts/run-web-tests.mjs`); it now takes 5 processes, 7 to 15 s idle and 72 s at a
load average of 40 to 54.

### Contract tests that hold on every Node version and operating system

What the contract suite once got wrong only on some CI cells (fixed
2026-11-01):

- **A process started in a group of its own** (`detached`, so a hung test
  file can be stopped with everything it started) no longer hears the
  terminal's Ctrl-C or a canceled job's SIGTERM: whoever started it forwards
  them (`scripts/run-web-tests.mjs` stops every file it started, 2026-11-07).
- **Starting a process is measured, never assumed**: on a runner five times
  slower than its peers, starting Node and the network guard took longer than
  a margin chosen on a fast machine; a test of a time bound times a short run
  first and scales from it (`tests/web-test-runner-contract.test.mjs`).
- **A rate is judged only on a watchdog tick that came on time**: when the
  whole process paused past a measuring window, the timer ran before the
  bytes that arrived were read, and a healthy transfer looked empty
  (`scripts/install-linux-engine-libraries.mjs`).
- **A file path from a module URL** is
  `fileURLToPath(new URL('…', import.meta.url))`, never `.pathname`: on
  Windows that is `/D:/a/…`, which Node opens as `D:\D:\a\…`.
- **The repository's `package.json` declares `"type": "module"`**, so every
  `.js` in it is an ES module (the token validators in `scripts/` included).
  Without it Node warned about each file it had to reparse, through
  `console.error`, and on Node 22 the warning landed inside a test that
  counts `console.error` calls.
- **An exact power of ten** is parsed from its literal, ``Number(`1e${n}`)``:
  `Math.pow(10, -4)` is `9.999999999999999e-5` on Node 22, which made
  figure export's tick step for 0.001–0.003 jump from 5e-4 to 1e-3.
- **A test waits for the event it needs**, never a fixed pause. Windows
  timers fire on a 15.6 ms tick and its message delivery is slower, so a
  pause that covered a delivery on Linux and macOS did not there:
  `BroadcastChannel` delivery outlasted the 20 ms the memory census test
  slept, and a range the page cache was handed had not arrived 10 ms later
  in `formats-io-contract` (run 41292953076, Node 22 and 24). Every wait in
  the suite is one of four kinds:
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

  A link the test gives a reader is measured as the test means it, round
  trip included: the page cache is given one round trip of the link to
  answer, so a held range "read around" on a 1 ms link was fetched again on
  a runner whose page answered slower than 1 ms; the test now prices bytes
  dear on a 1 s link instead. To check the suite against slow delivery and
  coarse timers before CI does, run it under
  `NODE_OPTIONS="--import ./scripts/slow-delivery-preload.mjs" npm test`:
  every port and channel message arrives up to 25 ms late, in order, and
  every timer runs on a 16 ms tick, in order. At 2026-11-01's HEAD it failed
  the range cache, notebook widget, served-store and reader-channel tests
  (a fixed pause or a bounded count of `setImmediate` turns waiting for a
  port's message) and left `viewer-embed-contract` hanging on a port its
  failed test never closed. Once those were fixed it left
  `gene-read-every-source` running forever after all its tests passed: a
  read the reader started in the background (the row pointers of a Zarr
  folder over HTTP, for a gene index it was keeping on disk) carried no
  signal and outlived the source's disconnect, retrying its host. Closing
  a store now stops every read made through it (`close-stops-reads-contract`).
- **A test whose assertions follow from a plan holds this machine's speed.**
  An opening's plan weighs the speed of the machine besides the link: how
  fast chunks decode (`progressive-open.js` `decodingRate`, averaged over
  everything the process decoded), how long making every cell's order takes
  (`subsample-plan.js` `geneSampleOf`), how fast a kept gene index was read
  and what a gene index pass's work cost. Measured, each follows how busy
  the machine is, so a test asserting which opening a link plans passed alone
  and failed beside pytest and R checks (`gene-sample-contract`,
  2026-11-04: "a fast link: file order" opened every cell in a `blocks`
  order). Every such rate is read through one module,
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
  and holds the link too (`holdMachineAtRest({link: LINK_AT_REST})`):
  `open-cells-folder-reuse` failed on a loaded Windows runner because its
  900-cell opening took over no index of the 300 cells open before it
  (run 42180511574), and does so on every run with every interval timed 40
  times longer -- the first gene of the 900 cells is then priced at 14.7 s,
  genes are read for 8 of them, and the share opens in the file's cell
  order. A byte reader whose requests a test holds at a gate takes a clock
  of its own (`createByteReader(source, {clock: virtualTime().now})`): timed,
  a request held for as long as a busy machine takes to open the gate was a
  slow link's, and the reader read further ahead than the test's reads
  expected (`shared-reads-cancellation`, HDF5: 12 of 12 runs failed with
  intervals timed 400 times longer, none with the clock).
  **To find such a test**, run the suite with every measured interval made
  longer: preload a module that replaces `performance.now` by
  `base + (now - base) * N` (`node --import ./slow-clock.mjs --test …`, N =
  40 and 400). A test that passes at 1 and fails at 40 asserts something the
  machine's speed decides. (A deadline a test measures on the same clock
  expires N times sooner: read the failure before believing it.)
- **A test undoes what it set up in reverse order.** `node:test` runs a
  test's `after` hooks first registered first, so a test that fakes
  `fetch` (or `Worker`, `window`, `document`) and then makes a source with
  it put the real global back before the source was torn down, and the
  leftover read above went to the real network. Such a test registers its
  undoings with `undoLater(t, undo)` (`tests/helpers/teardown.mjs`), which
  runs them last registered first, each awaited.
- **No contract test reaches the network.** `npm test` preloads
  `scripts/network-guard-preload.mjs` into every test file's process: a
  request to any host but this machine's loopback (`fetch`, and every TCP
  connection beneath it) is refused before a DNS lookup, written to stderr
  with the stack that asked, and the file that made it fails when its
  process exits, naming each URL -- also when the code caught the refusal
  and retried, or asked after its test had ended (the leftover reads above
  were exactly that, and passed). A test fakes every host it reads. The
  one opt-in is for a test whose subject is a real host (none is today):
  `allowRealHost(t, host)` (`tests/helpers/real-network.mjs`) lets that
  host through for that test only. A file run on its own takes the guard
  with `node --import ./scripts/network-guard-preload.mjs --test
  tests/<file>.test.mjs`; with the slow-delivery preload, `NODE_OPTIONS`
  adds to it (`network-guard-contract` checks the guard itself).
- **A contract about cost counts, it does not time.** A threshold on a clock
  is flaky by construction on a shared runner: a 1 MB publish of the
  community-annotation Worker measured 12.1 parses of CPU on one hosted job
  and 5 on a Mac (Web contracts run 42129386525). `blocksRun`
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
- **A test waits for the work, not for a number of turns of the event loop.**
  Crypto, zlib and file work run on other threads, so no count of `setImmediate`
  turns is enough for them. A fake schedule hands the write's own promise back
  (`manualSchedule` in `tests/helpers/memory-history-database.mjs`); a test with
  mocked timers ticks them once the code has set the timer it means to end.
- **Virtual time counts zlib's thread-pool work**
  (`tests/helpers/virtual-time.mjs`). Node's `DecompressionStream` inflates
  off the event loop, where `process.getActiveResourcesInfo()` lists
  nothing, so the clock used to move while a chunk was still inflating. A
  macOS runner that inflated slowly then counted a global heap read in 4
  rounds instead of 3.

This page focuses on “what you can do today” to validate changes.

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
(`../cellucid-python`, `../cellucid-r`, `../cellucid-datasets`) unless an
environment variable names another place (a worktree, a CI checkout
elsewhere):

| Variable | Default | Read by | When the sibling is absent |
|---|---|---|---|
| `CELLUCID_PYTHON_ROOT` | `../cellucid-python` | `scripts/generate-history-golden.mjs` and `tests/history-golden-fixture.test.mjs` (the shared golden history, `tests/fixtures/history_golden/`), `tests/notebook-widget-requests-contract.test.mjs` (the Python widget copy), `scripts/generate-action-catalog.mjs` and `tests/action-catalog-document-contract.test.mjs` (`src/cellucid/action_catalog.json`), the Python viewer specs (`viewer-channel`, `viewer-round-trip`, `viewer-host-selection-delivery`, `notebook-widget-jupyterlab`: the source they run) | the fixture tests skip; the catalog copy is named and skipped (`t.diagnostic`) |
| `CELLUCID_R_ROOT` | `../cellucid-r` | `scripts/generate-history-golden.mjs` and `tests/history-golden-fixture.test.mjs` (the shared golden history, `tests/testthat/fixtures/history_golden/`), `tests/notebook-widget-requests-contract.test.mjs` (the R widget copy), `scripts/generate-action-catalog.mjs` and `tests/action-catalog-document-contract.test.mjs` (`inst/action_catalog.json`), `tests/browser/notebook-widget-jupyterlab.spec.mjs` (the R package source), `scripts/efficiency-floor/served.R` | the golden-history copy test skips; the catalog copy is named and skipped |
| `CELLUCID_DATASETS_ROOT` | `../cellucid-datasets` | `tests/session-preset-control-inventory-contract.test.mjs`: the starting views the sample catalog advertises (`anndata/datasets.json`), read the way the viewer reads them, held to the control inventory | skipped when the sibling is not checked out; a directory named by the variable that is not a `cellucid-datasets` checkout fails |
| `CELLUCID_PYTHON_DOCS` | `../cellucid-python/docs` | `tests/info-tips-contract.test.mjs` (every info tip's link resolves to a heading of these docs), `tests/citation-contract.test.mjs`, `tests/analysis-options-schema-contract.test.mjs`, `tests/action-reference-contract.test.mjs` and `scripts/generate-action-reference.mjs` (the action reference page) | the documentation checks skip |
| `CELLUCID_PYTHON` | none | the Python viewer specs: an interpreter with cellucid's dependencies (JupyterLab and anywidget too for `notebook-widget-jupyterlab`); the package itself comes from `CELLUCID_PYTHON_ROOT`, never an installed copy | the spec skips |
| `CELLUCID_RSCRIPT` | none | `tests/browser/r-viewer-channel.spec.mjs`: an `Rscript` whose library has the package's dependencies | the spec skips |
| `CELLUCID_R` | none | `tests/browser/notebook-widget-jupyterlab.spec.mjs`, R part: an R with IRkernel, cellucid-r's dependencies and pkgload | that test skips |

In CI the Node job sets none of these: it checks out, beside this one, only the
files of each sibling its contracts read (`.github/workflows/test.yml`, sparse:
the two widget copies; both packages' copies of the action catalog, the MCP
tools and the history model; both golden history folders; this site's Markdown
pages; and the sample catalog `anndata/datasets.json` with its starting views).
`tests/ci-sibling-checkouts-contract.test.mjs` holds those lists to what the
generators write, so a new shared file cannot pass on a whole checkout and fail
on every CI job. Only `tests/anndata-copy-cross-language.test.mjs` reads more
(a whole package and a Python with its dependencies); it runs where those are.

(browser-spec-stability)=
### Browser specs that wait for state, not time

A spec that pauses a fixed time for work it cannot see (a redraw, a read, a
notification, a toggle event) passes where that work is quick and fails where
it is not: a hosted runner draws a frame of software GL in a second, macOS
Firefox in seconds. It then fails on another engine and shard from run to run,
and a red gate that teaches nothing but to re-run it slows every change and
every upgrade. The browser suite keeps `retries: 0` for the same reason: a
retry would hide exactly these failures.

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
  whose expiry fails the test. `PENDING_RACE_FIXES` named the race-covering
  sleeps found on 2026-11-01 with the fix each needed; each was replaced by a
  wait on the product's state, a sleep is never added to it, and its `todo`
  test reports it whenever it is not empty. A timer of no delay (a turn of
  the event loop) and a timer that rejects are not fixed waits.
- **No spec reaches the network.** `tests/browser/helpers/outside-requests.mjs`
  watches every context of every spec: a response that names a server address
  outside this machine fails the test, naming it (a route's answer names none).
  It only watches, because a catch-all route would turn off the HTTP cache for
  every spec. The test server serves `index.html` with the suite's own sample
  catalog in place of the published one (`scripts/serve-browser-tests.mjs`), so
  opening the app never depends on the published site, whatever the spec's
  address. A spec plays the hosts it reads from with `page.route`
  (`playAnthropic` for the AI provider); a made-up API key typed in a spec really
  went to api.anthropic.com, and its 401 reached WebKit's console only when the
  answer beat the end of the test (run 42129386525).
- **One cell of the matrix is not a fair test of pixels or GPU-heavy work.**
  Firefox on GitHub's macOS runner takes about a minute for a 1024² fragment
  pass, and has drawn the same bytes held in two GPU objects a few pixels
  apart with the cameras asserted equal (3 of 1.4 million off by 1, 13 of
  921,600 off by at most 2; the same objects drawn twice were identical:
  `tests/browser/helpers/hosted-gl.mjs`). There a spec still holds, exactly,
  every input to a draw (cameras, order, the bytes of the stores, cells drawn)
  and, where it records them, what the compared frames submitted to the GPU
  (below). Its pixel comparison of two such objects is bounded, not skipped: no
  more than one lit pixel in a thousand may differ (`pixelParityVerdict`,
  `lod-rank-order-resize`, `segmented-point-stores`), which still fails on what
  moves whole points (a wrong order, color or camera, cells missing) and
  cannot show a difference of fewer pixels that only that renderer would
  draw. The figure export's point pass drawn in a worker and on the page is
  compared by its two GL contexts only, and `segmented-point-stores` draws
  12,000 points there, still three segments. A test whose budget follows the
  GPU measures it (`smoke-dataset-switch`, `renderer-smoke`,
  `segmented-point-stores`). That those few pixels are that renderer's
  rounding and nothing else has not been established: they were measured on
  the hosted cell only. No other cell has such a treatment.
- **What two frames submitted is asserted before their pixels.**
  `tests/browser/helpers/gl-submissions.mjs` records every draw call of a
  frame: program, mode, vertex count, viewport, whether it rasterizes, and
  every uniform of its program at the call; calls that follow one another
  under the same state count as one, so a store drawn one call a segment
  submits what a store drawn in one call does. `segmented-point-stores` and
  `lod-rank-order-resize` assert the two sides equal on every cell. Two pixels
  that differed on hosted Windows Firefox (`segmented-point-stores`, nightly
  runs 41269083996 and 42180511574) had been read from two kinds of frame: the
  first frame after a level change on one side (it also ran the
  transform-feedback pass that writes the new level's highlight words) and a
  later pulse frame on the other. Four animation frames after a change had
  stood for "a pulse frame has been drawn", which holds where frames are
  fast; ambient frames are paced by what frames cost the main thread and how
  late they leave the page's frames, so on a slow renderer
  the first frame was still the last one drawn. Whether that frame's two
  pixels are a defect of that renderer is not known; the spec now reads
  steady frames (below) and compares them for equality there as everywhere.
- **Under a highlight a scene never settles; its frames repeat.** A
  highlight's pulse draws a frame whenever one is due, so `settled` never
  holds. A spec that only acts next waits with `waitForPageAtRest`, for which
  a scene drawn for the pulse alone is at rest. A spec that reads pixels
  waits for a steady pulse frame (`steadyPulseFrame` in `gl-submissions.mjs`):
  nothing invalidated, and the last two frames drawn since the change
  submitted the same, with the pulse's phase pinned
  (`setHighlightAnimationTime`).
- **A wait is as long as the machine, and the test's timeout is its fence.**
  A loop that asks the product's state once an animation frame does not count
  its frames, and sets itself no clock. `segmented-point-stores` waited up to
  900 frames a scene for a condition that could not hold (it compared a
  `drawOrder` the store order does not have) and went on whatever they
  showed: a quarter of a minute on a developer machine, and the whole
  ten-minute budget where a frame takes seconds (hosted macOS and Windows
  Firefox, "Test timeout of 600000ms"). Counts that ran out are in the
  history of this workflow: "The live store never adopted its tree order"
  (`frustum-culling-parity`, Linux Firefox), a frame read before any had been
  drawn under its rings and 126,568 pixels off (`highlight-pulse-layer-parity`,
  Linux Firefox, a loop that went on in silence after 600 frames), "the
  sidebar toggle started no slide" (`presented-layout`, Linux WebKit: a 300 ms
  slide that had started and ended between two looks; the slide sampled now
  does not end by itself). `tests/browser-spec-counted-waits-contract.test.mjs`
  (`npm test`) finds every `for` loop that counts to a number, awaits a frame
  or a timer and leaves early, and every loop bounded by a deadline of its
  own; each must be listed with why its count is the subject. The same test
  holds the time a spec gives an assertion or a wait: the suite's own is a
  minute for an assertion and the test's timeout for everything else, and a
  spec that gave one less (`{ timeout: 10_000 }`, 110 of them) said how fast
  the machine must be.
- **A selection gesture is made once, when the page can take it.** A press on
  the plot is read by the tool in use against the frame on screen; the tool's
  mode is switched by an action, and a pick is of the frame that was drawn.
  Specs repeated a gesture until its tool armed and called it dropped after
  three seconds: on a runner that armed later a second gesture landed on top
  of the first (`select-bar`, macOS Chromium; `uint32-categories-end-to-end`,
  Linux WebKit and Firefox). `tests/browser/helpers/selection-gesture.mjs`
  waits for every action asked so far (`getActions().idle()`) and for the
  page to be at rest, makes the gesture once, and waits for the tool's step
  controls.
- **A slow machine on a developer machine.**
  `CELLUCID_BROWSER_TEST_CPU_THROTTLE=6` slows the main thread of every page
  six times (Chromium only, the DevTools protocol's CPU throttling;
  `tests/browser/helpers/test.mjs`), so work is ordered as on a hosted
  runner: what a spec read right after an action is then not there yet. It
  slows a page's own thread, not its workers or the GPU.
- **A page is closed at a defined point.** A change of the history is handed
  to the browser's storage some tens of milliseconds after it is made, and is
  then kept whole whatever becomes of the page
  ({doc}`29_action_layer_domains_and_history`, "All of a write or none"). A
  spec about what a closed or reloaded page kept waits for that first
  (`window._cellucidHistories.flush()`, or the saved state listed with what
  it holds); a spec about a page that goes at once asserts what the product
  promises then: all of the last write or none of it, never a part
  (`saved-state-close.spec.mjs`, which closes a tab at several moments of one
  write in every engine).
- **A page is parsed only once it is styled.** WebKit does not hold
  DOMContentLoaded for a stylesheet; Chromium and Firefox do. The page boot
  starts the code at DOMContentLoaded, so with the stylesheet late (held 4 s:
  the code started 0.3 s in) WebKit ran the bootstrap in a page with no
  styles: design tokens had no value (the analysis tools failed to load on
  one) and markup hidden by a class was displayed. A spec that navigated and
  waited for the welcome dialog to be visible saw the markup, pressed Escape
  before any listener existed, and the page showed the dialog later over a
  page the spec believed free of it (Linux WebKit, run 42180485733 and
  earlier: `history-graph`, `highlight-pulse-layer-parity`,
  `knn-missing-connectivity-warning`, `my-datasets`). `index.html` now carries
  a style gate (`CELLUCID_STYLE_GATE`): a script straight after the stylesheet
  link, which every engine runs only once the stylesheets before it are in,
  so nothing of `<body>` is parsed and DOMContentLoaded does not come before
  them; a stylesheet that did not load is reported on the code-load failure
  screen instead of starting the code without it. `style-gate.spec.mjs` holds
  the stylesheet, and one it imports, for 2.5 s in every engine;
  `tests/startup-code-graph-contract.test.mjs` holds the gate's place in the
  page. `dismissWelcome` waits for the state the product keeps
  (`#welcome-modal[open]`: the welcome is a `<dialog>`, never displayed
  until opened, styled or not).
- **Offline only through `goOffline`.** A page still loading something of its
  own when the connection goes reports that failure too: the analysis tools it
  imports as it starts (hundreds of files, still arriving on a hosted runner:
  `community-annotation-startup-lifecycle`, Linux WebKit, run 42180485733), a
  read of the opening still under way. `tests/browser/helpers/offline.mjs`
  waits on every page of the context for the page to be at rest and for
  `window._comparisonModule`, which the page sets in the continuation of that
  import, then takes the context offline;
  `tests/browser-spec-offline-contract.test.mjs` allows `setOffline(true)`
  nowhere else.
- **A click whose press and release a layout change split is reported.**
  Playwright checks what a press lands on and releases a moment later at the
  same point. The Session panel's Dataset box fills in while the start
  dataset opens and moves everything under it (the foot of the panel by 718 px
  in six steps): a release then lands on something else, the browser sends
  the `click` to their common ancestor, and the control pressed hears
  nothing, with no error (`ui-visual-consistency.spec.mjs`, the license link,
  Linux WebKit, run 42180485733). `tests/browser/helpers/split-clicks.mjs`
  watches every context of every test for a trusted mouse press and release
  at one point whose pressed element is no longer under it, and fails the
  test naming both elements (`split-clicks.spec.mjs`). The cure is the
  spec's: click once the page has stopped moving (`waitForPageAtRest`).
  `openSessionSection` (`panels.mjs`), which opens a header under that box
  whatever a spec is about, waits for the interface to be built
  (`window._cellucidUi`), clicks, and clicks again only for a split click
  reported on its own header.
  A drag has the same hazard one step earlier: `keyboard-move-lifecycle` read
  a panel summary's box and pressed two seconds later, on the Agents panel,
  which the history rows of the test's own steps had moved there, and tore
  that panel off (hosted macOS Firefox, run 42202865024). It now waits for
  the steps to be recorded, moves the pointer onto the summary itself and
  checks what the press landed on before it drags.
- **A change made through the interface is done when the dispatcher says
  so.** A control's change is an action, applied when those before it have
  finished and with what it reads on the way. A spec that goes on through
  the page's internals (a domain restored directly, a pick) waits for
  `getActions().idle()` first (`whenActionsDone`, `dataset-ready.mjs`):
  `state-domains-roundtrip` restored a domain while the field chosen just
  before was still being applied, and was overtaken by it (WebKit, 1 run in
  14).
- **Repeat before pushing.** `npm run test:browser:repeat` runs specs many
  times in each engine, Chromium, Firefox and WebKit one after another, with
  one worker and no retries as CI runs them, in browser processes no
  longer-lived than CI's (six files and forty tests, the bounded runner's
  limits: a WebKit run past that stopped answering a new page's first
  navigation, a failure no CI cell can have), and fails on any single failed
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

### Long, sharded and opt-in suites

- **Every browser spec file declares tests.** Playwright says nothing about a
  spec file that holds none, and a spec cut to its imports once lost thirteen
  tests while the run stayed green; `tests/browser-specs-hold-tests-contract.test.mjs`
  (`npm test`) fails on a spec file under `tests/browser/` that calls no
  `test(…)` / `test.describe(…)` (or hands `test` to a helper that declares
  them, `declareWalkTests`) or does not end with a newline.
- **The control walk** (`tests/browser/helpers/control-walk.mjs`) walks every
  control the page shows and fails on a control that changes state with no
  step, an outside change, more than one step after the step shown, or an
  undo that does not put back exactly
  ({doc}`29_action_layer_domains_and_history`, "Drift guards"); new steps on
  another branch none of which follows the step shown (Clear history keeping
  the state, Branch without this step replaying the later steps) are judged
  as an edit of the history.
  It is one test per sidebar section, one for everything outside them and the
  graph-and-vectors test, spread over `control-walk-1-of-4.spec.mjs` …
  `4-of-4`, and runs on Chromium only. Run one part with
  `npx playwright test tests/browser/control-walk- --project=chromium --grep "<section id>"`,
  or the four files side by side with `--shard=<i>/4`; each part prints
  where its time went, by section, its slowest controls and phase (about 15
  minutes for the whole walk on a GPU). A control is judged by what it does,
  never by how long it takes: under the runners' software GL a control with
  the velocity overlay flowing or smoke animating takes 15 to 25 s, nearly
  all of it the page's frames. What bounds the walk is its progress: a wait
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

### What the hosted CI runners differ in

A spec that passes on a developer machine can fail only in CI, for reasons no
change of the spec's subject explains. These were measured on the Web
contracts workflow (`cellucid/.github/workflows/test.yml`), run 41213981681
(push of 2026-10-31) and the nightly runs before it:

- **Software GL.** Hosted Linux runners draw on a software rasterizer:
  Chromium on SwiftShader, Firefox on llvmpipe, WebKit behind the masked name
  "Apple GPU" at llvmpipe's speed (262,144 points in 30 ms a draw). The shard
  running `webgl2-runtime.spec.mjs` logs the renderer of each cell. macOS
  WebKit draws on the virtual machine's GPU; macOS Firefox's GL is about 2,500
  times slower than a developer Mac's.
- **The page measures its device and acts on what it measures.** Level of
  detail and frustum culling turn on when the GPU's share of a frame of every
  cell is measured over the interactive budget, and stay on only where their
  frames are measured faster (`rendering/full-detail-frame-cost.js`,
  `render-controls.js`); an opening is progressive when showing its cells is
  predicted slow from what showing a cell cost earlier openings
  (`data/opening-cost.js`); worker pools follow the cores the browser
  reports, memory budgets the memory it reports. On software GL those
  measures came out otherwise, and at other times, than on a developer
  machine, so specs about none of those choices saw the page change state
  under them, at a time the runner decided. Five CI rounds patched one spec
  each (a person's click held the choice still, a dataset's file was held
  back until it was, a camera zoomed in first, an "at once" opening waited
  for the start sample's own record); each round found the next spec
  (2026-10-28 to 11-01: the verdict reached before `analysis-cancel` watched,
  a failed publication in `camera-defaults` turned the page's own culling
  off, the start sample's record made `text-matrix-stream`'s "at once"
  opening progressive, Web contracts run 41307458156). The cause was one:
  a test ran on whatever device the runner was. **Every spec now runs on a
  fixed device**, the `device` option of the `test` fixture
  (`tests/browser/helpers/device.mjs`, installed on every context the test
  opens, `browser.newContext` included):
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
      about a slow link gives its link here (`{latencyMs: 400, bytesPerMs:
      32000}`) instead of delaying answers.
    - `decodeBytesPerMs` (1 GB/s), `orderMadeMs` (0),
      `keptReadBytesPerMs` and `indexPassMs` (null: none measured).
    - `openingMsPerCell`: what showing a cell costs (the machine's
      `shownCellMs`): nothing measured by default, so every opening is
      priced by its reads alone, and the page's own records are not made.
      `browser.newContext({ device: { openingMsPerCell: 2 } })` makes an
      opening progressive; `'measured'` lets the page record (the record's
      own spec in `open-cells-in-place`).

    Why the link and the decoding rate are held too (they were left
    measured until 2026-11-06, "every fixture is served from the machine
    running the spec, far inside the opening target"): a request a loaded
    machine answers itself late is a slow link to the reader. In
    `open-cells-in-place`, beside five other runs, the 1,500-cell opening
    read genes for 864 of its cells, which opens the share in the file's
    cell order, so Apply showed it as a new reading instead of in place.
    `held-machine.spec.mjs` holds the file's answers back 370 ms: measured,
    the raised share's genes are read for a sample and it is read anew, in
    every engine; on the fixture's device it opens in place with genes for
    every cell. In `kept-gene-index` a hosted macOS WebKit runner read
    genes for a sample on one of the two visits, twice: the first visit's
    gene on 606 of its 1,000 cells (run 42180511574), and the later visit's
    on 811 of them in the file's cell order, where the first had read every
    cell in the share's own order (run 42202865024). Every value of both
    visits was its cell's (checked against the fixture's formula, 1,000 of
    1,000 and 811 of 811); the spec compares the visits position by
    position.

    WebKit only: a spec of `persistentStorage` runs with one worker. WebKit
    keeps a persistent context's origin storage outside its profile folder
    (`helpers/test.mjs` `clearOriginStorage`), and two workers running such
    specs at once lost what they kept (`kept-gene-index` and `kept-parts`:
    6 of 12 runs; CI runs one worker).

    The holds reach a page that runs its sources (this machine). A page
    opened as a deployed site runs its bundle, whose module instances an
    init script cannot name, and measures everything: a spec opened there
    asserts nothing that follows from a plan.

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
  | The GPU's time for a level of a backdrop's sky (`rendering/backdrops/backdrop-look-renderer.js`) | How far the sky is refined (its faces' resolution) | Measured | not verified for every spec: the look shows from its first level on any GPU, and a slow one stops refining sooner (a coarser sky), so a spec may assert a look only against the same page |
  | What a neighbor row cost from the file, and a pass over the held edges (`app/main.js` `knnRowReader`) | Which of the two gives a cell's neighbors | Measured | both give the same rows |
  | The cores | Pool sizes | Held: `hardwareConcurrency` (page `navigator`, the machine's `cores` in workers) | held |
  | The memory the browser reports | Budgets, the cells an auto-fit opens | Held where reported (`deviceMemory`, Chromium); Firefox and WebKit report none and take 8 GiB | determined |
  | An allocation the engine refused; a tab that crashed while opening | The memory allowance halved | Measured | only where the runner really cannot allocate, which the page must then follow |
  | Free storage (`storage.estimate()`) | A gene index on disk, a copy kept or refused | Measured | the fixtures are megabytes; a runner without room fails by a named storage error |
  | Elapsed time in loops that yield (index passes, selection queries, noise generation, inflate) | When a loop gives the thread back | Measured | scheduling only: the slices are bounded by time, so a slow runner does less a slice, not longer slices |

  The per-spec holds this replaced (`drawEveryCellByHand`,
  `drawEveryCellFromTheStart`, a camera zoomed in before a simulated slow
  GPU, an "at once" opening that waited for the start sample's record) are
  gone. A spec that sets level of detail or culling through the viewer
  (`setAdaptiveLOD`, `setFrustumCulling`) still makes that opening's choice,
  and the render controls leave the rest of its reports alone. Specs that
  count or classify draws tell a point-pass draw by its program (one of the
  renderer's `programs`), not by its count: the GPU gather of a tree's draw
  order (`high-perf/gpu-permute.js`) draws `pointCount` points too, with
  blending off, into its own target, and skip transform-feedback passes
  (`RASTERIZER_DISCARD`); and they count per animation frame
  (`document.timeline.currentTime`), not per pointer move: a slow host draws
  several frames per move. `full-detail-frame-cost.spec.mjs` asserts each
  verdict exactly, number for number: before the timer, fences that
  signaled at a simulated time still let the runner decide (a frame
  measured 150 ms faster was called no faster on Linux Chromium; 40,000
  cells under a flat 60 ms measured 124 ms against 96 ms with level of
  detail and culling, and were kept, run 41255961977).
- **Controls that move under the pointer.** Playwright clicks where an
  element was when it last checked it stable; a control that moves in
  between is missed. The Dataset box above the Session panel's controls
  changed height on its own: its read line wraps "Reading… 7.5 MB of 22 MB
  read in 8 requests" over two lines and its plain total over one, 1.5 s
  after the last bytes, and the box gave its facts up for "Reading … again…"
  while Apply read the dataset again. The second Apply of a share round trip
  then landed 15 px below the button and did nothing (macOS WebKit, run
  41255961977). The box now keeps its room (`dataset-card.js`), and
  `session-restore-gene-sample.spec.mjs` asserts that Apply stays put from
  the click until the read line settles. The highlight count above the
  Analysis panel moves too: with level of detail on it counts the highlighted
  cells the level draws ("333 of 1,000 highlighted cells visible"), a line
  that wraps, so the automatic choice turning level of detail on and off
  moved Run 42 px between the press and the release of its click, which then
  reached no button (`analysis-cancel.spec.mjs`, Linux WebKit, run
  41292953076), and the `Force LOD level` row appeared and went with it
  (40 px). No control may move because a line the page changes by itself
  changes length: the count holds one line of its own (cut, its title says
  it whole), the selection bar's status lines hold one line, and the force
  row stays, grayed out while level of detail is off;
  `analysis-cancel.spec.mjs` asserts that Run and the count's header stay
  put while the page turns level of detail on and off (it failed before,
  Run moving 40 px). A sidebar scan that records every control's place each
  frame found no other control moved by the page alone after an opening.
- **Work the page starts on its own after an opening.** The cell info card
  asks for a large pane's lookup tree once the page is idle, so "no tree
  before the first hover" holds only while the card is off
  (`cell-info.spec.mjs` turns it off while 300,000 cells open, then switches
  it on and hovers in one task). The gene list is read after the first
  view, and listing it refreshes the interface, which publishes the active
  view again: a spec that counts publications starts once
  `data-startup-phase` is `ready` and `state.whenGenesListed()` has settled
  (`renderer-smoke.spec.mjs` smoke settlement, Windows Firefox). The
  opening's frame-cost measurement asks for its next sample frame once the
  last one's fence passes, seconds later on macOS Firefox:
  `getFrameScheduleState().settled` is false while a sample's fence is
  awaited (`frameCostPhase` names the measurement's phase), so
  `sceneQuietness` no longer called a scene quiet that the measurement was
  about to draw (`idle-frame-scheduling.spec.mjs`, three tests, nightly run
  41269083996). The measurement now asks for frames only where an answer
  could change what is drawn: the first frame of every cell it would sample
  is looked at (`presentedDrawsFewer`), and where no pane has a reduced
  level for its cells and every pane's cells lie in its view, it stops
  (`nothing-fewer`) and the next framing is looked at afresh; culling is
  tried only where a cell may lie outside the view. On macOS Firefox a frame
  of 8 to 300 cells took 1.5 to 5 s (`smoke-dataset-switch.spec.mjs` trace,
  nightly run 41269083996), so the 3 to 6 frames each opening's measurement
  asked for were tens of seconds a test. The scene is not settled either
  while a measurement is armed: arming asks for a frame at once (a
  measurement armed for the next framing on a still scene once waited for
  whatever drew next, a hover or a pulse, and could turn level of detail on
  under it), and `settled` holds until it has answered. One helper waits
  for all of it, `waitForPageAtRest` (`tests/browser/helpers/dataset-ready.mjs`):
  startup `ready`, the gene list read, the frame loop settled (the
  measurement included; a scene drawing for a highlight's pulse alone is at
  rest too) and every view on its final geometry. A spec that sets the
  camera or the GPU waits for it first: a camera set while the sample's own
  published view was still to come was undone by it, and the zoomed framing
  of `full-detail-frame-cost` and `graphics-context-restore` was never
  measured (WebKit, full suite, 2026-11-01).
- **Where a 2,500-times-slower GPU spends a spec's time** (macOS Firefox,
  `smoke-dataset-switch.spec.mjs` trace): every frame, whatever it draws,
  1.5 to 5 s at 1280 × 720 (a test's two-frame waits 3 to 10 s each); each
  `getError()` read 3 to 8 s (it waits for the GPU); a smoke density splat
  of 8 to 300 cells 2 to 6 s; the cloud noise at its smallest size (32³)
  31 s once, already paced by fences and sized by the GPU's speed
  (`gpu-noise-generator.js`). The per-frame cost is per pixel, not per
  cell: no automatic choice the page has (level of detail, culling, the
  noise batch size) can lower it; only drawing fewer pixels could (a render
  scale or antialiasing chosen from the measured frame cost), which is a
  decision about how the picture looks, left to the owner.
- **Teardown while work is still in flight.** On a slow runner an opening's
  first frame can still be awaited when the page is closed: the frame-phase
  callback that records it ran after `dispose()` and asked the disposed
  viewer for its statistics (`open-cells-in-place.spec.mjs`, Linux WebKit).
  Application teardown now retires the live dataset publication
  (`createLatestDatasetPublicationContinuationOwner().retire()`), whose
  signal ends that wait and every read of it, and viewer disposal drops the
  deferred WebGL error checks still waiting for their fence, whose polls
  otherwise asked for animation frames after the disposal fence
  (`viewer-terminal-disposal.spec.mjs`, macOS Firefox). A failed in-place
  Apply of Cells to open is announced once ("… could not be read again with
  these cells …"), not left to an unhandled rejection
  (`tests/dataset-reopen-failure-contract.test.mjs`). The failure itself
  (Linux WebKit, "WebKit encountered an internal error", then "The selected
  H5AD file is empty"): Apply sent the picked `File` to the reader worker
  again, the worker's key of a file (name, size, last modified) did not
  match the new copy's (it read 0 bytes), so the store was not taken over
  and the copy was opened anew, as empty. Reading again with other cells
  now sends no source to a live worker, which reads the one it holds
  (`heldSource`, `tests/reader-reopen-held-source-contract.test.mjs`), and a
  picked file of 0 bytes that cannot be read is said to be unreadable, not
  empty.
- **What the product takes back.** A box the page unticks itself (a failed
  velocity allocation) is clicked, not `check()`ed: `check()` fails when the
  page's answer lands before it looks again (Windows Chromium and Firefox).
- **Core counts.** The analysis worker pool has a worker per core less the
  page's and the reader worker's, so a 3-core macOS runner has one and a
  density curve is summed in one slice there;
  `analysis-density-workers-identity.spec.mjs` requires the slice count the
  pool's size plans, and the bit-for-bit identity at any count.
- **Browser-reported metrics.** A browser records no largest contentful
  paint if the first input (the Escape that closes the welcome dialog) came
  before a text or image paint, so `analytics-production.spec.mjs` requires
  LCP only when the browser recorded one; the page takes each performance
  observer's undelivered entries (`takeRecords`) before it disconnects them
  on the first hide. A routed response the page let go of by navigating
  away is disposed, and fulfilling it fails (Windows Chromium): the route
  handler lets that one error pass.
- **Pixel-parity evidence.** A difference in `unplaced-cells.spec.mjs`
  attaches both frames, their difference and what each was drawn from
  (camera, level of detail, cells drawn, draw order, GL renderer) to the
  report. That evidence found its cause: the reference was drawn in a cell
  lookup tree's order and the second frame in source order (Linux Firefox,
  221,584 bytes in runs 41213981681 and 41255961977; Linux Chromium too in
  run 41292953076), because the cell info card asks for that tree when the
  browser is idle, which on a runner came after the spec had changed the
  positions. Translucent points overlap in draw order, so two frames are
  comparable only when drawn in the same one: the spec turns the card off
  and asserts both frames' order before comparing pixels, and
  `isPresentedGeometrySettled` waits for any tree that will reorder a store.
  Frames compared across a reload, a restore or an opening are drawn alike
  because the device is fixed (the automatic choice followed a measured
  frame time, `graphics-context-restore.spec.mjs`, macOS Firefox, nightly run
  41269083996). A spec that waits for "only the highlight pulse is left"
  asks `animationCause === 'highlight'` with nothing invalidated, never the
  cause of the last frame drawn: right after a restore that was still the
  pulse from before the loss, the restored context had drawn nothing, and
  Linux WebKit read every pixel back as 0 (push run 41292953076). That says
  the frame of the change has been drawn, not that a later one has: a spec
  that compares the pixels of two such scenes reads steady pulse frames
  (`gl-submissions.mjs`, above).
- **A one-shot GL failure** injected into an allocation (`createFramebuffer`,
  `createTransformFeedback`) is aimed at the caller under test by its stack
  (`velocity-multiview-lifecycle.spec.mjs`): the spatial tree's GPU gather
  allocates both at a time the GPU's speed decides.
- **WebKit's navigation wording.** A `page.goto` or `page.reload` stops the
  loads the left document still had going, and WebKit reports each as that
  document's failure ("Cannot load … due to access control checks.",
  "Importing a module script failed."; Playwright passes the fetch one as a
  page error), a worker's module graph included: a page that had turned
  culling on 190 ms before a spec's `page.goto` was still loading its
  spatial-tree worker's modules (push run 41307458156). A worker the page
  ends itself (`terminate()`) reports nothing. Every page of every test
  context, a spec's own `browser.newContext()` included, is prepared by
  `tests/browser/helpers/test.mjs` (`setAsideLeavingReports`,
  `tests/browser/helpers/leaving-document.mjs`): while its `goto`, `reload`,
  `goBack` or `goForward` from a document runs, those wordings, and only
  those, reach none of its listeners; no spec handles them itself
  (`tests/browser-leaving-document-contract.test.mjs`). A read of a kept
  file through a blob is stopped and reported the same way ("Cannot load
  blob:http://… due to access control checks.", a page error): a page left
  while its Stored listing read a kept dataset's last-used time failed
  `kept-gene-index` in WebKit (3 of 120 runs with one worker, on main as
  well; a page that reads such files and is left made 1,169 reports in 40
  navigations). Playwright's WebKit cuts that sentence at the colon of
  `blob:` (the error's name "Cannot load blob", its message "ttp://…"), so
  the helper matches a page error by the first line of its stack, the
  sentence as WebKit said it.
- **Open: a Chromium page that stops answering for about a minute on a
  loaded machine.** With the machine's load average between 85 and 316 (ten
  cores, other work beside three Playwright workers), a Chromium page now
  and then answered no `page.evaluate` for about 60 s (the spec's 60 s poll
  for its open cells then ended at 0) and went on afterwards: 2 of 198 runs
  of `open-cells-in-place`, `kept-gene-index` and `kept-parts` on main
  (e0aaaebc5), 13 of 426 on the branch that holds the machine, none of 140
  in Firefox or 176 in WebKit. The page had started (its viewer made, in
  one trace its first frame drawn) and the reads resumed by themselves.
  Not traced; a hosted runner runs one worker, and no CI run of these specs
  on Chromium failed so in the forty runs looked at.
- **WebKit refuses a read of a file being rewritten.** A page's read of a
  kept file through a blob, while a worker rewrites the file through a
  synchronous access handle, is refused now and then and logged ("Failed to
  load resource: … (WebKitBlobResource error 1.)", 5 of 2,645 reads;
  Chromium and Firefox refuse none): the Stored listing read a kept
  dataset's last-used time while a later visit's reader wrote it
  (`kept-parts`, one run of 120). The write now takes the lock the listing
  reads under (`transcode/space.js` `filesOf`): none of 3,001 reads refused.
- **Shared network addresses.** GitHub's API answers 60 unauthenticated
  requests an hour per address, and a runner's address is shared: its
  refusal (403 with `x-ratelimit-remaining: 0`), which the GitHub provider
  handles by reading at the ref, is set aside by
  `data-source-coverage-live-servers.spec.mjs`; nothing else from that host is.
- **Sibling checkouts.** The Node job checks out `cellucid-python`,
  `cellucid-r` and `cellucid-datasets` at their default branch when it starts.
  A cross-repository change whose sibling push lands later fails until the next
  run; `notebook-widget-requests-contract.test.mjs` names each copy's commit
  and runs its behavior contracts on both copies, so the failing copy is the
  one behind. Push the siblings first.
- **No WebKit on Windows.** The nightly and dispatch matrix runs every engine
  on every operating system except WebKit on Windows (removed 2026-10-31).
  Playwright's Windows WebKit is the WinCairo port, not Safari, which does not
  exist on Windows: it has no origin-private file system
  (`navigator.storage.getDirectory`) and no `OffscreenCanvas`, so the specs
  that need them failed there on every run. Safari is covered by WebKit on
  macOS (every push) and Linux.

Playwright's headless Chromium (1.62, and 1.65 on 2026-11-07) draws with
SwiftShader on an Apple M1 Pro too (`ANGLE (Google, Vulkan 1.3.0 (SwiftShader
Device …))`, 2026-11-01) unless launched with `--use-angle=metal`, and it
composites the page in software (`SystemInfo.getInfo`: `gpu_compositing:
disabled_software`, `rasterization: disabled_software`): the configuration
Chromium chooses itself when SwiftShader is its only GL. The plain `chromium`
project is therefore the local reproduction of software GL. It is slower than a
GPU but faster than a hosted runner, so a race that needs a runner's speed may
still not show.

Do not force SwiftShader onto the whole GPU process instead (`--use-gl=angle
--use-angle=swiftshader --enable-unsafe-swiftshader`): that is a machine no
runner and no reader has. It also turns on GPU compositing and GPU
rasterization, with Skia Graphite on Dawn over SwiftShader, so every
composited frame is drawn on the CPU by the GPU process's main thread, the
thread every WebGL context's commands run on, and on macOS that thread holds
each frame about 100 ms until it is finished
(`SkiaOutputSurfaceImplOnGpu::EndAccessImages`). An animation the compositor
thread runs on its own (the opening's spinner and indeterminate bar) asks for
a frame at every vsync, so the thread is busy with frames all the time and
WebGL commands run only in the gaps: through the opening, each synchronous
WebGL call of the page waited 74-263 ms (Chrome trace, 2026-11-07). When the
figure-export worker created its WebGL2 context at the same time
(`figure-export.spec.mjs`, "export cancellation settles worker work and
cleans download staging", which asks the worker for `spriteZoomLimit` as the
page loads), the page's next `getError` was not run for 88 s to over 150 s:
the opening could not end, so its animation never stopped. All 12 runs hung;
all 8 with `reducedMotion: 'reduce'` (which stops every animation) passed, as
did all 9 in the three engines' default configurations. A spec that fails only
under the forced configuration measures that configuration, not the app. To
run the
automatic level of detail and culling in every engine as a runner does, run
the specs on a copy of the tree with `INTERACTIVE_FRAME_BUDGET_MS` set to 0 in
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

## Documentation build sanity check (optional but recommended)

If you changed docs pages, you can build the `cellucid-python` docs locally:

```bash
make -C cellucid-python/docs html
```

This ensures MyST links/toctrees are valid.

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
