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
runner allows; `tests/browser-shard-capacity-contract.test.mjs` plans every
shard of the matrix for every engine and runner from one Playwright listing
(`tests/helpers/shard-plan-reporter.mjs` replicates Playwright 1.62.0's
`createTestGroups` and `filterForShard`, refuses another version, and is
compared test for test with real `--shard` listings of every other shard
and the last, which hold every boundary). It listed each of the 21 shards
in a process of its own before: 26 s idle and 166 s under load, past the
suite's 120 s a test; it now takes 5 processes, 7 to 15 s idle and 72 s at a
load average of 40 to 54.

### Contract tests that hold on every Node version and operating system

What the contract suite once got wrong only on some CI cells (fixed
2026-11-01):

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
  a deferred the test releases (`holdDatasetOpening` and `releaseAndDrain` in
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
  - `openingMsPerCell`: what showing a cell costs, held in the page's module
    (`holdShownCells`): nothing measured by default, so every opening is
    priced by its reads alone; the page's own records leave it as it is.
    `browser.newContext({ device: { openingMsPerCell: 2 } })` makes an
    opening progressive; `'measured'` lets the page record (the record's own
    spec in `open-cells-in-place`).
  - `hardwareConcurrency` (4) and, where an engine reports it (Chromium),
    `deviceMemory` (8), so the page's pools and budgets are sized alike
    everywhere. Workers read their own `navigator`, which an init script
    cannot reach: a pool a worker sizes still follows the machine, and its
    results do not depend on its size. The link's speed and the decoding
    rate are measured too and left so: every fixture is served from the
    machine running the spec, far inside the opening target.

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
  signalled at a simulated time still let the runner decide (a frame
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
  row stays, greyed out while level of detail is off;
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
  Linux WebKit read every pixel back as 0 (push run 41292953076).
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
  (`tests/browser-leaving-document-contract.test.mjs`).
- **Shared network addresses.** GitHub's API answers 60 unauthenticated
  requests an hour per address, and a runner's address is shared: its
  refusal (403 with `x-ratelimit-remaining: 0`), which the GitHub provider
  handles by reading at the ref, is set aside by
  `data-source-coverage-live-servers.spec.mjs`; nothing else from that host is.
- **Sibling checkouts.** The Node job checks out `cellucid-python`,
  `cellucid-r` and `cellucid-datasets` at their default branch when it starts.
  A cross-repository change whose sibling push lands later fails until the next
  run; `notebook-widget-requests-contract.test.mjs` names each copy's commit
  and runs its behaviour contracts on both copies, so the failing copy is the
  one behind. Push the siblings first.
- **No WebKit on Windows.** The nightly and dispatch matrix runs every engine
  on every operating system except WebKit on Windows (removed 2026-10-31).
  Playwright's Windows WebKit is the WinCairo port, not Safari, which does not
  exist on Windows: it has no origin-private file system
  (`navigator.storage.getDirectory`) and no `OffscreenCanvas`, so the specs
  that need them failed there on every run. Safari is covered by WebKit on
  macOS (every push) and Linux.

Playwright's headless Chromium (1.62) drew with SwiftShader on an Apple M1
Pro too (`ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device …))`, 2026-11-01)
unless launched with `--use-angle=metal`. Reproduce software GL locally with Chromium under SwiftShader: a config that
extends `playwright.config.mjs` with
`launchOptions: { args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] }`
for the `chromium` project. It is slower than a GPU but faster than a hosted
runner, so a race that needs a runner's speed may still not show. To run the
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

### F) Sessions

1) Save a session (`.cellucid-session`).
2) Reload the page and load the session.
3) Confirm:
   - camera/layout restored
   - active fields and filters restored
   - highlight pages/groups restored (or lazy-restored shortly after)

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
  colour buffer) and bounded draw submissions;
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
