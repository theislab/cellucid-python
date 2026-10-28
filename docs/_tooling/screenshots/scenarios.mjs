// The named scenarios `capture.mjs` can reproduce.
//
// A scenario is data, not a script. Everything a reviewer needs in order to
// judge whether the image shows the right thing — which dataset, which clicks,
// which region, where the pointer is — is visible in one object, and the same
// object is what regenerates the image later. Adding a figure to the
// documentation means adding an entry here, never hand-cropping a screenshot.
//
// FIELDS
//   id          Stable handle used on the command line and in `captures.json`.
//   topic       Directory under `docs/_static/screenshots/`.
//   name        File stem. The emitted file is `<topic>/<name>.<format>`.
//   describes   One sentence. Becomes the `:alt:` text, so it must describe
//               what is visible, not what the reader should conclude.
//   dataset     Dataset id in `cellucid-datasets/anndata/datasets.json`.
//   query       Extra launch parameters, merged into the URL.
//   keepWelcome Leave the first-run overlay up (default: dismiss it).
//   steps       Ordered actions. See `applyStep` in `capture.mjs` for the
//               vocabulary; `{ run: page => ... }` is the escape hatch.
//   shot        { mode: 'window' }
//               { mode: 'element', selector, padding }
//               { mode: 'region', selectors: [...], padding }
//   cursor      { on, at, dx, dy, state } — a real hover plus a drawn pointer.
//   rings       [{ on, label }] — optional highlight, off unless asked for.
//   settle      { attempts, gapMs } or { frames } for animated scenes.
//   emitWidth   Override the width rule. Never larger than the native capture.
//   format      'png' (default) or 'webp'. WebP is for images a page places
//               with its own `<img>` markup (the home page); every `{figure}`
//               is PNG.
//   quality     WebP only: the encoder quality in (0, 1], fixed per scenario.
//   servedAs    Show the application under this origin instead of the run's
//               own (answered by the run's server through a route): the
//               Google and Microsoft sign-in steps exist only on the origins
//               those providers registered, `http://127.0.0.1:8765`.
//   entry       async urls => address. Open this address (a Viewer URL
//               `cellucid serve` or `cellucid_view()` printed) instead of the
//               run's application with its sample catalog.
//   awaitStartup false: do not wait for startup to finish or for the welcome
//               overlay; the subject is a startup still under way.
//
// SIZE RULE (enforced in `capture.mjs`, restated here because this is the file
// people edit): the emitted width is the captured region's CSS width times two,
// capped at 1440. Do not set `emitWidth` to make an image "look right" on one
// page — the figures on a page are sized by what they contain, and a panel that
// is 246 CSS px wide is a 492 px file.
//
// READY SIGNALS that are known to work on this application:
//   PANCREAS_READY     both of the two below; every pancreas scenario starts here
//   #filter-count      "Showing all 3,696 cells"   — data loaded and counted
//   html[data-startup-phase="ready"]                — the sample's published
//                      default session has been applied and the viewer runs
//   #stats             "Points: … • Field: … • Centroids: …"  — colouring live
//   #legend            present once a field paints a legend
// `#dataset-cells` is NOT a ready signal: it abbreviates to "4K".
//
// The field selectors are rebuilt once the obs field list arrives, which
// is after `#filter-count` settles. Always put `{ waitForStable:
// '#categorical-field' }` before selecting a field; without it the selection is
// silently discarded and the image shows "Field: None".
//
// The sidebar accordions are `<details class="accordion-section">`. On first
// load these are open — session, coloring-filtering, highlighted-cells,
// page-analysis — and these are closed: visualization, navigation,
// compare-views, community-annotation, cinematic-camera, figure-export,
// shortcuts, benchmark. A control inside a closed section is not visible, so
// a step that uses one opens its section first. Use `openSection` /
// `closeSection` (or `onlySections(...)`) rather than clicking a summary, so a
// scenario is idempotent.

import { OVERLAY_ROOT_ID, paintOverlay } from './overlay.mjs';

// The point count is published before startup has finished: the sample's
// published default session (`startup.step('session', …)` in the web app's
// `main.js`) is applied *after* it, and that restore sets the navigation mode,
// the open accordions and the camera. A step taken in between is undone when
// the restore lands — `core-navigation-freefly` chose Free-fly and then timed
// out waiting for its controls in one run of three, and the annotation panel
// crops resolved to an empty region because Session and Visualization were
// re-opened above them. The root element's `data-startup-phase` reaches
// `ready` only after that restore, so waiting for it is what makes "ready"
// mean ready.
const PANCREAS_READY = {
  run: async page => {
    await page.waitForFunction(() => (
      document.querySelector('#filter-count')?.textContent
        ?.includes('Showing all 3,696 cells') === true
    ));
    await page.waitForFunction(
      () => document.documentElement.dataset.startupPhase === 'ready'
    );
    // The Dataset box counts what is read while the first views are read
    // ("Reading… 81 KB read in 25 requests") and settles 1.5 s after the last
    // bytes; a frame taken before that shows a count that depends on timing.
    // A scenario that holds the page's long timers (`holdLongTimersInPage`)
    // also holds that 1.5 s window, so its line stays as it was; such a
    // scenario does not show the Dataset box.
    await page.waitForFunction(() => {
      if (window.__captureHoldsLongTimers === true) return true;
      const line = document.getElementById('dataset-read');
      return line === null || line.hidden || !/Reading/.test(line.textContent);
    }, undefined, { timeout: 120_000 });
  },
};

// Every notification raised so far has expired and left the DOM.
//
// Loading the sample raises its own — the restored default session, the
// opening of its AnnData store, every field that gets coloured — and each expires on its
// own four- or five-second timer. Left alone they do two kinds of damage: in a
// window capture they appear or not depending on how fast the host ran the
// steps, so the same scenario produces different bytes; and under a
// `#notification-center` crop they drain while the frames settle, so the
// `column-reverse` stack changes height under the crop and
// `assertCropStillFramed` fails it. Waiting for an empty center removes the
// churn instead of racing it. An ERROR carries no timer, so a scenario whose
// subject is an error must drain *before* raising it.
const NOTIFICATIONS_DRAINED = {
  run: page => page.waitForFunction(
    () => document.querySelector('#notification-center')
      ?.childElementCount === 0,
    undefined,
    { timeout: 60_000 }
  ),
};

/**
 * Stop the clock that is about to delete a timed notification.
 *
 * A SUCCESS or INFO notification dismisses itself four seconds after it is
 * shown (`DEFAULT_DURATIONS_MS` in the web app's `notification-center.js`).
 * Resolving the crop costs about a second on a loaded machine and each
 * screenshot about seven hundred milliseconds, so a capture of one races its
 * subject's lifetime, and the outcome depends on the capture host's load —
 * the same scenario captured cleanly at load average 4 and failed at load
 * average 15 with no code changed.
 *
 * So the deadline is removed rather than raced: timers of a second or longer
 * scheduled from this point on are pushed past the end of the capture. Nothing
 * is dropped and no callback is skipped; they simply do not come due while the
 * shutter is open. Shorter timers — yields, the 300 ms exit removal — run
 * normally, so the work that raises the notification is not stalled. It is the
 * same kind of intervention as the seeded `Math.random`, the pinned time zone
 * and Playwright's `animations: 'disabled'`; nothing about the notification's
 * content, wording or appearance is touched, and `assertCropStillFramed` still
 * proves the real element was in frame.
 *
 * @param {import('@playwright/test').Page} page
 */
function holdNotificationTimers(page) {
  return page.evaluate(holdLongTimersInPage);
}

/**
 * The page-side half of `holdNotificationTimers`. Serialised into the browser,
 * so it closes over nothing. A scenario whose subject is raised during
 * startup, before any step can run, installs it with `addInitScript` instead.
 */
function holdLongTimersInPage() {
  window.__captureHoldsLongTimers = true;
  const realSetTimeout = window.setTimeout;
  const CAPTURE_HORIZON_MS = 600_000;
  window.setTimeout = function (handler, delay, ...args) {
    const requested = Number(delay) || 0;
    return realSetTimeout.call(
      window,
      handler,
      requested >= 1000 ? Math.max(requested, CAPTURE_HORIZON_MS) : requested,
      ...args
    );
  };
}

/** @type {ReadonlyArray<Object>} */
export const SCENARIOS = Object.freeze([
  {
    id: 'orientation-window-overview',
    topic: '_proof',
    name: 'window-pancreas-overview',
    describes:
      'The Cellucid web app with the Pancreatic endocrinogenesis sample ' +
      'loaded in its 3D Orbit view, coloured by the clusters field, with the ' +
      'sidebar and legend visible.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { waitForStable: '#categorical-field' },
      { select: '#categorical-field', label: 'clusters' },
      { waitForText: '#stats', text: 'Field: clusters (category)' },
      { waitFor: '#legend' },
    ],
    shot: { mode: 'window' },
  },

  {
    id: 'loading-panel-dataset-picker',
    topic: '_proof',
    name: 'panel-session-dataset-picker',
    describes:
      'The Session panel showing the Dataset box for the loaded Pancreatic ' +
      'endocrinogenesis sample above the Sample dataset dropdown, with the ' +
      'pointer on the dropdown.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { openSection: '#session-section' },
      { scrollIntoView: '#dataset-select' },
    ],
    shot: { mode: 'element', selector: '#session-section', padding: 6 },
    cursor: { on: '#dataset-select', at: 'right', dx: -16, state: 'hover' },
  },

  {
    id: 'coloring-panel-categorical-field',
    topic: '_proof',
    name: 'panel-coloring-categorical',
    describes:
      'The Coloring and Filtering panel with the clusters field selected as ' +
      'the categorical colouring source.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { openSection: '#coloring-filtering-section' },
      { waitForStable: '#categorical-field' },
      { select: '#categorical-field', label: 'clusters' },
      { waitForText: '#stats', text: 'Field: clusters (category)' },
      { scrollIntoView: '#categorical-field' },
    ],
    shot: { mode: 'element', selector: '#coloring-filtering-section', padding: 6 },
  },

  {
    id: 'coloring-window-annotated-field-selector',
    topic: '_proof',
    name: 'window-annotated-color-by',
    describes:
      'The Cellucid web app with the categorical field selector ringed and ' +
      'the pointer pressing it, showing where the colouring field is chosen.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { openSection: '#coloring-filtering-section' },
      { waitForStable: '#categorical-field' },
      { select: '#categorical-field', label: 'clusters' },
      { waitForText: '#stats', text: 'Field: clusters (category)' },
      { scrollIntoView: '#categorical-field' },
    ],
    shot: { mode: 'window' },
    cursor: { on: '#categorical-field', at: 'right', dx: -20, state: 'press' },
    rings: [{ on: '#categorical-field-row', label: '1' }],
  },

  {
    id: 'legend-region-categories',
    topic: '_proof',
    name: 'region-legend-categories',
    describes:
      'The legend listing every category of the selected clusters field ' +
      'beside the reported visible-point count.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { waitForStable: '#categorical-field' },
      { select: '#categorical-field', label: 'clusters' },
      { waitForText: '#stats', text: 'Field: clusters (category)' },
      { waitFor: '#legend' },
      { scrollIntoView: '#legend' },
    ],
    shot: { mode: 'region', selectors: ['#legend', '#filter-count'], padding: 10 },
  },

  {
    id: 'orientation-welcome-first-run',
    topic: '_proof',
    name: 'window-welcome-first-run',
    describes:
      'The Cellucid welcome overlay as it appears on first load, over the ' +
      'sample dataset behind it.',
    dataset: 'pancreas',
    keepWelcome: true,
    steps: [PANCREAS_READY],
    shot: { mode: 'window' },
  },

  // --- benchmarking_performance -------------------------------------------
  // The Performance Benchmark accordion is the only interface the benchmarking
  // section documents. Opening it starts live performance monitoring
  // (`activateBenchmarkingPanel` in the web app's `main.js`), which un-hides
  // `#benchmark-stats` and fills it with the FPS, frame time and Timing
  // Details min / p95 / max of whatever is loaded — here, the docs page's own
  // render loop on the capture host. That readout is a measurement of the
  // machine, not of Cellucid (see the note under this scenario), and it also
  // grows the panel about a second after it opens, which used to push Copy
  // Situation Report and Analyze Performance out of the published frame. So
  // the scenario hides the readout with a style rule rather than racing it:
  // the image shows every control and no number. The metric names live in a
  // table in `n_benchmarking_performance/05_benchmark_tools.md`. Collapsing
  // the six sections that are open on first load is what keeps the accordion
  // inside the 1000 px viewport, so the element crop is never silently
  // clamped.

  {
    id: 'benchmark-panel-controls',
    topic: 'benchmarking_performance',
    name: 'benchmark-panel-controls',
    describes:
      'The Benchmark panel before any run: Synthetic data with six presets '
      + 'from 500K to 40M, Point count, Pattern and Embeddings, the Load '
      + 'synthetic data button, and Diagnostics with Analyze performance and '
      + 'Copy situation report; the live performance readout is omitted.',
    dataset: 'pancreas',
    steps: [
      PANCREAS_READY,
      { closeSection: '#session-section' },
      { closeSection: '#visualization-section' },
      { closeSection: '#compare-views-section' },
      { closeSection: '#coloring-filtering-section' },
      { closeSection: '#highlighted-cells-section' },
      { closeSection: '#page-analysis-section' },
      { openSection: '#benchmark-section' },
      // `!important` outranks the inline `display: block` the monitor sets.
      {
        run: page => page.addStyleTag({
          content: '#benchmark-stats { display: none !important; }',
        }),
      },
      { waitFor: '#bottleneck-analyze-btn' },
      { scrollIntoView: '#benchmark-section' },
    ],
    shot: { mode: 'element', selector: '#benchmark-section', padding: 6 },
  },

  // There is deliberately NO scenario for the metric grid after a run. It was
  // written, captured, and withdrawn: `.benchmark-stats` shows a live FPS and
  // frame time, so the image carries a measurement, and a measurement taken in
  // headless Chromium on a shared machine is a number about the capture host,
  // not about Cellucid. The trial run emitted "FPS 2 / 547.37 ms" at 100,000
  // points. Publishing that would read as a performance claim. The metric
  // *names* are documented as a table in
  // `n_benchmarking_performance/05_benchmark_tools.md`, which is verifiable
  // against `index.html` and does not rot with the host's load average.

  // --- developer: attempted, withdrawn ---------------------------------------
  // `p_developer_docs/08_ui_modules_map.md` is the one page in that section
  // where a picture of the running app would beat a diagram: its tables map
  // module file -> what it owns, and a reader's real question is the inverse
  // ("I clicked this, which file is it?").
  //
  // A scenario for it was written and run, and it does not work with this tool
  // as it stands. Two independent blockers, both recorded here so the next
  // person does not rediscover them:
  //
  //   1. HEIGHT. The sidebar with its accordions expanded is several thousand
  //      CSS pixels tall. VIEWPORT is fixed at 1000 and `resolveClip` clamps
  //      into it, so the crop silently loses everything below the fold — the
  //      run emitted 288x988 css and captured two of eleven sections.
  //   2. CLOSED STATE DOES NOT HOLD. Driving every section shut with
  //      `closeSection` and capturing produced an image with Session and
  //      Visualization still open, so something re-opens them after the steps
  //      run. That something is the sample's published default session, which
  //      was applied after the point count — the old ready signal — and so
  //      after the steps; `PANCREAS_READY` now waits for it, which removes
  //      this blocker. Blocker 1 stands.
  //
  // A third, smaller point: ring labels render as a fixed 22 px circular badge,
  // so file paths cannot go on the image; they would need a numbered key on the
  // page. That part is fine — the key is more readable and it translates.
  //
  // A scenario may now name a taller `viewport`, which lifts blocker 1; the
  // figure has not been made because the page reads better as a table: it
  // carries corrected prose tables instead, which is the defect that
  // actually mattered there: it documented 19 of 29 modules and named the wrong
  // initializer for three of them.

  // ==========================================================================
  // a_orientation · c_core_interactions · d_fields_coloring_legends ·
  // e_filtering
  //
  // All on `pancreas`: 3,696 cells, 3,753 genes, categorical
  // `clusters_coarse` / `clusters` / `cell_type`, continuous `S_score` /
  // `G2M_score`, and `X_pca` as its latent space for the outlier filter. It is the
  // smallest published sample that carries all four filter surfaces, so a
  // reader can reproduce any of these frames in under a minute.
  //
  // The obs selects carry the field's index in `obsData.fields` as its option
  // value — 0 S_score, 1 G2M_score, 2 clusters_coarse, 3 clusters,
  // 4 cell_type — and the legend sorts categories alphabetically, so
  // `clusters` renders as Alpha, Beta, Delta, Ductal, Epsilon, Ngn3 high EP,
  // Ngn3 low EP, Pre-endocrine.
  ...buildFilteringScenarios(),

  // Community annotation and sessions. Appended as one builder for the same
  // reason as the block above: it can be reviewed and re-run without touching a
  // line another unit owns. The annotation scenarios run entirely against the
  // synthetic Worker in `community-annotation-mock.mjs` — no GitHub account is
  // used and no repository is written to.
  ...buildAnnotationAndSessionScenarios(),

  // The analysis panel. Same reason for the separate builder: the five figures
  // here share one long prelude — a dataset, a categorical field, and one
  // confirmed highlight page — and nothing in it is owned by a sibling unit.
  ...buildAnalysisScenarios(),

  // The local-file refusals on b_data_loading. Separate builder for the same
  // reason as the blocks above, and because these are the only scenarios that
  // hand the application a file from disk rather than driving its controls.
  ...buildDataLoadingRefusalScenarios(),

  // The Session panel on b_data_loading: the Dataset box, Open… and the Link
  // or server bar. Separate builder for the same reason as the blocks above;
  // every link points at a loopback host of this run, never the internet, so
  // a capture is reproducible offline.
  ...buildDataLinkScenarios(),

  // Opening every format, the opening itself (progress, a choice, a
  // failure), the Google Drive and OneDrive sign-in lines, and the viewers
  // `cellucid serve` and `cellucid_view()` open. Separate builder because these
  // are the only scenarios that build their own files and start servers.
  ...buildFormatAndViewerScenarios(),

  // Orientation, 2D navigation, connectivity, the camera path, the vector
  // field overlay, figure export and the benchmark: figures earlier sweeps
  // captured by hand, now each a scenario.
  ...buildViewAndExportScenarios(),

  // f_highlighting_selection. Separate builder for the same reason as the
  // blocks above: every figure there is one stage of the same few real
  // gestures on the same 2D view, so they share one prelude and one set of
  // fixed gesture coordinates that nothing else in this file uses.
  ...buildHighlightingScenarios(),

  // Multiview (Grid compare). Separate builder for the same reason as the
  // blocks above: every figure here is the same two-panel grid — the live view
  // plus one kept snapshot under locked cameras — and they share one prelude,
  // one zoom gesture and one layout assertion that nothing else uses.
  ...buildMultiviewScenarios(),

  // The documentation home page. Separate builder because these are the only
  // scenarios that emit WebP at named widths, in a light and a dark variant of
  // each scene, for a page that places them with its own `<img>` markup.
  ...buildHomeScenarios(),

  // The phone layout (a_orientation/02_system_requirements): a 390x844 touch
  // screen, at or under the 720 px width where the sidebar becomes a
  // full-screen sheet. The only scenarios with their own viewport.
  ...buildPhoneScenarios(),
]);

// ---------------------------------------------------------------------------
// Orientation, core interactions, fields/legends and filtering.
//
// Kept in one hoisted builder rather than inline so this block can be appended
// to, and reviewed, without touching a line any sibling unit owns.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildFilteringScenarios() {
  const READY = PANCREAS_READY;

  // The field selectors are rebuilt when the obs field list arrives,
  // which is after the point count settles; selecting before that is silently
  // discarded. `waitForStable` is what makes these two lines truthful.
  const clusters = [
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label: 'clusters' },
    { waitForText: '#stats', text: 'Field: clusters (category)' },
    { waitFor: '#legend' },
  ];
  /**
   * Navigation and Compare Views open, the sidebar scrolled so Navigation's
   * Mode select and the Shows picker below it are both in the window.
   */
  const SHOWS_PICKER_IN_VIEW = [
    { openSection: '#navigation-section' },
    { openSection: '#compare-views-section' },
    { waitFor: '#representation-controls:not([hidden])' },
    {
      run: page => page.locator('#navigation-section').evaluate(
        element => element.scrollIntoView({ block: 'start' })
      ),
    },
  ];
  /**
   * Show one representation in the selected view with the Shows picker, and
   * wait for the view's navigation mode to follow it.
   */
  // The picker is rebuilt when the view takes the representation, so the
  // choice is asserted where it lands — the picker reading it and the view's
  // navigation mode following it — rather than at the instant of the click.
  const show = (label, mode) => [
    { waitForStable: '#view-representation-select' },
    { run: async page => {
      const [value] = await page.selectOption('#view-representation-select', { label });
      await page.waitForFunction(
        ([expected, navigation]) => (
          document.querySelector('#view-representation-select')?.value === expected
            && document.querySelector('#navigation-mode')?.value === navigation
        ),
        [value, mode],
      );
    } },
    { wait: 2000 },
  ];
  const sScore = [
    { waitForStable: '#continuous-field' },
    { select: '#continuous-field', label: 'S_score' },
    { waitForText: '#stats', text: 'Field: S_score (continuous)' },
    { waitFor: '#legend' },
  ];
  const ins1 = [
    { waitForStable: '#gene-expression-search' },
    { fill: '#gene-expression-search', value: 'Ins1' },
    { waitFor: '#gene-expression-dropdown .dropdown-item' },
    {
      run: page =>
        page.locator('#gene-expression-dropdown .dropdown-item').first().click(),
    },
    { waitForText: '#stats', text: 'Field: Ins1 (continuous)' },
    { waitFor: '#legend' },
  ];

  /** Uncheck one categorical legend row by its visible label. */
  const hide = label => ({
    run: page => page.uncheck(`#legend input[aria-label="Show category ${label}"]`),
  });

  /** Move a range input to an exact value and let the app commit it. */
  const setRange = (selector, value, nth = 0) => ({
    run: async page => {
      const input = page.locator(selector).nth(nth);
      await input.fill(String(value));
      await input.dispatchEvent('change');
      await page.waitForTimeout(500);
    },
  });

  const MIN_SLIDER = '#legend .legend-filter input[type="range"]';
  const settle = { attempts: 20, gapMs: 200 };

  return [
    // ---------------------------------------------------------------- filtering
    {
      id: 'filtering-panel-unfiltered',
      topic: 'filtering',
      name: 'panel-clusters-unfiltered',
      describes:
        'The Coloring & filtering panel with the clusters field selected: Color'
        + ' by with its three field selectors, Show cell info, the Legend group '
        + 'with Show centroids, Show labels and the Outlier filter at 100 '
        + 'percent, the Show all and Hide all buttons above the eight-category '
        + 'legend, and Filters reading Showing all 3,696 cells above No filters '
        + 'active.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#coloring-filtering-section' }],
      shot: { mode: 'element', selector: '#coloring-filtering-section', padding: 6 },
      settle,
    },
    {
      id: 'filtering-active-filters-empty',
      topic: 'filtering',
      name: 'active-filters-empty',
      describes:
        'The Filters block with its (i), reading Showing all 3,696 cells above '
        + 'the words No filters active.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#active-filters-container' }],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      settle,
    },
    {
      id: 'filtering-scope-tooltip',
      topic: 'filtering',
      name: 'active-filters-scope-tooltip',
      describes:
        'The (i) of the Filters block pressed, and the information box it '
        + 'opens: every filter hiding cells in the selected view, from the '
        + 'legend, the range sliders and the outlier filter, that unticking one '
        + 'pauses it and removing it drops it, that filters belong to the view '
        + 'they were set in, and Learn more links.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#active-filters-container' },
        { click: '#active-filter-scope-btn' },
        { waitFor: '#active-filter-scope-tooltip' },
      ],
      shot: {
        mode: 'region',
        selectors: ['#active-filters-container', '#active-filter-scope-tooltip'],
        padding: 8,
      },
      cursor: { on: '#active-filter-scope-btn', at: 'center', state: 'press' },
      settle,
    },
    {
      id: 'filtering-window-all-visible',
      topic: 'filtering',
      name: 'window-clusters-all-visible',
      describes:
        'The whole application window with the Pancreas embedding coloured by '
        + 'clusters, every legend row checked, and the Filters block reading '
        + 'Showing all 3,696 cells above No filters active.',
      dataset: 'pancreas',
      steps: [READY, { closeSection: '#session-section' }, ...clusters, { scrollIntoView: '#active-filters-container' }, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-categories-hidden',
      topic: 'filtering',
      name: 'window-clusters-three-hidden',
      describes:
        'The same window after unchecking Ductal, Ngn3 high EP and Ngn3 low EP:'
        + ' those points are absent from the canvas, their legend rows are '
        + 'unchecked and read 0 of their cell counts, and the Filters block reads'
        + ' Showing 1,876 of 3,696 cells above one clusters row.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Ngn3 low EP'),
        hide('Ngn3 high EP'),
        { scrollIntoView: '#active-filters-container' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-active-filters-one-row',
      topic: 'filtering',
      name: 'active-filters-one-row',
      describes:
        'The Filters block reading Showing 1,876 of 3,696 cells above one row: '
        + 'a ticked checkbox, the text clusters: hiding Ductal, cut off with an '
        + 'ellipsis, and a remove cross at its right.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Ngn3 low EP'),
        hide('Ngn3 high EP'),
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      settle,
    },
    {
      id: 'filtering-type-category-visibility',
      topic: 'filtering',
      name: 'type-category-visibility',
      describes:
        'The categorical legend for clusters: Show all and Hide all above eight'
        + ' rows, three of them (Ductal, Ngn3 high EP, Ngn3 low EP) unchecked and'
        + ' reading 0 of their cell counts, each row carrying a colour swatch, a '
        + 'label, a cell count, a rename pencil and a delete bin.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Ngn3 low EP'),
        hide('Ngn3 high EP'),
        { scrollIntoView: '#legend' },
      ],
      shot: { mode: 'element', selector: '#legend', padding: 6 },
      settle,
    },
    {
      id: 'filtering-type-numeric-range',
      topic: 'filtering',
      name: 'type-numeric-range',
      describes:
        'The Visible range block of the continuous S_score legend with Live '
        + 'filtering switched off: Minimum and Maximum sliders with their values,'
        + ' and the Filter button beside Reset enabled, with the pointer on '
        + 'Filter.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        { click: '#legend .legend-filter button.switch' },
        { waitFor: '#legend .legend-filter button.switch[aria-checked="false"]' },
        { waitFor: '#legend .legend-filter .btn-primary:not([disabled])' },
      ],
      shot: { mode: 'element', selector: '#legend .legend-filter', padding: 6 },
      cursor: { on: '#legend .legend-filter .btn-primary', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'filtering-button-disabled-live-on',
      topic: 'filtering',
      name: 'filter-button-disabled-live-on',
      describes:
        'The same Visible range block with Live filtering on: the Filter button'
        + ' is greyed out under the pointer, beside an enabled Reset button.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        { waitFor: '#legend .legend-filter button.switch[aria-checked="true"]' },
        { waitFor: '#legend .legend-filter .btn-primary[disabled]' },
      ],
      shot: { mode: 'element', selector: '#legend .legend-filter', padding: 6 },
      cursor: { on: '#legend .legend-filter .btn-primary', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'filtering-window-continuous-range',
      topic: 'filtering',
      name: 'window-continuous-range-applied',
      describes:
        'The whole window coloured by S_score with the Minimum slider raised: '
        + 'only the higher-scoring cells are drawn, and the Filters block reads a'
        + ' reduced count above one S_score range row.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...sScore,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 20),
        { scrollIntoView: '#active-filters-container' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-outlier-slider-95',
      topic: 'filtering',
      name: 'outlier-slider-95',
      describes:
        'The Outlier filter with its (i), its slider moved left of the end stop'
        + ' and its readout showing 95 percent, with the pointer on the slider.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#outlier-filter-container' },
        setRange('#outlier-filter', 95),
      ],
      shot: { mode: 'element', selector: '#outlier-filter-container', padding: 8 },
      cursor: { on: '#outlier-filter', at: 'right', dx: -40, state: 'hover' },
      settle,
    },
    {
      id: 'filtering-window-outlier-100',
      topic: 'filtering',
      name: 'window-outlier-100-off',
      describes:
        'The whole window coloured by clusters with the Outlier filter at its '
        + 'default 100 percent: scattered cells sit between and around the '
        + 'clusters, and the Filters block reads Showing all 3,696 cells.',
      dataset: 'pancreas',
      steps: [READY, { closeSection: '#session-section' }, ...clusters, { scrollIntoView: '#outlier-filter-container' }, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-outlier-95',
      topic: 'filtering',
      name: 'window-outlier-95-applied',
      describes:
        'The same window with the Outlier filter at 95 percent: fewer scattered'
        + ' cells between the clusters, every legend count reduced, and a '
        + 'clusters: outlier ≤ 95% row in the Filters block above Showing 3,506 '
        + 'of 3,696 cells.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...clusters,
        { scrollIntoView: '#outlier-filter-container' },
        setRange('#outlier-filter', 95),
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      cursor: { on: '#outlier-filter', at: 'right', dx: -40, state: 'hover' },
      settle,
    },
    {
      id: 'filtering-active-filters-three-rows',
      topic: 'filtering',
      name: 'active-filters-three-rows',
      describes:
        'The Filters block holding three rows of different kinds at once: an '
        + 'S_score range row ending in its cell count, cut off with an ellipsis, '
        + 'a clusters hiding row, and a clusters: outlier ≤ 95% row with no '
        + 'count.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 15),
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Epsilon'),
        { scrollIntoView: '#outlier-filter-container' },
        setRange('#outlier-filter', 95),
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      settle,
    },
    {
      id: 'filtering-active-filters-row-disabled',
      topic: 'filtering',
      name: 'active-filters-row-disabled',
      describes:
        'The same three-row Filters list with the clusters hiding row unchecked'
        + ' under the pointer: that row is greyed and struck through while the '
        + 'other two stay black, and the count line above has risen.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 15),
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Epsilon'),
        { scrollIntoView: '#outlier-filter-container' },
        setRange('#outlier-filter', 95),
        { scrollIntoView: '#active-filters-container' },
        {
          run: page =>
            page.uncheck(
              '#active-filters [data-filter-id="obs-category-3"]  input[type="checkbox"]'
            ),
        },
        { waitFor: '#active-filters [data-filter-id="obs-category-3"][data-disabled="true"]' },
      ],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      cursor: {
        on: '#active-filters [data-filter-id="obs-category-3"] input[type="checkbox"]',
        at: 'center',
        state: 'hover',
      },
      settle,
    },
    {
      id: 'filtering-active-filters-remove',
      topic: 'filtering',
      name: 'active-filters-remove-button',
      describes:
        'A one-row Filters list with the pointer resting on the remove cross at'
        + ' the right edge of the row.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#legend' },
        hide('Ductal'),
        hide('Epsilon'),
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      cursor: {
        on: '#active-filters [data-filter-id="obs-category-3"] .list-item-remove',
        at: 'center',
        state: 'hover',
      },
      settle,
    },
    {
      id: 'filtering-window-zero-visible',
      topic: 'filtering',
      name: 'window-zero-visible-cells',
      describes:
        'The whole window with a completely empty canvas after pressing Hide '
        + 'all in the clusters legend: every row unchecked and reading 0 of its '
        + 'count, and the Filters block reading Showing 0 of 3,696 cells above '
        + 'one clusters hiding row.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...clusters,
        { scrollIntoView: '#legend' },
        { run: page => page.locator('#legend button', { hasText: 'Hide All' }).click() },
        { waitForText: '#filter-count', text: 'Showing 0 of 3,696 cells' },
        { scrollIntoView: '#active-filters-container' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-legend-category-unavailable',
      topic: 'filtering',
      name: 'legend-category-unavailable',
      describes:
        'The clusters legend under a narrow S_score range filter: Alpha, '
        + 'Ductal, Ngn3 high EP and Ngn3 low EP stay ticked with small counts, '
        + 'while Beta, Delta, Epsilon and Pre-endocrine are greyed out with '
        + 'disabled checkboxes and read 0 cells.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 62),
        ...clusters,
        { scrollIntoView: '#legend' },
      ],
      shot: { mode: 'element', selector: '#legend', padding: 6 },
      settle,
    },
    {
      id: 'filtering-window-gene-range-active',
      topic: 'filtering',
      name: 'window-gene-range-active',
      describes:
        'The whole window coloured by the gene Ins1 with the Minimum slider '
        + 'raised, so only the cells expressing it are drawn and an Ins1 range '
        + 'row is listed in the Filters block.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...ins1,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 3),
        { scrollIntoView: '#active-filters-container' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-gene-range-scope-lost',
      topic: 'filtering',
      name: 'window-gene-range-scope-lost',
      describes:
        'The same window after switching the active field to clusters: every '
        + 'cell is drawn again, the Ins1 row has left the Filters block, and the '
        + 'count line reads Showing all 3,696 cells.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...ins1,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 3),
        ...clusters,
        { waitForText: '#filter-count', text: 'Showing all 3,696 cells' },
        { scrollIntoView: '#active-filters-container' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },

    // ------------------------------------------------- fields, coloring, legends
    {
      id: 'fields-selector-rows',
      topic: 'web_app',
      name: 'field-selectors-three-routes',
      describes:
        'The Color by heading with its (i), above the three field selectors '
        + 'stacked in order — Categorical obs holding clusters, Continuous obs '
        + 'holding None, and an empty Gene expression search box — each followed '
        + 'by four small icon buttons, with the pointer pressing the first.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#color-by-controls' }],
      shot: {
        mode: 'region',
        selectors: ['#color-by-controls .control-group-header', '#gene-expression-row'],
        padding: 10,
      },
      cursor: { on: '#categorical-field', at: 'right', dx: -18, state: 'press' },
      settle,
    },
    {
      id: 'fields-legend-categorical',
      topic: 'web_app',
      name: 'legend-categorical-clusters',
      describes:
        'The categorical legend for clusters: the Show all and Hide all buttons'
        + ' above eight rows, each holding a checkbox, a colour swatch, a label, '
        + 'a cell count, a rename pencil and a delete bin.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#legend' }],
      shot: { mode: 'element', selector: '#legend', padding: 6 },
      settle,
    },
    {
      id: 'fields-legend-continuous',
      topic: 'web_app',
      name: 'legend-continuous-sscore',
      describes:
        'The continuous legend for S_score: the Color scale (Viridis) bar with '
        + 'its two end values, Log color scale off and Rescale to range on, each '
        + 'with its (i), and below a rule the Visible range block with Live '
        + 'filtering on, Minimum and Maximum sliders, and Filter and Reset.',
      dataset: 'pancreas',
      steps: [READY, ...sScore, { scrollIntoView: '#legend' }],
      shot: { mode: 'element', selector: '#legend', padding: 6 },
      settle,
    },
    {
      id: 'fields-colormap-menu',
      topic: 'web_app',
      name: 'legend-colormap-menu-open',
      describes:
        'The colormap list open under the S_score colour bar: seven named '
        + 'gradients, Viridis, Plasma, Inferno, Magma, Cividis, Turbo and '
        + 'Twilight, with Viridis marked as the active one.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        { run: page => page.locator('#legend .colorbar-gradient').click() },
        { waitFor: '#legend .colormap-option' },
        { scrollIntoView: '#legend' },
      ],
      shot: { mode: 'element', selector: '#legend', padding: 6 },
      settle,
    },
    {
      id: 'fields-gene-search-open',
      topic: 'web_app',
      name: 'gene-search-dropdown-open',
      describes:
        'The Gene expression search box holding the letters Rp with its list '
        + 'open below, naming the genes that contain them.',
      dataset: 'pancreas',
      steps: [
        READY,
        { waitForStable: '#gene-expression-search' },
        { scrollIntoView: '#gene-expression-row' },
        { fill: '#gene-expression-search', value: 'Rp' },
        { waitFor: '#gene-expression-dropdown .dropdown-item' },
      ],
      shot: {
        mode: 'region',
        selectors: ['label[for="gene-expression-search"]', '#gene-expression-row', '#gene-expression-dropdown'],
        padding: 8,
      },
      settle,
    },
    {
      id: 'fields-gene-search-no-match',
      topic: 'web_app',
      name: 'gene-search-no-match',
      describes:
        'The Gene expression list after typing a symbol this sample does not '
        + 'publish: No gene matches "CD19", a sentence giving the 3,753 gene '
        + 'names the dataset publishes, and a link Why a gene may be missing.',
      dataset: 'pancreas',
      steps: [
        READY,
        { waitForStable: '#gene-expression-search' },
        { scrollIntoView: '#gene-expression-row' },
        { fill: '#gene-expression-search', value: 'CD19' },
        { waitForText: '#gene-expression-dropdown', text: 'No gene matches' },
      ],
      shot: {
        mode: 'region',
        selectors: ['label[for="gene-expression-search"]', '#gene-expression-row', '#gene-expression-dropdown'],
        padding: 8,
      },
      settle,
    },
    {
      id: 'fields-window-gene-coloring',
      topic: 'web_app',
      name: 'window-color-by-gene-ins1',
      describes:
        'The whole window coloured by expression of the gene Ins1 on a Viridis '
        + 'scale, dark where the gene is absent and bright in one lobe of the '
        + 'embedding, with the continuous legend in the sidebar.',
      dataset: 'pancreas',
      steps: [READY, { closeSection: '#session-section' }, ...ins1, { scrollIntoView: '#legend' }, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'fields-window-centroids',
      topic: 'web_app',
      name: 'window-centroids-and-labels',
      describes:
        'The Pancreas embedding coloured by clusters with Show centroids and '
        + 'Show labels both ticked, so each cluster carries a marker at its '
        + 'centre with its category name printed on it.',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        ...clusters,
        { scrollIntoView: '#centroid-controls' },
        { run: page => page.check('#toggle-centroid-points') },
        { run: page => page.check('#toggle-centroid-labels') },
        { wait: 1200 },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },

    // ------------------------------------------------------- core interactions
    {
      id: 'core-navigation-orbit',
      topic: 'web_app',
      name: 'navigation-controls-orbit',
      describes:
        'The Navigation block with Mode set to Orbit and the pointer pressing '
        + 'the select: a Keyboard speed slider, Google Earth-style drag and Show '
        + 'orbit anchor.',
      dataset: 'pancreas',
      steps: [
        READY,
        // Alone, so the whole block fits above the sidebar's bottom edge.
        onlySections('#navigation-section'),
        { scrollIntoView: '#navigation-controls' },
        { select: '#navigation-mode', value: 'orbit' },
        { waitFor: '#orbit-controls' },
      ],
      shot: { mode: 'element', selector: '#navigation-controls', padding: 6 },
      cursor: { on: '#navigation-mode', at: 'right', dx: -18, state: 'press' },
      settle,
    },
    {
      id: 'core-navigation-planar',
      topic: 'web_app',
      name: 'navigation-controls-planar',
      describes:
        'The Navigation block with Mode set to Planar: a Keyboard speed slider,'
        + ' Zoom to cursor and Invert axes.',
      dataset: 'pancreas',
      steps: [
        READY,
        // Alone, so the whole block fits above the sidebar's bottom edge.
        onlySections('#navigation-section'),
        { scrollIntoView: '#navigation-controls' },
        { select: '#navigation-mode', value: 'planar' },
        { waitFor: '#planar-controls' },
      ],
      shot: { mode: 'element', selector: '#navigation-controls', padding: 6 },
      settle,
    },
    {
      id: 'core-navigation-freefly',
      topic: 'web_app',
      name: 'navigation-controls-freefly',
      describes:
        'The Navigation block with Mode set to Free-fly: Look sensitivity and '
        + 'Move speed sliders with their readouts, Invert look axes and '
        + 'Projectile shooting, and Capture pointer with its (i).',
      dataset: 'pancreas',
      steps: [
        READY,
        // Alone, so the whole block fits above the sidebar's bottom edge.
        onlySections('#navigation-section'),
        { scrollIntoView: '#navigation-controls' },
        { select: '#navigation-mode', value: 'free' },
        { waitFor: '#freefly-controls' },
      ],
      shot: { mode: 'element', selector: '#navigation-controls', padding: 6 },
      settle,
    },
    {
      id: 'core-window-dimension-3d',
      topic: 'web_app',
      name: 'window-dimension-3d',
      describes:
        'The whole window showing the Pancreas 3D embedding on the Studio '
        + 'backdrop, with the Navigation Mode select reading Orbit and, under '
        + 'Compare views, the Shows picker reading UMAP · 3D.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, ...SHOWS_PICKER_IN_VIEW, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'core-window-dimension-2d',
      topic: 'web_app',
      name: 'window-dimension-2d',
      describes:
        'The same window after choosing UMAP · 2D in the Shows picker, with the'
        + ' pointer on it: the cells now lie on a plane seen face-on over a '
        + 'dotted backdrop, and the Navigation Mode select has followed the view '
        + 'to Planar.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        ...SHOWS_PICKER_IN_VIEW,
        ...show('UMAP · 2D', 'planar'),
        ...SHOWS_PICKER_IN_VIEW,
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      cursor: { on: '#view-representation-select', at: 'right', dx: -18, state: 'press' },
      settle,
    },
    {
      id: 'core-window-dimension-1d',
      topic: 'web_app',
      name: 'window-dimension-1d',
      describes:
        'The same window after choosing UMAP · 1D in the Shows picker: the '
        + 'cells are spread along one axis, and the Navigation Mode select '
        + 'has followed the view to Planar.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        ...SHOWS_PICKER_IN_VIEW,
        ...show('UMAP · 1D', 'planar'),
        ...SHOWS_PICKER_IN_VIEW,
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'core-render-mode-select',
      topic: 'web_app',
      name: 'render-mode-select',
      describes:
        'The Render mode select reading Points, ringed and with the pointer '
        + 'pressing it, above Shader quality, Point size and the Depth cues '
        + 'sliders (3D lighting, Atmospheric fog, Perspective size scaling).',
      dataset: 'pancreas',
      steps: [READY, { openSection: '#visualization-section' }, { scrollIntoView: '#render-mode' }],
      shot: {
        mode: 'region',
        selectors: ['label[for="render-mode"]', '#depth-controls'],
        padding: 8,
      },
      cursor: { on: '#render-mode', at: 'right', dx: -18, state: 'press' },
      settle,
    },
    {
      // Smoke animates by default, so a capture of it would never settle.
      // Animation speed 0 reads `Still`, and a still smoke view is not
      // redrawn at all once its noise exists (`needsContinuousFrames()` in the
      // web app's `smoke-renderer.js`): two consecutive frames are then the
      // same bytes, which is what the settle loop waits for. The ray-start
      // jitter is a fixed blue-noise tile, not a per-frame offset, so nothing
      // else moves.
      id: 'core-window-smoke-cloud',
      topic: 'web_app',
      name: 'window-smoke-cloud',
      describes:
        'The whole window in Volumetric smoke cloud mode: the Pancreas cells '
        + 'drawn as a soft cloud tinted with the colours of their clusters, and '
        + 'the sidebar showing Render mode set to Volumetric smoke cloud (beta) '
        + 'with its Beta tag above the Smoke quality, Smoke shape, Smoke lighting'
        + ' and Performance controls, Animation speed at Still.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { openSection: '#visualization-section' },
        { select: '#render-mode', value: 'smoke' },
        { waitFor: '#smoke-speed' },
        {
          run: async page => {
            const speed = page.locator('#smoke-speed');
            await speed.fill('0');
            await speed.dispatchEvent('input');
            await speed.blur();
          },
        },
        { waitForText: '#smoke-speed-display', text: 'Still' },
        { scrollIntoView: '#smoke-controls' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { attempts: 40, gapMs: 250 },
    },

    // ------------------------------------------------------------- orientation
    {
      id: 'orientation-window-map',
      topic: 'web_app',
      name: 'window-orientation-map',
      describes:
        'The whole application window with the Pancreas sample loaded, the '
        + 'sidebar showing Coloring & filtering with the embedding coloured by '
        + 'clusters, and the eight-row legend beside the matching colours on the '
        + 'canvas.',
      dataset: 'pancreas',
      steps: [READY, { closeSection: '#session-section' }, ...clusters, { scrollIntoView: '#legend' }, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },

    // ------------------------------------------------------------ accessibility
    // The Keyboard Shortcuts accordion is the app's own authoritative shortcut
    // list, and the accessibility page restates it row for row. Capturing the
    // whole `<details>` keeps its summary in frame, so the reader can see which
    // panel to open, and the six sections that are open on first load are
    // collapsed because an element crop is clamped to the 1000 px viewport.
    //
    // `#shortcuts-section` carries no `data-state-serializer-skip`, so its
    // open/closed state is part of the session UI inventory
    // (`state-serializer/ui-controls.js`) and the built-in sample's published
    // starting state re-closes it. Opening before that lands produces a capture
    // of the collapsed summary and nothing else, so the state has to be waited
    // for first and the panel checked afterwards.
    {
      id: 'accessibility-panel-keyboard-shortcuts',
      topic: 'web_app',
      name: 'keyboard-shortcuts',
      describes:
        'The Keyboard shortcuts panel open in the Cellucid sidebar, listing '
        + 'five groups, each with its (i): General (R, F, H, 1 2 3, O P G and ?),'
        + ' Navigation (W A S D and Shift), Zoom (= ] and - [), Free-fly (E Q and'
        + ' Esc) and Selection (Alt+click, Alt+drag, Shift+Alt+click, '
        + 'Ctrl+Alt+click, Esc to cancel a selection, Esc to leave Select mode, '
        + 'and X to clear highlights).',
      dataset: 'pancreas',
      steps: [
        READY,
        { closeSection: '#session-section' },
        { closeSection: '#visualization-section' },
        { closeSection: '#navigation-section' },
        { closeSection: '#compare-views-section' },
        { closeSection: '#coloring-filtering-section' },
        { closeSection: '#highlighted-cells-section' },
        { closeSection: '#page-analysis-section' },
        { openSection: '#shortcuts-section' },
        { scrollIntoView: '#shortcuts-section' },
        { run: async page => {
          const open = await page.evaluate(
            () => document.querySelector('#shortcuts-section')?.open ?? null
          );
          if (open !== true) {
            throw new Error(
              '#shortcuts-section closed again after it was opened; something '
              + 'restored the serialized accordion state after this step.'
            );
          }
        } },
      ],
      shot: { mode: 'element', selector: '#shortcuts-section', padding: 6 },
      settle,
    },
  ];
}

// ---------------------------------------------------------------------------
// Community annotation and sessions.
//
// The annotation feature is the only documented surface that signs in to a
// third party. Every scenario below reaches it through the synthetic Worker in
// `./community-annotation-mock.mjs`, installed in `setup` before navigation, so
// the screens are the product's real screens while no account is used and no
// repository is written to. The vote tallies and the consensus arithmetic in
// these images are computed by the application from the mocked annotator files
// on a real Pull — they are not drawn by this tool.
//
// Dataset: `pancreas`, annotatable column `clusters`, target category `Beta`.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildAnnotationAndSessionScenarios() {
  const READY = PANCREAS_READY;
  const CA = '#community-annotation-section';
  const PANEL = '#community-annotation-controls';
  const MODAL = '.community-annotation-modal-overlay [role="document"]';
  const DETAIL = '.community-annotation-voting-detail';

  // Collapsing the sections that are open on first load is what keeps the
  // annotation accordion inside the 1000 px viewport; an element crop is
  // clamped to the viewport, so a panel scrolled out of it crops to nothing.
  const openPanel = [
    READY,
    { closeSection: '#session-section' },
    { closeSection: '#visualization-section' },
    { closeSection: '#compare-views-section' },
    { closeSection: '#coloring-filtering-section' },
    { closeSection: '#highlighted-cells-section' },
    { closeSection: '#page-analysis-section' },
    { openSection: CA },
    { scrollIntoView: PANEL },
  ];

  /** Make `clusters` annotatable and pin the consensus rule under capture. */
  const annotate = (minAnnotators, threshold) => ({
    run: async page => {
      await page.evaluate(async ({ min, th }) => {
        const { getCommunityAnnotationSession } =
          await import('/assets/js/app/community-annotations/session.js');
        const session = getCommunityAnnotationSession();
        session.setFieldAnnotated('clusters', true);
        session.setAnnotatableConsensusSettings('clusters', {
          minAnnotators: min, threshold: th,
        });
      }, { min: minAnnotators, th: threshold });
    },
  });

  /** Pull, so the merged view is built by the product from the mocked files. */
  const pull = [
    { click: `${PANEL} button:has-text("GitHub sync")` },
    { waitFor: MODAL },
    { click: '.community-annotation-modal-overlay button:has-text("Pull latest")' },
    { waitForText: '.community-annotation-modal-overlay', text: 'Pulled latest annotations.' },
  ];

  const closeModal = [
    { click: '.community-annotation-modal-overlay button:has-text("Close")' },
    { waitForHidden: MODAL },
  ];

  const selectClusters = [
    { openSection: '#coloring-filtering-section' },
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label: '🗳️ clusters' },
    { waitFor: '#legend' },
  ];

  const votingOnBeta = [
    ...selectClusters,
    { click: '[aria-label="Vote on category Beta"]' },
    { waitFor: DETAIL },
  ];

  const mock = options => async ({ page }) => {
    const module = await import('./community-annotation-mock.mjs');
    await module.installGitHubMock({ page, ...options });
  };

  /** The four renderer controls, set to values no default holds. */
  const RENDERER_CHANGES = [
    { select: '#hp-shader-quality', value: 'light' },
    { run: async page => {
      // Adaptive LOD gates the forced-level slider, so it is enabled first.
      await page.setChecked('#hp-lod-enabled', true);
      await page.fill('#hp-lod-force', '4');
      await page.dispatchEvent('#hp-lod-force', 'input');
      await page.setChecked('#hp-frustum-culling', true);
    } },
  ];

  /**
   * Save the current state to a file, reload the application so every control
   * returns to its default, then restore that file through the real Load
   * picker. One scenario, one context: a bundle saved here is not visible to
   * any other scenario.
   */
  //
  // The reload starts the application over, so it is waited for exactly as a
  // fresh scenario is: `PANCREAS_READY` (the sample's own default session has
  // been applied — otherwise "Session fully restored." below can match *that*
  // restore instead of this one) and `NOTIFICATIONS_DRAINED`. With
  // `holdNotifications`, the restore's own notifications are kept on screen
  // for a crop of them; see `holdNotificationTimers`.
  const saveThenRestore = ({ holdNotifications = false } = {}) => ({
    run: async page => {
      const download = page.waitForEvent('download');
      await page.click('#save-state-btn');
      const file = await download;
      const saved = await file.path();
      await page.reload({ waitUntil: 'domcontentloaded' });
      const welcome = page.locator('#welcome-modal');
      await welcome.waitFor({ state: 'visible' });
      await page.keyboard.press('Escape');
      await welcome.waitFor({ state: 'hidden' });
      await PANCREAS_READY.run(page);
      await NOTIFICATIONS_DRAINED.run(page);
      if (holdNotifications) await holdNotificationTimers(page);
      const chooser = page.waitForEvent('filechooser');
      await page.click('#load-state-btn');
      await (await chooser).setFiles(saved);
      await page.waitForFunction(() => (
        document.querySelector('#notification-center')?.textContent
          ?.includes('Session fully restored.') === true
      ));
    },
  });

  const mockUsers = async which => {
    const module = await import('./community-annotation-mock.mjs');
    return module[which]();
  };

  return [
    // ------------------------------------------------ 1. the entry point
    {
      id: 'annotation-lifecycle-01-disconnected',
      topic: 'community_annotation',
      name: 'lifecycle-01-disconnected',
      describes:
        'The Community annotation panel before any repository is connected: '
        + 'Repository with its (i), the line Not connected. and a Connect repo '
        + 'button, with the pointer on that button.',
      dataset: 'pancreas',
      setup: mock({ signedIn: false, bindRepo: false }),
      steps: openPanel,
      shot: { mode: 'element', selector: PANEL, padding: 10 },
      cursor: { on: `${PANEL} button:has-text("Connect repo")`, at: 'center', state: 'hover' },
    },

    // ------------------------------------------- 2. wizard step 1, sign in
    {
      id: 'annotation-lifecycle-02-signin',
      topic: 'community_annotation',
      name: 'lifecycle-02-signin',
      describes:
        'Step 1 of 4 of the GitHub sync wizard: a header showing the dataset '
        + 'with GitHub and Repo both not connected, the four steps Sign in, '
        + 'Install app, Choose repository and Sync with Sign in current, the line '
        + '"Cellucid gets no repository access until you install its GitHub App '
        + 'on repositories you choose.", and a Continue with GitHub button.',
      dataset: 'pancreas',
      setup: mock({ signedIn: false, bindRepo: false }),
      steps: [
        ...openPanel,
        { click: `${PANEL} button:has-text("Connect repo")` },
        { waitFor: MODAL },
        { waitForText: MODAL, text: 'Step 1 of 4' },
      ],
      shot: { mode: 'element', selector: MODAL, padding: 0 },
      cursor: {
        on: '.community-annotation-modal-overlay button:has-text("Continue with GitHub")',
        at: 'center', state: 'hover',
      },
    },

    // -------------------------------------- 3. wizard step 2, install app
    {
      id: 'annotation-lifecycle-03-install-app',
      topic: 'community_annotation',
      name: 'lifecycle-03-install-app',
      describes:
        'Step 2 of 4 of the GitHub sync wizard, Install app: the header shows '
        + 'the dataset, the signed-in GitHub account and Repo not connected; '
        + 'the line "Add repo opens GitHub to choose repositories; then '
        + 'Reload.", Add repo (under the pointer) and Reload, three repository '
        + 'cards marked Public or Private, and "Showing 1–3 of 3 '
        + 'repositories. Page 1 of 1."',
      dataset: 'pancreas',
      setup: mock({ bindRepo: false }),
      steps: [
        ...openPanel,
        { click: `${PANEL} button:has-text("Connect repo"), ${PANEL} button:has-text("Choose repo"), ${PANEL} button:has-text("GitHub sync")` },
        { waitFor: MODAL },
        { waitFor: '.community-annotation-repo-card' },
      ],
      shot: { mode: 'element', selector: MODAL, padding: 0 },
      cursor: {
        on: '.community-annotation-modal-overlay button:has-text("Add repo")',
        at: 'center', state: 'hover',
      },
    },

    // ---------------------------- 4. wizard step 3, no repository matches
    {
      id: 'annotation-state-no-repository-matches-filter',
      topic: 'community_annotation',
      name: 'state-no-repository-matches-filter',
      describes:
        'Step 3 of 4 of the GitHub sync wizard, Choose repository, with a '
        + 'filter typed that matches no repository and the line No repositories '
        + 'match this filter. under it.',
      dataset: 'pancreas',
      setup: mock({ bindRepo: false, repoCount: 240 }),
      steps: [
        ...openPanel,
        { click: `${PANEL} button:has-text("Connect repo"), ${PANEL} button:has-text("Choose repo"), ${PANEL} button:has-text("GitHub sync")` },
        { waitFor: '.community-annotation-repo-card' },
        { click: '.community-annotation-wizard-next' },
        { waitFor: '[aria-label="Selectable annotation repositories"]' },
        { fill: '[aria-label="Filter repositories"]', value: 'no-such-repository' },
        { waitForText: MODAL, text: 'No repositories match this filter.' },
      ],
      shot: { mode: 'element', selector: MODAL, padding: 0 },
    },

    // ------------------------------------------ 5. wizard step 4, sync
    {
      id: 'annotation-lifecycle-05-sync',
      topic: 'community_annotation',
      name: 'lifecycle-05-sync',
      describes:
        'Step 4 of 4 of the GitHub sync wizard, Sync: the dataset, GitHub '
        + 'account and repository in the header, the line "Pull fetches the '
        + 'latest votes; Publish writes them for everyone.", Pull latest and '
        + 'Publish, an Auto pull checkbox with its interval, and the status line '
        + '"Pulled latest annotations."',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page, users: module.usersConsensus(), config: module.configDocument(),
        });
      },
      steps: [...openPanel, ...pull],
      shot: { mode: 'element', selector: MODAL, padding: 0 },
    },

    // ------------------------------------------ 6. connected status panel
    {
      id: 'annotation-lifecycle-06-connected-status',
      topic: 'community_annotation',
      name: 'lifecycle-06-connected-status',
      describes:
        'The connected Community annotation panel: GitHub sync with the '
        + 'Dataset, GitHub and Repo rows and a GitHub sync… button, a Profile '
        + 'block naming the signed-in user with Edit and Clear, and the '
        + 'collapsed Annotated columns, Consensus column and Snapshot and '
        + 'cache groups.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page, users: module.usersConsensus(), config: module.configDocument(),
        });
      },
      steps: [...openPanel, ...pull, ...closeModal],
      shot: { mode: 'element', selector: PANEL, padding: 10 },
    },

    // ------------------------------------------- 7. MANAGE ANNOTATION
    {
      id: 'annotation-lifecycle-07-manage-annotation',
      topic: 'community_annotation',
      name: 'lifecycle-07-manage-annotation',
      describes:
        'Annotated columns open in the Community annotation panel: the '
        + 'Categorical obs picker on the voting column clusters with Add, Remove '
        + 'and Close, and its consensus settings with a Threshold slider and a '
        + 'Min annotators field.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page, users: module.usersConsensus(), config: module.configDocument(),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        { click: `${PANEL} .disclosure-header:has-text("Annotated columns")` },
        // The consensus settings belong to whichever column is chosen in the
        // picker, so they do not appear until one is chosen.
        { waitForStable: `${PANEL} select[id^="community-annotation-categorical-obs"]` },
        { select: `${PANEL} select[id^="community-annotation-categorical-obs"]`, label: '\u{1F5F3}\uFE0F clusters' },
        { waitForText: PANEL, text: 'Min annotators' },
        // The open group makes the panel taller than the window, so the
        // crop is the group itself, down to the next group's header.
        { scrollIntoView: `${PANEL} .disclosure-header:has-text("Consensus column")` },
      ],
      shot: {
        mode: 'region',
        // The group's whole box (its header row carries the (i) beside the
        // header), down to the next group's header row.
        selectors: [
          `${PANEL} .disclosure:has(> .disclosure-header-row .disclosure-header:has-text("Annotated columns"))`,
          `${PANEL} .disclosure-header-row:has(.disclosure-header:has-text("Consensus column"))`,
        ],
        padding: 10,
      },
    },

    // ------------------------------------------------ 8. the votable badge
    {
      id: 'annotation-lifecycle-08-votable-badge',
      topic: 'community_annotation',
      name: 'lifecycle-08-votable-badge',
      describes:
        'The Categorical obs select holding clusters, prefixed by a ballot-box '
        + 'badge marking it open for voting, and the legend beneath it: each '
        + 'category with its consensus label (Beta cell under Beta) and its up '
        + 'and down vote counts.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page, users: module.usersConsensus(), config: module.configDocument(),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...selectClusters,
      ],
      shot: { mode: 'region', selectors: ['#color-by-controls .control-group-header', '#categorical-field', '#legend'], padding: 8 },
    },

    // -------------------------------------------------- 9. Pending
    {
      id: 'annotation-lifecycle-09-pending',
      topic: 'community_annotation',
      name: 'lifecycle-09-pending',
      describes:
        'The voting view for the Beta cluster with a banner reading "Pending '
        + '(1 voter)", one suggestion card (Beta cell, its ontology id, marker '
        + 'and evidence, one upvote), and the New suggestion form below.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page,
          users: module.usersPending(),
          config: module.configDocument({ minAnnotators: 3 }),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...votingOnBeta,
      ],
      shot: { mode: 'element', selector: DETAIL, padding: 0 },
    },

    // ------------------------------------------------- 11. Consensus
    {
      id: 'annotation-lifecycle-11-consensus',
      topic: 'community_annotation',
      name: 'lifecycle-11-consensus',
      describes:
        'The voting view for the Beta cluster with a banner reading '
        + '"Consensus: "Beta cell" (100% • 4 voters)", the suggestion card with '
        + 'its ontology id, marker, evidence and four upvotes, and the New '
        + 'suggestion form below.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page,
          users: module.usersConsensus(),
          config: module.configDocument({ minAnnotators: 3 }),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...votingOnBeta,
      ],
      shot: { mode: 'element', selector: DETAIL, padding: 0 },
    },

    // -------------------------------------------------- 12. Disputed
    {
      id: 'annotation-lifecycle-12-disputed',
      topic: 'community_annotation',
      name: 'lifecycle-12-disputed',
      describes:
        'The voting view for the Beta cluster with a banner reading '
        + '"Disputed: "Beta (immature), Beta cell" (0% • 4 voters)", two '
        + 'suggestion cards each at net 0 with two upvotes and two downvotes, '
        + 'and the New suggestion form below.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page,
          users: module.usersDisputed(),
          config: module.configDocument({ minAnnotators: 3 }),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...votingOnBeta,
      ],
      shot: { mode: 'element', selector: DETAIL, padding: 0 },
    },

    // ------------------------------------------ 10. casting a vote (hover)
    {
      id: 'annotation-lifecycle-10-vote',
      topic: 'community_annotation',
      name: 'lifecycle-10-vote',
      describes:
        'One suggestion card showing its label, net score, ontology, marker '
        + 'and evidence lines, its up and down vote buttons, the proposer '
        + 'attribution and the comment box, with the pointer on the upvote '
        + 'button.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page,
          users: module.usersConsensus(),
          config: module.configDocument({ minAnnotators: 3 }),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...votingOnBeta,
      ],
      shot: {
        mode: 'element',
        selector: '.community-annotation-suggestion-card',
        padding: 8,
      },
      cursor: {
        on: '.community-annotation-suggestion-card .community-annotation-vote-btn',
        at: 'right', dx: -8, state: 'hover',
      },
    },

    // --------------------------------------------- 13. propose a suggestion
    {
      id: 'annotation-lifecycle-13-new-suggestion',
      topic: 'community_annotation',
      name: 'lifecycle-13-new-suggestion',
      describes:
        'The New suggestion form filled in with a label, an ontology id and '
        + 'marker genes, with Search CAP, Search Ontology and Search Markers '
        + 'buttons beside the fields and Add and Clear buttons below, and the '
        + 'pointer on Add.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page,
          users: module.usersConsensus(),
          config: module.configDocument({ minAnnotators: 3 }),
        });
      },
      steps: [
        ...openPanel, ...pull, ...closeModal,
        annotate(3, 0.5),
        ...votingOnBeta,
        // The form sits at the foot of the modal's own scroller, below the
        // suggestion list, so it has to be scrolled to before it is fillable.
        { scrollIntoView: '.community-annotation-modal-overlay [placeholder="Label (required)"] >> nth=-1' },
        { fill: '.community-annotation-modal-overlay [placeholder="Label (required)"] >> nth=-1', value: 'Delta cell' },
        { fill: '.community-annotation-modal-overlay [placeholder="Ontology id (optional)"] >> nth=-1', value: 'CL:0000173' },
        { fill: '.community-annotation-modal-overlay [placeholder="Marker genes (optional)"] >> nth=-1', value: 'Sst, Rbp4' },
        // The whole form, heading included, inside the modal's scroller.
        { run: page => page.locator('.community-annotation-new').last()
          .evaluate(element => element.scrollIntoView({ block: 'center' })) },
      ],
      shot: {
        mode: 'element',
        selector: '.community-annotation-new >> nth=-1',
        // The form's heading sits above its box.
        padding: 22,
      },
      cursor: {
        on: '.community-annotation-new >> nth=-1 >> button:has-text("Add")',
        at: 'center', state: 'hover',
      },
    },

    // ---------------------------------------------------- 17. publish
    {
      id: 'annotation-lifecycle-17-publish',
      topic: 'community_annotation',
      name: 'lifecycle-17-publish',
      describes:
        'The same Sync step after Publish, its status line reading "Publish '
        + 'complete." under the Pull latest and Publish buttons.',
      dataset: 'pancreas',
      setup: async ({ page }) => {
        const module = await import('./community-annotation-mock.mjs');
        await module.installGitHubMock({
          page, users: module.usersConsensus(), config: module.configDocument(),
        });
      },
      steps: [
        ...openPanel, ...pull,
        { run: async page => {
          await page.evaluate(async () => {
            const { getCommunityAnnotationSession } =
              await import('/assets/js/app/community-annotations/session.js');
            getCommunityAnnotationSession().setProfile({
              username: 'ghid_4242', githubUserId: 4242, login: 'r-okafor',
              displayName: 'R. Okafor', title: 'Postdoc', orcid: '', linkedin: '',
            });
          });
        } },
        { click: '.community-annotation-modal-overlay button:has-text("Publish")' },
        { waitForText: '.community-annotation-modal-overlay', text: 'Publish complete.' },
      ],
      shot: { mode: 'element', selector: MODAL, padding: 0 },
    },

    // ======================================================== sessions
    // A session bundle is a file, so these scenarios drive the real file
    // controls: `Save` produces a download, and `Load` opens a file
    // chooser that Playwright answers. The whole round trip runs inside one
    // scenario, because the tool gives each scenario a fresh browser context
    // and a bundle saved in one context is not visible to another.

    {
      id: 'sessions-save-restore-01-session-controls',
      topic: 'sessions_sharing',
      name: 'save-restore-01-session-controls',
      describes:
        'The Session accordion showing the Dataset box for the loaded sample, '
        + 'the Open data group with the Sample dataset dropdown, the Open… '
        + 'button and the Link or server bar, the Session file group with Save '
        + 'and Load below them, and the collapsed Data mapping group.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
      ],
      shot: { mode: 'element', selector: '#session-section', padding: 8 },
    },

    {
      id: 'sessions-save-restore-02-save-state-action',
      topic: 'sessions_sharing',
      name: 'save-restore-02-save-state-action',
      describes:
        'The Session file group with its (i), Save and Load side by side, and '
        + 'the pointer on Save.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
      ],
      shot: { mode: 'element', selector: '#session-state-controls', padding: 10 },
      cursor: { on: '#save-state-btn', at: 'center', state: 'hover' },
    },

    {
      id: 'sessions-renderer-controls-default',
      topic: 'sessions_sharing',
      name: 'save-restore-05-renderer-controls-default',
      describes:
        'The Cells and Performance controls at their defaults: Shader quality '
        + 'Full (lighting + fog), Point size, the Depth cues sliders, and under '
        + 'Performance Level of detail (LOD) and Frustum culling unchecked and '
        + 'Antialiasing checked.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
      ],
      shot: { mode: 'region', selectors: ['#renderer-controls', '#render-performance-controls'], padding: 8 },
    },

    {
      id: 'sessions-renderer-controls-changed',
      topic: 'sessions_sharing',
      name: 'save-restore-06-renderer-controls-changed',
      describes:
        'The same controls changed: Shader quality Light, Level of detail '
        + '(LOD) checked with Force LOD level at 4, and Frustum culling '
        + 'checked.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        ...RENDERER_CHANGES,
      ],
      shot: { mode: 'region', selectors: ['#renderer-controls', '#render-performance-controls'], padding: 8 },
    },

    {
      // The point of this image is that the four values below came back from a
      // file, not from the defaults the app starts at. The scenario therefore
      // changes them, saves, reloads the page, and restores — all in one
      // context, because the saved bundle only exists inside it.
      id: 'sessions-renderer-controls-restored',
      topic: 'sessions_sharing',
      name: 'save-restore-08-renderer-controls-restored',
      describes:
        'The same controls after a saved session was loaded into a fresh '
        + 'page: Shader quality Light, Level of detail (LOD) checked with Force '
        + 'LOD level 4, and Frustum culling checked, as they were saved.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        ...RENDERER_CHANGES,
        saveThenRestore(),
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        { run: async page => {
          const state = await page.evaluate(() => ({
            quality: document.querySelector('#hp-shader-quality').value,
            lod: document.querySelector('#hp-lod-enabled').checked,
            force: document.querySelector('#hp-lod-force').value,
            frustum: document.querySelector('#hp-frustum-culling').checked,
          }));
          // Assert the claim the caption makes. A restored control that came
          // back holding a default would make this image a lie, and a silent
          // one, because the panel looks the same either way.
          const expected = {
            quality: 'light', lod: true, force: '4', frustum: true,
          };
          for (const [key, value] of Object.entries(expected)) {
            if (state[key] !== value) {
              throw new Error(
                `restored renderer control ${key} is ${JSON.stringify(state[key])}, `
                + `expected ${JSON.stringify(value)}`
              );
            }
          }
        } },
      ],
      shot: { mode: 'region', selectors: ['#renderer-controls', '#render-performance-controls'], padding: 8 },
    },

    {
      id: 'sessions-save-confirmation',
      topic: 'sessions_sharing',
      name: 'save-restore-03-save-confirmation',
      describes:
        'A notification reading "Session saved successfully" after Save '
        + 'produced a session bundle.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
        NOTIFICATIONS_DRAINED,
        { run: async page => {
          await holdNotificationTimers(page);
          const download = page.waitForEvent('download');
          await page.click('#save-state-btn');
          await download;
          await page.waitForFunction(() => (
            document.querySelector('#notification-center')?.textContent
              ?.includes('Session saved successfully') === true
          ));
        } },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
    },

    {
      id: 'sessions-restore-confirmation',
      topic: 'sessions_sharing',
      name: 'save-restore-04-restore-confirmation',
      describes:
        'The one notification a completed restore publishes: the restore\'s '
        + 'progress entry, now titled "Loaded" and reading "Session fully '
        + 'restored.".',
      dataset: 'pancreas',
      steps: [
        READY,
        // No RENDERER_CHANGES here, unlike the scenarios either side. A
        // restored state with LOD on raises "LOD node mappings ready (N ms)"
        // as well, and that notification races the restore's own: across
        // three otherwise identical runs it landed above or below the
        // restore's entry, so the stack — and the bytes — changed from run to
        // run, and its N is a timing of the capture host besides.
        // `saveThenRestore` already waits for "Session fully restored.", the
        // only notification a completed restore publishes.
        { openSection: '#session-section' },
        saveThenRestore({ holdNotifications: true }),
        // Under 1000 ms on purpose: `holdNotificationTimers` defers every
        // in-page timer of a second or more, this step's own included.
        { waitForStable: '#notification-center', forMs: 700 },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
    },

    {
      id: 'sessions-restored-view',
      topic: 'sessions_sharing',
      name: 'save-restore-07-restored-view',
      describes:
        'The whole application window immediately after a session restore, '
        + 'with the dataset drawn in the viewer and the sidebar showing the '
        + 'restored control values.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        ...RENDERER_CHANGES,
        { openSection: '#session-section' },
        saveThenRestore(),
        { openSection: '#coloring-filtering-section' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
    },

    {
      id: 'sessions-official-sample-state-applied',
      topic: 'sessions_sharing',
      name: 'official-01-sample-state-applied',
      describes:
        'The whole application window shortly after the Pancreas built-in '
        + 'sample was chosen, showing the embedding coloured by cell type in a '
        + '3-D orbit view on the light Studio backdrop, with the Session '
        + 'accordion open beside it.',
      dataset: 'pancreas',
      steps: [
        READY,
        { run: async page => {
          // Every claim the caption makes, checked against the running app
          // rather than against the prose that describes the published state.
          const applied = await page.evaluate(() => ({
            shows: document.querySelector('#view-representation-select')
              ?.selectedOptions[0]?.textContent ?? null,
            navigation: document.querySelector('#navigation-mode')?.value ?? null,
            theme: document.querySelector('#theme-select')?.value ?? null,
            background: document.querySelector('#background-select')?.value ?? null,
          }));
          const expected = {
            shows: 'UMAP · 3D', navigation: 'orbit', theme: 'light', background: 'studio',
          };
          for (const [key, value] of Object.entries(expected)) {
            if (applied[key] !== value) {
              throw new Error(
                `published default state ${key} is ${JSON.stringify(applied[key])}, `
                + `expected ${JSON.stringify(value)}`
              );
            }
          }
        } },
        { openSection: '#session-section' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
    },

    {
      id: 'sessions-refusal-not-a-session-file',
      topic: 'sessions_sharing',
      name: 'refusal-02-not-a-session-file',
      describes:
        'A Cellucid notification reading "This file is not a Cellucid session '
        + 'file of format 2 (it does not start with CELLUCID_SESSION/2).", with a'
        + ' dismiss button at its right, after a file that is not a session '
        + 'bundle was chosen in Load.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
        // Drained first: triggered while the sample's own notifications were
        // still expiring, the capture failed with "#notification-center moved
        // while the image settled (380x346 became 380x292)", and the run that
        // passed showed two unrelated startup messages under the refusal.
        NOTIFICATIONS_DRAINED,
        { run: async page => {
          // A real file on disk, because the picker is the product's own and a
          // reader reaches this state by choosing the wrong file, not by
          // synthesising one. A 1x1 PNG named like a figure is exactly the
          // mistake the caption describes.
          const { mkdtemp, writeFile } = await import('node:fs/promises');
          const { tmpdir } = await import('node:os');
          const { join } = await import('node:path');
          const directory = await mkdtemp(join(tmpdir(), 'cellucid-docs-'));
          const decoy = join(directory, 'pancreas-figure.png');
          await writeFile(decoy, Buffer.from(
            'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQ'
            + 'DwAEhQGAhKmMIQAAAABJRU5ErkJggg==',
            'base64'
          ));
          const chooser = page.waitForEvent('filechooser');
          await page.click('#load-state-btn');
          await (await chooser).setFiles(decoy);
          await page.waitForFunction(() => (
            document.querySelector('#notification-center')?.textContent
              ?.includes('not a Cellucid session file') === true
          ));
        } },
        // The refusal is the restore's own progress row turning into its
        // failure (title, icon and message all change in that transition), so
        // the crop is measured only once the stack has stopped changing. An
        // ERROR notification carries no dismissal timer, so nothing expires
        // after this.
        { waitForStable: '#notification-center', forMs: 1000 },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
    },

    {
      id: 'sessions-refusal-different-dataset',
      topic: 'sessions_sharing',
      name: 'refusal-01-different-dataset',
      describes:
        'A notification titled "Loading failed" reading "Session file: This '
        + 'session was saved on a different dataset than the one that is open '
        + 'now (3,696 cells and 3,753 genes when it was saved, 219,731 cells and '
        + '5,754 genes now). Open the dataset the session was saved on, then '
        + 'load the session again."',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
        // A session records where its dataset came from, and Load opens a
        // sample's, a link's or a server's dataset itself before applying
        // it. A file from this computer is the one source a page cannot open
        // again, so the refusal is reached as a reader reaches it: the session
        // is saved on the Pancreas sample opened as an .h5ad from disk, and
        // loaded while another dataset (the Garcia sample) is open.
        { run: async page => {
          const before = await openFromDisk(page, 'file', await formatFixture('pancreas.h5ad'));
          await datasetReplaced(before).run(page);
          await NOTIFICATIONS_DRAINED.run(page);
          const download = page.waitForEvent('download');
          await page.click('#save-state-btn');
          const saved = await (await download).path();
          await page.selectOption('#dataset-select', 'dataset:local-demo:garcia');
          await page.waitForFunction(() => (
            document.querySelector('#filter-count')?.textContent
              ?.includes('Showing all 219,731 cells') === true
          ), undefined, { timeout: 180_000 });
          // The switch raises its own download and field notifications, and
          // any work still in flight keeps a progress row in the center, so
          // an empty center is the signal that the new dataset has settled.
          // A fixed pause here used to let them drain under the crop.
          await NOTIFICATIONS_DRAINED.run(page);
          const chooser = page.waitForEvent('filechooser');
          await page.click('#load-state-btn');
          await (await chooser).setFiles(saved);
          await page.waitForFunction(() => (
            document.querySelector('#notification-center .notification[data-type="error"]') !== null
          ), undefined, { timeout: 60_000 });
        } },
        { waitForStable: '#notification-center', forMs: 1000 },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
    },

    // -------------------------------- role cannot be determined (refusal)
    {
      id: 'annotation-state-repository-unreadable',
      topic: 'community_annotation',
      name: 'state-repository-unreadable',
      describes:
        'The whole application window with the Community Annotation panel back '
        + 'in its disconnected state and a notification explaining that the '
        + "user's role for the repository could not be determined and that the "
        + 'annotation repo was disconnected.',
      dataset: 'pancreas',
      // The refusal is raised during startup — the mock pre-binds the
      // repository, so the role lookup runs before any step can — and it is
      // an ERROR published with an 8-second timer. The image published before
      // this scenario held its timers showed the panel and no notification:
      // the refusal had expired while the frames settled. So the clock is
      // held from the first line of page code, and the startup notifications
      // that therefore also stay are dismissed with their own × buttons,
      // leaving the refusal as the only one, as a reader would.
      setup: async context => {
        await mock({ repoStatus: 404 })(context);
        await context.page.addInitScript(holdLongTimersInPage);
      },
      steps: [
        ...openPanel,
        { waitForText: '#notification-center', text: 'Unable to determine your role' },
        {
          run: async page => {
            const others = page.locator(
              '#notification-center > .notification',
              { hasNotText: 'Unable to determine your role' }
            );
            while ((await others.count()) > 0) {
              await others.first().locator('.notification-dismiss').click();
              await page.waitForTimeout(400);
            }
            await page.waitForFunction(() => {
              const center = document.querySelector('#notification-center');
              return center?.childElementCount === 1
                && center.textContent.includes('Unable to determine your role');
            });
          },
        },
      ],
      shot: { mode: 'window' },
      rings: [{ on: '#community-annotation-controls', label: '1' }],
    },
  ];
}

// ---------------------------------------------------------------------------
// The analysis panel: h_analysis.
//
// Every figure here needs the same three things before the analysis section is
// worth photographing — the pancreas sample, `cell_type` as the active
// categorical field, and one confirmed highlight page — so they are built
// together and share one prelude.
//
// The page is made the way the documentation tells a reader to make it:
// Annotation based highlighting, Alt+click, Confirm. There is no selector for
// a point inside a WebGL canvas, so the seed cell is found by scanning a fixed
// grid of Alt+clicks and reading the panel's own step readout, which is the
// same signal a user reads. A hard-coded pixel would silently photograph the
// wrong cell type the first time the default camera moves; the scan asserts
// the group it found instead, and fails loudly when no click reaches it.
// ---------------------------------------------------------------------------

/** Sidebar sections that are open on first load and are not analysis. */
const NON_ANALYSIS_SECTIONS = Object.freeze([
  '#session-section',
  '#visualization-section',
  '#compare-views-section',
  '#coloring-filtering-section',
  '#highlighted-cells-section',
]);

/**
 * Make `Page 1` out of the 262 `Ngn3 low EP` cells, exactly as a reader would.
 *
 * The scan skips the leftmost 340 CSS px: the sidebar overlays the canvas
 * there, and `#sidebar-toggle` sits at its top-left corner, so a click in that
 * band hides the whole sidebar instead of selecting anything.
 *
 * Alt without a modifier *replaces* the candidate set for a categorical field,
 * but the standing preview changes which cells are pickable, so each miss is
 * cancelled before the next attempt. That keeps every attempt independent and
 * makes the successful one `Step 1`, which is what the page label records.
 *
 * @param {import('@playwright/test').Page} page
 */
async function confirmNgn3LowEndocrineProgenitorPage(page) {
  await stageNgn3LowEndocrineProgenitorStep(page);
  await page.evaluate(() =>
    document.querySelector('#annotation-confirm-btn').click()
  );
  await page.waitForFunction(() => (
    document.querySelector('.highlight-page-tab-count')?.textContent === '262'
  ));
}

/**
 * Stage the 262 `Ngn3 low EP` cells as `Step 1` of Page 1, unconfirmed: the
 * Alt+click scan of `confirmNgn3LowEndocrineProgenitorPage`, stopped before
 * Confirm.
 *
 * @param {import('@playwright/test').Page} page
 */
async function stageNgn3LowEndocrineProgenitorStep(page) {
  const readout = page.locator('#highlight-mode-description');
  const left = 340;
  const columns = 14;
  const rows = 10;
  for (let row = 0; row < rows; row += 1) {
    for (let column = 0; column < columns; column += 1) {
      await page.mouse.move(
        left + ((column + 0.5) * (1440 - left)) / columns,
        ((row + 0.5) * 1000) / rows
      );
      await page.keyboard.down('Alt');
      await page.mouse.down();
      await page.mouse.up();
      await page.keyboard.up('Alt');
      if (/Step 1: 262 cells/.test((await readout.textContent()) ?? '')) {
        await page.mouse.move(4, 995);
        return;
      }
      await page.evaluate(() =>
        document.querySelector('#annotation-cancel-btn')?.click()
      );
    }
  }
  throw new Error(
    'No Alt+click in the scan reached the 262-cell Ngn3 low EP group; the '
      + 'default camera or the sample changed.'
  );
}

/**
 * Leave the analysis section as the only open one, so a panel crop starts at
 * the top of the sidebar's scroller rather than halfway down it.
 *
 * @param {import('@playwright/test').Page} page
 */
async function showOnlyAnalysisSection(page) {
  await page.evaluate(sections => {
    for (const selector of sections) {
      document.querySelector(selector).open = false;
    }
    document.querySelector('#page-analysis-section').open = true;
  }, NON_ANALYSIS_SECTIONS);
}

/**
 * Open one analysis mode and bring its header to the top of the sidebar.
 *
 * The accordion header is a `<div role="button">` outside the viewport until
 * the sidebar is scrolled, so it is activated through its own click handler
 * rather than through the mouse.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} mode
 */
async function openAnalysisMode(page, mode) {
  await page.evaluate(
    id => document.querySelector(`#analysis-header-${id}`).click(),
    mode
  );
  await page.waitForSelector(`#analysis-panel-${mode}`, { state: 'visible' });
  await scrollAnalysisModeToTop(page, mode);
}

/**
 * @param {import('@playwright/test').Page} page
 * @param {string} mode
 */
async function scrollAnalysisModeToTop(page, mode) {
  await page.evaluate(
    id => document.querySelector(`#analysis-header-${id}`).scrollIntoView({ block: 'start' }),
    mode
  );
}

/**
 * Press a mode's Run button and wait for the result the panel publishes.
 *
 * Differential expression and marker discovery are the two modes that do not
 * auto-run, and both are minutes-scale on a large sample; the wait is given its
 * own budget rather than the context default, so a slow machine reports a real
 * failure instead of a timeout at 45 seconds.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} mode
 */
async function runAnalysisAndWaitForResult(page, mode) {
  await page.evaluate(
    id => document.querySelector(`#analysis-panel-${id} .analysis-run-btn`).click(),
    mode
  );
  await page.waitForFunction(
    id => (
      document.querySelector(`#analysis-panel-${id} .analysis-expand-btn`) !== null
        && document.querySelector(`#analysis-panel-${id} .analysis-preview-plot .main-svg`) !== null
    ),
    mode,
    { timeout: 300_000 }
  );
}

/**
 * Open the expanded view from a mode's own Expand button and wait for its plot.
 *
 * The modal grows into place, and `capture.mjs` measures the crop rectangle
 * once, before the settle loop. A crop read while the overlay is still opening
 * is smaller than the overlay that is eventually drawn, so the image comes out
 * with its own edges cut off — waiting for the plot to exist is not enough, the
 * geometry has to stop moving too.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} mode
 */
async function expandAnalysisResult(page, mode) {
  // A real pointer click, as a reader makes it: a scripted `.click()` counts
  // as keyboard activation, and the modal would then draw its keyboard focus
  // ring on the close button, which a mouse user never sees.
  await page.click(`#analysis-panel-${mode} .analysis-expand-btn`);
  await page.mouse.move(4, 995);
  await page.waitForSelector('.analysis-modal-content', { state: 'visible' });
  await page.waitForFunction(
    () => document.querySelector('.analysis-modal-plot .main-svg') !== null,
    undefined,
    { timeout: 300_000 }
  );
  await page.waitForFunction(async () => {
    const read = () => {
      const box = document.querySelector('.analysis-modal-content')?.getBoundingClientRect();
      return box === undefined ? null : `${box.x},${box.y},${box.width},${box.height}`;
    };
    const first = read();
    if (first === null) return false;
    await new Promise(resolve => setTimeout(resolve, 400));
    return read() === first;
  }, undefined, { polling: 400, timeout: 60_000 });
}

/** @returns {Array<Object>} */
function buildAnalysisScenarios() {
  // Pancreas, `cell_type` active, one confirmed page and its free complement.
  const ANALYSIS_READY = [
    PANCREAS_READY,
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label: 'cell_type' },
    { waitForText: '#stats', text: 'Field: cell_type (category)' },
    { run: confirmNgn3LowEndocrineProgenitorPage },
    { run: showOnlyAnalysisSection },
  ];

  // `cell_type` across both pages as a bar plot, auto-run, left with its
  // sidebar result on screen.
  const DETAILED_CATEGORICAL_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'detailed') },
    { select: '#detailed-type', value: 'categorical' },
    { waitFor: '#detailed-categorical' },
    { select: '#detailed-categorical', value: 'cell_type' },
    { waitFor: '#analysis-plot-type' },
    { run: page => page.evaluate(() =>
      document.querySelector('#analysis-panel-detailed .analysis-page-tabs .chip[data-page-id="restof__page_1"]').click()
    ) },
    // Both pages among the drawn plot's traces is the only honest signal that
    // the second page reached the analysis rather than only the tab strip.
    // (The plot's legend is not it: a legend that would squeeze the small
    // sidebar plot is left out.)
    { run: page => page.waitForFunction(() => (
      (document.querySelector('#analysis-panel-detailed .js-plotly-plot')?.data ?? [])
        .map(trace => trace.name)
        .join('|') === 'Page 1|Rest of Page 1'
    )) },
    { run: page => scrollAnalysisModeToTop(page, 'detailed') },
  ];

  // Differential expression of Page 1 against its complement, run from the
  // form and left with its sidebar result on screen.
  const DIFFERENTIAL_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'differential') },
    // Compare is Pages by default and the test Wilcoxon rank-sum, so the two
    // groups are the whole choice.
    { waitForStable: '#analysis-panel-differential select[name="pageA"]' },
    { select: '#analysis-panel-differential select[name="pageA"]', value: 'page_1' },
    { waitForStable: '#analysis-panel-differential select[name="pageB"]' },
    { select: '#analysis-panel-differential select[name="pageB"]', value: 'restof__page_1' },
    { run: page => runAnalysisAndWaitForResult(page, 'differential') },
    { run: page => scrollAnalysisModeToTop(page, 'differential') },
  ];

  /**
   * Marker genes on `cell_type` (Group by), shown as the result's **Show**
   * choice: the heatmap a run opens on, or the clustered heatmap with its
   * dendrograms. The test is the form's default, Wilcoxon rank-sum.
   *
   * @param {'heatmap'|'clustered'} show
   * @returns {Array<Object>}
   */
  const markerGenes = show => [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'genesPanel') },
    { waitForStable: '#analysis-panel-genesPanel select[name="groupBy"]' },
    { select: '#analysis-panel-genesPanel select[name="groupBy"]', value: 'cell_type' },
    { run: page => runAnalysisAndWaitForResult(page, 'genesPanel') },
    ...(show === 'heatmap' ? [] : [
      { click: `#analysis-panel-genesPanel .segmented[aria-label="Show"] .segmented-btn[data-value="${show}"]` },
      // The clustered heatmap is drawn once the workers return the trees: the
      // dendrograms' own subplots (x2y2 for the groups, x3y3 for the genes)
      // beside the heatmap are the signal it is on screen.
      { run: page => page.waitForFunction(() => document.querySelector(
        '#analysis-panel-genesPanel .analysis-preview-plot .cartesianlayer .x2y2, '
        + '#analysis-panel-genesPanel .analysis-preview-plot .cartesianlayer .x3y3'
      ) !== null, undefined, { timeout: 120_000 }) },
    ]),
    { run: page => scrollAnalysisModeToTop(page, 'genesPanel') },
  ];

  // Differential expression, Page 1 against its complement, not run yet.
  const DIFFERENTIAL_FORM = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'differential') },
    { waitForStable: '#analysis-panel-differential select[name="pageA"]' },
    { select: '#analysis-panel-differential select[name="pageA"]', value: 'page_1' },
    { waitForStable: '#analysis-panel-differential select[name="pageB"]' },
    { select: '#analysis-panel-differential select[name="pageB"]', value: 'restof__page_1' },
    { run: page => scrollAnalysisModeToTop(page, 'differential') },
  ];

  /** Choose one variable of a two-step variable control (type, then field). */
  const chooseVariable = (mode, option, type, field) => ({
    run: async page => {
      const root = `#analysis-panel-${mode} [data-option="${option}"]`;
      await page.selectOption(`${root} select >> nth=0`, type);
      await page.waitForFunction(
        selector => document.querySelectorAll(`${selector} select`).length > 1,
        root
      );
      await page.selectOption(`${root} select >> nth=1`, { label: field });
    },
  });

  // Correlation of the two cell-cycle scores across both pages.
  const CORRELATION_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'correlation') },
    chooseVariable('correlation', 'variableX', 'continuous', 'S_score'),
    chooseVariable('correlation', 'variableY', 'continuous', 'G2M_score'),
    // Both pages, coloured by cell type: Colour by sits under More options.
    { click: '#analysis-panel-correlation .chip[data-page-id="restof__page_1"]' },
    { run: async page => {
      const more = page.locator('#analysis-panel-correlation .disclosure-header', { hasText: 'More options' });
      if ((await more.getAttribute('aria-expanded')) !== 'true') await more.click();
      const display = page.locator('#analysis-panel-correlation .disclosure-header', { hasText: 'Display' });
      if ((await display.count()) > 0 && (await display.first().getAttribute('aria-expanded')) !== 'true') {
        await display.first().click();
      }
      await page.selectOption('#analysis-panel-correlation select[name="colorBy"]', { label: 'cell_type' });
      await page.mouse.move(4, 995);
    } },
    { run: page => runAnalysisAndWaitForResult(page, 'correlation') },
    { run: page => scrollAnalysisModeToTop(page, 'correlation') },
  ];

  /** Type genes into a mode's first gene set, as a paste does. */
  const addGenes = (mode, genes) => ({
    run: async page => {
      const input = page.locator(`#analysis-panel-${mode} .gene-set-input input`).first();
      await input.fill(genes);
      await input.press('Enter');
      await page.waitForFunction(
        ([selector, count]) => document.querySelectorAll(`${selector} .gene-chip`).length >= count,
        [`#analysis-panel-${mode}`, genes.split(/[\s,]+/).filter(Boolean).length]
      );
      // Focus leaves the box as the next click (Run) takes it.
      await input.blur();
    },
  });

  // Gene signature of three beta-cell genes, both pages.
  const SIGNATURE_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'signature') },
    // A mode opens on Page 1 alone; the complement is added for the comparison.
    { click: '#analysis-panel-signature .chip[data-page-id="restof__page_1"]' },
    addGenes('signature', 'Ins1, Ins2, Iapp'),
    { run: page => runAnalysisAndWaitForResult(page, 'signature') },
    { run: page => scrollAnalysisModeToTop(page, 'signature') },
  ];

  // Dot plot of eight pancreatic markers by cluster.
  const DOTPLOT_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'dotplot') },
    { waitForStable: '#analysis-panel-dotplot select[name="groupBy"]' },
    { select: '#analysis-panel-dotplot select[name="groupBy"]', value: 'clusters' },
    addGenes('dotplot', 'Krt19, Sox9, Neurog3, Fev, Ins1, Gcg, Sst, Ghrl'),
    { run: page => runAnalysisAndWaitForResult(page, 'dotplot') },
    { run: page => scrollAnalysisModeToTop(page, 'dotplot') },
  ];

  // Composition of the Kanemaru heart cells by sample, tested between sexes.
  const COMPOSITION_RESULT = [
    { run: page => page.waitForFunction(
      () => document.documentElement.dataset.startupPhase === 'ready',
      undefined,
      { timeout: 180_000 }
    ) },
    { run: showOnlyAnalysisSection },
    { run: page => openAnalysisMode(page, 'composition') },
    { waitForStable: '#analysis-panel-composition select[name="cellTypes"]' },
    { select: '#analysis-panel-composition select[name="cellTypes"]', value: 'LVL1' },
    { select: '#analysis-panel-composition select[name="samples"]', value: 'sample_ID' },
    { select: '#analysis-panel-composition select[name="condition"]', value: 'sex' },
    { run: page => runAnalysisAndWaitForResult(page, 'composition') },
    { run: page => scrollAnalysisModeToTop(page, 'composition') },
  ];

  return [
    // ----------------------------------------- Detailed, categorical variable
    {
      id: 'analysis-detailed-categorical-preview',
      topic: 'analysis',
      name: 'detailed-categorical',
      describes:
        'The Detailed panel with Variable set to Categorical obs then '
        + 'cell_type, a Compare pages block with Page 1 and Rest of Page 1 '
        + 'both selected, a Plot Type select reading Bar Plot under a pointer, '
        + 'and a grouped bar plot of cell-type counts beneath it whose '
        + 'horizontal axis carries two category names, Ductal and Delta, at '
        + 'its two ends.',
      dataset: 'pancreas',
      steps: DETAILED_CATEGORICAL_RESULT,
      shot: { mode: 'region', selectors: ['#analysis-header-detailed', '#analysis-panel-detailed'], padding: 6 },
      cursor: { on: '#analysis-plot-type', at: 'center', state: 'hover' },
    },

    // ------------------------------------ Differential expression, in the bar
    {
      id: 'analysis-de-sidebar-result',
      topic: 'analysis',
      name: 'de-sidebar-result',
      describes:
        'The sidebar result of a differential expression run: a small volcano '
        + 'plot with the genes up in Page 1 in red on the right, those down in '
        + 'blue on the left and the rest grey, one of its gene labels shown, and '
        + 'the Expand button below it under the pointer.',
      dataset: 'pancreas',
      steps: DIFFERENTIAL_RESULT,
      shot: {
        mode: 'region',
        selectors: [
          '#analysis-panel-differential .analysis-preview-container',
          '#analysis-panel-differential .analysis-expand-btn',
        ],
        padding: 4,
      },
      cursor: {
        on: '#analysis-panel-differential .analysis-expand-btn',
        at: 'center',
        state: 'hover',
      },
    },

    // ------------------------------------ Differential expression, full view
    {
      id: 'analysis-de-volcano-expanded',
      topic: 'analysis',
      name: 'de-volcano-expanded',
      describes:
        'A large modal headed Differential expression: Page 1 vs Rest of Page '
        + '1: a volcano plot with seven labelled genes top-left over a legend '
        + 'reading Up (395), Down (698) and Not significant (2660), an Export row'
        + ' of PNG, SVG and CSV buttons above the Plot options on the right, a '
        + 'Summary statistics table bottom-left, and a Statistical analysis table'
        + ' of the top genes bottom-right.',
      dataset: 'pancreas',
      steps: [
        ...DIFFERENTIAL_RESULT,
        { run: page => expandAnalysisResult(page, 'differential') },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
    },

    // ----------------------------------------------- Marker genes, Ranked
    {
      id: 'analysis-marker-genes-ranked',
      topic: 'analysis',
      name: 'marker-genes-ranked',
      describes:
        'The Marker genes panel grouped by cell_type, with its Run button, '
        + 'and beneath it the result header, the Show choice (List, Heatmap, '
        + 'Clustered) on Heatmap, and a heatmap of eight cell-type rows '
        + 'against their top marker genes with a Z-score colour bar.',
      dataset: 'pancreas',
      steps: markerGenes('heatmap'),
      shot: { mode: 'region', selectors: ['#analysis-header-genesPanel', '#analysis-panel-genesPanel'], padding: 6 },
    },

    // --------------------------------------------- Marker genes, Clustered
    {
      id: 'analysis-marker-genes-clustered',
      topic: 'analysis',
      name: 'marker-genes-clustered',
      describes:
        'A heatmap with cell-type groups down one axis and gene names across '
        + 'the other, dendrograms attached to both axes, and a colour bar for '
        + 'the Z-scored expression.',
      dataset: 'pancreas',
      steps: [
        ...markerGenes('clustered'),
        { run: page => expandAnalysisResult(page, 'genesPanel') },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-plot' },
    },

    // ------------------------------------ Marker genes, the whole full view
    {
      id: 'analysis-marker-genes-expanded',
      topic: 'analysis',
      name: 'marker-genes-expanded',
      describes:
        'The expanded Marker genes view of cell_type with a clustered '
        + 'heatmap and dendrograms in the plot area, an Export row of PNG, '
        + 'SVG and CSV buttons above the plot options, the summary statistics '
        + 'naming the field, method and counts, and the table of top marker '
        + 'genes, with a pointer on the CSV button; the run-time Duration row '
        + 'is omitted.',
      dataset: 'pancreas',
      steps: [
        ...markerGenes('clustered'),
        { run: page => expandAnalysisResult(page, 'genesPanel') },
        // The summary table ends with a Duration row: the wall-clock time the
        // discovery took on the capture host ("82.1s", then "1m 32s" a day
        // later). That is a measurement of the machine, not of Cellucid —
        // the reason the benchmark scenario hides its live readout — and it
        // made every recapture differ. The row is removed from the frame and
        // the alt text says so.
        {
          run: page => page.evaluate(() => {
            const rows = [...document.querySelectorAll(
              '.analysis-modal-content .analysis-stats-table tr'
            )].filter(row => row.cells[0]?.textContent === 'Duration');
            if (rows.length !== 1) {
              throw new Error(`expected one Duration row, found ${rows.length}`);
            }
            rows[0].style.display = 'none';
          }),
        },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
      cursor: {
        on: '.analysis-export-toolbar button[title="Download data as CSV"]',
        at: 'center',
        state: 'hover',
      },
    },

    // ------------------------------------------- Detailed, the full view
    {
      id: 'analysis-detailed-expanded',
      topic: 'analysis',
      name: 'detailed-expanded',
      describes:
        'The expanded Detailed modal headed COMPARING: CELL_TYPE with a '
        + 'grouped bar plot top-left whose horizontal axis prints five of the '
        + 'eight category names, an Export row and PLOT OPTIONS column on the '
        + 'right, a SUMMARY STATISTICS table bottom-left, and a STATISTICAL '
        + 'ANALYSIS panel bottom-right holding a chi-squared card reading N/A, '
        + 'with a pointer on the CSV button.',
      dataset: 'pancreas',
      steps: [
        ...DETAILED_CATEGORICAL_RESULT,
        { run: page => expandAnalysisResult(page, 'detailed') },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
      cursor: {
        on: '.analysis-export-toolbar button[title="Download data as CSV"]',
        at: 'center',
        state: 'hover',
      },
    },

    // --------------------------- Detailed preview, as the auto-run example
    {
      id: 'analysis-autorun-detailed-preview',
      topic: 'analysis',
      name: 'autorun-vs-run-a',
      describes:
        'The Detailed panel holding a finished result — Variable set to '
        + 'Categorical obs then cell_type, two pages selected, Bar Plot chosen '
        + 'and a grouped bar plot drawn beneath — with no Run button anywhere '
        + 'in the panel.',
      dataset: 'pancreas',
      steps: DETAILED_CATEGORICAL_RESULT,
      shot: { mode: 'region', selectors: ['#analysis-header-detailed', '#analysis-panel-detailed'], padding: 6 },
    },

    // ------------------------------------ Volcano thresholds, the A/B pair
    {
      id: 'analysis-volcano-threshold-permissive',
      topic: 'analysis',
      name: 'volcano-threshold-a',
      describes:
        'A volcano plot at a log2 fold-change threshold of 0.5, with broad red '
        + 'and blue clouds of coloured points either side of a narrow grey '
        + 'band, and a legend below the axis title counting each colour.',
      dataset: 'pancreas',
      steps: [
        ...DIFFERENTIAL_RESULT,
        { run: page => expandAnalysisResult(page, 'differential') },
        { run: page => setVolcanoFoldChangeThreshold(page, 0.5) },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-plot' },
    },
    {
      id: 'analysis-volcano-threshold-strict',
      topic: 'analysis',
      name: 'volcano-threshold-b',
      describes:
        'The same volcano plot at a log2 fold-change threshold of 3, with the '
        + 'threshold lines pushed far apart and almost every point now grey, '
        + 'and a legend below the axis title counting each colour.',
      dataset: 'pancreas',
      steps: [
        ...DIFFERENTIAL_RESULT,
        { run: page => expandAnalysisResult(page, 'differential') },
        { run: page => setVolcanoFoldChangeThreshold(page, 3) },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-plot' },
    },

    // ------------------------------------------- the panel, every mode shut
    {
      id: 'analysis-panel-modes',
      topic: 'analysis',
      name: 'analysis-panel-modes',
      describes:
        'The Analysis panel with every analysis collapsed, in three groups: '
        + 'Genes (Marker genes, Differential expression, Dot plot, Gene '
        + 'signature), Cells (Composition) and Explore (Quick, Detailed, '
        + 'Correlation); each row carries a chevron, an (i) and a copy-window '
        + 'button.',
      dataset: 'pancreas',
      steps: [...ANALYSIS_READY, { scrollIntoView: '#page-analysis-section' }],
      shot: { mode: 'element', selector: '#page-analysis-section', padding: 6 },
    },

    // ------------------------------------ Differential expression, not run
    {
      id: 'analysis-de-form-before-run',
      topic: 'analysis',
      name: 'autorun-vs-run-b',
      describes:
        'The Differential expression form before a run: Compare set to Pages, '
        + 'Group A Page 1 and Group B Rest of Page 1, the Swap A and B link, More'
        + ' options collapsed, and the Run button under the pointer; no result '
        + 'yet.',
      dataset: 'pancreas',
      steps: DIFFERENTIAL_FORM,
      shot: { mode: 'region', selectors: ['#analysis-header-differential', '#analysis-panel-differential'], padding: 6 },
      cursor: { on: '#analysis-panel-differential .analysis-run-btn', at: 'center', state: 'hover' },
    },

    // -------------------------------------- the page chips, Rest of Page 1
    {
      id: 'analysis-rest-of-page-chips',
      topic: 'analysis',
      name: 'rest-of-page-selector',
      describes:
        'The page chips of the Detailed panel: Page 1 with its colour swatch '
        + 'and the count 262, and Rest of Page 1 with its own swatch and count, '
        + 'both selected, with the pointer on Rest of Page 1.',
      dataset: 'pancreas',
      steps: DETAILED_CATEGORICAL_RESULT,
      shot: { mode: 'element', selector: '#analysis-panel-detailed .analysis-page-tabs', padding: 6 },
      cursor: {
        on: '#analysis-panel-detailed .analysis-page-tabs .chip[data-page-id="restof__page_1"]',
        at: 'center',
        state: 'hover',
      },
    },

    // ----------------------------------------------- the Expand button
    {
      id: 'analysis-expand-button',
      topic: 'analysis',
      name: 'expand-button',
      describes:
        'A Detailed result in the sidebar: the small bar plot of cell_type for '
        + 'the two pages and, under it, the Expand button with the pointer on it.',
      dataset: 'pancreas',
      steps: DETAILED_CATEGORICAL_RESULT,
      shot: {
        mode: 'region',
        selectors: ['#analysis-panel-detailed .analysis-preview-plot', '#analysis-panel-detailed .analysis-expand-btn'],
        padding: 6,
      },
      cursor: { on: '#analysis-panel-detailed .analysis-expand-btn', at: 'center', state: 'hover' },
    },

    // --------------------------------------- a mode in a floating window
    {
      id: 'analysis-floating-window',
      topic: 'analysis',
      name: 'floating-analysis-window',
      setup: pinHighlightGlow,
      describes:
        'The whole Cellucid window with the Correlation analysis copied into a '
        + 'floating window over the point cloud: its own header with copy, '
        + 'maximise, close and collapse buttons, and its own copy of the '
        + 'Correlation form.',
      dataset: 'pancreas',
      steps: [
        ...ANALYSIS_READY,
        { run: page => page.evaluate(() => document.querySelector(
          '.analysis-mode[data-mode="correlation"] .analysis-mode-copy-btn'
        ).click()) },
        { waitFor: '.analysis-window-panel' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
    },

    // ------------------------------------------------------------- Quick
    {
      id: 'analysis-quick-insights',
      topic: 'analysis',
      name: 'quick-insights',
      describes:
        'The Quick panel for Page 1, 262 cells: a Composition section with a '
        + 'bar for clusters_coarse and one for clusters, each a single Ngn3 low '
        + 'EP segment named beneath at 100.0%, a Statistics section with the '
        + 'mean, median and standard deviation of S_score and G2M_score, and the '
        + 'collapsed Page selection.',
      dataset: 'pancreas',
      steps: [
        ...ANALYSIS_READY,
        { run: page => openAnalysisMode(page, 'simple') },
        { waitFor: '#analysis-panel-simple .insights-section' },
        { run: page => scrollAnalysisModeToTop(page, 'simple') },
      ],
      shot: { mode: 'region', selectors: ['#analysis-header-simple', '#analysis-panel-simple'], padding: 6 },
    },
    {
      id: 'analysis-quick-no-pages',
      topic: 'analysis',
      name: 'analysis-empty-no-pages',
      describes:
        'The Quick panel before any cells were confirmed into a page: it '
        + 'reports Page 1 with 0 cells and has no composition or statistics to '
        + 'show.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'cell_type' },
        { waitForText: '#stats', text: 'Field: cell_type (category)' },
        { run: showOnlyAnalysisSection },
        { run: page => openAnalysisMode(page, 'simple') },
        { wait: 1500 },
        { run: page => scrollAnalysisModeToTop(page, 'simple') },
      ],
      shot: { mode: 'region', selectors: ['#analysis-header-simple', '#analysis-panel-simple'], padding: 6 },
    },

    // -------------------------------------------------------- Correlation
    {
      id: 'analysis-correlation-result',
      topic: 'analysis',
      name: 'correlation-results',
      describes:
        'The Correlation form with Page 1 and Rest of Page 1 selected, X axis '
        + 'S_score and Y axis G2M_score (Continuous obs), More options open on '
        + 'Display with Colour by cell_type, and the Run button.',
      dataset: 'pancreas',
      steps: CORRELATION_RESULT,
      shot: {
        mode: 'region',
        selectors: ['#analysis-header-correlation', '#analysis-panel-correlation .analysis-run-btn'],
        padding: 6,
      },
    },
    {
      id: 'analysis-correlation-expanded',
      topic: 'analysis',
      name: 'correlation-expanded',
      describes:
        'The expanded Correlation view of S_score against G2M_score: the '
        + 'scatter with its per-page trend lines, the options beside it, and the '
        + 'summary and statistics tables with r, r², p and n per page.',
      dataset: 'pancreas',
      steps: [...CORRELATION_RESULT, { run: page => expandAnalysisResult(page, 'correlation') }],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
    },

    // --------------------------------------------------- Gene signature
    {
      id: 'analysis-signature-result',
      topic: 'analysis',
      name: 'gene-signature',
      describes:
        'The Gene signature panel with one gene set holding Ins1, Ins2 and Iapp, '
        + 'both pages selected, and the signature score of each page drawn '
        + 'beneath.',
      dataset: 'pancreas',
      steps: SIGNATURE_RESULT,
      shot: { mode: 'region', selectors: ['#analysis-header-signature', '#analysis-panel-signature'], padding: 6 },
    },
    {
      id: 'analysis-signature-expanded',
      topic: 'analysis',
      name: 'gene-signature-expanded',
      describes:
        'The expanded Gene signature view of Ins1, Ins2 and Iapp: the score '
        + 'distribution of each page, the options beside it, a summary table '
        + 'per page and the two-group tests with their effect sizes.',
      dataset: 'pancreas',
      steps: [...SIGNATURE_RESULT, { run: page => expandAnalysisResult(page, 'signature') }],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
    },

    // ---------------------------------------------------------- Dot plot
    {
      id: 'analysis-dotplot-result',
      topic: 'analysis',
      name: 'dotplot-result',
      describes:
        'The Dot plot form with Group by clusters and one gene set, Set 1, '
        + 'holding eight pancreatic markers as chips (Krt19, Sox9, Neurog3, Fev, '
        + 'Ins1, Gcg, Sst, Ghrl), the box to add or paste more, Clear, Add set '
        + 'and Import, More options, the Run button, and under it the line the '
        + 'result opens with: 8 genes, 8 groups of clusters and the cell count.',
      dataset: 'pancreas',
      steps: DOTPLOT_RESULT,
      shot: {
        mode: 'region',
        selectors: ['#analysis-header-dotplot', '#analysis-panel-dotplot .analysis-result-header'],
        padding: 6,
      },
    },
    {
      id: 'analysis-dotplot-expanded',
      topic: 'analysis',
      name: 'dotplot-expanded',
      describes:
        'The expanded Dot plot of the eight pancreatic markers across the '
        + 'eight clusters: one row per cluster, dot colour the mean expression '
        + 'and dot size the fraction of cells expressing it, with its colour '
        + 'bar and size legend, Export, and a summary table below.',
      dataset: 'pancreas',
      steps: [...DOTPLOT_RESULT, { run: page => expandAnalysisResult(page, 'dotplot') }],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
    },

    // ------------------------------------------------------- Composition
    // Kanemaru: 131,636 heart cells from 10 samples of two sexes, so a
    // composition by sample with a condition has replicates to test.
    {
      id: 'analysis-composition-result',
      topic: 'analysis',
      name: 'composition-result',
      describes:
        'The Composition panel on the Kanemaru heart sample with Cell types '
        + 'LVL1, Samples sample_ID and Condition sex, and the stacked bars '
        + 'beneath: one bar per sample, grouped under its condition, split by '
        + 'cell type.',
      dataset: 'kanemaru',
      steps: COMPOSITION_RESULT,
      shot: { mode: 'region', selectors: ['#analysis-header-composition', '#analysis-panel-composition'], padding: 6 },
    },
    {
      id: 'analysis-composition-expanded',
      topic: 'analysis',
      name: 'composition-expanded',
      describes:
        'The expanded Composition view on the Kanemaru heart sample: the '
        + 'stacked bars per sample by sex, the options beside them, and the '
        + 'table of each cell type\'s proportions and test result.',
      dataset: 'kanemaru',
      steps: [...COMPOSITION_RESULT, { run: page => expandAnalysisResult(page, 'composition') }],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
    },

    // ----------------------------- Differential expression, running, cancelled
    // The sample is read from a host of this run that stops answering when
    // Run is pressed, so the job holds at its first step for the photograph;
    // Cancel then ends it as a reader's would.
    {
      id: 'analysis-de-running',
      topic: 'analysis',
      name: 'de-running',
      describes:
        'The Differential expression panel while it runs: the Run button '
        + 'disabled, and under it a box naming the comparison (Page 1 against '
        + 'Rest of Page 1, with their cell counts), the step under way with how '
        + 'many genes of how many are done, a progress bar and a Cancel button.',
      dataset: 'pancreas',
      query: async () => ({ sampleCatalog: `${await gatedSampleCatalog()}/` }),
      steps: [
        ...DIFFERENTIAL_FORM,
        { run: () => { SAMPLE_CATALOG_GATE.hold = true; } },
        { run: page => page.evaluate(() => document.querySelector('#analysis-panel-differential .analysis-run-btn').click()) },
        { waitFor: '#analysis-panel-differential .job-status' },
        { waitForText: '#analysis-panel-differential .job-status', text: 'genes' },
        { wait: 1500 },
        { run: page => scrollAnalysisModeToTop(page, 'differential') },
      ],
      shot: { mode: 'region', selectors: ['#analysis-header-differential', '#analysis-panel-differential'], padding: 6 },
      settle: { attempts: 6, gapMs: 300 },
    },
    {
      id: 'analysis-de-cancelled',
      topic: 'analysis',
      name: 'de-cancelled',
      describes:
        'The Differential expression panel after Cancel: the Run button '
        + 'enabled again, and under it a box saying the run was cancelled, '
        + 'after how long, that nothing was kept, and a Run again button.',
      dataset: 'pancreas',
      query: async () => ({ sampleCatalog: `${await gatedSampleCatalog()}/` }),
      steps: [
        ...DIFFERENTIAL_FORM,
        { run: () => { SAMPLE_CATALOG_GATE.hold = true; } },
        { run: page => page.evaluate(() => document.querySelector('#analysis-panel-differential .analysis-run-btn').click()) },
        { waitForText: '#analysis-panel-differential .job-status', text: 'genes' },
        { wait: 1500 },
        { click: '#analysis-panel-differential .job-status-action' },
        { waitForText: '#analysis-panel-differential .job-status', text: 'Run again' },
        { run: page => scrollAnalysisModeToTop(page, 'differential') },
      ],
      shot: { mode: 'region', selectors: ['#analysis-header-differential', '#analysis-panel-differential'], padding: 6 },
    },

    // ------------------------------ the selection lifecycle, before analysis
    {
      id: 'analysis-life-01-pending-step',
      topic: 'analysis',
      name: 'life-01-confirm-a-selection',
      describes:
        'The Highlighting panel with the Annotation tool chosen and a pending '
        + 'step reading "Step 1: 262 cells" with Confirm, undo, redo and Cancel, '
        + 'above Page 1, which still holds no cells. The pointer is on Confirm.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'cell_type' },
        { waitForText: '#stats', text: 'Field: cell_type (category)' },
        { run: stageNgn3LowEndocrineProgenitorStep },
        onlySections('#highlighted-cells-section'),
      ],
      shot: { mode: 'element', selector: '#highlighted-cells-section', padding: 6 },
      cursor: { on: '#annotation-confirm-btn', at: 'right', dx: -12, state: 'hover' },
    },
    {
      id: 'analysis-life-02-page-ready',
      topic: 'analysis',
      name: 'life-02-two-pages-ready',
      describes:
        'The same Highlighting panel after Confirm: Page 1 holds 262 cells, and'
        + ' the list below it reads 262 cells highlighted with Clear, above the '
        + 'confirmed group with its checkbox and remove button.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'cell_type' },
        { waitForText: '#stats', text: 'Field: cell_type (category)' },
        { run: confirmNgn3LowEndocrineProgenitorPage },
        onlySections('#highlighted-cells-section'),
        { run: page => page.mouse.move(4, 995) },
      ],
      shot: { mode: 'element', selector: '#highlighted-cells-section', padding: 6 },
    },

    // ------------------------------- Differential expression, its two groups
    {
      id: 'analysis-de-page-selection',
      topic: 'analysis',
      name: 'de-page-selection',
      describes:
        'The top of the Differential expression form: Compare set to Pages, '
        + 'Group A reading Page 1 and Group B reading Rest of Page 1.',
      dataset: 'pancreas',
      steps: DIFFERENTIAL_FORM,
      shot: {
        mode: 'region',
        selectors: [
          '#analysis-panel-differential [data-option="compare"]',
          '#analysis-panel-differential [data-option="pageB"]',
        ],
        padding: 10,
      },
    },
    {
      id: 'analysis-de-fdr-denominator',
      topic: 'analysis',
      name: 'de-fdr-denominator',
      describes:
        'The statistics table of an expanded differential expression result: '
        + 'Genes tested, Not tested, Significant at the adjusted-p and |log2FC| '
        + 'thresholds, Up in A, Down in A, the p-value convention, Multiple '
        + 'testing (Benjamini–Hochberg), the expression matrix and the fold '
        + 'change.',
      dataset: 'pancreas',
      steps: [
        ...DIFFERENTIAL_RESULT,
        { run: page => expandAnalysisResult(page, 'differential') },
        { waitFor: '.de-modal-stats' },
        // The statistics pane is short by default; its splitter is dragged up,
        // as a reader would, until the whole table shows.
        { run: async page => {
          const box = await page.locator('.analysis-modal-resizer-horizontal').boundingBox();
          await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
          await page.mouse.down();
          await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2 - 300, { steps: 12 });
          await page.mouse.up();
          await page.mouse.move(4, 995);
          await page.waitForTimeout(500);
        } },
        { scrollIntoView: '.de-modal-stats' },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-stats-panel', padding: 0 },
    },

    // -------------------------------------------- the export row, any mode
    {
      id: 'analysis-export-toolbar',
      topic: 'analysis',
      name: 'export-toolbar',
      describes:
        'The export row of an expanded analysis result: an Export label '
        + 'followed by PNG, SVG and CSV buttons.',
      dataset: 'pancreas',
      steps: [...DETAILED_CATEGORICAL_RESULT, { run: page => expandAnalysisResult(page, 'detailed') }],
      shot: { mode: 'element', selector: '.analysis-modal-content .analysis-options-export-row', padding: 3 },
    },
  ];
}

/**
 * Set the expanded volcano's **Min |log2FC|** filter (the Differential
 * expression display option `minLog2FC`, shown in the expanded view's options).
 *
 * The field commits on `change`, so both `input` and `change` are sent. The
 * wait is on the counts the legend prints rather than on the field's own
 * value, because the field holds the number that was asked for whether or not
 * the plot ever redrew.
 *
 * @param {import('@playwright/test').Page} page
 * @param {number} threshold
 */
async function setVolcanoFoldChangeThreshold(page, threshold) {
  const before = await page.evaluate(() => (
    [...document.querySelectorAll('.analysis-modal-plot .legend text')]
      .map(node => node.textContent).join('|')
  ));
  await page.evaluate(value => {
    const field = document.querySelector('.analysis-modal-options input[name="minLog2FC"]');
    if (field === null) {
      throw new Error('The expanded volcano has no Min |log2FC| field.');
    }
    field.valueAsNumber = value;
    field.dispatchEvent(new Event('input', { bubbles: true }));
    field.dispatchEvent(new Event('change', { bubbles: true }));
  }, threshold);
  await page.waitForFunction(previous => (
    [...document.querySelectorAll('.analysis-modal-plot .legend text')]
      .map(node => node.textContent).join('|') !== previous
  ), before, { timeout: 60_000 });
}

// ---------------------------------------------------------------------------
// The local-file refusals on `b_data_loading`.
//
// These are the only scenarios that hand the application a file from disk. The
// CSV renamed `.h5ad` is built here, because building it is the evidence that
// the refusal is the product's own and not a doctored screenshot; the
// current-schema `.h5ad` whose `obsm` holds only `X_pca` is committed beside
// `fixtures/generate-refusal-fixtures.py`, which regenerates it.
//
// Both figures on that page shipped empty — 776 px of background scatter under
// an `:alt:` describing a message that was never in the frame: the refusal was
// then a timed state that deleted itself four seconds after it appeared, and
// the capture was still working when it went. The web app now keeps an ERROR
// until it is dismissed, and `assertCropStillFramed` in `capture.mjs` is why a
// future variant of that failure would fail the run instead of shipping.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildDataLoadingRefusalScenarios() {
  // Two screenshots back to back. Everything the settle loop would otherwise
  // wait for has already been waited for by the steps below — the dataset is
  // loaded and counted, the canvas has stopped drawing, and the notification
  // center is empty until the refusal lands in it — so a long gap buys
  // nothing.
  const REFUSAL_SETTLE = Object.freeze({ attempts: 2, gapMs: 0 });

  // A refusal is reported with `notifications.fail(...)`, which turns the
  // loading row into an ERROR without naming a duration, and an ERROR stays
  // until it is dismissed (`DEFAULT_DURATIONS_MS` in the web app's
  // `notification-center.js`). The subject therefore holds still once it has
  // landed; what does not hold still is everything raised before it, so the
  // center is drained first (see `NOTIFICATIONS_DRAINED`). That is also what
  // makes the figure show one message rather than four.
  //
  // Selecting a file is a `change` on a hidden input, so the Session panel is
  // opened and scrolled to first: the notification is captured, but a reader
  // who reruns this should see the control that raised it.
  const AT_LOCAL_CONTROLS = [
    PANCREAS_READY,
    { openSection: '#session-section' },
    { scrollIntoView: '#user-data-open-btn' },
    NOTIFICATIONS_DRAINED,
  ];

  /**
   * Offer one file to the picker behind **Open… › File…** (it takes an
   * `.h5ad` or a `.zarr.zip`) and wait for the notification text.
   *
   * @param {(directory: string) => Promise<string>} write - writes the fixture
   * @param {string} text - the substring the notification must carry
   */
  const offerH5ad = (write, text) => ({
    run: async page => {
      const { mkdtemp } = await import('node:fs/promises');
      const { tmpdir } = await import('node:os');
      const { join } = await import('node:path');
      const directory = await mkdtemp(join(tmpdir(), 'cellucid-docs-'));
      const fixture = await write(directory);
      await page.locator('#user-data-file-input').setInputFiles(fixture);
      await page.waitForFunction(
        needle => (
          document.querySelector('#notification-center')?.textContent
            ?.includes(needle) === true
        ),
        text,
        { timeout: 120_000 }
      );
    },
  });

  return [
    // -------------------------------------------- a CSV renamed to `.h5ad`
    {
      id: 'loading-refusal-h5ad-wrong-file',
      topic: 'data_loading',
      name: 'fail-h5ad-wrong-file',
      describes:
        'A Cellucid notification reading "The selected file is not a valid '
        + 'HDF5 file, so it cannot be an AnnData (.h5ad) file. Save it again '
        + 'with anndata…", with a dismiss button at its right.',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        offerH5ad(
          async directory => {
            // A real CSV, not random bytes. The extension is what the reader
            // got wrong; the leading bytes are what Cellucid checks.
            const { writeFile } = await import('node:fs/promises');
            const { join } = await import('node:path');
            const decoy = join(directory, 'not-really.h5ad');
            await writeFile(
              decoy,
              'cell_id,cell_type\nAAACCTG-1,T cell\nAAACGGG-1,B cell\n',
              'utf8'
            );
            return decoy;
          },
          'not a valid HDF5 file, so it cannot be an AnnData (.h5ad) file'
        ),
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: REFUSAL_SETTLE,
    },

    // ------------------------------------------ a valid `.h5ad` with no UMAP
    // The file opens on a best guess. The notification stays until it is
    // dismissed, so it is photographed on its own, Show button included.
    {
      id: 'loading-best-guess-notification',
      topic: 'data_loading',
      name: 'best-guess-notification',
      describes:
        'A Cellucid notification reading "This file has no default UMAP '
        + '(X_umap_1d, X_umap_2d, X_umap_3d, X_umap, X_UMAP, X_wnn.umap or X_ref.umap). '
        + 'Showing \'X_pca\' as 2D. Change '
        + 'it under Session › Data mapping.", with a Show button under the '
        + 'text and a dismiss button at its right.',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        offerH5ad(pcaOnlyFixture, 'This file has no default UMAP'),
      ],
      shot: {
        mode: 'element',
        selector: '#notification-center .notification:has(.notification-action)',
        padding: 10,
      },
      settle: REFUSAL_SETTLE,
    },

    // ------------------------------------------------- Open… with its menu open
    // A browser's file dialog picks files or one folder, never both, so the
    // one button asks first. The menu is a sibling list under the button, so
    // the crop spans both.
    {
      id: 'loading-open-menu',
      topic: 'data_loading',
      name: 'choose-local-data-control',
      describes:
        'Close-up of File or folder: the full-width Open… button with its menu '
        + 'open, showing File… (AnnData, Seurat/SCE .rds/.rda, Loom, 10x .h5, '
        + '.h5mu) and Folder… (Zarr, 10x outs, MatrixMarket), with the pointer on'
        + ' File….',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        { scrollIntoView: '#user-data-block' },
        { click: '#user-data-open-btn' },
        { waitFor: '#user-data-folder-item' },
      ],
      shot: {
        mode: 'region',
        selectors: ['#user-data-block', '#user-data-open-menu'],
        padding: 8,
      },
      cursor: { on: '#user-data-file-item', at: 'left', dx: 44, state: 'hover' },
      settle: REFUSAL_SETTLE,
    },

    // ---------------------------- the best guess's Show, opening the checklist
    // Data mapping never opens on its own; Show on the best-guess
    // notification is what opens its Embeddings checklist.
    {
      id: 'loading-best-guess-show',
      topic: 'data_loading',
      name: 'best-guess-show-checklist',
      describes:
        'The Data mapping group after Show on the best-guess notification: '
        + 'Embeddings open and summarised as 1 on · PCA 2D · guessed, its one '
        + 'row, PCA in 2D, on and marked Default and guessed, above the collapsed'
        + ' Genes, Neighbour graph, Latent space and Vector fields groups and the'
        + ' Revert and Apply buttons.',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        offerH5ad(pcaOnlyFixture, 'This file has no default UMAP'),
        { click: '#notification-center .notification-action' },
        { waitFor: '#data-mapping-embeddings' },
        { scrollIntoView: '#data-mapping-controls' },
      ],
      shot: { mode: 'element', selector: '#data-mapping-controls', padding: 8 },
      settle: REFUSAL_SETTLE,
    },

    // ------------------------------- the Embeddings list and "Other arrays"
    // Arrays that are off and wider than three columns are listed last,
    // collapsed. Garcia's 24-column `scvi` (its latent space) is one; the
    // capture fails at the waitForText below if the sample's own session
    // turns it on, rather than photographing an empty disclosure.
    {
      id: 'loading-data-mapping-other-arrays',
      topic: 'data_loading',
      name: 'data-mapping-other-arrays',
      describes:
        'The Data mapping group of the Garcia sample with its Embeddings '
        + 'checklist open as one flat list of rows, each with its checkbox, name,'
        + ' dimension and Default mark, and at its end Other arrays (1) open, '
        + 'listing scVI 1,2, which is off, with its width, 24 columns.',
      dataset: 'garcia',
      steps: [
        { run: async page => {
          await page.waitForFunction(() => (
            document.querySelector('#filter-count')?.textContent
              ?.includes('Showing all 219,731 cells') === true
          ), undefined, { timeout: 180_000 });
          await page.waitForFunction(
            () => document.documentElement.dataset.startupPhase === 'ready',
            undefined,
            { timeout: 180_000 }
          );
        } },
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
        { run: page => page.locator('#data-mapping-disclosure .disclosure-header').first().click() },
        { waitFor: '#data-mapping-embeddings' },
        { run: page => page.locator('#data-mapping-more .disclosure-header').first().click() },
        { waitForText: '#data-mapping-more', text: '24 columns' },
        { scrollIntoView: '#data-mapping-more' },
      ],
      shot: { mode: 'region', selectors: ['#data-mapping-controls .control-group-header', '#data-mapping-disclosure'], padding: 6 },
      settle: REFUSAL_SETTLE,
    },

    // --------------------------- the outlier filter without a latent space
    // `pca-only.h5ad` holds one two-column array, so nothing in it looks
    // like a latent space: the filter is its label, its (i) and one button.
    {
      id: 'filtering-outlier-no-latent',
      topic: 'filtering',
      name: 'outlier-no-latent',
      describes:
        'The Outlier filter in the legend for a file without a latent space: '
        + 'the label, its (i) and one Choose latent space button, with no '
        + 'slider, and the pointer on the button.',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        offerH5ad(pcaOnlyFixture, 'This file has no default UMAP'),
        { openSection: '#coloring-filtering-section' },
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'cell_type' },
        { waitFor: '#outlier-filter-action' },
        { waitForText: '#outlier-filter-action', text: 'Choose latent space' },
        { scrollIntoView: '#outlier-filter-container' },
      ],
      shot: { mode: 'element', selector: '#outlier-filter-container', padding: 8 },
      cursor: { on: '#outlier-filter-action', at: 'right', dx: -16, state: 'hover' },
      settle: REFUSAL_SETTLE,
    },
  ];
}

/**
 * The committed `.h5ad` whose `obsm` holds only a two-column `X_pca`
 * (`fixtures/generate-refusal-fixtures.py`).
 *
 * @returns {Promise<string>}
 */
async function pcaOnlyFixture() {
  const { fileURLToPath } = await import('node:url');
  const { dirname, join } = await import('node:path');
  return join(dirname(fileURLToPath(import.meta.url)), 'fixtures', 'pca-only.h5ad');
}

// ---------------------------------------------------------------------------
// Opening data from a link: b_data_loading/14_open_from_link and the pages
// that show the Session panel.
//
// A catalog link points at the docs catalog host (the published sample
// catalog, `datasets.json` and the AnnData stores it lists, read as a linked
// folder); a single-file link points at the web repository's permissive-CORS
// mirror, which answers range requests the way a bucket or a repository's file
// host does.

function buildDataLinkScenarios() {
  const SETTLE = Object.freeze({ attempts: 2, gapMs: 0 });
  const startupReady = {
    run: page => page.waitForFunction(
      () => document.documentElement.dataset.startupPhase === 'ready',
      undefined,
      { timeout: 120_000 }
    ),
  };
  const AT_ADDRESS_BAR = [
    { openSection: '#session-section' },
    { scrollIntoView: '#data-address-block' },
  ];
  /** An `.h5ad` read in place from the permissive-CORS mirror. */
  /**
   * The Pancreas sample as an `.h5ad`, read in place from a loopback host of
   * this run that answers range requests as a bucket does.
   */
  const linkedFile = async () => {
    const { dirname, basename } = await import('node:path');
    const file = await formatFixture('pancreas.h5ad');
    return { data: `${await serveFolder(dirname(file))}/${basename(file)}` };
  };
  const LINKED_FILE_OPEN = [
    startupReady,
    { run: page => page.waitForFunction(
      () => document.getElementById('dataset-info')?.dataset.state === 'ready'
    ) },
    // The read line counts on while the first views are read; the box is
    // photographed once it says what was read, not "Reading…".
    { run: page => page.waitForFunction(() => {
      const line = document.getElementById('dataset-read');
      return line !== null && !line.hidden && !/Reading/.test(line.textContent);
    }, undefined, { timeout: 120_000 }) },
  ];

  return [
    {
      id: 'loading-session-panel',
      topic: 'data_loading',
      name: 'data-loading-session-panel',
      describes:
        'The Session panel with the Pancreatic endocrinogenesis sample loaded: '
        + 'the Dataset box with its name, cells and genes, format and host, '
        + 'description and a collapsed Details; Open data with the Sample '
        + 'dataset dropdown, the Open… button and the Link or server bar; Save '
        + 'and Load under Session file; and the collapsed Data mapping group.',
      dataset: 'pancreas',
      steps: [PANCREAS_READY, NOTIFICATIONS_DRAINED, { openSection: '#session-section' }],
      shot: { mode: 'element', selector: '#session-section', padding: 6 },
    },
    {
      id: 'loading-session-panel-dark',
      topic: 'data_loading',
      name: 'data-loading-session-panel-dark',
      describes:
        'The same Session panel in the dark theme: the Dataset box, Open data '
        + 'with the Sample dataset dropdown, Open… and the Link or server bar, '
        + 'Session file with Save and Load, and the collapsed Data mapping group.',
      dataset: 'pancreas',
      theme: 'dark',
      steps: [
        PANCREAS_READY,
        { openSection: '#visualization-section' },
        { select: '#theme-select', value: 'dark' },
        { run: page => page.waitForFunction(() => document.documentElement.dataset.theme === 'dark') },
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
      ],
      shot: { mode: 'element', selector: '#session-section', padding: 6 },
    },
    {
      id: 'loading-dataset-box-details',
      topic: 'data_loading',
      name: 'dataset-box-details',
      describes:
        'The Dataset box after an .h5ad file was opened from a web address: its '
        + 'name, the cells and genes line, "h5ad file" and the host, how much '
        + 'has been read in how many requests, and Details open, listing '
        + 'File, Format, Size, Read, Cells, Genes, Obs fields, Embeddings and '
        + 'the rest.',
      query: linkedFile,
      steps: [
        ...LINKED_FILE_OPEN,
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
        { waitFor: '#dataset-read' },
        { click: '#dataset-details .disclosure-header' },
        { waitFor: '#dataset-details-list' },
        { scrollIntoView: '#dataset-facts' },
      ],
      shot: { mode: 'element', selector: '#dataset-facts', padding: 6 },
      settle: SETTLE,
    },
    {
      id: 'loading-link-control-typed',
      topic: 'data_loading',
      name: 'open-from-link-control',
      describes:
        'The Link or server bar with "zenodo:23288289" typed into it, the '
        + 'pointer on its Open button.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        { fill: '#data-address', value: 'zenodo:23288289' },
      ],
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      cursor: { on: '#data-address-open-btn', state: 'hover' },
    },
    {
      // A catalog of several datasets lists them in the dataset dropdown,
      // whose label then names the catalog's host, so the crop is the whole
      // Open data group: the dropdown above, the Link line below.
      id: 'loading-link-catalog-connected',
      topic: 'data_loading',
      name: 'connect-data-link',
      describes:
        'The Open data group connected to a dataset catalog folder: the '
        + 'dropdown labelled "Dataset from" the catalog\'s host and holding its '
        + 'datasets, Open…, and the Link or server bar holding the address, '
        + 'with one line under it reading Link, the address and a compact '
        + 'Disconnect button.',
      dataset: 'pancreas',
      query: urls => ({ data: `${urls.catalogOrigin}/` }),
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        ...AT_ADDRESS_BAR,
        { waitFor: '#data-link-disconnect-container' },
        { waitForText: 'label[for="dataset-select"]', text: 'Dataset from' },
      ],
      shot: { mode: 'element', selector: '#open-data-controls', padding: 8 },
      settle: SETTLE,
    },
    {
      // `cellucid serve` on the Pancreas sample's own store, on a free
      // loopback port, opened by typing its address. The server is the
      // product's own; it is stopped when the capture process exits.
      id: 'loading-address-server-connected',
      topic: 'data_loading',
      name: 'connect-remote-server',
      describes:
        'The Link or server bar holding the address of a cellucid serve '
        + 'server on this computer, with one line under it reading Server, the '
        + 'server\'s address and a compact Disconnect button.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        ...AT_ADDRESS_BAR,
        { run: async page => {
          const address = await serveSampleStore('pancreas.zarr');
          await page.fill('#data-address', address);
          await page.click('#data-address-open-btn');
          await page.waitForSelector('#remote-disconnect-container', { state: 'visible', timeout: 120_000 });
        } },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      settle: SETTLE,
    },
    {
      id: 'loading-link-refusal-unrecognised',
      topic: 'data_loading',
      name: 'fail-data-link-unrecognised',
      describes:
        'A Cellucid notification reading "Cellucid does not recognise "not a '
        + 'link"", followed by the kinds of link it opens.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        NOTIFICATIONS_DRAINED,
        { fill: '#data-address', value: 'not a link' },
        { click: '#data-address-open-btn' },
        { waitForText: '#notification-center', text: 'does not recognise' },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: SETTLE,
    },
  ];
}

/**
 * Serve one store of the sample catalog with `cellucid serve` on a free
 * loopback port and return the address to type into the bar
 * (`127.0.0.1:<port>`) once its health endpoint answers. The executable is
 * `$CELLUCID_CLI`, else `cellucid` on the PATH; the catalog is
 * `$CELLUCID_DATASETS_REPO/anndata`, as in `capture.mjs`. The server is
 * stopped when this process exits.
 *
 * @param {string} storeName - a store in the catalog folder
 * @returns {Promise<string>}
 */
async function serveSampleStore(storeName) {
  const { spawn } = await import('node:child_process');
  const { fileURLToPath } = await import('node:url');
  const { dirname, join, resolve } = await import('node:path');
  const here = dirname(fileURLToPath(import.meta.url));
  const datasets = process.env.CELLUCID_DATASETS_REPO
    ?? resolve(here, '..', '..', '..', '..', 'cellucid-datasets');
  const port = await freePort(`cellucid serve ${storeName}`);
  const child = spawn(process.env.CELLUCID_CLI ?? 'cellucid', [
    'serve', join(datasets, 'anndata', storeName),
    '--host', '127.0.0.1', '--port', String(port),
    '--no-browser', '--no-web-ui', '--quiet',
  ], { stdio: 'ignore' });
  let spawnError = null;
  child.once('error', error => { spawnError = error; });
  process.once('exit', () => child.kill());
  const deadline = Date.now() + 60_000;
  for (;;) {
    if (spawnError !== null) {
      throw new Error(`cellucid serve could not start: ${spawnError.message}`);
    }
    if (child.exitCode !== null) {
      throw new Error(`cellucid serve exited with status ${child.exitCode}`);
    }
    try {
      const health = await fetch(`http://127.0.0.1:${port}/_cellucid/health`);
      if (health.ok) return `127.0.0.1:${port}`;
    } catch {
      // Not listening yet.
    }
    if (Date.now() > deadline) {
      child.kill();
      throw new Error(`cellucid serve did not answer on port ${port} within 60 s`);
    }
    await new Promise(next => setTimeout(next, 250));
  }
}

// ---------------------------------------------------------------------------
// Highlighting and selection: f_highlighting_selection.
//
// Every figure is a real gesture on the Pancreas sample, made the way the
// pages tell a reader to make it: a tool button, then `Alt` held while the
// mouse presses, drags and releases over the plot. The 2D embedding is used
// because a flat projection makes the geometry of a gesture legible, and the
// gestures are fixed CSS coordinates in the fixed 1440x1000 viewport over that
// fixed view. A coordinate is only as good as the camera it was measured
// against, so every gesture that completes asserts the step readout the panel
// writes for it (`Step 1: 826 cells`, `Step 1: 591 cells`, …) and every
// gesture that is photographed mid-drag asserts the drag class the renderer
// puts on the canvas. A moved default camera therefore fails the capture
// instead of photographing a different selection under an unchanged caption.
// ---------------------------------------------------------------------------

/**
 * The one lasso every lasso figure draws, as viewport CSS points: an ellipse
 * of 36 segments around the upper part of the left lobe of the 2D embedding
 * (mostly Ductal cells), starting at its lower-left and running clockwise.
 * The whole loop encloses 826 cells.
 *
 * @param {number} fraction how much of the loop has been drawn, 0..1
 * @returns {Array<{x: number, y: number}>}
 */
function highlightLassoPath(fraction) {
  const centerX = 535;
  const centerY = 400;
  const radiusX = 95;
  const radiusY = 140;
  const segments = 36;
  const start = 0.75 * Math.PI;
  const points = [];
  for (let index = 0; index <= Math.round(segments * fraction); index += 1) {
    const angle = start + (2 * Math.PI * index) / segments;
    points.push({
      x: Math.round(centerX + radiusX * Math.cos(angle)),
      y: Math.round(centerY + radiusY * Math.sin(angle)),
    });
  }
  return points;
}

/**
 * A straight drag as twelve evenly spaced pointer positions.
 *
 * @param {{x: number, y: number}} from
 * @param {{x: number, y: number}} to
 */
function highlightDragPath(from, to) {
  const points = [];
  for (let index = 0; index <= 12; index += 1) {
    points.push({
      x: Math.round(from.x + ((to.x - from.x) * index) / 12),
      y: Math.round(from.y + ((to.y - from.y) * index) / 12),
    });
  }
  return points;
}

/**
 * Hold `Alt`, press at the first point and move through the rest — the
 * gesture every highlight tool reads. `release: false` leaves the button and
 * `Alt` held, which is how a figure shows a gesture in progress: the preview
 * the renderer paints while dragging exists only while the button is down.
 *
 * @param {import('@playwright/test').Page} page
 * @param {Array<{x: number, y: number}>} path
 * @param {{release: boolean}} options
 */
async function altGesture(page, path, { release }) {
  const [first, ...rest] = path;
  await page.mouse.move(first.x, first.y);
  await page.keyboard.down('Alt');
  await page.mouse.down();
  for (const point of rest) {
    await page.mouse.move(point.x, point.y, { steps: 3 });
  }
  if (release) {
    await page.mouse.up();
    await page.keyboard.up('Alt');
  }
}

/**
 * Draw the pointer glyph where a canvas gesture is, plus any rings.
 *
 * `capture.mjs`'s own `cursor` field anchors to an element and, for `press`,
 * issues a real `mouse.down()`. A gesture on the canvas is neither: its point
 * is a coordinate in the plot, and the button is either already held (a
 * second press would restart the gesture) or already released (a fresh
 * press would start a new one). So the glyph is painted here, with the same
 * overlay code, at the coordinate the gesture actually used, and the real
 * pointer is left exactly where the gesture left it.
 *
 * @param {import('@playwright/test').Page} page
 * @param {{x: number, y: number}} point
 * @param {Array<{on: string, label?: string}>} [rings]
 */
async function paintGesturePointer(page, point, rings = []) {
  const ringRects = [];
  for (const ring of rings) {
    const box = await page.locator(ring.on).first().boundingBox();
    if (box === null || box.width <= 0 || box.height <= 0) {
      throw new Error(`Ring selector has no painted area: ${ring.on}`);
    }
    ringRects.push({ ...box, label: ring.label ?? null });
  }
  await page.evaluate(paintOverlay, {
    rootId: OVERLAY_ROOT_ID,
    cursor: { x: point.x, y: point.y, state: 'press' },
    rings: ringRects,
  });
}

/**
 * Pin the highlight glow, the one animation a highlight figure cannot settle.
 *
 * Highlighted and preview cells are drawn with a glow that pulses on the
 * shader uniform `u_time` (`0.5 + 0.5 * sin(u_time * 4.0)` in the web app's
 * `high-perf-shaders.js`), fed from `performance.now()` minus the moment the
 * highlight renderer was built. No two frames agree, so the settle loop never
 * converges, and a frame-pinned capture is no better: the phase depends on
 * how long startup took on the capture host. So the two uniforms the pulse
 * is drawn from are held at time zero: `u_time` at 0 and `u_glow` at its value
 * there, `0.5 + 0.5 · sin(0)` = 0.5 (`highlightGlowAt` in the web app's
 * `high-perf-shaders.js`) — the picture the renderer itself draws for a still
 * highlight style. It is the same kind of intervention as the seeded
 * `Math.random`: the clock the app reads for everything else, its timers, and
 * every other uniform are untouched.
 *
 * @param {{context: import('@playwright/test').BrowserContext}} options
 */
function pinHighlightGlow({ context }) {
  return context.addInitScript(pinShaderTimeInPage);
}

/** The page-side half of `pinHighlightGlow`; serialised, closes over nothing. */
function pinShaderTimeInPage() {
  const timeLocations = new WeakSet();
  const glowLocations = new WeakSet();
  for (const Context of [window.WebGLRenderingContext, window.WebGL2RenderingContext]) {
    if (Context === undefined) continue;
    const getUniformLocation = Context.prototype.getUniformLocation;
    const uniform1f = Context.prototype.uniform1f;
    Context.prototype.getUniformLocation = function (program, name) {
      const location = getUniformLocation.call(this, program, name);
      if (location !== null && name === 'u_time') timeLocations.add(location);
      if (location !== null && name === 'u_glow') glowLocations.add(location);
      return location;
    };
    Context.prototype.uniform1f = function (location, value) {
      if (timeLocations.has(location)) return uniform1f.call(this, location, 0);
      if (glowLocations.has(location)) return uniform1f.call(this, location, 0.5);
      return uniform1f.call(this, location, value);
    };
  }
}

/** @returns {Array<Object>} */
function buildHighlightingScenarios() {
  // Open on first load and not the subject; collapsed so the panel fits.
  // Declared here, not beside the helpers: `SCENARIOS` calls this builder
  // while the module is still initialising, before a later top-level
  // `const` exists.
  const NON_HIGHLIGHT_SECTIONS = [
    '#session-section',
    '#visualization-section',
    '#compare-views-section',
    '#coloring-filtering-section',
    '#page-analysis-section',
  ];
  const PANEL = { mode: 'element', selector: '#highlighted-cells-section', padding: 6 };
  const WINDOW = { mode: 'window' };
  const settle = { attempts: 20, gapMs: 220 };
  const readout = text => ({ waitForText: '#highlight-mode-description', text });

  const LASSO_END = highlightLassoPath(0.8).at(-1);
  const BETA_CELL = { x: 1060, y: 400 };
  const EMPTY_SPACE = { x: 860, y: 250 };
  const SEED = { x: 535, y: 440 };
  // Proximity drags outward and KNN drags up, by the same 112 px, from the
  // same seed cell, so the two figures compare like with like.
  const PROXIMITY_END = { x: 635, y: 490 };
  const KNN_END = { x: 535, y: 328 };

  /**
   * Pancreas, flat 2D view in Planar navigation, the given colouring, and
   * every other open section collapsed so Highlighting sits whole in the
   * window. Showing UMAP · 2D moves the view to Planar; waiting for that
   * and then for the notification center to drain is what makes the view
   * final before the first coordinate is used.
   */
  const ready2d = field => [
    PANCREAS_READY,
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label: field },
    {
      waitForText: '#stats',
      text: field === 'None' ? 'Field: None' : `Field: ${field} (category)`,
    },
    { openSection: '#compare-views-section' },
    { waitForStable: '#view-representation-select' },
    { select: '#view-representation-select', label: 'UMAP · 2D' },
    { run: page => page.waitForFunction(
      () => document.querySelector('#navigation-mode')?.value === 'planar'
    ) },
    ...NON_HIGHLIGHT_SECTIONS.map(selector => ({ closeSection: selector })),
    NOTIFICATIONS_DRAINED,
  ];
  const READY = ready2d('cell_type');

  const LASSO_RELEASED = [
    ...READY,
    { click: '#highlight-mode-lasso' },
    { run: page => altGesture(page, highlightLassoPath(1), { release: true }) },
    readout('Step 1: 826 cells'),
  ];
  const LASSO_CONFIRMED = [
    ...LASSO_RELEASED,
    { click: '#lasso-confirm-btn' },
    { waitForText: '#highlight-count', text: '826 cells highlighted' },
    // Off the panel, so no control keeps the hover styling of the click.
    { run: page => page.mouse.move(1300, 940) },
  ];
  const ANNOTATION_CLICKED = [
    ...READY,
    { run: page => altGesture(page, [BETA_CELL], { release: true }) },
    readout('Step 1: 591 cells'),
  ];

  return [
    {
      id: 'highlighting-choose-mode',
      topic: 'highlighting_selection',
      name: 'choose-highlight-mode',
      describes:
        'The Highlighting panel before anything is selected: the Selection (i) '
        + 'and the four Tool buttons — Annotation (active), KNN drag, Proximity '
        + 'drag and Lasso — ringed, with the pointer on Lasso; the Select mode '
        + 'switch, off; Highlights with one page tab reading Page 1 (0) and the '
        + 'empty list reading No cells highlighted; and the collapsed Create '
        + 'categorical.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [...READY, { scrollIntoView: '#highlighted-cells-section' }],
      shot: PANEL,
      cursor: { on: '#highlight-mode-lasso', at: 'center', state: 'hover' },
      rings: [{ on: '#highlight-mode', label: '1' }],
      settle,
    },
    {
      id: 'highlighting-lasso-draw',
      topic: 'highlighting_selection',
      name: '01-draw-lasso',
      describes:
        'The 2D planar Pancreas embedding mid-lasso: Alt is held, the pointer '
        + 'is pressed at the end of a dotted path drawn most of the way round the'
        + ' upper part of the left lobe, and the cells the path already encloses '
        + 'are drawn in the highlight colour while the panel shows no step yet.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...READY,
        { click: '#highlight-mode-lasso' },
        { run: page => altGesture(page, highlightLassoPath(0.8), { release: false }) },
        { waitFor: '#glcanvas.lassoing' },
        { run: page => paintGesturePointer(page, LASSO_END) },
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-lasso-step',
      topic: 'highlighting_selection',
      name: '02-review-lasso-step',
      describes:
        'The Highlighting panel after the lasso was released, showing the step '
        + 'readout "Step 1: 826 cells" with the help line "Alt to intersect, '
        + 'Shift+Alt to add, Ctrl+Alt to subtract", and the ringed Confirm, undo, '
        + 'redo and Cancel step controls with the pointer on Confirm; the page '
        + 'tab still reads Page 1 (0) and the list still reads "No cells '
        + 'highlighted".',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [...LASSO_RELEASED, { scrollIntoView: '#highlighted-cells-section' }],
      shot: PANEL,
      cursor: { on: '#lasso-confirm-btn', at: 'right', dx: -12, state: 'hover' },
      rings: [{ on: '#lasso-step-controls', label: '2' }],
      settle,
    },
    {
      id: 'highlighting-lasso-confirmed',
      topic: 'highlighting_selection',
      name: '03-confirm-lasso-group',
      describes:
        'The Highlighting panel after Confirm: the page tab reads Page 1 826, '
        + 'and the ringed group list shows the summary "826 cells highlighted", a '
        + 'Clear button, and one group row labelled "Lasso (826 cells) (826)" '
        + 'with a checkbox and a remove control.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [...LASSO_CONFIRMED, { scrollIntoView: '#highlighted-cells-section' }],
      shot: PANEL,
      rings: [{ on: '#highlighted-groups', label: '3' }],
      settle,
    },
    {
      id: 'highlighting-lasso-result-window',
      topic: 'highlighting_selection',
      name: '04-view-highlighted-cells',
      describes:
        'The whole window after confirming the lasso group: the 826 enclosed '
        + 'cells in the upper part of the left lobe are drawn in the highlight '
        + 'colour, and the Highlighting panel reports 826 cells highlighted with '
        + 'one Lasso group.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: LASSO_CONFIRMED,
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-annotation-click',
      topic: 'highlighting_selection',
      name: 'select-category-by-click',
      describes:
        'Annotation selection after Alt+click on a cell of the green Beta arm: '
        + 'all 591 Beta cells are drawn in the highlight colour, a tooltip by the'
        + ' pointer names Beta, 591 cells and cell_type, and the panel reads Step'
        + ' 1: 591 cells with the help line "Alt to replace, Shift+Alt to add, '
        + 'Ctrl+Alt to subtract" above the Confirm, undo, redo and Cancel '
        + 'controls.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...ANNOTATION_CLICKED,
        { run: page => paintGesturePointer(page, BETA_CELL) },
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-proximity-drag',
      topic: 'highlighting_selection',
      name: 'drag-proximity-radius',
      describes:
        'Proximity drag mid-gesture on the 2D embedding: a dashed circle '
        + 'centred on a dot at the seed cell in the left lobe, the cells inside '
        + 'it drawn in the highlight colour, and the pointer pressed on the '
        + 'circle\'s edge.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...READY,
        { click: '#highlight-mode-proximity' },
        {
          run: page => altGesture(
            page, highlightDragPath(SEED, PROXIMITY_END), { release: false }
          ),
        },
        { waitFor: '#glcanvas.proximity-dragging' },
        { run: page => paintGesturePointer(page, PROXIMITY_END) },
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-knn-drag',
      topic: 'highlighting_selection',
      name: 'drag-knn-degree',
      describes:
        'KNN drag mid-gesture from the same seed cell, dragged upward by the '
        + 'same distance: a small dashed circle labelled 3° around the seed, and '
        + 'the highlighted cells spreading through the left lobe along the '
        + 'neighbour graph rather than filling a disc, with the pointer pressed '
        + 'above the seed.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...READY,
        { click: '#highlight-mode-knn' },
        // Entering KNN mode loads the neighbour graph if it is not loaded yet,
        // under a "Loading neighbor graph" notification, and a press made
        // before it lands only starts that load. Draining the notification
        // center waits for the load (and its completion notice) to finish.
        NOTIFICATIONS_DRAINED,
        {
          run: page => altGesture(
            page, highlightDragPath(SEED, KNN_END), { release: false }
          ),
        },
        { waitFor: '#glcanvas.knn-dragging' },
        { run: page => paintGesturePointer(page, KNN_END) },
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-two-pages',
      topic: 'highlighting_selection',
      name: 'compare-highlight-pages',
      describes:
        'The Highlighting panel with two page tabs ringed, their names '
        + 'shortened — Page 1 with a blue swatch carrying 826, and the active '
        + 'Page 2 with a red swatch carrying 267, each with a combine (⋈) and a '
        + 'delete (×) control — the list reading 267 cells highlighted with one '
        + 'Proximity (267 cells) (267) group, and the pointer on the button that '
        + 'adds another page.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...LASSO_CONFIRMED,
        { click: '#add-highlight-page' },
        { waitForText: '#highlight-pages-tabs', text: 'Page 2' },
        { click: '#highlight-mode-proximity' },
        {
          run: page => altGesture(
            page,
            highlightDragPath({ x: 1090, y: 640 }, { x: 1140, y: 670 }),
            { release: true }
          ),
        },
        readout('Step 1: 267 cells'),
        { click: '#proximity-confirm-btn' },
        { waitForText: '#highlight-count', text: '267 cells highlighted' },
        { scrollIntoView: '#highlighted-cells-section' },
      ],
      shot: PANEL,
      cursor: { on: '#add-highlight-page', at: 'center', state: 'hover' },
      rings: [{ on: '#highlight-pages-tabs', label: '1' }],
      settle,
    },
    {
      id: 'highlighting-across-views',
      topic: 'highlighting_selection',
      name: 'share-highlight-across-views',
      describes:
        'Two 3D panels of the same dataset side by side with the sidebar '
        + 'closed, each under a title chip (1 · UMAP 3D · cell_type, 2 · UMAP 3D '
        + '· cell_type): the left panel at the default close-up camera, the right'
        + ' panel zoomed out and rotated; the cells confirmed by one lasso drawn '
        + 'in the left panel are drawn in the highlight colour in both.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        PANCREAS_READY,
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'cell_type' },
        { waitForText: '#stats', text: 'Field: cell_type (category)' },
        { openSection: '#compare-views-section' },
        { click: '#split-keep-view-btn' },
        {
          run: page => page.waitForFunction(() => (
            document.querySelectorAll('#split-view-badges-list > .split-badge').length === 2
          )),
        },
        // Independent cameras are what let the second panel look from
        // another angle; unlink them if the kept view started linked.
        { run: async page => {
          if (await page.getAttribute('#camera-lock-btn', 'aria-pressed') === 'true') {
            await page.click('#camera-lock-btn');
          }
        } },
        { waitFor: '#camera-lock-btn[aria-pressed="false"]' },
        { click: '#highlight-mode-lasso' },
        { click: '#sidebar-toggle' },
        { waitForText: '#sidebar-toggle', text: '☰' },
        NOTIFICATIONS_DRAINED,
        {
          // Zoom out and orbit the focused right panel only. The pointer
          // rests before release so no drag velocity is carried past it.
          run: async page => {
            await page.mouse.move(1080, 500);
            for (let notch = 0; notch < 5; notch += 1) {
              await page.mouse.wheel(0, 300);
              await page.waitForTimeout(80);
            }
            await page.mouse.down();
            await page.mouse.move(1260, 420, { steps: 15 });
            await page.waitForTimeout(300);
            await page.mouse.up();
          },
        },
        // One lasso, in the left panel, around the Ductal flap.
        {
          run: page => altGesture(page, [
            { x: 50, y: 370 }, { x: 200, y: 300 }, { x: 330, y: 285 },
            { x: 425, y: 330 }, { x: 420, y: 420 }, { x: 300, y: 475 },
            { x: 180, y: 515 }, { x: 60, y: 430 }, { x: 50, y: 370 },
          ], { release: true }),
        },
        readout('Step 1: 1,204 cells'),
        { click: '#sidebar-toggle' },
        { waitForText: '#sidebar-toggle', text: '✕' },
        { click: '#lasso-confirm-btn' },
        { waitForText: '#highlight-count', text: '1,204 cells highlighted' },
        { click: '#sidebar-toggle' },
        { waitForText: '#sidebar-toggle', text: '☰' },
        { run: page => page.mouse.move(720, 990) },
        NOTIFICATIONS_DRAINED,
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-missed-click-notice',
      topic: 'highlighting_selection',
      name: 'read-missed-click-notice',
      describes:
        'A proximity Alt+click on empty space: the pointer is pressed above the'
        + ' embedding where no cell is drawn, and the ringed line in the '
        + 'Highlighting panel, under the tool buttons, reads "That click did not '
        + 'land on a cell, so there is nothing to select. Zoom in or raise Point '
        + 'size to make cells easier to hit: only the cells the viewer is '
        + 'currently drawing can be clicked."',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...READY,
        { click: '#highlight-mode-proximity' },
        { run: page => altGesture(page, [EMPTY_SPACE], { release: true }) },
        readout('That click did not land on a cell'),
        {
          run: page => paintGesturePointer(page, EMPTY_SPACE, [
            { on: '#highlight-mode-description', label: '1' },
          ]),
        },
      ],
      shot: WINDOW,
      settle,
    },
    {
      id: 'highlighting-no-active-field-notice',
      topic: 'highlighting_selection',
      name: 'read-no-active-field-notice',
      describes:
        'The Highlighting panel with the Annotation tool chosen and no field '
        + 'coloured, and the Selection (i) tip open, its first line ringed: '
        + '"Color by a field first (Coloring & filtering), then Alt+click a '
        + 'cell."',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...ready2d('None'),
        { scrollIntoView: '#highlighted-cells-section' },
        { click: '#highlight-selection-info-btn' },
        { waitForText: '#highlight-tool-hint', text: 'Color by a field first' },
      ],
      shot: { mode: 'region', selectors: ['#highlighted-cells-section', '#highlight-selection-info-tooltip'], padding: 6 },
      rings: [{ on: '#highlight-tool-hint', label: '1' }],
      settle,
    },
    {
      id: 'highlighting-annotation-step',
      topic: 'highlighting_selection',
      name: 'read-annotation-step',
      describes:
        'The Highlighting panel after an annotation Alt+click on a Beta cell, '
        + 'with the ringed step readout "Step 1: 591 cells" and the help line '
        + '"Alt to replace, Shift+Alt to add, Ctrl+Alt to subtract" above '
        + 'Confirm, undo, redo and Cancel.',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [...ANNOTATION_CLICKED, { scrollIntoView: '#highlighted-cells-section' }],
      shot: PANEL,
      rings: [{ on: '#highlight-mode-description', label: '1' }],
      settle,
    },
    {
      id: 'highlighting-named-page',
      topic: 'highlighting_selection',
      name: 'highlighting-selected-page',
      describes:
        'The Highlighting panel with one page renamed "Left lobe" carrying 826 '
        + 'cells, the summary "826 cells highlighted", and one group row labelled '
        + '"Lasso (826 cells) (826)".',
      dataset: 'pancreas',
      setup: pinHighlightGlow,
      steps: [
        ...LASSO_CONFIRMED,
        // Rename by double-clicking the tab name, as the UI page documents.
        { run: page => page.dblclick('.highlight-page-tab-name') },
        { fill: '.highlight-page-tab-input', value: 'Left lobe' },
        { press: 'Enter' },
        { waitForText: '#highlight-pages-tabs', text: 'Left lobe' },
        // Enter returns focus to the tab; a click on the summary line moves
        // it off, so the tab is drawn without a keyboard focus ring.
        { click: '#highlight-count' },
        { run: page => page.mouse.move(1300, 940) },
        { scrollIntoView: '#highlighted-cells-section' },
      ],
      shot: PANEL,
      settle,
    },
  ];
}


// ---------------------------------------------------------------------------
// Multiview: c_core_interactions and n_benchmarking_performance.
//
// Grid compare lays its panels out in the presented area — the part of the
// canvas the sidebar does not cover — as equal cells (`resolvePresentedLayout`
// in the web app's `rendering/viewer/presented-view-layout.js`). These figures
// exist to show that layout, so each one asserts it before the shutter rather
// than trusting it: every panel's title chip is centred on the panel it
// labels, so the chips' centres pin down where the panels are.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildMultiviewScenarios() {
  const settle = { attempts: 20, gapMs: 220 };

  /** Colour the selected view and wait for the stats line to say so. */
  const colourBy = (select, field, kind) => [
    { waitForStable: select },
    { select, label: field },
    { waitForText: '#stats', text: `Field: ${field} (${kind})` },
  ];

  /**
   * Keep view: the live view becomes badge 1 and the new snapshot, badge 2,
   * is the selected view. The grid exists once both badges and both panel
   * title chips do.
   */
  const KEEP_VIEW = [
    { openSection: '#compare-views-section' },
    { click: '#split-keep-view-btn' },
    {
      run: page => page.waitForFunction(() => (
        document.querySelectorAll('#split-view-badges-list > .split-badge').length === 2
        && document.querySelectorAll('#view-title-layer .view-title-chip').length === 2
      )),
    },
  ];

  /**
   * The published sample's camera is framed for one full-window view; a panel
   * a little under half as wide shows the embedding too close. Four wheel
   * notches over the right panel zoom both panels out (the cameras are locked),
   * and the pointer then rests in the bottom-left corner of the window, clear of
   * the embedding and of every control, so nothing in the frame carries hover
   * styling.
   */
  const ZOOM_OUT = {
    run: async page => {
      await page.mouse.move(1150, 500);
      for (let notch = 0; notch < 4; notch += 1) {
        await page.mouse.wheel(0, 300);
        await page.waitForTimeout(80);
      }
      await page.mouse.move(150, 995);
    },
  };

  /**
   * The layout under test: two equal panels side by side, starting at the
   * sidebar's right edge when it is open and at the window's left edge when it
   * is closed. Waiting for it, rather than reading it once, also waits out the
   * sidebar's open/close slide.
   */
  const TWO_EQUAL_PANELS_IN_PRESENTED_AREA = {
    run: page => page.waitForFunction(() => {
      const sidebar = document.querySelector('#sidebar');
      const left = sidebar.classList.contains('hidden')
        ? 0
        : Math.ceil(sidebar.getBoundingClientRect().right);
      const chips = [...document.querySelectorAll('#view-title-layer .view-title-chip')];
      const panelWidth = (window.innerWidth - left) / chips.length;
      return chips.length === 2 && chips.every((chip, index) => {
        const rect = chip.getBoundingClientRect();
        const centre = rect.left + rect.width / 2;
        return Math.abs(centre - (left + panelWidth * (index + 0.5))) < 1;
      });
    }),
  };

  /** Scroll the sidebar so one section's header sits at its top edge. */
  const scrollSidebarTo = selector => ({
    run: page => page.locator(selector).evaluate(
      element => element.scrollIntoView({ block: 'start' })
    ),
  });

  return [
    {
      id: 'core-multiview-two-panels',
      topic: 'web_app',
      name: 'multiview-two-panels',
      describes:
        'The whole window in Grid compare with the sidebar open: two equal '
        + 'panels to the right of the sidebar show the Pancreas 3D embedding from'
        + ' the same linked camera — view 1 coloured by cell_type, view 2 '
        + 'recoloured by S_score on a Viridis scale — each under a title chip '
        + 'naming its embedding and field. The sidebar shows Compare views with '
        + 'the two view rows, each marked linked and Orbit, View 2 shows UMAP · '
        + '3D, Keep view, Link cameras, Layout Grid compare and Clear, and '
        + 'Coloring & filtering with S_score in the Continuous obs select.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        ...KEEP_VIEW,
        // The snapshot is the selected view, so this recolours it alone.
        ...colourBy('#continuous-field', 'S_score', 'continuous'),
        { waitForText: '#view-title-layer .view-title-chip.active', text: 'S_score' },
        ZOOM_OUT,
        TWO_EQUAL_PANELS_IN_PRESENTED_AREA,
        scrollSidebarTo('#compare-views-section'),
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'performance-dark-theme-multiview',
      topic: 'web_app',
      name: 'dark-theme-multiview',
      describes:
        'The whole window in the dark theme on the Grid (dark) background, in '
        + 'Grid compare with the sidebar open: two equal panels show the Pancreas'
        + ' 3D embedding coloured by cell_type from the same linked camera, each '
        + 'under a title chip. The sidebar shows Visualization with Theme Dark '
        + 'and Background Grid (dark), and Compare views below it with the two '
        + 'view rows.',
      dataset: 'pancreas',
      theme: 'dark',
      steps: [
        PANCREAS_READY,
        { openSection: '#visualization-section' },
        { select: '#theme-select', value: 'dark' },
        { run: page => page.waitForFunction(() => document.documentElement.dataset.theme === 'dark') },
        { select: '#background-select', value: 'grid-dark' },
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        ...KEEP_VIEW,
        ZOOM_OUT,
        TWO_EQUAL_PANELS_IN_PRESENTED_AREA,
        scrollSidebarTo('#visualization-section'),
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'core-connectivity-multiview',
      topic: 'web_app',
      name: 'connectivity-multiview',
      describes:
        'Two equal panels filling the whole window with the sidebar closed (the'
        + ' ☰ button at the top left), both showing the Pancreas 3D embedding '
        + 'coloured by cell_type from the same linked camera, each under a title '
        + 'chip, and in both the neighbour-graph edges, drawn at line width 2.0, '
        + 'form a grey web between the cells.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        ...KEEP_VIEW,
        // Turned on after the snapshot exists and while it is the selected
        // view, so edges in the live panel show the setting reaching it too.
        { openSection: '#visualization-section' },
        { run: page => page.check('#toggle-connectivity') },
        { waitForText: '#connectivity-info', text: '33k total · 33k visible · 33k shown' },
        // Thicker than the 1.0 default so the edges survive the page's
        // downscale; the alt text names the width.
        {
          run: async page => {
            await page.fill('#connectivity-width', '20');
            await page.dispatchEvent('#connectivity-width', 'change');
          },
        },
        { waitForText: '#connectivity-width-display', text: '2.0' },
        { click: '#sidebar-toggle' },
        { waitForText: '#sidebar-toggle', text: '☰' },
        ZOOM_OUT,
        TWO_EQUAL_PANELS_IN_PRESENTED_AREA,
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
  ];
}

// ---------------------------------------------------------------------------
// The documentation home page (docs/index.md).
//
// The home page places its images with its own `<img>` markup rather than a
// `{figure}`: each picture exists in a light and a dark variant, shown by the
// theme the reader picked, and the hero also at two widths for `srcset`. So
// these scenarios emit WebP at a fixed quality instead of PNG — a full-window
// capture of hundreds of thousands of cells is megabytes as PNG and a few
// hundred kilobytes as WebP — and they name their width:
//
//   hero      1280 and 2400 px: the page column is at most 70rem (1120 CSS
//             px), so 1280 serves a 1x screen and a phone, 2400 a 2x screen.
//   features  1080 px: a feature card is at most about 540 CSS px wide, so
//             1080 is 2x there and about 3x on a phone, where it spans the
//             screen.
//
// Every scene is the app's own look, on the backdrop each theme is shown on
// (Studio in light, Nebula in dark), with the published sample data.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildHomeScenarios() {
  const settle = { attempts: 24, gapMs: 250 };
  const WEBP = { format: 'webp', quality: 0.86 };
  const THEMES = [
    { theme: 'light', background: 'studio', backdropName: 'Studio' },
    { theme: 'dark', background: 'nebula', backdropName: 'Nebula' },
  ];

  const SUO_READY = {
    run: async page => {
      await page.waitForFunction(
        () => document.documentElement.dataset.startupPhase === 'ready'
      );
      await page.waitForFunction(() => (
        document.querySelector('#filter-count')?.textContent
          ?.includes('Showing all 561,947 cells') === true
      ));
    },
  };

  /** The app's theme and backdrop, set the way a reader sets them. */
  const look = ({ theme, background }) => [
    { openSection: '#visualization-section' },
    { waitForStable: '#theme-select' },
    { select: '#theme-select', value: theme },
    {
      run: page => page.waitForFunction(
        expected => document.documentElement.dataset.theme === expected,
        theme
      ),
    },
    { select: '#background-select', value: background },
  ];

  const colourBy = (select, field, kind) => [
    { waitForStable: select },
    { select, label: field },
    { waitForText: '#stats', text: `Field: ${field} (${kind})` },
  ];

  /** Only these accordion sections open, the first at the top of the sidebar. */
  const onlySections = (...open) => ({
    run: page => page.evaluate(wanted => {
      for (const section of document.querySelectorAll('details.accordion-section')) {
        section.open = wanted.includes(`#${section.id}`);
      }
      document.querySelector(wanted[0])?.scrollIntoView({ block: 'start' });
    }, open),
  });

  /**
   * Turn the published camera a little and come closer, so the 3D embedding
   * reads as a solid and fills the space beside the sidebar; then rest the
   * pointer below the sidebar, clear of the canvas (whose hover shows cell
   * info) and of every control.
   */
  const TURN_AND_CLOSER = {
    run: async page => {
      await page.mouse.move(900, 500);
      await page.mouse.down();
      await page.mouse.move(960, 470, { steps: 12 });
      await page.waitForTimeout(250);
      await page.mouse.up();
      await page.mouse.wheel(0, -200);
      await page.mouse.move(150, 995);
    },
  };

  /** Show the 2D UMAP in the view, as the Shows picker does. */
  const TO_2D = [
    { openSection: '#navigation-section' },
    { openSection: '#compare-views-section' },
    { waitFor: '#representation-controls:not([hidden])' },
    { waitForStable: '#view-representation-select' },
    // The picker shows the view's embedding, so it reads 3D again until the
    // switch has landed; the switch is done when the view and the navigation
    // mode both say 2D, which is what is waited for instead of the control.
    {
      run: async page => {
        await page.selectOption('#view-representation-select', { label: 'UMAP · 2D' });
        await page.waitForFunction(() => (
          document.querySelector('#navigation-mode')?.value === 'planar'
            && document.querySelector('#view-representation-select')?.value === 'X_umap'
        ));
      },
    },
  ];

  /** Two wheel notches closer over the canvas, then the pointer at rest. */
  const CLOSER_2D = {
    run: async page => {
      await page.mouse.move(880, 500);
      for (let notch = 0; notch < 2; notch += 1) {
        await page.mouse.wheel(0, -200);
        await page.waitForTimeout(80);
      }
      await page.mouse.move(150, 995);
    },
  };

  const CLOSE_SIDEBAR = [
    { click: '#sidebar-toggle' },
    { waitForText: '#sidebar-toggle', text: '☰' },
  ];

  const hero = (variant, width) => ({
    id: `home-hero-${variant.theme}-${width}`,
    topic: 'home',
    name: `hero-${variant.theme}-${width}`,
    describes:
      `Cellucid in the ${variant.theme} theme on the ${variant.backdropName} `
      + 'backdrop: a 3D UMAP of the developing human immune system atlas '
      + '(Suo), all 561,947 cells coloured by cell type, beside the Coloring '
      + 'and filtering panel, whose legend lists each cell type with its '
      + 'number of cells.',
    dataset: 'suo',
    theme: variant.theme,
    steps: [
      SUO_READY,
      ...look(variant),
      ...colourBy('#categorical-field', 'cell_type', 'category'),
      onlySections('#coloring-filtering-section', '#highlighted-cells-section'),
      TURN_AND_CLOSER,
      NOTIFICATIONS_DRAINED,
    ],
    shot: { mode: 'window' },
    settle,
    emitWidth: width,
    ...WEBP,
  });

  const feature = (variant, { key, describes, dataset, steps, shot, settle: own }) => ({
    id: `home-feature-${key}-${variant.theme}`,
    topic: 'home',
    name: `feature-${key}-${variant.theme}`,
    describes: describes(variant),
    dataset,
    theme: variant.theme,
    steps: steps(variant),
    shot: shot ?? { mode: 'window' },
    settle: own ?? settle,
    emitWidth: 1080,
    ...WEBP,
  });

  const FEATURES = [
    {
      key: 'genes',
      dataset: 'pancreas',
      describes: variant =>
        `The Pancreas 2D UMAP in the ${variant.theme} theme, coloured by `
        + 'expression of the gene Ttr on a log colour scale (Viridis): '
        + 'brightest in the endocrine cells on the right, a darker gradient '
        + 'along the trajectory, and pale grey for the few cells with none; '
        + 'the colour bar is open in the sidebar.',
      steps: variant => [
        PANCREAS_READY,
        ...look(variant),
        ...TO_2D,
        onlySections('#coloring-filtering-section'),
        { waitForStable: '#gene-expression-search' },
        { fill: '#gene-expression-search', value: 'Ttr' },
        { waitFor: '#gene-expression-dropdown .dropdown-item' },
        {
          run: page =>
            page.locator('#gene-expression-dropdown .dropdown-item').first().click(),
        },
        { waitForText: '#stats', text: 'Field: Ttr (continuous)' },
        { run: page => page.locator('#legend button.switch:has-text("Log color scale")').click() },
        {
          run: page => page.waitForFunction(() => [
            ...document.querySelectorAll('#legend button.switch[aria-checked="true"]'),
          ].some(node => node.textContent.includes('Log color scale'))),
        },
        CLOSER_2D,
        NOTIFICATIONS_DRAINED,
      ],
    },
    {
      key: 'analysis',
      dataset: 'pancreas',
      describes: variant =>
        `The expanded Marker genes window for cell_type in the ${variant.theme} `
        + 'theme, over the blurred Pancreas embedding: a heatmap of z-scored '
        + 'expression with the eight cell types down the side and their top '
        + 'marker genes across, Export and Plot options beside it, and below, '
        + 'the summary statistics and a table of the top genes of the '
        + 'selected group.',
      steps: variant => [
        PANCREAS_READY,
        ...look(variant),
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        { run: showOnlyAnalysisSection },
        { run: page => openAnalysisMode(page, 'genesPanel') },
        { waitForStable: '#analysis-panel-genesPanel select[name="groupBy"]' },
        { select: '#analysis-panel-genesPanel select[name="groupBy"]', value: 'cell_type' },
        { run: page => runAnalysisAndWaitForResult(page, 'genesPanel') },
        { run: page => expandAnalysisResult(page, 'genesPanel') },
        { run: page => page.mouse.move(4, 995) },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'element', selector: '.analysis-modal-content', padding: 32 },
    },
    {
      key: 'velocity',
      dataset: 'pancreas',
      describes: variant =>
        `The Pancreas 2D UMAP in the ${variant.theme} theme with the sidebar `
        + 'closed, coloured by cell type, and the RNA velocity overlay on: '
        + 'streaks of particles flowing along the embedding, brightest where '
        + 'the velocity is largest.',
      // The overlay is a particle animation that never stops changing, so the
      // frame is pinned rather than settled; `captures.json` records that.
      settle: { frames: 120 },
      steps: variant => [
        PANCREAS_READY,
        ...look(variant),
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        ...TO_2D,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#velocity-overlay-enabled' },
        { run: page => page.check('#velocity-overlay-enabled') },
        { waitFor: '#velocity-overlay-settings' },
        ...CLOSE_SIDEBAR,
        CLOSER_2D,
        NOTIFICATIONS_DRAINED,
      ],
    },
    {
      key: 'compare',
      dataset: 'pancreas',
      describes: variant =>
        `Two equal panels in the ${variant.theme} theme with the sidebar `
        + 'closed, both showing the Pancreas 3D embedding from one linked '
        + 'camera: on the left coloured by cell type, on the right by the '
        + 'S-phase score on a Viridis scale, each under a title chip naming '
        + 'its field.',
      steps: variant => [
        PANCREAS_READY,
        ...look(variant),
        { openSection: '#coloring-filtering-section' },
        ...colourBy('#categorical-field', 'cell_type', 'category'),
        { openSection: '#compare-views-section' },
        { openSection: '#compare-views-section' },
        { click: '#split-keep-view-btn' },
        {
          run: page => page.waitForFunction(() => (
            document.querySelectorAll('#view-title-layer .view-title-chip').length === 2
          )),
        },
        ...colourBy('#continuous-field', 'S_score', 'continuous'),
        { waitForText: '#view-title-layer .view-title-chip.active', text: 'S_score' },
        ...CLOSE_SIDEBAR,
        {
          run: async page => {
            await page.mouse.move(1080, 500);
            for (let notch = 0; notch < 3; notch += 1) {
              await page.mouse.wheel(0, 300);
              await page.waitForTimeout(80);
            }
            await page.mouse.move(4, 995);
          },
        },
        NOTIFICATIONS_DRAINED,
      ],
    },
  ];

  return [
    ...THEMES.flatMap(variant => [hero(variant, 1280), hero(variant, 2400)]),
    ...FEATURES.flatMap(spec => THEMES.map(variant => feature(variant, spec))),
  ];
}

// ---------------------------------------------------------------------------
// Opening each format, the opening itself, the providers' sign-in lines, and
// the Python and R viewers: b_data_loading 01–20 and the package guides.
//
// Every file a scenario opens is derived from public data by
// `fixtures/build-format-fixtures.py` (the Pancreas sample as an `.h5ad`, a
// Zarr ZIP, a Loom file, a Cell Ranger `outs/` folder, a MatrixMarket folder
// and Seurat objects; the web repository's CITE-seq corpus file), built once
// per run into a cache folder outside every repository. Every address is a
// loopback host of this run; the Google and Microsoft hosts are answered by
// routes, so no capture reaches a provider.
// ---------------------------------------------------------------------------

let formatFixturesReady = null;

/**
 * The path of one format fixture, building the set on first use.
 *
 * The cache folder is named by the digest of the two generator scripts, so an
 * edited recipe is rebuilt and an unchanged one is reused across runs.
 * `$CELLUCID_DOCS_PYTHON` (anndata, mudata, h5py, scipy) and
 * `$CELLUCID_DOCS_RSCRIPT` (SeuratObject) run them; both default to the
 * names on the PATH.
 *
 * @param {string} name - a file or folder the generator writes
 * @returns {Promise<string>}
 */
async function formatFixture(name) {
  formatFixturesReady ??= (async () => {
    const { createHash } = await import('node:crypto');
    const { access, readFile, writeFile } = await import('node:fs/promises');
    const { execFile } = await import('node:child_process');
    const { tmpdir } = await import('node:os');
    const { dirname, join } = await import('node:path');
    const { fileURLToPath } = await import('node:url');
    const here = join(dirname(fileURLToPath(import.meta.url)), 'fixtures');
    const hash = createHash('sha256');
    for (const script of ['build-format-fixtures.py', 'build-format-fixtures.R']) {
      hash.update(await readFile(join(here, script)));
    }
    const out = join(tmpdir(), `cellucid-docs-formats-${hash.digest('hex').slice(0, 12)}`);
    const complete = join(out, '.complete');
    try {
      await access(complete);
      return out;
    } catch {
      // Not built yet.
    }
    await new Promise((resolve, reject) => execFile(
      process.env.CELLUCID_DOCS_PYTHON ?? 'python3',
      [join(here, 'build-format-fixtures.py'), out],
      { maxBuffer: 1 << 24 },
      (error, _stdout, stderr) => (error
        ? reject(new Error(`build-format-fixtures.py failed: ${stderr || error.message}`))
        : resolve())
    ));
    await writeFile(complete, '');
    return out;
  })();
  const { join } = await import('node:path');
  return join(await formatFixturesReady, name);
}

/**
 * Serve a folder over loopback the way a bucket or a file host does: CORS,
 * HEAD, single ranges. With `stallAt`, no byte past that offset is sent: a
 * read across it sends what comes before and holds the connection open, so an
 * opening that streams the file stops part-way and stays there for the
 * photograph. With `gate`, nothing is
 * answered while `gate.hold` is true. Stopped when this process exits.
 *
 * @param {string} root
 * @param {{stallAt?: number|null, gate?: {hold: boolean}|null,
 *   contentTypeOf?: ((path: string) => string|undefined)|null,
 *   extra?: Map<string, {body: Buffer, type: string}>|null, ranges?: boolean}} [options]
 * @returns {Promise<string>} the origin
 */
async function serveFolder(root, { stallAt = null, gate = null, contentTypeOf = null, extra = null, ranges = true } = {}) {
  const http = await import('node:http');
  const { createReadStream } = await import('node:fs');
  const { stat } = await import('node:fs/promises');
  const { resolve, sep } = await import('node:path');
  const base = resolve(root);
  const server = http.createServer(async (request, response) => {
    response.setHeader('Access-Control-Allow-Origin', '*');
    response.setHeader('Access-Control-Allow-Headers', 'Range');
    response.setHeader(
      'Access-Control-Expose-Headers',
      'Content-Length, Content-Range, Accept-Ranges, ETag, Last-Modified'
    );
    if (request.method === 'OPTIONS') {
      response.writeHead(204);
      response.end();
      return;
    }
    // A held request is never answered: whatever asked for it waits.
    if (gate?.hold === true) return;
    const own = extra?.get(decodeURIComponent(new URL(request.url, 'http://x').pathname).slice(1));
    if (own !== undefined) {
      response.writeHead(200, { 'Content-Type': own.type, 'Content-Length': String(own.body.length) });
      response.end(request.method === 'HEAD' ? undefined : own.body);
      return;
    }
    const file = resolve(base, `.${decodeURIComponent(new URL(request.url, 'http://x').pathname)}`);
    let info = null;
    if (file.startsWith(base + sep)) {
      info = await stat(file).catch(() => null);
    }
    if (info === null || !info.isFile()) {
      response.writeHead(404, { 'Content-Type': 'text/plain' });
      response.end('Not found');
      return;
    }
    let start = 0;
    let end = info.size - 1;
    let status = 200;
    const range = ranges ? /^bytes=(\d*)-(\d*)$/.exec(request.headers.range ?? '') : null;
    if (range !== null) {
      if (range[1] === '') {
        start = Math.max(0, info.size - Number(range[2]));
      } else {
        start = Number(range[1]);
        if (range[2] !== '') end = Math.min(end, Number(range[2]));
      }
      status = 206;
      response.setHeader('Content-Range', `bytes ${start}-${end}/${info.size}`);
    }
    response.writeHead(status, {
      'Accept-Ranges': ranges ? 'bytes' : 'none',
      'Content-Length': String(end - start + 1),
      'Content-Type': contentTypeOf?.(file.slice(base.length + 1).split(sep).join('/'))
        ?? 'application/octet-stream',
      ETag: `"${info.size}-${Math.round(info.mtimeMs)}"`,
      'Last-Modified': info.mtime.toUTCString(),
    });
    if (request.method === 'HEAD' || info.size === 0) {
      response.end();
      return;
    }
    // Bytes past `stallAt` are never sent, whether asked for whole or by
    // range: the reader is left waiting at that point of the file.
    if (stallAt !== null && start >= stallAt) return;
    const stalls = stallAt !== null && end >= stallAt;
    const stream = createReadStream(file, {
      start,
      end: stalls ? stallAt - 1 : end,
    });
    stream.pipe(response, { end: !stalls });
  });
  // A client that reads hundreds of files in a row over one connection (the
  // web build a server downloads) must not have it closed between two of
  // them: Node's default 5 s keep-alive races the next request.
  server.keepAliveTimeout = 120_000;
  server.headersTimeout = 125_000;
  const key = `serveFolder ${root} ${JSON.stringify({ stallAt, ranges, gate: gate !== null, extra: extra === null ? null : [...extra.keys()] })}`;
  const port = await freePort(key);
  await new Promise(done => server.listen(port, '127.0.0.1', done));
  process.once('exit', () => server.close());
  return `http://127.0.0.1:${server.address().port}`;
}

/** A free loopback port. */
/**
 * A free loopback port that is the same on every run for the same `key`.
 *
 * Several captures print a port: a server's address in the Link or server
 * bar, a link's host in the Dataset box or the opening panel. A port the
 * system hands out (`listen(0)`) differs on every run, so those images could
 * never be captured twice with the same bytes. The port is instead derived
 * from the key (`20000`–`39999`), and the next one is taken only when that
 * one is busy on this machine.
 *
 * @param {string} key - What the port is for, unique to its scenario
 * @returns {Promise<number>}
 */
async function freePort(key) {
  if (typeof key !== 'string' || key === '') {
    throw new TypeError('freePort needs a key, so the port is the same on every run.');
  }
  const { createServer } = await import('node:net');
  const { createHash } = await import('node:crypto');
  const start = createHash('sha256').update(key).digest().readUInt32BE(0) % 20000;
  for (let attempt = 0; attempt < 200; attempt += 1) {
    const port = 20000 + ((start + attempt) % 20000);
    const free = await new Promise(done => {
      const probe = createServer();
      probe.once('error', () => done(false));
      probe.listen(port, '127.0.0.1', () => probe.close(() => done(true)));
    });
    if (free) return port;
  }
  throw new Error(`No free loopback port near the one for "${key}".`);
}

/**
 * Open a file or folder the way a reader does: **Open…**, then **File…** or
 * **Folder…**, then the browser's chooser.
 *
 * @param {import('@playwright/test').Page} page
 * @param {'file'|'folder'} kind
 * @param {string} target
 * @returns {Promise<string>} the dataset name shown before
 */
async function openFromDisk(page, kind, target) {
  // The name and the origin line together: a file of the sample opened
  // from disk keeps the sample's name, but not its origin.
  const before = await page.evaluate(() => (
    `${document.getElementById('dataset-name')?.textContent}|`
      + `${document.getElementById('dataset-origin')?.textContent}`
  ));
  await page.click('#user-data-open-btn');
  await page.waitForSelector('#user-data-open-menu', { state: 'visible' });
  const [chooser] = await Promise.all([
    page.waitForEvent('filechooser'),
    page.click(kind === 'file' ? '#user-data-file-item' : '#user-data-folder-item'),
  ]);
  await chooser.setFiles(target);
  return before;
}

/** Wait until another dataset than `before` is open and drawn. */
function datasetReplaced(before) {
  return {
    run: page => page.waitForFunction(
      previous => (
        `${document.getElementById('dataset-name')?.textContent}|`
          + `${document.getElementById('dataset-origin')?.textContent}` !== previous
          && document.getElementById('dataset-info')?.dataset.state === 'ready'
          && document.documentElement.dataset.startupPhase === 'ready'
      ),
      before,
      { timeout: 180_000 }
    ),
  };
}

/** Open a fixture from disk and wait for it to replace the Pancreas sample. */
function openFixture(kind, name) {
  return {
    run: async page => {
      const before = await openFromDisk(page, kind, await formatFixture(name));
      await datasetReplaced(before).run(page);
    },
  };
}

/**
 * Google's hosts, answered on this machine: the Drive API says 404 for the
 * file (what it says to anyone the file is not shared with), and Google's two
 * scripts are stubs. `tokenError` is the OAuth error the stubbed sign-in
 * answers a click with; without it the click is never answered.
 */
async function stubGoogle(context, { tokenError = null } = {}) {
  await context.route('https://www.googleapis.com/**', route => route.fulfill({
    status: 404,
    contentType: 'application/json',
    headers: { 'Access-Control-Allow-Origin': '*' },
    body: JSON.stringify({ error: { code: 404, message: 'File not found: .', errors: [{ domain: 'global', reason: 'notFound', message: 'File not found: .' }] } }),
  }));
  const answer = tokenError === null
    ? ''
    : `const client = this; setTimeout(() => client.callback({ error: ${JSON.stringify(tokenError)} }), 0);`;
  await context.route('https://accounts.google.com/**', route => route.fulfill({
    contentType: 'text/javascript',
    body: 'window.google = window.google || {}; google.accounts = { oauth2: { '
      + `initTokenClient: () => ({ requestAccessToken() { ${answer} } }) } };`,
  }));
  await context.route('https://apis.google.com/**', route => route.fulfill({
    contentType: 'text/javascript',
    body: 'window.gapi = { load(name, options) { window.google = window.google || {}; google.picker = {}; (options.callback || options)(); } };',
  }));
}

function buildFormatAndViewerScenarios() {
  const DATASET_BOX_WITH_DETAILS = [
    NOTIFICATIONS_DRAINED,
    { openSection: '#session-section' },
    { click: '#dataset-details .disclosure-header' },
    { waitFor: '#dataset-details-list' },
    { scrollIntoView: '#dataset-facts' },
  ];
  const SETTLE = Object.freeze({ attempts: 2, gapMs: 0 });
  const AT_ADDRESS_BAR = [
    { openSection: '#session-section' },
    { scrollIntoView: '#data-address-block' },
  ];
  // A private Google Drive file pasted into the bar, up to its sign-in line.
  const GOOGLE_PRIVATE_FILE_STEPS = [
    PANCREAS_READY,
    NOTIFICATIONS_DRAINED,
    ...AT_ADDRESS_BAR,
    { fill: '#data-address', value: 'https://drive.google.com/file/d/1AbCdEfGhIjKlMnOpQrStUvWxYz0123456/view' },
    { click: '#data-address-open-btn' },
    { waitFor: '#data-sign-in-step' },
    { waitForText: '#data-sign-in-step', text: 'Private file on Google Drive' },
    { run: page => page.waitForFunction(() => !document.getElementById('data-sign-in-btn')?.disabled) },
  ];
  const datasetBox = ({ id, name, describes, kind, file }) => ({
    id,
    topic: 'data_loading',
    name,
    describes,
    dataset: 'pancreas',
    steps: [PANCREAS_READY, NOTIFICATIONS_DRAINED, openFixture(kind, file), ...DATASET_BOX_WITH_DETAILS],
    shot: { mode: 'element', selector: '#dataset-facts', padding: 6 },
    settle: SETTLE,
  });

  return [
    // ------------------------------------------------ an .h5ad, whole window
    {
      id: 'loading-open-h5ad-window',
      topic: 'data_loading',
      name: 'h5ad-current-loaded',
      describes:
        'The whole Cellucid window after the Pancreas sample was opened as an '
        + '.h5ad file with Open… › File…: the Session panel\'s Dataset box names '
        + 'the file, its cells and genes and "h5ad file · this computer", and '
        + 'the canvas shows its UMAP.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        openFixture('file', 'pancreas.h5ad'),
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
      ],
      shot: { mode: 'window' },
    },
    datasetBox({
      id: 'loading-open-zarr-zip-dataset-box',
      name: 'zarr-zip-loaded',
      describes:
        'The Dataset box after a Zarr ZIP archive was opened from this '
        + 'computer: the dataset\'s name, its cells and genes, "Zarr v2 ZIP · '
        + 'this computer", its description, and Details open with its source and '
        + 'citation, the file, its format and size, how it is read, and what it '
        + 'holds.',
      kind: 'file',
      file: 'pancreas.zarr.zip',
    }),
    datasetBox({
      id: 'loading-open-rds-dataset-box',
      name: 'rds-dataset-box',
      describes:
        'The Dataset box after a Seurat object saved with saveRDS() was opened: '
        + 'its cells and genes, "R object (.rds) · this computer", the line '
        + '"Read whole once (it cannot be read in parts), kept in this browser\'s '
        + 'storage for next time", and Details open with Read (whole, once, from '
        + 'start to end, nothing uploaded), Kept (the origin-private file system '
        + 'until it is cleared), Read in parts (a Cellucid-ready AnnData) and what '
        + 'it holds.',
      kind: 'file',
      file: 'pancreas.rds',
    }),
    datasetBox({
      id: 'loading-open-loom-dataset-box',
      name: 'loom-dataset-box',
      describes:
        'The Dataset box after a Loom file was opened: the title stored in the '
        + 'file, its cells and genes, "Loom (.loom) · this computer", and '
        + 'Details open with the file, its size, how it is read, and what it '
        + 'holds: four obs fields and a UMAP.',
      kind: 'file',
      file: 'pancreas.loom',
    }),
    datasetBox({
      id: 'loading-open-matrix-market-dataset-box',
      name: 'matrix-market-dataset-box',
      describes:
        'The Dataset box after a MatrixMarket folder with a cell table was '
        + 'opened with Open… › Folder…: the folder\'s name, its cells and genes, '
        + '"MatrixMarket folder · this computer", the line saying where the '
        + 'converted data is kept, and Details open.',
      kind: 'folder',
      file: 'pancreas_matrix',
    }),
    {
      id: 'loading-open-tenx-outs-window',
      topic: 'data_loading',
      name: 'tenx-outs-window',
      describes:
        'The whole Cellucid window after a Cell Ranger outs/ folder was opened '
        + 'with Open… › Folder…: the Dataset box reads "10x Cell Ranger outs '
        + 'folder · this computer", and the canvas shows Cell Ranger\'s UMAP '
        + 'coloured by its clustering, gene_expression_clusters.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        openFixture('folder', 'outs'),
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'gene_expression_clusters' },
        { waitForText: '#stats', text: 'gene_expression_clusters' },
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
      ],
      shot: { mode: 'window' },
    },

    // ----------------------------------------- an .rda of two Seurat objects
    {
      id: 'loading-open-rda-choose-object',
      topic: 'data_loading',
      name: 'rda-choose-object',
      describes:
        'The opening panel over the viewer reading "Opening '
        + 'pancreas-objects.rda" and "pancreas-objects.rda holds 2 single-cell '
        + 'objects. Open which one?", with one button per Seurat object, '
        + '"pancreas (Seurat, 3,696 cells)" and "pancreas_ductal (Seurat, 916 '
        + 'cells)", and Cancel.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        { run: holdNotificationTimers },
        { run: async page => openFromDisk(page, 'file', await formatFixture('pancreas-objects.rda')) },
        { waitForText: '.opening-overlay .opening-actions', text: 'pancreas_ductal' },
        { run: page => page.mouse.move(4, 995) },
      ],
      shot: { mode: 'element', selector: '.opening-overlay .opening-panel', padding: 16 },
      settle: SETTLE,
    },

    // ------------------------------------------------- the Modality select
    {
      id: 'loading-h5mu-modality-select',
      topic: 'data_loading',
      name: 'data-mapping-modality',
      describes:
        'The Data mapping group of a CITE-seq .h5mu with Genes open: a Modality '
        + 'select reading "rna · 28 features" above Expression, Gene names, '
        + 'Values and Highly variable genes, and Revert and Apply below.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        openFixture('file', 'citeseq-subset.h5mu'),
        // The file has no X_umap, so it opens on a best guess and says so in a
        // notification that stays until dismissed; it is outside this crop.
        { openSection: '#session-section' },
        { run: async page => {
          const genes = page.locator('[data-mapping-section="genes"] .disclosure-header');
          if ((await genes.getAttribute('aria-expanded')) !== 'true') await genes.click();
          await page.waitForSelector('#data-mapping-modality', { state: 'visible' });
        } },
        { scrollIntoView: '#data-mapping-controls' },
      ],
      shot: { mode: 'element', selector: '#data-mapping-controls', padding: 6 },
      cursor: { on: '#data-mapping-modality', at: 'right', dx: -12, state: 'hover' },
      settle: SETTLE,
    },

    // ------------------------------------------------ the opening, under way
    // A Seurat `.rds` linked at startup is read once from the start. The host
    // serves whole files only (as some web servers do), sends 60% and holds
    // the rest, so the panel stops part-way with its byte count, speed and
    // time left on screen. (From a host that serves ranges the count stays at
    // 0 B until the read completes — reported as a defect, not pictured.)
    {
      id: 'loading-opening-progress',
      topic: 'data_loading',
      name: 'opening-progress',
      describes:
        'The opening panel centred over the empty viewer while a linked Seurat '
        + '.rds is downloaded: "Opening" and the file\'s address, the step '
        + '"Downloading the whole file (9.4 MB): … does not serve parts of '
        + 'files", a progress bar part-way, how much has arrived with the speed '
        + 'and time left, the note saying the host serves the file only whole, '
        + 'how much is downloaded and that it is kept, and a Cancel button.',
      awaitStartup: false,
      query: async () => {
        const file = await formatFixture('pancreas.rds');
        const { stat } = await import('node:fs/promises');
        const { dirname, basename } = await import('node:path');
        const { size } = await stat(file);
        const origin = await serveFolder(dirname(file), { stallAt: Math.round(size * 0.6), ranges: false });
        return { data: `${origin}/${basename(file)}` };
      },
      steps: [
        { waitFor: '.opening-overlay .opening-panel' },
        { run: page => page.waitForFunction(() => (
          /of/.test(document.querySelector('.opening-overlay .opening-detail')?.textContent ?? '')
        ), undefined, { timeout: 120_000 }) },
        { wait: 2500 },
      ],
      shot: { mode: 'element', selector: '.opening-overlay .opening-panel', padding: 24 },
      settle: { attempts: 3, gapMs: 0 },
    },

    // ------------------------------------------------- the opening, failed
    {
      id: 'loading-opening-failed',
      topic: 'data_loading',
      name: 'opening-failed',
      describes:
        'The panel over the empty viewer after a linked file could not be '
        + 'opened: "The data could not be opened", one sentence saying why (the '
        + 'file was not found at its host, HTTP 404) and what to do, a folded '
        + 'Details, and the buttons Open without this link, focused, and Try '
        + 'again.',
      awaitStartup: false,
      query: async () => {
        const origin = await serveFolder(await formatFixture('.'));
        return { data: `${origin}/cells.h5ad` };
      },
      steps: [
        { waitForText: '.opening-overlay .opening-title', text: 'could not be opened' },
        { waitForText: '.opening-overlay .opening-actions', text: 'Open without this link' },
      ],
      shot: { mode: 'element', selector: '.opening-overlay .opening-panel', padding: 24 },
      settle: SETTLE,
    },

    // -------------------------------------------- Google Drive, private file
    // Served as `http://127.0.0.1:8765`, an origin Cellucid registered with
    // Google, because sign-in exists only there. The Drive API answers 404 for
    // a file the visitor may not read, as it does for any private file, and
    // Google's two scripts are answered by stubs: nothing reaches Google.
    {
      id: 'loading-gdrive-private-sign-in',
      topic: 'data_loading',
      name: 'gdrive-private-sign-in',
      describes:
        'The Link or server bar holding a Google Drive file link, with one line '
        + 'under it reading "Private file on Google Drive.", a Sign in to Google '
        + 'button and an (i).',
      dataset: 'pancreas',
      servedAs: 'http://127.0.0.1:8765',
      setup: ({ context }) => stubGoogle(context),
      steps: GOOGLE_PRIVATE_FILE_STEPS,
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      settle: SETTLE,
    },

    // ------------------- Google Drive, an account its administrator blocks
    // The same private file; Google's token answer is the one a Google
    // Workspace account gets when its administrator has not allowed
    // Cellucid (`admin_policy_enforced`), from the stubbed sign-in script.
    {
      id: 'loading-gdrive-admin-blocked',
      topic: 'data_loading',
      name: 'gdrive-admin-blocked',
      describes:
        'The Link or server bar after Sign in to Google, for an account whose '
        + 'Google Workspace administrator does not allow Cellucid: the line '
        + 'under it says so and what to do, with a link "How an administrator '
        + 'allows Cellucid" and the button Use another account.',
      dataset: 'pancreas',
      servedAs: 'http://127.0.0.1:8765',
      setup: ({ context }) => stubGoogle(context, { tokenError: 'admin_policy_enforced' }),
      steps: [
        ...GOOGLE_PRIVATE_FILE_STEPS,
        { click: '#data-sign-in-btn' },
        { waitForText: '#data-sign-in-btn', text: 'Use another account' },
        { waitFor: '#data-sign-in-help' },
      ],
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      settle: SETTLE,
    },

    // ------------------------------------------------------ OneDrive file
    // OneDrive serves no file anonymously, so every link asks for a Microsoft
    // sign-in first. Microsoft's hosts are refused by routes: the picture is
    // the line Cellucid shows before anything is asked of Microsoft.
    {
      id: 'loading-onedrive-sign-in',
      topic: 'data_loading',
      name: 'onedrive-sign-in',
      describes:
        'The Link or server bar holding a OneDrive share link, with one line '
        + 'under it reading "OneDrive file: sign in with Microsoft to open it.", '
        + 'a Sign in with Microsoft button and an (i).',
      dataset: 'pancreas',
      servedAs: 'http://127.0.0.1:8765',
      setup: async ({ context }) => {
        for (const host of ['https://login.microsoftonline.com/**', 'https://graph.microsoft.com/**']) {
          await context.route(host, route => route.abort());
        }
      },
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        ...AT_ADDRESS_BAR,
        { fill: '#data-address', value: 'https://1drv.ms/u/c/0a1b2c3d4e5f6a7b/EaBcDeFgHiJkLmNoPqRsTuVw' },
        { click: '#data-address-open-btn' },
        { waitFor: '#data-sign-in-step' },
        { waitForText: '#data-sign-in-step', text: 'sign in with Microsoft' },
        { run: page => page.waitForFunction(() => !document.getElementById('data-sign-in-btn')?.disabled) },
      ],
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      settle: SETTLE,
    },

    // ---------------------------------------------- a GitHub path, typed
    {
      id: 'loading-github-path-typed',
      topic: 'data_loading',
      name: 'connect-github-catalog',
      describes:
        'The Link or server bar with the GitHub path '
        + '"theislab/cellucid-demo-custom-datasets/catalog" typed into it and the '
        + 'pointer on its Open button.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        { fill: '#data-address', value: 'theislab/cellucid-demo-custom-datasets/catalog' },
      ],
      shot: { mode: 'element', selector: '#data-address-block', padding: 8 },
      cursor: { on: '#data-address-open-btn', state: 'hover' },
    },

    // --------------------------- a GitHub repository with no datasets.json
    // GitHub's hosts are answered by routes with the 404 they give for a
    // file a repository does not hold.
    {
      id: 'loading-link-refusal-github-not-found',
      topic: 'data_loading',
      name: 'fail-github-catalog-not-found',
      describes:
        'A Cellucid notification saying the linked data could not be opened '
        + 'because the GitHub folder holds no datasets.json, after a repository '
        + 'path rather than its catalog folder was opened.',
      dataset: 'pancreas',
      setup: async ({ context }) => {
        for (const host of ['https://raw.githubusercontent.com/**', 'https://api.github.com/**', 'https://cdn.jsdelivr.net/**']) {
          await context.route(host, route => route.fulfill({
            status: 404,
            contentType: 'text/plain',
            headers: { 'Access-Control-Allow-Origin': '*' },
            body: '404: Not Found',
          }));
        }
      },
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        NOTIFICATIONS_DRAINED,
        { fill: '#data-address', value: 'theislab/cellucid-annotation' },
        { click: '#data-address-open-btn' },
        { run: page => page.waitForFunction(() => (
          document.querySelector('#notification-center .notification[data-type="error"]') !== null
        ), undefined, { timeout: 120_000 }) },
        { waitForStable: '#notification-center', forMs: 1000 },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: SETTLE,
    },

    // ---------------------------------- a server address with nothing there
    {
      id: 'loading-address-refusal-server-unreachable',
      topic: 'data_loading',
      name: 'fail-remote-server-unreachable',
      describes:
        'A Cellucid notification saying the Cellucid server could not be '
        + 'reached at a loopback address, and what to check.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        NOTIFICATIONS_DRAINED,
        { run: async page => page.fill('#data-address', `127.0.0.1:${await freePort('nothing listening (server unreachable)')}`) },
        { click: '#data-address-open-btn' },
        { run: page => page.waitForFunction(() => (
          document.querySelector('#notification-center .notification[data-type="error"]') !== null
        ), undefined, { timeout: 120_000 }) },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: SETTLE,
    },
    {
      id: 'loading-keep-dataset-after-failure',
      topic: 'data_loading',
      name: 'keep-dataset-after-failed-connect',
      describes:
        'The whole Cellucid window after a server address with nothing behind it '
        + 'was opened: the Pancreas sample is still drawn and still named in '
        + 'the Dataset box, and one error notification sits in the lower right.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...AT_ADDRESS_BAR,
        NOTIFICATIONS_DRAINED,
        { run: async page => page.fill('#data-address', `127.0.0.1:${await freePort('nothing listening (keep the dataset)')}`) },
        { click: '#data-address-open-btn' },
        { run: page => page.waitForFunction(() => (
          document.querySelector('#notification-center .notification[data-type="error"]') !== null
        ), undefined, { timeout: 120_000 }) },
        { scrollIntoView: '#dataset-facts' },
      ],
      shot: { mode: 'window' },
      settle: SETTLE,
    },

    // --------------------------------------------------- cellucid serve
    // The product's own server on the Pancreas sample's store, serving this
    // run's application as its web build (`--web-source-url`), opened at the
    // Viewer URL it prints.
    {
      id: 'server-cellucid-serve-window',
      topic: 'server',
      name: 'open-served-dataset',
      describes:
        'Cellucid opened at the Viewer URL of cellucid serve, coloured by '
        + 'clusters: the Dataset box names the Pancreas sample with its cells '
        + 'and genes, "Zarr v2 folder" and the server\'s loopback address, and '
        + 'how much was read; the Link or server bar holds that address with '
        + 'one line under it reading Server, the address and Disconnect.',
      entry: async urls => startPythonServe(urls),
      steps: [
        { run: page => page.waitForFunction(
          () => document.documentElement.dataset.startupPhase === 'ready',
          undefined,
          { timeout: 180_000 }
        ) },
        { waitFor: '#remote-disconnect-container' },
        { waitForStable: '#categorical-field' },
        { select: '#categorical-field', label: 'clusters' },
        { waitForText: '#stats', text: 'Field: clusters (category)' },
        { run: page => page.waitForFunction(() => {
          const line = document.getElementById('dataset-read');
          return line !== null && !line.hidden && !/Reading/.test(line.textContent);
        }, undefined, { timeout: 120_000 }) },
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
      ],
      shot: { mode: 'window' },
    },

    // --------------------------------- show_anndata() from a Python process
    {
      id: 'python-viewer-driven-from-python',
      topic: 'jupyter',
      name: 'show-anndata-driven-from-python',
      describes:
        'The viewer cellucid.show_anndata() serves, opened at its viewer URL, '
        + 'after Python called viewer.set_color_by("clusters") and '
        + 'viewer.highlight_cells() on the Ductal cells: the Pancreas sample '
        + 'coloured by clusters with the Ductal cells highlighted.',
      setup: pinHighlightGlow,
      entry: async urls => startPythonViewer(urls),
      steps: [
        { run: page => page.waitForFunction(
          () => document.documentElement.dataset.startupPhase === 'ready',
          undefined,
          { timeout: 180_000 }
        ) },
        { run: () => pythonViewerDriven },
        { waitForText: '#stats', text: 'Field: clusters' },
        NOTIFICATIONS_DRAINED,
        { openSection: '#highlighted-cells-section' },
      ],
      shot: { mode: 'window' },
    },

    // ------------------------------------------- cellucid_view() from R
    {
      id: 'r-viewer-driven-from-r',
      topic: 'r_package',
      name: 'cellucid-view-driven-from-r',
      describes:
        'Cellucid opened by cellucid_view() on a Seurat object in an R session, '
        + 'after R called cellucid_color_by(viewer, "clusters") and '
        + 'cellucid_highlight() on the Ductal cells: the points are coloured by '
        + 'clusters with the Ductal cells highlighted, and the Dataset box names '
        + 'the object.',
      setup: pinHighlightGlow,
      entry: async urls => startRViewer(urls),
      steps: [
        { run: page => page.waitForFunction(
          () => document.documentElement.dataset.startupPhase === 'ready',
          undefined,
          { timeout: 180_000 }
        ) },
        { run: () => rViewerDriven },
        { waitForText: '#stats', text: 'Field: clusters' },
        NOTIFICATIONS_DRAINED,
        { openSection: '#session-section' },
      ],
      shot: { mode: 'window' },
    },
  ];
}

/**
 * `cellucid serve` on the Pancreas sample's store, with this run's
 * application as its web build. Returns the Viewer URL once the health
 * endpoint answers; stopped when this process exits.
 *
 * @param {{appOrigin: string}} urls
 * @returns {Promise<string>}
 */
async function startPythonServe(urls) {
  const { spawn } = await import('node:child_process');
  const { mkdtemp } = await import('node:fs/promises');
  const { tmpdir } = await import('node:os');
  const { dirname, join, resolve } = await import('node:path');
  const { fileURLToPath } = await import('node:url');
  const here = dirname(fileURLToPath(import.meta.url));
  const datasets = process.env.CELLUCID_DATASETS_REPO
    ?? resolve(here, '..', '..', '..', '..', 'cellucid-datasets');
  const port = await freePort('cellucid serve (Viewer URL)');
  const cache = await mkdtemp(join(tmpdir(), 'cellucid-docs-web-'));
  const output = [];
  const child = spawn(process.env.CELLUCID_CLI ?? 'cellucid', [
    'serve', join(datasets, 'anndata', 'pancreas.zarr'),
    '--host', '127.0.0.1', '--port', String(port), '--no-browser',
    '--web-source-url', await webBuildOrigin(), '--web-cache-dir', cache,
  ], {
    stdio: ['ignore', 'pipe', 'pipe'],
    // Into a pipe Python buffers its output by block, so the banner with the
    // Viewer URL would arrive only when the server exits.
    env: { ...process.env, PYTHONUNBUFFERED: '1' },
  });
  child.stdout.on('data', chunk => output.push(String(chunk)));
  child.stderr.on('data', chunk => output.push(String(chunk)));
  process.once('exit', () => child.kill());
  const deadline = Date.now() + 300_000;
  for (;;) {
    if (child.exitCode !== null) {
      throw new Error(`cellucid serve exited with status ${child.exitCode}:\n${output.join('')}`);
    }
    const viewer = /Viewer URL:\s+(\S+)/.exec(output.join(''));
    if (viewer !== null) return viewer[1];
    if (Date.now() > deadline) {
      child.kill();
      throw new Error(`cellucid serve printed no Viewer URL within 300 s:\n${output.join('')}`);
    }
    await new Promise(next => setTimeout(next, 250));
  }
}

/** Settles once the R session has coloured and highlighted the viewer. */
let rViewerDriven = null;

/**
 * An R session (`$CELLUCID_DOCS_RSCRIPT`, else `Rscript`; the cellucid package
 * from `$CELLUCID_DOCS_R_LIBS` when set) that reads the
 * Pancreas Seurat object, shows it with `cellucid_view(pane = "none")`, and,
 * once the browser connects, colours it and highlights the Ductal cells from
 * R. Returns the viewer URL R printed; `rViewerDriven` settles when R has sent
 * both commands. Stopped when this process exits.
 *
 * @param {{appOrigin: string}} urls
 * @returns {Promise<string>}
 */
async function startRViewer(urls) {
  const { spawn } = await import('node:child_process');
  const { mkdtemp } = await import('node:fs/promises');
  const { tmpdir } = await import('node:os');
  const { join } = await import('node:path');
  const rds = await formatFixture('pancreas.rds');
  const cache = await mkdtemp(join(tmpdir(), 'cellucid-docs-web-'));
  const port = await freePort('cellucid_view() (R viewer)');
  // `$CELLUCID_DOCS_R_LIBS` puts a library first on the search path, so the
  // viewer is driven by the cellucid package of the tree being documented
  // (built from it with `R CMD build` and `R CMD INSTALL -l`) rather than an
  // older one installed in the environment.
  const libraries = process.env.CELLUCID_DOCS_R_LIBS ?? '';
  const script = [
    `if (nzchar(${JSON.stringify(libraries)})) .libPaths(c(${JSON.stringify(libraries)}, .libPaths()))`,
    'suppressPackageStartupMessages({ library(cellucid); library(SeuratObject) })',
    `pancreas <- readRDS(${JSON.stringify(rds)})`,
    `viewer <- cellucid_view(pancreas, pane = "none", port = ${port}L, quiet = TRUE,`,
    `  web_source_url = ${JSON.stringify(await webBuildOrigin())}, web_cache_dir = ${JSON.stringify(cache)})`,
    'cat("VIEWER_URL", viewer$viewer_url, "\\n"); flush(stdout())',
    'cellucid_wait_ready(viewer, timeout = 300)',
    'cellucid_color_by(viewer, "clusters")',
    'cellucid_highlight(viewer, colnames(pancreas)[pancreas$clusters == "Ductal"])',
    'cat("DRIVEN\\n"); flush(stdout())',
    'repeat httpuv::service(100)',
  ].join('\n');
  const output = [];
  const child = spawn(process.env.CELLUCID_DOCS_RSCRIPT ?? 'Rscript', ['-e', script], {
    stdio: ['ignore', 'pipe', 'pipe'],
  });
  child.stdout.on('data', chunk => output.push(String(chunk)));
  child.stderr.on('data', chunk => output.push(String(chunk)));
  process.once('exit', () => child.kill());
  rViewerDriven = (async () => {
    const deadline = Date.now() + 300_000;
    while (!output.join('').includes('DRIVEN')) {
      if (child.exitCode !== null || Date.now() > deadline) {
        throw new Error(`The R session did not drive the viewer:\n${output.join('')}`);
      }
      await new Promise(next => setTimeout(next, 250));
    }
  })();
  // Awaited by a step; observed here too, so a session that dies before that
  // step runs fails the scenario instead of the whole process.
  rViewerDriven.catch(() => {});
  const deadline = Date.now() + 300_000;
  for (;;) {
    if (child.exitCode !== null) {
      throw new Error(`The R session exited with status ${child.exitCode}:\n${output.join('')}`);
    }
    const viewer = /VIEWER_URL (\S+)/.exec(output.join(''));
    if (viewer !== null) return viewer[1];
    if (Date.now() > deadline) {
      child.kill();
      throw new Error(`The R session printed no viewer URL within 300 s:\n${output.join('')}`);
    }
    await new Promise(next => setTimeout(next, 250));
  }
}

/** While `hold` is true the gated sample catalog answers nothing. */
const SAMPLE_CATALOG_GATE = { hold: false };
let gatedCatalogOrigin = null;

/**
 * The sample catalog (`cellucid-datasets/anndata/`) on a host of this run
 * that stops answering while `SAMPLE_CATALOG_GATE.hold` is set: an analysis
 * started then holds at its first read, so its running state can be
 * photographed. Reset to answering for every scenario that asks for it.
 *
 * @returns {Promise<string>}
 */
async function gatedSampleCatalog() {
  SAMPLE_CATALOG_GATE.hold = false;
  if (gatedCatalogOrigin === null) {
    const { dirname, join, resolve } = await import('node:path');
    const { fileURLToPath } = await import('node:url');
    const here = dirname(fileURLToPath(import.meta.url));
    const datasets = process.env.CELLUCID_DATASETS_REPO
      ?? resolve(here, '..', '..', '..', '..', 'cellucid-datasets');
    gatedCatalogOrigin = await serveFolder(join(datasets, 'anndata'), { gate: SAMPLE_CATALOG_GATE });
  }
  return gatedCatalogOrigin;
}

// ---------------------------------------------------------------------------
// Orientation, navigation and connectivity, the camera path, the vector field
// overlay, figure export and the benchmark: the figures of a_orientation,
// c_core_interactions, i_vector_field_velocity, k_figure_export and
// n_benchmarking_performance that earlier sweeps captured by hand.
// ---------------------------------------------------------------------------

const ALL_SECTIONS = Object.freeze([
  '#session-section',
  '#visualization-section',
  '#navigation-section',
  '#compare-views-section',
  '#coloring-filtering-section',
  '#community-annotation-section',
  '#highlighted-cells-section',
  '#page-analysis-section',
  '#cinematic-camera-section',
  '#figure-export-section',
  '#shortcuts-section',
  '#benchmark-section',
]);

/** Leave only these sidebar sections open, so a crop starts at the top. */
function onlySections(...open) {
  return {
    run: page => page.evaluate(([all, keep]) => {
      for (const selector of all) {
        document.querySelector(selector).open = keep.includes(selector);
      }
    }, [ALL_SECTIONS, open]),
  };
}

function buildViewAndExportScenarios() {
  const settle = { attempts: 20, gapMs: 200 };
  const colourBy = label => [
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label },
    { waitForText: '#stats', text: `Field: ${label} (category)` },
  ];
  /** Show one representation in the selected view and wait for it to land. */
  const show = (label, mode) => [
    { openSection: '#compare-views-section' },
    { waitFor: '#representation-controls:not([hidden])' },
    { waitForStable: '#view-representation-select' },
    { run: async page => {
      const [value] = await page.selectOption('#view-representation-select', { label });
      await page.waitForFunction(
        ([expected, navigation]) => (
          document.querySelector('#view-representation-select')?.value === expected
            && document.querySelector('#navigation-mode')?.value === navigation
        ),
        [value, mode],
      );
    } },
    { wait: 1500 },
  ];
  const SUO_READY = {
    run: page => page.waitForFunction(
      () => document.documentElement.dataset.startupPhase === 'ready'
        && document.getElementById('dataset-info')?.dataset.state === 'ready',
      undefined,
      { timeout: 300_000 }
    ),
  };
  const KANEMARU_READY = SUO_READY;
  /** Orbit the camera by a real drag on the canvas, away from the sidebar. */
  const orbitBy = dx => ({
    run: async page => {
      await page.mouse.move(900, 500);
      await page.mouse.down();
      for (let step = 1; step <= 10; step += 1) {
        await page.mouse.move(900 + (dx * step) / 10, 500);
      }
      await page.mouse.up();
      await page.waitForTimeout(800);
    },
  });
  const TWO_KEYFRAMES = [
    onlySections('#cinematic-camera-section'),
    { click: '#cinematic-save-btn' },
    orbitBy(120),
    { click: '#cinematic-save-btn' },
    { waitFor: '.cinematic-transport-bar' },
  ];
  const VELOCITY_ON = [
    { openSection: '#visualization-section' },
    { waitFor: '#velocity-overlay-controls:not([hidden])' },
    { run: page => page.check('#velocity-overlay-enabled') },
    { waitFor: '#velocity-overlay-settings:not([hidden])' },
  ];
  const EXPORT_PANEL = [
    PANCREAS_READY,
    ...colourBy('clusters'),
    onlySections('#figure-export-section'),
    { waitFor: '#figure-export-btn' },
  ];
  const openExportGroup = title => ({
    run: async page => {
      const header = page.locator('#figure-export-section .figure-export-accordion .disclosure-header', { hasText: title });
      if ((await header.getAttribute('aria-expanded')) !== 'true') await header.click();
      await page.waitForTimeout(400);
    },
  });

  return [
    // -------------------------------------------------- the app, at a glance
    {
      id: 'orientation-app-overview-cell-type',
      topic: 'web_app',
      name: 'app-overview-cell-type',
      describes:
        'The whole Cellucid window with the Pancreas sample coloured by '
        + 'cell_type: the sidebar shows Coloring and filtering with the '
        + 'categorical field selected and its eight-row legend, and the 3D '
        + 'embedding fills the canvas in the matching colours.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('cell_type'),
        onlySections('#coloring-filtering-section'),
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'orientation-welcome-startup',
      topic: 'web_app',
      name: 'welcome-startup',
      describes:
        'The Cellucid welcome overlay as a first visit opens it, with Choose a '
        + 'dataset and Learn more, over the Suo atlas loaded behind it.',
      dataset: 'suo',
      keepWelcome: true,
      steps: [SUO_READY, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'orientation-startup-loaded-build',
      topic: 'web_app',
      name: 'startup-loaded-build',
      describes:
        'The Suo atlas loaded after the welcome overlay is dismissed, with '
        + 'every sidebar section collapsed so the footer shows, ending in the '
        + 'Version and Build line.',
      dataset: 'suo',
      steps: [SUO_READY, onlySections(), NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },

    // --------------------------------------------- 2D: picker and navigation
    {
      id: 'core-dimension-navigation-controls-2d',
      topic: 'web_app',
      name: 'dimension-navigation-controls-2d',
      describes:
        'The Navigation block reading Planar above the Compare views block '
        + 'whose Shows picker reads UMAP · 2D: choosing a 2D embedding for a '
        + 'view switched its navigation to Planar.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        onlySections('#navigation-section', '#compare-views-section'),
        ...show('UMAP · 2D', 'planar'),
      ],
      shot: { mode: 'region', selectors: ['#navigation-controls', '#view-representation-select'], padding: 4 },
      settle,
    },

    // ---------------------------------------------------- connectivity edges
    {
      id: 'core-connectivity-edges-window',
      topic: 'web_app',
      name: 'connectivity-edges-controls',
      describes:
        'The whole window with Show connectivity edges on: the neighbour graph '
        + 'drawn as thin lines between the Pancreas cells, and the sidebar '
        + 'showing the Connectivity block with Edge color, Edge opacity, Line '
        + 'width and Max edges.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('clusters'),
        onlySections('#visualization-section'),
        { waitFor: '#connectivity-controls:not([hidden])' },
        { openSection: '#visualization-section' },
        { run: page => page.check('#toggle-connectivity') },
        { waitFor: '#connectivity-sliders:not([hidden])' },
        { scrollIntoView: '#connectivity-controls' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { attempts: 30, gapMs: 300 },
    },

    // ------------------------------------------------------------ camera path
    {
      id: 'core-camera-path-unconfigured',
      topic: 'web_app',
      name: 'camera-path-unconfigured',
      describes:
        'The Camera path panel with no saved positions: Keyframes with Save '
        + 'position and Clear all, and the Playback settings below.',
      dataset: 'pancreas',
      steps: [PANCREAS_READY, onlySections('#cinematic-camera-section')],
      shot: { mode: 'element', selector: '#cinematic-camera-section', padding: 6 },
      settle,
    },
    {
      id: 'core-camera-path-configured',
      topic: 'web_app',
      name: 'camera-path-configured-panel',
      describes:
        'The Camera path panel after two positions were saved: keyframes KF 1 '
        + 'and KF 2 with the KF 1 → KF 2 duration on auto, End with a return '
        + 'to the first keyframe, Set all with Apply, and Playback with Default '
        + 'speed, Position interpolation Catmull-Rom spline, Rotation '
        + 'interpolation SLERP (smooth), Easing Linear, Loop playback and '
        + 'Autoplay.',
      dataset: 'pancreas',
      steps: [PANCREAS_READY, ...TWO_KEYFRAMES, NOTIFICATIONS_DRAINED],
      shot: { mode: 'element', selector: '#cinematic-camera-section', padding: 6 },
      settle,
    },
    {
      id: 'core-camera-path-transport',
      topic: 'web_app',
      name: 'camera-path-transport-visible',
      describes:
        'The whole window after two camera positions were saved: the playback '
        + 'bar across the top of the canvas with its Play button and timeline '
        + 'markers, and the Camera path panel in the sidebar.',
      dataset: 'pancreas',
      steps: [PANCREAS_READY, ...colourBy('clusters'), ...TWO_KEYFRAMES, NOTIFICATIONS_DRAINED],
      shot: { mode: 'window' },
      settle,
    },

    // ------------------------------------------------ the vector field overlay
    {
      id: 'velocity-dataset-without-field',
      topic: 'vector_field_velocity',
      name: 'read-dataset-without-vector-field',
      describes:
        'The Connectivity block of the Visualization panel with the Kanemaru '
        + 'heart sample loaded, and no Vector field block under it, because the '
        + 'dataset has none.',
      dataset: 'kanemaru',
      steps: [KANEMARU_READY, onlySections('#visualization-section'), { scrollIntoView: '#points-controls' }],
      shot: { mode: 'element', selector: '#points-controls', padding: 8 },
      settle,
    },
    {
      id: 'velocity-dataset-with-field',
      topic: 'vector_field_velocity',
      name: 'read-dataset-with-vector-field',
      describes:
        'The same part of the Visualization panel with the Pancreas sample '
        + 'loaded: the Connectivity block followed by a Vector field block '
        + 'whose Show overlay box is unchecked.',
      dataset: 'pancreas',
      steps: [PANCREAS_READY, onlySections('#visualization-section'), { scrollIntoView: '#points-controls' }],
      shot: { mode: 'element', selector: '#points-controls', padding: 8 },
      rings: [{ on: '#velocity-overlay-controls' }],
      settle,
    },
    {
      id: 'velocity-window-without-overlay',
      topic: 'vector_field_velocity',
      name: '01-embedding-without-overlay',
      describes:
        'The Pancreas 2D embedding coloured by cell_type, with the Vector field '
        + 'block ringed in the sidebar and its Show overlay box unchecked.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('cell_type'),
        ...show('UMAP · 2D', 'planar'),
        onlySections('#visualization-section'),
        { scrollIntoView: '#velocity-overlay-controls' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      rings: [{ on: '#velocity-overlay-controls' }],
      settle,
    },
    {
      // The particles move for as long as the overlay is on, so the image is
      // pinned to a number of frames from a seeded start, not settled.
      id: 'velocity-window-with-overlay',
      topic: 'vector_field_velocity',
      name: '02-embedding-with-overlay',
      describes:
        'The same 2D embedding with Show overlay checked: bright particle '
        + 'trails stream across the cells, and the Vector field controls are open'
        + ' in the sidebar, ending in the line that names the embeddings the '
        + 'field is available on.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('cell_type'),
        ...show('UMAP · 2D', 'planar'),
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        { scrollIntoView: '#velocity-overlay-controls' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { frames: 180 },
    },
    {
      id: 'velocity-core-controls',
      topic: 'vector_field_velocity',
      name: 'set-core-overlay-controls',
      describes:
        'The Vector field block with Show overlay checked: the Vector field '
        + 'select on velocity, the Particle density, Flow speed, Trail length, '
        + 'Particle size and Opacity sliders with their values, Color scheme, '
        + 'Sync with LOD, and the collapsed Advanced visual settings.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        { scrollIntoView: '#velocity-overlay-controls' },
      ],
      shot: { mode: 'element', selector: '#velocity-overlay-controls', padding: 8 },
      settle: { attempts: 6, gapMs: 200 },
    },
    {
      id: 'velocity-window-overlay-2d',
      topic: 'vector_field_velocity',
      name: 'overlay-controls',
      describes:
        'The Pancreas sample in 2D with the velocity overlay running over the '
        + 'cells and the Vector field controls open in the sidebar.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('clusters'),
        ...show('UMAP · 2D', 'planar'),
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        { scrollIntoView: '#velocity-overlay-settings' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { frames: 180 },
    },
    {
      id: 'velocity-window-overlay-3d',
      topic: 'vector_field_velocity',
      name: 'pancreas-velocity-3d',
      describes:
        'The Pancreas sample in its 3D view coloured by clusters with the '
        + 'velocity overlay running, and the Vector field controls open in the '
        + 'sidebar.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('clusters'),
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        { scrollIntoView: '#velocity-overlay-controls' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { frames: 180 },
    },
    {
      id: 'velocity-window-advanced',
      topic: 'vector_field_velocity',
      name: 'advanced-settings',
      describes:
        'The 2D velocity overlay with Advanced visual settings open in the sidebar: '
        + 'Particle rendering, Trail settings, HDR and bloom, Color grading '
        + 'and Cinematic effects.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        ...colourBy('clusters'),
        ...show('UMAP · 2D', 'planar'),
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        { click: '[aria-controls="velocity-advanced-settings"]' },
        { waitFor: '#velocity-advanced-settings:not([hidden])' },
        { scrollIntoView: '#velocity-advanced-settings' },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { frames: 180 },
    },

    // ----------------------------------------------------------- figure export
    {
      id: 'export-panel-overview',
      topic: 'figure_export',
      name: 'panel-overview',
      describes:
        'The Figure export panel with its four collapsed groups, Framing, '
        + 'Labels and annotations, Style and Download, above the Export button, '
        + 'with the pointer on Export.',
      dataset: 'pancreas',
      steps: EXPORT_PANEL,
      shot: { mode: 'element', selector: '#figure-export-section', padding: 6 },
      cursor: { on: '#figure-export-btn', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'export-framing',
      topic: 'figure_export',
      name: 'framing',
      describes:
        'The Framing group open: Plot size on Screen (half) with its W and H '
        + 'boxes, and Preview with its Show preview checkbox, above the Export '
        + 'button.',
      dataset: 'pancreas',
      steps: [...EXPORT_PANEL, openExportGroup('Framing')],
      shot: { mode: 'element', selector: '#figure-export-section', padding: 6 },
      settle,
    },
    {
      id: 'export-labels',
      topic: 'figure_export',
      name: 'labels-annotations',
      describes:
        'The Labels and annotations group open: the Title box pre-filled, '
        + 'starting with the dataset name; the Annotations checkboxes Axes, '
        + 'Legend, Reference grid, 3D orientation and Depth sort; the X and Y '
        + 'axis labels; Legend: Right; and Emphasize selection.',
      dataset: 'pancreas',
      steps: [...EXPORT_PANEL, openExportGroup('Labels and annotations')],
      shot: { mode: 'element', selector: '#figure-export-section', padding: 6 },
      settle,
    },
    {
      id: 'export-style',
      topic: 'figure_export',
      name: 'style',
      describes:
        'The Style group open: background, font family and the text sizes for '
        + 'the base text, legend, ticks, axis labels, title and centroid labels.',
      dataset: 'pancreas',
      steps: [
        ...EXPORT_PANEL,
        openExportGroup('Style'),
        // The open group is taller than the room below the collapsed sections
        // above it, so the section is scrolled to the top of the sidebar.
        { run: page => page.evaluate(() => document.querySelector('#figure-export-section')
          .scrollIntoView({ block: 'start' })) },
      ],
      shot: { mode: 'element', selector: '#figure-export-section', padding: 6 },
      settle,
    },
    {
      id: 'export-download',
      topic: 'figure_export',
      name: 'download',
      describes:
        'The Download group open with the format set to SVG: the SVG point '
        + 'strategy select on its placeholder, shortened to "Choose S…", and '
        + 'Export all views (split-view) with the pointer on its checkbox, above '
        + 'the Export button.',
      dataset: 'pancreas',
      steps: [
        ...EXPORT_PANEL,
        openExportGroup('Download'),
        // SVG, so the point-strategy select (which PNG hides) is in the frame.
        { select: '#figure-export-section select[aria-label="Download format"]', value: 'svg' },
        { waitFor: '#figure-export-section select[aria-label="SVG point strategy"]' },
      ],
      shot: { mode: 'element', selector: '#figure-export-section', padding: 6 },
      cursor: { on: '#figure-export-all-views', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'export-preview-window',
      topic: 'figure_export',
      name: 'preview',
      describes:
        'The whole window with Preview on in the Framing group: the sidebar '
        + 'shows the figure as it will be written, beside the live view on the '
        + 'canvas.',
      dataset: 'pancreas',
      steps: [
        ...EXPORT_PANEL,
        openExportGroup('Framing'),
        { run: page => page.check('#figure-export-preview-enabled') },
        { wait: 3000 },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { attempts: 30, gapMs: 300 },
    },
    {
      id: 'export-blocked-dialog',
      topic: 'figure_export',
      name: 'read-export-blocked-dialog',
      describes:
        'The Export blocked dialog over the app: "Cellucid will not create an '
        + 'export that differs from the active view", Velocity overlay not '
        + 'exported with its reason, and a single Back button.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        onlySections('#visualization-section'),
        ...VELOCITY_ON,
        onlySections('#figure-export-section'),
        { click: '#figure-export-btn' },
        { waitFor: '.figure-export-warning' },
      ],
      shot: { mode: 'window' },
      settle: { frames: 120 },
    },

    // -------------------------------------------------------------- benchmark
    // The live readout is hidden for the reason `benchmark-panel-controls`
    // gives: it measures the capture host, not Cellucid.
    {
      id: 'benchmark-window-synthetic',
      topic: 'benchmarking_performance',
      name: 'benchmark-panel',
      describes:
        'The whole window after the Benchmark panel generated one million '
        + 'points of the Model surface pattern: the synthetic point cloud on '
        + 'the canvas and the panel\'s presets, point count, pattern and '
        + 'embeddings controls in the sidebar; the live readout is omitted.',
      dataset: 'pancreas',
      steps: [
        PANCREAS_READY,
        onlySections('#benchmark-section'),
        { run: page => page.addStyleTag({ content: '#benchmark-stats { display: none !important; }' }) },
        { click: '#benchmark-section .benchmark-preset[data-count="1000000"]' },
        { select: '#benchmark-pattern', value: 'glb' },
        { click: '#benchmark-run' },
        { run: page => page.waitForFunction(
          () => /1,000,000/.test(document.getElementById('filter-count')?.textContent ?? ''),
          undefined,
          { timeout: 180_000 }
        ) },
        NOTIFICATIONS_DRAINED,
      ],
      shot: { mode: 'window' },
      settle: { attempts: 30, gapMs: 300 },
    },
  ];
}

/** Settles once the Python process has coloured and highlighted the viewer. */
let pythonViewerDriven = null;

/**
 * A Python process (`$CELLUCID_DOCS_PYTHON`, else `python3`) that shows the
 * Pancreas sample's store with `cellucid.show_anndata()` — outside a notebook
 * it prints the viewer URL instead of an iframe — and, once the browser
 * connects, colours it and highlights the Ductal cells from Python. Returns
 * the viewer URL; `pythonViewerDriven` settles when both commands were sent.
 * Stopped when this process exits.
 *
 * @param {{appOrigin: string}} urls
 * @returns {Promise<string>}
 */
async function startPythonViewer(urls) {
  const { spawn } = await import('node:child_process');
  const { mkdtemp } = await import('node:fs/promises');
  const { tmpdir } = await import('node:os');
  const { dirname, join, resolve } = await import('node:path');
  const { fileURLToPath } = await import('node:url');
  const here = dirname(fileURLToPath(import.meta.url));
  const datasets = process.env.CELLUCID_DATASETS_REPO
    ?? resolve(here, '..', '..', '..', '..', 'cellucid-datasets');
  const store = join(datasets, 'anndata', 'pancreas.zarr');
  const cache = await mkdtemp(join(tmpdir(), 'cellucid-docs-web-'));
  const script = [
    'import time',
    'import anndata as ad',
    'import numpy as np',
    'import cellucid',
    `viewer = cellucid.show_anndata(${JSON.stringify(store)}, web_source_url=${JSON.stringify(await webBuildOrigin())}, web_cache_dir=${JSON.stringify(cache)})`,
    'print("VIEWER_URL", viewer.viewer_url, flush=True)',
    'viewer.wait_for_ready(timeout=300)',
    'viewer.set_color_by("clusters")',
    `clusters = ad.read_zarr(${JSON.stringify(store)}).obs["clusters"]`,
    'viewer.highlight_cells(np.flatnonzero(clusters == "Ductal"), color="#e41a1c")',
    'print("DRIVEN", flush=True)',
    'while True:',
    '    time.sleep(1)',
  ].join('\n');
  const output = [];
  const child = spawn(process.env.CELLUCID_DOCS_PYTHON ?? 'python3', ['-c', script], {
    stdio: ['ignore', 'pipe', 'pipe'],
  });
  child.stdout.on('data', chunk => output.push(String(chunk)));
  child.stderr.on('data', chunk => output.push(String(chunk)));
  process.once('exit', () => child.kill());
  pythonViewerDriven = (async () => {
    const deadline = Date.now() + 300_000;
    while (!output.join('').includes('DRIVEN')) {
      if (child.exitCode !== null || Date.now() > deadline) {
        throw new Error(`The Python process did not drive the viewer:\n${output.join('')}`);
      }
      await new Promise(next => setTimeout(next, 250));
    }
  })();
  pythonViewerDriven.catch(() => {});
  const deadline = Date.now() + 300_000;
  for (;;) {
    if (child.exitCode !== null) {
      throw new Error(`The Python process exited with status ${child.exitCode}:\n${output.join('')}`);
    }
    const viewer = /VIEWER_URL (\S+)/.exec(output.join(''));
    if (viewer !== null) return viewer[1];
    if (Date.now() > deadline) {
      child.kill();
      throw new Error(`The Python process printed no viewer URL within 300 s:\n${output.join('')}`);
    }
    await new Promise(next => setTimeout(next, 250));
  }
}

let webBuildHost = null;

/**
 * This run's web repository as a web build host (`cellucid-web-assets.json`
 * and the files it lists), for the servers that download the app they serve:
 * `cellucid serve`, `show_anndata()` and `cellucid_view()` then show this
 * tree's application.
 *
 * @returns {Promise<string>}
 */
async function webBuildOrigin() {
  if (webBuildHost === null) {
    const { dirname, resolve } = await import('node:path');
    const { fileURLToPath } = await import('node:url');
    const here = dirname(fileURLToPath(import.meta.url));
    const { readFile } = await import('node:fs/promises');
    const { join } = await import('node:path');
    const { pathToFileURL } = await import('node:url');
    const web = process.env.CELLUCID_WEB_REPO ?? resolve(here, '..', '..', '..', '..', 'cellucid');
    // The inventory is computed over the files this host serves, with the
    // web repository's own generator and the asset list its committed
    // inventory names: on a clean tree it is that committed file byte for
    // byte, and on a tree with uncommitted edits it still describes what is
    // served, which the downloaders verify file by file.
    const generator = await import(
      pathToFileURL(join(web, 'scripts', 'generate-web-asset-inventory.mjs')).href
    );
    const committed = JSON.parse(await readFile(join(web, 'cellucid-web-assets.json'), 'utf8'));
    const inventory = await generator.buildWebAssetInventory({
      rootDir: web,
      paths: committed.assets.map(asset => asset.path),
    });
    // Each file goes out with the content type its inventory entry names,
    // which the downloaders check, as the published site sends them.
    const types = new Map(inventory.assets.map(asset => [asset.path, asset.content_type]));
    webBuildHost = await serveFolder(web, {
      contentTypeOf: path => types.get(path),
      extra: new Map([['cellucid-web-assets.json', {
        body: generator.serializeWebAssetInventory(inventory),
        type: 'application/json; charset=utf-8',
      }]]),
    });
  }
  return webBuildHost;
}

// ---------------------------------------------------------------------------
// The phone layout
// ---------------------------------------------------------------------------

// A function declaration, like every builder here: `SCENARIOS` is built when
// the module loads, before any `const` below it exists.
function buildPhoneScenarios() {
  // A phone held upright: 390x844 CSS px, a touch screen.
  const PHONE_VIEWPORT = Object.freeze({ width: 390, height: 844, touch: true });
  // The sidebar's open or close slide has finished.
  const SIDEBAR_SETTLED = {
    run: page => page.waitForFunction(
      () => document.getElementById('sidebar')?.getAnimations().length === 0
    ),
  };
  return [
    {
      id: 'orientation-phone-sidebar',
      topic: 'web_app',
      name: 'phone-sidebar-open',
      describes:
        'Cellucid on a phone-sized touch screen with the Pancreas sample open: '
        + 'the sidebar fills the screen, with its close button (✕) in a strip '
        + 'at the top right above the Session section.',
      dataset: 'pancreas',
      viewport: PHONE_VIEWPORT,
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        { run: page => page.waitForFunction(
          () => !document.getElementById('sidebar')?.classList.contains('hidden')
        ) },
        SIDEBAR_SETTLED,
      ],
      shot: { mode: 'window' },
    },
    {
      id: 'orientation-phone-plot',
      topic: 'web_app',
      name: 'phone-plot',
      describes:
        'The same phone screen after tapping ✕: the whole Pancreas sample in '
        + 'view on the plot, the menu button (☰) at the top left, and the '
        + 'selection bar with its Select mode switch at the bottom.',
      dataset: 'pancreas',
      viewport: PHONE_VIEWPORT,
      steps: [
        PANCREAS_READY,
        NOTIFICATIONS_DRAINED,
        SIDEBAR_SETTLED,
        { run: page => page.locator('#sidebar-toggle').tap() },
        { run: page => page.waitForFunction(
          () => document.getElementById('sidebar')?.classList.contains('hidden')
        ) },
        SIDEBAR_SETTLED,
      ],
      shot: { mode: 'window' },
    },
  ];
}
