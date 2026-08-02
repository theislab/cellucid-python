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
//   name        File stem. The emitted file is `<topic>/<name>.png`.
//   describes   One sentence. Becomes the `:alt:` text, so it must describe
//               what is visible, not what the reader should conclude.
//   dataset     Dataset id in `cellucid-datasets/exports/datasets.json`.
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
//
// SIZE RULE (enforced in `capture.mjs`, restated here because this is the file
// people edit): the emitted width is the captured region's CSS width times two,
// capped at 1440. Do not set `emitWidth` to make an image "look right" on one
// page — the figures on a page are sized by what they contain, and a panel that
// is 246 CSS px wide is a 492 px file.
//
// READY SIGNALS that are known to work on this application:
//   #filter-count      "Showing all 3,696 points"   — data loaded and counted
//   #stats             "Points: … • Field: … • Centroids: …"  — colouring live
//   #legend            present once a field paints a legend
// `#dataset-cells` is NOT a ready signal: it abbreviates to "4K".
//
// The field selectors are rebuilt once the observation manifest arrives, which
// is after `#filter-count` settles. Always put `{ waitForStable:
// '#categorical-field' }` before selecting a field; without it the selection is
// silently discarded and the image shows "Field: None".
//
// The sidebar accordions are `<details class="accordion-section">`. On first
// load these are open — session, visualization, compare-views,
// coloring-filtering, highlighted-cells, page-analysis — and these are closed:
// community-annotation, cinematic-camera, figure-export, shortcuts, benchmark.
// Use `openSection` / `closeSection` rather than clicking a summary, so a
// scenario is idempotent.

const PANCREAS_READY = { waitForText: '#filter-count', text: 'Showing all 3,696 points' };

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
      'The Session panel showing the loaded Pancreatic endocrinogenesis ' +
      'dataset summary above the sample picker, with the pointer on the ' +
      'sample dropdown.',
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
  // section documents, and `.benchmark-stats` is `display: none` until a run
  // finishes, so the panel has two genuinely different states. Both are
  // captured: the controls a reader has to find, and the metric grid whose
  // *labels* the prose refers to. Collapsing the six sections that are open on
  // first load is what keeps the accordion inside the 1000 px viewport, so the
  // element crop is never silently clamped.

  {
    id: 'benchmark-panel-controls',
    topic: 'benchmarking_performance',
    name: 'benchmark-panel-controls',
    describes:
      'The Performance Benchmark panel before any run, showing the six ' +
      'point-count presets from 100K to 20M, the custom point count field, ' +
      'the data pattern dropdown, and the Load Synthetic Data, Copy ' +
      'Situation Report and Analyze Performance buttons.',
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
  //      run. Until that is understood, a "collapsed sidebar" scenario is not
  //      reproducible.
  //
  // A third, smaller point: ring labels render as a fixed 22 px circular badge,
  // so file paths cannot go on the image; they would need a numbered key on the
  // page. That part is fine — the key is more readable and it translates.
  //
  // Fixing this means a per-scenario viewport override plus understanding the
  // accordion re-open, both of which are changes to shared tooling. Until then
  // the page carries corrected prose tables instead, which is the defect that
  // actually mattered there: it documented 19 of 29 modules and named the wrong
  // initializer for three of them.

  // ==========================================================================
  // a_orientation · c_core_interactions · d_fields_coloring_legends ·
  // e_filtering
  //
  // All on `pancreas`: 3,696 cells, 3,753 genes, categorical
  // `clusters_coarse` / `clusters` / `cell_type`, continuous `S_score` /
  // `G2M_score`, latent outlier quantiles generated at 0.95. It is the
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

  // The field selectors are rebuilt when the observation manifest arrives,
  // which is after the point count settles; selecting before that is silently
  // discarded. `waitForStable` is what makes these two lines truthful.
  const clusters = [
    { waitForStable: '#categorical-field' },
    { select: '#categorical-field', label: 'clusters' },
    { waitForText: '#stats', text: 'Field: clusters (category)' },
    { waitFor: '#legend' },
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
        'The Coloring and Filtering panel with the clusters field selected: the '
        + 'three field selectors, the outlier slider at 100 percent, the '
        + 'eight-category legend, and the line Showing all 3,696 points above an '
        + 'empty filter list.',
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
        'The block labelled Active filters (selected view only), reading Showing '
        + 'all 3,696 points above the words No filters active.',
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
        'The information dialog beside the Active filters label, open and reading '
        + 'that filters stay with the selected panel when its coloring changes '
        + 'while other panels remain unchanged, with the pointer on the small i '
        + 'button that opened it.',
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
        + 'clusters, every legend row checked, and the sidebar reporting Showing '
        + 'all 3,696 points.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#active-filters-container' }],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-categories-hidden',
      topic: 'filtering',
      name: 'window-clusters-three-hidden',
      describes:
        'The same window after unchecking Ductal, Ngn3 low EP and Ngn3 high EP: '
        + 'those points are absent from the canvas, their legend rows are '
        + 'unchecked with grey swatches, and the sidebar reports a reduced count '
        + 'above one filter row.',
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
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-active-filters-one-row',
      topic: 'filtering',
      name: 'active-filters-one-row',
      describes:
        'The filter list holding one row that names the clusters field and the '
        + 'three categories it hides, followed by a visible-over-available cell '
        + 'count, with a checkbox at its left edge and a remove cross at its '
        + 'right.',
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
        'The categorical legend for clusters: a Show All and a Hide All button '
        + 'above eight rows, three of them unchecked with grey swatches, each row '
        + 'carrying a colour swatch, a label and a cell count.',
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
        'The Filtering block of the continuous legend with Live filtering turned '
        + 'off: the helper line reads Drag sliders, then click Filter, and the '
        + 'FILTER button beside RESET is enabled with the pointer on it.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        { run: page => page.click('#legend .legend-filter .toggle-switch') },
        { waitForText: '#legend', text: 'Drag sliders, then click Filter' },
      ],
      shot: { mode: 'element', selector: '#legend .legend-filter', padding: 6 },
      cursor: { on: '#legend .legend-filter button.primary', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'filtering-button-disabled-live-on',
      topic: 'filtering',
      name: 'filter-button-disabled-live-on',
      describes:
        'The same Filtering block with Live filtering on: the helper line reads '
        + 'Changes apply as you drag and the FILTER button is greyed out under '
        + 'the pointer.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        { waitForText: '#legend', text: 'Changes apply as you drag' },
      ],
      shot: { mode: 'element', selector: '#legend .legend-filter', padding: 6 },
      cursor: { on: '#legend .legend-filter button.primary', at: 'center', state: 'hover' },
      settle,
    },
    {
      id: 'filtering-window-continuous-range',
      topic: 'filtering',
      name: 'window-continuous-range-applied',
      describes:
        'The whole window coloured by S_score with the Min slider raised: only '
        + 'the higher-scoring cells are drawn, and the sidebar reports a reduced '
        + 'count above one S_score range row.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...sScore,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 20),
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-outlier-slider-95',
      topic: 'filtering',
      name: 'outlier-slider-95',
      describes:
        'The control labelled Outlier filter (latent space) inside Display '
        + 'options, its slider moved left of the end stop and its readout showing '
        + '95 percent, with the pointer on the slider.',
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
        'The whole window coloured by clusters with the outlier slider at its '
        + 'default 100 percent: scattered cells sit between and around the '
        + 'clusters, and the sidebar reads Showing all 3,696 points.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#outlier-filter-container' }],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-outlier-95',
      topic: 'filtering',
      name: 'window-outlier-95-applied',
      describes:
        'The same window with the outlier slider at 95 percent: the scattered '
        + 'cells between the clusters are gone, the readout beside the slider '
        + 'says 95 percent, and a clusters outlier row has appeared in the filter '
        + 'list above a reduced count.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#outlier-filter-container' },
        setRange('#outlier-filter', 95),
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
        'The filter list holding three rows of different kinds at once: a '
        + 'clusters hiding row, an S_score numeric range row and a clusters '
        + 'outlier row, where only the first two carry a cell count.',
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
        'The same three-row filter list with the clusters hiding row unchecked: '
        + 'that row is greyed and struck through while the other two stay black, '
        + 'and the count line above has risen.',
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
              '#active-filters [data-filter-id="obs-category-3"] .filter-checkbox'
            ),
        },
        { waitForText: '#active-filters', text: '(disabled)' },
      ],
      shot: { mode: 'element', selector: '#active-filters-container', padding: 8 },
      cursor: {
        on: '#active-filters [data-filter-id="obs-category-3"] .filter-checkbox',
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
        'A one-row filter list with the pointer resting on the remove cross at '
        + 'the right edge of the row.',
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
        on: '#active-filters [data-filter-id="obs-category-3"] .filter-remove-btn',
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
        'The whole window with a completely empty canvas after pressing Hide All '
        + 'in the clusters legend, the sidebar reading Showing 0 of 3,696 points '
        + 'above one clusters hiding row.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#legend' },
        { run: page => page.locator('#legend button', { hasText: 'Hide All' }).click() },
        { waitForText: '#filter-count', text: 'Showing 0 of 3,696 points' },
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-legend-category-unavailable',
      topic: 'filtering',
      name: 'legend-category-unavailable',
      describes:
        'The clusters legend under a narrow S_score range filter: most rows now '
        + 'print a visible-over-available count, and one row is greyed out with '
        + 'its checkbox and colour well disabled.',
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
        'The whole window coloured by the gene Ins1 with the Min slider raised, '
        + 'so only expressing cells are drawn and an Ins1 range row is listed in '
        + 'the filter list.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...ins1,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 3),
        { scrollIntoView: '#active-filters-container' },
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'filtering-window-gene-range-scope-lost',
      topic: 'filtering',
      name: 'window-gene-range-scope-lost',
      describes:
        'The same window after switching the active field to clusters: every cell '
        + 'is drawn again, the Ins1 row has left the filter list, and the count '
        + 'line reads Showing all 3,696 points.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...ins1,
        { scrollIntoView: '#legend' },
        setRange(MIN_SLIDER, 3),
        ...clusters,
        { waitForText: '#filter-count', text: 'Showing all 3,696 points' },
        { scrollIntoView: '#active-filters-container' },
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
        'The three field selectors stacked in order — Categorical obs holding '
        + 'clusters, Continuous obs holding None, and an empty Gene Expression '
        + 'search box — each followed by four small icon buttons.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#categorical-field-row' }],
      shot: {
        mode: 'region',
        selectors: ['#categorical-field-row', '#gene-expression-row'],
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
        'The categorical legend for clusters: the hint Click a swatch to pick a '
        + 'color, the Show All and Hide All buttons, and eight rows each holding '
        + 'a checkbox, a colour swatch, a label, a cell count, a rename pencil '
        + 'and a delete bin.',
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
        'The continuous legend for S_score: a Viridis colour bar with its two '
        + 'numeric end points, a Log color scale toggle reading Off, a Rescale '
        + 'colorbar to slider range toggle reading On, and the whole Filtering '
        + 'block beneath them.',
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
        'The colormap picker expanded under the S_score colour bar, a grid of '
        + 'named gradient swatches with Viridis marked as the active one.',
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
        'The Gene Expression search box holding the letters Rp with its dropdown '
        + 'open below, listing the gene names that contain them.',
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
        selectors: ['#gene-expression-row', '#gene-expression-dropdown'],
        padding: 8,
      },
      settle,
    },
    {
      id: 'fields-gene-search-no-match',
      topic: 'web_app',
      name: 'gene-search-no-match',
      describes:
        'The Gene Expression dropdown after typing a symbol this sample does not '
        + 'publish: a no-match line, a sentence giving the number of gene names '
        + 'the dataset publishes, and a link about why a gene may be missing.',
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
        selectors: ['#gene-expression-row', '#gene-expression-dropdown'],
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
        + 'embedding, with the continuous legend open in the sidebar.',
      dataset: 'pancreas',
      steps: [READY, ...ins1, { scrollIntoView: '#legend' }],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'fields-window-centroids',
      topic: 'web_app',
      name: 'window-centroids-and-labels',
      describes:
        'The Pancreas embedding coloured by clusters with Show centroids and Show '
        + 'labels both ticked, so each cluster carries a marker at its centre and '
        + 'its category name printed beside it.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#centroid-controls' },
        { run: page => page.check('#toggle-centroid-points') },
        { run: page => page.check('#toggle-centroid-labels') },
        { wait: 1200 },
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
        'The Navigation block with Mode set to Orbit, exposing a keyboard speed '
        + 'slider, a Google Earth-style drag checkbox and a Show orbit anchor '
        + 'checkbox.',
      dataset: 'pancreas',
      steps: [
        READY,
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
        'The Navigation block with Mode set to Planar, exposing a pan speed '
        + 'slider, a Zoom to cursor checkbox and an Invert axes checkbox.',
      dataset: 'pancreas',
      steps: [
        READY,
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
        'The Navigation block with Mode set to Free-fly, exposing look '
        + 'sensitivity and move speed sliders with numeric readouts, invert look '
        + 'and projectile checkboxes, and a Capture pointer button.',
      dataset: 'pancreas',
      steps: [
        READY,
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
        'The whole window showing the Pancreas 3D embedding inside its grid box, '
        + 'with the Dimension select reading 3D and the Navigation Mode select '
        + 'reading Orbit.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#dimension-controls' }],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'core-window-dimension-2d',
      topic: 'web_app',
      name: 'window-dimension-2d',
      describes:
        'The same window after choosing 2D: the cells now lie on a plane seen '
        + 'face-on, the Dimension select reads 2D, and the Navigation Mode '
        + 'select has followed it to Planar.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#dimension-controls' },
        { select: '#dimension-select', value: '2' },
        { wait: 2000 },
        { scrollIntoView: '#dimension-controls' },
      ],
      shot: { mode: 'window' },
      cursor: { on: '#dimension-select', at: 'right', dx: -18, state: 'press' },
      settle,
    },
    {
      id: 'core-window-dimension-1d',
      topic: 'web_app',
      name: 'window-dimension-1d',
      describes:
        'The same window after choosing 1D: the cells are spread along one '
        + 'axis, the Dimension select reads 1D, and the Navigation Mode select '
        + 'has followed it to Planar.',
      dataset: 'pancreas',
      steps: [
        READY,
        ...clusters,
        { scrollIntoView: '#dimension-controls' },
        { select: '#dimension-select', value: '1' },
        { wait: 2000 },
        { scrollIntoView: '#dimension-controls' },
      ],
      shot: { mode: 'window' },
      settle,
    },
    {
      id: 'core-render-mode-select',
      topic: 'web_app',
      name: 'render-mode-select',
      describes:
        'A control labelled RENDER MODE holding a dropdown that reads Points, '
        + 'ringed and with the pointer pressing it, above a dashed box labelled '
        + 'DEPTH PERCEPTION whose first slider, POINT SIZE (LOG), reads 0.75.',
      dataset: 'pancreas',
      steps: [READY, { scrollIntoView: '#render-mode' }],
      shot: {
        mode: 'region',
        selectors: ['label[for="render-mode"]', '#depth-controls'],
        padding: 8,
      },
      cursor: { on: '#render-mode', at: 'right', dx: -18, state: 'press' },
      settle,
    },

    // ------------------------------------------------------------- orientation
    {
      id: 'orientation-window-map',
      topic: 'web_app',
      name: 'window-orientation-map',
      describes:
        'The whole application window with the Pancreas sample loaded, the '
        + 'sidebar scrolled to Coloring and Filtering, the embedding coloured by '
        + 'clusters, and the eight-row legend beside the matching colours on the '
        + 'canvas.',
      dataset: 'pancreas',
      steps: [READY, ...clusters, { scrollIntoView: '#legend' }],
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
        'The Keyboard Shortcuts panel open in the Cellucid sidebar, listing '
        + 'five groups of shortcuts: the ungrouped global keys R, F, H, 1/2/3, '
        + 'O/P/G and ?; Navigation (Orbit / Planar / Free-fly) with W A S D and '
        + 'Shift; Orbit and Planar only with the zoom keys; Free-fly specific '
        + 'with E/Q and Esc; and Highlighting with Alt+click, Alt+drag, '
        + 'Shift+Alt+click, Ctrl+Alt+click, Esc and X.',
      dataset: 'pancreas',
      steps: [
        READY,
        // `background: grid` is not the built-in default, so it is only true
        // once the sample's published starting state has been applied.
        { run: page => page.waitForFunction(
          () => document.querySelector('#background-select')?.value === 'grid'
        ) },
        { closeSection: '#session-section' },
        { closeSection: '#visualization-section' },
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
   * returns to its default, then restore that file through the real Load State
   * picker. One scenario, one context: a bundle saved here is not visible to
   * any other scenario.
   */
  const saveThenRestore = () => ({
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
      await page.waitForFunction(() => (
        document.querySelector('#filter-count')?.textContent
          ?.includes('Showing all 3,696 points') === true
      ));
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
        'The Community Annotation sidebar panel before any repository is '
        + 'connected, reading "No annotation repo connected." above a Connect '
        + 'repo button, with the pointer on that button.',
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
        'Step 1 of 4 of the GitHub sync wizard, headed "Sign in with GitHub", '
        + 'with a Continue with GitHub button and a header showing the dataset '
        + 'plus GitHub and Repo both reading not connected.',
      dataset: 'pancreas',
      setup: mock({ signedIn: false, bindRepo: false }),
      steps: [
        ...openPanel,
        { click: `${PANEL} button:has-text("Connect repo")` },
        { waitFor: MODAL },
        { waitForText: MODAL, text: 'Sign in with GitHub' },
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
        'Step 2 of 4 of the wizard, headed "Install the GitHub App", listing '
        + 'repository cards marked Public or Private with Add repo and Reload '
        + 'buttons above and a pagination status line below.',
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
        'Step 3 of 4 of the wizard with a filter typed that matches no '
        + 'repository, showing "No repositories match this filter." above a '
        + 'status line reading "No repositories to display."',
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
        'Step 4 of 4 of the wizard, headed "Sync (pull / publish)", with Pull '
        + 'latest and Publish buttons, an Auto pull selector, and a status line '
        + 'reading "Pulled latest annotations."',
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
        'The connected Community Annotation panel showing Dataset, Github and '
        + 'Repo rows, a GitHub sync button, a Profile block naming the '
        + 'signed-in user, and the Manage annotation, Derived consensus column '
        + 'and Consensus snapshot sections.',
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
        'The Manage annotation section expanded, showing the categorical obs '
        + 'picker with Add, Remove and Close buttons and the Annotatable '
        + 'consensus settings block with a Threshold slider and a Min '
        + 'annotators field.',
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
        { click: `${PANEL} .analysis-accordion-header:has(.analysis-accordion-title:text-is("MANAGE ANNOTATION"))` },
        // The consensus settings block belongs to whichever column is chosen
        // in the manage picker, so it does not appear until one is chosen.
        { waitForStable: `${PANEL} select[aria-label="Categorical obs"]` },
        { select: `${PANEL} select[aria-label="Categorical obs"]`, label: '\u{1F5F3}\uFE0F clusters' },
        { waitForText: PANEL, text: 'Annotatable consensus settings' },
      ],
      shot: { mode: 'element', selector: PANEL, padding: 10 },
    },

    // ------------------------------------------------ 8. the votable badge
    {
      id: 'annotation-lifecycle-08-votable-badge',
      topic: 'community_annotation',
      name: 'lifecycle-08-votable-badge',
      describes:
        'The Coloring and Filtering categorical obs dropdown with the clusters '
        + 'entry prefixed by a ballot-box badge marking it open for voting.',
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
      shot: { mode: 'region', selectors: ['#categorical-field', '#legend'], padding: 8 },
    },

    // -------------------------------------------------- 9. Pending
    {
      id: 'annotation-lifecycle-09-pending',
      topic: 'community_annotation',
      name: 'lifecycle-09-pending',
      describes:
        'The community voting modal for the Beta cluster with a status banner '
        + 'reading Pending and a voter count below the column minimum, above a '
        + 'single suggestion card.',
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
        'The community voting modal for the Beta cluster with a status banner '
        + 'reading Consensus, the winning label, 100 percent and four voters, '
        + 'above the suggestion card carrying its ontology id, marker and '
        + 'evidence lines.',
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
        'The community voting modal for the Beta cluster with a status banner '
        + 'reading Disputed and naming two tied labels, above two suggestion '
        + 'cards carrying identical net scores and identical up and down vote '
        + 'counts.',
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
        at: 'center', state: 'hover',
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
        { fill: '.community-annotation-modal-overlay [placeholder="Ontology id (optional, e.g. CL:0000625)"] >> nth=-1', value: 'CL:0000173' },
        { fill: '.community-annotation-modal-overlay [placeholder="Marker genes (optional, comma-separated)"] >> nth=-1', value: 'Sst, Rbp4' },
      ],
      shot: {
        mode: 'element',
        selector: '.community-annotation-new >> nth=-1',
        padding: 8,
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
        'Step 4 of 4 of the wizard after publishing, with the status line '
        + 'reading "Publish complete." beside the Pull latest and Publish '
        + 'buttons.',
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
    // controls: `Save State` produces a download, and `Load State` opens a file
    // chooser that Playwright answers. The whole round trip runs inside one
    // scenario, because the tool gives each scenario a fresh browser context
    // and a bundle saved in one context is not visible to another.

    {
      id: 'sessions-save-restore-01-session-controls',
      topic: 'sessions_sharing',
      name: 'save-restore-01-session-controls',
      describes:
        'The Session accordion showing the loaded dataset summary, the sample '
        + 'picker, the local, remote-server and GitHub loading controls, and '
        + 'the Session state block with Save State and Load State at the bottom.',
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
        'The Session state block with Save State and Load State side by side '
        + 'and the pointer on Save State.',
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
        'The renderer settings block at its defaults, with Shader quality set '
        + 'to Full, Level-of-Detail unchecked and Frustum culling unchecked.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
      ],
      shot: { mode: 'element', selector: '#renderer-controls', padding: 8 },
    },

    {
      id: 'sessions-renderer-controls-changed',
      topic: 'sessions_sharing',
      name: 'save-restore-06-renderer-controls-changed',
      describes:
        'The renderer settings block with Shader quality set to Light, '
        + 'Level-of-Detail checked, Force LOD level set to 4 and Frustum '
        + 'culling checked.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        ...RENDERER_CHANGES,
      ],
      shot: { mode: 'element', selector: '#renderer-controls', padding: 8 },
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
        'The renderer settings block after restoring a saved session, showing '
        + 'Shader quality Light, Level-of-Detail checked, Force LOD level 4 and '
        + 'Frustum culling checked, matching the values that were saved.',
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
      shot: { mode: 'element', selector: '#renderer-controls', padding: 8 },
    },

    {
      id: 'sessions-save-confirmation',
      topic: 'sessions_sharing',
      name: 'save-restore-03-save-confirmation',
      describes:
        'A notification reading "Session saved successfully" after Save State '
        + 'produced a session bundle.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
        { run: async page => {
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
        'Notifications after a completed restore, ending with "Session fully '
        + 'restored." and "Session loaded successfully".',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#visualization-section' },
        { scrollIntoView: '#renderer-controls' },
        ...RENDERER_CHANGES,
        { openSection: '#session-section' },
        saveThenRestore(),
        { waitForText: '#notification-center', text: 'Session loaded successfully' },
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
        + '3-D orbit view on a light grid background, with the Session '
        + 'accordion open beside it.',
      dataset: 'pancreas',
      steps: [
        READY,
        { run: async page => {
          // Every claim the caption makes, checked against the running app
          // rather than against the prose that describes the published state.
          const applied = await page.evaluate(() => ({
            dimension: document.querySelector('#dimension-select')?.value ?? null,
            navigation: document.querySelector('#navigation-mode')?.value ?? null,
            theme: document.querySelector('#theme-select')?.value ?? null,
            background: document.querySelector('#background-select')?.value ?? null,
          }));
          const expected = {
            dimension: '3', navigation: 'orbit', theme: 'light', background: 'grid',
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
      ],
      shot: { mode: 'window' },
    },

    {
      id: 'sessions-refusal-not-a-session-file',
      topic: 'sessions_sharing',
      name: 'refusal-02-not-a-session-file',
      describes:
        'Two notifications reading "Invalid session file (bad MAGIC header)." '
        + 'after a file that is not a session bundle was chosen in Load State.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
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
              ?.includes('Invalid session file') === true
          ));
        } },
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
    },

    {
      id: 'sessions-refusal-different-dataset',
      topic: 'sessions_sharing',
      name: 'refusal-01-different-dataset',
      describes:
        'Notifications reporting that the session was saved on a different '
        + 'dataset than the one open now, naming the cell and gene counts the '
        + 'file was saved against and the counts of the dataset currently open.',
      dataset: 'pancreas',
      steps: [
        READY,
        { openSection: '#session-section' },
        { scrollIntoView: '#session-state-controls' },
        { run: async page => {
          // Save against pancreas, switch to a genuinely different published
          // dataset, then load the pancreas bundle. The refusal is the
          // product's own identity check, not a doctored file.
          const download = page.waitForEvent('download');
          await page.click('#save-state-btn');
          const saved = await (await download).path();
          await page.selectOption('#dataset-select', 'dataset:local-demo:he');
          await page.waitForFunction(() => (
            document.querySelector('#filter-count')?.textContent
              ?.includes('Showing all 3,696 points') === false
          ), undefined, { timeout: 120_000 });
          await page.waitForTimeout(1500);
          const chooser = page.waitForEvent('filechooser');
          await page.click('#load-state-btn');
          await (await chooser).setFiles(saved);
          await page.waitForFunction(() => (
            document.querySelector('#notification-center')?.textContent
              ?.includes('saved on a different dataset') === true
          ), undefined, { timeout: 60_000 });
        } },
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
      setup: mock({ repoStatus: 404 }),
      steps: [
        ...openPanel,
        { waitForText: '#notification-center', text: 'Unable to determine your role' },
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
        await page.evaluate(() =>
          document.querySelector('#annotation-confirm-btn').click()
        );
        await page.waitForFunction(() => (
          document.querySelector('.highlight-page-tab-count')?.textContent === '262'
        ));
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
  await page.evaluate(
    id => document.querySelector(`#analysis-panel-${id} .analysis-expand-btn`).click(),
    mode
  );
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
      document.querySelector('.analysis-page-tab[data-page-id="restof__page_1"]').click()
    ) },
    // Both pages in the plot's own legend is the only honest signal that the
    // second page reached the analysis rather than only the tab strip.
    { run: page => page.waitForFunction(() => (
      [...document.querySelectorAll('#analysis-panel-detailed .legend text')]
        .map(node => node.textContent)
        .join('|') === 'Page 1|Rest of Page 1'
    )) },
    { run: page => scrollAnalysisModeToTop(page, 'detailed') },
  ];

  // Differential expression of Page 1 against its complement, run from the
  // form and left with its sidebar result on screen.
  const DIFFERENTIAL_RESULT = [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'differential') },
    { select: '#analysis-panel-differential select[aria-label="Page A"]', value: 'page_1' },
    { select: '#analysis-panel-differential select[aria-label="Page B"]', value: 'restof__page_1' },
    { select: '#analysis-panel-differential select[name="method"]', value: 'wilcox' },
    { run: page => runAnalysisAndWaitForResult(page, 'differential') },
    { run: page => scrollAnalysisModeToTop(page, 'differential') },
  ];

  /**
   * Marker discovery on `cell_type` in one of its three display modes.
   *
   * @param {'ranked'|'clustered'|'custom'} mode
   * @returns {Array<Object>}
   */
  const markerGenes = mode => [
    ...ANALYSIS_READY,
    { run: page => openAnalysisMode(page, 'genesPanel') },
    { select: '#analysis-panel-genesPanel select[name="obsCategory"]', value: 'cell_type' },
    { select: '#analysis-panel-genesPanel select[name="mode"]', value: mode },
    { select: '#analysis-panel-genesPanel select[name="method"]', value: 'wilcox' },
    { run: page => runAnalysisAndWaitForResult(page, 'genesPanel') },
    { run: page => scrollAnalysisModeToTop(page, 'genesPanel') },
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
        + 'horizontal axis carries a single category name.',
      dataset: 'pancreas',
      steps: DETAILED_CATEGORICAL_RESULT,
      shot: { mode: 'element', selector: '#analysis-panel-detailed', padding: 10 },
      cursor: { on: '#analysis-plot-type', at: 'center', state: 'hover' },
    },

    // ------------------------------------ Differential expression, in the bar
    {
      id: 'analysis-de-sidebar-result',
      topic: 'analysis',
      name: 'de-sidebar-result',
      describes:
        'A small volcano plot in the sidebar with a red cloud on the right, a '
        + 'blue cloud on the left and a grey cloud in the middle, its axes '
        + 'labelled and no legend, and an Expand button below it with a '
        + 'pointer on it.',
      dataset: 'pancreas',
      steps: DIFFERENTIAL_RESULT,
      shot: {
        mode: 'region',
        selectors: [
          '#analysis-panel-differential .analysis-preview-container',
          '#analysis-panel-differential .analysis-expand-btn',
        ],
        padding: 12,
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
        'A large modal headed DIFFERENTIAL EXPRESSION: PAGE 1 VS REST OF PAGE '
        + '1, divided into four regions: a volcano plot with labelled genes '
        + 'top-left over a legend reading Up (396), Down (701) and Not '
        + 'significant (2656) drawn below its axis title, an Export row of '
        + 'PNG, SVG and CSV buttons and a PLOT OPTIONS column on the right, a '
        + 'SUMMARY STATISTICS table bottom-left, and a STATISTICAL ANALYSIS '
        + 'panel bottom-right holding the Top Differentially Expressed Genes '
        + 'table.',
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
        'The Marker Genes panel with Mode set to Ranked Genes and the help '
        + 'line "Show top markers sorted by significance", a Wilcoxon method '
        + 'select, a Discover Markers button, and beneath it a heatmap of '
        + 'eight cell-type rows against many gene columns with a Z-score '
        + 'colour bar, its gene axis carrying evenly spaced, separated names.',
      dataset: 'pancreas',
      steps: markerGenes('ranked'),
      shot: { mode: 'element', selector: '#analysis-panel-genesPanel', padding: 10 },
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
        'The expanded Marker Genes modal headed MARKER GENES: CELL_TYPE with a '
        + 'clustered heatmap and dendrograms in the plot area, an Export row '
        + 'of PNG, SVG and CSV buttons above a PLOT OPTIONS column, a SUMMARY '
        + 'STATISTICS panel naming the category, mode, method and counts, and '
        + 'a STATISTICAL ANALYSIS panel holding the Top Marker Genes table, '
        + 'with a pointer on the CSV button.',
      dataset: 'pancreas',
      steps: [
        ...markerGenes('clustered'),
        { run: page => expandAnalysisResult(page, 'genesPanel') },
      ],
      shot: { mode: 'element', selector: '.analysis-modal-content' },
      cursor: {
        on: '.analysis-export-toolbar button[title="Download data"]',
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
        on: '.analysis-export-toolbar button[title="Download data"]',
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
      shot: { mode: 'element', selector: '#analysis-panel-detailed', padding: 10 },
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
  ];
}

/**
 * Move the expanded volcano's |log2FC| threshold slider.
 *
 * The control commits on `change`, not on `input`, so both events are sent. The
 * wait is on the counts the legend prints rather than on the slider's own
 * value, because the slider holds the number that was asked for whether or not
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
    const slider = document.querySelector(
      '.analysis-modal-options input[id$="-analysis-opt-foldChangeThreshold"]'
    );
    if (slider === null) {
      throw new Error('The expanded volcano has no |log2FC| threshold slider.');
    }
    slider.valueAsNumber = value;
    slider.dispatchEvent(new Event('input', { bubbles: true }));
    slider.dispatchEvent(new Event('change', { bubbles: true }));
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
// an `:alt:` describing a message that was never in the frame. The refusal is a
// timed state that deletes itself four seconds after it appears, and the
// capture was still working when it went. `holdTimers` below is why that cannot
// happen again, and `assertCropStillFramed` in `capture.mjs` is why a future
// variant of it would fail the run instead of shipping.
// ---------------------------------------------------------------------------

/** @returns {Array<Object>} */
function buildDataLoadingRefusalScenarios() {
  // Two screenshots back to back. Everything the settle loop would otherwise
  // wait for has already been waited for by the steps below — the dataset is
  // loaded and counted, the canvas has stopped drawing, and the notification
  // center is empty until the refusal lands in it — so a long gap buys nothing
  // and only spends the subject's lifetime.
  const REFUSAL_SETTLE = Object.freeze({ attempts: 2, gapMs: 0 });

  // Stop the clock that is about to delete the subject.
  //
  // A refusal is a *timed* state: `dataset-connections.js` reports it with
  // `notifications.fail(...)`, which reaches `_updateNotification` with the
  // type changed to `error` and no duration, so the notification center falls
  // back to `defaultDuration` — 4000 ms — and schedules its own dismissal
  // (`notification-center.js`: `defaultDuration = 4000`, `fail()`,
  // `_scheduleDismiss`). The web repository's own browser suite says the same
  // thing in `data-source-failure-modes.spec.mjs`: "Notifications auto-dismiss,
  // so a failure message has to be sampled while it is on screen."
  //
  // Budgeting against that deadline does not work, and was measured not to
  // work. Resolving the crop costs about a second on a loaded machine and each
  // screenshot about seven hundred milliseconds, so the tool needs roughly half
  // the window on an idle host and more than all of it on a busy one: the same
  // scenario captured cleanly at load average 4 and failed at load average 15,
  // with no code changed between the two runs. A capture whose correctness
  // depends on the capture host's load is not reproducible, and it is exactly
  // how `fail-h5ad-wrong-file.png` and `fail-missing-umap-embedding.png` came
  // to ship as photographs of the empty canvas.
  //
  // So the deadline is removed rather than raced. Timers scheduled from this
  // point on are pushed past the end of the capture — nothing is dropped and no
  // callback is skipped, they simply do not come due while the shutter is open.
  // It is the same kind of intervention as the seeded `Math.random`, the pinned
  // time zone and Playwright's `animations: 'disabled'`: the environment is
  // held still so that what the application draws can be photographed. Nothing
  // about the notification's content, wording or appearance is touched, and
  // `assertCropStillFramed` still proves the real element was in frame.
  const holdTimers = page => page.evaluate(() => {
    const realSetTimeout = window.setTimeout;
    const CAPTURE_HORIZON_MS = 600_000;
    window.setTimeout = function (handler, delay, ...args) {
      return realSetTimeout.call(
        window,
        handler,
        Math.max(Number(delay) || 0, CAPTURE_HORIZON_MS),
        ...args
      );
    };
  });

  // Selecting a file is a `change` on a hidden input, so the Session panel is
  // opened and scrolled to first: the notification is captured, but a reader
  // who reruns this should see the control that raised it.
  //
  // Then wait for the notification center to empty. Loading a dataset raises a
  // handful of its own — the restored session, each artifact download, the
  // field that got coloured — and every one of them expires on its own timer.
  // While they are draining, the center's height changes under a
  // `column-reverse` stack, the crop region never holds still, and the settle
  // loop is still sampling when the refusal itself expires. Starting from an
  // empty center removes the churn instead of racing it, and it is also what
  // makes the figure show one message rather than four.
  const AT_LOCAL_CONTROLS = [
    PANCREAS_READY,
    { openSection: '#session-section' },
    { scrollIntoView: '#user-data-h5ad-btn' },
    {
      run: page => page.waitForFunction(
        () => document.querySelector('#notification-center')
          ?.childElementCount === 0,
        undefined,
        { timeout: 60_000 }
      ),
    },
  ];

  /**
   * Offer one file to the `.h5ad` picker and wait for the failure text.
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
      await holdTimers(page);
      await page.locator('#user-data-h5ad-input').setInputFiles(fixture);
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
        + 'HDF5/H5AD file. Choose an AnnData .h5ad file or regenerate it with '
        + 'AnnData.", with a dismiss button at its right.',
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
          'not a valid HDF5/H5AD file'
        ),
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: REFUSAL_SETTLE,
    },

    // ------------------------------------- a valid `.h5ad` with no embedding
    {
      id: 'loading-refusal-missing-umap-embedding',
      topic: 'data_loading',
      name: 'fail-missing-umap-embedding',
      describes:
        'A Cellucid notification reading "No exact UMAP embedding found in '
        + 'obsm. Expected one or more of X_umap_1d, X_umap_2d, or X_umap_3d. '
        + 'Available obsm keys: X_pca.", with a dismiss button at its right.',
      dataset: 'pancreas',
      steps: [
        ...AT_LOCAL_CONTROLS,
        offerH5ad(
          async () => {
            const { fileURLToPath } = await import('node:url');
            const { dirname, join } = await import('node:path');
            return join(
              dirname(fileURLToPath(import.meta.url)),
              'fixtures',
              'pca-only.h5ad'
            );
          },
          'No exact UMAP embedding found in obsm'
        ),
      ],
      shot: { mode: 'element', selector: '#notification-center', padding: 10 },
      settle: REFUSAL_SETTLE,
    },
  ];
}
