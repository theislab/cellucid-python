# Analysis architecture

This page documents the Cellucid **Page Analysis** subsystem: how analysis UIs are registered, how data is queried, which compute backends exist, and how the module avoids memory/performance degradation during long sessions.

It is written for contributors adding:
- a new analysis mode,
- a new plot type,
- a new compute operation,
- or changes that affect analysis ↔ highlights ↔ sessions interactions.

## At a glance

**Audience**
- Computational users: read “What analysis depends on” + “Troubleshooting”.
- Developers: read fully and then follow {doc}`17_extension_point_add_analysis_mode` for implementation steps.

**Time**
- 30–60 minutes

**Prerequisites**
- {doc}`06_state_datastate_and_events` (highlights/pages and events)
- {doc}`09_data_loading_pipeline_and_caching` (var fields and lazy loading context)

---

## Where analysis is initialized (not a normal UI module)

Unlike most sidebar modules, Page Analysis is created directly in:
- `cellucid/assets/js/app/main.js`

It looks for a DOM container:
- `#page-analysis-section`

Then constructs:
- `createComparisonModule({ state, container })`
  - implementation: `cellucid/assets/js/app/analysis/comparison-module.js`

For debugging convenience, the instance is stored as:
- `window._comparisonModule`

:::{note}
Page Analysis has its own UI architecture under `cellucid/assets/js/app/analysis/ui/`.
It is not orchestrated by `cellucid/assets/js/app/ui/core/ui-coordinator.js`.
:::

---

## Mental model: “analysis runs on pages + highlights”

Most Cellucid analyses are framed around:

- **Highlight pages** (independent groups of highlighted cells)
- **Visibility** (filters/outliers)
- **Variables** (categorical obs, continuous obs, gene expression)

The comparison module listens to state events:
- `page:changed` (pages added/removed/renamed/switched)
- `highlight:changed` (membership changes)

Code:
- `cellucid/assets/js/app/analysis/comparison-module.js`

This is why highlight pages are treated as first-class state: they are analysis inputs, not just UI decoration.

---

## Code map (analysis directory)

Analysis lives under:
- `cellucid/assets/js/app/analysis/`

Key entry point:

- `analysis/comparison-module.js`
  - orchestrates the analysis UI and binds it to `DataState`. Modules import
    each other directly; there is no barrel `index.js` for the analysis tree.

Subdirectories (high-level):
- `analysis/core/`: plot plugin contract, validation, the plot registry
- `analysis/compute/`: compute operations and their Worker and CPU backends
- `analysis/data/`: data layer and streaming gene loader
- `analysis/stats/`: statistical tests, corrections, result formatting
- `analysis/plots/`: plot infrastructure and types
- `analysis/ui/`: UI manager, analysis-type UIs, floating windows, shared UI components
- `analysis/shared/`: memory monitor, plot theme, debug utils, cleanup utilities

---

## Core architectural pieces

### 1) `ComparisonModule` (orchestrator)

Responsibilities:
- registers analysis “types” (Quick/Detailed/Correlation/DE/Signature/Marker Genes, etc.)
- manages shared configuration (`currentConfig`)
- coordinates page selection and mode switching
- triggers plot creation and manages plot lifetimes
- restores analysis windows from their exact session descriptors

Code:
- `cellucid/assets/js/app/analysis/comparison-module.js`

### 2) `DataLayer` (unified access to state + datasets)

The analysis module needs a safe, consistent way to:
- extract per-cell values for obs/var fields
- handle missing values
- respect filtering/visibility/highlights
- cache expensive derived arrays (e.g. bulk gene extraction)

This is handled by:
- `cellucid/assets/js/app/analysis/data/data-layer.js`

### 3) Compute backends (Worker / CPU)

Compute is structured around “operations” with capability metadata: which
operations can run in a Worker and which can run on the main thread. The
operations are correlation (`COMPUTE_CORRELATION`) and differential expression
(`COMPUTE_DIFFERENTIAL`); marker-gene discovery dispatches its own
`MARKERS_*` messages to the same worker pool. There is no GPU analysis backend.

Key files:
- `cellucid/assets/js/app/analysis/compute/operations.js`
- `cellucid/assets/js/app/analysis/compute/operation-handlers.js`
- `cellucid/assets/js/app/analysis/compute/compute-manager.js`
- `cellucid/assets/js/app/analysis/compute/worker-pool.js`
- `cellucid/assets/js/app/analysis/compute/data-worker.js`

Design goal:
- select one declared, available backend before an operation starts; execution
  errors are reported directly and never switch to another backend.

### 4) Plot infrastructure (Plotly + custom helpers)

Plots are managed via:
- a plot registry / factory (`analysis/plots/plot-factory.js`)
- a Plotly loader (`analysis/plots/plotly-loader.js`) to keep initial load light
- plot types registered as modules under `analysis/plots/types/`

Theme integration:
- CSS variables must be resolved into concrete colors for Plotly.
- Centralized in `analysis/shared/plot-theme.js` and updated on `cellucid:theme-change`.

### 5) Memory monitor + cleanup

Long analysis sessions can degrade memory if plots/caches are not cleaned up.

The module includes:
- a memory monitor that triggers periodic and pressure-based cleanup
- explicit cleanup hooks for components

Code:
- `cellucid/assets/js/app/analysis/shared/memory-monitor.js`
- `cellucid/assets/js/app/analysis/shared/resource-cleanup.js`

---

## Session interaction (analysis windows + artifacts)

The session serializer can be given references to analysis managers once the module exists:
- `sessionSerializer.setAnalysisRefs({ comparisonModule, analysisWindowManager })`

This enables:
- reopening floating analysis windows after session restore
- restoring analysis artifacts represented by registered session chunks

See:
- {doc}`10_sessions_persistence_and_serialization`

---

## Performance and correctness edge cases

### 1) Variable selection + mutual exclusion

Many analysis UIs treat variable selection as mutually exclusive:
- categorical obs vs continuous obs vs gene expression

When adding new modes, be explicit about:
- which variable types are accepted
- what happens when the selection becomes invalid mid-session (field deleted, dataset switched)

### 2) Visibility vs highlight membership

Common pitfall:
- “highlight group contains 10k cells” but filters hide 9k; analysis must decide which notion it uses.

Recommendation:
- expose the choice in UI (visible-only vs all highlighted) when it affects interpretation.

### 3) Dataset reload

On dataset reload (especially “in-place” reloads for local-user), analysis must:
- clear caches tied to the old dataset
- purge old plots and detach event handlers
- reset page selection safely if page ids changed

There are dev-only helpers to sanity check this:
- see `window._comparisonModule` and local-user switch self-test wiring in `ui/modules/dataset-controls.js`.

### 4) The per-gene streaming path

Differential expression and marker discovery repeat the same steps once per
gene, so anything in that loop is multiplied by the gene count. The loop keeps
these properties, and each one is pinned by a test:

- **Decoding happens in the decode worker.** A quantized gene payload is read on
  the main thread as transferred (still gzip when the URL ends in `.gz`) and
  handed to the quantization worker, which inflates it under the exact declared
  length and dequantizes it. Both the worker and the main-thread backend decode
  through one table-driven kernel, `decodeQuantizedToFloat32` in
  `cellucid/assets/js/data/quantization-contract.js`, which is bit identical to
  evaluating `minValue + code * scale` per cell.
- **Per-run work is done once per run.** `StreamingGeneLoader` validates both
  group specifications when a run starts, not once per gene.
- **Rank tests tabulate values without a `Map`.** Wilcoxon (DE and markers),
  Spearman ranks and order statistics number distinct values through
  `DistinctValueIndex` in `compute/math-utils.js` up to 4,096 distinct
  values and sort natively above that. Results are bit identical to the
  sort-based definition (`tests/analysis-engine-rank-kernel-identity.test.mjs`,
  `tests/statistical-tests.test.mjs`).
- **Page membership is recomputed only when the page changes.**
  `DataLayer` keeps each page's membership and digest beside the structure it
  came from (the page, its group list, each group, its `enabled` flag, its
  `cellIndices` array and length) and reuses them while all of those are the
  same objects. Highlight `cellIndices` arrays are therefore never edited in
  place once published; replace the array instead.
- **Yields are time-sliced.** Long loops hand the event loop back every 12 ms
  through a posted message (`createSlicedYield` in
  `shared/concurrency-utils.js`), not after a fixed item count through
  `setTimeout(0)`, whose 4 ms clamp turned a 20,000-gene run into seconds of
  idle waiting.

---

## Troubleshooting (analysis)

### Symptom: “Analysis panel is blank”

Likely causes (ordered):
1) `#page-analysis-section` DOM element missing (HTML change).
2) Analysis module threw during init (console error).
3) Plotly failed to load (blocked network or CSP).

How to confirm:
- Console: search for the first error after `[Main]` logs.
- Check whether `window._comparisonModule` exists.

Fix:
- Restore the expected DOM container id.
- For Plotly issues: confirm `analysis/plots/plotly-loader.js` can fetch its assets under your CSP.

### Symptom: “Plots get slower over time / tab memory grows”

Likely causes:
- caches not being cleared on dataset reload
- Plotly plots not being purged (old DOM nodes retained)
- missing cleanup handler registration

How to confirm:
- Use the memory monitor debug logs (enable analysis debug via `?debug=1` or localStorage `debug=1`).
- Take a heap snapshot and look for retained plot DOM nodes.

Fix:
- Ensure `purgePlot(...)` is called when replacing plots.
- Register cleanup handlers with the memory monitor for new components.

---

Next: {doc}`12_figure_export_architecture` or {doc}`17_extension_point_add_analysis_mode`.
