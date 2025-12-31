# Extension point: add an analysis mode

This guide explains how to extend Cellucid’s **Page Analysis** subsystem.

There are multiple “kinds” of analysis extensions; pick the right one before you start:

1) Add a **new analysis UI mode** (new accordion “type” like Quick/Detailed/etc.)
2) Add a **new plot type**
3) Add a **new compute operation / transform / stat test**

## Prerequisites

Read first:
- {doc}`11_analysis_architecture`

Then decide what you’re actually adding.

---

## Option A: Add a new analysis UI mode (most common)

### Step A1: Create a new UI factory

Add a new file under:
- `cellucid/assets/js/app/analysis/ui/analysis-types/`

The factory should return an object that follows the expected UI interface.
At minimum, analysis UIs typically implement:
- `onPageSelectionChange(pageIds)`
- `destroy()` (strongly recommended)

### Step A2: Register it in `ComparisonModule`

In:
- `cellucid/assets/js/app/analysis/comparison-module.js`

Add a new registration in `_registerAnalysisTypes()`:

- pick a stable `id` (used as the mode key)
- provide a user-facing `name`
- set `factory` to your UI factory
- optionally provide `factoryOptions` if you need callbacks

Reference implementation:
- existing registrations for `simple`, `detailed`, `correlation`, `differential`, etc.

### Step A3: Ensure the container exists

The analysis UI manager uses a `containerMap` created from DOM elements with `data-mode="..."`.
Confirm your new mode has a matching container in the analysis accordion markup.

Code pointer:
- analysis accordion HTML is built inside `ComparisonModule` (see `comparison-module.js`).
- UI manager contract: `cellucid/assets/js/app/analysis/ui/analysis-ui-manager.js`

### Step A4: Wire data access through `DataLayer`

Do not read raw `DataState` internals from the UI directly.
Use:
- `cellucid/assets/js/app/analysis/data/data-layer.js`

This ensures:
- consistent missing-value handling
- consistent caching
- consistent page/highlight semantics

### Step A5: Implement cleanup (mandatory for large additions)

If your UI creates:
- Plotly plots
- WebGL resources
- large cached arrays

…you must provide a cleanup path:
- call `purgePlot(...)` when replacing plots
- unregister memory monitor handlers if you register them
- clear timers and event listeners

Memory monitor:
- `cellucid/assets/js/app/analysis/shared/memory-monitor.js`

---

## Option B: Add a new plot type

Plot types are usually registered by importing the module:
- `cellucid/assets/js/app/analysis/plots/types/*.js`

Steps:
1) Add your plot module under `analysis/plots/types/`.
2) Ensure it registers with the plot registry (see existing plot type files).
3) Import it from `comparison-module.js` (the file currently imports plot types to register them).
4) Add UI affordances if users need to choose it.

Remember:
- Plotly theming requires resolved colors (see `analysis/shared/plot-theme.js`).

---

## Option C: Add a compute operation / transform / stat test

This is the deepest extension point and the easiest to get wrong.

Steps (high level):
1) Define the operation metadata in:
   - `analysis/compute/operations.js`
2) Implement the handler in:
   - `analysis/compute/operation-handlers.js` (or a specialized module)
3) Decide where it can run:
   - GPU-capable? Worker-capable? CPU-only fallback?
4) Add validation for payload schemas (avoid silent corruption).
5) Update UI factories to use the new operation through the compute manager.

Stat tests typically live under:
- `analysis/stats/statistical-tests.js`

Transforms typically live under:
- `analysis/data/transform-pipeline.js` and registry modules.

---

## Testing and validation

Because analysis interacts with highlights, pages, datasets, and caching, validate:

- dataset switch resets caches cleanly
- page add/remove/rename does not corrupt selection state
- plots are purged when replaced
- memory monitor cleanup does not break active UI

See:
- {doc}`14_testing_ci_and_release_process`

---

## Troubleshooting

### Symptom: “My new mode never appears”

Likely causes:
- not registered in `_registerAnalysisTypes()`
- container missing in `containerMap`

Fix:
- ensure id matches `data-mode="..."` container

### Symptom: “Mode works once then gets slower”

Likely causes:
- plots not purged
- caches accumulate with no cleanup

Fix:
- add `destroy()` and purge/clear resources
- register cleanup with the memory monitor if appropriate

---

Next: {doc}`18_extension_point_add_export_renderer`.
