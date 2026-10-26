# State (`DataState`) and events

This page is the **source-of-truth developer guide** for how Cellucid models app state.

If you are changing anything that affects:
- fields (color-by, legends, deleted/renamed/user-defined),
- filters/visibility/outliers,
- highlights or highlight pages,
- multiview snapshots,
- the embedding each view shows (and so its dimension, 1D/2D/3D),
- vector-field overlays,

…you should read this page before implementing changes.

## At a glance

**Audience**
- Computational users: skim “Events” and “How to debug state”.
- Developers: read fully (this is where most subtle bugs come from).

**Time**
- 30–60 minutes (longer if you follow code pointers)

**Prerequisites**
- Familiarity with the architecture overview: {doc}`05_app_architecture_overview`

---

## `DataState` in one sentence

`DataState` is the **state coordinator** that owns core typed arrays and exposes the public state API by mixing in method surfaces from focused manager modules.

Code:
- `cellucid/assets/js/app/state/core/data-state.js`

---

## What `DataState` owns (core data + invariants)

### Core typed arrays (hot-path state)

`DataState` owns the arrays that back rendering and visibility:

- `pointCount`: number of points/cells
- `positionsArray`: the live view's normalized positions, with exactly
  `positionsDimension` components per cell (x, xy or xyz — never padded to
  XYZ); both are updated together whenever the live view shows another
  embedding, and read by smoke
  density and the session fingerprint
- `colorsArray`: RGBA color buffer (uint8 packed); smoke density reads it too,
  so the cloud carries the active field's colours
- `categoryTransparency`: per-point visibility, a `Uint8Array` holding exactly 0
  or 1 (every filter keeps or drops a cell; highlight dimming and fog are done
  in the shaders), validated by `requirePointTransparency`
  (`rendering/alpha-visibility.js`)
- outlier quantiles are not copied into the state: the outlier filter reads the
  active field's own `outlierQuantiles` through `getActiveOutlierQuantiles()`

Key invariant:
- **indices are stable**: point `i` refers to the same cell across all arrays (positions/colors/alpha/highlights).

### Field data + caches

Loaded field values live in `_fieldDataCache` (obs fields) and
`_varFieldDataCache` (gene-expression fields), shared across views. They are
bounded in **bytes**, not entries (`field/residency.js`): once inactive fields
hold more than their budget, the least recently used are released. The budget
is `residencyBudgetBytes('inactiveFields')` (`data/memory-budget.js`): a
sixteenth of the estimated device memory (see {ref}`memory-budgets`). The active field, any field with filter state, and fields still
loading are always kept, and release runs only when the active field changes, so
a session restore or an analysis that loads several fields cannot lose one.
(An entry-count bound did not bound memory: browsing genes kept every one at
4 bytes per cell.)

(memory-budgets)=
### Memory budgets

Every cache that decides what stays resident takes its budget from
`data/memory-budget.js`. A budget is a residency policy: it decides what is
released, spilled to disk or read again, never what can be opened or shown.

- **The device** is estimated from what the browser reports:
  `navigator.deviceMemory` (Chromium; earlier releases report at most 8 GiB,
  Chromium 151 reports 16 on a 16 GiB Mac) and the JavaScript heap limit
  (`performance.memory.jsHeapSizeLimit`), which V8 sets at about a quarter of
  physical memory (3.76 GB on that Mac). The heap limit raises an 8 GiB report
  that may be a cap and stands in when deviceMemory is absent. A browser that
  reports neither (Firefox, Safari) is taken to be 8 GiB.
- **The allowance** is half the device: every budget together. The other half
  is the system, the browser, other tabs and what the views show.
- **The shares** of the allowance sum to exactly one, so no two caches count
  the same bytes: the gene-column index 3/8 (it spills to disk beyond it),
  positions of layouts no view shows 1/8, inactive fields 1/8, an analysis
  run's default budget 1/8, the AnnData adapter's cache 1/16, a latent space
  held for outlier quantiles 1/16, decoded chunk blocks 1/16, decode workers
  1/32 (about 40 MB each while they decode; with the cores they size the
  pool) and selection undo history 1/32.
- **A refused allocation** of fewer bytes than the allowance halves the
  allowance of the thread that saw it for every budget read afterwards.

| Device | Allowance | Gene index | Positions, inactive fields, analysis (each) | Decode workers |
|---|---|---|---|---|
| 4 GiB | 2 GiB | 768 MiB | 256 MiB | 1 |
| 8 GiB (or unreported) | 4 GiB | 1.5 GiB | 512 MiB | 3 |
| 16 GiB | 8 GiB | 3 GiB | 1 GiB | 6 |

### Multi-view contexts (live + snapshots)

Cellucid supports:
- one **live** view (`viewId = "live"`)
- any number of snapshot views (“Keep view”)

There is no count cap: GPU memory decides what stays resident. When an
allocation is refused, the colour stores and alpha textures of kept views that
are not presented are released first and rebuilt when they are shown again
(`tests/browser/kept-views-memory.spec.mjs`). The grid layout fits any number
of panes. `createSnapshotView` refuses only in smoke render mode.

`DataState` maintains a `viewContexts` map (`viewId → context`) so each view can have:
- independent active field selections (obs/var)
- independent filters/outlier settings
- the embedding it shows, and whether it follows its dimension's default
  (`followsDefault`); its dimension is that embedding's
- its own camera and navigation mode; **Link cameras** shares one camera
  between views that show the same embedding (viewer-side:
  `viewer.setCamerasLocked()`, `viewer.getCameraLinkGroups()`)

Per-cell arrays (colors, category transparency, the visibility mask, outlier
quantiles) are owned, not copied, as views change. The **active** view's context
shares the arrays `DataState` itself holds, so syncing it after every coloring,
filtering or centroid change costs nothing per cell. `setActiveView` moves the
arrays: the outgoing context keeps the ones the state was using, and the state
adopts the incoming context's. A new snapshot (“Keep view”) **shares** the
transparency array and visibility mask with the view it was kept from, copy on
write (`state/managers/per-cell-array-ownership.js`): every writer goes through
`writableTransparency`, `writableCellVisibilityMask`, `writableColors` or
`writableColorAlphaBytes`, which copy first if another view still holds the
array. A kept view's colours are one copy that the viewer and renderer keep by
reference; once published, an array is never written again. A context handed
to `restoreContext` is copied, because its caller may still hold it. So never
write a per-cell array directly — always take it through a `writable*` helper.
A kept view now adds about 6–7.4 bytes per cell instead of 26.

### Field switch cost

Switching the colored field recolors every cell on the main thread, so its
per-cell loops are kept to table reads. Category colors come from one RGBA word
per category (`_writeCategoryColors`, shared by coloring and the color refresh
after a filter edit); continuous colors sample the colormap without allocating.
The value scans that only depend on a field's data are proved once per array:
a continuous field's range and non-finite check (`scanContinuousValues`), and a
categorical field's code-range check (`requireCategoryCodesInInventory`).
Category counts resolve the active filters once per count instead of once per
cell. Field arrays are never written in place once published — a changed
field gets a new array — which is what makes those per-array proofs valid.
At 5,000,000 cells (Chromium, Apple M1 Pro) a cached switch to a continuous
field blocked the main thread for about 480 ms and now for about 140 ms; to a
categorical field, about 400 ms and now about 180 ms.

### Highlights (multi-page highlight system)

Highlights are modeled as:
- **highlight pages**: independent collections of highlight groups
- **highlight groups**: named groups with colors, enabled flags, and cell index membership
- one derived per-point overlay buffer: `highlightArray` (`Uint8Array`)

The key point:
- **highlights are layered on top of visibility** (filtered-out points can still “be in a highlight group”).

Selections in progress (the annotation, lasso, KNN and proximity candidates)
are `CellSet`s (`app/utils/cell-set.js`): one bit a cell, with word-wise union,
intersection and difference. Never hold a set of cells in a JavaScript `Set`
or `Map`: V8 refuses more than 16,777,216 entries, so a selection of a large
cluster in a 30M-cell dataset would fail. Cell lists are checked with
`requireCellList` (the same module), which refuses duplicates without a `Set`,
for the same reason. `tests/cell-set-contract.test.mjs`
enforces this in the selection modules.

### Representations + vector fields

A view shows one *representation*: one embedding, at its columns, of any
dimension. Cellucid tracks them via:
- `dimensionManager` (from the data layer): the layout registry. Every
  position layout has a `layoutKey`; positions, normalization and cache
  eviction are all keyed by layout (`loadLayout`, `getPositions`,
  `getNormTransform`). The data mapping publishes which layouts are *enabled*
  — any number, of any mix of dimensions — and at most one *default* per
  dimension (`setRepresentations`, `getEnabledLayoutsForDimension`,
  `getDefaultLayout`). The layout a view shows is
  `dimensionManager.getViewLayout(viewId)`; nothing stores a view's dimension
  beside it (`state.getViewDimensionLevel(viewId)` reads the layout's).
- `state.getEnabledRepresentations()` lists what the picker offers, and
  `state.getDataMapping()` / `state.applyDataMapping()`
  (`state/managers/mapping-manager.js`) read and change the mapping.

Coordinates of an enabled layout are read the first time a view shows it and
kept while a view shows it; afterwards they stay in a cache whose budget is a
sixteenth of the estimated device memory
(`residencyBudgetBytes('positions')`, see {ref}`memory-budgets`), released least
recently shown first, and entirely under memory pressure. The budget never
refuses what the views show: the layouts shown plus the one being shown are
always admitted, even above it.

Vector fields (velocity/drift overlays) are optional and use:
- `vectorFieldManager` (from the data layer): lists available fields and loads
  vectors per layout (`loadField(fieldId, layoutKey)`); each field is bound to
  the embeddings it was computed on (`hasFieldLayout(fieldId, layoutKey)`), and
  is drawn in every view that shows one of them

---

## How the public API is assembled (mixins + managers)

`DataState` is a single class, but its API surface is assembled from multiple modules:

- `cellucid/assets/js/app/state/core/data-state.js`
  - constructs `DataState`
  - mixes in method prototypes:
    - `DataStateViewMethods`
    - `DataStateFieldMethods`
    - `DataStateColorMethods`
    - `DataStateFilterMethods`
    - `highlightStateMethods`
  - defines “manager getters”:
    - `state.views` → `ViewManager`
    - `state.fields` → `FieldManager`
    - `state.filters` → `FilterManager`
    - `state.colors` → `ColorManager`
    - `state.highlights` → `HighlightManager`

### Why this structure exists

- Keeps concerns separated (filters vs coloring vs view orchestration).
- Makes it easier to reason about “what should be serialized”.
- Avoids accidental cyclic imports between UI and state.

---

## Events (what the UI and analysis listen to)

`DataState` extends a small `EventEmitter` (`cellucid/assets/js/app/utils/event-emitter.js`).
Events are the main way UI modules stay decoupled from state internals.

### Core events (stable names)

| Event | Payload | Meaning |
|---|---|---|
| `visibility:changed` | none | Visibility/transparency changed (filters/outliers updated) |
| `field:changed` | `{ source, fieldIndex, changeType, detail }` | Field metadata/state changed (rename/delete/load/filter/color changes) |
| `highlight:changed` | none | Highlight groups changed (membership/colors/enabled) |
| `page:changed` | none | Highlight page set changed (add/remove/rename/switch) |
| `view:representation-changed` | `{ viewId, layoutKey, dimension, previousLayoutKey }` | The representation a view shows changed (any view, active or not), or focus moved to a view showing another one |
| `mapping:changed` | `{ mapping, cause, origin }` | The data mapping was applied; `origin` is `user` (Session › Data mapping), `session` (a session restore) or `notebook` (a host program's command) |
| `representations:changed` | `{ enabled }` | The enabled representations changed (fires with `mapping:changed`) |
| `vectorFields:changed` | list of available fields | Vector-field availability changed (dataset load / manager attached) |

Code pointers:
- Visibility + field events are emitted from view-context sync code:
  - `cellucid/assets/js/app/state/managers/view-context-viewer-sync.js`
- Highlight/page events:
  - `cellucid/assets/js/app/state/managers/highlight-manager.js`
- Representation + vectorFields events:
  - `cellucid/assets/js/app/state/managers/view-manager.js`
- Mapping events:
  - `cellucid/assets/js/app/state/managers/mapping-manager.js`

### Event ordering and “why did this fire twice?”

Common patterns:
- One user action (e.g. “switch active field”) can cause **both**:
  - `field:changed` (active field selection)
  - `visibility:changed` (if outlier defaults/filters change)
- Bulk restores (session load) may intentionally fire fewer events by using batch mode.

If you are debugging a double-render / double-recompute:
- Confirm whether the UI is subscribed in multiple modules.
- Confirm whether the state is emitting from both a “setter” and a derived recompute.

---

## Batch mode (bulk changes without repeated recompute)

When applying many state changes (especially during session restore), you want to avoid:
- re-computing visibility N times
- re-uploading buffers N times

`DataState` provides:
- `beginBatch()`
- `endBatch()`

Semantics:
- In batch mode, managers mark “dirty” flags instead of immediately recomputing.
- On `endBatch()`, a single consolidated recomputation occurs.

Code:
- `cellucid/assets/js/app/state/managers/view-manager.js`

---

## State ↔ viewer synchronization (how changes become pixels)

The renderer is updated explicitly by state methods (not by UI code directly).
Common flows:

- Active field changes:
  - loads field data (if needed)
  - updates `colorsArray`
  - calls `viewer.updateColors(colorsArray)`
  - triggers `field:changed`

- Filter/outlier changes:
  - updates `categoryTransparency`
  - calls `viewer.updateTransparency(categoryTransparency)`
  - triggers `visibility:changed`

:::{important}
The two halves of a visibility change have very different costs, and code
written as if they were the same has been a recurring source of regressions.

- **Deriving `categoryTransparency` is `O(pointCount)` every time**, no matter
  how few cells moved. Do not add another full-dataset pass beside it; fold new
  conditions into the existing loop.
- **Publishing it to the GPU is proportional to what changed.** The viewer
  compares the new alpha bytes against the uploaded texture and re-uploads only
  the dirty rows, coalesced into a bounded number of contiguous regions. A
  republication in which nothing moved uploads nothing *and* reports that fact
  back, so dependent work — notably repacking the highlight geometry — is
  skipped entirely.

Two consequences for contributors:

1. Never call `updateTransparency` with a freshly allocated array that happens to
   hold identical values expecting it to be free; it is free because the *bytes*
   match, and it is your job not to perturb them.
2. If you add a consumer that must react to visibility, hang it off the
   “generation moved” signal rather than off the call itself, or you reintroduce
   the cost the comparison exists to avoid.
:::

- Representation changes all go through
  `setViewRepresentation(viewId, layoutKey, { followsDefault })`. The callers
  are `setViewEmbedding(viewId, ref)` (an explicit choice: the picker, the
  API), `restoreViewEmbedding(viewId, ref)` (a session: the view follows the
  default again when it shows it), `showDefaultRepresentation(viewId, level)`
  and `cycleRepresentation(viewId, level)` (keys 1, 2 and 3). It:
  - refuses a layout that is not enabled; a request for the layout the view
    already shows only records `followsDefault`
  - loads that layout's positions, pinned while the change is staged, and
    stages its centroids, with rollback for every owner on failure
  - updates `positionsArray` and `positionsDimension` (live view)
  - calls `viewer.updatePositions(newPositions, dimensionLevel)` (or
    `viewer.setViewPositions` for a kept view), then
    `dimensionManager.setViewLayout(viewId, layoutKey)` and
    `viewer.setViewEmbedding(viewId, { dimension, layoutKey })`, which drops
    that view's LOD, frustum, picking and vector state for the old layout
  - triggers `view:representation-changed`

The bulk of this wiring lives under:
- `cellucid/assets/js/app/state/managers/view-context-*`
- `cellucid/assets/js/app/state/managers/*-manager.js`

---

## How to debug `DataState` quickly

### 1) Use the dev globals

In the browser console:
- `window._cellucidState` gives you the active state.
- `window._cellucidViewer` gives you the viewer.

### 2) Subscribe to events interactively

Example pattern:

```js
const state = window._cellucidState;
state.on('field:changed', (e) => console.log('field:changed', e));
state.on('visibility:changed', () => console.log('visibility:changed'));
state.on('highlight:changed', () => console.log('highlight:changed'));
```

### 3) Confirm invariants

Useful checks:
- `state.pointCount` matches `state.colorsArray.length / 4` and `state.positionsArray.length / state.positionsDimension`.
- `state.getActiveViewId()` is what you think it is (live vs snapshot).

If these invariants break, most UI behavior becomes misleading.

---

## Edge cases (common sources of subtle bugs)

### Dataset reload / switching

Risks:
- listeners are still subscribed to old state/viewer objects
- caches carry over incorrectly
- highlights/snapshot views persist unexpectedly

Mitigation patterns:
- centralize reload in `main.js`
- ensure state reset paths clear view contexts, highlights, and caches as appropriate

### Empty datasets

Cellucid supports “no dataset selected” and “empty data” states.
Make sure new code handles:
- `pointCount === 0`
- a dataset without genes (an empty `var`)

### Async representation changes

Representation changes are async (may involve fetching embeddings).
If you add code that awaits during one, watch for:
- “cross-await state corruption” (the active view changed while awaiting)

The code uses a `_dimensionChangeLock` to serialize changes; preserve that property when refactoring.

### Highlights vs filters

Users often expect highlights to disappear when filtered; Cellucid treats them as independent.
If you change highlight behavior, document:
- whether highlight membership is preserved when hidden
- whether analysis uses “visible highlights” or “all highlights”

---

## Troubleshooting

### Symptom: “Filters don’t update the viewer” / “count changed but points didn’t”

Likely causes:
- `categoryTransparency` updated but `viewer.updateTransparency` not called.
- A view context is active but you updated live-view state (or vice versa).

How to confirm:
- Log `state.getActiveViewId()` and inspect which view you are editing.
- Put a breakpoint in `viewer.updateTransparency`.

Fix:
- Ensure state methods (not UI code) own the viewer sync calls.
- Use the view-context-aware APIs when operating on snapshot views.

### Symptom: “The picker changes, but the embedding does not”

Likely causes:
- the representation is not enabled (`dimensionManager.isEnabled(layoutKey)`)
- async load failed (network / field list mismatch / values refused)
- the change targeted another view than the one on screen

Fix:
- Confirm `state.getEnabledRepresentations()` and
  `state.getViewLayoutKey(viewId)` for the view you expect.
- Check Network tab for embedding fetch failures.

---

Next: {doc}`07_rendering_pipeline_webgl_and_performance_notes`.
