# State (`DataState`) and events

This page is the **source-of-truth developer guide** for how Cellucid models app state.

If you are changing anything that affects:
- fields (color-by, legends, deleted/renamed/user-defined),
- filters/visibility/outliers,
- highlights or highlight pages,
- multiview snapshots,
- dimension switching (1D/2D/3D),
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
- `positionsArray`: normalized positions (used for smoke density; updated on dimension switches)
- `colorsArray`: RGBA color buffer (uint8 packed)
- `categoryTransparency`: per-point alpha (float32), used for filtering/visibility
- `outlierQuantilesArray`: per-point outlier quantiles (float32), used by outlier filter logic

Key invariant:
- **indices are stable**: point `i` refers to the same cell across all arrays (positions/colors/alpha/highlights).

### Field data + caches

Loaded field values are cached in bounded LRU caches to prevent unbounded memory growth:

- `_fieldDataCache`: obs fields (defaults to max 50 entries)
- `_varFieldDataCache`: gene-expression fields (defaults to max 20 entries)

These caches are shared across views and are critical for performance on large datasets.

### Multi-view contexts (live + snapshots)

Cellucid supports:
- one **live** view (`viewId = "live"`)
- many **snapshot** views (“Keep view”)

`DataState` maintains a `viewContexts` map (`viewId → context`) so each view can have:
- independent active field selections (obs/var)
- independent filters/outlier settings
- independent dimension level (embedding)
- synchronized or independent camera (camera sync is viewer-side; state tracks view id)

### Highlights (multi-page highlight system)

Highlights are modeled as:
- **highlight pages**: independent collections of highlight groups
- **highlight groups**: named groups with colors, enabled flags, and cell index membership
- one derived per-point overlay buffer: `highlightArray` (`Uint8Array`)

The key point:
- **highlights are layered on top of visibility** (filtered-out points can still “be in a highlight group”).

### Dimension switching + vector fields

Cellucid supports multi-dimensional embeddings via:
- `dimensionManager` (from the data layer): knows which embeddings exist and how to load them
- `activeDimensionLevel` (live view)
- per-view dimension levels (stored in view contexts)

Vector fields (velocity/drift overlays) are optional and use:
- `vectorFieldManager` (from the data layer): lists available fields and loads vectors per dimension

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
| `dimension:changed` | `level` (number) | Active view dimension changed |
| `vectorFields:changed` | list of available fields | Vector-field availability changed (dataset load / manager attached) |

Code pointers:
- Visibility + field events are emitted from view-context sync code:
  - `cellucid/assets/js/app/state/managers/view-context-viewer-sync.js`
- Highlight/page events:
  - `cellucid/assets/js/app/state/managers/highlight-manager.js`
- Dimension + vectorFields events:
  - `cellucid/assets/js/app/state/managers/view-manager.js`

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

- Dimension changes:
  - loads new positions for the requested dimension
  - updates `positionsArray`
  - calls `viewer.updatePositions(newPositions)`
  - triggers `dimension:changed`

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
- `state.pointCount` matches `state.colorsArray.length / 4` and `state.positionsArray.length / 3`.
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
- missing manifests (e.g. no `var_manifest.json`)

### Async dimension switching

Dimension changes are async (may involve fetching embeddings).
If you add code that awaits during a dimension change, watch for:
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

### Symptom: “Dimension selector changes, but embedding does not”

Likely causes:
- embedding for the target dimension is missing
- async load failed (network / manifest mismatch)
- viewer positions cache not updated for the active view

Fix:
- Confirm dimension metadata in `dimensionManager` and the dataset identity.
- Check Network tab for embedding fetch failures.

---

Next: {doc}`07_rendering_pipeline_webgl_and_performance_notes`.
