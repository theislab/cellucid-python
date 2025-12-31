# Rendering pipeline (WebGL) and performance notes

This page explains how Cellucid renders millions of points efficiently, and what changes are likely to cause performance regressions or WebGL failures.

It is written for contributors who might touch:
- `cellucid/assets/js/rendering/*`
- state→viewer synchronization code (`cellucid/assets/js/app/state/*`)
- overlay features (smoke, connectivity edges, vector fields)

## At a glance

**Audience**
- Computational users: read “Common failure modes” + “Troubleshooting”.
- Developers: read the whole page before changing renderer APIs or buffer formats.

**Time**
- 30–60 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`
- {doc}`06_state_datastate_and_events`

---

## Renderer entry point

The viewer is created in `cellucid/assets/js/rendering/viewer.js` via:
- `createViewer({ canvas, labelLayer, viewTitleLayer, sidebar, onViewFocus })`

Hard constraint:
- **WebGL2 only**. If `canvas.getContext('webgl2')` fails, the app throws early.

---

## Major render subsystems (what gets drawn)

`viewer.js` orchestrates multiple render passes/subsystems:

- **Scatter points (main cloud)**
  - Backend: `HighPerfRenderer` (`cellucid/assets/js/rendering/high-perf-renderer.js`)
  - Receives: positions (float32), colors (RGBA uint8 packed), and alpha/transparency.
  - Notes:
    - Spatial indexing is dimension-aware (binary tree / quadtree / octree), used for picking/LOD.

- **Smoke / volumetric density**
  - Backend: `SmokeRenderer` (`cellucid/assets/js/rendering/smoke-cloud/`)
  - Uses a density volume built from visible points (GPU splatting path exists via `viewer.buildSmokeVolumeGPU`).

- **Connectivity edges**
  - Backend: instanced lines + edge textures (in `viewer.js`)
  - Visibility is derived from filter alpha and LOD; edges can be capped for UI performance.

- **Highlight tools**
  - Backend: `HighlightTools` (`cellucid/assets/js/rendering/highlight-renderer.js`)
  - Handles interactive selection and highlight overlay rendering.

- **Centroids**
  - Small, separate centroid shader program (centroid count is typically small).

- **Overlays**
  - Overlay framework: `OverlayManager` + overlay context
  - Vector field / velocity overlay: `VelocityOverlay` (`cellucid/assets/js/rendering/overlays/velocity/velocity-overlay.js`)

---

## Viewer public API (how state becomes pixels)

The viewer exposes a small API surface used by `DataState` and `main.js`.

### Core dataset load

- `viewer.setData({ positions, colors, outlierQuantiles, transparency, dimensionLevel })`
  - Initializes point buffers and builds spatial indices.
  - `colors` are RGBA uint8 packed; `transparency` is a separate float array used for filtering.

### Incremental updates (critical for performance)

- `viewer.updateColors(colors)`
  - Updates color buffer without reloading positions.

- `viewer.updateTransparency(alphaArray)`
  - Updates alpha/visibility texture/buffer.
  - Also triggers highlight/overlay rebuild hooks (e.g., selection needs the new visibility mask).

- `viewer.updatePositions(positions)`
  - Used for dimension switching (same point count, different embedding).
  - Must preserve the same indices; mismatched lengths are rejected.

Design rule:
- Prefer incremental updates (`updateColors`, `updateTransparency`, `updatePositions`) over full `setData` reloads.
Full reloads are the fastest way to introduce jank on large datasets.

---

## Multiview rendering (live + snapshots)

Cellucid can render:
- one live view
- multiple snapshot views in a grid layout

Key idea:
- A snapshot view is not “a second dataset”; it is a second **view context** with its own camera/dimension/filter/field choices, rendered from the same underlying point identity index space.

The viewer:
- tracks per-view dimension levels (`viewer.setViewDimension(...)` is called from state when available)
- caches per-view positions for snapshot views
- uses optimized paths for alpha/transparency sharing so it does not re-upload N full buffers per view

---

## Overlays (vector fields) and why they’re tricky

Vector field overlays are conceptually simple (“animate particles along vectors”), but interact with:
- dimension switching (vectors are dimension-specific),
- filtering (particles should respect visibility),
- multiview (each view has its own dimension level and visibility),
- performance (particle systems can be expensive).

Implementation notes:
- The overlay is “opt-in”: it is initialized only when enabled/needed.
- When visibility changes, `viewer.updateTransparency` marks the overlay “visibility dirty” so it can lazily rebuild spawn sources.

Developer advice:
- Avoid coupling overlay state to UI controls directly; route it through state events and viewer methods.

---

## Performance footguns (common mistakes)

### 1) Accidental hot-path allocations

Avoid:
- allocating new arrays/objects inside per-frame loops
- creating new typed arrays on every slider change when a reusable scratch buffer would work

Prefer:
- reusing scratch buffers (float32/uint8) sized to `pointCount`
- debouncing expensive work in UI modules

### 2) Full buffer re-uploads when only alpha changed

Filtering often only needs alpha changes.
Do not rebuild/re-upload positions or full color buffers for filtering.

### 3) Doing DOM work in render-critical flows

The renderer should not depend on DOM reads/writes at render frequency.
UI changes should be event-driven and decoupled from the render loop.

---

## Common failure modes

### WebGL2 unavailable

Symptom:
- immediate error: “WebGL2 is required but not supported…”

Causes:
- old browser
- corporate GPU policy / disabled WebGL
- remote desktop environments

Mitigation:
- test with a modern browser first
- provide a clear error message (already does)

### WebGL context lost

Symptoms:
- canvas freezes
- DevTools shows “WebGL context lost”

Causes:
- GPU memory pressure (large datasets + high settings)
- browser tab throttling
- driver instability

Mitigation:
- reduce smoke settings / disable heavy overlays
- lower dataset size for reproduction
- use the performance troubleshooting docs: {doc}`../n_benchmarking_performance/index`

---

## Troubleshooting (renderer-level)

### Symptom: “Points disappear / everything is invisible”

Likely causes:
- `categoryTransparency` got set to ~0 for all points (filters/outlier threshold).
- Alpha array length mismatch caused viewer to reject update.

How to confirm:
- In console: check `window._cellucidState.filteredCount` and `window._cellucidState.categoryTransparency`.
- Put a breakpoint in `viewer.updateTransparency`.

Fix:
- Ensure filter logic updates transparency correctly and recomputes counts.

### Symptom: “Dimension switch breaks rendering”

Likely causes:
- `updatePositions` received wrong length (`pointCount * 3` required).
- Embedding fetch failed and returned empty/invalid buffer.

How to confirm:
- DevTools → Network: check embedding requests.
- Console: look for `[Viewer] updatePositions: position count mismatch`.

Fix:
- Validate embedding buffers before calling viewer updates.
- Ensure the dimension manager returns consistent point ordering.

---

Next: {doc}`08_ui_modules_map`.
