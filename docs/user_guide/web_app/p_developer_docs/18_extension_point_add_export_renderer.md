# Extension point: add an export renderer

This guide explains how to extend the **Figure Export** subsystem by adding a new renderer or export format option.

Before you start, read:
- {doc}`12_figure_export_architecture`

---

## Step 0: Decide what you’re adding

Common additions:

1) A new **SVG mode** (variant of full/optimized/hybrid)
2) A new **PNG export option** (metadata, sizing, sampling strategy)
3) A new **format** (only if the browser/export stack supports it cleanly)

Design constraints:
- export must be **WYSIWYG** relative to the current interactive view
- export must not regress interactive performance (work should run only on preview/export)
- large datasets must have explicit safety controls (avoid freezing the tab)

---

## Step 1: Locate the export module entry points

Start at:
- `cellucid/assets/js/app/ui/modules/figure-export/index.js`
- `cellucid/assets/js/app/ui/modules/figure-export/figure-export-engine.js`

Renderer implementations live under:
- `cellucid/assets/js/app/ui/modules/figure-export/renderers/`

Shared utilities live under:
- `cellucid/assets/js/app/ui/modules/figure-export/utils/`

---

## Step 2: Implement the renderer

### SVG renderers

If you add a new SVG renderer mode:
- implement it alongside existing SVG renderer logic
- reuse shared layout and legend builders so exported figures remain consistent

Be explicit about:
- whether points are emitted as vector elements or rasterized
- whether you rely on density reduction / sampling
- which metadata/provenance fields you include in the output

### PNG renderer changes

PNG export uses:
- Canvas2D for composition
- WebGL2 rasterization for point rendering (shader-accurate appearance)

If you change PNG behavior:
- ensure you preserve high-DPI export correctness
- ensure metadata injection remains valid (PNG `tEXt` chunks)

---

## Step 3: Add UI wiring

Update the export UI so users can select the new option:
- `cellucid/assets/js/app/ui/modules/figure-export/figure-export-ui.js`

Make sure:
- defaults are safe for large datasets
- “dangerous” modes require explicit confirmation (large dataset dialog)

---

## Step 4: Confirm state/view consistency (WYSIWYG)

Export must reflect:
- current camera (zoom/rotation/planar/orbit)
- current visibility (filters/outliers)
- current colors and legends
- current view id (live vs snapshot)
- current dimension level

If any of these are wrong, users will lose trust in exports.

---

## Step 5: Add performance safeguards

For large datasets:
- avoid generating enormous SVG strings without user confirmation
- use reduction/hybrid strategies
- debounce preview re-renders
- consider progressive rendering for previews (small sample first)

---

## Testing checklist

At minimum:
- small dataset SVG full vector works
- large dataset triggers large-export UI and hybrid/optimized modes complete
- PNG metadata is present (open exported PNG and inspect chunks if needed)
- exports match the interactive view for:
  - different themes
  - different colormaps
  - filtered vs unfiltered views
  - snapshot views

See:
- {doc}`14_testing_ci_and_release_process`

---

## Troubleshooting

### Symptom: “Export doesn’t match what I see”

Likely causes:
- wrong view id snapshotting (live vs snapshot)
- stale legend model
- visibility mask mismatch

Fix:
- ensure the export engine snapshots state and uses the same field/visibility logic as the viewer.

### Symptom: “Hybrid SVG points look different than the viewer”

Likely cause:
- rasterization path not using the same shader variant as interactive viewer.

Fix:
- use the WebGL point rasterizer utility (see `figure-export/README.md`).
