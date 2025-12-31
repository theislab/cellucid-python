# Figure export architecture

This page documents the Cellucid **Figure Export** subsystem: how SVG/PNG export reproduces the current view (WYSIWYG), how large datasets are handled safely, and where to add new export renderers.

## At a glance

**Audience**
- Computational users: read “Export modes” + “Troubleshooting”.
- Developers: read fully before changing rendering buffers, legend models, or export code paths.

**Time**
- 20–40 minutes

**Prerequisites**
- {doc}`06_state_datastate_and_events` (state needed for WYSIWYG)
- {doc}`07_rendering_pipeline_webgl_and_performance_notes` (render buffers and projection)

---

## Code map (where export lives)

The figure export module is a dedicated UI module under:
- `cellucid/assets/js/app/ui/modules/figure-export/`

Start with:
- `cellucid/assets/js/app/ui/modules/figure-export/README.md`

Key files (high-level):
- `ui/modules/figure-export/index.js`: module entry point (wires UI + engine)
- `ui/modules/figure-export/figure-export-ui.js`: sidebar panel + controls + preview UX
- `ui/modules/figure-export/figure-export-engine.js`: snapshots the current view state and orchestrates rendering/download

Renderers:
- `ui/modules/figure-export/renderers/svg-renderer.js`
- `ui/modules/figure-export/renderers/png-renderer.js`

---

## Design goals (why the module is structured this way)

- **WYSIWYG**: exports reproduce the current interactive view (camera/zoom/rotation/filters/colors).
- **Performance-safe**: heavy work runs only when the user explicitly previews/exports (not on every render frame).
- **Extensible**: shared layout/components are reused across SVG and PNG.

---

## Export modes (SVG)

SVG export supports multiple modes to balance fidelity vs scalability:

- **Full vector** (`full-vector`)
  - Every point is a vector element.
  - Best fidelity, worst scalability on very large datasets.

- **Optimized vector** (`optimized-vector`)
  - Uses density-preserving reduction so the SVG stays editable but smaller.

- **Hybrid** (`hybrid`)
  - Points are rasterized into an `<image>`; annotations remain vector.
  - Best for huge datasets while keeping labels/legends editable.

Large dataset safety:
- The module prompts for an explicit choice on large datasets to avoid freezing the browser.

---

## PNG export (and why it uses WebGL)

PNG export uses Canvas2D for composition, but rasterizes points via WebGL2:
- this preserves shader-accurate appearance (e.g. 3D sphere shading)
- avoids maintaining two independent point renderers

The PNG exporter also embeds metadata via PNG `tEXt` chunks:
- `Software`, `Source`, `Creation Time`, `Description`, etc.

This is useful for provenance (“how was this figure made?”) and reproducibility.

---

## Shared components (legends, axes, orientation)

The module includes shared builders:
- axes builder (1D/2D; planar camera mode)
- legend builder (categorical + continuous)
- orientation indicator (3D widget)
- centroid overlay builder

Key design constraint:
- Legends should be sourced from the same model as interactive legends so colors match exactly.

---

## Session serialization interaction

Figure Export UI state is intentionally excluded from generic UI control serialization.

Why:
- export UI includes transient “preview results” state and modals
- persisting it would make sessions fragile and confusing (“why did an old preview reopen?”)

Mechanisms:
- DOM subtree exclusion using `data-state-serializer-skip="true"` in `cellucid/index.html`
- additional skip logic inside the figure export UI module

See:
- `cellucid/assets/js/app/state-serializer/README.md`
- {doc}`10_sessions_persistence_and_serialization`

---

## Troubleshooting (export)

### Symptom: “Export is blank / missing points”

Likely causes (ordered):
1) Visibility mask was not applied consistently (filters/outliers).
2) Projection mismatch (dimension switch or view mismatch).
3) Large dataset reduction produced an unexpected downsample.

How to confirm:
- Compare export preview vs the interactive view at the moment of export.
- Confirm which export mode is selected (full/optimized/hybrid).

Fix:
- Ensure export engine snapshots the correct view id and uses the same visibility/transparency logic as the viewer.

### Symptom: “SVG export freezes the tab”

Likely cause:
- Full vector export on too many points.

Fix:
- Use optimized-vector or hybrid mode.
- Consider filtering to a subset before exporting.

---

Next: {doc}`18_extension_point_add_export_renderer` (how to add a new export renderer safely).
