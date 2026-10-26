# Figure export architecture

This page documents the Cellucid **Figure export** subsystem: how SVG/PNG export reproduces the current view (WYSIWYG), how large datasets are handled safely, and where to add new export renderers.

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

Every per-cell pass and every encode runs in the export worker
(`utils/figure-export-worker.js`, reached through
`utils/figure-export-worker-client.js`), never on the page's main thread:

- the point passes draw on OffscreenCanvas WebGL2 contexts the worker keeps for
  the page's lifetime (`webgl-point-rasterizer.js`), one per antialiasing
  setting; a pass returns its frame with `transferToImageBitmap()` and the
  canvas shrinks to one pixel. Creating and losing a context per pass cost a
  context and a shader compile every time, and made WebKitGTK/WPE abort the
  page on a GL call made without a current context during the teardown;
- the page copies each frame onto the plot's device pixels, then hands the
  composed canvas's pixels to the worker as an ImageData whose buffer is
  transferred (`pagePixels` in `figure-export-worker-client.js`; never an
  ImageBitmap of a page canvas, which Linux WebKit encoded without the
  canvas's pixels), where `OffscreenCanvas.convertToBlob` encodes it (the encoder
  `toBlob` uses, so the bytes are the same) and the metadata chunks are added;
- axis bounds and the vector strategies' circle markup are computed there too;
  the markup returns as one Blob, so the SVG is assembled from Blob parts and
  never exists as one string.

Each rendered file opens a worker session: a typed array is copied to the
worker once per session, and closing the session releases it. An export that is
cancelled rejects its pending requests at once; the worker's late answers are
discarded.

### Where the point pass is drawn

The point passes (`spriteZoomLimit`, `rasterizePoints`) need a WebGL2 context
on an OffscreenCanvas, and not every browser offers one in a worker:

- WebKitGTK and WPE (GNOME Web/Epiphany, and Playwright's Linux WebKit) have no
  WebGL in workers: [WebKit bug 209809](https://bugs.webkit.org/show_bug.cgi?id=209809).
- Safari 16.4–16.6, and 17.0–17.1 on macOS Ventura, returned no WebGL context
  to an OffscreenCanvas in a worker while the same context worked on the page:
  [WebKit bug 253431](https://bugs.webkit.org/show_bug.cgi?id=253431),
  [WebKit bug 183720](https://bugs.webkit.org/show_bug.cgi?id=183720).

So before the first point pass the client asks the worker once whether it can
create such a context (`webgl2`, answered by `webgl2RasterAvailable()`), and
keeps the answer for the page's lifetime. Where the worker can, every pass is
drawn there. Where it cannot, the client runs the same handler
(`createFigureExportWorkerHandler`) on the page for the point passes only; the
bounds, the vector markup and every encode stay in the worker. There is one
rasterizer and no setting: the placement follows the measured capability.
`setFigureExportWorkerWebgl2(false)` states the answer for tests, and
`CELLUCID_FIGURE_EXPORT_POINT_PASS=page` runs a browser suite with every
point pass on the page.

The PNG exporter also embeds metadata via UTF-8 PNG `iTXt` chunks:
- `Software`, `Website`, `Creation Time`, `Dataset`, `Dataset ID`, `Color Field`,
  `Source File`, `Description` and a compact JSON `Comment` blob.

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
2) Projection mismatch (the view shows another embedding, or the wrong view was exported).
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
