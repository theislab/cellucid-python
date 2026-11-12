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
canceled (**Cancel export**, the button's label while an export runs; or
closing the panel) rejects its pending requests at once. The worker hears the
session's `cancel` ahead of its queue and skips every request of that session
still queued. The pass already running ends (one synchronous draw) and its late
answer is discarded.

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

Both placements draw the same pixels only where their two contexts are the
same: `rasterContext` (a point-pass operation, so it runs where the passes
run) describes the context of a placement (`describeRasterContext`: the
renderer and vendor, the attributes the browser granted, the drawing buffer's
samples). `figure-export.spec.mjs` compares the two descriptions before the
two rasters, and attaches both with the differing pixels when either
differs: macOS Firefox once drew 7 of 15,346 pixels up to 3 levels apart on
antialiased rims (nightly run 41269083996), with nothing recorded to say why.

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

Key design constraints:
- Legends should be sourced from the same model as interactive legends so colors match exactly.
- **A legend entry is never cut** (owner decision, 2026-11-10). A legend is as
  large as its entries need: `planLegend()` (`components/legend-builder.js`)
  plans each entry's whole text, its columns (each as wide as its longest
  entry, `(hidden)` marker included) and the box, with the one text
  measurement of `utils/text-measure.js` (the canvas's `measureText` for the
  figure's font); `utils/layout.js` reserves exactly that box and the figure
  grows around it; the SVG and Canvas builders draw the plan and measure
  nothing. Before, the box was capped at the plot's size and the builders cut
  names to their columns with "…" (Suo `cell_type`, 69 categories, at
  720 × 500; and a grid panel's legend of 69 entries drawn as swatches and
  "…").
- **Axes are laid out once, from their text** (F158, 2026-11-11).
  `planFigureAxes()` (`components/axes-builder.js`) places every axis line,
  tick, tick label and axis label of a plot and gives the margins they need on
  each side, from the measured tick labels and the fonts: the y-axis label past
  the widest y tick label, a corner tick label or an over-long axis label
  reserved where it reaches. The layout plans each plot's axes and reserves
  those margins; `renderSvgAxes()` and `drawCanvasAxes()` draw that plan and
  measure nothing. Before, the y label stood a fixed distance from the axis
  and ran into wide tick labels (17 px ticks, 20 px label: "UMAP 2" over
  "-0.5"), the margins were a fixed 62 px at any text size, and the PNG put
  the y tick labels 2 px from their ticks where the SVG put them 8 px.
- **One layout call per export.** `computeFigureLayout(payload)`
  (`utils/layout.js`) reads each panel's axes, legend and the highlight's
  legend entry from the payload and lays the figure out; the PNG and SVG
  renderers both call it, so they cannot place anything differently.
- **Highlighted cells are in the legend.** A figure draws the active highlight
  page's cells in the highlight's color (`HIGHLIGHT_COLOR`,
  `rendering/shaders/high-perf-shaders.js`, the one color the viewer, the point
  raster and the vector rings use); `utils/highlight-legend.js` gives the
  legend entry for them (that color, the page's name from
  `payload.selection.pageName`, its cell count) and which panels draw them,
  and `planLegend({ highlight })` plans it as one more legend row.
- **Grid geometry has one owner.** `computeGridFigureLayout()` lays out the
  whole multi-panel figure for both renderers. A panel legend that cannot name
  its entries in half its cell beside the plot moves below the plot; rows grow
  to the tallest need, so the figure can come out larger than the requested
  size.
- **The only size limit is the browser's canvas.** A million-category legend
  once asked for a canvas millions of pixels wide and the PNG never came; the
  cap that answered it cut every legend. Now the PNG canvas is allocated at the
  figure's size and checked by drawing a pixel in its far corner and reading
  it back (`utils/figure-canvas.js`); where the browser cannot draw it, the
  export is refused before rendering with one sentence naming the size and the
  remedy (turn the legend off or show fewer categories; or a smaller size or
  DPI when the figure is too large without its legend).

---

## State, history and sessions

The Figure export section is excluded from the sidebar's control inventory
(`data-state-serializer-skip="true"` in `cellucid/index.html`): its controls
are not settings one by one, and its transient state (a preview, an open
dialog) is never kept. Its settings are instead one record of the `settings`
state domain (`figureExport`: `exportSessionSettings()` captures it,
`importSessionSettings()` puts it back, both validated by
`figure-export/utils/session-settings.js`), so a history step and a session
file hold them with every other setting
({doc}`29_action_layer_domains_and_history`).

No action changes the panel's settings, and nothing observes a change made in
the panel: such a change is not a step of its own. The next step that captures
the `settings` domain holds it (a person's settings change recorded by
observation folds it into that step; a dispatched `setting.set` keeps it as
"Changed outside the history"), and a Save records it as "Changed outside the
history".

The figure renderer also draws the read `view.picture` (a PNG of a view,
returned as data and never downloaded), which the printable history report
uses for its pictures (`app/history/exports/report.js`). What a figure cannot
draw (the vector field's flow, the connectivity edges) is in a picture with
`source: "screen"`: the view's pane of the viewer's canvas as drawn now
(`app/reads/views.js`; the drawing buffer is kept between frames), so an
agent tuning an overlay sees what you see. A figure picture taken while an
overlay is on says so and names the screen picture.

See:
- `cellucid/assets/js/app/state-serializer/README.md` (the control inventory)
- `cellucid/assets/js/app/state/domains/settings.js`
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
