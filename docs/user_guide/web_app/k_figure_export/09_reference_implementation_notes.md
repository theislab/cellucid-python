# Reference (implementation notes)

**Audience:** developers and power users  
**Time:** 5–15 minutes  
**What you’ll learn:**
- Where the figure export implementation lives
- Which current architecture and authoring references define the behavior
- What to check first when debugging export mismatches

**Prerequisites:**
- Comfort reading implementation docs and browser console logs

---

## Source and documentation references

- Figure export module README: `cellucid/assets/js/app/ui/modules/figure-export/README.md`
- Figure export architecture: {doc}`../p_developer_docs/12_figure_export_architecture`
- Screenshot authoring: {doc}`../../python_package/h_developer_docs/15_docs_development_and_style_guide`

---

## Architecture (what calls what)

Figure export is split into three layers so it stays fast and maintainable:

1) **UI (sidebar panel)**  
   Collects user settings and handles preview/framing interactions.
   - `cellucid/assets/js/app/ui/modules/figure-export/figure-export-ui.js`

2) **Engine (state snapshot + orchestration)**  
   Captures view buffers/render state, builds a payload with metadata, chooses
   filenames, stages every renderer result, and triggers one download.
   - `cellucid/assets/js/app/ui/modules/figure-export/figure-export-engine.js`

3) **Renderers (SVG/PNG)**  
   Turn the payload into a Blob (no DOM manipulation for SVG; a composed
   Canvas for PNG).
   - `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-renderer.js`
   - `cellucid/assets/js/app/ui/modules/figure-export/renderers/png-renderer.js`

4) **Export worker**  
   Runs every per-cell pass and every encode off the page's main thread:
   whether a figure draws any cell, the shader-accurate WebGL2 point pass (on OffscreenCanvas), the
   vector point markup (returned as one Blob), and PNG encoding with
   `OffscreenCanvas.convertToBlob` plus the `iTXt` metadata. The page only
   composes, so exporting a large dataset does not freeze the tab. In a
   browser whose workers have no WebGL2 (WebKitGTK and WPE,
   [WebKit bug 209809](https://bugs.webkit.org/show_bug.cgi?id=209809)) the
   same handler draws the point pass on the page.
   - `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-export-worker.js`
   - `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-export-worker-client.js`

Key utilities/components:
- Layout (`computeSingleViewLayout()` and `computeGridFigureLayout()`, the whole multi-panel figure both renderers draw, with each cell's label row, plot rect and panel legend rect): `cellucid/assets/js/app/ui/modules/figure-export/utils/layout.js`
- Legends, including `planLegend()` (each entry's whole text, the columns and the box a legend needs; the builders draw the plan and measure nothing), `resolveSharedGridLegend()` (may one legend stand for the whole grid?), the `(hidden)` entry marking, and the log-scale colorbar note: `cellucid/assets/js/app/ui/modules/figure-export/components/legend-builder.js`
- The one text measurement legends are laid out with (the canvas's `measureText` for the figure's font): `cellucid/assets/js/app/ui/modules/figure-export/utils/text-measure.js`
- The PNG canvas, allocated at the figure's size and checked to really draw, or the one-sentence refusal: `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-canvas.js`
- Per-panel provenance record shared by both formats: `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-provenance.js`
- Panel letters (`A`…`Z`, `AA`…) shared by the drawn label and the metadata, and the panel title (`panelTitle`: letter and view name) that the grid layout measures and both renderers draw: `cellucid/assets/js/app/ui/modules/figure-export/utils/panel-label.js`
- `n = N selected` badge, counted per panel: `cellucid/assets/js/app/ui/modules/figure-export/utils/selection-badge.js`
- Render-mode fidelity gate: `cellucid/assets/js/app/ui/modules/figure-export/utils/render-mode.js`
- Viewer-background reader and reference-layer appearance (mirrors `viewer.setBackground()`; the background enum itself is `cellucid/assets/js/rendering/viewer-backgrounds.js`): `cellucid/assets/js/app/ui/modules/figure-export/utils/viewer-background.js`
- Reference-grid geometry, shared by SVG rules and Canvas strokes, and the entry point for either reference layer: `cellucid/assets/js/app/ui/modules/figure-export/components/reference-grid.js`
- Studio backdrop as vectors, from the constants and functions its shader is built from (`cellucid/assets/js/rendering/studio-backdrop.js`): `cellucid/assets/js/app/ui/modules/figure-export/components/studio-backdrop-export.js`
- The one placement of a view on its plot, used by every layer (points, centroids, labels, reference layers, backdrop look, preview, raster plan): the plot takes the view's or the crop's aspect (`fitPlotSize()`, used by the layout) and the view fills it (`createPlotPlacement()`, which refuses any other aspect, so nothing is letterboxed): `cellucid/assets/js/app/ui/modules/figure-export/utils/plot-fit.js`
- Axes as rulers of the plot (`computeSourceAxisBounds()`: the coordinates the plot's edges show; `planAxisTicks()`: the ticks between them whose labels, at the figure's font, never overlap; `planFigureAxes()`: every line and text of a plot's axes and the margins they need, which the layout reserves and both formats draw): `cellucid/assets/js/app/ui/modules/figure-export/components/axes-builder.js`
- The legend entry for highlighted cells (the highlight's color, the page's name, its cell count) and which panels draw them: `cellucid/assets/js/app/ui/modules/figure-export/utils/highlight-legend.js`
- Model-to-plot projection shared by both reference layers: `cellucid/assets/js/app/ui/modules/figure-export/utils/plot-projection.js`
- Data-overlay fidelity gates (connectivity, velocity) and the disclosed viewer chrome: `cellucid/assets/js/app/ui/modules/figure-export/utils/overlay-fidelity.js`
- Atmospheric fog for the vector strategies: `cellucid/assets/js/app/ui/modules/figure-export/utils/fog.js`
- Background-derived annotation ink: `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-ink.js`
- Density reduction: `cellucid/assets/js/app/ui/modules/figure-export/utils/density-reducer.js`
- Point projection: `cellucid/assets/js/app/ui/modules/figure-export/utils/point-projector.js`
- Vector point markup (written in the export worker): `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-point-markup.js`
- Session record of the panel's settings: `cellucid/assets/js/app/ui/modules/figure-export/utils/session-settings.js`
- WebGL point rasterizer (PNG + Hybrid SVG): `cellucid/assets/js/app/ui/modules/figure-export/utils/webgl-point-rasterizer.js`
- PNG UTF-8 iTXt metadata: `cellucid/assets/js/app/ui/modules/figure-export/utils/png-metadata.js`
- Deterministic stored ZIP batches: `cellucid/assets/js/app/ui/modules/figure-export/utils/zip-archive.js`
- Fidelity warnings dialog: `cellucid/assets/js/app/ui/modules/figure-export/components/fidelity-warning-dialog.js`

---

## Payload + metadata (what gets embedded)

The engine constructs a payload containing:
- dataset identity (name/id/source type/URLs/local folder or file name),
- `meta.views`: **one provenance record per exported panel** — view id/label,
  that panel's field key/kind, and that panel's filter summary — each read from
  the panel's own atomic snapshot,
- export options (format/size/DPI/strategy/legend/axes/reference grid/background/crop),
- per-view render state + camera state snapshots.

There is no top-level view/field/filters value. A grid stamped with the active
view's field and filters would attribute one panel's scientific state to the
whole figure, which is a false claim in a published figure.

This payload is embedded into exports:
- **PNG**: as UTF-8 `iTXt` chunks including a JSON `Comment` blob (see `buildPngTextMetadata()` in `cellucid/assets/js/app/ui/modules/figure-export/renderers/png-renderer.js`).
- **SVG**: as RDF/Dublin Core metadata plus a `cellucid:json` field (see `buildSvgMetadata()` in `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-renderer.js`).

Both formats build their description and JSON from the one builder in
`utils/figure-provenance.js`, so the SVG and PNG of a single export cannot
disagree. Figure-level claims — the PNG `Color Field` chunk, the SVG
`cellucid:colorField` tag, and the filename's color-field segment — are written
only when `unanimousFieldKey()` finds every panel carrying the same field.

For one job, the engine downloads the native SVG/PNG Blob. For two or more
jobs, it validates every native Blob, writes a stable-order stored ZIP with
fixed ZIP timestamps, and downloads that one `application/zip` Blob. A render
or archive failure occurs before the sole download call.

:::{important}
The metadata can include the name you gave a local folder or file
(`datasetUserPath`), and a dataset URL for remote sources. It is never a
filesystem path: `datasetUserPath` is set only from `file.name` or from the
first path segment of a chosen directory, because a browser never hands the page
an absolute location. Treat the *name* as a privacy concern when sharing figures
publicly — see {ref}`figure-export-source-file-privacy`.
:::

---

## Known limitations (current)

- Connectivity/edges overlay is not exported (point layer only).
- Only the **Points** render mode is exported; any other mode is a fidelity
  blocker, not a degraded render (`utils/render-mode.js`).
- A multi-panel export uses **one shared legend only when every panel publishes
  the same field and a semantically equal legend model, hidden categories
  included** (`resolveSharedGridLegend()`). Panels that disagree each get their
  own legend inside their cell; no colored panel is exported without a legend.
- No legend is ever capped, dropped or cut. Every legend is as large as its
  entries need, each column as wide as its longest entry; the figure grows
  around it. A grid's plots share one height (`computeGridFigureLayout`):
  the largest whole height whose composition — plots of their views' aspects,
  axis margins where a panel has axes (`figurePanelAxes`), labels, legends,
  title — fits the requested box, and the figure is exactly that composition;
  the panel's figure fills a width instead (`fillGridWidth`: no height to keep
  to), widened to it, and a single view's (`fillSingleViewWidth`) is the
  widest plot whose figure keeps to the width, each step narrowing the plot by
  what the figure is too wide (the rest never shrinks as the plot narrows, so
  it never steps past the answer). The units are one model
  (`utils/figure-size.js`): a figure pixel is 1/96 inch, the width is set in
  millimetres and text in points, a PNG at D dpi is D/96 device pixels a
  figure pixel, an SVG says millimetres. A grid keeps the on-screen
  arrangement (`payload.gridColumns`, the viewer's `getPresentedGridColumns`).
  Legends take room from the plots down to `PANEL_LEGEND_SHARE` (0.5) of the
  height they have without legends, and never below `PANEL_PLOT_MIN_WIDTH` ×
  `PANEL_PLOT_MIN_HEIGHT` (80 × 64 px); past that the figure grows. A panel
  legend stands beside its plot (Legend: Right) unless, beside the plot the
  grid has without legends, it would be wider than that plot; then, and
  always for Legend: Bottom, below it at its column's width. A legend's box
  is its plan's (`planLegend()`), never stretched to the plot. A legend
  beside the plot or the grid is at least `LEGEND_MIN_EMS_BESIDE` (10) ems of
  its font wide (120 px at 12 px, 25 mm at 7 pt), a panel legend at least `PANEL_LEGEND_MIN_WIDTH` (96 px). These live
  in `utils/layout.js`.
- Axes are drawn only where they are rulers (`figureViewHasAxes()` in
  `components/axes-builder.js`): a 1-D or 2-D view in planar navigation. A
  3-D view, or a flat view orbited, has no axes and no axis margins.
- The vector strategies and every format's centroid discs draw each point at
  the sprite size the shader gives it at its depth
  (`createViewportPointDiameter()` in `utils/point-size.js`, which is
  `pointSpriteSize` of `rendering/point-sprite.js`: `u_pointSize` mixed with
  its perspective size by the size attenuation, clamped to 0.5 px and the
  pane's cap (`renderState.pointSizeMax`: 160 px, or the GPU's largest point)),
  times the plot's scale. The raster passes run the viewer's own shaders,
  with the cap times the raster's scale. A cell near the camera is drawn
  whole and opaque in every path, as on screen (no fade by distance).
- The one limit on a figure's size is the browser's canvas (PNG): a figure it
  cannot draw is refused before rendering with one sentence naming its size
  and the remedy (`utils/figure-canvas.js`).
- SVG “full vector” can be unusable for large point clouds; Hybrid/PNG are the intended practical paths.

---

## What to capture when filing a bug (minimum repro package)

Include enough context to reproduce without guessing:

- **App version/commit** (or deployment date) if available
- **Browser + version**, **OS**, and (if relevant) **GPU model**
- **Dataset identity**:
  - dataset name/id,
  - how it was loaded (local file / server / link / Jupyter),
  - approximate point count and #categories in the active field
- **Exact steps** (click-by-click) to reach the failing state
- **Export settings**:
  - figure width (mm, or the journal column) and text sizes (pt),
  - format (SVG/PNG) and DPI,
  - strategy (Full vector, Optimized vector, or Hybrid for SVG),
  - include axes/legend/title/background/crop,
  - export all views on/off (and whether you were in grid compare)
- **The exported artifact** attached (native SVG/PNG or the batch ZIP)
- **For ZIP batches**: the exact ordered entry names
- **Console logs**:
  - any `[FigureExport]` warnings/errors,
  - screenshot of fidelity warning dialog text if it appeared
- **Deterministic delivery check**:
  - requested job count,
  - observed native filename or ZIP filename,
  - ZIP entry names for a batch.

---

## Next steps

- {doc}`07_troubleshooting_figure_export` (user-facing symptom mapping)
- {doc}`05_metadata_and_provenance` (what should be recorded/embedded)
