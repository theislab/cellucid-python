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
   Turn the payload into a Blob (no DOM manipulation for SVG; Canvas/WebGL for PNG).
   - `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-renderer.js`
   - `cellucid/assets/js/app/ui/modules/figure-export/renderers/png-renderer.js`

Key utilities/components:
- Layout (single-view layout plus `computeGridPaneLayout()`, which owns one grid cell: label row, plot rect, panel legend rect): `cellucid/assets/js/app/ui/modules/figure-export/utils/layout.js`
- Legends, including `resolveSharedGridLegend()` (may one legend stand for the whole grid?), the `(hidden)` entry marking, and the log-scale colorbar note: `cellucid/assets/js/app/ui/modules/figure-export/components/legend-builder.js`
- Per-panel provenance record shared by both formats: `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-provenance.js`
- Panel letters (`A`…`Z`, `AA`…) shared by the drawn label and the metadata: `cellucid/assets/js/app/ui/modules/figure-export/utils/panel-label.js`
- `n = N selected` badge, counted per panel: `cellucid/assets/js/app/ui/modules/figure-export/utils/selection-badge.js`
- Render-mode fidelity gate: `cellucid/assets/js/app/ui/modules/figure-export/utils/render-mode.js`
- Viewer-background reader and reference-grid appearance (mirrors `viewer.setBackground()`): `cellucid/assets/js/app/ui/modules/figure-export/utils/viewer-background.js`
- Reference-grid geometry, shared by SVG rules and Canvas strokes: `cellucid/assets/js/app/ui/modules/figure-export/components/reference-grid.js`
- Data-overlay fidelity gates (connectivity, velocity) and the disclosed viewer chrome: `cellucid/assets/js/app/ui/modules/figure-export/utils/overlay-fidelity.js`
- Atmospheric fog for the vector strategies: `cellucid/assets/js/app/ui/modules/figure-export/utils/fog.js`
- Background-derived annotation ink: `cellucid/assets/js/app/ui/modules/figure-export/utils/figure-ink.js`
- Density reduction: `cellucid/assets/js/app/ui/modules/figure-export/utils/density-reducer.js`
- Point projection: `cellucid/assets/js/app/ui/modules/figure-export/utils/point-projector.js`
- WebGL point rasterizer (PNG + Hybrid SVG): `cellucid/assets/js/app/ui/modules/figure-export/utils/webgl-point-rasterizer.js`
- PNG UTF-8 iTXt metadata: `cellucid/assets/js/app/ui/modules/figure-export/utils/png-metadata.js`
- Deterministic stored ZIP batches: `cellucid/assets/js/app/ui/modules/figure-export/utils/zip-archive.js`
- Fidelity warnings dialog: `cellucid/assets/js/app/ui/modules/figure-export/components/fidelity-warning-dialog.js`

---

## Payload + metadata (what gets embedded)

The engine constructs a payload containing:
- dataset identity (name/id/source type/URLs/user path),
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
The metadata can include a local dataset path (`datasetUserPath`) when using local file/server workflows. Treat this as a privacy concern when sharing figures publicly.
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
- A panel legend is skipped when its cell cannot host a readable one:
  `PANEL_LEGEND_MAX_SHARE` (0.5 of the cell's inner width or body height) minus
  `PANEL_LEGEND_GAP` (8 px) must reach `PANEL_LEGEND_MIN_WIDTH` (96 px, Legend:
  Right) or `PANEL_LEGEND_MIN_HEIGHT` (44 px, Legend: Bottom). These live in
  `utils/layout.js` and are exported as `LAYOUT_CONSTANTS`.
- A shared legend beside the grid is dropped, in favor of per-panel legends,
  when the grid leaves under 96 px of spare width (right) or 72 px of spare
  height (bottom).
- SVG “full vector” can be unusable for large point clouds; Hybrid/PNG are the intended practical paths.

---

## What to capture when filing a bug (minimum repro package)

Include enough context to reproduce without guessing:

- **App version/commit** (or deployment date) if available
- **Browser + version**, **OS**, and (if relevant) **GPU model**
- **Dataset identity**:
  - dataset name/id,
  - how it was loaded (export folder / server / URL / Jupyter),
  - approximate point count and #categories in the active field
- **Exact steps** (click-by-click) to reach the failing state
- **Export settings**:
  - plot size (W×H),
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
