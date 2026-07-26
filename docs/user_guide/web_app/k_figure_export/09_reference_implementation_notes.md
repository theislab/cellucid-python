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
- Layout: `cellucid/assets/js/app/ui/modules/figure-export/utils/layout.js`
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
- active view id/label,
- active field key/kind,
- filter summary,
- export options (format/size/DPI/strategy/legend/axes/background/crop),
- per-view render state + camera state snapshots.

This payload is embedded into exports:
- **PNG**: as UTF-8 `iTXt` chunks including a JSON `Comment` blob (see `buildPngTextMetadata()` in `cellucid/assets/js/app/ui/modules/figure-export/renderers/png-renderer.js`).
- **SVG**: as RDF/Dublin Core metadata plus a `cellucid:json` field (see `buildSvgMetadata()` in `cellucid/assets/js/app/ui/modules/figure-export/renderers/svg-renderer.js`).

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
- Multi-panel export uses a shared legend only when all panels share the same active field.
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
