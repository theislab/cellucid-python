# Figure Export (Publication-Grade SVG/PNG)

Cellucid’s Figure Export tools help you turn an interactive view into a **shareable, publication-grade artifact**:
vector (`.svg`) and raster (`.png`) outputs that preserve key state (camera, visibility, color-by, legends) and can be used outside the web app.

These pages are written for mixed audiences:
- **Wet lab / non-technical**: click-by-click workflows and “what success looks like”.
- **Computational users**: format choices, performance implications, reproducibility and provenance.
- **Power users / developers**: implementation notes, limitations, and debugging playbooks.

:::{tip}
If you just need a good-looking figure quickly, start with {doc}`02_export_ui_walkthrough`, then skim {doc}`04_quality_knobs_and_best_practices`.
:::

---

## Fast path (pick your goal)

| You want to… | Start with | Why | Then read |
|---|---|---|---|
| Export something for slides fast | {doc}`02_export_ui_walkthrough` | Minimal decisions, good defaults | {doc}`04_quality_knobs_and_best_practices` |
| Choose the right SVG/PNG mode | {doc}`03_export_formats_and_renderers` | Understand tradeoffs | {doc}`04_quality_knobs_and_best_practices` |
| Make a figure you can reproduce later | {doc}`01_figure_export_goals_wysiwyg_and_reproducibility` | Know what’s guaranteed | {doc}`05_metadata_and_provenance` |
| Debug “export looks wrong / too big / fails” | {doc}`07_troubleshooting_figure_export` | Symptom → diagnosis → fix | {doc}`06_edge_cases` |

---

## Related concepts

- Exporting analysis results (tables/plots): {doc}`../h_analysis/09_exporting_analysis_results`
- Saving/restoring state before exporting: {doc}`../l_sessions_sharing/index`
- Color/legend choices affect exports: {doc}`../d_fields_coloring_legends/index`
- Filtering/visibility affects what gets exported: {doc}`../e_filtering/index`

---

## What Figure Export captures (quick)

Figure export is designed to be **WYSIWYG**: it snapshots the *current view state* at export time and renders an artifact you can use outside the app.

In typical usage, exports include:

- **The active view** (live view or the focused snapshot)  
  If you enable “Export all views (split-view)”, exports include the full grid (live + snapshots) as a multi-panel figure.
- **Camera + dimension (1D/2D/3D)** so the exported framing matches what you curated.
- **Visibility** (filters and hide/show state): invisible points are not exported.
- **Color-by + legend mapping** sourced from the same legend model used by the viewer.
- **Point appearance** based on the current viewer render state:
  - PNG exports use shader-accurate rendering (so 3D “spheres” match what you see).
  - Hybrid SVG exports rasterize points but keep annotations vector.
- **Annotations you choose**: title, axes, legend, 3D orientation indicator.
- **Overlays tied to the view buffers** such as highlights (rings) and centroid overlays (if enabled in the viewer).

:::{important}
Figure export settings are intentionally **not persisted in session bundles**. If you want to reproduce a specific export later, record the export settings (or rely on embedded metadata; see {doc}`05_metadata_and_provenance`).
:::

---

## What Figure Export does not capture (quick)

Some things are intentionally not exported:

- **Connectivity/edges overlay** (if enabled in the viewer, export currently includes points only).
- **Transient UI**: tooltips, hover labels, toasts, selection lasso UI, etc.
- **Byte-for-byte identical files** across exports: exports include timestamps and embedded metadata, and SVG rendering can vary by font availability.

---

## Output size and DPI (important for PNG)

Cellucid treats **Plot size** as the size of the **plot content area**.

- The exported file can be **larger than your plot size** because legend/axes/title are placed *around* the plot (without shrinking the plot).
- For PNG, **DPI scales the final pixel dimensions**:

`output_pixels ≈ plot_pixels × (DPI / 96)` (plus extra for legend/axes/title)

Example: a 1200×900 plot at 300 DPI will produce an image roughly `1200×(300/96) ≈ 3750` pixels wide (plus annotation padding).

---

## Downloads, filenames, and privacy

- Exports download via your browser (usually to your default Downloads folder).
- Filenames are conservative and include dataset/field/view info, e.g.:
  - `<dataset>_<color-field>_<view>_<timestamp>.svg`
  - `<dataset>_<color-field>_<view>_dpi300_<timestamp>.png` (when exporting multiple DPIs)

:::{important}
Exports can embed provenance metadata that may include **dataset names/ids and source paths/URLs**. If you are sharing figures publicly, skim {doc}`05_metadata_and_provenance` for how to inspect or strip metadata.
:::

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-panel-overview
Suggested filename: figure_export/01_export-panel-overview.png
Where it appears: User Guide → Web App → Figure Export → index.md
Capture:
  - UI location: Figure Export panel open with key sections visible (format/mode, size, preview/export)
  - State prerequisites: dataset loaded; a color-by field set; legend visible
Crop:
  - Include: export panel + enough of the canvas to orient the reader
  - Exclude: private dataset identifiers, file paths, personal avatars
Alt text:
  - Figure Export panel open in the sidebar.
Caption:
  - The Figure Export panel controls output format, size, and rendering mode; preview before exporting large figures.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Figure Export panel overview.
:width: 100%

The Figure Export panel is the entry point for publication-grade SVG/PNG exports.
```

---

## Recommended reading order

1) {doc}`01_figure_export_goals_wysiwyg_and_reproducibility`
2) {doc}`02_export_ui_walkthrough`
3) {doc}`03_export_formats_and_renderers`
4) {doc}`04_quality_knobs_and_best_practices`
5) {doc}`05_metadata_and_provenance`
6) {doc}`06_edge_cases`
7) {doc}`07_troubleshooting_figure_export`
8) {doc}`08_screenshots` (capture checklist)
9) {doc}`09_reference_implementation_notes`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
