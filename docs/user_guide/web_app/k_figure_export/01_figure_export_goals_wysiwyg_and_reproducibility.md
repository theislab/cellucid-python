# Figure export goals (WYSIWYG and reproducibility)

**Audience:** everyone (wet lab + computational + developers writing papers)  
**Time:** 5–10 minutes  
**What you’ll learn:**
- What “WYSIWYG” means for figure export
- What “reproducible export” does (and does not) promise
- How export relates to sessions, dataset identity, and provenance

**Prerequisites:**
- A dataset loaded in the web app

---

## Why figure export exists

Cellucid is built for **interactive exploration**: you rotate/zoom, change color-by fields, apply filters, and curate the exact view you want to communicate.

When you move from exploration to communication, screenshots are rarely enough:
- they are limited by screen resolution,
- they often miss provenance details (what dataset/version? what filters?),
- and they can be painful to edit for papers (labels, legends, layout).

Figure export exists to produce **stable, shareable artifacts**:

- **Publication-ready SVG** for editability (text, legend layout, panel labels).
- **High-DPI PNG** for pixel-perfect appearance (especially for 3D/shader-rendered points).
- **WYSIWYG state capture**: export uses the current view’s camera, visibility, color mapping, and relevant overlays.
- **Provenance**: exports embed structured metadata so you can later answer “what exactly did we export?” (see {doc}`05_metadata_and_provenance`).

:::{note}
Figure export is intentionally designed to be **performance-safe**: it does not add overhead to the normal render loop. Heavy work happens only when you explicitly preview/export.
:::

---

## WYSIWYG vs reproducible (two different promises)

It helps to separate two promises that people often mix:

- **WYSIWYG** (“what you see is what you get”):  
  when you press **Export**, the exported figure should match the state you are currently looking at (camera, visibility, colors, legends, highlights).
- **Reproducible export**:  
  if you (or a collaborator) restore the same dataset/version and the same saved state later, exporting again should produce the *same figure* (within stated limits).

There are also different levels of “reproducible”:

1) **Scientifically reproducible**: clusters/gradients/highlights are the same and support the same interpretation.  
2) **Visually reproducible**: looks indistinguishable to a human viewer.  
3) **Byte-for-byte reproducible**: the exported file bytes are identical. This is rarely realistic because exports include timestamps and embedded metadata.

Common sources of variation include:

- **Browser differences**: SVG and font rendering can differ slightly across browsers.
- **Font availability**: exported SVGs reference font families; if a font is missing on another machine, a substitute font can shift text widths/line breaks.
- **GPU/WebGL differences**: PNG and Hybrid SVG point rendering can depend on WebGL2 availability and GPU/driver behavior.
- **Export strategy choices**:
  - *Optimized vector* intentionally reduces points to preserve density and keep files usable.
  - *Hybrid SVG* rasterizes points (by design) while keeping annotations vector.
  - *Raster strategy* outputs PNG even if you requested SVG.
  These are still valid exports, but they change what “exact match” means.

---

## What figure export should reproduce (checklist)

When you export a figure, you should expect the following to be captured (subject to the options you choose in the export panel):

**View scope**
- The **active view** (live view or the currently focused snapshot view).  
  If you enable **Export all views (split-view)** while in grid mode, the export includes a multi-panel grid of views.

**Geometry**
- Camera framing and projection (rotation/zoom/pan; orbit vs planar mode).
- Active **dimension** (1D/2D/3D).
- Cropping via **Frame export** (if enabled and confirmed).

**Visibility and appearance**
- Visibility after filters and hide/show logic (invisible points are not drawn).
- Color-by field + legend mapping (categorical palettes or continuous colormaps).
- Current point size and (for PNG/Hybrid) shader-accurate point appearance.

**Annotations and overlays**
- Title (if set).
- Axes and axis labels (if enabled).
- Legend (if enabled; position right/bottom).
- 3D orientation indicator (if enabled in 3D orbit mode).
- Highlights and selection emphasis (if enabled; non-selected points can be muted).
- Centroid overlay (if centroids are enabled in the viewer).

:::{tip}
If you need a “methods-ready” record of what was exported, use {doc}`05_metadata_and_provenance`. PNGs and SVGs embed structured metadata including dataset identity, field key, filter summary, and export settings.
:::

---

## What figure export does NOT reproduce (must be explicit)

Figure export is intentionally scoped. Common “expected but missing” items include:

- **Connectivity/edges overlay**: even if connectivity lines are enabled in the viewer, export currently exports the **point layer only** (no edges).
- **Transient UI**: hover tooltips, lasso outlines, selection handles, notifications/toasts, menus.
- **Exact point-for-point identity when using reduction**:
  - Optimized vector uses density-preserving reduction; the exported point set is *not* identical to the on-screen set.
  - Hybrid SVG rasterizes points; individual points are not editable as separate SVG circles.
- **Guaranteed byte-for-byte identical files**: exports include timestamps and embedded provenance; re-exporting later will produce a different filename and metadata timestamps even if the figure looks identical.
- **Downstream tool quirks**: some SVG editors have limitations (very large SVGs can be slow or crash; embedded images in Hybrid SVG can be handled differently across tools).

---

## Practical reproducibility recipe (paper + collaboration)

Use this workflow when you want to be able to recreate a figure weeks later (or hand it to a collaborator):

1) **Freeze the data identity**
   - Record the dataset identifier/version you loaded (export folder version tag, dataset id, or hosting URL).
   - If you are iterating on exported data, decide which version is “figure-final”.

2) **Freeze the app state**
   - Save a **session bundle** that captures the view (camera, filters, active field, highlights, etc.).  
     See: {doc}`../l_sessions_sharing/index`

3) **Record export settings**
   - Plot size (W×H)
   - Format (SVG/PNG) and DPI (for PNG)
   - Large dataset strategy (full/optimized/hybrid/raster)
   - Legend/axes/title/background choices
   - Frame export crop settings (if used)

   :::{important}
   Figure export UI state is not persisted in sessions. If you want to re-export later with the same settings, either record them manually or rely on the embedded metadata in the exported file.
   :::

4) **Export and verify**
   - Export the figure.
   - Open it **outside Cellucid** (browser for SVG; image viewer for PNG; Illustrator/Inkscape if you plan to edit).
   - Confirm that what you intended is present: correct field/legend, correct view, correct filters/visibility.

5) **Package for collaboration**
   - Share: dataset export (or dataset URL) + session bundle + exported figure(s).
   - If you post-process in Illustrator, keep both:
     - the original exported artifact (SVG/PNG), and
     - the edited final figure (AI/SVG/PDF), so provenance isn’t lost.

Related pages:
- {doc}`../l_sessions_sharing/index`
- {doc}`05_metadata_and_provenance`

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-wysiwyg-vs-output
Suggested filename: figure_export/02_wysiwyg-vs-exported-output.png
Where it appears: User Guide → Web App → Figure Export → 01_figure_export_goals_wysiwyg_and_reproducibility.md
Capture:
  - Show the viewer state and the exported artifact opened outside the app (side-by-side or two windows)
  - Prefer a state with a legend + title + at least one overlay (highlight, vector field, etc.) if supported
Redact:
  - Remove: private dataset ids/names and local file paths
Alt text:
  - Side-by-side view of the Cellucid viewer and the exported figure.
Caption:
  - Compare the viewer state to the exported artifact to validate what “WYSIWYG” and “reproducible export” mean in practice.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for WYSIWYG vs exported output comparison.
:width: 100%

Compare the in-app view to the exported artifact to validate which parts of state are reproduced.
```

---

## Next steps

- {doc}`02_export_ui_walkthrough` (click-by-click export)
- {doc}`03_export_formats_and_renderers` (choose SVG/PNG mode)
- {doc}`07_troubleshooting_figure_export` (if exports look different)
