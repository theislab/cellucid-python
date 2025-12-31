# Export UI walkthrough

**Audience:** wet lab + computational users (figures for slides, reports, and papers)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Where the Figure Export panel lives and how to open it
- The difference between preview and export
- How to set size/aspect ratio and use the framing overlay
- How export behaves for live views vs snapshots (and multi-view layouts, if supported)

**Prerequisites:**
- A dataset loaded
- (Recommended) A defined color-by field so you can sanity-check legends

---

## Fast path (export a decent figure in <2 minutes)

This is the “no excuses” workflow that works for most people.

1) **Get the view into the exact state you want to publish**
   - Set the embedding dimension (2D vs 3D).
   - Choose your **color-by** field and confirm the legend looks right.
   - Apply filters (if any) and confirm you didn’t accidentally hide important cells.

2) **Make sure the correct view is active**
   - If you have snapshot views: click inside the panel you want to export (in grid compare) or switch to “edit selected view” to make it obvious.

3) **Open** the **Figure Export** panel (left sidebar).

4) **Choose plot size**
   - Start with **Screen (half)** or **1200 × 900**.
   - If you want a larger, “final” export, increase plot size (but note PNG + high DPI can get huge; see below).

5) **Choose download format**
   - **PNG (300 DPI)** is a strong default for slides and manuscripts.
   - **SVG** is the right choice if you want to edit labels/legend layout in Illustrator/Inkscape.

6) **(Optional) Crop without changing camera**
   - Enable **Show preview** → enable **Frame export** → drag the frame → click **Confirm** to apply.

7) Click **Export**.

8) **Open the file outside Cellucid** and do a 10-second sanity check:
   - correct view (live vs snapshot),
   - correct legend/colors,
   - correct filters/visibility,
   - readable text at your intended size.

---

## Where the Figure Export panel lives

Figure export is controlled from the **Figure Export** panel in the left sidebar.

**Which view gets exported?**

- By default, Cellucid exports the **active view**:
  - the live view, or
  - the snapshot view you currently have focused/selected.
- In multiview **grid compare**, click inside the panel you want to export first (focus matters).
- If you enable **Export all views (split-view)** and you are in grid compare mode, Cellucid exports a **multi-panel grid** (live + snapshots).

**Are export settings saved?**

- Figure export settings are **not stored in session bundles** (by design).
- Treat export settings as “ephemeral UI”: if you care about reproducing a figure later, record the settings or rely on embedded metadata in the exported artifact (see {doc}`05_metadata_and_provenance`).

---

## Preview vs export

**Preview** is a lightweight, export-style approximation:

- It is meant to help you confirm **layout and framing** (especially cropping).
- It is intentionally **downsampled and debounced** so the UI stays responsive on large datasets.
- It supports **colorblindness simulation** (Normal / Deuteranopia / Protanopia / Tritanopia) to help you catch problematic palettes.

**Export** is the real thing:

- It renders from the current view buffers and camera matrices.
- It applies your chosen format and strategy (SVG full/optimized/hybrid or PNG).
- It can show dialogs/warnings:
  - a **Large Dataset Export** chooser when exporting large SVGs with “Ask” strategy,
  - **Export Fidelity Warnings** when something cannot match the on-screen view exactly (e.g., 3D shader-accurate points cannot be represented as pure SVG circles, or WebGL2 is unavailable).

:::{tip}
If your preview looks wrong, click **Refresh** in the preview controls (especially after changing views, fields, or filters).
:::

---

## Setting size and aspect ratio

Cellucid exposes **Plot size** controls as width × height in pixels. In the UI you’ll see:

- A size preset dropdown:
  - **Screen (half)** *(default)*
  - **Screen (current)**
  - **1200 × 900**
  - **1600 × 1200**
  - **1920 × 1080**
- Editable numeric inputs for **W** and **H**.

### Plot size vs exported file size

The numbers you enter represent the desired **plot content size**, not the entire exported figure.

- If you include axes/legend/title, the final file will be **larger** because those elements are laid out *around* the plot (the plot is not shrunk to make room).

### PNG DPI scaling (the part that surprises people)

For PNG export, Cellucid scales the final pixel dimensions by DPI:

`output_pixels ≈ plot_pixels × (DPI / 96)` (plus annotation padding)

Example:
- Plot width 1200 at 300 DPI → roughly `1200 × (300/96) ≈ 3750` pixels wide for the plot area.

Practical advice:
- If PNG exports are “too huge”, either:
  - reduce plot size, or
  - export at 150 DPI for drafts/slides.

### Auto text sizing

By default, **Auto text** is enabled:
- text sizes are automatically scaled to the plot size,
- and you can still override by disabling Auto text and setting sizes manually.

This is usually what you want for multi-panel figures: pick a plot size first, then let the UI scale labels/legend appropriately.

---

## Titles, captions, legends, axes

### Title

- The **Title** field auto-fills based on your dataset name and active field/visible categories (when available).
- If you type your own title, auto-fill turns off.
- To return to auto-fill behavior, clear the title field completely.

### Legend

- **Legend** is enabled by default.
- Position can be:
  - **Legend: Right**
  - **Legend: Bottom**

Legend contents are sourced from the same legend model as the viewer, so categorical colors and continuous colormaps should match what you see.

### Axes

- **Axes** is enabled by default.
- Axes ticks are computed from the **currently visible region** (after filtering and cropping).
- Axis labels default to **X** and **Y**.
  - Clear X and/or Y to hide axis labels entirely.

How to interpret axes:
- In 2D/planar views, axes aim to correspond to embedding coordinates (when available).
- In 3D orbit views, axes show **camera-space** coordinates (useful for orientation, but not the same as “UMAP_1/UMAP_2” in a strict sense). If you want an explicit 3D orientation cue, keep **3D orientation** enabled.

### 3D orientation + depth sort

- **3D orientation** adds a small orientation widget in 3D orbit mode.
- **Depth sort** is enabled by default and helps 3D point clouds export with the correct near/far ordering.
  - Turning it off can be faster, but may reduce fidelity in dense 3D views.

### Style controls that affect export readability

In the **Style** section you can choose:
- Background: **Match viewer**, **White BG**, **Transparent**, **Custom…**
- Font family: Arial / Helvetica / Times (safe cross-platform defaults)
- Text sizing:
  - Auto text on/off
  - Base / Legend / Ticks / Axis / Title / Centroids font sizes

:::{tip}
If you plan to edit the SVG in Illustrator/Inkscape, choose a common font (Arial/Helvetica) to minimize font-substitution surprises.
:::

---

## Selection emphasis (optional)

If you have an active selection/highlight and want the exported figure to communicate it clearly, enable:

- **Emphasize selection**

This makes non-selected points:
- turn gray, and
- drop in opacity (controlled by the **α** input; default ~0.15).

Practical use cases:
- A “selected cluster vs background” figure for a talk.
- Highlighting a rare population without hiding the rest of the dataset (hiding can be misleading).

:::{important}
Emphasize selection does nothing if there are no highlighted/selected cells.
:::

---

## Download options and large dataset strategies

### Download format options

In the **Download** section you can choose:

- **SVG**: one SVG export.
- **PNG**: one PNG export at the selected DPI.
- **SVG + PNG**: exports both formats (useful when you want an editable SVG *and* a WYSIWYG raster).
- **PNG (150/300/600)**: exports three PNGs at common DPIs.
- **All**: exports SVG plus PNGs at 150/300/600.

### DPI (for PNG)

When exporting a single PNG, you can choose:
- 150 / 300 / 600 DPI

Higher DPI means:
- more pixels (larger files, more time),
- but better raster fidelity for printing and zooming.

### Large dataset strategy (for SVG exports)

SVG export can become impractical when you have tens/hundreds of thousands of visible points.
Cellucid gives you control via **Large dataset strategy**:

- **Large data: Ask** *(recommended default)*  
  If visible points exceed the threshold (default ~50k), Cellucid prompts you to choose.
- **Full vector**  
  Exports every point as an SVG circle (maximum editability; can be huge/slow).
- **Optimized vector**  
  Uses density-preserving reduction to keep approximately the **Keep** point count (default ~100k).  
  Good for medium datasets when you need a fully-vector plot but can tolerate approximation.
- **Hybrid** *(recommended for 3D and very large datasets)*  
  Rasterizes points (WYSIWYG, including 3D shading) but keeps annotations vector.
- **Raster (PNG)**  
  Exports PNG even if you requested SVG (maximum compatibility; no editable SVG points).

:::{note}
If the current view uses shader-accurate 3D point rendering, pure SVG circles cannot reproduce it. Cellucid may automatically switch SVG exports to **Hybrid** and show an export fidelity warning.
:::

---

## Framing overlay (crop) — mental model

**Frame export** is a “photography-style crop” tool:

- You keep the same camera framing in the viewer,
- but choose a smaller sub-region to export.

This is useful when:
- the camera framing you like includes extra whitespace,
- you want to crop to a specific aspect ratio for a figure panel,
- you want to “zoom” for export without disturbing the interactive view.

### How to use Frame export (step-by-step)

1) Enable **Show preview** (preview is required for framing).
2) Enable **Frame export**.
3) Drag the frame in the preview:
   - drag inside the frame to move it,
   - drag edges/corners to resize.
4) (Optional) Click **Lock aspect** to keep the crop aspect synced to the plot size.
5) Click **Confirm** to apply the crop (the preview will switch to a framed preview).
6) Export.

Useful actions:
- **Reset** resets the crop to the full frame.
- **Fit plot** adjusts your plot height to match the current crop aspect.
- Double-clicking the preview resets the crop.

:::{important}
**Confirm** toggles the preview between:
- “edit framing” (full preview with a draggable frame), and
- “framed preview” (shows what the crop will export).

As long as **Frame export** is enabled, the current frame is used for export; use **Confirm** to validate the crop before you export.
:::

---

## Live view vs snapshot export

### Exporting a single view (default)

By default, export produces a **single-panel figure** from the **active view**:

- **Live view**: export the current working view.
- **Snapshot view**: export the currently focused snapshot.

In multiview workflows:
- In **Grid compare**, click inside the view panel you intend to export first.
- In **Edit selected view**, the selected view is unambiguous.

### Exporting all views (split-view)

If you have snapshots and want a multi-panel figure:

1) Switch to multiview **Grid compare** layout.
2) Enable **Export all views (split-view)**.
3) Export.

What to expect:
- The export is arranged as an automatically computed grid.
- Panels are labeled (e.g., A, B, C…) and include per-panel view labels.
- A legend is shared only when all panels use the same active color-by field; otherwise legend behavior may be limited (export individual panels if you need separate legends/fields).

:::{tip}
If you want all panels to be directly comparable, lock cameras across views (see {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`).
:::

---

## Screenshot placeholders (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-preview-vs-export-controls
Suggested filename: figure_export/05_preview-vs-export-controls.png
Where it appears: User Guide → Web App → Figure Export → 02_export_ui_walkthrough.md
Capture:
  - UI location: Figure Export panel with both Preview and Export controls visible
  - State prerequisites: dataset loaded; non-empty canvas; legend visible
Crop:
  - Include: preview/export buttons + relevant size/mode controls
Alt text:
  - Figure Export panel showing preview and export controls.
Caption:
  - Preview is for sanity-checking crop and layout; Export produces the final artifact and may show warnings for large figures.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for preview vs export controls.
:width: 100%

Use preview to confirm crop/layout, then export the final SVG/PNG artifact.
```

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-framing-overlay
Suggested filename: figure_export/06_framing-overlay.png
Where it appears: User Guide → Web App → Figure Export → 02_export_ui_walkthrough.md
Capture:
  - Enable the framing/crop overlay and show resize handles on top of the canvas
  - Prefer a state where the crop visibly excludes some content (so the effect is obvious)
Alt text:
  - Framing overlay used to crop the exported region.
Caption:
  - Use the framing overlay to crop the exported region without changing the camera view.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the framing overlay.
:width: 100%

The framing overlay lets you control the exported crop independently of camera navigation.
```

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-exported-artifact-example
Suggested filename: figure_export/07_exported-artifact-opened.png
Where it appears: User Guide → Web App → Figure Export → 02_export_ui_walkthrough.md
Capture:
  - Export one SVG and one PNG, then open at least one outside the app (Preview, browser, Inkscape, Illustrator)
Redact:
  - Remove: private dataset ids/names and local file paths
Alt text:
  - An exported figure opened outside the web app.
Caption:
  - Always open the exported SVG/PNG outside Cellucid to confirm it is shareable and renders correctly in downstream tools.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for an exported artifact opened outside the app.
:width: 100%

Open the exported file outside the app to confirm it renders correctly in the tools you’ll actually use.
```

---

## Next steps

- {doc}`03_export_formats_and_renderers` (what each SVG/PNG mode means)
- {doc}`04_quality_knobs_and_best_practices` (publication-grade settings)
- {doc}`07_troubleshooting_figure_export` (export button does nothing, huge SVG, etc.)
