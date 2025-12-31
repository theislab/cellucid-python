# Quality knobs and best practices

**Audience:** wet lab + computational users preparing “final” figures  
**Time:** 15–30 minutes  
**What you’ll learn:**
- Practical rules-of-thumb for choosing SVG vs PNG
- How to get readable text/legends at the target figure size
- How to handle huge point clouds without misleading downsampling
- Common pitfalls that make “publication figures” look unprofessional

**Prerequisites:**
- A dataset loaded
- You can already export a basic figure (see {doc}`02_export_ui_walkthrough`)

---

## Decide first: vector vs raster (what matters for your use case)

Start by deciding what matters most **after** export:

- Choose **SVG** when you want to edit:
  - text labels and title,
  - legend layout,
  - panel arrangement (multi-panel figure assembly).

- Choose **PNG** when you want:
  - maximum compatibility (everyone can open it),
  - pixel-identical appearance (especially for 3D/shader-accurate point rendering),
  - minimal surprises across tools.

Then choose the point strategy:

- **Full vector**: maximum editability of points, but only practical for small datasets.
- **Optimized vector**: density-preserving reduction for medium datasets (approximate but still vector).
- **Hybrid**: raster points + vector annotations (recommended for very large datasets and most 3D exports).

:::{tip}
If you are unsure, export **PNG (300 DPI)** first. If you later need editability, export **SVG + PNG** and use the SVG only for text/legend/layout edits.
:::

---

## Size, resolution, and readability

### The two knobs that matter most

Most “this doesn’t look publication-grade” complaints come down to:

1) **Plot size** (W×H)
2) **DPI** (for PNG)

Everything else is secondary.

### Plot size vs output pixels (PNG)

Cellucid scales PNG pixel dimensions by DPI:

`output_pixels ≈ plot_pixels × (DPI / 96)` (plus extra for legend/axes/title)

Practical implications:
- A large plot size combined with 300–600 DPI can create *very large* images (slow to export, huge file size).
- If you want a smaller output, reduce plot size or export at 150 DPI for drafts.

### Rule-of-thumb sizing recipes

Use these as starting points; then adjust based on your legends and number of panels.

**Slides / lab meeting**
- Format: PNG
- DPI: 150–300
- Plot size: 1200×900 or 1600×1200
- Auto text: ON

**Manuscript single-panel**
- Format: PNG (300 DPI) or SVG (Hybrid/Full)
- Start with plot size: 900×700 to 1400×1000
- If text is too small, increase plot size rather than scaling up in PowerPoint/Illustrator.

**Multi-panel (2×2 or 3×2 grid)**
- Use **Export all views (split-view)** in grid compare mode.
- Increase overall plot size so each panel has enough pixels (e.g., 2000×1500+ as a starting point).

### Readability checklist

- Turn **Auto text** on unless you have a specific typography target.
- Use **Legend: Bottom** when the legend is wide and would otherwise squeeze the plot.
- Prefer **White BG** or **Match viewer** with high contrast for print.
- Verify in the tool you will actually use (PowerPoint, Illustrator, Inkscape, or a PDF workflow).

---

## Large point clouds: quality without lying

Large datasets are where figure export can become misleading if you’re not careful.

### Pick the right strategy

- If you see a **Large Dataset Export** dialog, it means you have many visible points and full vector SVG may be impractical.
- Recommended choices:
  - **Optimized vector** for medium datasets when you still want a fully vector plot.
  - **Hybrid** for very large datasets or when exporting 3D views.
  - **PNG** when you need maximum compatibility or the SVG is too heavy.

### Be explicit when you downsample

If you export with **Optimized vector**, you are intentionally reducing points.
That can be scientifically legitimate (density preservation), but it should be disclosed when relevant.

Suggested caption language (adapt as needed):
- “Points were exported with density-preserving reduction for visualization.”

### Protect rare populations

If a rare population matters:
- verify it is still visible in the optimized export,
- consider exporting a second panel focused on that population (Frame export crop + emphasize selection),
- or export a PNG/hybrid which preserves the full point set appearance more faithfully.

### Interpret export warnings correctly

Cellucid may show **Export Fidelity Warnings** such as:
- “Point reduction enabled” (optimized vector is approximate)
- “SVG strategy adjusted for WYSIWYG” (3D/shader point appearance requires Hybrid)
- “Shader-accurate export may degrade” (WebGL2/matrices missing; fallback may look flatter)

These warnings are not necessarily “errors”; they are telling you which promise changed (pixel-perfect vs density-preserving vs rasterized points).

---

## Color, contrast, and accessibility

### Categorical legends (many categories)

If you have many categories:
- prefer **Legend: Bottom** (more horizontal space),
- consider hiding/merging rarely used categories (if your workflow supports it),
- avoid exporting unreadable “wall of legend” figures (it looks unprofessional and is hard to interpret).

### Continuous colormaps

- Avoid colormaps that compress contrast (the figure can look “flat” in print).
- Prefer palettes that remain interpretable when printed in grayscale (when possible).

### Colorblindness preview (high value, low effort)

In the export preview, use the colorblind simulation dropdown:
- Deuteranopia
- Protanopia
- Tritanopia

This is a fast way to catch:
- categories that collapse to the same perceived color,
- low-contrast colormaps that become ambiguous.

:::{note}
Colorblindness simulation is **preview-only**. It does not change the exported artifact.
:::

### Background choices

- **White BG**: safest for papers and screenshots.
- **Transparent**: useful for overlaying in slides, but be careful:
  - some viewers show transparency as a checkerboard,
  - anti-aliased text/halo effects can look different depending on what you place behind it.

---

## Post-processing (Illustrator/Inkscape) — safe edits vs dangerous edits

### Safe edits (usually fine)

- Edit title text, axis labels, legend title, and annotation text.
- Reposition legend (without changing the plot scaling).
- Add panel labels (A, B, C…) or align panels.
- Add arrows/callouts that point to regions (without altering the data layer).

### Dangerous edits (can mislead)

- Moving points, deleting points, or manually recoloring subsets outside the documented color mapping.
- Non-uniformly scaling the plot region (changes relative distances).
- Cropping in a way that hides context without disclosure.
- Changing legend mapping after export (e.g., swapping category colors) without documenting.

### Avoid accidental rasterization

In vector editors it’s easy to accidentally rasterize everything (especially when applying effects).
Best practices:

- Keep a copy of the original exported SVG/PNG untouched.
- If you need complex effects, apply them to a duplicate layer or in a raster workflow, and keep the original as provenance.
- If using Hybrid SVG, remember: points are already raster; treat the SVG as “vector annotations + embedded image”.

---

## Screenshot placeholders (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-quality-knobs
Suggested filename: figure_export/09_quality-knobs.png
Where it appears: User Guide → Web App → Figure Export → 04_quality_knobs_and_best_practices.md
Capture:
  - UI location: quality-related controls (mode, size, legend/text options, any downsampling toggles)
Alt text:
  - Figure Export quality settings controls.
Caption:
  - Adjust size and mode first; then tune legend/text options so the exported figure is readable at its final intended scale.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for figure export quality controls.
:width: 100%

Most “publication quality” issues come from size and readability choices, not from minute rendering tweaks.
```

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-large-export-warning
Suggested filename: figure_export/10_large-export-warning.png
Where it appears: User Guide → Web App → Figure Export → 04_quality_knobs_and_best_practices.md
Capture:
  - Trigger any warning dialog shown for very large SVG exports (or a banner warning)
  - Include the recommended action (switch mode, reduce size, etc.) if shown
Alt text:
  - Warning dialog about a large export.
Caption:
  - Large SVG exports can be impractical to open/edit; switch to optimized or hybrid modes when prompted.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for large export warning.
:width: 100%

For huge datasets, expect warnings and use optimized/hybrid modes to keep outputs usable in downstream tools.
```

---

## Next steps

- {doc}`03_export_formats_and_renderers` (understand the modes)
- {doc}`05_metadata_and_provenance` (how to record provenance for papers)
- {doc}`07_troubleshooting_figure_export` (SVG too big, exports look wrong)
