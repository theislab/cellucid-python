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

Then the point strategy (**Full vector** until you choose another):

- **Full vector**: maximum editability of points, but only practical for small datasets.
- **Optimized vector**: density-preserving reduction for medium datasets (approximate but still vector).
- **Hybrid**: raster points + vector annotations (recommended for very large datasets and most 3D exports).

:::{tip}
If you are unsure, export **PNG** at 300 DPI first. If you later need editability, export **SVG + PNG** and use the SVG only for text/legend/layout edits.
:::

---

## Size, resolution, and readability

### The two knobs that matter most

Most “this doesn’t look publication-grade” complaints come down to:

1) **Figure width** (millimetres, or the journal column) and **text sizes** (points)
2) **DPI** (for PNG)

Everything else is secondary.

### Figure width vs output pixels (PNG)

The figure is exactly as wide as you set it, at every DPI; a PNG is that
width in inches times the DPI:

`PNG width (px) = width (mm) / 25.4 × DPI` — a double column (183 mm) at 300 DPI is 2163 px.

Practical implications:
- Set the width the figure will be printed at, not a larger one to scale down
  later: its text is set in points (Auto text: 7 pt, the size journals ask
  for) and keeps that size only at the width it was made for.
- A wide figure at 600 DPI can be a *very large* image (slow to export, huge
  file size); export at 150 DPI for drafts.

### Rule-of-thumb sizing recipes

Use these as starting points; then adjust based on your legends and number of panels.

**Slides / lab meeting**
- Format: PNG
- DPI: 150–300
- Width: double column (183 mm), or a custom width for a 16:9 slide (about 250 mm)
- Auto text: ON (or larger text, in points, for the back of the room)

**Manuscript single-panel**
- Format: PNG (300 DPI) or SVG (Hybrid/Full)
- Width: the column the panel goes in (single 89 mm, 1.5 column 120 mm, double 183 mm)
- Auto text sets 7 pt (ticks 6 pt, title 8 pt); never scale the file up or down afterwards, which changes the text's size.

**Multi-panel (2×2 or 3×2 grid)**
- Use **Export all views (split-view)** in grid compare mode.
- Use a double column (183 mm): the panels keep their on-screen arrangement and their views' shapes, and the figure fills the width at the height they make.
- The panels' shape is the views' shape on screen: several views side by side
  are tall and narrow. To choose it, use **Frame export** with **Lock aspect**.
- Size for the legends you will actually get: panels that disagree (different
  field, or the same field with different hidden categories) each carry their
  own legend. A panel legend stands beside its plot when it is no wider than
  the plot, and otherwise below it; legends take room from the plots down to
  half their size without legends, and past that the figure comes out larger
  than the size you typed (nothing is dropped or cut).

### Readability checklist

- Turn **Auto text** on unless you have a specific typography target.
- Use **Legend: Bottom** when a long legend would otherwise make the figure much wider than the plot: below the plot it fills the plot's width and grows downwards instead.
- Prefer **White** or **Match viewer** with high contrast for print.
- Verify in the tool you will actually use (PowerPoint, Illustrator, Inkscape, or a PDF workflow).

---

## Large point clouds: quality without lying

Large datasets are where figure export can become misleading if you’re not careful.

### Pick the right strategy

- SVG has no default strategy and no size-triggered chooser. Select one explicit
  strategy before export.
- Choose according to the artifact you need:
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

Cellucid may block export with an exact requirement:
- **Exact point export unavailable**: PNG or Hybrid lacks WebGL2 or required
  camera matrices.
- **Connectivity overlay not exported**: connectivity is visible, but the
  current figure exporter owns the point layer only.
- **Volumetric smoke cloud render mode not exported** (or **Active render mode
  unavailable**): figure export reproduces the **Points** render mode only.

These blockers never change the selected format, SVG strategy, backend, or
fidelity. Resolve the named condition and submit the same request again.

---

## Color, contrast, and accessibility

### Categorical legends (many categories)

If you have many categories:
- prefer **Legend: Bottom** (more horizontal space),
- consider merging rarely used categories (if your workflow supports it),
- avoid exporting unreadable “wall of legend” figures (it looks unprofessional and is hard to interpret).

:::{note}
Hiding categories does not shorten the exported legend. A category you hid
keeps its entry, marked `(hidden)` with a hollow gray swatch, so a filtered
figure cannot pass as a complete one. A filter on another field does shorten
it: a category left with no cells (the row the on-screen legend shows
disabled) is not listed, and the title (a grid panel's label) names that
filter. Merge or
drop categories upstream if you need a shorter legend.
:::

### Continuous colormaps

- Avoid colormaps that compress contrast (the figure can look “flat” in print).
- Prefer palettes that remain interpretable when printed in grayscale (when possible).
- If **Log color scale** is on for the colored field, the exported colorbar is
  labeled `Log10 color scale` with its true midpoint (the geometric mean of the
  endpoints). Keep that line in the figure: without it a reader interpolates the
  two endpoint labels linearly and reads off values that are wrong by orders of
  magnitude.
- A gene whose zeros are drawn gray (**Gray for zero** on in its legend, or a
  log scale) is exported the same way, and the colorbar says `Gray: 0` (with
  a gene read for a sample of the cells, `Gene read for … cells. Faded: not
  read · gray: 0`); with Gray for zero off its zeros take the scale's lowest
  color, and a gene sample's note says only `Faded: not read`.

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

- **White**: safest for papers and screenshots.
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

## Interface reference

```{figure} ../../../_static/screenshots/figure_export/style.png
:alt: The Style group open: background, font family and the text sizes for the base text, legend, ticks, axis labels, title and centroid labels.
:width: 532px

Style controls expose background, typography, and text sizing.
```



---

## Next steps

- {doc}`03_export_formats_and_renderers` (understand the modes)
- {doc}`05_metadata_and_provenance` (how to record provenance for papers)
- {doc}`07_troubleshooting_figure_export` (SVG too big, exports look wrong)
