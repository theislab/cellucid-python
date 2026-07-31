# Export formats and renderers

**Audience:** computational users + power users (format tradeoffs, editability, file size)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What Cellucid’s export formats are (SVG vs PNG)
- What the SVG renderers mean (full / optimized / hybrid)
- How to choose a mode based on dataset size and downstream editing
- Common incompatibilities with external tools (Illustrator/Inkscape, browsers)

**Prerequisites:**
- A dataset loaded
- Basic familiarity with “vector vs raster” is helpful but not required

---

## Overview: what you can export

Cellucid’s figure export supports:

**Formats**
- **SVG** (`.svg`) — vector container (best for editing text/legend/layout).
- **PNG** (`.png`) — high-DPI raster (best for pixel-fidelity and maximum compatibility).
- Optional batch exports:
  - **SVG + PNG**
  - **PNG (150/300/600)**
  - **All** (SVG + PNGs at 150/300/600)

Each batch is delivered as one deterministic, uncompressed ZIP. Its entries
retain the exact native SVG/PNG bytes and filenames produced by the renderers.

**What gets exported**
- By default: the **active view** (live view or focused snapshot).
- Optionally: **all views** in split-view grid mode (multi-panel export).

**Where files go**
- A single export downloads as `.svg` or `.png`; a batch downloads once as
  `.zip`.
- Filenames include dataset/field/view and a timestamp (and DPI tag when relevant). See {doc}`index` for a quick summary.

---

## Decision table (choose the right mode)

Use this table as a first pass; then adjust based on your dataset size and what you plan to do after export.

| Your goal | Recommended export | Why |
|---|---|---|
| Put a figure in slides quickly | PNG (300 DPI) | Good defaults, “looks like the viewer”, widely compatible |
| Submit a figure to a paper and do minor edits | SVG (Hybrid or Full) + optional PNG | SVG keeps text/legend editable; hybrid keeps point appearance sane for large/3D |
| Edit the layout/labels extensively | SVG (Full or Hybrid) | Vector annotations are easy to tweak; hybrid avoids gigantic point-count SVGs |
| You have >50k–200k visible points | SVG (Optimized) or Hybrid | Full vector SVG may be huge/slow; optimized/hybrid are practical |
| You have >200k visible points or 3D shading | Hybrid SVG or PNG | Pure SVG circles cannot reproduce shader-accurate 3D point appearance |
| Illustrator/Inkscape crashes on the SVG | Hybrid SVG or PNG | Points become rasterized; annotations stay readable; file becomes manageable |
| You need maximum compatibility, no editability | PNG (150/300/600) | Raster always opens; choose DPI based on target use |

---

## SVG export (vector)

SVG export is a vector-first path designed for editability.

### What is vector in Cellucid SVG exports?

- **Always vector** (when enabled):
  - title text
  - axes and tick labels
  - legends (categorical swatches / continuous colorbars) — one shared legend
    for a multi-panel figure whose panels match, otherwise one legend per panel
  - panel labels and the `n = N selected` badge
  - 3D orientation indicator
  - centroid labels (when enabled in the viewer)
- **Points depend on strategy**:
  - **Full vector**: points are SVG circles
  - **Optimized vector**: fewer SVG circles (density-preserving reduction)
  - **Hybrid**: points are rasterized into an embedded PNG `<image>` (annotations remain vector)

### Fonts in SVG

SVG exports reference a font family (Arial/Helvetica/Times). They do not pre-convert text to paths.

- If the same font is available in your editor (Illustrator/Inkscape), text widths should look consistent.
- If a font is missing, editors will substitute, which can shift text and legends.

### Why SVGs get huge (and how to avoid it)

SVG size scales roughly with:
- number of visible points (each point can become a `<circle>`),
- legend complexity (number of categories),
- and number of panels (multiview export).

If you hit size/performance limits, switch to **Optimized** or **Hybrid** (or export PNG).

### Full vector SVG

**What it is:** Every visible point is exported as an SVG `<circle>`.

**Best for:**
- small datasets (roughly <50k visible points),
- maximum editability (you can select/edit individual points in a vector editor),
- workflows where downstream tooling expects fully vector points.

**Tradeoffs:**
- files can become enormous and slow to open for large point clouds,
- editors may crash or become unusable.

### Optimized vector SVG

**What it is:** Cellucid applies density-preserving reduction in viewport space to keep approximately the configured number of points (default ~100k).

**What is preserved:**
- cluster shapes and density patterns (the goal is “same visual story”),
- overall color distribution.

**What is approximated:**
- the exact set of points (it is not pixel-for-pixel identical to the full view),
- fine-scale “which exact cells are drawn” at the point level.

**How to sanity-check:**
- export once in Full vector (for a small subset or smaller size) and compare,
- confirm that rare populations are still visible if they are important,
- if you need exact point appearance (especially in 3D), use Hybrid or PNG instead.

### Hybrid SVG

**What it is:** Points are rasterized into an embedded PNG image (inside the SVG), while annotations remain vector.

This is often the best default when:
- you have a large dataset,
- you want to keep legends/axes/text crisp and editable,
- and you need the exported points to match the viewer’s appearance (including 3D shading).

**What happens when you zoom in later?**
- text/legend/axes stay crisp (vector),
- points are raster (they will eventually pixelate if you zoom extremely far in).

:::{note}
Hybrid SVG requires WebGL2 and the current camera matrices for shader-accurate
points. If either is unavailable, export is blocked until the capability is
restored.
:::

---

## PNG export (raster)

PNG export produces a high-DPI raster image that is easy to paste into slides, manuscripts, and reports.

**When PNG is the right choice**
- You want “what I saw” fidelity with minimal surprises.
- You are exporting 3D views with shader-accurate “sphere” point rendering.
- You need maximum compatibility (everyone can open PNG).

**Common pitfalls (and fixes)**
- *“Text looks blurry”*: export at a larger plot size and/or higher DPI; don’t export tiny and scale up later.
- *“Files are enormous”*: reduce plot size and/or export at 150 DPI for drafts/slides.

**Provenance**
- PNG exports embed UTF-8 `iTXt` metadata, including a structured JSON blob with dataset identity, export settings, and one record per exported panel (view, field, filters). See {doc}`05_metadata_and_provenance`.

---

## Browser contract and downstream tool notes

### Browsers

- Chromium, Firefox, and Safari/WebKit use the same native single-file and ZIP
  batch contracts.
- PNG and Hybrid require WebGL2 plus the published camera matrices. If either is
  absent, **Export blocked** names the missing requirement and no file is
  created.
- Every format additionally requires the **Points** render mode; a volumetric
  smoke view is blocked rather than exported as a point cloud.
- Diagnose unexpected output from the embedded metadata, exact export settings,
  and any `[FigureExport]` console error. Changing browsers is not an export
  strategy.

### Illustrator/Inkscape tips

- If an SVG is slow/crashes your editor, open it in a web browser first to confirm it’s valid.
- For huge datasets, prefer **Hybrid** or **PNG** to avoid editor crashes.
- Hybrid SVGs contain an embedded raster image for points; most editors handle this fine, but operations like “select individual points” are not possible (by design).

### “Open SVG in a browser” sanity check

Yes: it’s a fast way to confirm the export is not corrupted and that legends/axes/title are present before you send the file to collaborators or open it in Illustrator.

---

## Interface reference

```{figure} ../../../_static/screenshots/figure_export/download.png
:alt: Figure Export Download controls with output format, DPI, and split-view choice.
:width: 224px

Download settings choose output format, resolution, and whether all views are exported.
```

---

## Next steps

- {doc}`04_quality_knobs_and_best_practices` (how to get the best output)
- {doc}`06_edge_cases` (huge legends, 0 visible points, smoke/fog quirks)
- {doc}`07_troubleshooting_figure_export` (SVG too big, PNG looks different)
