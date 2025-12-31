# Edge cases

**Audience:** computational users + power users  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Edge cases that commonly break or mislead figure exports
- How to recognize “this export is not trustworthy”
- Practical workarounds for pathological states (0 visible points, huge legends, smoke/fog)

**Prerequisites:**
- A dataset loaded

---

## Data/state edge cases

### Exporting with 0 visible points (all cells filtered out)

What you may see:
- an empty plot region (no points),
- axes/legend/title may still render,
- the exported artifact is “valid” but scientifically meaningless unless intentional.

How to handle:
- If you expected points: check filters and visibility first (see {doc}`../e_filtering/index`).
- If this is intentional (e.g., demonstrating a filter failure mode): consider adding an explicit caption note.

### Exporting tiny groups (rare populations, single-cell highlights)

Common failure mode: the “important” points disappear in an optimized export or become visually insignificant.

Best practices:
- Use **Emphasize selection** (mute non-selected points) rather than relying on the rare points being obvious.
- Avoid **Optimized vector** if the exact presence of rare points matters; prefer **Hybrid** or **PNG**.
- Increase plot size so single points are visible and anti-aliasing is not lost.

### NaN/Inf and invalid values

If positions/colors contain NaN/Inf:
- some points can be skipped during rendering,
- axis bounds can be computed incorrectly (or fall back to a default range),
- downstream tools may behave unpredictably for SVG.

If you see odd axes or missing regions:
- verify your exported dataset does not contain NaN/Inf in embeddings,
- try exporting a smaller subset or a different view to isolate the issue.

### Missing/renamed/deleted fields

Exports rely on the **active color-by field** and its legend model.
If the field was renamed/removed (or differs across snapshot views):
- the legend can be missing or not shared in multi-panel exports,
- your export may not match what you thought you were exporting.

Workarounds:
- confirm the active field immediately before exporting,
- for multiview exports, keep the same active field across panels when you need a shared legend,
- export panels individually when each panel uses a different field.

---

## Scale edge cases (big dataset, huge legends)

### Millions of points (SVG size explosion)

If you have hundreds of thousands to millions of visible points:
- **Full vector SVG** can become enormous and slow or impossible to open.

Recommended approach:
- Use **Hybrid SVG** or **PNG**.
- If you need a fully vector file, use **Optimized vector** and explicitly disclose point reduction when appropriate.

### Category explosion (legend overflow)

Categorical legends with hundreds/thousands of categories often produce:
- unreadable figures,
- huge file sizes,
- legends that dominate the layout.

Strategies:
- Switch legend position to **Bottom** (more space).
- Disable legend for the figure and mention in the caption how to interpret colors.
- Reduce categories upstream (merge, filter, or focus on a subset) if scientifically appropriate.

### Many views/snapshots (multi-panel export too large/slow)

Multi-panel export can become huge when you combine:
- many panels,
- large plot sizes,
- and high DPI PNG export.

Strategies:
- Export fewer panels at a time.
- Reduce DPI for drafts (150 DPI), then re-export final at 300 DPI.
- Increase overall plot size only as needed so each panel remains readable.

---

## Rendering edge cases (smoke/fog/overlays)

### 3D shader-accurate points and SVG limitations

Pure SVG circles cannot represent certain “shader-accurate” point appearances (e.g., 3D sphere shading).
Cellucid may:
- automatically switch SVG export to **Hybrid** (raster points + vector annotations), and/or
- warn that fidelity may degrade if WebGL2 is unavailable.

If you need the exported points to match the on-screen appearance:
- export **PNG** or **Hybrid SVG**,
- export from a browser/environment with WebGL2 available.

### Depth ordering in dense 3D views

If points look “inside out” or layering seems wrong:
- ensure **Depth sort** is enabled for export.

### Axes interpretation in 3D

In 3D orbit mode, axes are reported in **camera-space** coordinates:
- useful for orientation,
- not the same as “UMAP_1/UMAP_2” in the strict embedding coordinate sense.

If you need axes that correspond to embedding coordinates, export a 2D/planar view.

### Connectivity overlay not exported

Even if connectivity lines are visible in the viewer, exports currently include:
- the **point layer** (and related point-based overlays),
- **not** the connectivity edges.

If edges are essential for your figure, capture them via another workflow and document how they were generated.

---

## Environment edge cases

### WebGL2 / GPU differences

PNG and Hybrid SVG exports try to use WebGL2 for shader-accurate point rendering.
If WebGL2 is unavailable or restricted (common in locked-down environments):
- exports can fall back to simpler rendering (flatter dots),
- you may see fidelity warnings.

### Fonts and text layout differences

SVG exports reference font families. If the font is missing on another machine:
- the editor substitutes a font,
- text widths change and legends may shift.

Mitigations:
- use common fonts (Arial/Helvetica),
- or convert text to outlines during final figure production (after you are confident in the layout).

### Browser download restrictions

If clicking Export does nothing:
- your browser may be blocking downloads (popup/download settings),
- a corporate browser policy may restrict file creation,
- extensions can interfere.

Try:
- an incognito window,
- a different browser,
- exporting a small PNG first as a “download sanity check”.

### Privacy/provenance surprises

Exports can include:
- dataset names/ids in filenames,
- embedded metadata that may include source URLs or local dataset paths.

If you plan to share figures publicly, review {doc}`05_metadata_and_provenance` for inspection/stripping workflows.

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-zero-visible-points
Suggested filename: figure_export/12_zero-visible-points.png
Where it appears: User Guide → Web App → Figure Export → 06_edge_cases.md
Capture:
  - Apply filters so no points are visible, then attempt export (or show export preview)
  - Capture any warning/empty-state messaging the export UI provides
Alt text:
  - Export state when no points are visible due to filtering.
Caption:
  - Exports should make “0 visible points” obvious; confirm whether the export is empty by design or indicates a filtering mistake.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for exporting with 0 visible points.
:width: 100%

Edge cases like “all cells filtered out” should be explicitly visible in export preview/output.
```

---

## Next steps

- {doc}`07_troubleshooting_figure_export` (symptom → diagnosis → fix)
- {doc}`04_quality_knobs_and_best_practices` (large dataset strategies)
