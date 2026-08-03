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
- exact export validation terminates before rendering.

If you see odd axes or missing regions:
- verify your exported dataset does not contain NaN/Inf in embeddings,
- try exporting a smaller subset or a different view to isolate the issue.

### Missing/renamed/deleted fields

Exports rely on the **active color-by field** and its legend model.
If the field was renamed/removed:
- a panel with no active field is drawn with no legend (there is nothing to
  explain) and its provenance record carries `field: {key: null, kind: null}`,
- your export may not match what you thought you were exporting.

If the field merely *differs across panels*, nothing is lost: each panel is
drawn with its own legend instead of one shared legend (see
{doc}`02_export_ui_walkthrough`). What changes is the figure-level metadata —
`Color Field` and the filename's color-field segment are omitted, because no
single field describes the figure.

Workarounds:
- confirm the active field immediately before exporting,
- for multiview exports, keep the same field *and* the same hidden categories
  across panels when you want one shared legend instead of per-panel legends,
- give the grid enough plot size that each panel legend has room to be drawn
  (below roughly 208 px of usable cell width for a right legend, a panel is
  exported without its legend).

---

## Scale edge cases (big dataset, huge legends)

### An export that is simply too large

There are two ceilings, and they fail at different moments.

**The field's own limit.** `W` and `H` accept a whole number between **100** and
**20000**. A value outside that never starts an export; the input rejects it in
place with `Width must be between 100 and 20000.` (or the height equivalent),
or `Width must be a whole number.` for a non-integer.

**The GPU's limit.** The point layer is re-rendered through WebGL at export
resolution, so the *scaled* plot rectangle has to fit in a drawing buffer the
graphics driver will allocate. This is the one that catches people, because the
number that must fit is not the number you typed — it is the plot width times
`DPI / 96`. A 6000 px plot at 600 DPI asks for 37,500 px, and a typical limit is
16,384. When it does not fit you get a toast reading:

> Export failed: Figure export 37500×… exceeds this GPU's exact 16384×16384 raster limit.

and, in the rarer case where the buffer is accepted but comes back short:

> Export failed: Figure export could not allocate the exact …×… WebGL2 drawing buffer.

The limit is the machine's, not Cellucid's, so the same settings can succeed on
one computer and fail on another.

**What to do.** Lower DPI before lowering plot size: 300 DPI is the usual journal
requirement and 600 DPI doubles every dimension for no gain in a raster figure
that is already oversampled. Or export **SVG**, which has no DPI stage and no
raster ceiling for the vector strategies.

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

:::{note}
Hiding categories does **not** shrink the legend: hidden categories keep their
entry, marked `(hidden)` with a hollow swatch, so the figure cannot be mistaken
for an unfiltered one. Merging or deleting categories reduces legend size;
hiding them does not.
:::

### Many views/snapshots (multi-panel export too large/slow)

Multi-panel export can become huge when you combine:
- many panels,
- large plot sizes,
- and high DPI PNG export.

Strategies:
- Export fewer panels at a time.
- Reduce DPI for drafts (150 DPI), then re-export final at 300 DPI.
- Increase overall plot size only as needed so each panel remains readable.
- Remember that panels which disagree each carry their own legend inside their
  cell, which costs plot area and adds vector content per panel. Making the
  panels agree (same field, same hidden categories) collapses them back to one
  shared legend beside the grid.

---

## Rendering edge cases (smoke/fog/overlays)

### Export blocked: the view contains something the figure cannot

Some things the viewer draws have no representation in an exported figure. When
one of them is on, **Export** does not produce a degraded figure and does not
warn afterwards — it stops before rendering anything and says what is in the
way.

```{figure} ../../../_static/screenshots/figure_export/read-export-blocked-dialog.png
:alt: The Export blocked dialog over the app, headed "Export blocked", explaining that Cellucid will not create an export that differs from the active view, listing "Velocity overlay not exported" with its reason, and offering a single Back button.
:width: 1440px

Pressing **Export** with the velocity overlay running. One heading, one reason
per blocker, one **Back** button — there is no “export anyway”.
```

The blockers, and what each one wants you to do:

| Blocker | What to do |
|---|---|
| **Velocity overlay not exported** | The figure would show the point layer without the flow the screen shows. Turn `Show overlay` off, export, turn it back on. |
| **Connectivity overlay not exported** | Same reasoning for the neighbour-graph edges. Turn `Show connectivity edges` off. |
| **_<Mode>_ render mode not exported** | The active render mode has no figure renderer. Switch back to `Points`. |
| **Velocity overlay state unavailable** / **Active render mode unavailable** / **Active viewer background unavailable** | The exporter could not read the state it needs in order to reproduce the view faithfully, so it refuses rather than guessing. Reload before exporting. |

This is a different mechanism from the **Never exported** line in the Labels &
Annotations group. That line names four interaction aids — orbit compass, lasso
and selection radius indicators, multiview panel title chips, projectile sandbox
— which are dropped silently by design because they are aids, not measurements.
A blocker is something that would change what the figure appears to *say*.

### Reference grid

The viewer's **Background** control decides two things at once: the colour the
canvas is cleared to, and whether the matplotlib-style reference grid box is
drawn behind the data. `Grid (light)` — the default — and `Grid (dark)` draw
one; `White` and `Black` do not.

Exports reproduce that grid. It is drawn as **vector rules** in SVG (editable
`<line>` elements you can recolour or delete in Illustrator or Inkscape) and as
the same geometry on the canvas in PNG, so both formats carry the identical
grid. The rules take their colour, spacing, line weight, and per-plane fade
from the viewer's active background, so a light-grid view exports a light grid
and a dark-grid view exports a dark one.

**Reference grid** in the Annotations block turns it off when you want a
grid-free figure. It is disabled when the viewer background draws no grid.

If the background control cannot be read at all, export is blocked with
**Active viewer background unavailable**: an export that cannot tell whether
the viewer draws a grid can neither include nor omit it honestly.

### 3D shader-accurate points and SVG limitations

Pure SVG circles cannot represent certain “shader-accurate” point appearances (e.g., 3D sphere shading).
Cellucid never changes the selected strategy:
- **Full vector** and **Optimized vector** encode the explicitly requested SVG
  circles.
- **Hybrid** and **PNG** use the shader-rendered point pass and require WebGL2
  plus complete camera matrices.

If you need the exported points to match the on-screen appearance:
- export **PNG** or **Hybrid SVG**,
- ensure the current environment exposes WebGL2 before exporting.

### Atmospheric fog and 3D lighting in vector exports

**Atmospheric fog** is a depth cue, not decoration: it fades and thins points
with distance from the camera, and it is on by default. Every export carries
it — PNG and Hybrid through the shader, Full Vector and Optimized Vector
through the same Beer-Lambert arithmetic applied to each circle — so a vector
figure reads its depth axis the same way the screen does.

**3D lighting** (the sphere-impostor shading) has no flat-disc equivalent and is
reproduced only by **PNG** and **Hybrid**. It is applied equally to every point
regardless of depth, so a vector figure differs in overall shading but never in
the relative reading of one cell against another.

### Depth ordering in dense 3D views

If points look “inside out” or layering seems wrong:
- ensure **Depth sort** is enabled for export.

### Axes interpretation in 3D

In 3D orbit mode, axes are reported in **camera-space** coordinates:
- useful for orientation,
- not the same as “UMAP_1/UMAP_2” in the strict embedding coordinate sense.

If you need axes that correspond to embedding coordinates, export a 2D/planar view.

### Volumetric smoke render mode not exported

Figure export renders the **point layer** only. If the viewer's render mode is
**Volumetric smoke cloud**, export is blocked with **Volumetric smoke cloud
render mode not exported** rather than producing a point-cloud image stamped
with metadata claiming it is the current view — that would be a different
figure, not a degraded one.

Switch the render mode back to **Points** before exporting. If the render-mode
control cannot be read at all, export is blocked the same way (**Active render
mode unavailable**), because an export that cannot confirm what the viewer is
drawing cannot promise to reproduce it.

### Connectivity overlay not exported

Even if connectivity lines are visible in the viewer, exports currently include:
- the **point layer** (and related point-based overlays),
- **not** the connectivity edges.

If edges are essential for your figure, capture them via another workflow and document how they were generated.

### Velocity vector field not exported

The velocity overlay is a data layer for the same reason connectivity is: it
shows a measured field. Figure export renders the point layer only, so an
enabled velocity overlay blocks the export with **Velocity overlay not
exported** rather than silently producing a figure without the flow the active
view shows. Turn the overlay off before exporting.

### Viewer chrome is never exported

Navigation aids are not measurements and no figure carries them. The orbit
compass, lasso and selection-radius indicators, multiview panel title chips, and
the projectile sandbox are always absent from an export. The Annotations block
states this so it is not discovered only after the file is opened. These do
**not** block an export — the compass in particular is on by default in orbit
mode, and blocking on it would make ordinary 3D exports impossible.

### Dark-background figures

Titles, tick labels, axis labels, panel letters, legend text, the plot frame,
and the selection badge take their colour from the figure's own background. A
figure exported with **Background: Match viewer** from `Grid (dark)` or `Black`, or
with a dark **Custom** colour, is written in light ink and stays readable. A
transparent figure keeps the light ink, because the page it will be placed on
is unknown and paper is the safe assumption.

---

## Environment edge cases

### WebGL2 / GPU differences

PNG and Hybrid SVG exports try to use WebGL2 for shader-accurate point rendering.
If WebGL2 is unavailable or restricted (common in locked-down environments):
- export is blocked with **Exact point export unavailable**.

### Fonts and text layout differences

SVG exports reference font families. If the font is missing on another machine:
- the editor substitutes a font,
- text widths change and legends may shift.

Mitigations:
- use common fonts (Arial/Helvetica),
- or convert text to outlines during final figure production (after you are confident in the layout).

### Browser download restrictions

If clicking Export produces no file:
- your browser may be blocking downloads (popup/download settings),
- a corporate browser policy may restrict file creation,
- extensions can interfere.

Use deterministic diagnostics:
- check the site download permission and the browser download list,
- inspect `[FigureExport]` console errors,
- submit one small PNG request to separate rendering capacity from file size,
- for a batch, expect exactly one ZIP download rather than multiple automatic
  downloads.

### Privacy/provenance surprises

Exports can include:
- dataset names/ids in filenames,
- embedded metadata that may include source URLs or the name of the local folder
  or file you loaded. Never a filesystem path — see
  {ref}`figure-export-source-file-privacy` for the exact inventory.

If you plan to share figures publicly, review {doc}`05_metadata_and_provenance` for inspection/stripping workflows.

---

## Next steps

- {doc}`07_troubleshooting_figure_export` (symptom → diagnosis → fix)
- {doc}`04_quality_knobs_and_best_practices` (large dataset strategies)
