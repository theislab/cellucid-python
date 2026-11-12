# Figure export edge cases

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
- an empty plot region (no points) with the notice "No visible cells" in place of axes,
- the legend and title still render,
- the exported artifact is “valid” but scientifically meaningless unless intentional.

How to handle:
- If you expected points: check filters and visibility first (see {doc}`../e_filtering/index`).
- If this is intentional (e.g., demonstrating a filter failure mode): consider adding an explicit caption note.

### Exporting tiny groups (rare populations, single-cell highlights)

Common failure mode: the “important” points disappear in an optimized export or become visually insignificant.

Best practices:
- Use **Emphasize selection** (mute non-selected points) rather than relying on the rare points being obvious.
- Avoid **Optimized vector** if the exact presence of rare points matters; prefer **Hybrid** or **PNG**.
- Export at a higher DPI, or frame the region (**Frame export**) so single points are larger in the plot, and anti-aliasing is not lost.

### NaN/Inf and invalid values

If positions/colors contain NaN/Inf:
- exact export validation terminates before rendering.

If you see odd axes or missing regions:
- verify your dataset file does not contain NaN/Inf in embeddings,
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
- give the grid enough size for its panel legends, or the figure comes out
  larger than the size you typed: a panel legend is never dropped or cut, and
  one wider than its plot (or every one, when beside their plots they would
  not fit the size) stands below it.

---

## Scale edge cases (big dataset, huge legends)

### An export that is simply too large

There are three ceilings, and they fail at different moments.

**The field's own limit.** `W` takes a width in millimetres of at least
**10** (a plot and a line of text), with no widest. A width below that never
starts an export; the box rejects it in place with `Figure width (mm) must be
at least 10.`, or `Figure width (mm) must be a decimal number.` for one that
is not a number.

**The GPU's limit.** The point layer is re-rendered through WebGL at export
resolution, so the *scaled* plot rectangle has to fit in a drawing buffer the
graphics driver will allocate. This is the one that catches people, because the
number that must fit is not the width you typed — it is the plot's width in
inches times the DPI. A 1600 mm plot at 600 DPI asks for 37,795 px, and a
typical limit is 16,384. When it does not fit you get a toast reading:

> Export failed: Figure export 37500×… exceeds this GPU's exact 16384×16384 raster limit.

and, in the rarer case where the buffer is accepted but comes back short:

> Export failed: Figure export could not allocate the exact …×… WebGL2 drawing buffer.

The limit is the machine's, not Cellucid's, so the same settings can succeed on
one computer and fail on another.

**The browser's canvas limit (PNG).** A PNG is composed on one browser canvas
the size of the whole figure, legend included, at `DPI / 96`. Every browser
refuses a canvas past a side length or an area of its own, lower on a machine
short of memory, and none of them says so. Cellucid asks the browser directly
— it allocates the canvas and checks that a pixel drawn in its far corner is
really there — and refuses the export before anything is rendered when it is
not. When the figure would fit without its legend, the toast reads:

> Export failed: The legend of 4,000 categories needs a …×… px figure, more than this browser can draw; turn the legend off or show fewer categories.

and otherwise:

> Export failed: This figure is …×… px, more than this browser can draw; choose a smaller size or a lower DPI.

with the figure's size in device pixels. The legend is never shrunk or cut to fit instead. An
SVG is written as text and has no such limit.

**What to do.** Lower DPI before lowering the width: 300 DPI is the usual journal
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
- huge file sizes,
- legends that dominate the layout: every entry is drawn whole, so the legend
  (and the figure) grows to hold them all, and past the browser's canvas limit
  a PNG export is refused (see above).

Strategies:
- Switch legend position to **Bottom**: the legend then fills the plot's width
  and grows downwards instead of sideways.
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
- wide figures,
- and high DPI PNG export.

Strategies:
- Export fewer panels at a time.
- Reduce DPI for drafts (150 DPI), then re-export final at 300 DPI.
- Widen the figure only as needed so each panel remains readable (a double column holds a 2 × 2 grid well).
- Remember that panels which disagree each carry their own legend inside their
  cell, which costs plot area (or makes the figure larger) and adds vector
  content per panel. Making the panels agree (same field, same hidden
  categories) collapses them back to one shared legend beside the grid.

---

## Rendering edge cases (smoke/fog/overlays)

### Export blocked: the view contains something the figure cannot

Some things the viewer draws have no representation in an exported figure. When
one of them is on, **Export** does not produce a degraded figure and does not
warn afterwards — it stops before rendering anything and says what is in the
way.

```{figure} ../../../_static/screenshots/figure_export/read-export-blocked-dialog.png
:alt: The Export blocked dialog over the app: "Cellucid will not create an export that differs from the active view", Connectivity overlay not exported with its reason, and the buttons Back and Turn off edges and export.
:width: 1440px

Pressing **Export** with the connectivity edges shown. One heading, one reason
per blocker, **Back**, and **Turn off edges and export**, which turns the edges
off and exports the points as shown — there is no “export anyway”.
```

The blockers, and what each one wants you to do:

| Blocker | What to do |
|---|---|
| **Connectivity overlay not exported** | Same reasoning for the neighbor-graph edges. Turn `Show connectivity edges` off. |
| **_<Mode>_ render mode not exported** | The active render mode has no figure renderer. Switch back to `Points`. |
| **Velocity overlay state unavailable** (whether the vector field is shown cannot be read, so whether to draw its arrows) / **Active render mode unavailable** / **Active viewer background unavailable** | The exporter could not read the state it needs in order to reproduce the view faithfully, so it refuses rather than guessing. Reload before exporting. |

This is a different mechanism from the **Never exported** line in the Labels &
Annotations group. That line names four interaction aids — orbit compass, lasso
and selection radius indicators, multiview panel title chips, projectile sandbox
— which are dropped silently by design because they are aids, not measurements.
A blocker is something that would change what the figure appears to *say*.

### Reference grid

The viewer's **Background** control decides two things at once: the color the
canvas is cleared to, and which reference layer is drawn behind the data: the
matplotlib-style grid box (`Grid (light)` and `Grid (dark)`), the studio floor
and dot surface (`Studio (light)` — the default — and `Studio (dark)`), a
painted look (`Atelier`, `Nebula`, `Abyss` and `Aurora`, light and dark), or
nothing (`White` and `Black`).

Exports reproduce that layer. The grid box is drawn as **vector rules** in SVG
(editable `<line>` elements you can recolor or delete in Illustrator or
Inkscape) and as the same geometry on the canvas in PNG, so both formats carry
the identical grid. The rules take their color, spacing, line weight, and
per-plane fade from the viewer's active background, so a light-grid view
exports a light grid and a dark-grid view exports a dark one.

The studio backdrop is reproduced from the same constants and functions its
shader is built from, with the floor height, stage center and screen density of
the view it was drawn in:

- floor rules are `<line>` elements stroked with a gradient sampled along each
  rule, because a rule fades with distance from the stage, with the viewing
  angle, and as its spacing family crossfades;
- dots, of the dot surface (also the one a level 3D view sees behind the
  data) or of a floor seen from straight above, are circles;
- the soft glow at the center of the view and the pool under the data are
  radial gradients. The glow lightens the viewer's own surface color, so it is
  drawn only with **Background: Match viewer** (or a custom color equal to
  it); the pool is drawn over any background.

The known differences are below one 8-bit color step or a quarter of a CSS
pixel: the gradient between a rule's 24 color stops is interpolated linearly,
a rule's width is the one it has where it is strongest, and the screen's
sub-step dither is not drawn. After a dimension or navigation change the screen
crossfades between floor and dots for 0.3 s; an export shows the settled
surface. The handover between floor and dot surface as a 3D view levels
depends on the camera alone, so an export taken mid-handover matches the
screen.

**Reference grid** in the Annotations block turns the layer off when you want a
figure without it. It is disabled when the viewer background draws none.

If the background control cannot be read at all, export is blocked with
**Active viewer background unavailable**: an export that cannot tell whether
the viewer draws a reference layer can neither include nor omit it honestly.


#### Background looks

A look (`Atelier`, `Nebula`, `Abyss`, `Aurora`) is a painted picture, not line
geometry, so it cannot become vector elements. The export draws it again with
the viewer's own renderer and shaders, for the figure's plot area, from the
exact state the view shows: the same look, the same colors taken from the
cells, the same contact shadow (Atelier) and the same camera, with each pixel
placed exactly as the points are placed, so the look stays registered with the
cells under any crop. PNG carries it at the file's resolution; SVG embeds it as
one image under the points, at twice the plot size. Star cores, paper grain and
marine snow keep the size they have on screen.

The look is drawn only when the figure's **Background** is **Match viewer**: on
**White**, **Transparent** or a custom color the figure keeps the paper you
asked for. The fade after choosing a look is not part of the figure.

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

### Axes on a 3D view

A 3D view has no axes, and no axis margins, whatever **Axes** says; nor does a
1-D or 2-D view tilted by orbiting. Axes are drawn only where they are rulers:
on a 1-D or 2-D view seen flat (planar navigation), whose cells all lie on the
plane the camera looks straight down at. A perspective picture shows cells at
many depths, so one length on the plot is a different distance at every depth,
and ticks would be read as UMAP coordinates the picture does not have. Keep
**3D orientation** on to show a 3D view's directions, and export a 2D view in
planar navigation when the figure needs embedding coordinates.

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

### The vector field: arrows on a grid

The velocity overlay animates particles, which no figure can hold. With
`Show overlay` on, a figure (PNG, SVG, the preview and an agent's export
picture alike) draws a still picture of the same flow instead: arrows on a
grid over the plot, as scVelo's `velocity_embedding_grid` draws them.

- **Same field, cells and camera.** The arrows come from the field the overlay
  shows, for every panel whose embedding it is bound to, from the cells the
  panel shows (its filters and hidden categories) with a vector, carried
  through the panel's own camera onto the plot. A panel showing another
  embedding draws none, as on screen.
- **How each arrow is made.** The plot is covered by a grid of 40 nodes along
  its longer side. Each node takes the Gaussian-weighted mean of nearby cells'
  vectors (a width of half a grid step), divided by the larger of 1 and the
  weight it gathered, so a node few cells reach draws a short arrow; a node
  that gathered under 1% of the 99th percentile weight draws none. Arrows
  are scaled by their mean length (half a grid step), as matplotlib's quiver
  autoscale that scVelo uses, so a fast region draws arrows longer than a step.
- **Colors.** An arrow is colored as the overlay colors a particle: its mean
  speed over the field's largest, through the overlay's `Color scheme`.
  Arrows are outlined in the paper's color so they read on any background.
- **Key.** The key names the field and how it is drawn (`velocity: mean
  arrows of nearby cells`) and gives the colors' speed range in the
  embedding's units (`speed 0 – 1.23`). With **Legend** on it is two rows of
  the legend — an arrow and what it is, the colors' ramp and their speed —
  under the entries and the highlight's row, laid out with the legend so it
  never covers a cell. In a grid, panels that draw one field with one speed
  range share it in the shared legend; panels whose keys differ each list
  theirs in their own legend. With **Legend** off the key stands in the plot's
  bottom-right corner, its line wrapped at its words to a narrow plot.
- **Provenance.** The figure's description and its JSON record name the field
  each panel draws (`Vector field: velocity (arrows: …)`).

The overlay's look (particle size, trails, bloom, color grading) is the
animation's and does not change the figure.

### Viewer chrome is never exported

Navigation aids are not measurements and no figure carries them. The orbit
compass, lasso and selection-radius indicators, multiview panel title chips, and
the projectile sandbox are always absent from an export. The Annotations block
states this so it is not discovered only after the file is opened. These do
**not** block an export — the compass in particular is on by default in orbit
mode, and blocking on it would make ordinary 3D exports impossible.

### Dark-background figures

Titles, tick labels, axis labels, panel letters, legend text, the plot frame,
and the selection badge take their color from the figure's own background. A
figure exported with **Background: Match viewer** from `Grid (dark)` or `Black`, or
with a dark **Custom** color, is written in light ink and stays readable. A
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
