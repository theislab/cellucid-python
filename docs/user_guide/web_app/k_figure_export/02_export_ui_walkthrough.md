# Export UI walkthrough

**Audience:** wet lab + computational users (figures for slides, reports, and papers)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Where the **Figure export** panel lives and how to open it
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

3) **Open** the **Figure export** panel (left sidebar).

4) **Choose the figure's width**
   - Pick the journal column it will be printed in — **Single column · 89 mm**, **1.5 column · 120 mm**, **Double column · 183 mm** *(default)* — or type a width in millimetres.
   - The figure is exactly that wide at every resolution; its height follows from the view's shape, and its text keeps its point size (see below).

5) **Choose download format**
   - **PNG** at **300** DPI (the DPI it starts at) is a strong default for slides and manuscripts.
   - **SVG** is the right choice if you want to edit labels/legend layout in Illustrator/Inkscape.

6) **(Optional) Crop without changing camera**
   - Enable **Show preview** → enable **Frame export** → drag the frame → click **Confirm** to apply.

7) Click **Export**. While the export runs, the button reads **Cancel export**: press it to stop the export where it is (nothing is downloaded, and nothing is reported as failed). A double click on **Export** (or a held Enter) exports once: its second click never presses **Cancel export**.

8) **Open the file outside Cellucid** and do a 10-second sanity check:
   - correct view (live vs snapshot),
   - correct legend/colors,
   - correct filters/visibility,
   - readable text at your intended size.

---

## Where the Figure export panel lives

Figure export is controlled from the **Figure export** panel in the left sidebar.

```{figure} ../../../_static/screenshots/figure_export/panel-overview.png
:alt: The Figure export panel with its four collapsed groups, Framing, Labels and annotations, Style and Download, above the Export button, with the pointer on Export.
:width: 532px

The panel opens with all four groups collapsed. They are mutually exclusive:
opening one closes the others, so work through them in order rather than
expecting to see them side by side.
```

**Which view gets exported?**

- By default, Cellucid exports the **active view**:
  - the live view, or
  - the snapshot view you currently have focused/selected.
- In multiview **grid compare**, click inside the panel you want to export first (focus matters).
- If you enable **Export all views (split-view)** and you are in grid compare mode, Cellucid exports a **multi-panel grid** (live + snapshots).

**Are export settings saved?**

- Yes: every setting in this panel (size, frame, background, fonts, labels, legend, preview, format, DPI, SVG strategy, **Export all views**) is a step of the history: undo takes it back, and a session file keeps it and puts it back when loaded.
- The exported artifact also embeds metadata (see {doc}`05_metadata_and_provenance`).

---

## Preview vs export

**Preview** is a lightweight, export-style approximation:

- It is meant to help you confirm **layout and framing** (especially cropping).
- It is intentionally **downsampled and debounced** so the UI stays responsive on large datasets.
- It supports **colorblindness simulation** (Normal / Deuteranopia / Protanopia / Tritanopia) to help you catch problematic palettes.

**Export** is the real thing:

- It renders from the current view buffers and camera matrices.
- It applies your chosen format and strategy (SVG full/optimized/hybrid or PNG).
- SVG has no preselected strategy. You must explicitly choose **Full vector**,
  **Optimized vector**, or **Hybrid**; Cellucid never changes that choice.
- **Export blocked** lists an exact missing requirement, such as WebGL2/camera
  matrices for PNG or Hybrid, a visible connectivity layer that the current
  exporter does not include, or a render mode other than **Points** (figure
  export draws the point layer only, so a volumetric smoke view would export as
  a different image, not a degraded one). When the only blocker is the
  connectivity edges, the dialog unblocks it in one step: **Turn off edges and
  export** turns them off, as their switch would, and exports. The velocity
  overlay does not block: the figure draws its field as arrows on a grid.

:::{tip}
If your preview looks wrong, click **Refresh** in the preview controls (especially after changing views, fields, or filters).
:::

---

## Setting the figure's size

A figure is sized as it will be printed: you give its **width**, and the
figure fills it exactly, at the height its plots make it; its text is set in
**points** and keeps that size at any width and resolution. In the UI you'll
see, under **Figure width**:

- A width preset — the column widths journals print figures at (Nature's;
  Cell, Science and PLOS print within a few millimetres of them):
  - **Single column · 89 mm**
  - **1.5 column · 120 mm**
  - **Double column · 183 mm** *(default)*
  - **Custom width** — any other width, typed in the **W** box (millimetres; a
    point or a comma, at least 10). Typing a preset's width selects the preset.
- One line under it with the figure's printed size and, when Export writes a
  PNG, its pixels at the chosen DPI: `183 × 130 mm · 2163 × 1541 px` (its pointer says the DPI). It reads the preview's layout; with the preview off, or with
  **Export all views**, it says the width and that the plots (or the panels)
  make the height.

```{figure} ../../../_static/screenshots/figure_export/framing.png
:alt: The Framing group open: Figure width on its default, Double column · 183 mm, with its W box at 183 mm and the printed size under it (183 mm wide · height from its plot), and Preview with its Show preview checkbox, above the Export button.
:width: 532px

The **Framing** group. The preset writes its width into the **W** box you can
type into; the line under it says the figure's printed size.
```

An agent, a notebook or a replay sets the same two settings with
`figure-export.set` (`size_preset`, `width` in millimetres, `font_sizes` in
points): "180 mm wide" is `width: 180`, never a conversion to pixels.

### The figure's height

- Each plot takes the **shape of the view on screen** — the pane's width to
  height, or the frame's when **Frame export** is on — and the figure is as wide
  as asked, so its height follows: a 16:10 pane in a double column gives a
  plot about 130 × 81 mm. The view fills its plot exactly, with the screen's
  own projection: there is no blank band inside or beside the plot, no cell
  the screen shows is cut off, and points keep the size they have relative to
  the view. To change the figure's shape, frame it (**Frame export**; **Lock
  aspect** keeps the frame's shape while you resize it) or use **Fit the
  cells**, which frames the cells.
- Axes, title and legend stand **inside** the width: the plot is as wide as the
  width leaves it, never shrunk below a readable 80 figure pixels (21 mm).
- The legend is as large as its entries need. Beside the plot (**Legend:
  Right**) it fills the plot's height in as many columns as its entries take,
  each column as wide as its longest name, and it is at least 10 ems of its
  text wide (25 mm at 7 pt, a color bar's length); below the plot (**Legend:
  Bottom**) it fills the plot's width in as many rows as they take. No category
  name is ever shortened with "…". A legend asked beside the plot that would
  make the figure wider than asked (eight long names, or the vector field's
  key, beside a short plot in a single column) stands below the plot instead,
  as a grid's legends do, and the export picture (`view.picture` with source
  "export") says so. A figure still wider than asked (a title, or text this
  large) is said there too. A title wider than the figure widens it.
- With **Export all views**, the panels keep the arrangement they have on
  screen (columns and rows, see
  {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`)
  and the figure fills the width: every panel's plot has one height and its
  own view's shape (views of one shape give panels of one size), the plots of a
  row and of a column line up, and the figure is its panels, labels and legends
  at the tallest plots that keep to the width. Legends take their room from the
  plots down to half the size the plots have without legends; past that the
  figure grows for them (see {ref}`legends-in-a-multi-panel-export`).
- When a plot is far from the shape of the cells it draws — a tall pane of a
  wide embedding (several views side by side, or the chat's side panel beside
  them) gives a tall figure with the cells in a band — the preview says so
  under the picture (`Plot narrower than the cells:`, why under its pointer) and
  offers **Fit the cells** (below), which frames the cells so each plot takes
  their shape.

A width below 10 mm is refused by its box before the export starts, with
`Figure width (mm) must be at least 10.`; a text size is a whole or half point
from 1 to 375.

### How big the file actually comes out

Worked example, single view of a 4:3 pane, title on, axes on, a short legend
(four categories) on the right, Auto text:

| Width | Figure (printed) | PNG at 300 DPI | The plot |
|---|---:|---:|---:|
| Single column, 89 mm | 89 × 60 mm | 1050 × 706 px | 35 × 27 mm |
| Double column, 183 mm | 183 × 130 mm | 2163 × 1541 px | 130 × 97 mm |

Around the plot stand its axis space (what the axes' text needs: here about 16
mm left of the plot, for two-digit ticks at 6 pt and 7 pt axis labels; see
{ref}`figure-export-axes`), the title row, the legend and the outer padding.
They are the same size at every width — they are text and its room — so they
take most of a single column and a third of a double one. An SVG is the same
size, written in millimetres (its `viewBox` in figure pixels, 1/96 inch), so
Illustrator, Inkscape and a journal's checks open it at its printed size.

A longer legend adds its own width: beside the plot, its entries take one row
per line of the plot's height and as many columns as that leaves them, each
column as wide as its longest name in the export's font.

### PNG DPI

A PNG is the figure's width in inches times the DPI:

`PNG width (px) = width (mm) / 25.4 × DPI` — 183 mm at 300 DPI is 2163 px.

Its `pHYs` chunk says the DPI, so a layout program places it at its printed
size. For drafts and slides, 150 DPI halves each side; 600 DPI doubles it.

### Auto text sizing

By default, **Auto text** is on: the text is set at the sizes journals print
figure text at — 7 pt (panel titles, legend, axis labels), ticks and centroid
labels 6 pt, the figure's title 8 pt — at any width. Turn it off to set your
own sizes, in points (whole or half points). Points keep their printed size:
the same 7 pt is 7 pt in a single column and in a double one, at 150 or 600 DPI.

---

## Titles, captions, legends, axes

### Title

- The **Title** field auto-fills with your dataset name and active field, and, when you have hidden some of its categories, which ones (as the Filters list says it). For one view it then names the filters on other fields as the view's label does ("Developing human lung (He) • LVL2 • LVL0: showing only Epithelial"): the legend lists no category those filters leave without cells. In a grid, each panel's label names its own filters.
- A view you have named (double-click its row in **Compare views**) is titled by its name instead: "Developing human lung (He) • Early epithelium" (see {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`).
- If you type your own title, auto-fill turns off.
- To return to auto-fill behavior, clear the title field completely.

### Legend

- **Legend** is enabled by default.
- Position can be:
  - **Legend: Right**
  - **Legend: Bottom**

Legend contents are sourced from the same legend model as the viewer, so categorical colors and continuous colormaps should match what you see.

Two things the exported legend spells out that the on-screen legend does not
have to, because a figure travels without the app around it:

- **Categories you hid keep their entry and are marked `(hidden)`**, with a
  hollow gray outline instead of their color. The entry stays so a filtered
  figure cannot be mistaken for a complete one; the color goes so nobody hunts
  the plot for points that were never drawn. A category that a filter on
  another field leaves with no cells (the row the on-screen legend shows
  disabled, "No cells available in this category after other filters") is
  not listed at all, ticked or not: the figure has none of its cells, and the
  title (a grid panel's label) names that filter ("LVL2 • LVL0: showing only
  Epithelial").
- **A logarithmic colorbar is labeled `Log10 color scale`** together with its
  true midpoint (the geometric mean of the two endpoints), so a reader cannot
  interpolate the endpoint labels linearly and quote a value that is wrong by
  orders of magnitude.
- **Highlighted cells get an entry.** The figure draws the cells of the
  highlight page you are on over the coloring in the highlight's color (cyan),
  each dot inside a black ring, as the screen does (in every format: PNG,
  Hybrid and the vector SVGs), so the legend lists them: a cyan swatch and
  `Highlighted: <page> (<N> cells)`, the page's name and its highlighted cell
  count. Without it, a reader of the figure takes that color for the coloring's
  (cyan cells over a gene's colormap read as "high expression"). The entry is
  one more row under the categories or the color bar, and the legend is laid
  out with it; a view colored by no field gets a legend of the entry alone. It
  appears only where highlighted cells are drawn: a panel whose filters hide
  them all has none, and with **Legend** off there is none. In a multi-panel
  export the shared legend carries it once, or each panel's own legend where
  that panel draws highlighted cells. To leave the highlight out of a figure,
  clear the page or switch to an empty one before exporting.

(figure-export-axes)=
### Axes

- **Axes** is enabled by default, and draws axes only where they are
  rulers: on a 1-D or 2-D view seen flat (planar navigation, the default for
  1-D and 2-D). A 3-D view, or a flat view tilted by orbiting, has no axes
  and no axis margins, whatever the checkbox says (its tooltip says so): a
  perspective picture shows cells at many depths, so one length on the plot
  is a different distance at every depth, and ticks would invite reading
  coordinates the picture does not have. Keep **3D orientation** on to show
  a 3-D view's directions.
- The axes are **rulers of the plot**: their ends are the coordinates the plot's
  edges show (the view's, or the frame's), so a tick sits over the cells at its
  value and the axes frame exactly what is drawn. Only the round values that
  fall between the ends get a tick, and a short axis (a small grid panel) gets
  fewer of them at a coarser round step, so its labels, measured in the
  figure's font, never run into each other.
- The axes' room follows their text: the tick labels are measured, the y-axis
  label stands past the widest of them and the x-axis label under them, and the
  margins grow with the **Ticks** and **Axis** sizes (any size from 1 to 375
  pt), so an axis label never runs into a tick label and nothing is cut at the
  figure's edge or covered by the legend or a panel's title. An axis label
  longer than its side starts at the plot's corner and runs on past the far
  end. The PNG and the SVG place every tick and label at the same place: both
  draw one layout of the axes.
- A plot with no visible cell shows "No visible cells" instead of axes.
- Axis labels default to **X** and **Y**.
  - Clear X and/or Y to hide axis labels entirely.

How to interpret axes: they are embedding coordinates, measured on the
embedding's plane, which the camera looks straight down at; every cell lies on
it, so a tick sits over the cells at its value.

### Reference grid

- **Reference grid** is enabled by default and exports the reference layer the
  viewer background draws: the matplotlib-style grid box of **Grid (light)**
  and **Grid (dark)**, or the floor and dot surface of **Studio (light)** and
  **Studio (dark)**. A painted look (**Atelier**, **Nebula**, **Abyss**,
  **Aurora**) is redrawn as an image under the points when the figure's
  background is **Match viewer**.
- It follows the viewer: a light-grid view exports a light grid, a dark-studio
  view exports a dark floor. Nothing is hard-coded.
- In SVG the layer arrives as ordinary vector elements — `<line>` rules (the
  studio floor's with gradient strokes), `<circle>` dots and gradient fills —
  so you can recolor, thin, or delete it in Illustrator or Inkscape. PNG
  carries the identical geometry.
- The checkbox is disabled when the viewer background (**White** or **Black**)
  draws no reference layer.

### 3D orientation + depth sort

- **3D orientation** adds a small orientation widget in 3D orbit mode. It
  takes the same share of the plot as it would of the on-screen view: its size
  is scaled by the plot's pixels over the view's, as the points are, so a small
  grid panel is not covered by it (its angle text shrinks with it).
- **Depth sort** is enabled by default and helps 3D point clouds export with the correct near/far ordering.
  - Turning it off can be faster, but may reduce fidelity in dense 3D views.

### Style controls that affect export readability

In the **Style** section you can choose:
- Background: **Match viewer**, **White**, **Transparent**, **Custom…**
  - Annotation ink follows the background you choose: a dark figure gets light
    titles, ticks, axis labels, and frame instead of near-black ones.
- Font family: Arial / Helvetica / Times (safe cross-platform defaults)
- Text sizing:
  - Auto text on/off
  - Base / Legend / Ticks / Axis / Title / Centroids text sizes, in points

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

It also stamps an **`n = 1,234 selected`** badge in the top-left corner of the
plot area, so the count travels with the figure. In a multi-panel export **every
panel carries its own badge**, counted from the cells that panel actually draws:
panels carry their own filters, so one number could not be true for all of them.

Practical use cases:
- A “selected cluster vs background” figure for a talk.
- Highlighting a rare population without hiding the rest of the dataset (hiding can be misleading).

:::{important}
Emphasize selection does nothing if there are no highlighted/selected cells.
:::

---

## Download options and large dataset strategies

```{figure} ../../../_static/screenshots/figure_export/download.png
:alt: The Download group open with the format set to SVG: the SVG point strategy select on Full vector, and Export all views (split-view) with the pointer on its checkbox, above the Export button.
:width: 532px

The **Download** group. The SVG point strategy opens on **Full vector**
(every point an editable circle); choose **Optimized vector** or **Hybrid** for
a large dataset. Cellucid never changes the strategy you chose.
```

### Download format options

In the **Download** section you can choose:

- **SVG**: one SVG export.
- **PNG**: one PNG export at the selected DPI.
- **SVG + PNG**: one ZIP containing both native files.
- **PNG (150/300/600)**: one ZIP containing three PNGs at the named DPIs.
- **All**: one ZIP containing the SVG plus PNGs at 150/300/600.

Cellucid renders and validates every batch member before creating the ZIP or
starting the single browser download. If any member fails, no archive is
downloaded.

### DPI (for PNG)

When exporting a single PNG, you can choose:
- 150 / 300 / 600 DPI

Higher DPI means:
- more pixels (larger files, more time),
- but better raster fidelity for printing and zooming.

### Explicit SVG point strategy

SVG export can become impractical when you have tens/hundreds of thousands of visible points.
When the selected download includes SVG, one strategy draws its points:
**Full vector** until you choose another (an agent's `svg_strategy` is the
same setting):

- **Full vector**  
  Exports every point as an SVG circle (maximum editability; can be huge/slow).
- **Optimized vector**  
  Uses deterministic density-preserving reduction with the explicit **Keep**
  target (initially 100,000; if fewer points are visible, all visible points
  remain available).
  Good for medium datasets when you need a fully-vector plot but can tolerate approximation.
- **Hybrid**
  Rasterizes points (WYSIWYG, including 3D shading) but keeps annotations vector.

:::{note}
Full and Optimized vector deliberately encode points as SVG circles. Hybrid
deliberately embeds the shader-rendered point pass. Cellucid executes the
strategy you selected and never substitutes another strategy or format.
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
4) (Optional) Click **Lock aspect** to keep the frame's shape while you resize it.
5) Click **Confirm** to apply the crop (the preview will switch to a framed preview).
6) Export.

The plot takes the frame's shape, filled by the frame exactly, and the figure
keeps its width, so the frame decides the figure's height — identically in
PNG, every SVG mode, and the preview.

A frame is drawn at the export's full resolution, not enlarged from the
whole view: its points are rasterized at the pixels they land on, so a
framed PNG (or the point layer of an SVG in Hybrid mode) is as sharp as an
unframed one of the same size, while every point keeps the size it has
relative to the frame. The one exception is a frame so small that
the GPU's largest point sprite cannot reach the size a point must have
there: the points are then rasterized at that largest zoom and the frame is
scaled the rest of the way.

Useful actions:
- **Fit the cells** (beside **Frame export**) frames every cell the figure
  draws, with a small margin (3% of their extent, and a point's radius): the
  plot takes the cells' shape instead of the pane's, so a narrow pane — the
  chat open, a grid of views — no longer gives a narrow plot with empty bands
  above and below a wide embedding. It turns the preview and **Frame export**
  on, leaves **Lock aspect** off, and never widens the view (only cells the
  pane shows count). With **Export all views**, one frame holds every panel's
  cells. When a plot is far (1.5 times or more) from the shape of the cells
  it draws, the preview says so in one line under the picture, with the same
  button. An agent does the same with `figure_export_set` and
  `frame: "cells"`; `figure_export_get` reads back the frame it made.
- **Reset** resets the crop to the full frame.
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
- Panels are labeled (e.g., A, B, C…) and titled by their view: a view you
  named is titled by its name (`A. Early`, `B. Mid`, `C. Late`), any other by
  its field and filters (`A. LVL2 • LVL0: showing only Epithelial`). To label
  the panels of a time course or a comparison, name each view first
  (double-click its row in **Compare views**).
- Every colored panel is explained by a legend — one shared legend beside the
  grid when the panels agree, otherwise a legend inside each panel (see below).
- With **Emphasize selection** on, each panel carries its own `n = N selected`
  badge, counted from that panel's own cells.
- Embedded metadata describes **every** panel, one record each (see
  {doc}`05_metadata_and_provenance`).

(legends-in-a-multi-panel-export)=
#### Legends in a multi-panel export

Cellucid draws **one shared legend beside the grid** only when every panel
really says the same thing: the same color-by field *and* a legend model equal
in every respect, **including which categories are hidden**. Two panels colored
by `cell_type` with different categories filtered out do not qualify — a single
legend would tell the reader both panels draw the same set of categories, which
is exactly the claim the `(hidden)` marking exists to prevent.

Whenever the panels do not qualify, **each panel gets its own legend drawn
inside its cell**. A panel legend never covers the plot, the axes, or the
`A. Live` panel label, and it never shortens a name.

Each panel is titled by its letter and its view's name (`B. Sox9 (ductal)`;
a view without a name by its field and filters). The letters come from the
figure, in panel order, so a view's name needs no letter of its own. A title
wider than its column widens that column, and the plots shrink to keep the
figure to its width: a title is never cut off.

How a grid is sized:

- **The figure fills its width, and the figure is its composition.** The
  panels are arranged as on screen; every panel's plot has one height and its
  view's shape; the figure is its panels, labels and legends at the tallest
  plots whose composition keeps to the width, and its height is theirs. Every
  legend — shared or per panel — is as large as its entries need, every entry
  whole, and takes its room from the plots down to half the height they have
  without legends; past that, or to keep every plot at least 80 × 64 figure
  pixels (21 × 17 mm), the figure comes out wider than asked rather than
  dropping or cutting anything.
- **The shared legend** stands beside the grid, its top level with the first
  row's plots (**Legend: Right**, at least 10 ems of its text wide, in as many columns as
  its entries take), or below it (**Legend: Bottom**, in as many rows as they
  take).
- **A panel legend stands beside its plot when it is no wider than the
  plot.** With **Legend: Right**, a panel legend fills its plot's height in as
  many columns as its entries take. One that would be wider than its plot (at
  the plot size the width gives the grid without legends) — a long list of long
  names, such as 69 cell types — stands **below its plot**, at its column's
  width, in as many rows as its entries take. When that arrangement would not
  fit the width (the legends would take more than half the plots' height, or
  widen the figure past it), every legend beside its plot is
  tried, then every legend below its plot, and the first that fits is drawn:
  two wide views with long legends in a narrow figure get their legends below
  their plots rather than plots shrunk to strips beside many-columned
  legends. With **Legend: Bottom**, every panel legend stands below its plot.
- **The plots line up.** The plots of a row start at one height and those of
  a column at one x; a column has axis margins only when one of its panels
  draws axes. A legend's frame is the box its entries need.

:::{tip}
If you want all panels to be directly comparable, show the same embedding in them and keep **Link cameras** on, so they share one camera (see {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`). Matching fields, filters, and hidden categories is also what buys you a single shared legend instead of one legend per panel.
:::

---

## Interface reference

```{figure} ../../../_static/screenshots/figure_export/preview.png
:alt: The whole window with Preview on in the Framing group: the sidebar shows the figure as it will be written, beside the live view on the canvas.
:width: 1440px

Preview exposes the current framing and visual state before a file is written.
```

---

## Next steps

- {doc}`03_export_formats_and_renderers` (what each SVG/PNG mode means)
- {doc}`04_quality_knobs_and_best_practices` (publication-grade settings)
- {doc}`07_troubleshooting_figure_export` (export button does nothing, huge SVG, etc.)
