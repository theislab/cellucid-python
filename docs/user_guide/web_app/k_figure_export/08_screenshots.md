# Verified figure-export captures

Captured from the running web app at a 1440×1000 viewport with a device pixel
ratio of 2, light theme, against the **Pancreatic endocrinogenesis (scVelo)**
sample.

## Panel structure

```{figure} ../../../_static/screenshots/figure_export/panel-overview.png
:alt: The Figure export panel with its four collapsed groups, Framing, Labels and annotations, Style and Download, above the Export button, with the pointer on Export.
:width: 532px

The panel as it opens: **Export** with its (i), four collapsed groups, and the
**Export** button. The groups are mutually exclusive — opening one closes the others —
so you never see more than one of the captures below on screen at a time.
```

## Preview and framing

```{figure} ../../../_static/screenshots/figure_export/preview.png
:alt: The whole window with Preview on in the Framing group: the sidebar shows the figure as it will be written, beside the live view on the canvas.
:width: 1440px

Preview exposes the current framing and visual state before any file is written.
```

```{figure} ../../../_static/screenshots/figure_export/framing.png
:alt: The Framing group open: Plot size on Screen (half) with its W and H boxes, and Preview with its Show preview checkbox, above the Export button.
:width: 532px

**Framing.** The width and height boxes are the **plot content** size, not the
size of the file — see {doc}`04_quality_knobs_and_best_practices` for the
arithmetic that turns one into the other. **Frame export**, the drag-to-crop
frame, appears only once `Show preview` is on (see the window capture above).
```

## Labels, style, and output

```{figure} ../../../_static/screenshots/figure_export/labels-annotations.png
:alt: The Labels and annotations group open: the Title box pre-filled, starting with the dataset name; the Annotations checkboxes Axes, Legend, Reference grid, 3D orientation and Depth sort; the X and Y axis labels; Legend: Right; and Emphasize selection.
:width: 532px

**Labels and annotations.** The **Title** box arrives pre-filled from the
dataset (see {doc}`05_metadata_and_provenance`); below it the **Annotations**
checkboxes, the axis labels, the legend position, and **Emphasize
selection**.
```

```{figure} ../../../_static/screenshots/figure_export/style.png
:alt: The Style group open: background, font family and the text sizes for the base text, legend, ticks, axis labels, title and centroid labels.
:width: 532px

Style controls expose background, typography, and text sizing.
```

```{figure} ../../../_static/screenshots/figure_export/download.png
:alt: The Download group open with the format set to SVG: the SVG point strategy select on its placeholder, shortened to "Choose S…", and Export all views (split-view) with the pointer on its checkbox, above the Export button.
:width: 532px

**Download.** Output format, then an explicit SVG point strategy for SVG (or
the DPI for PNG), and
whether all views are exported. A single output downloads in its native format;
a multi-file request downloads once as a ZIP containing the native files.
```

## When the app refuses

```{figure} ../../../_static/screenshots/figure_export/read-export-blocked-dialog.png
:alt: The Export blocked dialog over the app: "Cellucid will not create an export that differs from the active view", Velocity overlay not exported with its reason, and a single Back button.
:width: 1440px

Pressing **Export** with the velocity overlay running raises this. It is a
modal, not a toast, and there is no “export anyway”. Besides **Back**, a view
blocked only by the edges or the vector field offers **Turn off … and
export**, which turns those layers off and exports the points as shown
(the screenshot predates that button). The rule behind it is that a figure must match the view it claims to
be, so anything the exporter cannot reproduce blocks the export rather than
being silently dropped. See {doc}`06_edge_cases`.
```
