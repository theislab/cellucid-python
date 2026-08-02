# Verified figure-export captures

Captured from the running web app at a 1440×1000 viewport with a device pixel
ratio of 2, light theme, against the **Pancreatic endocrinogenesis (scVelo)**
sample.

## Panel structure

```{figure} ../../../_static/screenshots/figure_export/panel-overview.png
:alt: The Figure Export panel with its four collapsed groups — Framing, Labels and Annotations, Style, Download — above the Export button, with the pointer on Export.
:width: 516px

The panel as it opens: one line of purpose, four collapsed groups, and
**Export**. The groups are mutually exclusive — opening one closes the others —
so you never see more than one of the captures below on screen at a time.
```

## Preview and framing

```{figure} ../../../_static/screenshots/figure_export/preview.png
:alt: Cellucid Figure Export preview for the deterministic 120-cell fixture with framing controls visible.
:width: 1440px

Preview exposes the current framing and visual state before any file is written.
```

```{figure} ../../../_static/screenshots/figure_export/framing.png
:alt: The Figure Export Framing group open, showing the plot size preset with its width and height boxes, the preview toggle with a colourblind-simulation selector and Refresh, and the frame-export crop controls with the aspect readout.
:width: 472px

**Framing.** The width and height boxes are the **plot content** size, not the
size of the file — see {doc}`04_quality_knobs_and_best_practices` for the
arithmetic that turns one into the other. `Show preview` has to be on before the
drag-to-crop frame can be used, which is why the crop hint sits beneath it.
```

## Labels, style, and output

```{figure} ../../../_static/screenshots/figure_export/labels-annotations.png
:alt: The Figure Export Labels and Annotations group with the Title box pre-filled with the dataset name followed by the coloured field and its visible categories, above the annotation checkboxes and the line listing elements that are never exported.
:width: 472px

**Labels & Annotations.** Two things worth noticing: the **Title:** box arrives
pre-filled from the dataset (see {doc}`05_metadata_and_provenance`), and the
grey line at the bottom names the on-screen decorations that are deliberately
never drawn into a figure.
```

```{figure} ../../../_static/screenshots/figure_export/style.png
:alt: Figure Export Style controls for background, typography, and text sizes.
:width: 224px

Style controls expose background, typography, and text sizing.
```

```{figure} ../../../_static/screenshots/figure_export/download.png
:alt: The Figure Export Download group open, showing the download format selector, the PNG DPI selector, the SVG point-strategy selector and the export-all-views checkbox, with the pointer on that checkbox.
:width: 472px

**Download.** Output format, resolution, an explicit SVG point strategy, and
whether all views are exported. A single output downloads in its native format;
a multi-file request downloads once as a ZIP containing the native files.
```

## When the app refuses

```{figure} ../../../_static/screenshots/figure_export/read-export-blocked-dialog.png
:alt: The Export blocked dialog over the app, headed "Export blocked", explaining that Cellucid will not create an export that differs from the active view, listing "Velocity overlay not exported" with its reason, and offering a single Back button.
:width: 1440px

Pressing **Export** with the velocity overlay running raises this. It is a
modal, not a toast, and its only button is **Back** — there is no “export
anyway”. The rule behind it is that a figure must match the view it claims to
be, so anything the exporter cannot reproduce blocks the export rather than
being silently dropped. See {doc}`06_edge_cases`.
```
