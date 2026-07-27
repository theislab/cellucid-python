# Verified core-interaction captures

## 2D dimension and navigation

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: Cellucid Compare Views controls showing 2D dimension and Planar navigation selected.
:width: 246px

A 2D dataset selects Planar navigation by default. Dimension and navigation
remain explicit controls if the user wants another valid configuration.
```

```{figure} ../../../_static/screenshots/web_app/dimension-2d-planar-default.png
:alt: Cellucid displaying a two-dimensional embedding with Planar navigation.
:width: 1440px

The corresponding 2D view preserves a flat, screen-aligned exploration surface.
```

## Comparing views

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: Cellucid showing two side-by-side views of the same dataset.
:width: 1440px

Two kept views support side-by-side comparison while retaining one dataset
identity.
```

## Connectivity edges

```{figure} ../../../_static/screenshots/web_app/connectivity-edges-controls.png
:alt: Cellucid KNN Connectivity controls beside a 3D dataset with visible cyan graph edges.
:width: 1440px

When a connectivity graph is present, the Render controls expose edge color,
opacity, line width, and the maximum rendered edge count.
```

```{figure} ../../../_static/screenshots/web_app/connectivity-multiview.png
:alt: Two Cellucid views showing the same point cloud and connectivity edges.
:width: 1440px

Connectivity rendering remains synchronized across kept views of the dataset.
```

## Camera paths

The **Camera Path** accordion starts collapsed. Opening it does not start camera
motion.

```{figure} ../../../_static/screenshots/web_app/camera-path-unconfigured.png
:alt: Camera Path panel with no saved keyframes and no playback bar.
:width: 246px

With fewer than two saved positions, the top playback bar is absent and camera
motion remains inactive.
```

```{figure} ../../../_static/screenshots/web_app/camera-path-configured-panel.png
:alt: Camera Path panel with two saved keyframes and playback settings.
:width: 246px

Two valid keyframes expose timing and interpolation settings in the sidebar.
Loop playback and Autoplay follow the saved path settings.
```

```{figure} ../../../_static/screenshots/web_app/camera-path-transport-visible.png
:alt: Cellucid with two camera keyframes and the playback transport visible above the canvas.
:width: 1440px

The top transport appears only after two valid keyframes exist. Playback waits
for **Play** unless the user explicitly enables **Autoplay**.
```
