# Verified core-interaction captures

## Navigation modes

The **Mode** dropdown is followed by a box of options belonging to that mode
alone, so the panel has three different shapes. All three are captured from the
Pancreas sample at twice the device pixel ratio.

```{figure} ../../../_static/screenshots/web_app/navigation-controls-orbit.png
:alt: A panel headed NAVIGATION with a MODE dropdown reading Orbit under the mouse pointer, and a bordered box below holding a KEYBOARD SPEED slider reading 0.40x, a ticked GOOGLE EARTH-STYLE DRAG checkbox and a ticked SHOW ORBIT ANCHOR checkbox.
:width: 472px

**Orbit.** Scenario `core-navigation-orbit`.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-planar.png
:alt: The same panel with the MODE dropdown reading Planar and a bordered box holding a KEYBOARD SPEED slider reading 0.0075x, a ZOOM TO CURSOR (PINCH-STYLE) checkbox and an INVERT AXES checkbox.
:width: 472px

**Planar.** Scenario `core-navigation-planar`.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-freefly.png
:alt: The same panel with the MODE dropdown reading Free-fly and a taller bordered box holding LOOK SENSITIVITY and MOVE SPEED sliders reading 0.05x and 1.00 u/s, INVERT LOOK AXES and PROJECTILE SHOOTING checkboxes, and a CAPTURE POINTER checkbox with a small circled i beside it.
:width: 472px

**Free-fly.** Scenario `core-navigation-freefly`.
```

## Render mode

```{figure} ../../../_static/screenshots/web_app/render-mode-select.png
:alt: A control labelled RENDER MODE holding a dropdown that reads Points, ringed in blue with the mouse pointer pressing it, above a dashed box labelled DEPTH PERCEPTION whose first slider, POINT SIZE (LOG), reads 0.75.
:width: 480px

`Render mode` sits at the top of **Visualization** and has exactly two options,
`Points` and `Volumetric smoke cloud (alpha)`. Almost everything below it —
starting with the **Depth perception** box shown here — belongs to the points
renderer and is hidden in smoke mode; **Image quality** is the exception and
stays in both modes. Choosing smoke also shows an `Alpha` pill beside the
`Render mode:` label, which is why no pill is visible here. Scenario
`core-render-mode-select`.
```

## Dimension switching

```{figure} ../../../_static/screenshots/web_app/window-dimension-3d.png
:alt: The whole application window with a three-dimensional coloured point cloud inside a light grid box, a DIMENSION dropdown reading 3D and a MODE dropdown reading Orbit.
:width: 1440px

**3D**, the Pancreas default. Scenario `core-window-dimension-3d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-2d.png
:alt: The same window after choosing 2D; the cells lie on one plane seen face-on, the DIMENSION dropdown reads 2D and the MODE dropdown reads Planar with the planar options beneath it, with the mouse pointer pressing the dimension dropdown.
:width: 1440px

**2D.** The coordinates changed and the navigation mode followed them to
`Planar`. Scenario `core-window-dimension-2d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-1d.png
:alt: The same window after choosing 1D; the cells form a single straight line of coloured segments across the grid box, the DIMENSION dropdown reads 1D and the MODE dropdown reads Planar with the planar options beneath it.
:width: 1440px

**1D**, also driven with `Planar`. Scenario `core-window-dimension-1d`.
```

## The older 2D captures

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: A Compare Views control column showing a Dimension dropdown set to 2D above a Navigation Mode dropdown set to Planar with planar options beneath it.
:width: 246px

Both controls set for a 2D embedding: `2D` and `Planar`. Unless you have picked
a mode yourself, choosing the dimension is what puts `Planar` there — see
{doc}`05_dimension_switching_1d_2d_3d`.
```

```{figure} ../../../_static/screenshots/web_app/dimension-2d-planar-default.png
:alt: Cellucid displaying a two-dimensional embedding face-on with Planar navigation.
:width: 1440px

The same configuration on the canvas: with Planar navigation, a 2D embedding is
seen face-on as a flat, screen-aligned surface.
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
