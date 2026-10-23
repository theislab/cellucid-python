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
:alt: A control labelled RENDER MODE holding a dropdown that reads Points, ringed in blue with the mouse pointer pressing it, above a dashed box labelled DEPTH PERCEPTION whose first slider, POINT SIZE (LOG), reads 2.87.
:width: 480px

`Render mode` sits at the top of **Visualization** and has exactly two options,
`Points` and `Volumetric smoke cloud (beta)`. Almost everything below it —
starting with the **Depth perception** box shown here — belongs to the points
renderer and is hidden in smoke mode; **Image quality** is the exception and
stays in both modes. Choosing smoke also shows a `Beta` pill beside the
`Render mode` label, which is why no pill is visible here. Scenario
`core-render-mode-select`.
```

```{figure} ../../../_static/screenshots/web_app/window-smoke-cloud.png
:alt: The whole window in Volumetric smoke cloud mode: the Pancreas cells drawn as a soft, lit cloud tinted with the colours of their clusters inside the grid box, and the sidebar showing Render mode set to Volumetric smoke cloud (beta) with its Beta tag, then the Volumetric smoke controls with Grid density 128³, Ray quality High, Render resolution 0.50x, Noise detail 128³, Cloud density 1.00 and Animation speed Still.
:width: 1440px

**Volumetric smoke cloud**, coloured by `clusters`. The cloud carries each
cluster's legend colour, the `Beta` pill is showing, and every smoke control is
at its initial value except `Animation speed`, held at 0 (`Still`) so the
frame stops changing and the capture is reproducible. Scenario
`core-window-smoke-cloud`.
```

## Switching embeddings (earlier interface)

These captures show the earlier **Dimension** dropdown, which is now the
**View N shows** picker in **Compare views**; the coordinates and the
navigation mode change the same way.

```{figure} ../../../_static/screenshots/web_app/window-dimension-3d.png
:alt: The whole application window with a three-dimensional coloured point cloud inside a light grid box, a DIMENSION dropdown reading 3D and a MODE dropdown reading Orbit.
:width: 1440px

**3D**, the Pancreas default (earlier interface). Scenario `core-window-dimension-3d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-2d.png
:alt: The same window after choosing 2D; the cells lie on one plane seen face-on, the DIMENSION dropdown reads 2D and the MODE dropdown reads Planar with the planar options beneath it, with the mouse pointer pressing the dimension dropdown.
:width: 1440px

**2D** (earlier interface). The coordinates changed and the navigation mode
followed them to `Planar`. Scenario `core-window-dimension-2d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-1d.png
:alt: The same window after choosing 1D; the cells form a single straight line of coloured segments across the grid box, the DIMENSION dropdown reads 1D and the MODE dropdown reads Planar with the planar options beneath it.
:width: 1440px

**1D**, also driven with `Planar` (earlier interface). Scenario `core-window-dimension-1d`.
```

## The older 2D captures

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: A Compare Views control column showing a Dimension dropdown set to 2D above a Navigation Mode dropdown set to Planar with planar options beneath it.
:width: 246px

Both controls set for a 2D embedding: `2D` and `Planar` (earlier interface:
the Dimension dropdown is now the embedding picker). Unless you have picked a
mode for the view yourself, showing a 2D embedding is what puts `Planar` there
— see {doc}`05_dimension_switching_1d_2d_3d`.
```

```{figure} ../../../_static/screenshots/web_app/dimension-2d-planar-default.png
:alt: Cellucid displaying a two-dimensional embedding face-on with Planar navigation.
:width: 1440px

The same configuration on the canvas: with Planar navigation, a 2D embedding is
seen face-on as a flat, screen-aligned surface.
```

## Comparing views

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels fill the area to the right of the sidebar and show the Pancreas 3D embedding from one shared camera — the left panel, the live view, coloured by cell_type, and the right panel, a kept view, recoloured by S_score on a Viridis scale — each under a title chip naming its field. The sidebar, scrolled to Compare Views in an earlier interface, lists badge 1 cell_type and the selected badge 2 S_score above the Layout select set to Grid compare, and Coloring & Filtering shows S_score in the Continuous obs select.
:width: 1440px

Two views compare one dataset side by side (earlier interface: the rows now
also name the embedding each view shows). The panels split the area to the
right of the sidebar equally, and each keeps its own coloring: the live view by
cell_type, the selected kept view by S_score.
```

## Connectivity edges

```{figure} ../../../_static/screenshots/web_app/connectivity-edges-controls.png
:alt: Cellucid KNN Connectivity controls beside a 3D dataset with visible cyan graph edges.
:width: 1440px

When a connectivity graph is present, the Render controls expose edge color,
opacity, line width, and the maximum rendered edge count.
```

```{figure} ../../../_static/screenshots/web_app/connectivity-multiview.png
:alt: Two equal panels filling the whole window with the sidebar closed, both showing the Pancreas 3D embedding coloured by cell_type from one shared camera, and in both the KNN connectivity edges, drawn at line width 2.0, form a grey web between the cells.
:width: 1440px

Connectivity settings apply to every panel: the edges were switched on, and
widened to 2.0, while the kept view was the selected view, and the live
panel draws them too. With the sidebar closed the two panels split the whole
window.
```

## Camera paths

The **Camera path** accordion starts collapsed. Opening it does not start camera
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
