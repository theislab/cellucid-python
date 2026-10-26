# Verified core-interaction captures

## Navigation modes

The **Mode** dropdown is followed by a box of options belonging to that mode
alone, so the panel has three different shapes. All three are captured from the
Pancreas sample at twice the device pixel ratio.

```{figure} ../../../_static/screenshots/web_app/navigation-controls-orbit.png
:alt: The Navigation block with Mode set to Orbit and the pointer pressing the select: a Keyboard speed slider, Google Earth-style drag and Show orbit anchor.
:width: 488px

**Orbit.** Scenario `core-navigation-orbit`.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-planar.png
:alt: The Navigation block with Mode set to Planar: a Keyboard speed slider, Zoom to cursor and Invert axes.
:width: 488px

**Planar.** Scenario `core-navigation-planar`.
```

```{figure} ../../../_static/screenshots/web_app/navigation-controls-freefly.png
:alt: The Navigation block with Mode set to Free-fly: Look sensitivity and Move speed sliders with their readouts, Invert look axes and Projectile shooting, and Capture pointer with its (i).
:width: 488px

**Free-fly.** Scenario `core-navigation-freefly`.
```

## Render mode

```{figure} ../../../_static/screenshots/web_app/render-mode-select.png
:alt: The Render mode select reading Points, ringed and with the pointer pressing it, above Shader quality, Point size and the Depth cues sliders (3D lighting, Atmospheric fog, Perspective size scaling).
:width: 496px

`Render mode` sits under **Cells** in **Visualization** and has exactly two
options, `Points` and `Volumetric smoke cloud (beta)`. Almost everything below
it — `Shader quality`, `Point size` and the **Depth cues** sliders shown here —
belongs to the points renderer and is hidden in smoke mode; `Antialiasing` in **Performance**
is the exception and stays in both modes. Choosing smoke also shows a `Beta` pill beside the
`Render mode` label, which is why no pill is visible here. Scenario
`core-render-mode-select`.
```

```{figure} ../../../_static/screenshots/web_app/window-smoke-cloud.png
:alt: The whole window in Volumetric smoke cloud mode: the Pancreas cells drawn as a soft cloud tinted with the colours of their clusters, and the sidebar showing Render mode set to Volumetric smoke cloud (beta) with its Beta tag above the Smoke quality, Smoke shape, Smoke lighting and Performance controls, Animation speed at Still.
:width: 1440px

**Volumetric smoke cloud**, coloured by `clusters`. The cloud carries each
cluster's legend colour, the `Beta` pill is showing, and every smoke control is
at its initial value except `Animation speed`, held at 0 (`Still`) so the
frame stops changing and the capture is reproducible. Scenario
`core-window-smoke-cloud`.
```

## Switching embeddings

Each view shows one embedding, chosen with the **Shows** picker in **Compare
views**; the navigation mode follows the embedding's dimension unless you set
it yourself.

```{figure} ../../../_static/screenshots/web_app/window-dimension-3d.png
:alt: The whole window showing the Pancreas 3D embedding on the Studio backdrop, with the Navigation Mode select reading Orbit and, under Compare views, the Shows picker reading UMAP · 3D.
:width: 1440px

**3D**, the Pancreas default. Scenario `core-window-dimension-3d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-2d.png
:alt: The same window after choosing UMAP · 2D in the Shows picker, with the pointer on it: the cells now lie on a plane seen face-on over a dotted backdrop, and the Navigation Mode select has followed the view to Planar.
:width: 1440px

**2D**. The coordinates changed and the navigation mode
followed them to `Planar`. Scenario `core-window-dimension-2d`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-1d.png
:alt: The same window after choosing UMAP · 1D in the Shows picker: the cells are spread along one axis, and the Navigation Mode select has followed the view to Planar.
:width: 1440px

**1D**, also driven with `Planar`. Scenario `core-window-dimension-1d`.
```

## The 2D controls

```{figure} ../../../_static/screenshots/web_app/dimension-navigation-controls-2d.png
:alt: The Navigation block reading Planar above the Compare views block whose Shows picker reads UMAP · 2D: choosing a 2D embedding for a view switched its navigation to Planar.
:width: 480px

Both controls for a 2D embedding: **Navigation** on `Planar` and **Shows** on
*UMAP · 2D*. Unless you have picked a mode for the view yourself, showing a 2D
embedding is what puts `Planar` there — see
{doc}`05_dimension_switching_1d_2d_3d`. The canvas for the same configuration
is the **2D** capture above.
```

## Comparing views

```{figure} ../../../_static/screenshots/web_app/multiview-two-panels.png
:alt: The whole window in Grid compare with the sidebar open: two equal panels to the right of the sidebar show the Pancreas 3D embedding from the same linked camera — view 1 coloured by cell_type, view 2 recoloured by S_score on a Viridis scale — each under a title chip naming its embedding and field. The sidebar shows Compare views with the two view rows, each marked linked and Orbit, View 2 shows UMAP · 3D, Keep view, Link cameras, Layout Grid compare and Clear, and Coloring & filtering with S_score in the Continuous obs select.
:width: 1440px

Two views compare one dataset side by side; each row under **Views** names
the embedding the view shows and its field. The panels split the area to the
right of the sidebar equally, and each keeps its own coloring: view 1 by
cell_type, the selected view 2 by S_score.
```

## Connectivity edges

```{figure} ../../../_static/screenshots/web_app/connectivity-edges-controls.png
:alt: The whole window with Show connectivity edges on: the neighbour graph drawn as thin lines between the Pancreas cells, and the sidebar showing the Connectivity block with Edge color, Edge opacity, Line width and Max edges.
:width: 1440px

When a connectivity graph is present, the **Connectivity** block of
**Visualization** exposes edge color, opacity, line width, and the maximum
rendered edge count, with the counts of edges in total, visible and shown.
```

```{figure} ../../../_static/screenshots/web_app/connectivity-multiview.png
:alt: Two equal panels filling the whole window with the sidebar closed (the ☰ button at the top left), both showing the Pancreas 3D embedding coloured by cell_type from the same linked camera, each under a title chip, and in both the neighbour-graph edges, drawn at line width 2.0, form a grey web between the cells.
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
:alt: The Camera path panel with no saved positions: Keyframes with Save position and Clear all, and the Playback settings below.
:width: 532px

With fewer than two saved positions, the top playback bar is absent and camera
motion remains inactive.
```

```{figure} ../../../_static/screenshots/web_app/camera-path-configured-panel.png
:alt: The Camera path panel after two positions were saved: keyframes KF 1 and KF 2 with the KF 1 → KF 2 duration on auto, End with a return to the first keyframe, Set all with Apply, and Playback with Default speed, Position interpolation Catmull-Rom spline, Rotation interpolation SLERP (smooth), Easing Linear, Loop playback and Autoplay.
:width: 532px

Two valid keyframes expose timing and interpolation settings in the sidebar.
Loop playback and Autoplay follow the saved path settings.
```

```{figure} ../../../_static/screenshots/web_app/camera-path-transport-visible.png
:alt: The whole window after two camera positions were saved: the playback bar across the top of the canvas with its Play button and timeline markers, and the Camera path panel in the sidebar.
:width: 1440px

The top transport appears only after two valid keyframes exist. Playback waits
for **Play** unless the user explicitly enables **Autoplay**.
```
