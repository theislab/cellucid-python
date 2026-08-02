# Core Interactions (Viewer Basics)

These pages explain how to *drive the viewer* confidently: navigation, camera, render modes, multiview, and dimension switching.

**Audience:** everyone (novices get click-by-click; power users get exact controls + edge cases)  
**Time:** 15–30 minutes for full coverage

**Recommended reading order**

1) {doc}`01_navigation_modes_orbit_planar_free_fly`
2) {doc}`04_view_layout_live_snapshots_small_multiples`
3) {doc}`05_dimension_switching_1d_2d_3d`
4) {doc}`03_render_modes_points_vs_volumetric_smoke`
5) {doc}`02_camera_controls_advanced` (when you want to tune)
6) {doc}`06_troubleshooting_core_interactions` (when something feels “stuck”)
7) {doc}`07_screenshots` (every control on this page, photographed)

## Where these controls live

Everything in this section is in two sidebar accordions, **Visualization** and
**Compare Views**. Inside **Compare Views**, the **Navigation** block is the one
you will use most, and it changes shape with the mode you pick:

```{figure} ../../../_static/screenshots/web_app/navigation-controls-orbit.png
:alt: A panel headed NAVIGATION with a small circled i button, a MODE dropdown reading Orbit with the mouse pointer pressing it, and a bordered box below holding a KEYBOARD SPEED slider reading 0.40x, a ticked GOOGLE EARTH-STYLE DRAG checkbox and a ticked SHOW ORBIT ANCHOR checkbox.
:width: 472px

The **Navigation** block in Orbit mode, the default for a 3D embedding. The
`Mode` dropdown is the control this whole section revolves around.
```

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`move-to-end;1.5em;sd-mr-1` Navigation Modes
:link: 01_navigation_modes_orbit_planar_free_fly
:link-type: doc

Orbit vs planar vs free‑fly, with exact mouse/keyboard mappings and focus rules.
:::

:::{grid-item-card} {octicon}`sliders;1.5em;sd-mr-1` Camera Controls
:link: 02_camera_controls_advanced
:link-type: doc

Look sensitivity, movement speed, pointer lock, orbit anchor, and reset behavior.
:::

:::{grid-item-card} {octicon}`eye;1.5em;sd-mr-1` Render Modes
:link: 03_render_modes_points_vs_volumetric_smoke
:link-type: doc

Points vs volumetric smoke, including performance knobs and failure modes.
:::

:::{grid-item-card} {octicon}`columns;1.5em;sd-mr-1` Live + Snapshots
:link: 04_view_layout_live_snapshots_small_multiples
:link-type: doc

Keep view, grid compare, camera locking, badges, and per‑view state.
:::

:::{grid-item-card} {octicon}`versions;1.5em;sd-mr-1` Dimension Switching
:link: 05_dimension_switching_1d_2d_3d
:link-type: doc

How 1D/2D/3D switching works (per view), what changes, and what doesn’t.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 06_troubleshooting_core_interactions
:link-type: doc

Symptom → diagnosis → fix for “can’t rotate”, “controls disappear”, “smoke blank”, “context lost”, and more.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Verified captures
:link: 07_screenshots
:link-type: doc

The navigation, render-mode, dimension, multiview and camera-path controls as
they appear in the running app.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
