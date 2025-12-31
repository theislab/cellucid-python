# Core Interactions (Viewer Basics)

These pages explain how to *drive the viewer* confidently: navigation, camera, render modes, multiview, and dimension switching.

**Audience:** everyone (novices get click-by-click; power users get exact controls + edge cases)  
**Time:** 15–30 minutes for full coverage

**Recommended reading order**

1) `c_core_interactions/01_navigation_modes_orbit_planar_free_fly`
2) `c_core_interactions/04_view_layout_live_snapshots_small_multiples`
3) `c_core_interactions/05_dimension_switching_1d_2d_3d`
4) `c_core_interactions/03_render_modes_points_vs_volumetric_smoke`
5) `c_core_interactions/02_camera_controls_advanced` (when you want to tune)
6) `c_core_interactions/06_troubleshooting_core_interactions` (when something feels “stuck”)

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

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
