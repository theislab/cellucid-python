# Screenshots

This page is a **screenshot capture checklist** for the Core Interactions section.

It exists so you (or a collaborator) can capture screenshots once, systematically, without hunting through every page.

---

## How placeholders work

All pages in this section currently use:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Each placeholder is preceded by an HTML comment with:

- what state to capture,
- what to crop/redact,
- suggested caption/alt text,
- suggested filename conventions.

General guidance lives in:

- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## Recommended screenshot set (core interactions)

### Navigation modes (3 screenshots)

1) Orbit mode (with orbit anchor visible)
2) Planar mode (with zoom-to-cursor toggle visible)
3) Free-fly mode (with Capture pointer visible)

You can capture these from:

- `c_core_interactions/01_navigation_modes_orbit_planar_free_fly`

### Multiview badges + indicators (1–2 screenshots)

Capture:

- at least 2 snapshots,
- grid compare layout,
- cameras unlocked (so Orb/Pan/Fly badges show),
- a visible dimension badge (e.g. 3D).

Page:

- `c_core_interactions/04_view_layout_live_snapshots_small_multiples`

### Dimension switching (2 screenshots)

1) Dimension dropdown open showing available dims
2) Before/after switching dimensions (same view)

Page:

- `c_core_interactions/05_dimension_switching_1d_2d_3d`

### Render modes (2 screenshots)

1) Points mode (with render mode dropdown visible)
2) Smoke mode (with smoke controls visible)

Page:

- `c_core_interactions/03_render_modes_points_vs_volumetric_smoke`

---

## Highly recommended “failure mode” screenshots (debugging gold)

These save huge amounts of time in support and onboarding:

- WebGL context lost overlay (requires reload)
- Smoke mode blank/slow (with grid + ray quality visible)
- Pointer lock active (cursor hidden) and how to exit (`Esc`)

System requirements page includes a placeholder for the context-lost overlay:

- `a_orientation/02_system_requirements`
