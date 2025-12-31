# Screenshots

This page is a **production checklist** for screenshots needed by the vector field / velocity overlay documentation.

Use the screenshot playbook:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

:::{important}
All screenshots here are placeholders/specs for now. Replace the `{figure}` targets later with real images under `cellucid-python/docs/_static/screenshots/`.
:::

---

## Minimum required set (recommended)

If you only capture a small set at first, capture these:

1) Overlay controls location (where the block is in the sidebar)
2) Overlay panel expanded (toggle + dropdown + core sliders)
3) Off vs on (same camera)
4) Advanced settings example (show one of the “cinematic” knobs)
5) A common failure state (missing fields / disabled toggle / error toast)

---

## Overlay panel (orientation)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-panel-expanded
Suggested filename: vector_field_velocity/01_overlay-panel-expanded.png
Where it appears: User Guide → Web App → Vector Field / Velocity Overlay
Capture:
  - UI location: left sidebar → overlay module/panel (expanded)
  - State prerequisites: dataset loaded with at least one 2D vector field
  - Action to reach state: enable overlay, open advanced settings (if applicable)
Crop:
  - Include: overlay enable toggle + field dropdown + 2–3 core sliders
  - Exclude: personal bookmarks, private dataset names if sensitive
Annotations:
  - Callouts: toggle, dropdown, density slider
Alt text:
  - Velocity overlay settings panel expanded with enable toggle and field selector.
Caption:
  - Where to enable the overlay and select a vector field.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the expanded velocity overlay panel.
:width: 100%

Where to enable the overlay and select a vector field.
```

---

## UI location (where this feature lives)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-ui-location
Suggested filename: vector_field_velocity/00_ui-location.png
Where it appears: Vector Field / Velocity Overlay → enabling overlay
Capture:
  - UI location: left sidebar → “Visualization”
  - State prerequisites: dataset loaded; Render mode set to “Points”
  - Action to reach state:
    - Open “Visualization”
    - Ensure “Render mode: Points” is visible
    - Scroll so “Vector Field Overlay:” block is visible in-frame
Crop:
  - Include: “Visualization” header, “Render mode” control, and the “Vector Field Overlay:” block
  - Exclude: personal bookmarks, private dataset names if sensitive
Annotations:
  - Callouts: Render mode selector, Vector Field Overlay block title
Alt text:
  - Visualization section of the sidebar showing where the Vector Field Overlay controls are located.
Caption:
  - The Vector Field Overlay controls live under Visualization and are currently available in Points render mode.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Vector Field Overlay UI location.
:width: 100%

The Vector Field Overlay controls live under Visualization and are currently available in Points render mode.
```

---

## Off vs on (success state)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-off-vs-on
Suggested filename: vector_field_velocity/02_off-vs-on.png
Where it appears: Vector Field / Velocity Overlay → enabling overlay
Capture:
  - Same camera, same view; left image overlay off; right image overlay on
  - Use a moderately sparse region so particles are clearly visible
Alt text:
  - Side-by-side comparison of the embedding with the overlay off and on.
Caption:
  - The overlay adds animated particle flow on top of the embedding without changing the underlying points.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for overlay off vs on.
:width: 100%

The overlay adds animated particle flow on top of the embedding without changing the underlying points.
```

---

## Trail settings impact (quality knob)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-trail-settings
Suggested filename: vector_field_velocity/03_trail-settings.png
Where it appears: Vector Field / Velocity Overlay → advanced parameters
Capture:
  - Same field, same camera
  - Show two states: short trails vs long trails
Alt text:
  - Comparison showing how trail length changes the appearance of particle flow.
Caption:
  - Trail settings trade clarity for motion persistence; longer trails can become visually dense on large datasets.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for trail settings comparison.
:width: 100%

Trail settings trade clarity for motion persistence; longer trails can become visually dense on large datasets.
```

---

## Advanced settings expanded (what “Advanced Visual Settings” contains)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-advanced-settings-expanded
Suggested filename: vector_field_velocity/05_advanced-settings-expanded.png
Where it appears: Vector Field / Velocity Overlay → advanced parameters
Capture:
  - UI location: Vector Field Overlay → Advanced Visual Settings (expanded)
  - State prerequisites: overlay enabled on a dataset with vector fields
  - Action to reach state: expand the “Advanced Visual Settings” <details>
Crop:
  - Include: the “Advanced Visual Settings” summary and at least 2 sections (e.g., Particle Rendering + HDR & Bloom)
  - Exclude: unrelated sidebar panels
Annotations:
  - Callouts: Bloom strength slider, Turbulence slider (or any two key knobs)
Alt text:
  - Advanced Visual Settings expanded for the velocity overlay, showing additional sliders like bloom and turbulence.
Caption:
  - Advanced Visual Settings contain post-processing and motion-tuning knobs; use conservatively for scientific figures.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for expanded Advanced Visual Settings.
:width: 100%

Advanced Visual Settings contain post-processing and motion-tuning knobs; use conservatively for scientific figures.
```

---

## Common failure state (dropdown empty / toggle disabled)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-no-fields
Suggested filename: vector_field_velocity/04_no-fields.png
Where it appears: Vector Field / Velocity Overlay → troubleshooting
Capture:
  - Dataset loaded without vector fields (or in wrong dimension)
  - Overlay panel open showing disabled toggle or empty dropdown
Alt text:
  - Velocity overlay panel showing no available fields.
Caption:
  - If no vector field exists for the current dimension, the overlay cannot be enabled.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for no available velocity fields.
:width: 100%

If no vector field exists for the current dimension, the overlay cannot be enabled.
```

---

## Error toast example (failed load)

<!-- SCREENSHOT PLACEHOLDER
ID: velocity-overlay-failed-load-toast
Suggested filename: vector_field_velocity/06_failed-load-toast.png
Where it appears: Vector Field / Velocity Overlay → troubleshooting
Capture:
  - Trigger a realistic failure:
    - load a dataset that advertises vector fields but is missing the file, or
    - select a field/dimension combination that is invalid
  - Ensure the toast/notification is visible and readable
Crop:
  - Include: the toast message and enough context to know it’s the vector overlay
Redact:
  - Remove: URLs, local paths, private dataset ids if sensitive
Alt text:
  - Notification toast showing a failed vector field load for the overlay.
Caption:
  - When the vector field fails to load (missing data or shape mismatch), Cellucid shows an error toast and disables the overlay.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a failed vector field load toast.
:width: 100%

When the vector field fails to load (missing data or shape mismatch), Cellucid shows an error toast and disables the overlay.
```
