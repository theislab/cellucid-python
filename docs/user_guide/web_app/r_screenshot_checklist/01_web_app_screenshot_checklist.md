# Web App Screenshot Checklist (Docs Authoring)

This checklist helps you capture screenshots for the **Cellucid web app** docs in a consistent, reproducible way.

For global screenshot rules (captions, alt text, redaction), also see `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`.

---

## Folder + naming convention

Recommended location:

- `cellucid-python/docs/_static/screenshots/web_app/`

Recommended naming:

- prefix with ordering within a page: `01_...`, `02_...`
- include the feature + state: `core-interactions-orbit-mode.png`, `filters-all-cells-filtered-out.png`

---

## What every screenshot must include

Use this as a “pre-flight” check before committing images:

- **Alt text**: one sentence describing the key UI element/state (not “Screenshot of…”).
- **Caption**: one sentence telling the reader what to learn/confirm.
- **Inline capture note**: an HTML comment above the `{figure}` block describing how to reproduce it later.
- **Redaction**: remove dataset ids/names, repo/org names, usernames, sample ids, tokens, local paths (when sensitive).

---

## Minimal coverage per page (recommended)

For most web app pages, aim for:

1) **Orientation screenshot**: where this feature lives in the UI (sidebar accordion/panel visible).
2) **Success-state screenshot**: what “worked” looks like after the steps.
3) **Common failure screenshot** (optional but high value): the exact error/empty state users report.

---

## Coverage by web-app section

This is a checklist of screenshot *types* (not exact filenames) you should plan to capture.

### A) Orientation (`a_orientation/`)

- One full-screen UI “map” screenshot with labeled callouts (sidebar sections).
- System requirements: one example of a “WebGL2 required” message (or browser capability warning).

### B) Data loading (`b_data_loading/`)

- Loader panel empty state (no dataset loaded).
- Successful dataset loaded state showing dataset name + counts.
- One failing load with the visible on-screen error text (CORS / missing files / invalid export).

### C) Core interactions (`c_core_interactions/`)

- Orbit vs planar vs free-fly: one annotated screenshot per mode, showing the mode badge/control.
- Multiview: “live view only” vs “after Keep view” grid state with view badges visible.
- Dimension switching: before/after 1D→2D→3D (show the dimension badge/selector).
- WebGL context lost overlay (if you can capture it safely).

### D) Fields, coloring, legends (`d_fields_coloring_legends/`)

- Field dropdown open showing categorical vs continuous grouping.
- One categorical legend example (with category toggles if available).
- One continuous legend example (with min/max/clipping controls if exposed).

### E) Filtering (`e_filtering/`)

- Active filters panel showing a small filter stack (enabled + disabled items).
- “All cells filtered out” empty state (and the UI’s recovery affordance).

### F) Highlighting/selection (`f_highlighting_selection/`)

- Lasso selection in progress (cursor/selection path visible).
- Resulting highlight group list with color and counts visible.

### (Optional) Troubleshooting pages

- For each major module, capture the single most common failure state you see reported (and its on-screen message).

