# How to use this section

These pages are written to be **useful to multiple audiences**, even though they are “developer docs”.
You can stop early if you only need the top-level mental model.

## Audience layers (pick your depth)

1) **Fast path (non-technical / wet lab / collaborators)**
   - Goal: understand what’s happening well enough to follow a bug report or trust what “a session” means.
   - Read: {doc}`05_app_architecture_overview` + {doc}`10_sessions_persistence_and_serialization`.

2) **Practical path (computational users / power users)**
   - Goal: run Cellucid locally, understand loading modes, and provide actionable bug reports.
   - Read: {doc}`02_local_development_setup` + {doc}`09_data_loading_pipeline_and_caching` + {doc}`13_debugging_playbook`.

3) **Deep path (contributors / maintainers)**
   - Goal: add features without breaking performance, and understand state/events/rendering boundaries.
   - Read: {doc}`01_codebase_map_and_entry_points` → {doc}`06_state_datastate_and_events` → {doc}`07_rendering_pipeline_webgl_and_performance_notes` → {doc}`08_ui_modules_map` → {doc}`15_extension_points_overview`.

---

## Key repo separation (do not mix up)

Cellucid is intentionally split across repos:

- `cellucid/` (**web app**): the browser UI, state layer, WebGL rendering, figure export, community annotation UI, session bundles.
- `cellucid-python/` (**helper repo**): Python APIs + CLI (`prepare`, `serve`, `show_anndata`), plus these Sphinx docs.
- `cellucid-annotation/` (**helper repo/template**): GitHub repo template for community annotation voting + validation.

:::{important}
When debugging a web app issue, always record whether it happened in:
- the public hosted app (`https://www.cellucid.com`),
- a locally hosted copy of `cellucid/`, or
- the embedded viewer inside Jupyter (Cellucid loaded in an iframe).

These environments share code, but differ in **CORS**, **caching**, and **data-source behavior**.
:::

---

## “Definition of done” for changes

Most bugs in interactive visualization apps come from subtle state interactions.
When you change something in the web app, make sure you explicitly check:

- **State correctness**: do events fire exactly once? do caches invalidate? do counts match?
- **Persistence**: does save/load preserve the relevant settings? if not, is it intentionally excluded?
- **Performance**: did you add per-frame work or per-point allocations outside the renderer?
- **Cross-feature interactions**: filters ↔ highlights ↔ multiview snapshots ↔ analysis ↔ exports ↔ sessions.

If you’re adding a new feature, start by reading:
- {doc}`15_extension_points_overview`
- {doc}`16_extension_point_add_ui_module` (UI)
- {doc}`17_extension_point_add_analysis_mode` (analysis)
- {doc}`18_extension_point_add_export_renderer` (export)

---

## Documentation conventions used in these pages

- **“Code pointers”** are given as file paths (e.g. `cellucid/assets/js/app/main.js`) so you can jump directly to the source.
- **“Symptoms”** describe what a user sees; **“Likely causes”** are ordered; **“How to confirm”** is testable; **“Fix”** is step-by-step.
- “Dev-phase” notes are called out explicitly (things that may change quickly or lack migrations/back-compat).

Next: {doc}`01_codebase_map_and_entry_points`.
