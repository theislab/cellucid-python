# Developer Documentation (Web App Architecture)

This section is a **developer-facing map of how the Cellucid web app works**: where the important code lives, how state flows through the system, and how to debug/extend features without breaking performance.

It lives in the `cellucid-python` docs because most developers arrive from the Python side, but **the frontend code is in the `cellucid/` repo**.

:::{important}
If you are a wet-lab user or a general “how do I use the UI?” reader, this section is *not* the right starting point.
Start with the web app user guide landing page: {doc}`../index`.
:::

---

## Fast Path (Pick Your Goal)

| You want to… | Start here | Then read |
|---|---|---|
| Run the app locally and reproduce a bug | {doc}`02_local_development_setup` | {doc}`13_debugging_playbook` |
| Understand the big-picture architecture | {doc}`05_app_architecture_overview` | {doc}`06_state_datastate_and_events` → {doc}`07_rendering_pipeline_webgl_and_performance_notes` |
| Learn where the code lives (file map + entry points) | {doc}`01_codebase_map_and_entry_points` | {doc}`08_ui_modules_map` |
| Add a new sidebar feature/module | {doc}`16_extension_point_add_ui_module` | {doc}`15_extension_points_overview` |
| Add/extend analysis (plots, compute, UI modes) | {doc}`11_analysis_architecture` | {doc}`17_extension_point_add_analysis_mode` |
| Work on sessions (save/load) | {doc}`10_sessions_persistence_and_serialization` | {doc}`13_debugging_playbook` |
| Work on figure export | {doc}`12_figure_export_architecture` | {doc}`18_extension_point_add_export_renderer` |
| Understand loading sources (local/remote/GitHub/Jupyter) | {doc}`09_data_loading_pipeline_and_caching` | {doc}`../b_data_loading/index` (user-facing view) |

---

## What This Section Covers

- **Codebase map**: how `cellucid/index.html` boots the app and which folders correspond to state, UI, rendering, and data loading.
- **State model**: `DataState`, events, per-view contexts, multi-page highlights, caches, and batch updates.
- **Rendering pipeline**: the WebGL viewer, buffers, LOD, overlays (smoke + velocity), and performance footguns.
- **Persistence**: the `.cellucid-session` bundle format, what is saved/restored, and dataset mismatch behavior.
- **Debugging**: how to capture minimal repro steps + the exact logs that help fix issues quickly.
- **Extension points**: the safe, repeatable steps for adding a UI module, an analysis mode, or an export renderer.

:::{note}
Privacy/security and telemetry live in a separate user-facing section:
{doc}`../o_accessibility_privacy_security/index`.
This developer section links to it, but does not replace it.
:::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
