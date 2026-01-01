# Advanced notebooks (expert / developer)

These tutorials are for power users and developers who want to:
- treat the viewer as an interactive front-end for analysis pipelines
- round-trip state via session bundles (viewer → Python → AnnData)
- integrate custom messages/events
- validate and debug vector fields (velocity/drift overlays)

## Suggested order

1) {doc}`32_session_persistence_and_restoring_analysis_artifacts`
2) {doc}`23_programmatic_highlighting_and_selection_callbacks`
3) {doc}`33_vector_fields_and_velocity_overlay_end_to_end`
4) {doc}`31_custom_message_schemas_and_frontend_extensions`

## What you will learn

- Session bundle capture/restore and how it depends on dataset identity
- Robust event handling patterns (timeouts, deduping, multiple viewers)
- Vector field invariants (shape, dimension, scale) and failure modes
- Debugging the Python server + iframe + web UI cache stack

```{note}
For runnable hook/session notebooks, see {doc}`05_jupyter_embedding_hooks_sessions_gallery`.
```

```{toctree}
:maxdepth: 1
:hidden:

31_custom_message_schemas_and_frontend_extensions
32_session_persistence_and_restoring_analysis_artifacts
33_vector_fields_and_velocity_overlay_end_to_end
```
