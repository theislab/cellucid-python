# Python troubleshooting index

Start with the symptom you can observe. Preserve the first exception, HTTP
status, or browser-console error: it usually identifies the owning layer.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Installation and dependencies
:link: 01_installation_and_dependency_issues
:link-type: doc

Import errors, interpreter mismatches, and installation checks.
:::

:::{grid-item-card} Data preparation
:link: 02_data_preparation_issues
:link-type: doc

Shape, dtype, identity, embedding, metadata, expression, graph, and export
validation.
:::

:::{grid-item-card} Viewer embedding
:link: 03_viewer_embedding_issues
:link-type: doc

Notebook display, iframe reachability, asset delivery, and browser errors.
:::

:::{grid-item-card} Server mode
:link: 04_server_mode_issues
:link-type: doc

Bind failures, URLs, HTTP responses, remote access, and process lifecycle.
:::

:::{grid-item-card} Hooks and events
:link: 05_hooks_and_events_issues
:link-type: doc

Viewer IDs, callbacks, event validation, and session-bundle requests.
:::

:::{grid-item-card} Performance and memory
:link: 06_performance_and_memory_issues
:link-type: doc

Separate export time, I/O, browser rendering, and analysis costs.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
