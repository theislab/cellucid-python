---
orphan: true
---

# User Guide

Step-by-step tutorials to get you started with Cellucid visualization.

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web App
:link: web_app/index
:link-type: doc

Viewer basics, filtering, highlighting, analysis, figure export, sessions, performance, and troubleshooting.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Python Package
:link: python_package/index
:link-type: doc

Installation, data prep, `prepare()`/export, viewing modes, notebook integration, hooks/events, and troubleshooting.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` R Export (Preview)
:link: r_package/index
:link-type: doc

Exporting data with `cellucid-r` (Seurat/SCE recipes) into the same format used by the Cellucid web app.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:

web_app/index
python_package/index
r_package/index
```
