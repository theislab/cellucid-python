# Contributing

Cellucid is an ecosystem of repositories. Choose where you want to contribute:

::::{grid} 1 1 3 3
:gutter: 3

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web App (`cellucid`)
:link: contributing/web_app
:link-type: doc

UI workflows, rendering/performance, sessions, figure export, and community annotation frontend.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Python Package + CLI (`cellucid-python`)
:link: contributing/python
:link-type: doc

`prepare()` and the Cellucid-ready AnnData, `serve` backend, Jupyter integration, Python APIs, and this documentation site.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` R Package (`cellucid-r`)
:link: contributing/r
:link-type: doc

R writer (`cellucid_prepare()`), R-side recipes (Seurat/SCE), and R testing + docs.
:::

::::

## Quick routing guide

- **UI bug / rendering / figure export / sessions** → Web App
- **Python `prepare()` / CLI / server / notebooks / docs** → Python package
- **R `cellucid_prepare()` issues / Seurat or SCE recipes** → R package

```{toctree}
:maxdepth: 1
:hidden:

contributing/web_app
contributing/python
contributing/r
```

