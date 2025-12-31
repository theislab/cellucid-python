# Cellucid

**Interactive Single-Cell Data Visualization**

Cellucid is a **GPU-accelerated**, browser-first single-cell viewer for exploring massive datasets in real time — fly through 2D/3D embeddings, color by genes or metadata, filter and compare populations, and share reproducible annotations with collaborators.

::::{grid} 1 1 1 1
:gutter: 2

:::{grid-item}
```{button-link} https://cellucid.com
:color: primary
:expand:

{octicon}`browser;1em` Cellucid App
```
:::

::::

---

## Features

- **GPU-accelerated 2D/3D rendering** - Smooth navigation through millions of cells
- **Genes + metadata overlays** - Fast coloring for sparse/dense expression and `obs` fields
- **Filtering + selection workflows** - Highlight, subset, and compare populations interactively
- **Connectivity + vector overlays** - KNN edges and velocity/drift-style fields
- **Shareable exports** - Static folders you can host locally, on GitHub, or behind a server
- **Works everywhere** - Web app, local server mode, and Jupyter notebook embedding

---

## User Guide

Step-by-step tutorials to get you started with Cellucid visualization.

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web App
:link: user_guide/web_app/index
:link-type: doc

Viewer basics, filtering, highlighting, analysis, figure export, sessions, performance, and troubleshooting.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Python Package
:link: user_guide/python_package/index
:link-type: doc

Installation, data prep, `prepare()`/export, viewing modes, notebook integration, hooks/events, and troubleshooting.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` R Export (Preview)
:link: user_guide/r_package/index
:link-type: doc

Exporting data with `cellucid-r` (Seurat/SCE recipes) into the same format used by the Cellucid web app.
:::

::::

## Installation

Python:

```bash
pip install cellucid
```

R (export package):

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

```{toctree}
:maxdepth: 2
:hidden:
:caption: Sections

Web app <user_guide/web_app/index>
Python <user_guide/python_package/index>
R <user_guide/r_package/index>
contributing
changelog
```
