# API Reference

This page provides an overview of all public Cellucid objects, functions, and methods.

---

## Quick Reference

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {octicon}`play;1.5em;sd-mr-1` Jupyter
:link: jupyter
:link-type: doc

Interactive visualization in notebooks
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Server
:link: server
:link-type: doc

HTTP server for local viewing
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Export
:link: export
:link-type: doc

Prepare data for static hosting
:::

::::

---

## Functions

### Jupyter Display

Display visualizations directly in Jupyter notebooks.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.show_anndata
   cellucid.show
```

### Server

Run a local HTTP server to view data in browser.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.serve_anndata
   cellucid.serve
```

### Data Preparation

Export AnnData for static web hosting or sharing.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.prepare
```

---

## Classes

### Viewers

High-level classes for Jupyter integration.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.AnnDataViewer
   cellucid.CellucidViewer
```

### Servers

Server classes for fine-grained control.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.AnnDataServer
   cellucid.CellucidServer
```

### Adapters

Data adapters for AnnData handling.

```{eval-rst}
.. autosummary::
   :nosignatures:

   cellucid.AnnDataAdapter
```

---

```{toctree}
:maxdepth: 2
:hidden:

jupyter
server
export
viewers
adapters
```
