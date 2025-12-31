# Viewing and Loading Exports

This section answers: “I exported a folder — how do I actually view it?”

There are two main workflows:

- **Local viewing (fastest to start):** use the web app’s **browser folder picker** to load a single exported dataset folder.
- **Sharing (recommended):** create an “exports root” with `datasets.json` and host it (often on GitHub) so collaborators can load it easily.

```{note}
The web app has very detailed loading documentation in {doc}`../../web_app/b_data_loading/index`.
This section is the R-focused, export-folder-specific version.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`upload;1.5em;sd-mr-1` Open Exports in the Web App
:link: 01_open_exports_in_cellucid_web_app
:link-type: doc

Click-by-click guide for loading a local export folder (no server needed).
:::

:::{grid-item-card} {octicon}`share;1.5em;sd-mr-1` Host Exports for Sharing
:link: 02_host_exports_for_sharing
:link-type: doc

How to create an exports root + `datasets.json`, and how to share via GitHub or a static host.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Validate Exports & Debug Loading
:link: 03_validate_exports_and_debug_loading
:link-type: doc

Preflight checks (files, JSON, binary sanity) and common browser/CORS failure modes.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
