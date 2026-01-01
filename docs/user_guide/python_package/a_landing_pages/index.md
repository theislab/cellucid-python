# Landing Pages

**Audience:** everyone (choose your depth)  
**Time:** 10–30 minutes  
**Goal:** pick the right workflow and get your first Cellucid view running

This section is intentionally written for mixed audiences:
- **Wet lab / non-technical users**: step-by-step “do this, then that” instructions and clear success criteria.
- **Computational users**: explicit shapes/dtypes, configuration options, performance considerations, and edge cases.
- **Developers/maintainers**: architecture, file formats, and reproducibility notes.

Cellucid is a **web app** (the viewer UI) plus **helper packages** that bring your data into the viewer:
- **Cellucid web app**: what you see and interact with in the browser.
- **`cellucid` (cellucid-python)**: this package — export/serve data and embed the web app in notebooks.
- **`cellucid-annotation`**: helper repo for community annotation workflows.
- **`cellucid-r`**: planned; not ready yet.

## Choose your workflow (start here)

- “I have an AnnData and want to see it *now*” → {doc}`04_quick_start_3_levels` (Level 1: `show_anndata`)
- “I need a shareable, fast export folder (papers/collaboration)” → {doc}`04_quick_start_3_levels` (Level 2: `prepare`) + {doc}`../c_data_preparation_api/index`
- “I want the viewer in a browser (no notebook)” → {doc}`02_installation` (verify CLI) + {doc}`04_quick_start_3_levels` (use `cellucid serve`)
- “I need Python ↔ UI hooks (selection callbacks, programmatic highlights)” → {doc}`04_quick_start_3_levels` (Level 3) + {doc}`../e_jupyter_hooks/index`
- “Something failed” → {doc}`../i_troubleshooting_index/index`

```{note}
If you prefer notebook-style, highly verbose walkthroughs (with lots of edge cases and troubleshooting), start with {doc}`../f_notebooks_tutorials/index`.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`info;1.5em;sd-mr-1` What is cellucid-python?
:link: 01_what_is_cellucid_python
:link-type: doc

What the Python package does, how it relates to the Cellucid web app, and which workflows it supports.
:::

:::{grid-item-card} {octicon}`download;1.5em;sd-mr-1` Installation
:link: 02_installation
:link-type: doc

How to install, optional dependencies, platform notes, and installation troubleshooting.
:::

:::{grid-item-card} {octicon}`zap;1.5em;sd-mr-1` Quick start (3 levels)
:link: 04_quick_start_3_levels
:link-type: doc

Copy/paste: a minimal “show”, a practical export workflow, and an advanced server + hooks workflow.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Compatibility matrix (must be explicit)
:link: 03_compatibility_matrix_must_be_explicit
:link-type: doc

Exactly what works (and what doesn’t) across Jupyter, JupyterLab, VSCode notebooks, and Colab.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
