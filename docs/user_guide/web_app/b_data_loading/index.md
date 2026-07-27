# Data Loading in the Web App (All Paths)

These pages cover the end-to-end ways to get data into Cellucid (pre-exported folders, browser file picker, server mode, and Jupyter).

They are written for **mixed audiences**:
- **Wet lab / non-technical**: click-by-click, “what success looks like”, and safe defaults.
- **Computational users**: formats, performance tradeoffs, parameter choices, and reproducibility.
- **Power users / developers**: edge cases, failure modes, and how to debug loading issues.

:::{important}
If you are not sure which method to use, start with {doc}`01_loading_options_overview`.
:::

:::{note}
Vector fields (velocity/drift overlays): if your dataset includes per-cell vectors, Cellucid can visualize them as an animated overlay after loading.
You don’t need a special loading method—just make sure the vector field data is present and named correctly.
See {doc}`../i_vector_field_velocity/index` and {doc}`07_folder_file_format_expectations_high_level_link_to_spec`.
:::

---

## Fast Path (Choose Your Workflow)

Use this as a decision tree. You can always switch later.

| You have… | Best first choice | Why | Next page |
|---|---|---|---|
| A pre-exported folder from `prepare()` | Browser folder picker | Fastest, most reliable, no server | `03_browser_file_picker_tutorial` |
| A portable `.zarr.zip` or `.zip` | Browser Zarr ZIP picker | One validated file selection in every supported browser | `03_browser_file_picker_tutorial` |
| A `.h5ad` file | Server mode (recommended) or Jupyter | Python opens the file read-only-backed, reducing matrix memory pressure | `04_server_tutorial` or `05_jupyter_tutorial` |
| A `.zarr` directory | Server mode or Jupyter | Direct Python loading is supported and eager | `04_server_tutorial` or `05_jupyter_tutorial` |
| An in-memory `AnnData` in a notebook | Jupyter | Fastest way to iterate while analyzing | `05_jupyter_tutorial` |
| A dataset you want to share publicly | GitHub-hosted exports | Shareable URL, no server needed | `02_local_demo_tutorial` |

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 246px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

:::{note}
The Session panel keeps every loading action and connection state visible.
Its compact `i` buttons only open explanatory help: click one, or focus it with
{kbd}`Tab` and press {kbd}`Enter` or {kbd}`Space`. Press {kbd}`Escape`, click
elsewhere, or move focus outside to close it. See
{doc}`../a_orientation/04_ui_glossary_terminology` and
{doc}`../o_accessibility_privacy_security/01_accessibility`.
:::

---

## Getting Started

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Loading Options Overview
:link: 01_loading_options_overview
:link-type: doc

Understand all 14 loading options and choose the best approach for your workflow: in-memory, h5ad, zarr, or pre-exported.
:::

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Local Demo
:link: 02_local_demo_tutorial
:link-type: doc

Export once, then share via public GitHub-hosted exports (no server), or run a local demo viewer.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Standard Pancreas Sample
:link: 10_standard_pancreas_dataset
:link-type: doc

Load the scVelo pancreas sample and verify its 1D/2D/3D embeddings,
dimension-matched velocity, metadata, genes, connectivity, and provenance.
:::

::::

---

## Viewing Methods

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`upload;1.5em;sd-mr-1` Browser File Picker
:link: 03_browser_file_picker_tutorial
:link-type: doc

Load a prepared directory, one H5AD file, or one portable Zarr ZIP through the
three visible local-data controls.
:::

:::{grid-item-card} {octicon}`server;1.5em;sd-mr-1` Server Mode
:link: 04_server_tutorial
:link-type: doc

Run a local HTTP server for larger datasets with efficient lazy loading of gene expression.
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` Jupyter Integration
:link: 05_jupyter_tutorial
:link-type: doc

Visualize AnnData objects directly within Jupyter notebooks with embedded interactive widgets.
:::

::::

---

## Concepts & Troubleshooting

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`tag;1.5em;sd-mr-1` Dataset Identity
:link: 06_dataset_identity_why_it_matters
:link-type: doc

What makes a dataset “the same dataset” in Cellucid, and why it matters for sessions, sharing, and community annotation.
:::

:::{grid-item-card} {octicon}`file;1.5em;sd-mr-1` Format Expectations
:link: 07_folder_file_format_expectations_high_level_link_to_spec
:link-type: doc

What files/keys are required for exports, GitHub manifests, `.h5ad`, and `.zarr`.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting
:link: 08_troubleshooting_data_loading
:link-type: doc

Symptom → diagnosis → fix for file picker, server mode, GitHub exports, and common data issues.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Verified captures
:link: 09_screenshots
:link-type: doc

Current loading controls, a successful direct H5AD load, and browser-engine
acceptance captures.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Vector Fields Overlay
:link: ../i_vector_field_velocity/index
:link-type: doc

If your dataset includes vector fields (e.g. RNA velocity), enable the overlay after loading and pick a field for the current dimension.
:::

::::

---

## Quick Comparison

| Method | Best For | Requires Server |
|--------|----------|-----------------|
| Local Demo | Quick visualization, sharing static files | No |
| File Picker | Desktop viewing without setup | No |
| Server Mode | Large datasets, lazy gene loading | Yes |
| Jupyter | Interactive analysis workflows | No |

:::{note}
“Requires Server” here means “you run a small local Python process that serves data to the web viewer”.
This is still typically **local-only** (e.g. `127.0.0.1`) unless you explicitly bind to `0.0.0.0`.
:::

---

:::{tip}
Adding more pages:

- Add new pages into `cellucid-python/docs/user_guide/web_app/b_data_loading/`.
- Use a numeric prefix like `10_...` so they naturally sort.
- This page includes them automatically via a globbed toctree.
:::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
