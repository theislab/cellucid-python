# Data Loading in the Web App (All Paths)

These pages cover every current route into Cellucid: public samples, prepared
folders, browser-selected H5AD or Zarr ZIP files, local or remote Python
servers, Jupyter, and public GitHub catalogs.

They are written for **mixed audiences**:
- **Wet lab / non-technical**: click-by-click, “what success looks like”, and safe defaults.
- **Computational users**: formats, performance tradeoffs, parameter choices, and reproducibility.
- **Power users / developers**: edge cases, failure modes, and how to debug loading issues.

:::{important}
If you are not sure which method to use, start with {doc}`01_loading_options_overview`.
:::

:::{note}
Vector fields (velocity/drift overlays): if your dataset includes per-cell
vectors, Cellucid can draw them as an animated overlay after loading. Two
things decide whether the overlay appears, and only the first is about naming:

- **A prepared export** carries whatever `prepare(vector_fields=...)` wrote, and
  the **three browser pickers** discover `<field>_umap_<dim>d` keys in `obsm` by
  themselves. Nothing to switch on.
- **Every Python path** — `cellucid serve`, `serve_anndata()`, `show_anndata()` —
  leaves the vectors unread unless you ask, exactly like the neighbor graph.
  Pass `--vector-fields` or `serve_vector_fields=True`.

See {doc}`../i_vector_field_velocity/index` and {doc}`07_folder_file_format_expectations_high_level_link_to_spec`.
:::

---

## Fast Path (Choose Your Workflow)

Use this as a decision tree. You can always switch later.

| You have… | Best first choice | Why | Next page |
|---|---|---|---|
| A pre-exported folder from `prepare()` | Browser Prepared picker | Fastest, most reliable, no server | {doc}`03_browser_file_picker_tutorial` |
| A portable `.zarr.zip` or `.zip` | Browser Zarr ZIP picker | One validated file selection in every supported browser | {doc}`03_browser_file_picker_tutorial` |
| A `.h5ad` **under 512 MiB**, and no Python | Browser H5AD picker | One file selection, nothing installed, nothing uploaded | {doc}`03_browser_file_picker_tutorial` |
| A `.h5ad` **over 512 MiB**, or one that stalls the tab | Server mode (recommended) or Jupyter | Python opens the file read-only-backed, reducing matrix memory pressure | {doc}`04_server_tutorial` or {doc}`05_jupyter_tutorial` |
| A `.zarr` directory | Server mode or Jupyter | Direct Python loading is supported and eager | {doc}`04_server_tutorial` or {doc}`05_jupyter_tutorial` |
| An in-memory `AnnData` in a notebook | Jupyter | Fastest way to iterate while analyzing | {doc}`05_jupyter_tutorial` |
| A dataset collection you want to share publicly | GitHub-hosted exports | Exact shareable URLs, no running server | {doc}`11_custom_dataset_repository` |
| A static host, CDN, or intranet server you control | Self-hosted exports root | `?exportsBaseUrl=https://host/exports/` fills the app's own `Sample datasets:` dropdown from your catalog | {doc}`11_custom_dataset_repository` |

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Cellucid Session panel with a dataset loaded: a summary table of dataset, description, source, URL, cells, genes, obs fields and connectivity; a Sample datasets dropdown; H5AD, Zarr ZIP and Prepared buttons under Local data; Remote server and GitHub data text fields each with a Connect button; and Save State and Load State buttons.
:width: 524px

Every loading path lives in its own labelled group of the Session panel, and
Save State and Load State stay beside the dataset controls. The summary table at
the top always describes the dataset currently loaded, whichever path loaded it.
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

:::{grid-item-card} {octicon}`repo;1.5em;sd-mr-1` Custom Dataset Repository
:link: 11_custom_dataset_repository
:link-type: doc

Follow a complete three-dataset reference repository, validate its exact
catalog contract, and publish stable dataset-specific links.
:::

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Standard Pancreas Sample
:link: 10_standard_pancreas_dataset
:link-type: doc

Load the scVelo pancreas sample and verify its 1D/2D/3D embeddings,
dimension-matched velocity, metadata, genes, connectivity, and provenance.
:::

:::{grid-item-card} {octicon}`verified;1.5em;sd-mr-1` Sample Provenance
:link: 12_sample_dataset_provenance
:link-type: doc

What to cite for each built-in sample, which one can be rebuilt and checked
from its recipe, and what is recorded for the four that cannot.
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

| Method | Best For | Separate server command? |
|--------|----------|--------------------------|
| Public sample | Learning with known-good data | No |
| Browser file picker | Local viewing without Python | No |
| GitHub catalog | Public, dataset-specific share links | No |
| Self-hosted exports root | A catalog on your own static host or CDN | No — any HTTP server will do |
| Server mode | Large data and on-demand browser gene requests | Yes — `cellucid serve` |
| Jupyter | Interactive analysis workflows | No — the viewer manages a localhost server |

:::{note}
Server mode and Jupyter both use HTTP internally. The distinction above is
whether you start a separate terminal command. Both bind to `127.0.0.1` by
default; use the documented SSH-tunnel workflow for a remote machine.
:::

---

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
