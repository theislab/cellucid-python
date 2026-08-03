# Cellucid

**See every cell. Keep the biology in view.**

Cellucid is a GPU-accelerated, browser-first workspace for exploring
single-cell data interactively. Move through 1D, 2D, and 3D embeddings; color
by genes or cell metadata; filter and compare populations; inspect the neighbor
graph and velocity overlays; and carry the result into a figure, notebook,
saved session, or collaborative annotation round. Every term the app uses is
defined in the {doc}`UI glossary
<user_guide/web_app/a_orientation/04_ui_glossary_terminology>`.

::::{grid} 1 2 2 2
:gutter: 2

:::{grid-item}
```{button-link} https://www.cellucid.com
:color: primary
:expand:

{octicon}`browser;1em` Open Cellucid
```
:::

:::{grid-item}
```{button-ref} user_guide/web_app/a_orientation/03_quick_tour_60_seconds
:ref-type: doc
:color: secondary
:expand:

🚀 Take the 60-second tour
```
:::

::::

```{figure} _static/screenshots/web_app/suo-cell-type-close-up.jpg
:alt: A close-up of the Suo dataset in Cellucid, showing dense groups of categorically colored cells at varied point sizes.
:width: 70%
:align: center

The Suo atlas makes Cellucid's scale tangible: hundreds of thousands of cells
remain individually visible as you move from tissue-wide structure into a
close-up view.
```

## Start with the question in front of you

You do not need to learn the whole application before doing useful work. Pick
the route that matches today’s question.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Where are my populations?
:link: user_guide/web_app/d_fields_coloring_legends/index
:link-type: doc

Colour by a categorical `obs` field or a marker gene, then
{doc}`filter the noise away <user_guide/web_app/e_filtering/index>` and
{doc}`highlight a population <user_guide/web_app/f_highlighting_selection/index>`
without losing the surrounding context.
:::

:::{grid-item-card} Where are cells going?
:link: user_guide/web_app/b_data_loading/10_standard_pancreas_dataset
:link-type: doc

Open the real Pancreas sample, compare its 1D, 2D, and 3D embeddings, and
inspect the dimension-matched RNA-velocity overlay.
:::

:::{grid-item-card} How do I share my own data?
:link: user_guide/web_app/b_data_loading/11_custom_dataset_repository
:link-type: doc

Build a small public catalog with stable dataset identities, test it locally,
and share an exact URL that opens the intended dataset.
:::

:::{grid-item-card} Can the viewer stay beside my analysis?
:link: user_guide/web_app/b_data_loading/05_jupyter_tutorial
:link-type: doc

Embed Cellucid in Jupyter, inspect connection health, exchange selections and
state with Python, and stop the local viewer cleanly.
:::

::::

## Bring data in the way that fits the work

| What you have | Best first path | Why |
|---|---|---|
| A prepared Cellucid folder | {doc}`Browser Prepared picker <user_guide/web_app/b_data_loading/03_browser_file_picker_tutorial>` | Fast local exploration with no separate server |
| An H5AD under 512 MiB, or a portable Zarr ZIP | {doc}`Browser file picker <user_guide/web_app/b_data_loading/03_browser_file_picker_tutorial>` | One direct local selection, no Python installed |
| A large H5AD, Zarr directory, or prepared catalog | {doc}`Python server <user_guide/web_app/b_data_loading/04_server_tutorial>` | Exact server-side validation and on-demand browser gene requests |
| An in-memory `AnnData` object | {doc}`Jupyter integration <user_guide/web_app/b_data_loading/05_jupyter_tutorial>` | A tight analysis-to-visualization loop |
| A public collection of prepared datasets | {doc}`Custom repository guide <user_guide/web_app/b_data_loading/11_custom_dataset_repository>` | No running server; exact catalog and share links |
| A static host or CDN you already control | {doc}`Self-hosted exports root <user_guide/web_app/b_data_loading/11_custom_dataset_repository>` | `?exportsBaseUrl=` fills the app's own `Sample datasets:` dropdown from your catalog |

The {doc}`complete data-loading map
<user_guide/web_app/b_data_loading/01_loading_options_overview>` explains the
memory, network, privacy, and reproducibility tradeoffs before you commit to a
path.

## Learn at your own depth

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web app
:link: user_guide/web_app/index
:link-type: doc

Warm, task-focused guides for navigation, coloring, filtering, highlighting,
analysis, overlays, multiview, sessions, figures, accessibility, and
troubleshooting.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Python
:link: user_guide/python_package/index
:link-type: doc

Exact preparation, server, Jupyter, hooks, session, performance, and API
contracts—from a first notebook to a reproducible deployment.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` R
:link: user_guide/r_package/index
:link-type: doc

Prepare the same browser-ready format from matrices, Seurat, or
SingleCellExperiment, then validate and share it through the same viewer.
:::

::::

## One ecosystem, clear handoffs

- [cellucid](https://github.com/theislab/cellucid) is the web application.
- [cellucid-python](https://github.com/theislab/cellucid-python) prepares,
  serves, and embeds data.
- [cellucid-r](https://github.com/theislab/cellucid-r) prepares the same export
  contract from R.
- [cellucid-datasets](https://github.com/theislab/cellucid-datasets) hosts the
  public sample catalog, including the standard Pancreas dataset.
- [cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets)
  is a tiny, inspectable model for your own public catalog.
- [cellucid-annotation](https://github.com/theislab/cellucid-annotation)
  provides the repository contract for community annotation.

Install the Python tools with `pip install cellucid`. For R installation and
supported inputs, begin with the {doc}`R guide
<user_guide/r_package/installation>`.

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
