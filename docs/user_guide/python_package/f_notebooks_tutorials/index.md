# Notebooks / Tutorials (Very Detailed, Step-by-Step)

These are long-form, notebook-style tutorials for the **Cellucid Python package**. They are intentionally verbose:
- every step is explained,
- edge cases and “gotchas” are called out,
- and each tutorial ends with a large troubleshooting section.

```{note}
These pages are written in **MyST** and live in `docs/user_guide/python_package/f_notebooks_tutorials/`.

Some entries link to real `.ipynb` notebooks in the same folder. The docs site does **not** execute notebooks during the build (`nb_execution_mode = "off"`), so you should run notebooks locally to reproduce results.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Start Here
:link: 00_start_here
:link-type: doc

How to use this section, expected prerequisites, and how to run the notebooks reliably.
:::

:::{grid-item-card} {octicon}`person;1.5em;sd-mr-1` Beginner (Wet Lab Friendly)
:link: 01_beginner_notebooks_wet_lab_friendly
:link-type: doc

“I have an AnnData and I just want to see my cells.” Minimal choices, maximum clarity.
:::

:::{grid-item-card} {octicon}`graph;1.5em;sd-mr-1` Intermediate (Computational)
:link: 02_intermediate_notebooks_computational_workflows
:link-type: doc

Reproducible exports, scaling workflows, and practical “do this when…” guidance.
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Advanced (Expert / Developer)
:link: 03_advanced_notebooks_expert_developer
:link-type: doc

Hooks, session mechanics, vector fields, and format-level debugging.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Real-World Dataset Recipes
:link: 04_real_world_dataset_recipes_gallery
:link-type: doc

Dataset-specific preparation notebooks (end-to-end) you can adapt to your data.
:::

:::{grid-item-card} {octicon}`plug;1.5em;sd-mr-1` Jupyter Embedding + Hooks Sessions
:link: 05_jupyter_embedding_hooks_sessions_gallery
:link-type: doc

Notebooks focused on embedding, Python↔frontend events, and session bridging.
:::

::::

```{toctree}
:maxdepth: 2
:hidden:

00_start_here
01_beginner_notebooks_wet_lab_friendly
02_intermediate_notebooks_computational_workflows
03_advanced_notebooks_expert_developer
04_real_world_dataset_recipes_gallery
05_jupyter_embedding_hooks_sessions_gallery
```
