# Web App Guide

This guide documents the **Cellucid web app UI** end-to-end: how to load data, navigate, filter, analyze, export figures, share sessions, and troubleshoot.

If you are looking for Python-side APIs (`prepare`, `serve`, `show_anndata`, hooks), start with {doc}`../python_package/index` instead.

## Where to start (pick one)

- First time using Cellucid → {doc}`a_orientation/index` then {doc}`b_data_loading/index`
- “I just need to load my dataset” → {doc}`b_data_loading/index`
- “I need to understand filtering + selection” → {doc}`e_filtering/index` and {doc}`f_highlighting_selection/index`
- “I want RNA velocity / vector overlays” → {doc}`i_vector_field_velocity/index`
- “I want community annotation workflows” → {doc}`j_community_annotation/index`
- “Something is broken” → {doc}`q_troubleshooting_index/index`

## Chapters (what each folder contains)

### Getting started

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`info;1.5em;sd-mr-1` Orientation
:link: a_orientation/index
:link-type: doc

What Cellucid is, system requirements, a 60-second tour, and how to choose a workflow.
:::

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Data Loading
:link: b_data_loading/index
:link-type: doc

All supported loading paths (exports, file picker, server mode, Jupyter), plus dataset identity, format expectations, and troubleshooting (includes vector fields where relevant).
:::

::::

### Core UI workflows

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`zoom-in;1.5em;sd-mr-1` Core Interactions
:link: c_core_interactions/index
:link-type: doc

Navigation, selection primitives, and the “mental model” of interacting with points in the viewer.
:::

:::{grid-item-card} {octicon}`paintbrush;1.5em;sd-mr-1` Fields, Coloring, Legends
:link: d_fields_coloring_legends/index
:link-type: doc

How Cellucid treats `obs` fields, categorical vs continuous coloring, legends, and common field pitfalls.
:::

:::{grid-item-card} {octicon}`filter;1.5em;sd-mr-1` Filtering
:link: e_filtering/index
:link-type: doc

Filtering mental model, UI controls, and edge cases (including “why did my selection disappear?”).
:::

:::{grid-item-card} {octicon}`issue-opened;1.5em;sd-mr-1` Highlighting & Selection
:link: f_highlighting_selection/index
:link-type: doc

Selecting points, highlighting, synchronization rules, and how to reason about “active vs visible vs selected”.
:::

:::{grid-item-card} {octicon}`share;1.5em;sd-mr-1` Cross-highlighting
:link: g_cross_highlighting/index
:link-type: doc

How cross-highlighting works across viewers/embeddings and how to debug mismatches.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Analysis
:link: h_analysis/index
:link-type: doc

Analysis panel workflows (e.g. marker genes / differential expression), interpretation guidance, and troubleshooting.
:::

::::

### Specialized features

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`pulse;1.5em;sd-mr-1` Vector Fields / Velocity Overlay
:link: i_vector_field_velocity/index
:link-type: doc

Enable overlays, choose fields per dimension, tune parameters, and troubleshoot missing/incorrect velocity visuals.
:::

:::{grid-item-card} {octicon}`people;1.5em;sd-mr-1` Community Annotation
:link: j_community_annotation/index
:link-type: doc

How multi-user annotation voting works, author setup, annotator workflows, and UI reference.
:::

:::{grid-item-card} {octicon}`image;1.5em;sd-mr-1` Figure Export
:link: k_figure_export/index
:link-type: doc

Export UI walkthrough, formats, quality knobs, metadata/provenance, and edge cases for publication-ready figures.
:::

:::{grid-item-card} {octicon}`link;1.5em;sd-mr-1` Sessions & Sharing
:link: l_sessions_sharing/index
:link-type: doc

Saving/restoring UI state, share links, collaboration patterns, and “what makes a session compatible”.
:::

::::

### Performance, safety, and maintenance

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`meter;1.5em;sd-mr-1` Benchmarking & Performance
:link: n_benchmarking_performance/index
:link-type: doc

Performance mental model, best practices for large datasets, and how to diagnose GPU/rendering bottlenecks.
:::

:::{grid-item-card} {octicon}`shield;1.5em;sd-mr-1` Accessibility / Privacy / Security
:link: o_accessibility_privacy_security/index
:link-type: doc

Accessibility guidance plus a practical privacy/security model for “browser ↔ server ↔ data”.
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` Developer Docs
:link: p_developer_docs/index
:link-type: doc

Architecture notes, state/event model, debugging playbooks, and extension points for contributors.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Troubleshooting Index
:link: q_troubleshooting_index/index
:link-type: doc

Symptom → diagnosis → fix across installation, data loading, rendering, selection, analysis, and export.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Screenshot Checklist
:link: r_screenshot_checklist/index
:link-type: doc

One place to capture/track screenshots referenced across the web app docs.
:::

::::

```{toctree}
:maxdepth: 2
:hidden:

a_orientation/index
b_data_loading/index
c_core_interactions/index
d_fields_coloring_legends/index
e_filtering/index
f_highlighting_selection/index
g_cross_highlighting/index
h_analysis/index
i_vector_field_velocity/index
j_community_annotation/index
k_figure_export/index
l_sessions_sharing/index
n_benchmarking_performance/index
o_accessibility_privacy_security/index
p_developer_docs/index
q_troubleshooting_index/index
r_screenshot_checklist/index
```
