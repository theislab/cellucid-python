# Concepts & Mental Models

These pages help you build a correct mental model for the R workflow:

- What is an “export folder”?
- How does data move from R → disk → browser?
- What makes two exports “the same dataset” for sharing and reproducibility?

```{note}
If you just want to export something quickly, you can skip to {doc}`../a_landing_pages/04_quick_start_3_levels`.
If you are debugging confusing behavior, these concept pages are worth reading.
```

**Recommended reading order**

1) {doc}`01_what_is_an_export_folder`
2) {doc}`02_data_flow_r_to_browser`
3) {doc}`03_dataset_identity_and_reproducibility`

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Export Folder
:link: 01_what_is_an_export_folder
:link-type: doc

What gets written, why there are many small files, and how the viewer uses manifests.
:::

:::{grid-item-card} {octicon}`workflow;1.5em;sd-mr-1` Data Flow
:link: 02_data_flow_r_to_browser
:link-type: doc

The three viewing modes: local file picker, local static server, and hosted exports.
:::

:::{grid-item-card} {octicon}`shield-check;1.5em;sd-mr-1` Dataset Identity
:link: 03_dataset_identity_and_reproducibility
:link-type: doc

How `dataset_identity.json` affects sharing, sessions, and reproducibility.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
