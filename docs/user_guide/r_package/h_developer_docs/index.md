# Developer Docs

If you are only *using* the package, you can skip this section — start at
{doc}`../a_landing_pages/index` instead.

**Recommended reading order**

1) {doc}`01_codebase_architecture`
2) {doc}`02_testing_and_ci`
3) {doc}`03_release_process`

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`file-directory;1.5em;sd-mr-1` Codebase Architecture
:link: 01_codebase_architecture
:link-type: doc

What each file under `R/` is responsible for, and the order
`cellucid_prepare()` runs the export in.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Testing and CI
:link: 02_testing_and_ci
:link-type: doc

Running `testthat` and `R CMD check` locally, what the suite covers, and which
OS/R combinations CI runs.
:::

:::{grid-item-card} {octicon}`tag;1.5em;sd-mr-1` Release Process
:link: 03_release_process
:link-type: doc

The version sites that must move together, the tag the release workflow
requires, and the publishing targets.
:::

::::

```{note}
The export format itself is a contract with four other repositories. Before
changing anything a file layout depends on, read
{doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`.
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
