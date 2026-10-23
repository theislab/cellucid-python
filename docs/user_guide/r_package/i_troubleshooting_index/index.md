# Troubleshooting Index

Pick the row that matches the last thing that happened.

| Where you are stuck | Page |
|---|---|
| `install.packages()` / `library(cellucid)` fails | {doc}`01_installation_and_dependency_issues` |
| `cellucid_prepare()` stops, or the store looks wrong | {doc}`02_data_preparation_issues` |
| The store was written but looks incomplete or unreadable | {doc}`03_store_format_and_validation_issues` |
| The store exists but the web app will not open it | {doc}`04_web_app_loading_issues` |
| Writing is slow, or the store is too big to share | {doc}`05_performance_and_disk_usage_issues` |

```{note}
This index is the map. The full guide to `cellucid_prepare()` errors, symptom by
symptom, is {doc}`../c_data_preparation_api/11_troubleshooting_prepare`, and
every error message is listed in
{doc}`../g_api_reference_coverage/02_error_messages_and_exceptions_document_patterns`.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`download;1.5em;sd-mr-1` Installation and Dependencies
:link: 01_installation_and_dependency_issues
:link-type: doc

Missing package, wrong library path, blocked GitHub, no C toolchain, absent
`jsonlite` or `Matrix`, stale loaded version.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Data Preparation
:link: 02_data_preparation_issues
:link-type: doc

The root causes behind most `cellucid_prepare()` failures, in priority order.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Store Format and Validation
:link: 03_store_format_and_validation_issues
:link-type: doc

Checking a store with anndata, hidden metadata files, a `.zarr.zip` that will
not open, genes or defaults missing.
:::

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web App Loading
:link: 04_web_app_loading_issues
:link-type: doc

Folder picker leaving out files, links refused, blank canvas, fields that look
wrong.
:::

:::{grid-item-card} {octicon}`meter;1.5em;sd-mr-1` Performance and Disk Usage
:link: 05_performance_and_disk_usage_issues
:link-type: doc

Slow writes, large stores, many files, and hosts that serve single files.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
