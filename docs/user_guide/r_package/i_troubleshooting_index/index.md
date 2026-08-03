# Troubleshooting Index

Pick the row that matches the last thing that happened.

| Where you are stuck | Page |
|---|---|
| `install.packages()` / `library(cellucid)` fails | {doc}`01_installation_and_dependency_issues` |
| `cellucid_prepare()` stops, or the export looks wrong | {doc}`02_data_preparation_issues` |
| The export folder ran but looks incomplete or unreadable | {doc}`03_export_format_and_validation_issues` |
| The folder exists but the web app will not load it | {doc}`04_web_app_loading_issues` |
| The export is too slow, too large, or too big to share | {doc}`05_performance_and_disk_usage_issues` |

```{note}
This index is the map. The exhaustive export-time guide, symptom by symptom,
is {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`, and the
web app has its own loading guide at
{doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`.
```

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`download;1.5em;sd-mr-1` Installation and Dependencies
:link: 01_installation_and_dependency_issues
:link-type: doc

Missing package, wrong library path, blocked GitHub, absent `jsonlite` or
`Matrix`, stale loaded version.
:::

:::{grid-item-card} {octicon}`package;1.5em;sd-mr-1` Data Preparation
:link: 02_data_preparation_issues
:link-type: doc

The four root causes behind most `cellucid_prepare()` failures, in priority
order.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Export Format and Validation
:link: 03_export_format_and_validation_issues
:link-type: doc

Missing files, `.gz` confusion, empty `var/`, gene search finding nothing,
`datasets.json` not found.
:::

:::{grid-item-card} {octicon}`browser;1.5em;sd-mr-1` Web App Loading
:link: 04_web_app_loading_issues
:link-type: doc

Folder picker refusing the selection, endless spinner, blank canvas, CORS,
fields that look wrong.
:::

:::{grid-item-card} {octicon}`meter;1.5em;sd-mr-1` Performance and Disk Usage
:link: 05_performance_and_disk_usage_issues
:link-type: doc

Exports that take forever, repositories too big for GitHub, and filesystem
file-count limits.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
