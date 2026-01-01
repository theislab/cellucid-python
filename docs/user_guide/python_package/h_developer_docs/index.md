# Developer Docs (Python Package)

This section is the **developer-facing documentation** for the `cellucid` Python package (repo folder: `cellucid-python/`).

It explains how the Python side works (export format, servers, Jupyter embedding, hooks/events, sessions), how to contribute safely, and how to ship releases + docs.

:::{important}
If you are using Cellucid as an end user (wet lab / analyst / “I just want to view my data”), this is not the best starting point.
Start with the Python user guide landing page: {doc}`../index` and the Web App user guide: {doc}`../../web_app/index`.
:::

---

## Fast path (pick your goal)

| You want to… | Start here | Then read |
|---|---|---|
| Understand the big picture | {doc}`01_codebase_architecture` | {doc}`09_server_mode_architecture_endpoints_and_security` |
| Find “where in the repo is X?” | {doc}`02_repo_layout_and_entry_points` | {doc}`05_cli_architecture_and_commands` |
| Set up a local dev environment | {doc}`03_local_development_setup` | {doc}`13_testing_and_ci` |
| Change `prepare()` or the export format | {doc}`07_prepare_export_pipeline_architecture` | {doc}`08_export_format_spec_and_invariants` → {doc}`19_extension_point_add_new_export_feature` |
| Work on servers (`cellucid serve`, `serve_anndata`) | {doc}`09_server_mode_architecture_endpoints_and_security` | {doc}`16_performance_profiling_and_scaling` |
| Work on notebook embedding + hooks | {doc}`10_jupyter_embedding_architecture` | {doc}`11_hooks_events_protocol_and_schema` → {doc}`12_debugging_playbook` |
| Cut a release | {doc}`14_release_process` | {doc}`15_docs_development_and_style_guide` |

---

## Scope and mental model

Cellucid is an ecosystem:

- **Cellucid (web app)**: UI + rendering + interaction + sessions + exports. Code lives in the `cellucid/` repo.
- **cellucid-python (this repo)**: export format + servers + notebook embedding + hooks + Python/CLI APIs.
- **cellucid-annotation**: community annotation repo template + schema (web app feature, not primarily Python).
- **cellucid-datasets**: demo datasets (exports) used for examples/hosting.
- **cellucid-r (future)**: R exporter (not ready yet).

This developer section focuses on **the Python repo** but links to the **web app developer docs** when the behavior is frontend-driven:
{doc}`../../web_app/p_developer_docs/index`.

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
