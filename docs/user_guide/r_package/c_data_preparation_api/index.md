# Data Preparation API

This section is the **exhaustive documentation** for the R export API:

- `cellucid::cellucid_prepare(...)` (primary)
- `cellucid::prepare(...)` (alias)

It focuses on:
- exact input requirements (shapes, types, ordering),
- file outputs and on-disk formats,
- edge cases (missing values, huge datasets, weird column types),
- and troubleshooting (“why did this export fail?”).

```{note}
If you just want to export something quickly, start with {doc}`../a_landing_pages/04_quick_start_3_levels` and come back here when you need details.
```

**Recommended reading order**

1) {doc}`01_cellucid_prepare_overview`
2) {doc}`02_input_requirements_global`
3) {doc}`03_embeddings_and_coordinates`
4) {doc}`04_obs_cell_metadata`
5) {doc}`06_gene_expression_matrix` (+ {doc}`05_var_gene_metadata`)
6) {doc}`07_connectivities_knn_graph` (optional)
7) {doc}`08_vector_fields_velocity_displacement` (optional)
8) {doc}`09_output_format_specification_exports_directory`
9) {doc}`10_performance_tuning_prepare_export`
10) {doc}`11_troubleshooting_prepare_export`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
