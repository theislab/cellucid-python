# Data Preparation API (prepare/export) — The Big One

This section is the **exhaustive, no-surprises guide** to exporting data for Cellucid using:

- {func}`~cellucid.prepare` (Python API)

Cellucid is a **web app**. `cellucid-python` is a **helper package** you use from Python/CLI to:
- export datasets into a high-performance folder format (what this section covers),
- or serve/view data directly (see {doc}`../d_viewing_apis/index`).

This export format is designed to be:
- easy to host/share (static files + manifests),
- fast to load in the browser,
- reproducible (explicit manifests + stable dataset identity),
- and compatible with future helpers (e.g., a future `cellucid-R` exporter).

This section focuses on:
- exact input requirements (shapes, dtypes, ordering),
- on-disk outputs (what files get written and why),
- performance knobs (compression, quantization, subsetting),
- edge cases (NaN/Inf, huge categories, sparse vs dense),
- and troubleshooting (symptom → diagnosis → fix).

```{note}
If you just want to visualize something quickly (no export deep dive yet), start with:
- {doc}`../a_landing_pages/04_quick_start_3_levels` (choose your level), or
- {doc}`../d_viewing_apis/01_viewing_methods_overview` (all viewing modes).
```

## Recommended reading order

1) {doc}`01_prepare_export_overview`
2) {doc}`02_input_requirements_global`
3) {doc}`03_embeddings_and_coordinates`
4) {doc}`04_obs_cell_metadata`
5) {doc}`06_gene_expression_matrix` (+ {doc}`05_var_gene_metadata`)
6) {doc}`07_connectivities_knn_graph` (optional)
7) {doc}`08_vector_fields_velocity_displacement` (optional)
8) {doc}`09_output_format_specification_exports_directory`
9) {doc}`10_performance_tuning_guide_prepare_export`
10) {doc}`11_troubleshooting_prepare_export`

## “What do I need to export?” (at a glance)

Minimum viable export (fastest to produce):
- `latent_space` + `obs` + at least one embedding (`X_umap_2d` or `X_umap_3d`) → interactive viewer with metadata fields.

Common “most useful” export:
- add `gene_expression` + `var` → gene search + gene overlays.

Optional advanced features:
- add `connectivities` → graph-based features (KNN edges).
- add `vector_fields` → animated velocity/displacement overlays.

## API reference (when you need exact signatures)

- {doc}`../api/export` (includes {func}`~cellucid.prepare` + docstring/autodoc)

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
