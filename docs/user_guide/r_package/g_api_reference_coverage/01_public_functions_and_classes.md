# Public Functions and Objects

`cellucid-r` exports one user-facing function:

## `cellucid_prepare()`

Primary export function.

Core docs:
- overview: {doc}`../c_data_preparation_api/01_cellucid_prepare_overview`
- full API behavior: {doc}`../c_data_preparation_api/index`

Signature (high level):
- required, and carrying no default at all, so leaving one out is reported as a
  missing argument:
  - `obs_categorical_dtype`
  - `dataset_name`
  - `dataset_id`
- required data:
  - `latent_space`
  - `obs`
  - at least one of `X_umap_1d`, `X_umap_2d`, `X_umap_3d`
- optional:
  - `var`, `gene_expression`, `gene_identifiers`, `connectivities`, `vector_fields`
- tuning:
  - `compression`, `var_quantization`, `obs_continuous_quantization`, `obs_keys`, `force`
