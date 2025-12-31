# Public Functions and Objects

`cellucid-r` exports two user-facing functions (they are the same underlying implementation):

## `cellucid_prepare()`

Primary export function.

Core docs:
- overview: {doc}`../c_data_preparation_api/01_cellucid_prepare_overview`
- full API behavior: {doc}`../c_data_preparation_api/index`

Signature (high level):
- required:
  - `latent_space`
  - `obs`
  - at least one of `X_umap_1d`, `X_umap_2d`, `X_umap_3d`
- optional:
  - `var`, `gene_expression`, `gene_identifiers`, `connectivities`, `vector_fields`
- tuning:
  - `compression`, `var_quantization`, `obs_continuous_quantization`, `obs_keys`, `force`

## `prepare()`

Alias of `cellucid_prepare()`.

Why it exists:
- mirrors the Python package (`cellucid.prepare`) naming for cross-language consistency.

Docs:
- treat it exactly the same as `cellucid_prepare()`.
