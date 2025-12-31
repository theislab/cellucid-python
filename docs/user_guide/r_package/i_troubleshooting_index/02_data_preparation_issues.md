# Data Preparation Issues

Use this page when `cellucid_prepare()` fails or produces exports that look wrong.

```{important}
The exhaustive prepare/export troubleshooting guide is:
{doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`.
```

## Most common root causes (in priority order)

### 1) Cell order mismatch

Symptoms:
- export “succeeds” but clusters/metadata do not match the embedding
- gene expression appears on the wrong cells

Fix:
- choose a canonical `cells` vector and reorder *every* input to match it
- see {doc}`../c_data_preparation_api/02_input_requirements_global`

### 2) Expression orientation mismatch

Symptoms:
- `var has X rows, but gene_expression has Y genes`
- or export is slow/huge unexpectedly because you exported the wrong dimension

Fix:
- ensure `gene_expression` is **cells × genes**
- Seurat/SCE usually need `Matrix::t(...)`
- see {doc}`../c_data_preparation_api/06_gene_expression_matrix`

### 3) Too many genes exported

Symptoms:
- exports are enormous
- export time is “forever”

Fix:
- export a gene panel (`gene_identifiers`)
- use quantization (`var_quantization=8`) and compression (`compression=6`)
- see {doc}`../c_data_preparation_api/10_performance_tuning_prepare_export`

### 4) Metadata column types are wrong

Symptoms:
- numeric-looking columns become categorical
- legends explode into thousands of categories

Fix:
- coerce columns intentionally (factor vs numeric)
- export fewer columns (`obs_keys`)
- see {doc}`../c_data_preparation_api/04_obs_cell_metadata`

## Quick “what do I do next?” links

- Seurat export recipe: {doc}`../e_integrations_recipes/01_seurat_recipe`
- SCE export recipe: {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`
- Format/validation issues: {doc}`03_export_format_and_validation_issues`
