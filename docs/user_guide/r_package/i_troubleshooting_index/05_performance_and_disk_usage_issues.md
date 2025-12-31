# Performance and Disk Usage Issues

Use this page when exports are too slow, too large, or too hard to share.

## Symptom: export takes “forever”

**Likely causes**
- exporting too many genes
- writing to a slow filesystem (network drive, cloud-synced folder)

**How to confirm**
- count exported gene files:
  ```r
  length(list.files(file.path(out_dir, "var")))
  ```

**Fix**
- export a gene panel (`gene_identifiers`)
- set `var_quantization=8`
- set `compression=6`
- export to a local SSD

## Symptom: export folder is too big for GitHub

**Likely causes**
- full gene export (20k+ files) at large cell counts
- float32 export instead of quantized

**Fix**
- gene panel + quantization + compression (same as above)
- or use server mode for private sharing

## Symptom: the web app is slow after loading

**Likely causes**
- dataset too large for browser/GPU
- too many fields/categories

**Fix**
- reduce cells (subset)
- export fewer metadata columns
- prefer 2D embedding for weaker GPUs
- see: {doc}`../../web_app/n_benchmarking_performance/index`

## Symptom: “too many files” / filesystem errors

**Likely cause**
- exporting too many genes creates tens of thousands of files.

**Fix**
- gene panel
- consider splitting datasets by modality/cell type

## Deep guidance

- Performance tuning guide: {doc}`../c_data_preparation_api/10_performance_tuning_prepare_export`
- Gene expression pitfalls: {doc}`../c_data_preparation_api/06_gene_expression_matrix`
