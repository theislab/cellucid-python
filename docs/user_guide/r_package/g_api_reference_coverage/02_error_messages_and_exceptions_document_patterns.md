# Error Messages and Exceptions (Documentation Patterns)

This page is for users and maintainers who want to interpret `cellucid_prepare()` failures quickly.

## How `cellucid_prepare()` reports problems

- Fatal issues are raised with `stop(...)` (execution stops).
- Non-fatal issues (like missing requested genes) may be raised with `warning(...)`.

When reporting a bug, always include:
- the full error message,
- the call you used (arguments),
- and the shapes of your inputs (see the checklist in {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`).

## Common error families (and where to debug)

### “At least one embedding must be provided…”

Meaning:
- you didn’t pass any of `X_umap_1d/2d/3d`

Fix:
- provide at least one embedding

Docs:
- {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`

### “... must have exactly N columns”

Meaning:
- you passed an embedding/vector field with the wrong number of columns

Docs:
- embeddings: {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`
- vector fields: {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`

### “... has X cells, but embeddings have Y cells”

Meaning:
- row-order / subsetting mismatch across inputs

Docs:
- {doc}`../c_data_preparation_api/02_input_requirements_global`

### “var data.frame must be provided when gene_expression is given”

Meaning:
- you provided expression without gene metadata

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`
- {doc}`../c_data_preparation_api/06_gene_expression_matrix`

### “Matrix package is required …”

Meaning:
- you’re exporting sparse matrices or connectivities without Matrix installed

Fix:
```r
install.packages("Matrix")
```

## Common warnings

### “gene identifiers not found in var…”

Meaning:
- you requested genes via `gene_identifiers`, but some are missing from your gene ID list

Fix:
- intersect your panel with your gene IDs first, or accept that missing genes are skipped

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`
