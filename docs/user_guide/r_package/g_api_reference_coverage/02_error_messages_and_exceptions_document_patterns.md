# Error Messages and Exceptions (Documentation Patterns)

This page is for users and maintainers who want to interpret `cellucid_prepare()` failures quickly.

## How `cellucid_prepare()` reports problems

`cellucid_prepare()` reports every problem the same way: it calls `stop(...)` and
the export stops. There is no partial success and no severity ladder.

- The R sources contain no `warning()` call, so the package raises none of its
  own. Nothing is skipped, downgraded, or repaired for you.
- Rejection happens before publication. A rejected call publishes nothing: it
  creates no output directory, and with `force = TRUE` it leaves any existing
  generation exactly as it was.

When reporting a bug, always include:
- the full error message,
- the call you used (arguments),
- and the shapes of your inputs (see the checklist in {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`).

## Common error families (and where to debug)

### “At least one embedding must be provided: X_umap_1d, X_umap_2d, or X_umap_3d.”

Meaning:
- you didn’t pass any of `X_umap_1d/2d/3d`

Fix:
- provide at least one embedding

Docs:
- {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`

### “X_umap_2d must have exactly 2 columns, got 3.”

Meaning:
- you passed an embedding with the wrong number of columns

Vector fields report the related `Vector field '<name>' must have exactly 1, 2,
or 3 components.` when the array has no usable width at all, and
`Vector field '<name>' declares 2D but contains 1 components.` when its width
disagrees with the dimension its key declares.

Docs:
- embeddings: {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`
- vector fields: {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`

### “Vector field key 'velocity_umap' must exactly match …”

Meaning:
- a `vector_fields` list name is not a complete declaration. The full message
  quotes the required pattern, `'<field>_umap_<1|2|3>d'`. The dimension is
  declared by the key and never inferred, so an unsuffixed `velocity_umap` is
  rejected rather than sized from the array.

Fix:
- rename the list entry to `velocity_umap_2d` (or `_1d` / `_3d`) to match the
  array you are passing. The Python package enforces the identical rule.

Docs:
- {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`

### “gene_expression has 500 cells, but embeddings have 2000 cells.”

Meaning:
- row-order / subsetting mismatch across inputs

The same wording covers `Latent space has ... cells, but embeddings have ...
cells.` and `var has ... rows, but gene_expression has ... genes.`

Docs:
- {doc}`../c_data_preparation_api/02_input_requirements_global`

### “var must be a data.frame when gene_expression is supplied.”

Meaning:
- you provided expression without gene metadata

Fix:
- pass a `var` `data.frame` with one row per gene, in the column order of
  `gene_expression`

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`
- {doc}`../c_data_preparation_api/06_gene_expression_matrix`

### “gene_identifiers contains identifiers not found in var: MS4A1, CD3D”

Meaning:
- you requested genes via `gene_identifiers` that are absent from your gene ID
  list

This is a hard error, not a warning. Missing genes are **not** skipped and the
export does **not** proceed with the intersection.

Fix:
- intersect your panel with your gene IDs before the call:
  ```r
  markers <- intersect(markers, rownames(var))
  ```
- check `var_gene_id_column` first: the identifiers are read from
  `rownames(var)` when it is `NULL`, and from that exact column otherwise.

The message lists the first five missing identifiers and appends `...` when
there are more. Duplicated entries fail separately with
`gene_identifiers must not contain duplicate identifiers.`

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`

### “var_gene_id_column 'gene_symbol' not found in var. Available columns: …”

Meaning:
- you named a column that `var` does not have

Fix:
- pass `NULL` to use `rownames(var)`, or name an existing column exactly

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`

### “obs_categorical_dtype must be exactly one of "uint8" or "uint16".”

Meaning:
- `obs_categorical_dtype` is required and has no default. There is no automatic
  width selection, so `"auto"` and a missing argument both fail.

Fix:
- pass `"uint8"` for up to 255 categories per field, `"uint16"` for up to 65,535

Docs:
- {doc}`../c_data_preparation_api/04_obs_cell_metadata`

### “Gene identifiers 'HLA-DRB1/2' is not a portable identifier. …”

Meaning:
- an exported identifier cannot become a filename.

Every message of this family opens with the axis that raised it, so the first
words tell you which input to go and fix:

| Opening words | What was rejected |
|---|---|
| `Gene identifiers` | a gene ID selected for export |
| `Observation field keys` | a key read out of a `compact_v1` obs manifest |
| `obs_keys` | a name you passed in the `obs_keys` argument |
| `Vector field ids` | a `vector_fields` field id |
| `dataset_id` | the `dataset_id` argument |

`Observation field keys` and `obs_keys` are two different inputs, not two names
for one: the first is parsed out of a manifest you supplied, the second is the
argument you passed.

The full message states the rule: `Use 1-180 ASCII letters, numbers, '.', '_',
or '-', beginning with a letter or number and not ending with '.'.` Windows
device names get their own message under the same axis name, for example
`Vector field ids 'CON.velocity_umap' is reserved on Windows.`
Case-insensitive collisions are reported the same way and then list each
colliding group:

```text
Vector field ids contains names that collide on case-insensitive filesystems.
Rename the identifiers. Collisions:
  - 'velocity_umap' <- 'Velocity_UMAP', 'velocity_umap'
```

A vector-field ID is checked after its key is parsed, so the ID reported is the
part before the dimensional suffix.

Docs:
- {doc}`../c_data_preparation_api/05_var_gene_metadata`
- {doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`
- {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`

### “Categorical field 'organ' has 2 labels the viewer cannot show as written”

Meaning:
- a category label carries a character with no glyph, so the stored value is
  not the value a reader sees.

A label is drawn verbatim in the legend, the field selector, and every exported
figure. `cellucid_prepare()` rejects an empty label, a control character
(`U+0001`-`U+001F`, `U+007F`-`U+009F`), a zero-width character (`U+200B`,
`U+2060`, `U+FEFF`), and leading or trailing whitespace of any kind including
`U+00A0` NO-BREAK SPACE:

```text
Categorical field 'organ' has 2 labels the viewer cannot show as written:
  - 'Gut ' ends with U+0020 SPACE
  - 'Liver ' ends with U+0020 SPACE
Labels are drawn verbatim in the legend, the field selector, and exported
figures, so these characters are invisible on screen while still making the
label a different value from the one it looks like. Cellucid does not clean
them for you: trimming a label rewrites your annotation and can merge two
categories into one, moving cells between them. Clean the column and export
again, for example:
    obs[['organ']] <- factor(trimws(as.character(obs[['organ']]), whitespace = "[\\h\\v]"))
```

Every offending label in the field is listed in one message (up to ten, then a
count of the rest), so one edit fixes the column. Two labels that a
whitespace-collapsing renderer draws identically get their own message naming
both, for example `'T  cell'` and `'T cell'`.

`dataset_name`, `dataset_description`, `source_name`, `source_url`, and
`source_citation` use the same vocabulary and open with the argument name:
`source_name is displayed verbatim, so it must not carry characters that have
no glyph: 'Suo\u200b' contains U+200B ZERO WIDTH SPACE. ...`

The Python package raises the identical rule from `prepare()`.

Docs:
- {doc}`../c_data_preparation_api/04_obs_cell_metadata`
- {doc}`../b_concepts_mental_models/03_dataset_identity_and_reproducibility`

### “out_dir already exists; set force = TRUE to replace the complete generation.”

Meaning:
- the target directory is already there and `force` is `FALSE`

Fix:
- pass `force = TRUE` to replace the whole generation atomically, or choose a
  new `out_dir`

Docs:
- {doc}`../b_concepts_mental_models/01_what_is_an_export_folder`

### “An export generation is already active for /path/to/exports.”

Meaning:
- another live R or Python exporter owns that exact output directory

Fix:
- wait for the other export to finish; do not delete the hidden lock file

Docs:
- {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`

### “Matrix package is required …”

Meaning:
- you’re exporting sparse matrices or connectivities without Matrix installed

Fix:
```r
install.packages("Matrix")
```

### “Connectivity matrix must be exactly symmetric.”

Meaning:
- the graph is not symmetric in topology and weight

Neighbouring messages cover the rest of the graph contract:
`Connectivity diagonal values must all be exactly zero.`,
`Connectivity weights must all be non-negative.`,
`Connectivity values must all be finite.`, and
`Connectivity contains duplicate sparse coordinates at (…)`.

Fix:
- symmetrize the graph and drop stored zeros before the call

Docs:
- {doc}`../c_data_preparation_api/07_connectivities_knn_graph`

### “Every nonzero finite float32 value must be representable between 2^-149 and (2 - 2^-23) * 2^127 in absolute magnitude.”

Meaning:
- a value would become zero or infinity once written as float32

Related quantization messages read `Cannot quantize field '<name>': its
native-double variation collapses to one value in the viewer's float32 domain.`
and `Cannot quantize field '<name>': all values are constant, …`.

Docs:
- {doc}`../c_data_preparation_api/06_gene_expression_matrix`
- {doc}`../c_data_preparation_api/10_performance_tuning_prepare_export`
