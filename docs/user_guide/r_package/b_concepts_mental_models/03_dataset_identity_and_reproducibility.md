# Dataset Identity and Reproducibility

**Audience:** computational users, maintainers, anyone sharing datasets  
**Time:** 10 minutes  
**Goal:** understand what a store says about itself, and how to make stores reproducible.

## What the store says about itself: `uns['cellucid']`

Every store carries a mapping in `uns['cellucid']`, plus `uns['title']`:

| Key | From | Notes |
|---|---|---|
| `schema` | always `1` | the layout version of this mapping |
| `id` | `dataset_id` | the dataset's stable identifier |
| `name` | `dataset_name` | also written as `uns['title']` |
| `description` | `dataset_description` | only when supplied |
| `citation` | `citation_text`, `citation_url` | `{text, url}`; only when supplied, `url` optional |
| `source` | `source_name`, `source_url` | `{name, url}`; only when supplied, `url` optional |
| `defaults` | the launch-default arguments | only the ones you set |
| `outlierQuantiles` | `latent_key`, `centroid_min_points` | only with `latent_key`: `latentKey`, `minPoints`, and `fields`, one float32 quantile per cell for each categorical `obs` column |

`defaults` may hold `embedding1d` / `embedding2d` / `embedding3d`
(`embedding_defaults`; an obsm key, or `{key, columns}` with 0-based columns),
`geneIdColumn` (`var_gene_id_column`), `latentKey` (`latent_key`),
`vectorFieldDefault` (`vector_field_default`, the field name, which is also
the id the viewer lists, such as `velocity`), `obsKeys` (`obs_keys`), `normalizeEmbeddings`
(`normalize_embeddings`), and `centroidMinPoints` (`centroid_min_points`).
Every key you leave out stays the viewer's own choice.

## What counts as “the same dataset”

1) **The dataset id** (`uns['cellucid']['id']`). Sessions and community
   annotations are keyed on it, whichever way the store is opened (a local
   folder, a `.zarr.zip`, a link).
2) **The cells, in order.** Row order is cell identity, and the cell names in
   `obs` record it.
3) **The fields.** A session refers to obs columns, genes, and embeddings by name.

## The rules for `dataset_id` and `dataset_name`

Both are required; `cellucid_prepare()` never derives or rewrites them.

- `dataset_id`: 1–180 ASCII letters, digits, `.`, `_`, or `-`, starting with a
  letter or digit, not ending with `.`, and not a Windows device name — safe in
  every address and file name built from it.
- `dataset_name`, `dataset_description`, `citation_text`, `source_name`: shown
  verbatim, so no control or zero-width characters and no leading or trailing
  whitespace (an empty `dataset_description` is allowed).
- `citation_url`, `source_url`: absolute `http://` or `https://` addresses, each
  only together with its text (`citation_text`, `source_name`).

Use a stable, versioned id such as `pbmc3k_v1` or `mouse_brain_atlas_2024-08`.

## Reproducibility checklist

1) **Keep the script.** Store the exact R code that wrote the store.
2) **The bytes are reproducible.** The same inputs give a byte-identical store,
   folder or `.zarr.zip`: nothing time- or machine-dependent is written (ZIP
   entries carry a fixed 1980-01-01 timestamp).
3) **Record provenance in the call.** `citation_text`/`citation_url` and
   `source_name`/`source_url` travel with the store and are shown to whoever
   opens it.
4) **Freeze the inputs.** Note which assay/layer, which reductions, and which
   metadata columns you used.
5) **Make row order explicit.** Give `obs`, `gene_expression`, and the
   embeddings their cell names: `cellucid_prepare()` checks that they agree and
   stops at the first row where they do not.
6) **Version the dataset.** After a change that alters cells, fields, or
   embeddings, write a new store under a new `dataset_id` rather than
   overwriting the old one.

## Troubleshooting

### Symptom: “already exists; pass force = TRUE to replace it”

**Cause:** a store is already at `out_path`.

**Fix:** pass `force = TRUE` to replace it (the old store stays until the new
one is complete), or choose a new path. A folder that is not a Zarr store is
never replaced.

### Symptom: “My collaborator’s session doesn’t apply to my store”

**Likely causes:** a different `dataset_id`, a different cell order or number
of cells, or renamed columns.

**Fix:** agree on one `dataset_id` per dataset version, and write a new id when
the data changes.
