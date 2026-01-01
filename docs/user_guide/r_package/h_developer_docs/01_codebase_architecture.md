# Codebase Architecture

`cellucid-r` is intentionally small. Most of the implementation lives in one file.

## Repository layout (high level)

- `cellucid-r/R/cellucid_prepare.R`
  - exports `cellucid_prepare()` and `prepare()`
  - contains the exporter implementation and helper functions
- `cellucid-r/man/cellucid_prepare.Rd`
  - R help page generated/maintained for Bioconductor-style docs
- `cellucid-r/tests/testthat/`
  - unit tests validating core files, normalization, quantization, connectivity, vector fields
- `cellucid-r/vignettes/cellucid.Rmd`
  - minimal vignette showing a small export workflow
- `cellucid-r/PUBLISHING.md`
  - release/publishing checklist

## Data export pipeline (what happens in `cellucid_prepare()`)

At a high level:

1) **Validate embeddings** and infer `n_cells`.
2) **Normalize embeddings** (center + scale) and write `points_*d.bin`.
3) Validate/convert `latent_space` and `obs`.
4) Export optional **vector fields** (scaled with embedding normalization).
5) Export **obs**:
   - continuous values (float32 or quantized)
   - categorical codes (uint8/uint16) + outlier quantiles (latent-space)
   - centroids (embedding-space)
   - write `obs_manifest.json`
6) Export **gene expression** (optional):
   - validate `var` and `gene_expression`
   - write one dense vector per gene under `var/`
   - write `var_manifest.json`
7) Export **connectivities** (optional):
   - symmetrize and binarize
   - write edge pairs under `connectivity/`
   - write `connectivity_manifest.json`
8) Write `dataset_identity.json` (summary + pointers to files).

The user-guide docs mirror this structure:
- {doc}`../c_data_preparation_api/index`

## Design principles

- Minimal dependencies (only `jsonlite` required).
- Deterministic, file-based exports.
- Shared format with the Python exporter.

## Adding a new exported artifact (maintainer notes)

If you add a new feature that writes files:

1) Decide where it belongs:
   - `dataset_identity.json` (top-level discovery)
   - a dedicated manifest JSON
   - a new subdirectory of binaries
2) Add tests under `cellucid-r/tests/testthat/`.
3) Update the user guide format spec:
   - {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
