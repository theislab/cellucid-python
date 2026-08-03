# Codebase Architecture

**Audience:** contributors and maintainers  
**Time:** 15 minutes  
**What you'll get:** what each file under `R/` owns, and the order the export runs in

`cellucid-r` exposes one public function, `cellucid_prepare()`, and splits the
writer behind it across `R/` by responsibility: one module per axis, per encoding
concern, and per publication step.

## Repository layout (high level)

- `cellucid-r/R/cellucid_prepare.R`
  - exports `cellucid_prepare()`, the only public function
  - validates the arguments, then drives the export in order and delegates each
    step to the modules below
- `cellucid-r/R/` — argument and identity rules
  - `validate-arguments.R` — one function per argument value rule, plus the
    reporter for arguments that were never supplied at all
  - `validate-display-text.R` — the rule for text the viewer draws verbatim
    (names, descriptions, source lines, category labels, field identifiers);
    a cross-writer contract with `_display_text_defect()` in the Python package
  - `identifiers.R` — what a name may be: the portable-filename rule for
    `dataset_id`, and the display-text plus uniqueness rules for field
    identifiers
  - `dataset-identity.R` — the dataset's own name, id, description, creation
    time, and source
  - `output-path.R` — where the export may be written
  - `error-messages.R` — the one list syntax every message uses to show a set of
    values back to the caller
- `cellucid-r/R/` — input axes
  - `matrix-input.R` — base matrix, numeric data frame, or `Matrix::Matrix` into
    one internal storage
  - `embeddings.R` — the drawn coordinates, centered and scaled; keeps the scale
    for the vector fields drawn on top of them
  - `obs.R` — field-kind detection, categories and codes, and every obs payload
  - `var.R` — naming the exported genes and reading one gene's column
  - `connectivity.R` — the KNN graph, published as edge pairs from the strict
    upper triangle
  - `vector-fields.R` — per-cell vectors, scaled by their embedding's own scale
    and ordered by code point so both writers assign the same indices
  - `centroids.R` — category label positions and per-cell outlier distances
- `cellucid-r/R/` — encoding and manifests
  - `quantization.R` — continuous values narrowed to an integer code range,
    including the named constant-range case
  - `binary-payloads.R` — every payload byte, in the little-endian widths the
    browser contract names
  - `gzip-header.R` — the fixed ten-byte gzip header that makes a compressed
    payload byte-identical on every machine
  - `dtype-names.R` — the one place the bit width, the manifest dtype, and the
    payload file extension are derived from each other
  - `manifest.R` — the manifests, re-expanded and compared against the payloads
    actually written before the generation is published
- `cellucid-r/R/` — publication safety
  - `export-lock.R` — one export generation at a time per output directory:
    the native advisory lock plus the in-process registry
  - `export-transaction.R` — publishing a generation is one step or none, with a
    journal that resolves an interrupted export from what is on disk
  - `export-filesystem.R` — guarded filesystem operations shared by the lock and
    the transaction; a symlink or unexpected file type stops the export
  - `native.R` — the complete `.Call()` surface of `src/export_lock.c`
  - `zzz.R` — package unload hook that drains the native lock
- `cellucid-r/src/export_lock.c`
  - the native export lock reached through `.Call()`
- `cellucid-r/man/cellucid_prepare.Rd`
  - the R help page, hand-written and authoritative. This package runs no
    documentation generator and `R/` holds no roxygen comments, so edit the
    `.Rd` file directly.
- `cellucid-r/NAMESPACE`
  - hand-written too, and the only place the `useDynLib()` registration for
    `src/export_lock.c` exists
- `cellucid-r/tests/testthat/`
  - unit tests validating core files, normalization, quantization, connectivity, vector fields
- `cellucid-r/vignettes/cellucid.Rmd`
  - minimal vignette showing a small export workflow
- `cellucid-r/publishing.md`
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
   - validate exact finite binary values, symmetry, and zero diagonal
   - write edge pairs under `connectivity/`
   - write `connectivity_manifest.json`
8) Write `dataset_identity.json` (summary + pointers to files).

The user-guide docs mirror this structure:
- {doc}`../c_data_preparation_api/index`

## Design principles

- Minimal dependencies (only `jsonlite` required).
- Deterministic, file-based exports.
- Shared format with the Python exporter.
- One persistent, package-owned native lock per target prevents concurrent R
  and Python exporters from splitting a published generation.

## Adding a new exported artifact (maintainer notes)

If you add a new feature that writes files:

1) Decide where it belongs:
   - `dataset_identity.json` (top-level discovery)
   - a dedicated manifest JSON
   - a new subdirectory of binaries
2) Add tests under `cellucid-r/tests/testthat/`.
3) Update the user guide format spec:
   - {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
