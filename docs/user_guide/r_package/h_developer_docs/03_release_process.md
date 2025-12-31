# Release Process

This page summarizes the intended release/publishing flow for `cellucid-r`.

The full step-by-step checklist lives in:
- `cellucid-r/docs/publishing.md`

## Recommended release flow (high level)

1) Update `DESCRIPTION` version.
2) Update `NEWS.md` with user-visible changes.
3) Ensure tests pass (`devtools::test()` / `R CMD check`).
4) Merge to main.
5) Create a GitHub Release tag (e.g. `v0.99.1`).
6) Use CI to build a source tarball (`cellucid_<version>.tar.gz`).
7) Publish targets (depending on maturity):
   - GitHub (source of truth)
   - r-universe (easy binaries)
   - CRAN (manual submission)
   - Bioconductor (manual submission via contributions tracker)

## What contributors should keep in mind

- The R package is an exporter that must remain compatible with the viewer format.
- Changes to file formats should be documented and tested.
- Prefer backwards-compatible changes unless a format bump is intentional.

## Where to update docs

When export behavior changes, update:
- user guide API docs: {doc}`../c_data_preparation_api/index`
- format spec: {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
