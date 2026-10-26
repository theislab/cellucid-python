# Release Process

**Audience:** maintainers  
**Time:** 10 minutes  
**What you'll get:** the version sites that move together, and how a release is built and published

This page summarizes the intended release/publishing flow for `cellucid-r`.

The full step-by-step checklist lives in:
- `cellucid-r/publishing.md`

## Recommended release flow (high level)

1) Update `DESCRIPTION` version. It is the one source the rest are checked
   against.
2) Update `CITATION.cff`, `NEWS.md`, `README.md`, and
   `vignettes/installation.Rmd` to match, plus the pinned release literal in
   `tests/testthat/test-current-contract.R`.
3) Record user-visible changes in `NEWS.md`.
4) Ensure tests pass (`testthat::test_local()`, including the anndata round
   trip with `CELLUCID_ANNDATA_PYTHON` set, and `R CMD check --as-cran`).
5) Merge to main.
6) Create a GitHub Release tag. It must be exactly `v` followed by the
   `DESCRIPTION` version — `release.yaml` reads `DESCRIPTION`, computes
   `EXPECTED_TAG="v${VERSION}"`, and fails the build when the tag or the
   branch provenance does not match.
7) `release.yaml` then builds the source tarball, runs
   `R CMD check --as-cran` on it (a manual run can choose the default check
   instead), and uploads it with its SHA-256 as the `r-source-tarball`
   artifact.
8) Publish targets, in the order `cellucid-r/publishing.md` documents:
   - GitHub release (source of truth)
   - CRAN (manual upload of the tarball to the submission form)
   - r-universe (auto-builds from GitHub once the GitHub App is installed)
   - conda-forge (`r-cellucid`, built from the CRAN tarball)

`pkgdown.yaml` separately builds and deploys the package website at
<https://theislab.github.io/cellucid-r/>. That site is the R help pages; this
guide is the narrative documentation, and the two are maintained separately.

## Keeping this guide's version in step with the package

The page you are reading is part of the `cellucid-python` repository, and so is
{doc}`../installation`, which tells R users which version to install.
`cellucid-python/scripts/validate_release.py` checks those pages against the
**Python** package version in `cellucid-python/pyproject.toml`.

`cellucid-r` and `cellucid-python` are separate repositories with separate CI,
and neither checkout can read the other. Nothing therefore compares
`cellucid-r/DESCRIPTION` to this guide: bumping the Python version alone leaves
every gate green while the R installation page names a version the R package
never shipped, and bumping `cellucid-r/DESCRIPTION` alone leaves this guide
behind. The two versions are meant to move together; keeping them together is a
manual step in both release checklists.

- `cellucid-r`'s own test suite pins `DESCRIPTION` to every version site inside
  that repository, and `DESCRIPTION` carries a `Config/cellucid/version-marker`
  field so a `grep -rn CELLUCID_VERSION` sweep lists it.
- Releasing either package means bumping both to the same version and running
  `python scripts/validate_release.py` in `cellucid-python` plus
  `testthat::test_local()` in `cellucid-r`.

## What contributors should keep in mind

- The package writes a standard AnnData store plus `uns['cellucid']`; the web
  app and the Python package read and write the same contract.
- A change to what the store holds lands in the R writer, the Python writer,
  the web reader, their tests, and the documentation together; only the one
  current layout is kept.

## Where to update docs

When the store or the arguments change, update:
- the hand-written help pages under `cellucid-r/man/` (`cellucid_prepare.Rd`,
  `cellucid_view.Rd` and the other viewer pages), and the reference index in
  `cellucid-r/_pkgdown.yml` when a topic is added (the article index when a
  vignette is added; `.github/scripts/validate-pkgdown.R` refuses a site whose
  indexes miss one)
- `cellucid-r/README.md` and `cellucid-r/vignettes/cellucid.Rmd` (and
  `vignettes/viewer.Rmd` for the viewer)
- user guide API docs: {doc}`../c_data_preparation_api/index`
- the store specification: {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`
