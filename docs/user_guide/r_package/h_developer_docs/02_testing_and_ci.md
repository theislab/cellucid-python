# Testing and CI

This page is for contributors who want to validate changes to `cellucid-r`.

## Running tests locally

From within the `cellucid-r/` directory, you can run tests with `testthat`.

### Option A: `devtools` (recommended)

```r
install.packages("devtools")
devtools::test()
```

### Option B: `R CMD check` (closer to CRAN/Bioc behavior)

From a shell:

```sh
R CMD check cellucid-r
```

## What the tests cover (high level)

The existing tests validate:
- core files are written (`dataset_identity.json`, manifests, points)
- embedding normalization matches the Python implementation
- categorical codes and generated nullable outlier quantiles use their exact
  reserved markers, while gene and continuous-observation inputs are finite-only
- connectivity export writes edge pairs and chooses the right dtype
- vector fields are exported and scaled correctly
- `DESCRIPTION` agrees with every other in-repo version site (`NEWS.md`,
  `CITATION.cff`, `README.md`, `vignettes/installation.Rmd`)
- the hand-written help page under `man/` still matches the formals of
  `cellucid_prepare()`: the `\usage` block is read back out of the `.Rd` file
  and compared to the real signature, the same comparison `R CMD check` makes

## Optional dependencies in tests

`Matrix` is a `Suggests` dependency, and `tests/testthat/test-connectivity-contract.R`
calls `Matrix::` directly without guarding for it. Install `Matrix` before
running the suite; without it those tests error rather than skip.

## CI

The checked-in `R-CMD-check.yaml` workflow runs `R CMD check` on current
macOS, Windows, Ubuntu R-devel, Ubuntu R-release, and Ubuntu oldrel-1. Its
dependency setup installs the package requirements, including `Matrix`, before
the check.

See the publishing guide for workflow expectations:
- `cellucid-r/publishing.md`
