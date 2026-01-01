# Testing and CI (If Present)

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
- quantization uses reserved missing markers
- connectivity export writes edge pairs and chooses the right dtype
- vector fields are exported and scaled correctly

## Optional dependencies in tests

Some features/tests require `Matrix`. Tests typically skip gracefully if it’s missing.

## CI (if enabled)

If the repo has GitHub Actions workflows, typical checks include:
- `R CMD check` across OS/R versions
- `BiocCheck` for Bioconductor readiness

See the publishing guide for workflow expectations:
- `cellucid-r/PUBLISHING.md`
