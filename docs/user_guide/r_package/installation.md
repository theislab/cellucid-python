# Installation

**Audience:** everyone  
**Time:** 5–15 minutes (depending on your setup)
**Goal:** install Cellucid for R, verify version 0.9.1, and choose the correct <!-- CELLUCID_VERSION -->
source for current registry availability

:::{note}
**Active package version: 0.9.1.** Version 0.9.1 is the CRAN submission
release. The
[CRAN package index](https://cran.r-project.org/web/packages/available_packages_by_name.html)
is authoritative for registry availability. Use the CRAN command below when
the index lists `cellucid` 0.9.1; otherwise install the active source from the
official GitHub repository.
:::

## Requirements

- **R:** `cellucid` currently targets **R ≥ 4.3.0**.
- **Hard dependency:** `jsonlite`
- **Optional but recommended:** `Matrix` (sparse matrices + connectivities)
- **Source builds:** Rtools on Windows, the Xcode Command Line Tools on macOS,
  or a C toolchain plus R development headers on Linux

```{note}
`cellucid-r` is designed to be dependency-light and does not require Seurat or SingleCellExperiment.
Those packages are only needed if you want to follow the Seurat/SCE recipes in this guide.
```

## Install Cellucid

### From CRAN when version 0.9.1 is listed
<!-- CELLUCID_VERSION -->

When the CRAN package index lists `cellucid` 0.9.1, install it with: <!-- CELLUCID_VERSION -->

```r
install.packages("cellucid")
```

If the index does not list version 0.9.1 yet, use the official source <!-- CELLUCID_VERSION -->
installation below.

### From the official source repository

These routes compile Cellucid's small native export-lock primitive and
therefore require the source-build toolchain listed above. A platform CRAN
binary, when available, does not require a local compiler.

#### Option A: `remotes` (simple, common)

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

#### Option B: `pak` (often faster, better resolver)

```r
install.packages("pak")
pak::pak("theislab/cellucid-r")
```

## Verify the installation

```r
library(cellucid)
packageVersion("cellucid")
```

You should see:
- No error on `library(cellucid)`
- Version `0.9.1` from `packageVersion` <!-- CELLUCID_VERSION -->

## Optional dependencies (recommended)

### `Matrix`

If you plan to export:
- sparse gene expression matrices (common), or
- connectivity matrices (`connectivities=`),

install Matrix:

```r
install.packages("Matrix")
```

### Seurat / SingleCellExperiment (only for the recipes)

These are optional and only needed if you want to follow the extraction tutorials:
- {doc}`e_integrations_recipes/01_seurat_recipe`
- {doc}`e_integrations_recipes/02_singlecellexperiment_recipe`

## Troubleshooting installation

### Symptom: “package ‘jsonlite’ is not available”

**Likely causes**
- You are using an R environment without CRAN configured.
- You are offline or behind a restrictive proxy.

**How to confirm**
- Run `getOption("repos")` and verify you have a CRAN mirror set.

**Fix**
```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("jsonlite")
```

### Symptom: GitHub install fails (“HTTP error 401/403”, “rate limit”)

**Likely causes**
- Corporate proxy / blocked GitHub
- GitHub API rate limiting

**Fix options**
- Try `pak` (often more resilient).
- If needed, configure a GitHub PAT:
  - `Sys.setenv(GITHUB_PAT = "...")`

### Symptom: “there is no package called ‘Matrix’”

**Fix**
```r
install.packages("Matrix")
```

### Symptom: you installed but `library(cellucid)` loads the wrong version

**How to confirm**
```r
find.package("cellucid")
packageVersion("cellucid")
```

**Fix**
- Restart R (RStudio: Session → Restart R).
- Reinstall with:
  ```r
  remotes::install_github("theislab/cellucid-r", force = TRUE)
  ```

## Next steps

- Ready to export? Go to {doc}`a_landing_pages/04_quick_start_3_levels`.
