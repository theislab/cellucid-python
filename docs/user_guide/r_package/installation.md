# Installation

**Audience:** everyone  
**Time:** 5–15 minutes (depending on your setup)  
**Goal:** install Cellucid for R, verify version 0.9.5, and choose the correct <!-- CELLUCID_VERSION -->
source for current registry availability

:::{note}
**Active package version: 0.9.5.** Version 0.9.5 is the CRAN submission
release. The
[CRAN package index](https://cran.r-project.org/web/packages/available_packages_by_name.html)
is authoritative for registry availability. Use the CRAN command below when
the index lists `cellucid` 0.9.5; otherwise install the active source from the
official GitHub repository.
:::

## Requirements

- **R:** `cellucid` targets **R ≥ 4.1.0** (checked on R 4.1, 4.2, 4.3, oldrel, release and devel).
- **Hard dependencies:** `jsonlite`, and for the viewer `httpuv`, `later` and `curl`
- **Optional:** `Matrix` (sparse expression matrices and neighbour graphs; recommended),
  `SeuratObject` and `SingleCellExperiment` / `SummarizedExperiment` (each only
  to read its own objects), `DelayedArray` (assays kept on disk)
- **Source builds:** Rtools on Windows, the Xcode Command Line Tools on macOS,
  or a C toolchain plus the R and zlib development headers on Linux (for
  example `r-base-dev`, which brings `zlib1g-dev`, on Debian and Ubuntu)

```{note}
`cellucid-r` does not require Seurat or SingleCellExperiment: it reads their
objects when they are installed, and works from matrices and data frames when
they are not.
```

## Install Cellucid

### From CRAN when version 0.9.5 is listed
<!-- CELLUCID_VERSION -->

When the CRAN package index lists `cellucid` 0.9.5, install it with: <!-- CELLUCID_VERSION -->

```r
install.packages("cellucid")
```

If the index does not list version 0.9.5 yet, use the official source <!-- CELLUCID_VERSION -->
installation below.

### From the official source repository

These routes compile Cellucid's small C routines (the CRC-32 of the
`.zarr.zip` writer, the SHA-256 that checks the downloaded web app, the
latent-space outlier quantiles, one-pass matrix checks, and the decoders of a
saved session, linked against zlib) and therefore require the source-build
toolchain listed above. A platform CRAN
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
- Version `0.9.5` from `packageVersion` <!-- CELLUCID_VERSION -->

## Optional dependencies (recommended)

### `Matrix`

If you plan to write:
- sparse gene expression matrices (common), or
- a neighbour graph (`connectivities=`),

install Matrix:

```r
install.packages("Matrix")
```

### Seurat / SingleCellExperiment (only for their objects)

`cellucid_prepare(data = )` and `cellucid_view()` read a Seurat object with
`SeuratObject` and a SingleCellExperiment or SummarizedExperiment with
`SingleCellExperiment` / `SummarizedExperiment`; each is loaded only when an
object of its class arrives, and a missing one is named with the line that
installs it.

On R older than 4.4, `install.packages("SeuratObject")` fails because
SeuratObject 5 needs Matrix 1.6-4 or later and CRAN's current Matrix needs
R 4.4. Install Matrix 1.6-5, the last release for those R versions, first:

```r
remotes::install_version("Matrix", version = "1.6-5")
install.packages("SeuratObject")
```

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

- Ready to write your first store? Go to {doc}`a_landing_pages/04_quick_start_3_levels`.
