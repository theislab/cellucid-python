# Installation and Dependency Issues

**Audience:** everyone  
**Time:** 5–10 minutes  
**What you'll get:** a working `library(cellucid)`

Use this page when you can’t install or load `cellucid` in R.

## Symptom: `library(cellucid)` fails (“there is no package called 'cellucid'”)

**Likely causes**
- You didn’t install the package yet.
- You installed it in a different library than the one your current R session uses.

**How to confirm**
```r
.libPaths()
installed.packages()[,"Package"]
```

**Fix**
- Reinstall from GitHub:
  ```r
  install.packages("remotes")
  remotes::install_github("theislab/cellucid-r")
  ```
- Restart R and try again.

## Symptom: GitHub install fails (HTTP 401/403, rate limit)

**Likely causes**
- Corporate proxy / blocked GitHub
- GitHub API rate limiting

**Fix**
- Try `pak`:
  ```r
  install.packages("pak")
  pak::pak("theislab/cellucid-r")
  ```
- Configure a GitHub PAT:
  ```r
  Sys.setenv(GITHUB_PAT = "...")
  ```

## Symptom: “package ‘jsonlite’ is not available”

**Likely causes**
- No CRAN mirror configured.
- Offline/proxy issues.

**Fix**
```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("jsonlite")
```

## Symptom: “The Matrix package is required to store connectivities.”

You passed `connectivities` as a base R matrix without `Matrix` installed.
(`Matrix` is also what you need to hold a sparse `gene_expression` in the first
place.)

**Fix**
```r
install.packages("Matrix")
```

## Symptom: SeuratObject does not install on R older than 4.4

SeuratObject 5 needs Matrix 1.6-4 or later, and CRAN's current Matrix release
needs R 4.4, so the installer reports a dependency it cannot install.

**Fix**
```r
remotes::install_version("Matrix", version = "1.6-5")
install.packages("SeuratObject")
```

## Symptom: installing from GitHub fails with a compiler error

**Likely cause**
- The package compiles a few small C files (checksums, the outlier quantiles,
  one-pass matrix checks, the decoders of a saved session) and links zlib, and
  no C toolchain, or no zlib header (`zlib.h: No such file or directory`), is
  installed.

**Fix**
- Windows: install Rtools. macOS: `xcode-select --install`. Linux: a C
  compiler plus the R and zlib development headers (for example `r-base-dev`,
  which brings `zlib1g-dev`, on Debian and Ubuntu; `zlib-devel` on Fedora).
- Or install a platform binary from CRAN once the index lists the package; see
  {doc}`../installation`.

## Symptom: “data is a Seurat object from the SeuratObject package, which is not installed”

`cellucid_prepare(data = )` and `cellucid_view()` read an object with its own
package, loaded only when such an object arrives. Install the one the message
names: `install.packages("SeuratObject")`, or
`BiocManager::install("SingleCellExperiment")` (and `"SummarizedExperiment"`).
An object saved by Seurat v5 needs SeuratObject 5.

Without any of them, start from matrices and data frames:
- {doc}`../e_integrations_recipes/03_raw_matrices_and_data_frames_recipe`

## Symptom: “The viewer needs the httpuv package”

`httpuv`, `later` and `curl` are imports, installed with the package; this
message means the library holding them is not the one R loads. Reinstall them
with `install.packages(c("httpuv", "later", "curl"))`.

## Symptom: “I updated but R still loads an old version”

**How to confirm**
```r
packageVersion("cellucid")
find.package("cellucid")
```

**Fix**
- Restart R.
- Reinstall with:
  ```r
  remotes::install_github("theislab/cellucid-r", force = TRUE)
  ```
