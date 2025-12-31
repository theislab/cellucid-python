# Installation and Dependency Issues

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

## Symptom: “Matrix package is required …”

You’ll see this if you try to export:
- sparse gene expression matrices, or
- connectivities,
without Matrix installed.

**Fix**
```r
install.packages("Matrix")
```

## Symptom: “Seurat / SingleCellExperiment not installed”

These are optional and only required for the recipe/tutorial pages.

If you only want the minimal workflow, you can export from matrices/data.frames:
- {doc}`../e_integrations_recipes/03_raw_matrices_and_data_frames_recipe`

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
