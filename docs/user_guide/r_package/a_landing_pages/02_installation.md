# Installation

**Audience:** everyone  
**Time:** 5–15 minutes (depending on your setup)

## Requirements

- **R:** `cellucid` currently targets **R ≥ 4.3.0**.
- **Hard dependency:** `jsonlite`
- **Optional but recommended:** `Matrix` (sparse matrices + connectivities)

```{note}
`cellucid-r` is designed to be dependency-light and does not require Seurat or SingleCellExperiment.
Those packages are only needed if you want to follow the Seurat/SCE recipes in this guide.
```

## Install from GitHub (recommended for now)

### Option A: `remotes` (simple, common)

```r
install.packages("remotes")
remotes::install_github("theislab/cellucid-r")
```

### Option B: `pak` (often faster, better resolver)

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
- A version like `0.99.0` (or similar) from `packageVersion`

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
- {doc}`../e_integrations_recipes/01_seurat_recipe`
- {doc}`../e_integrations_recipes/02_singlecellexperiment_recipe`

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

- Ready to export? Go to {doc}`04_quick_start_3_levels`.
