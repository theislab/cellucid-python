# Host Exports for Sharing

**Audience:** anyone sharing datasets (lab leads, computational users, educators)  
**Time:** 15–60 minutes (depends on upload size)  
**Goal:** make your exported dataset easy for other people to load.

There are two broad sharing modes:

1) **Public sharing (no server):** publish exports to a public GitHub repo (recommended for demos/tutorials).
2) **Private sharing (server mode):** run a server on a trusted machine (recommended for sensitive data).

```{warning}
Publishing exports to a public GitHub repo is publishing your data. Treat it like making the dataset public.
```

## Option 1 — Public sharing via GitHub (recommended)

Cellucid supports loading datasets from a GitHub repo if you provide an **exports root**:

```
exports/
  datasets.json
  my_dataset/
    dataset_identity.json
    points_2d.bin.gz
    obs_manifest.json
    obs/
    ...
```

### Step 1: Put each dataset in its own folder

If you exported to `out_dir = ".../pbmc_demo"`, you can treat that as the dataset folder.

Recommended layout:

```
exports/
  pbmc_demo/
    <contents written by cellucid_prepare()>
```

### Step 2: Create `datasets.json`

`datasets.json` is a manifest that lets the web app discover datasets in the exports root.

Format (version 1):

```json
{
  "version": 1,
  "default": "pbmc_demo",
  "datasets": [
    {"id": "pbmc_demo", "path": "pbmc_demo/", "name": "PBMC demo", "n_cells": 10000, "n_genes": 2000}
  ]
}
```

#### Option A: Generate `datasets.json` in R (no Python required)

```r
library(jsonlite)

`%||%` <- function(a, b) if (!is.null(a)) a else b

exports_root <- file.path(getwd(), "exports")
dataset_dirs <- list.dirs(exports_root, full.names = TRUE, recursive = FALSE)

datasets <- list()
for (dir in dataset_dirs) {
  ident_path <- file.path(dir, "dataset_identity.json")
  if (!file.exists(ident_path)) next

  ident <- jsonlite::read_json(ident_path, simplifyVector = TRUE)

  datasets[[length(datasets) + 1L]] <- list(
    id = ident$id %||% basename(dir),
    path = paste0(basename(dir), "/"),
    name = ident$name %||% basename(dir),
    n_cells = ident$stats$n_cells %||% NULL,
    n_genes = ident$stats$n_genes %||% NULL
  )
}

stopifnot(length(datasets) > 0)

default_id <- datasets[[1]]$id

manifest <- list(version = 1, default = default_id, datasets = datasets)
jsonlite::write_json(manifest, file.path(exports_root, "datasets.json"), auto_unbox = TRUE, pretty = TRUE)
```

#### Option B: Generate `datasets.json` using `cellucid-python` (if you have Python)

```python
from cellucid.prepare_data import generate_datasets_manifest
generate_datasets_manifest("./exports", default_dataset="pbmc_demo")
```

The Python helper is documented in:
- {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`

### Step 3: Push `exports/` to GitHub

Create a public GitHub repo and commit:
- `exports/datasets.json`
- `exports/<dataset_folder>/...`

Then, in the Cellucid web app, use the GitHub loading option and point it at your repo/path.

The authoritative step-by-step is here:
- {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`

## Option 2 — Private sharing via server mode

If your dataset is private (patient IDs, unpublished experiments, etc.), do not publish to GitHub.

Instead:
- run a server on a trusted machine, and
- let collaborators access it over VPN/SSH tunnel.

Cellucid’s server mode is currently provided by the Python tooling:
- {doc}`../../web_app/b_data_loading/04_server_tutorial`

## Troubleshooting pointers

- “GitHub connect says `datasets.json not found`” → you uploaded only the dataset folder; you need an exports root + datasets.json.
- “My repo is too big for GitHub” → export fewer genes / quantize / compress, or use server mode.
- “CORS errors when hosting on a custom server” → use the supported server mode workflow or configure CORS headers.
