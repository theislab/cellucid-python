# Local & Remote Demo (Share Without Running a Server)

This tutorial shows how to **export a dataset once** and then view/share it **without** running a long-lived Python server.

You’ll learn two closely-related workflows:

1) **Remote demo (recommended): public GitHub repo**
   - You push an `exports/` folder to a **public** GitHub repo.
   - Anyone can open Cellucid and load your dataset from that repo.

2) **Local demo (optional): run the web app locally**
   - You put your exported dataset inside a local copy of the Cellucid web app.
   - You run a simple static file server (or just a local host).
   - This is useful for offline demos or if you want to host Cellucid yourself.

If you just want to load your own local files right now, use {doc}`03_browser_file_picker_tutorial`.

## At A Glance

**Audience**
- Wet lab / beginner: follow the step-by-step checklist; you do not need to understand file formats.
- Computational users: pay attention to dataset IDs, compression/quantization, and GitHub limits.
- Developers: this explains the exact directory structure (`datasets.json` + `dataset_identity.json`) used by the frontend.

**Time**
- Minimum viable share (one dataset, public repo): ~15–30 minutes (mostly uploading)
- Adding multiple datasets + polish: ~30–60 minutes

**Prerequisites**
- Python environment with `cellucid` installed
- A dataset in one of these forms:
  - AnnData object (`adata`)
  - a `.h5ad` file
  - a `.zarr` directory
- (For remote demo) a GitHub account + willingness to publish a **public** repo

## Important: Privacy and Sharing

When you publish exports to a public GitHub repo, you are effectively publishing your dataset.

Before you do this:
- Verify there is **no private metadata** in `obs` (e.g. patient IDs, sample identifiers, internal notes).
- Consider removing/renaming columns that should not be public.
- Consider publishing only a reduced dataset or a subset of fields.

If you need **private** sharing:
- Use **Server Mode** ({doc}`04_server_tutorial`) on a private machine + VPN/SSH tunnel.
- Or use the **Browser File Picker** on each collaborator’s machine.

## The Required Folder Layout (Do This Exactly)

Cellucid’s GitHub loader expects a top-level `datasets.json` **manifest** plus one folder per dataset.

Recommended layout:

```text
exports/                     # "exports root" (this is what you point Cellucid at)
├── datasets.json            # required
├── pbmc_demo/               # dataset folder (name can differ from dataset_id)
│   ├── dataset_identity.json
│   ├── points_2d.bin.gz
│   ├── points_3d.bin.gz
│   ├── obs_manifest.json
│   ├── var_manifest.json
│   ├── obs/
│   ├── var/
│   ├── vectors/             # optional (vector fields like velocity/drift)
│   └── connectivity/        # optional
└── another_dataset/
    └── ...
```

**Key rules**
- `datasets.json` is at the **exports root**, not inside the dataset folder.
- Each dataset folder must contain a `dataset_identity.json`.
- All paths in `datasets.json` are **relative** to the exports root.

## Step 1 — Export Your Dataset (Create the Dataset Folder)

You can export from:
- an in-memory `AnnData` (`adata`)
- a `.h5ad` file
- a `.zarr` directory

The central idea:
- You create one dataset folder, e.g. `./exports/pbmc_demo/`.
- That folder is what Cellucid will actually fetch.

```python
# Option A: export from an in-memory AnnData object
#
# This is the most common workflow for computational users.

from __future__ import annotations

from pathlib import Path

# NOTE: in your real workflow, you will already have `adata`.
# Here we only show the structure.

from cellucid import prepare

export_dir = Path("./exports/pbmc_demo")

# Optional: vector fields for the overlay (RNA velocity, drift, etc.).
# These must be in embedding coordinates and match the same cell order as your points.
# Naming convention: <field>_umap_<dim>d (e.g., velocity_umap_2d, velocity_umap_3d).
#
# vector_fields = {
#     "velocity_umap_2d": adata.obsm.get("velocity_umap_2d"),
#     "velocity_umap_3d": adata.obsm.get("velocity_umap_3d"),
# }

# Example sketch (fill in your own AnnData parts):
# X_umap = adata.obsm["X_umap"]
# prepare(
#     latent_space=adata.obsm.get("X_pca", X_umap),
#     obs=adata.obs,
#     var=adata.var,
#     gene_expression=adata.X,
#     connectivities=adata.obsp.get("connectivities"),
#     X_umap_2d=adata.obsm.get("X_umap_2d", X_umap if X_umap.shape[1] == 2 else None),
#     X_umap_3d=adata.obsm.get("X_umap_3d", X_umap if X_umap.shape[1] == 3 else None),
#     vector_fields=vector_fields,          # optional: velocity/drift overlays
#     out_dir=export_dir,
#     dataset_id="pbmc_demo",              # stable ID (important for sharing)
#     dataset_name="PBMC demo",            # human-readable
#     dataset_description="PBMC scRNA-seq (example)",
#     compression=6,                        # gzip level (tradeoff: size vs CPU)
#     var_quantization=8,                   # smaller, lossy; good default
#     obs_continuous_quantization=8,
# )
```

### Export knobs you should understand (even as a beginner)

- `dataset_id`:
  - A **stable** identifier that will appear in URLs and manifests.
  - Avoid changing it after you publish.
  - See {doc}`06_dataset_identity_why_it_matters` for what breaks when it changes.

- `dataset_name` / `dataset_description`:
  - Human-friendly labels shown in the UI.

- `compression` (gzip level):
  - Higher = smaller files, slower export.
  - A good starting point is `compression=6`.

- `var_quantization` and `obs_continuous_quantization`:
  - Controls how values are compressed/quantized.
  - Smaller bit-depth makes exports much smaller but loses precision.
  - For visualization, `8` is often a good default.

- `vector_fields` (optional):
  - Per-cell displacement vectors visualized by the vector field / velocity overlay.
  - Must be provided per-dimension (`*_umap_2d`, `*_umap_3d`) and match the same cell order as the embedding.
  - If you export them, you’ll see a vector-field dropdown in the UI after loading.
  - See {doc}`../i_vector_field_velocity/index` and the naming conventions in `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`.

If you want the full, parameter-by-parameter specification, see the Python package guide section:
- {doc}`../../python_package/c_data_preparation_api/index`

## Step 2 — Generate `datasets.json` (Create the Exports Manifest)

Once you have one (or many) dataset folders under `./exports/`, generate `datasets.json`.

This file is required for the GitHub loader.

```python
from __future__ import annotations

from cellucid.prepare_data import generate_datasets_manifest

# This scans ./exports/*/dataset_identity.json
# and writes ./exports/datasets.json
generate_datasets_manifest("./exports", default_dataset="pbmc_demo")
```

## Step 3 — Validate Locally (Before You Upload)

Before publishing, verify that your export loads:

- **Fastest validation:** open Cellucid and use the **Folder** file picker (Option #3).
- **Most realistic validation:** run the CLI server (Option #6) and open `?remote=...`.

This avoids debugging “GitHub problems” that are actually export problems.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-github-validate-local-file-picker
Suggested filename: data_loading/04_github-validate-export-folder.png
Where it appears: Data Loading → 02_local_demo_tutorial → Step 3 — Validate Locally
Capture:
  - UI location: left sidebar → file picker
  - State prerequisites: you have an exports folder on disk
  - Action to reach state: click “Browse local data…” and select the exported folder
Crop:
  - Include: the picker entry point + the loaded dataset name in the UI
  - Exclude: your home directory path in the picker UI (privacy)
Redact:
  - Remove: any sensitive dataset names
Alt text:
  - Cellucid showing a successfully loaded exported folder.
Caption:
  - Explain that this is how you validate exports before publishing.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for validating an export using the browser folder picker.
:width: 100%

Validate your export locally with the folder picker before publishing it to GitHub.
```

## Step 4 — Publish to a Public GitHub Repo (Remote Demo)

### 4.1 Create a repository

You can use either:
- a repo that contains only data (recommended), e.g. `my-lab/cellucid-exports`
- a general repo that contains code + an `exports/` folder

### 4.2 Commit the exports root

The repo must include:
- `exports/datasets.json`
- at least one dataset folder under `exports/`

### 4.3 GitHub constraints (do not skip)

- **Public repo only**: Cellucid’s GitHub loader fetches data via `raw.githubusercontent.com`.
- **Branch names**: the loader recognizes only these branch names if you include a branch in the path:
  - `main`, `master`, `develop`, `dev`, `gh-pages`
- **File size limits**:
  - GitHub blocks individual files > 100 MB in normal git.
  - Git LFS is often not usable for direct “raw file” fetch workflows.

If your export is too large for GitHub:
- Use **Server Mode** instead ({doc}`04_server_tutorial`).

### 4.4 What repo path to enter in Cellucid

In the Cellucid UI (GitHub connection), you will enter one of these:

- `owner/repo` (if `datasets.json` is at repo root)
- `owner/repo/exports` (if your exports root is in a folder)
- `owner/repo/gh-pages/exports` (if using `gh-pages` branch)

Cellucid will then:
- fetch `datasets.json`
- show the dataset list
- fetch only the files it needs on demand

**Shareable URL**

Once it works, you can share a direct link:

```text
https://www.cellucid.com?github=owner/repo/exports
```

(Replace `owner/repo/exports` with your chosen path.)

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-github-connect
Suggested filename: data_loading/05_github-connect.png
Where it appears: Data Loading → 02_local_demo_tutorial → Step 4.4
Capture:
  - UI location: left sidebar → GitHub repo connection
  - State prerequisites: none
  - Action to reach state: open Cellucid, locate GitHub repo input
Crop:
  - Include: GitHub repo input + Connect button + success message (if possible)
  - Exclude: any private owner/repo names
Redact:
  - Replace: repo path with a public example like `theislab/cellucid/exports`
Annotations:
  - Callouts: #1 the repo path format, #2 Connect button, #3 dataset dropdown after connect
Alt text:
  - GitHub repo connection input with a repository path entered.
Caption:
  - Explain the exact format expected.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the GitHub repo connection input.
:width: 100%

Connect Cellucid to a public GitHub repo by entering `owner/repo/path-to-exports`.
```

## Optional: “Local Demo” (Run the Web App Locally With Your Exports)

This is useful if you want:
- an offline demo
- a self-hosted instance of Cellucid
- a “known-good” environment where you control the exact viewer version

**High-level idea**
- You clone the Cellucid web app repo.
- You copy your exported dataset folder(s) into the web app’s `assets/exports/`.
- You ensure `assets/exports/datasets.json` is updated.
- You run a simple static server.

Because this mixes repositories (web + python), we keep it high-level here.
If you want, we can expand this into a dedicated, fully scripted tutorial.

## Edge Cases

- **Changed dataset IDs after publishing**: links break; collaborators load the wrong dataset.
- **Missing `datasets.json`**: GitHub loader cannot discover datasets.
- **Wrong repo path**: the loader fetches `datasets.json` from the wrong folder (404).
- **Non-standard branch name**: the GitHub loader currently recognizes only common branch names.
- **Large files**: GitHub may reject pushes, or raw fetch may be blocked by a corporate proxy.
- **Vector fields “missing” after publish**:
  - You didn’t export them (no `vectors/` directory and no `vector_fields` block in `dataset_identity.json`).
  - Or you exported only 2D vectors but you are viewing the dataset in 3D (dimension mismatch).

## Troubleshooting (Large, Practical)

Use this section like a lookup table.

---

### Symptom: “Connected to GitHub, but it says `datasets.json not found`”

**Likely causes (ordered)**
1) Your repo path points at the wrong folder.
2) You uploaded only a single dataset folder and forgot `datasets.json`.
3) Your exports root is in a different branch than you think.

**How to confirm**
- Open this URL in a browser (replace values):

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

- If that URL 404s, Cellucid will also fail.

**Fix**
- Regenerate `datasets.json` using `generate_datasets_manifest("./exports")`.
- Ensure the file is at the exports root you are pointing Cellucid at.

**Prevention**
- Keep a single exports root folder named `exports/` in your repo.
- Always validate the raw URL before sending the link to collaborators.

---

### Symptom: “Datasets list shows up, but loading a dataset fails (404 for a file)”

**Likely causes**
- `datasets.json` lists a dataset `path` that does not exist.
- You renamed the dataset folder after generating the manifest.

**How to confirm**
- Compare `datasets.json["datasets"][i]["path"]` to the actual repo folder names.

**Fix**
- Re-run `generate_datasets_manifest` after renaming folders.

---

### Symptom: “It loads, but gene search is extremely slow”

**Likely causes**
- Export not compressed/quantized (files too large)
- Using `.h5ad` browser mode (not an export)
- Huge dataset where each gene vector is large even after compression

**Fix**
- Re-export with `compression=6` and `var_quantization=8`.
- Prefer server mode for very large datasets.

---

### Symptom: “The vector field overlay toggle is disabled / no fields appear”

**Likely causes**
1) You didn’t export any vector fields (common).
2) Vector fields exist, but the naming convention is wrong (so they aren’t detected).
3) You exported only 2D vectors but you are looking at a 3D embedding (or vice versa).

**How to confirm**
- Open `exports/<dataset>/dataset_identity.json` and check for a `vector_fields` block.
- Confirm the repo contains a `vectors/` folder with files like:

  ```text
  vectors/<fieldId>_2d.bin.gz
  vectors/<fieldId>_3d.bin.gz
  ```

**Fix**
- Re-export with `prepare(..., vector_fields={...})` using the `*_umap_2d` / `*_umap_3d` naming convention.
- Make sure you commit/publish the `vectors/` directory.
- Switch the viewer to the dimension that actually has vectors.

For UI behavior and overlay-specific debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “GitHub blocks my push / files are too big”

**Likely causes**
- Individual files > 100 MB
- Total repo size becomes too large

**Fix**
- Reduce export size (quantization, fewer genes/fields).
- Switch to server mode.
- Host exports elsewhere and load via a server.

## Next Steps

- Want to load from your computer (no publishing)? → {doc}`03_browser_file_picker_tutorial`
- Want the most reliable workflow for big datasets? → {doc}`04_server_tutorial`
- Export format expectations: {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- If GitHub loading fails: {doc}`08_troubleshooting_data_loading`
