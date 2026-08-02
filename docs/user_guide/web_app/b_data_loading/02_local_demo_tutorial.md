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
For a complete, inspectable public repository with three tiny datasets, use
{doc}`11_custom_dataset_repository`.

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
from __future__ import annotations

from pathlib import Path

from cellucid import prepare

export_dir = Path("./exports/pbmc_demo")

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    X_umap_1d=adata.obsm["X_umap_1d"],
    X_umap_2d=adata.obsm["X_umap_2d"],
    X_umap_3d=adata.obsm["X_umap_3d"],
    vector_fields={
        "velocity_umap_1d": adata.obsm["velocity_umap_1d"],
        "velocity_umap_2d": adata.obsm["velocity_umap_2d"],
        "velocity_umap_3d": adata.obsm["velocity_umap_3d"],
    },
    vector_field_default="velocity_umap",
    out_dir=export_dir,
    dataset_id="pbmc_demo",
    dataset_name="PBMC demo",
    dataset_description="PBMC scRNA-seq with aligned velocity",
    obs_categorical_dtype="uint16",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)
```

This complete call assumes `adata` contains every named matrix. The executable
{doc}`../../python_package/f_notebooks_tutorials/prepare_pancreas` notebook
deliberately exports only the checked-in source file's real 2-D UMAP and
connectivity; it does not invent 1-D/3-D coordinates or velocity. For the full
standard Pancreas artifact with deterministic 1-D, 2-D, and 3-D embeddings and
aligned velocity, use {doc}`10_standard_pancreas_dataset` and the reproducible
{doc}`../../python_package/f_notebooks_tutorials/33_vector_fields_and_velocity_overlay_end_to_end`
workflow. For a dataset without a particular optional product, omit that
product deliberately—for example, omit both `X_umap_3d` and
`velocity_umap_3d` rather than inventing a third coordinate or passing `None`.

### Export knobs you should understand (even as a beginner)

- `dataset_id`:
  - A **stable** identifier that will appear in URLs and manifests.
  - Use 1–180 ASCII characters: start with a letter or digit, then use only
    letters, digits, `.`, `_`, or `-`; do not end with `.` or use a reserved
    Windows device name.
  - Avoid changing it after you publish.
  - See {doc}`06_dataset_identity_why_it_matters` for what breaks when it changes.

- `dataset_name` / `dataset_description`:
  - Human-friendly labels shown in the UI. `dataset_name` must be non-empty
    and unpadded, without control characters; Unicode is supported.

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
  - See {doc}`../i_vector_field_velocity/index` for naming and UI controls.

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

- **Fastest validation:** open Cellucid and use the **Prepared** file picker (Option #3).
- **Most realistic validation:** run the CLI server (Option #6) and open the
  exact printed prepared-data Viewer URL
  (`http://127.0.0.1:<port>/?source=remote`).

This avoids debugging “GitHub problems” that are actually export problems.

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: Cellucid web app with the sidebar open and a single-cell embedding colored by cell type.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
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
- **Branches**:
  - Shorthand without a branch always means `main`; Cellucid does not probe
    other branch names.
  - If your exports live on a different branch, name it explicitly:
    - `owner/repo@my-branch/exports`
    - or `https://github.com/owner/repo/tree/my-branch/exports`
- **File size limits**:
  - GitHub blocks ordinary Git files larger than 100 MiB.
  - Git LFS raw repository content is a pointer, not a transparent replacement
    for the prepared bytes Cellucid requests.

If your export is too large for GitHub:
- Use **Server Mode** instead ({doc}`04_server_tutorial`).

### 4.4 What repo path to enter in Cellucid

Expand **Session**, type the path into the **GitHub data:** field, and press
**Connect**. Enter one of these:

- `owner/repo` (if `datasets.json` is at repo root)
- `owner/repo/exports` (if your exports root is in a folder)
- `owner/repo@my-branch/exports` (if your exports root is on a custom branch)

```{figure} ../../../_static/screenshots/data_loading/connect-github-catalog.png
:alt: Close-up of the GitHub data control with an owner/repo/exports path typed into the text field and a mouse pointer resting on the Connect button beside it.
:width: 496px

The **GitHub data:** field takes the exports root, not a full URL and not a
single dataset directory. Press **Connect** beside it.
```

Cellucid will then:
- fetch and validate `datasets.json`
- fetch and validate every listed `dataset_identity.json`
- show the validated dataset list
- fetch the selected dataset's remaining files on demand

**Shareable URL**

Once it works, share the exact dataset. The `dataset=` field is required when
the catalog contains multiple datasets:

```text
https://www.cellucid.com/?github=owner/repo/exports&dataset=pbmc_demo
```

(Replace `owner/repo/exports` with your chosen path.)

## Optional: “Local Demo” (Run the Web App Locally With Demo Exports)

The web app’s `local-demo` source loads from an **exports base URL** (it does not have to be inside the web app repo).

In production, Cellucid’s sample datasets are intended to live in a separate repository/site (e.g. `cellucid-datasets`) and the web app is configured via:
- `<meta name="cellucid-exports-base-url" content="...">` in `cellucid/index.html`, or
- `?exportsBaseUrl=...` as a runtime override.

High-level local-dev workflow:
- Run the web app locally.
- Point `exportsBaseUrl` at any static host that serves `exports/datasets.json` + dataset folders.

```{note}
If your exports host is a different origin, it must authorize the Cellucid
origin with CORS headers so the browser can fetch the JSON and binary files.
There is no alternate browser-side transport for a host that rejects the
cross-origin requests.
```

## Edge Cases

- **Changed dataset IDs after publishing**: links break; collaborators load the wrong dataset.
- **Missing `datasets.json`**: GitHub loader cannot discover datasets.
- **Wrong repo path**: the loader fetches `datasets.json` from the wrong folder (404).
- **Wrong branch**: your exports root is on a different branch than you think (pin the branch in the repo path).
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
- Read the second notification. Cellucid prints the exact address it asked for,
  in the form:

  ```text
  Resource not found: https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

- Open that URL in a browser tab. If it 404s there too, the path or branch is
  wrong; if it loads there but not in Cellucid, suspect the network.

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: Two stacked Cellucid notifications, the upper naming the raw.githubusercontent.com URL that returned not found and the lower explaining that the GitHub repository could not be reached and suggesting a check for a typo.
:width: 776px

What a wrong repository path looks like. The upper notification hands you the
URL to test by hand.
```

**Fix**
- Regenerate `datasets.json` with an explicit catalog default:

  ```python
  generate_datasets_manifest(
      "./exports",
      default_dataset="pbmc_demo",
  )
  ```

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
- Want a complete repository to copy and adapt? → {doc}`11_custom_dataset_repository`
- Want the most reliable workflow for big datasets? → {doc}`04_server_tutorial`
- Export format expectations: {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- If GitHub loading fails: {doc}`08_troubleshooting_data_loading`
