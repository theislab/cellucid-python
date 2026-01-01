# `prepare()` / export overview

**Audience:** everyone (wet lab / non-technical, computational, and developers)  
**Time:** 15–30 minutes to get a working export; 60+ minutes for the full deep dive  
**Goal:** export a dataset folder that the Cellucid web app loads quickly and reproducibly

`cellucid.prepare(...)` writes an **export folder** containing:
- binary data files (points, fields, genes, edges, vectors),
- compact JSON manifests describing those files,
- and a required `dataset_identity.json` describing what the web app should expect.

That folder can then be loaded in the **Cellucid web app** via:
- the browser file picker,
- `cellucid serve ...` (local server),
- static hosting (e.g., GitHub Pages / any HTTP host),
- or other supported loading paths (see the Web App loading docs).

```{important}
Exports are for **fast, reproducible visualization**. They are not intended to be a lossless archival format for analysis.
Keep your original AnnData / source pipeline as the ground truth.
```

## What you’ll learn

- What `prepare()` writes (and how the web app uses it)
- What is **required** vs **optional** input
- How to choose **export** vs **direct viewing** (AnnData server / Jupyter)
- What “success looks like” (files on disk + what you should see in the web app)
- The biggest footguns (row order, NaNs, huge gene exports, “it didn’t update”)

## Prerequisites

- `pip install cellucid`
- You have either:
  - an AnnData (`adata`) with embeddings/obs/(optional) X/var/connectivities, or
  - the equivalent raw arrays + pandas DataFrames
- A writable output directory (prefer an SSD for large exports)

---

## Fast path (wet lab / non-technical)

This is the “I just want to open my dataset and explore it” path.

### Step 1 — Get an export folder

Ask a computational collaborator (or your pipeline) to create an export folder using `cellucid.prepare(...)`.

What you should receive is a directory that contains at least:
- `dataset_identity.json`
- `points_2d.bin` or `points_3d.bin` (optionally with `.gz`)

### Step 2 — Load the folder in the Cellucid web app

Use the browser file picker workflow:
- Web App docs: {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`

<!-- SCREENSHOT PLACEHOLDER
ID: prepare-overview-webapp-file-picker
Where it appears: User Guide → Python Package → Data Preparation API → prepare/export overview
Capture:
  - Open Cellucid web app
  - Navigate to the dataset loading panel
  - Choose the “folder picker / local export folder” option
  - Show the OS folder picker with an exported folder highlighted (do NOT select private directories)
Crop:
  - Include: the loader UI + the picker action button OR the OS picker dialog with the export folder selected
  - Exclude: personal usernames/home paths if sensitive
Redact:
  - Remove: any private dataset names, patient IDs, repo URLs, account info
Annotations:
  - Callouts: (1) the loader action/button, (2) the export folder name, (3) the “Select/Choose” confirmation button
Alt text:
  - Cellucid dataset loader showing the option to pick an exported dataset folder from the local computer.
Caption:
  - Load a pre-exported folder created by `cellucid.prepare(...)` using the browser file picker.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for loading an export folder via the web app file picker.
:width: 100%

Load a pre-exported folder created by `cellucid.prepare(...)` using the browser file picker.
```

### Step 3 — Confirm “success”

You should see:
- points rendered (2D or 3D),
- metadata fields available to color/filter (if exported),
- optional: gene search available (if gene expression was exported),
- optional: vector overlay controls available (if vector fields were exported).

If it fails, jump to:
- {doc}`11_troubleshooting_prepare_export`
- Web app loading troubleshooting: {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`

---

## Practical path (computational users)

### When to use export vs direct viewing

Use **exports** when you want:
- the smoothest browser interaction (especially for large `n_cells`),
- deterministic, shareable folders,
- static hosting (no Python server required once exported),
- vector overlays from precomputed per-cell vectors.

Use **direct viewing** when you want:
- zero export step while iterating,
- lazy gene expression loading from `.h5ad`/`.zarr` via server mode,
- notebook-first workflows where the dataset changes often.

Start here for non-export viewing paths:
- {doc}`../d_viewing_apis/index`

### Minimal call structure (inputs → outputs)

At minimum, you provide:
- `latent_space`: `(n_cells, n_latent_dims)` (required)
- `obs`: a DataFrame with `n_cells` rows (required)
- at least one embedding: `X_umap_2d` or `X_umap_3d` (required)

Optionally, you also provide:
- `gene_expression` + `var` (for gene overlays/search)
- `connectivities` (for KNN graph edges)
- `vector_fields` (for velocity/displacement overlays)

And you always choose:
- `out_dir` (where files are written)
- optional performance knobs (`compression`, `var_quantization`, `obs_continuous_quantization`)

### A complete AnnData-based example (recommended)

```python
from cellucid import prepare

# Required: at least one embedding
X_umap_3d = adata.obsm.get("X_umap_3d")
X_umap_2d = adata.obsm.get("X_umap_2d")
X_umap = adata.obsm.get("X_umap")  # common key (2D or 3D)

if X_umap_3d is None and X_umap is not None and X_umap.shape[1] == 3:
    X_umap_3d = X_umap
if X_umap_2d is None and X_umap is not None and X_umap.shape[1] == 2:
    X_umap_2d = X_umap

# Required: latent_space (used for outlier quantiles on categorical fields)
latent = adata.obsm.get("X_pca", X_umap_3d if X_umap_3d is not None else X_umap_2d)
if latent is None:
    raise ValueError("Provide adata.obsm['X_pca'] or an embedding to use as latent_space.")

# Optional: vector fields (velocity/drift) in embedding space
vector_fields = {
    "velocity_umap_2d": adata.obsm.get("velocity_umap_2d"),
    "velocity_umap_3d": adata.obsm.get("velocity_umap_3d"),
}

prepare(
    latent_space=latent,
    obs=adata.obs,

    # Optional but common (enables gene search/overlays)
    var=adata.var,
    gene_expression=adata.X,

    # Optional
    connectivities=adata.obsp.get("connectivities"),
    vector_fields=vector_fields,

    # Embeddings (provide any/all of 1D/2D/3D)
    X_umap_2d=X_umap_2d,
    X_umap_3d=X_umap_3d,

    # Output + performance knobs
    out_dir="./exports/pbmc_demo",
    dataset_name="PBMC demo",
    dataset_id="pbmc_demo",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,

    # Iteration safety
    force=True,
)
```

### What files get written (at a glance)

Minimum:

```text
out_dir/
├── dataset_identity.json
└── points_2d.bin(.gz)  or  points_3d.bin(.gz)
```

Typical full export:

```text
out_dir/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # if gene expression exported
├── connectivity_manifest.json        # if connectivities exported
├── points_1d.bin(.gz)                # if provided
├── points_2d.bin(.gz)                # if provided
├── points_3d.bin(.gz)                # if provided
├── obs/                              # field binaries (continuous + categorical)
├── var/                              # gene expression binaries (one per gene)
├── connectivity/                     # edge binaries
└── vectors/                          # vector field binaries
```

Full spec: {doc}`09_output_format_specification_exports_directory`

### “What success looks like” (export-time)

On a successful run, `prepare()` prints an “Export Settings” summary and then logs that it wrote:
- `points_<dim>d.bin(.gz)` for each embedding you provided
- `obs_manifest.json` + `obs/` files
- optional: `var_manifest.json` + `var/` files
- optional: `connectivity_manifest.json` + `connectivity/` files
- optional: `vectors/` files (and vector metadata in `dataset_identity.json`)
- always: `dataset_identity.json`

### How to open the export after it’s written

- In the browser (no Python running): {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`
- With a local server: run `cellucid serve ./exports/pbmc_demo` and open the printed URL
- In Jupyter (embedded): see {doc}`../api/jupyter` and {doc}`../d_viewing_apis/index`

---

## Deep path (expert / developer)

### Key invariants (the “contract”)

1) **Cell identity is row order.** There is no separate cell id table.
2) **All cell-aligned inputs must match `n_cells`.**
3) **Embeddings are normalized per-dimension** to a stable `[-1, 1]` coordinate range.
4) **Categorical fields produce two payloads**:
   - category codes (uint8/uint16 with a reserved missing marker),
   - per-cell outlier quantiles (float32 or quantized) computed from `latent_space`.
5) **Quantization is lossy** (by design) and meant for visualization performance.

Start with the global contract details:
- {doc}`02_input_requirements_global`

### Exports and the “14 loading paths”

The exported folder format is the common denominator across:
- local folder picking,
- static hosting (including GitHub exports roots),
- server-backed loading,
- and some notebook workflows.

Web app context:
- {doc}`../../web_app/b_data_loading/01_loading_options_overview`

### Versioning and compatibility

- `dataset_identity.json["version"]` is the identity schema version (current exports write `2`).
- `obs_manifest.json["_format"]` and `var_manifest.json["_format"]` are the manifest encoding versions (current exports write `compact_v1`).
- `dataset_identity.json["cellucid_data_version"]` records the `cellucid` Python package version used to export.

---

## Edge cases and common footguns

- **Row-order mismatch** (most common): embeddings/obs/gene_expression were not aligned after filtering/concatenation.
- **NaN/Inf in embeddings**: normalization will produce NaN outputs → viewer may fail or render nothing.
- **“It didn’t update”**: existing files were skipped because `force=False` (or you reused an `out_dir`).
- **Gene export explosion**: exporting many genes writes *one file per gene* and can create huge folders + slow file systems.
- **Too many categories**: categorical fields with very large category counts are unusable in the UI and may exceed dtype limits.
- **Vector dims don’t match**: providing a 3D vector field but only exporting 2D points (or vice versa) causes the field to be skipped.

## Troubleshooting (overview)

- Deep troubleshooting: {doc}`11_troubleshooting_prepare_export`
- Web app loading issues: {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`

## Next steps

- Read the global invariants and preflight checks: {doc}`02_input_requirements_global`
