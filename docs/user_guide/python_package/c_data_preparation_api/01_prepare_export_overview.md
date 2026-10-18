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
- `obs_manifest.json`
- an `embeddings/` folder holding at least one coordinate file, such as
  `embeddings/0_2d.bin` (optionally with `.gz`)

The export is in format 3, so it needs a viewer release that reads format 3.

### Step 2 — Load the folder in the Cellucid web app

Use the browser file picker workflow:
- Web App docs: {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 524px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

### Step 3 — Confirm “success”

You should see:
- points rendered (1D, 2D, or 3D),
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
- browser requests for gene expression on demand; `.h5ad` remains
  read-only-backed while `.zarr` is materialized eagerly in Python,
- notebook-first workflows where the dataset changes often.

Start here for non-export viewing paths:
- {doc}`../d_viewing_apis/index`

### Minimal call structure (inputs → outputs)

At minimum, you provide:
- `latent_space`: `(n_cells, n_latent_dims)` (required)
- `obs`: a DataFrame with `n_cells` rows (required)
- `embeddings`: at least one embedding, keyed by id, such as
  `embeddings={"X_umap": adata.obsm["X_umap"]}` (required)

Optionally, you also provide:
- `gene_expression` + `var` (for gene overlays/search)
- `connectivities` (for KNN graph edges)
- `vector_fields` (for velocity/displacement overlays)

And you always choose:
- `out_dir` (where files are written)
- optional performance knobs (`compression`, `var_quantization`, `obs_continuous_quantization`)

### A complete AnnData-based example (recommended)

```python
from cellucid import Embedding, prepare

# Required: at least one embedding. Every one you pass is exported, keyed by
# its id; a wide array is plotted at the columns you name.
embeddings = {"X_umap": Embedding(adata.obsm["X_umap"], label="UMAP", source_key="X_umap")}
if "X_tsne" in adata.obsm:
    embeddings["X_tsne"] = Embedding(adata.obsm["X_tsne"], label="t-SNE", source_key="X_tsne")
if "X_umap_3d" in adata.obsm:
    embeddings["X_umap_3d"] = adata.obsm["X_umap_3d"]

# Optional: vector fields (velocity/drift), keyed by field id and then by the
# embedding each array is drawn over.
vector_fields = {}
if "velocity_umap" in adata.obsm:
    vector_fields["velocity"] = {"X_umap": adata.obsm["velocity_umap"]}

prepare(
    # Required: latent_space (used for outlier quantiles on categorical fields)
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,

    # Optional but common (enables gene search/overlays)
    var=adata.var,
    gene_expression=adata.X,

    # Optional
    connectivities=adata.obsp.get("connectivities"),
    vector_fields=vector_fields or None,

    # Embeddings, and the one shown first when a dimension has several
    embeddings=embeddings,
    embedding_defaults={"2d": "X_umap"},

    # Output + performance knobs
    out_dir="./exports/pbmc_demo",
    dataset_name="PBMC demo",
    dataset_id="pbmc_demo",
    # Optional: fix this only in a canonical byte-reproducible builder.
    created_at="2026-07-26T12:34:56Z",
    obs_categorical_dtype="uint16",
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
├── obs_manifest.json
├── obs/
└── embeddings/                       # one file per embedding, e.g. 0_2d.bin(.gz)
```

Typical full export:

```text
out_dir/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # if gene expression exported
├── connectivity_manifest.json        # if connectivities exported
├── embeddings/                       # one file per embedding you provided
├── obs/                              # field binaries (continuous + categorical)
├── var/                              # gene expression binaries (one per gene)
├── connectivity/                     # edge binaries
└── vectors/                          # vector field binaries
```

Full spec: {doc}`09_output_format_specification_exports_directory`

### “What success looks like” (export-time)

On a successful run, `prepare()` prints an “Export Settings” summary and then logs that it wrote:
- `embeddings/<index>_<dim>d.bin(.gz)` for each embedding you provided
- `obs_manifest.json` + `obs/` files
- optional: `var_manifest.json` + `var/` files
- optional: `connectivity_manifest.json` + `connectivity/` files
- optional: `vectors/` files (and vector metadata in `dataset_identity.json`)
- always: `dataset_identity.json`

### How to open the export after it’s written

- In the browser (no Python running): {doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`
- With a local server: run `cellucid serve ./exports/pbmc_demo` and open the printed URL
- In Jupyter (embedded): see {doc}`../g_api_reference_coverage/api/jupyter` and {doc}`../d_viewing_apis/index`

---

## Deep path (expert / developer)

### Key invariants (the “contract”)

1) **Cell identity is row order.** There is no separate cell id table.
2) **All cell-aligned inputs must match `n_cells`.**
3) **Each embedding is normalized on its own** to a stable `[-1, 1]` coordinate range.
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

- `dataset_identity.json["version"]` is the identity schema version (current exports write `3`).
- `obs_manifest.json["_format"]` and `var_manifest.json["_format"]` are the manifest encoding versions (current exports write `compact_v2`).
- There is no compatibility path: a version-2 export (with `points_<dim>d.bin`
  files at its root) must be regenerated.
- `dataset_identity.json["cellucid_data_version"]` records the `cellucid` Python package version used to export.

---

## Edge cases and common footguns

- **Row-order mismatch** (most common): embeddings/obs/gene_expression were not aligned after filtering/concatenation.
- **NaN/Inf in embeddings**: the candidate is rejected before output
  publication.
- **“It didn’t update”**: a non-empty target with `force=False` raises
  `FileExistsError`; use `force=True` only when you intend a complete
  replacement.
- **Gene export explosion**: exporting many genes writes *one file per gene* and can create huge folders + slow file systems.
- **Too many categories**: categorical fields with very large category counts are unusable in the UI and may exceed dtype limits.
- **Vector dims don’t match**: the candidate is rejected; every vector array
  must have the dimension of the embedding it is bound to.

## Troubleshooting (overview)

- Deep troubleshooting: {doc}`11_troubleshooting_prepare_export`
- Web app loading issues: {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`

## Next steps

- Read the global invariants and preflight checks: {doc}`02_input_requirements_global`
