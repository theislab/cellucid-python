# Folder / file format expectations (high-level; link to spec)

**Audience:** everyone (beginners get “what should I have?”, experts get exact file/key expectations)  
**Time:** 10–15 minutes  
**What you’ll learn:** what Cellucid expects when you load an export folder, a GitHub exports repo, a `.h5ad`, or a `.zarr` store  
**Prerequisites:** none (helpful: you know which loading option you plan to use)

---

## The four things Cellucid can load (practical view)

Cellucid can load your data through different *paths*, but they all boil down to one of these **inputs**:

1) **Pre-exported dataset folder** (recommended)  
2) **GitHub “exports root”** (a folder that contains `datasets.json` + multiple exported datasets)  
3) **AnnData `.h5ad`**  
4) **AnnData `.zarr` directory**

Everything else (browser file picker, server mode, Jupyter) is *how you point the viewer at one of the above*.

---

## Exported dataset folder (what must be inside)

An exported dataset folder is created by `cellucid.prepare(...)`.

Minimum expected files:

```text
my_export/
├── dataset_identity.json          # required
└── points_2d.bin(.gz)             # at least one of 1d/2d/3d must exist
    points_3d.bin(.gz)             # optional (if you exported 3D)
```

Typical full layout (recommended for most datasets):

```text
my_export/
├── dataset_identity.json
├── obs_manifest.json
├── var_manifest.json                 # optional (if gene expression exported)
├── connectivity_manifest.json        # optional (if connectivities exported)
├── points_1d.bin(.gz)                # optional
├── points_2d.bin(.gz)                # optional
├── points_3d.bin(.gz)                # optional
├── obs/
├── var/
├── connectivity/
└── vectors/                          # optional (vector fields)
```

Notes:
- The `.gz` suffix appears when you export with compression enabled (gzip).
- The viewer reads `dataset_identity.json` first; it describes what dimensions/files exist.
- If you load a folder and Cellucid says “not a valid export”, it usually means `dataset_identity.json` is missing or malformed.

See also:
- Export API reference + high-level format: {doc}`../../python_package/api/export`
- Why identity matters: {doc}`06_dataset_identity_why_it_matters`

---

## GitHub exports root (multi-dataset sharing)

For GitHub-hosted sharing, you do **not** point Cellucid at a single dataset folder.

You point it at an **exports root** that contains a `datasets.json` manifest and one folder per dataset:

```text
exports/                     # this is what you point Cellucid at
├── datasets.json            # required for GitHub loader
├── pbmc_demo/               # dataset folder
│   ├── dataset_identity.json
│   └── ...
└── another_dataset/
    └── ...
```

Important distinctions:
- `datasets.json` is required for the GitHub loader, but **not** required for local folder picking.
- Dataset folder name and dataset id are different concepts:
  - folder name: `pbmc_demo/` (path on disk / in repo)
  - dataset id: `dataset_identity.json["id"]` (identity inside Cellucid)

See {doc}`02_local_demo_tutorial` for the exact publishing workflow.

---

## AnnData `.h5ad` and `.zarr` (what Cellucid reads)

When you load `.h5ad`/`.zarr` (directly in the browser, via server mode, or in Jupyter), Cellucid reads the dataset as **AnnData**.

### Minimum requirements (all AnnData modes)

You must have at least one embedding in `obsm`:
- `obsm["X_umap_3d"]` with shape `(n_cells, 3)` (recommended)
- `obsm["X_umap_2d"]` with shape `(n_cells, 2)`
- `obsm["X_umap"]` with shape `(n_cells, 2)` or `(n_cells, 3)`

### Optional (but commonly expected)

- **Obs fields** (metadata for coloring/filtering): `adata.obs[...]`
- **Gene metadata** (for gene search / naming): `adata.var[...]` and/or `adata.var.index`
- **Gene expression**: `adata.X` (dense or sparse)
- **Connectivity**: `adata.obsp["connectivities"]` (KNN graph)

### Vector fields (optional overlays)

Vector fields are optional per-cell displacement vectors (e.g. RNA velocity, drift) that Cellucid can render as an animated overlay.

If you expect vector fields, the key requirements are:
- **Correct dimension**: 2D vectors for 2D embeddings, 3D vectors for 3D, etc.
- **Correct row order**: vectors must match the same cell order as points/embeddings.
- **Discoverable naming / metadata**: so the viewer can find the fields.

See the overlay UI docs:
- {doc}`../i_vector_field_velocity/index`

#### AnnData naming convention (`obsm`)

Cellucid discovers vector fields from `adata.obsm` keys using the convention:

- `<field>_umap_<dim>d`

Examples:
- `velocity_umap_2d`, `velocity_umap_3d`
- `T_fwd_umap_2d`, `T_bwd_umap_2d`
- `drift_umap_3d`

Notes:
- Keys starting with `X_` are reserved for embeddings and ignored as vector fields.
- The overlay dropdown shows fields available for the **current** dimension only.

#### Export format (`vectors/` + `dataset_identity.json`)

Exports store vectors as binary files under `vectors/` and describe them in `dataset_identity.json["vector_fields"]`.

Typical filenames:

```text
vectors/<fieldId>_<dim>d.bin
vectors/<fieldId>_<dim>d.bin.gz      # if export compression is enabled
```

Typical identity schema (high-level):

```json
{
  "vector_fields": {
    "default_field": "velocity_umap",
    "fields": {
      "velocity_umap": {
        "available_dimensions": [2, 3],
        "files": {
          "2d": "vectors/velocity_umap_2d.bin.gz",
          "3d": "vectors/velocity_umap_3d.bin.gz"
        }
      }
    }
  }
}
```

#### How to confirm vector fields exist (quick checks)

- **Exports**: open `<export_dir>/dataset_identity.json` and search for `vector_fields`; also confirm `vectors/` exists.
- **AnnData**: print `adata.obsm.keys()` and look for keys like `velocity_umap_2d`.

Internal reference (naming and schema details):
- `cellucid/markdown/VECTOR_FIELD_OVERLAY_CONVENTIONS.md`

---

## Common “format mismatch” problems (and what they mean)

### “No embedding / no UMAP”

Meaning:
- Cellucid could not find a supported embedding key in `obsm`, or the shape is wrong.

Fix:
- Ensure you have one of the supported keys and it has shape `(n_cells, 2)` or `(n_cells, 3)`.

### “Fields list is empty”

Meaning:
- `obs` is empty, or fields were not exported correctly.

Fix:
- For exports: ensure `obs_manifest.json` exists and `dataset_identity.json["obs_fields"]` lists fields.
- For AnnData: ensure `adata.obs` actually contains columns and is not all missing.

### “Gene search returns nothing”

Meaning:
- There is no gene expression (`X`) and/or gene ids are not where Cellucid expects.

Fix:
- Provide `adata.X` and a usable gene id column.
- In Jupyter/server mode, pass `gene_id_column=...` if needed.

---

## Edge cases and limits (read if you’re debugging weirdness)

- **NaN/Inf in embeddings**: can produce missing points, huge bounding boxes, or WebGL crashes.
- **Duplicate IDs** (cells or genes): can cause confusing “wrong gene” or selection behavior.
- **Huge categorical fields**: legends/palettes don’t scale to 50k–100k categories; consider filtering/aggregating.
- **Browser `.h5ad`**: for large files, this often fails due to memory limits (prefer server mode or exports).

---

## Next steps

- Pick a workflow: {doc}`01_loading_options_overview`
- Load locally without Python: {doc}`03_browser_file_picker_tutorial`
- Load big `.h5ad`/`.zarr`: {doc}`04_server_tutorial`
- Debug loading failures: {doc}`08_troubleshooting_data_loading`
