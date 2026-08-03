# Folder / file format expectations (high-level; link to spec)

**Audience:** everyone (beginners get “what should I have?”, experts get exact file/key expectations)  
**Time:** 10–15 minutes  
**What you’ll learn:** what Cellucid expects when you load an export folder, a
GitHub exports repo, an H5AD file, a browser Zarr ZIP, or a Python-side Zarr
directory
**Prerequisites:** none (helpful: you know which loading option you plan to use)

---

## The four things Cellucid can load (practical view)

Cellucid can load your data through different *paths*, but they all boil down to one of these **inputs**:

1) **Pre-exported dataset folder** (recommended)  
2) **GitHub “exports root”** (a folder that contains `datasets.json` + multiple exported datasets)  
3) **AnnData `.h5ad`**  
4) **AnnData Zarr v2 store** — one portable ZIP in the browser UI, or the
   directory itself in Python server/Jupyter workflows

Everything else (browser file picker, server mode, Jupyter) is *how you point the viewer at one of the above*.

---

## Exported dataset folder (what must be inside)

An exported dataset folder is created by Python `cellucid.prepare(...)` or R
`cellucid_prepare(...)`.

Minimum expected files:

```text
my_export/
├── dataset_identity.json          # required
├── obs_manifest.json              # required
└── points_2d.bin(.gz)             # at least one non-empty 1d/2d/3d file
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
- Payload files are named by an **integer index**, never by a field name. Read
  the paths out of the manifests rather than assembling them from a label.

### Every binary payload is little-endian

Multi-byte values in `.bin` / `.bin.gz` payloads are stored least-significant
byte first, on every platform, by contract — not by whatever the exporting
machine happened to use. Both writers pin it: the Python exporter converts the
payload before writing, and `cellucid_prepare()` in R passes
`endian = "little"` at every `writeBin`. An export produced on a big-endian
machine is therefore byte-identical to one produced anywhere else, and the
browser can read it without inspecting a byte-order flag.

You only need to care about this if you are writing an exporter or reading the
payloads yourself. The exact rule, and what it means for each dtype, is in
{doc}`../../python_package/c_data_preparation_api/09_output_format_specification_exports_directory`.

See also:
- Export API reference + high-level format: {doc}`../../python_package/g_api_reference_coverage/api/export`
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
- The catalog's only top-level keys are `version`, `default`, and `datasets`.
  Version is exactly `1`, and `default` names an exact listed ID.
- Each dataset entry requires unique `id` and `path` fields. `path` is a safe
  relative directory ending in `/`.
- An entry may also contain `name`, `description`, `n_cells`, and `n_genes`;
  when present, each value must exactly match `dataset_identity.json`.
- Connecting validates the catalog and every listed identity before presenting
  the collection. Remaining payloads are fetched for the selected dataset as
  needed.

See {doc}`11_custom_dataset_repository` for the complete reference repository,
exact UI value, schema, validation sequence, and share links. The shorter
export-first workflow is {doc}`02_local_demo_tutorial`.

---

## AnnData H5AD and Zarr (what Cellucid reads)

The browser reads one H5AD file or one ZIP containing a complete Zarr v2 store.
Server mode and Jupyter accept H5AD paths, Zarr directories, or an in-memory
AnnData object. Python opens H5AD read-only-backed and materializes Zarr
eagerly with `anndata.read_zarr`.

### Minimum requirements (all AnnData modes)

You must have at least one embedding in `obsm`:
- `obsm["X_umap_1d"]` with shape `(n_cells, 1)`
- `obsm["X_umap_2d"]` with shape `(n_cells, 2)`
- `obsm["X_umap_3d"]` with shape `(n_cells, 3)`

A suffixed key and its column count must agree exactly, and a key that names its
dimension always decides that dimension. When a file declares none of the three
and carries the plain `X_umap` that `sc.tl.umap()` writes, that array is read at
the dimension its own column count states — 1, 2, or 3 columns — so an ordinary
Scanpy object loads as written. Any other width is refused. This holds for every
AnnData mode: the browser file picker, `cellucid serve`, and the notebook
viewers.

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
vectors/<index>_<dim>d.bin
vectors/<index>_<dim>d.bin.gz        # if export compression is enabled
```

`<index>` is the vector field's position, not its name. Payload files across
the whole export are named by an integer index and never by a field name; the
manifest beside them is what says which field an index belongs to. Always read
the paths out of `dataset_identity.json["vector_fields"]` rather than assembling
them from a field label.

Typical identity schema (high-level):

```json
{
  "vector_fields": {
    "default_field": "velocity_umap",
    "fields": {
      "velocity_umap": {
        "available_dimensions": [2, 3],
        "files": {
          "2d": "vectors/0_2d.bin.gz",
          "3d": "vectors/0_3d.bin.gz"
        }
      }
    }
  }
}
```

#### How to confirm vector fields exist (quick checks)

- **Exports**: open `<export_dir>/dataset_identity.json` and search for `vector_fields`; also confirm `vectors/` exists.
- **AnnData**: print `adata.obsm.keys()` and look for keys like `velocity_umap_2d`.

Naming and schema details:
- {doc}`../i_vector_field_velocity/01_what_vector_fields_are_user_facing`

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
- There is no gene expression (`X`), or the genes were named from a `var` column whose vocabulary nobody is searching in.

Fix:
- Provide `adata.X` and name the genes from the `var` column your readers will type.
- In Jupyter/server mode, pass `gene_id_column=...` if needed.

---

## Edge cases and limits (read if you’re debugging weirdness)

- **NaN/Inf in embeddings**: invalid current input; readers reject it before
  dataset adoption.
- **Duplicate gene names**: invalid current-format input; regenerate the
  dataset with one distinct name per gene. Names are never filenames, so any
  printable text is acceptable — but two genes sharing a name are
  indistinguishable in the field selector, the legend, and every exported
  figure.
- **Huge categorical fields**: legends/palettes don’t scale to 50k–100k categories; consider filtering/aggregating.
- **Browser byte ceilings**: the browser paths refuse anything that would need
  more than **512 MiB** for a single working set — an `.h5ad` larger than that on
  disk, a declared prepared payload larger than that, or a gzip stream that says
  it decompresses to more. JSON metadata has a separate, smaller ceiling of
  64 MiB. None of these apply to server mode or Jupyter, where Python holds the
  data and the browser receives only what it asks for.
- **Byte order**: never a compatibility problem. Prepared payloads are
  little-endian by contract; H5AD is byte-order-agnostic, because HDF5 converts
  to host order while reading; Zarr records its own byte order per array and the
  reader swaps when it has to.

---

## Next steps

- Pick a workflow: {doc}`01_loading_options_overview`
- Load locally without Python: {doc}`03_browser_file_picker_tutorial`
- Load big `.h5ad`/`.zarr`: {doc}`04_server_tutorial`
- Publish a public prepared catalog: {doc}`11_custom_dataset_repository`
- Debug loading failures: {doc}`08_troubleshooting_data_loading`
