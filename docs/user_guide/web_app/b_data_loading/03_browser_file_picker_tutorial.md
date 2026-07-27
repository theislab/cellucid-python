# Browser File Picker (No Python Setup)

This tutorial covers the **browser-only** ways to load data into Cellucid.

You will use the Cellucid web UI to pick one of:
- a **prepared export folder**
- a **`.h5ad` file** (works, but has hard limits)
- a portable **`.zarr.zip` or `.zip` archive** containing one Zarr v2 store

If your dataset is large or you want maximum reliability/performance, use **Server Mode** ({doc}`04_server_tutorial`).

## At A Glance

**Audience**
- Wet lab / beginner: follow the click-by-click path; you do not need Python.
- Computational users: pay attention to the `.h5ad` vs `.zarr` limitations and when to pre-export.

**Time**
- Prepared folder picker: ~2–5 minutes
- `.h5ad` / Zarr ZIP picker: ~5–10 minutes

**Prerequisites**
- A current stable desktop release of Chrome, Edge, Firefox, or Safari with
  WebGL2 and `DecompressionStream('gzip')`
- A dataset in one of these forms:
  - exported folder from `prepare()`
  - `.h5ad`
  - `.zarr.zip` or `.zip` containing one complete AnnData Zarr v2 store

**Privacy model**
- File picker modes load data from your computer into your browser.
- Your files are **not uploaded** to Cellucid servers.

## Fast Path (Wet Lab Friendly)

1) Open [Cellucid](https://www.cellucid.com/).
2) Expand **Session** in the left sidebar and find **Local data:**.
3) Choose the button matching exactly what you have:
   - **Prepared** — a folder created by `prepare()`
   - **H5AD** — one current-schema `.h5ad`, no larger than 512 MiB
   - **Zarr ZIP** — one `.zarr.zip` or `.zip` archive
4) Wait for the dataset to load.
5) Confirm success:
   - you see points rendered in the canvas
   - a field selector/legend has content

If you don’t see points after loading, jump to the troubleshooting section at the end.

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 246px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

## Option #3 — Load a Prepared Export Folder

### When to use this
- You (or a collaborator) already ran `cellucid.prepare(...)` in Python.
- You want the **fastest** and **most reliable** browser experience.

### Why this is recommended
Exported folders are designed for the viewer:
- embeddings and obs fields are compact
- gene expression is stored in a way Cellucid can fetch on demand
- the browser avoids loading a massive monolithic file

### What you click
Click **Prepared** and select the export directory.

### What success looks like
- Points appear quickly.
- Field lists are populated.
- Gene search responds quickly (first gene may take a moment).

**Vector fields (optional)**
- If the export includes vector fields (e.g. RNA velocity), you’ll be able to enable the overlay after loading.
- If the overlay toggle is disabled or the dropdown is empty, it usually means the dataset has no vectors for the current dimension (2D vs 3D) or they weren’t exported.
- See {doc}`../i_vector_field_velocity/index`.

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: Cellucid web app with the sidebar open and a single-cell embedding colored by cell type.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

## Option #4 — Load a `.h5ad` File Directly (Quick Preview)

### When to use this
- You have a **small** `.h5ad` and just want a quick look.
- You do not want to run Python locally.

### Important performance limitation (do not skip)

Browser `.h5ad` loading is **not truly lazy**.

Due to browser limitations, the viewer must load the entire `.h5ad` file into memory
before it can read data.

Practical consequence:
- large `.h5ad` files can freeze the tab or crash the browser

### Recommended alternative for larger `.h5ad`
Use server mode (recommended):

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

Then open the exact **Viewer URL** printed in the terminal. With the default
port it is:

```text
http://127.0.0.1:8765/?anndata=true
```

### `.h5ad` minimum requirements

Your `.h5ad` must include at least one embedding:

- `obsm['X_umap_1d']` with shape `(n_cells, 1)`
- `obsm['X_umap_2d']` with shape `(n_cells, 2)`
- `obsm['X_umap_3d']` with shape `(n_cells, 3)`

The key is part of the data contract: Cellucid does not infer the dimension
from an unsuffixed key or from the array shape.

Optional (but highly recommended):
- `obs` columns for coloring and filtering
- expression matrix `X` (dense or sparse) for gene coloring
- `obsp['connectivities']` if you want connectivity visualization

Optional (if you want the vector field / velocity overlay):
- Per-cell vectors in `obsm` using the naming convention `<field>_umap_<dim>d` (e.g. `velocity_umap_2d`, `velocity_umap_3d`).
- The overlay dropdown shows fields available for the current dimension only.

## Option #5 — Load a Zarr ZIP Archive Directly

### When to use this
- You have a complete AnnData Zarr v2 store packaged as `.zarr.zip` or `.zip`.
- You want a portable single-file selection in any supported browser.

### Practical expectations
- Cellucid validates the ZIP directory and Zarr metadata before adopting the
  dataset.
- Gene expression is read by Zarr chunk on demand.
- Archive entries, decoded chunks, and the active embedding still consume
  browser memory and are bounded by the documented reader limits.

For an extremely large store, prefer a prepared Cellucid export. Python server
mode loads `.zarr` eagerly, so use it only when the store fits server memory.

The archive may contain the Zarr store at its root or inside exactly one root
directory. It must preserve dotfiles such as `.zgroup`, `.zattrs`, and
`.zarray`. Vector fields use the same exact `obsm` convention as H5AD:
`velocity_umap_2d`, `T_fwd_umap_3d`, and other
`<field>_umap_<dim>d` keys.

## Common Failure Modes (and Why They Happen)

### “It worked for a demo dataset but not my data”
- Your file is missing required embeddings.
- Your `.h5ad` is too large for browser memory.
- Your Zarr ZIP is incomplete, contains multiple roots, or omits required
  `.zgroup`, `.zattrs`, or `.zarray` entries.

### “The file picker won’t let me select a folder”
- The **Prepared** control accepts directories only; H5AD and Zarr ZIP accept
  one file.
- In an embedded or managed browser, directory access may be disabled by the
  embedding policy. Open the standalone Cellucid page or choose the explicit
  Python server workflow.

### “It loads but fields are empty”
- `obs` is empty or not written correctly.
- For `.h5ad`/`.zarr`, field names may not be where you expect.

### “Gene search exists but everything is zero / blank”
- `adata.X` is empty or contains unexpected values.
- Or you are loading a dataset that lacks gene expression entirely.

### “Vector field overlay toggle is disabled / dropdown is empty”
- Your dataset does not contain any vector fields.
- Vector fields exist but are not named using the expected `*_umap_2d` / `*_umap_3d` convention.
- You’re in 3D but only 2D vectors exist (or vice versa).

## Browser and file-selection boundaries

- **Prepared** uses the directory-input capability available in supported
  current desktop browsers.
- **H5AD** and **Zarr ZIP** use ordinary single-file inputs and have the same
  validation contract in Chrome, Edge, Firefox, and Safari.
- **Permission prompts**: if you deny folder access, Cellucid cannot read files.
- **Large `.h5ad`**: the browser loads the whole file; memory spikes are normal.
- **Zarr ZIP root**: the archive must contain either root metadata at archive
  root or exactly one wrapper directory around the Zarr store.

## Troubleshooting (Massive)

Use this like a checklist. Most issues are diagnosable in < 2 minutes.

---

### Symptom: “I clicked Prepared, but nothing happens”

**Likely causes (ordered)**
1) You denied directory access or a managed-browser policy blocked it.
2) The UI is open in an embedded context that disallows directory access.
3) You clicked H5AD or Zarr ZIP instead of **Prepared**.

**How to confirm**
- Open Cellucid as a standalone page and click **Prepared**.
- Try selecting a different folder (a very small test export).

**Fix**
- Allow directory access and use the standalone page.
- If directory input is disabled by organizational policy, choose server mode
  explicitly ({doc}`04_server_tutorial`).

---

### Symptom: “`.h5ad` loads forever / browser freezes”

**Likely causes**
- The `.h5ad` file is too large for browser memory.

**How to confirm**
- Check file size on disk.
- Watch the browser tab’s memory usage.

**Fix**
- Use server mode:
  `cellucid serve data.h5ad --dataset-name "My dataset" --dataset-id my-dataset`.
- Or export once with `prepare()` and load the export folder.

---

### Symptom: “It says no embedding / no UMAP”

**Likely causes**
- None of the exact supported keys—`obsm['X_umap_1d']`,
  `obsm['X_umap_2d']`, or `obsm['X_umap_3d']`—is present.

**How to confirm**
- In Python:

  ```python
  print(adata.obsm.keys())
  ```

**Fix**
- Store each embedding under the key matching its exact column count. For
  example, a two-column UMAP belongs at `obsm['X_umap_2d']`.
- Compute UMAP and store it under one of the supported keys.

---

### Symptom: “Zarr ZIP selected, but it errors immediately”

**Likely causes**
- The archive does not contain exactly one valid AnnData Zarr v2 store.
- Required `.zgroup`, `.zattrs`, or `.zarray` entries are missing.
- The ZIP is encrypted, multi-disk, uses an unsupported compression method, or
  exceeds a reader limit.

**Fix**
- With the supported AnnData version, re-export with
  `adata.write_zarr("data.zarr")`; the resulting store must contain the Zarr
  v2 root markers `.zgroup` and `.zattrs`.
- Package the complete store while preserving its dotfiles, then select the
  resulting `.zarr.zip` or `.zip` file with **Zarr ZIP**.

---

### Symptom: “Fields list is empty, but I know I have metadata”

**Likely causes**
- `obs` was not saved in your file, or uses unexpected encodings.

**How to confirm**
- In Python: `print(adata.obs.head())`

**Fix**
- Save a clean copy of your AnnData (or export via `prepare()`).

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset truly has no vector fields.
2) Vector fields exist, but the naming convention is wrong (Cellucid can’t discover them).
3) Dimension mismatch: vectors exist for 2D but you’re viewing 3D (or vice versa).

**How to confirm**
- In Python, list `obsm` keys and look for entries like `velocity_umap_2d`:

  ```python
  print(adata.obsm.keys())
  ```

**Fix**
- Rename/regenerate vector fields to follow the `*_umap_2d` / `*_umap_3d` convention.
- Switch the viewer to the dimension that has vectors.
- If you’re using exports, re-export with `prepare(..., vector_fields={...})`.

For overlay UI behavior and deeper debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

## Next Steps

- Format expectations (exports vs `.h5ad` vs `.zarr`): {doc}`07_folder_file_format_expectations_high_level_link_to_spec`
- Full troubleshooting matrix: {doc}`08_troubleshooting_data_loading`
- If browser-only loading is too slow or unreliable for your dataset → {doc}`04_server_tutorial`
- If you are working in notebooks and want programmatic interaction → {doc}`05_jupyter_tutorial`
