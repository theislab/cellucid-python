# Browser File Picker (No Python Setup)

This tutorial covers the **browser-only** ways to load data into Cellucid.

You will use the Cellucid web UI to pick one of:
- a **pre-exported folder** (recommended)
- a **`.h5ad` file** (works, but has hard limits)
- a **`.zarr` directory** (often better than `.h5ad` in the browser)

If your dataset is large or you want maximum reliability/performance, use **Server Mode** ({doc}`04_server_tutorial`).

## At A Glance

**Audience**
- Wet lab / beginner: follow the click-by-click path; you do not need Python.
- Computational users: pay attention to the `.h5ad` vs `.zarr` limitations and when to pre-export.

**Time**
- Export folder picker: ~2–5 minutes
- `.h5ad` / `.zarr` picker: ~5–10 minutes

**Prerequisites**
- A modern desktop browser (Chrome/Edge/Firefox recommended)
- A dataset in one of these forms:
  - exported folder from `prepare()`
  - `.h5ad`
  - `.zarr` directory

**Privacy model**
- File picker modes load data from your computer into your browser.
- Your files are **not uploaded** to Cellucid servers.

## Fast Path (Wet Lab Friendly)

1) Open Cellucid: https://www.cellucid.com
2) In the left sidebar, find the data loading area (often labeled **Dataset Connections** or **Browse local data…**).
3) Choose the button matching what you have:
   - **Folder** (pre-exported folder) — recommended
   - **.h5ad** — small datasets only
   - **.zarr** — often better than `.h5ad`
4) Wait for the dataset to load.
5) Confirm success:
   - you see points rendered in the canvas
   - a field selector/legend has content

If you don’t see points after loading, jump to the troubleshooting section at the end.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-file-picker-entry-point
Suggested filename: data_loading/06_file-picker-buttons.png
Where it appears: Data Loading → 03_browser_file_picker_tutorial → Fast Path
Capture:
  - UI location: left sidebar → the three file picker buttons (folder / h5ad / zarr)
  - State prerequisites: Cellucid open; no dataset loaded
  - Action to reach state: open Cellucid in a new tab
Crop:
  - Include: the file picker buttons + any short explanatory text
  - Exclude: browser chrome and personal info
Redact:
  - Remove: any private dataset names
Annotations:
  - Callouts: #1 folder picker, #2 h5ad picker, #3 zarr picker
Alt text:
  - Data loading panel with buttons for folder, h5ad, and zarr.
Caption:
  - Tell readers which button to click depending on their file type.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the browser file picker buttons.
:width: 100%

Use the folder/h5ad/zarr buttons to load data directly from your computer.
```

## Option #3 — Load a Pre-exported Folder (Recommended)

### When to use this
- You (or a collaborator) already ran `cellucid.prepare(...)` in Python.
- You want the **fastest** and **most reliable** browser experience.

### Why this is recommended
Exported folders are designed for the viewer:
- embeddings and obs fields are compact
- gene expression is stored in a way Cellucid can fetch on demand
- the browser avoids loading a massive monolithic file

### What you click
Click the **Folder** / **Browse folder…** button and select the export directory.

### What success looks like
- Points appear quickly.
- Field lists are populated.
- Gene search responds quickly (first gene may take a moment).

**Vector fields (optional)**
- If the export includes vector fields (e.g. RNA velocity), you’ll be able to enable the overlay after loading.
- If the overlay toggle is disabled or the dropdown is empty, it usually means the dataset has no vectors for the current dimension (2D vs 3D) or they weren’t exported.
- See {doc}`../i_vector_field_velocity/index`.

<!-- SCREENSHOT PLACEHOLDER
ID: data-loading-export-folder-success
Suggested filename: data_loading/07_export-folder-success.png
Where it appears: Data Loading → 03_browser_file_picker_tutorial → Option #3
Capture:
  - UI location: left sidebar + canvas
  - State prerequisites: exported folder loaded successfully
  - Action to reach state: select an exported folder with the folder picker
Crop:
  - Include: dataset name, point count (if shown), one visible embedding, and the field selector
Redact:
  - Remove: private dataset ID/name if needed
Alt text:
  - Loaded embedding with a populated field selector.
Caption:
  - Describe the minimum indicators that confirm a successful load.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for a successful exported-folder load.
:width: 100%

After selecting an exported folder, you should see points and a populated field selector.
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
cellucid serve /path/to/data.h5ad
```

Then open:

```text
https://www.cellucid.com?remote=http://localhost:8765
```

### `.h5ad` minimum requirements

Your `.h5ad` must include at least one embedding:

- `obsm['X_umap_3d']` with shape `(n_cells, 3)` (recommended)
- `obsm['X_umap_2d']` with shape `(n_cells, 2)`
- `obsm['X_umap']` with shape `(n_cells, 2 or 3)`

Optional (but highly recommended):
- `obs` columns for coloring and filtering
- expression matrix `X` (dense or sparse) for gene coloring
- `obsp['connectivities']` if you want connectivity visualization

Optional (if you want the vector field / velocity overlay):
- Per-cell vectors in `obsm` using the naming convention `<field>_umap_<dim>d` (e.g. `velocity_umap_2d`, `velocity_umap_3d`).
- The overlay dropdown shows fields available for the current dimension only.

## Option #5 — Load a `.zarr` Directory Directly

### When to use this
- You have a `.zarr` export of AnnData and want browser-only viewing.
- You want better “chunked” behavior than `.h5ad`.

### Practical expectations
- The browser may still need to read a fair amount of metadata up front.
- Gene expression can often be fetched more lazily (chunk-by-chunk).

If your `.zarr` is extremely large, server mode is still the most reliable.

Vector fields are supported in `.zarr` as well (same convention as `.h5ad`): store them in `obsm` as `velocity_umap_2d`, `T_fwd_umap_3d`, etc.

## Common Failure Modes (and Why They Happen)

### “It worked for a demo dataset but not my data”
- Your file is missing required embeddings.
- Your `.h5ad` is too large for browser memory.
- Your `.zarr` directory is incomplete (missing `.zgroup`, `.zattrs`, etc.).

### “The file picker won’t let me select a folder”
- Some browsers restrict directory selection.
- Workaround: use server mode.

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

## Edge Cases (Browser-Specific)

- **Safari limitations**: directory picking and File System Access APIs vary by version.
- **Permission prompts**: if you deny folder access, Cellucid cannot read files.
- **Large `.h5ad`**: the browser loads the whole file; memory spikes are normal.
- **Nested `.zarr`**: some tools create nested stores; ensure you select the correct root.

## Troubleshooting (Massive)

Use this like a checklist. Most issues are diagnosable in < 2 minutes.

---

### Symptom: “I clicked Folder, but nothing happens”

**Likely causes (ordered)**
1) Your browser blocked the directory picker.
2) The UI is open in an embedded context that disallows directory access.
3) You clicked the wrong button (e.g. `.h5ad` picker for an export folder).

**How to confirm**
- Try a different browser (Chrome is the most reliable for folder access).
- Try selecting a different folder (a very small test export).

**Fix**
- Switch to Chrome/Edge.
- If you cannot use folder picking, run server mode instead ({doc}`04_server_tutorial`).

---

### Symptom: “`.h5ad` loads forever / browser freezes”

**Likely causes**
- The `.h5ad` file is too large for browser memory.

**How to confirm**
- Check file size on disk.
- Watch the browser tab’s memory usage.

**Fix**
- Use server mode: `cellucid serve data.h5ad`.
- Or export once with `prepare()` and load the export folder.

---

### Symptom: “It says no embedding / no UMAP”

**Likely causes**
- Missing `obsm['X_umap']` / `X_umap_2d` / `X_umap_3d`.

**How to confirm**
- In Python:

  ```python
  print(adata.obsm.keys())
  ```

**Fix**
- Compute UMAP and store it under one of the supported keys.

---

### Symptom: “Zarr directory selected, but it errors immediately”

**Likely causes**
- The directory is not a valid AnnData zarr store.
- Missing `.zgroup` / `.zattrs`.

**Fix**
- Re-export with `adata.write_zarr("data.zarr")`.
- Ensure you selected the top-level `.zarr` directory.

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
