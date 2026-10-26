# Browser File Picker (No Python Setup)

This tutorial covers the **browser-only** ways to load data into Cellucid.

You will use the Cellucid web UI to pick one of:
- a **`.h5ad` file**, of any size
- an **AnnData Zarr store folder** (`.zarr`, Zarr v2 or v3), as `cellucid.prepare()`
  and `cellucid_prepare()` write it by default
- a portable **`.zarr.zip` or `.zip` archive** containing one Zarr store (v2 or v3)

The same buttons open a Seurat or SingleCellExperiment `.rds`, a Loom file, a
10x Cell Ranger `.h5` or `outs/` folder, a MatrixMarket folder and a MuData
(`.h5mu`, MuData Zarr): {doc}`17_open_r_objects`, {doc}`18_open_loom_and_10x`
and {doc}`19_multimodal_mudata` say what each holds.

If your dataset lives on another machine, use **Server Mode** ({doc}`04_server_tutorial`): it hands the file to the same reader in the browser.

## At A Glance

**Audience**
- Wet lab / beginner: follow the click-by-click path; you do not need Python.
- Computational users: pay attention to how the expression matrix is stored, and when to write a Cellucid-ready AnnData with `prepare()`.

**Time**
- ~2–5 minutes

**Prerequisites**
- A current stable desktop release of Chrome, Edge, Firefox, or Safari with
  WebGL2 and `DecompressionStream('gzip')`
- A dataset in one of these forms:
  - `.h5ad`
  - a `.zarr` store folder
  - `.zarr.zip` or `.zip` containing one complete AnnData Zarr store (v2 or v3)
  - or any other format of {doc}`07_file_format_expectations` (R objects,
    Loom, 10x, MatrixMarket, MuData)

**Privacy model**
- File picker modes load data from your computer into your browser.
- Your files are **not uploaded** to Cellucid servers.

## Fast Path (Wet Lab Friendly)

1) Open [Cellucid](https://www.cellucid.com/).
2) Expand **Session** in the left sidebar and find **File or folder** under
   **Open data**.
3) Click **Open…** and choose what you have:
   - **File…** — one `.h5ad` (of any size), one `.zarr.zip` / `.zip` archive,
     one Seurat or SingleCellExperiment `.rds` / `.rda` / `.RData`
     ({doc}`17_open_r_objects`), or one `.loom`, 10x `.h5` or `.h5mu` file
     ({doc}`18_open_loom_and_10x`, {doc}`19_multimodal_mudata`)
   - **Folder…** — one `.zarr` store folder, one Cell Ranger `outs/` folder or
     one MatrixMarket folder (choose the folder itself)
   A browser's file dialog picks files *or* a folder, never both, which is why
   **Open…** asks first. What you pick decides how it is read. Or drag the file
   or folder from your file manager and drop it anywhere on the page; see
   {ref}`drop-to-open`.
4) Wait for the dataset to load.
5) Confirm success:
   - you see points rendered in the canvas
   - a field selector/legend has content

If you don’t see points after loading, jump to the troubleshooting section at the end.

```{figure} ../../../_static/screenshots/data_loading/choose-local-data-control.png
:alt: Close-up of File or folder: the full-width Open… button with its menu open, showing File… (AnnData, Seurat/SCE .rds/.rda, Loom, 10x .h5, .h5mu) and Folder… (Zarr, 10x outs, MatrixMarket), with the pointer on File….
:width: 496px

**Open…** under **File or folder**, with its menu open. **Folder…** asks
for a *folder*; **File…** asks for a single *file*.
```

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

Where **Open…** sits in the whole Session panel. Every loading path has its own
label under **Open data**, and **Save** / **Load** sit under **Session file**,
just below.
```

:::{important}
**If you cancel the file chooser, Cellucid says nothing at all.** No error, no
notification, no spinner — the app simply carries on with whatever was already
loaded. So “I chose **Folder…** and nothing happened” almost always means
the dialog was dismissed rather than confirmed. Open it again and complete the
selection.
:::

**Picking the file that is already open opens it again**, from the start: the
view returns to what the file says, and what you changed since (a selection,
a highlight, the coloring) is cleared, so save a session first to keep it
({doc}`../l_sessions_sharing/03_save_restore_ux`). The same holds for a newer copy of
the file that states the same dataset id. An `.rds` file or a MatrixMarket
folder opened again in the same tab reads what was kept the first time
({doc}`17_open_r_objects`).

## Load a Zarr Store Folder

### When to use this
- You have an AnnData Zarr store folder (`.zarr`), for example the one
  `cellucid.prepare()` or `cellucid_prepare()` writes by default, or
  `adata.write_zarr("data.zarr")`.

### What you click
Click **Open…**, then **Folder…**, and choose the `.zarr` folder itself
(not a file inside it, and not its parent folder). Your browser may ask you to confirm that the
site may read the folder's files ("upload" in its wording: nothing leaves your
computer). Dropping the folder onto the viewer does the same.

### What success looks like
- Points appear quickly: opening reads the store's metadata, the observation
  table and the embedding shown first.
- Field lists are populated.
- A gene reads its own chunks when you colour by it.

### Hidden metadata files
A Zarr v2 store keeps its metadata in files whose names start with a dot
(`.zgroup`, `.zarray`, `.zattrs`, `.zmetadata`); a Zarr v3 store uses
`zarr.json`, which is not hidden. The folder dialogs of Chrome, Edge and
Firefox include hidden files in the selection. **Safari's leaves them out**, so
a Zarr v2 folder chosen with **Folder…** in Safari arrives without its
metadata, and Cellucid says so: "This browser left out the folder's hidden
files … Drag the folder onto the page, or open its .zarr.zip." In Safari:

- **drag the `.zarr` folder** from Finder and drop it anywhere on the page, or
- open the store as a **`.zarr.zip`** archive (**Open… › File…**), or
- write it as Zarr v3 (`prepare(..., format="zarr", zarr_format=3)`), whose
  metadata file is not hidden, or serve it with `cellucid serve`, which every
  browser reads.

**Vector fields (optional)**
- `obsm` arrays named `<field>_<basis>` (scVelo's convention: `velocity_umap`,
  `velocity_umap_3d`, `velocity_branches`) are drawn by the velocity overlay
  over the embedding `X_<basis>` they are named for, when both have the same
  number of columns; nothing has to be switched on.
- See {doc}`../i_vector_field_velocity/index`.

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: The whole Cellucid window with the Pancreas sample coloured by cell_type: the sidebar shows Coloring and filtering with the categorical field selected and its eight-row legend, and the 3D embedding fills the canvas in the matching colours.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

## Load a `.h5ad` File Directly

### When to use this
- You have a `.h5ad`, of any size, and want to look at it without Python.
- You do not want to run Python locally.

(h5ad-direct-reading)=
### How a direct `.h5ad` is read

The file is read where it is, as each view needs it: nothing is uploaded and
the file is never copied into the browser's memory, so there is **no size
limit**. A background worker opens it, so the page stays responsive while it
does.

- Opening reads the metadata, the observation table and the embedding shown
  first, not the expression matrix.
- A gene is read when you colour by it. When `X` is stored by gene (CSC) or
  dense, a gene reads only its own column.
- When `X` is stored by cell (CSR, what Scanpy writes by default), a gene's
  values are spread through the whole matrix, so Cellucid reorganises it by
  gene in the background right after opening and shows its progress. A gene
  you choose meanwhile is answered by the same pass. The reorganised matrix is
  kept in memory when it fits (a share of the memory the browser reports) and
  otherwise in the browser's private file storage on your disk, removed when
  the file is closed.

Practical consequences:
- Memory holds what is shown: positions and the fields and genes in use, four
  bytes per cell each.
- A very large CSR matrix takes a while to reorganise the first time; the
  progress notice says how far it is. A column-major (CSC) matrix needs no
  reorganisation.
- A browser without private file storage (some private windows) scans the file
  for each gene once a CSR matrix is too large for memory, which is slow; use
  a window with site storage there.

### Server mode, for a file on another machine
Server mode hands the same file to this reader, from wherever it is stored:

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```

Then open the exact **Viewer URL** printed in the terminal.

### `.h5ad` minimum requirements

Your `.h5ad` needs at least one `obsm` array Cellucid can draw: a dense numeric
array with one row per cell and at least one column (a one-dimensional array
counts as one column). Which array each dimension shows is decided like this:

- A UMAP the file carries is used first. `obsm['X_umap_1d']`, `obsm['X_umap_2d']`
  and `obsm['X_umap_3d']` fill their own dimensions; when none of those three
  exists, the plain `X_umap` that `sc.tl.umap()` writes is read at the
  dimension its column count states (1, 2 or 3 columns). So an `.h5ad` straight out of Scanpy loads as
  written, with nothing renamed.
- A file with no UMAP is opened on a **best guess** (a t-SNE, a force layout,
  spatial coordinates, …, and last the first two columns of a wide array such
  as `X_pca`). The guess is marked *guessed* and a notification says what is
  shown.
- Any other array can be chosen afterwards under **Session › Data mapping**,
  at any columns.

`cellucid serve` hands a file to this same reader, so one file opens the same
way in this picker and through `cellucid serve`. See {doc}`13_data_mapping`.

Optional (but highly recommended):
- `obs` columns for coloring and filtering (numeric columns are read as
  32-bit floats: an integer above
  16,777,216 becomes its nearest 32-bit value, so keep identifier columns as
  strings if their exact values matter)
- an expression matrix (dense or sparse) for gene coloring: `X`, any layer or
  `raw.X`, named by the var index or any text var column, chosen under
  **Session › Data mapping › Genes** (see [Genes and expression](13_data_mapping.md#genes-and-expression))
- `obsp['connectivities']` if you want connectivity visualization

Optional (if you want the vector field / velocity overlay):
- Per-cell vectors in `obsm` using scVelo's naming convention `<field>_<basis>` for the embedding `X_<basis>` (e.g. `velocity_umap` for `X_umap`, `velocity_umap_3d` for `X_umap_3d`, `velocity_branches` for `X_branches`), with as many columns as that embedding.
- The overlay dropdown shows the fields bound to the embedding the view shows: `<field>_<basis>` is bound to `X_<basis>`, so it is drawn only while the view shows that embedding. The dropdown lists it by its `<field>` part (`velocity`).
- This picker binds them on its own — there is nothing to switch on, and the
  same holds for a file served by `cellucid serve` or `show_anndata()`. Any other array that fits
  an embedding (a 2-column `drift` on `X_pca`) is bound under **Session › Data
  mapping › Vector fields**; values are read only when the overlay shows them.

## Load a Zarr ZIP Archive Directly

### When to use this
- You have a complete AnnData Zarr store (v2 or v3) packaged as `.zarr.zip` or `.zip`.
- You want a portable single-file selection in any supported browser.

### Practical expectations
- Cellucid validates the ZIP directory and Zarr metadata before adopting the
  dataset.
- Gene expression is read by Zarr chunk on demand. Decoded chunks of a dense
  `X` are kept in a cache sized as a share of the memory the browser reports,
  so genes that share a chunk column strip decode those chunks once instead
  of once per gene.
- Decoded chunks and the active embedding consume browser memory; what is
  kept for reuse is released first when memory runs short.

A store folder needs no zipping: open it with **Open… › Folder…**, or serve it
with `cellucid serve`, which hands it to the browser key by key.

The archive may contain the Zarr store at its root or inside exactly one root
directory. It must preserve dotfiles such as `.zgroup`, `.zattrs`, and
`.zarray`. Vector fields use the same exact `obsm` convention as H5AD:
`velocity_umap`, `T_fwd_umap_3d`, and other
`<field>_<basis>` keys beside their `X_<basis>` embedding.

## What each refusal looks like on screen

Cellucid names the cause rather than saying “invalid file”, so the message is
usually the diagnosis. These are the three you are most likely to meet.

### You picked something that is not an H5AD

```{figure} ../../../_static/screenshots/data_loading/fail-h5ad-wrong-file.png
:alt: A Cellucid notification reading "The selected file is not a valid HDF5 file, so it cannot be an AnnData (.h5ad) file. Save it again with anndata…", with a dismiss button at its right.
:width: 760px

A CSV renamed to `.h5ad`. Cellucid checks the file's leading bytes, not its
name, so renaming a file never makes it loadable.
```

### The file has no UMAP, or nothing Cellucid can draw

A file without a UMAP is not refused: it opens on a best guess, and a
notification says what is shown, for example:

```text
This file has no default UMAP (X_umap_1d, X_umap_2d, X_umap_3d, X_umap, X_UMAP, X_wnn.umap or X_ref.umap).
Showing 'X_tsne' as 2D. Change it under Session › Data mapping.
```

Its **Show** button opens the Data mapping panel. See
{doc}`13_data_mapping`.

Only a file with no drawable `obsm` array at all is refused. The message lists
what `obsm` holds and why each entry was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell
(10,000) and at least one numeric column. obsm holds: X_df (table),
spatial (10,000 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example
with sc.tl.umap in Scanpy), then open the data again.
```

The list is the useful half: a table, a sparse array, or a dense array with the
cell that holds a missing or infinite value. Fix that array, or add an
embedding (for example `sc.tl.umap(adata)`), and open the file again.

## Common Failure Modes (and Why They Happen)

### “It worked for a demo dataset but not my data”
- Your file is missing required embeddings.
- Your Zarr ZIP is incomplete, contains multiple roots, or omits required
  `.zgroup`, `.zattrs`, or `.zarray` entries.

### “The file picker won’t let me select a folder”
- **Open… › Folder…** accepts folders only; **Open… › File…** accepts one
  file.
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
- Vector fields exist but are not named `<field>_<basis>` after an `X_<basis>` embedding in the file (for example `velocity_umap_2d` beside a bare `X_umap` with no `X_umap_2d`; scVelo's `velocity_umap` binds).
- You’re in 3D but only 2D vectors exist (or vice versa).

## Browser and file-selection boundaries

- **Folder…** uses the folder-input capability available in supported
  current desktop browsers; Safari's leaves out hidden files (see
  [Hidden metadata files](#hidden-metadata-files)).
- **File…** is an ordinary single-file input, with the same validation
  contract in Chrome, Edge, Firefox, and Safari.
- **Permission prompts**: if you deny folder access, Cellucid cannot read files.
- **Large `.h5ad`**: read from disk as each view needs it; a cell-major matrix is reorganised by gene in the background, with a progress notice.
- **Zarr ZIP root**: the archive must contain either root metadata at archive
  root or exactly one wrapper directory around the Zarr store.

## Troubleshooting (Massive)

Use this like a checklist. Most issues are diagnosable in < 2 minutes.

---

### Symptom: “I chose Folder…, but nothing happens”

**Likely causes (ordered)**
1) **You dismissed the file chooser.** A cancelled selection is deliberately
   silent — Cellucid shows no error and keeps the dataset you already had.
2) You denied folder access or a managed-browser policy blocked it.
3) The UI is open in an embedded context that disallows folder access.
4) You chose **File…** instead of **Folder…**.

**How to confirm**
- Choose **Open… › Folder…** again and complete the selection. If a
  dataset then loads, cause 1 was it.
- Open Cellucid as a standalone page and choose **Folder…**.

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
- For a cell-major (CSR) matrix, wait for the reorganisation progress notice
  to finish, or save the matrix column-major (CSC).
- Or write a Cellucid-ready AnnData once with `prepare()` and open that.

---

### Symptom: “It says Cellucid found nothing to plot”

**Likely causes**
- `obsm` holds no dense numeric array with one row per cell: only tables
  (DataFrames), sparse arrays, or arrays with another row count.
- Every array that could be drawn holds a value that cannot be: a NaN, an
  infinity, a value too large for 32-bit floats, or the same coordinates for
  every cell.

**How to confirm**
- Read the message itself: its `obsm holds: …` part lists every
  entry with its shape and, for a dense array, the reason it was not used
  (for example `NaN at cell 17`).
- In Python:

  ```python
  print({key: getattr(value, "shape", None) for key, value in adata.obsm.items()})
  ```

**Fix**
- Compute an embedding. `sc.tl.umap(adata)` writes the plain `X_umap`, which
  this picker reads as written.
- Or fix the array the message names (replace its missing values) and open the
  file again.

---

### Symptom: “It shows the wrong embedding”

**Likely cause**
- The file has no UMAP, so Cellucid opened it on a best guess (marked
  *guessed* in **Session › Data mapping**).
- Or the file has a UMAP and you wanted another array.

**Fix**
- Turn the array on, and for a wide array choose its columns, in **Session ›
  Data mapping**, mark it **Default** if every new view should show it, and
  press **Apply**. A view shows it once you choose it in **View N shows**
  (or it is the default). See {doc}`13_data_mapping`.

---

### Symptom: “Zarr ZIP selected, but it errors immediately”

**Likely causes**
- The archive does not contain exactly one valid AnnData Zarr store.
- Required `.zgroup`, `.zattrs`, or `.zarray` entries are missing.
- The ZIP is encrypted, multi-disk, uses an unsupported compression method, or
  exceeds a reader limit.

**Fix**
- With the supported AnnData version, write the store again with
  `adata.write_zarr("data.zarr")`; the resulting store must contain the Zarr
  v2 root marker `.zgroup` (stores from anndata before 0.8 have no root
  `.zattrs` and open as they are).
- Package the complete store while preserving its dotfiles, then select the
  resulting `.zarr.zip` or `.zip` file with **Open… › File…**.

---

### Symptom: “Fields list is empty, but I know I have metadata”

**Likely causes**
- `obs` was not saved in your file, or uses unexpected encodings.

**How to confirm**
- In Python: `print(adata.obs.head())`

**Fix**
- Save a clean copy of your AnnData (or write one with `prepare()`).

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset truly has no vector fields.
2) Vector fields exist, but their keys do not follow the `<field>_<basis>`
   default rule (no `X_<basis>` in the file matches, as for `velocity_umap_2d`
   beside a bare `X_umap`), or they do not have as many columns as their
   embedding, so nothing is bound to them by default.
3) The view shows an embedding the field is not bound to (for example the 3D
   UMAP while the field is bound to the 2D one).

**How to confirm**
- In Python, list `obsm` keys and look for entries like `velocity_umap` beside `X_umap`:

  ```python
  print(adata.obsm.keys())
  ```

**Fix**
- Bind the array yourself under **Session › Data mapping › Vector fields**, or
  rename it `<field>_<basis>` after its embedding `X_<basis>` (`velocity_umap`
  for `X_umap`).
- Show an embedding the field is bound to; the overlay info line names them.
- Or write the file again with `prepare(..., vector_fields={...})`.

For overlay UI behavior and deeper debugging, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

## Next Steps

Your dataset is open. The next thing most people want is to see something in it:

- **Colour the cells by a gene or a metadata column** →
  {doc}`../d_fields_coloring_legends/index`, or the three-minute version in
  {doc}`../a_orientation/03_quick_tour_60_seconds`
- **Hide the cells you are not looking at** → {doc}`../e_filtering/index`
- **Select a population and read out what is in it** →
  {doc}`../f_highlighting_selection/index` and {doc}`../h_analysis/index`

If loading itself is what you still need:

- Format expectations (`.h5ad`, Zarr, catalogs, `uns['cellucid']`): {doc}`07_file_format_expectations`
- Full troubleshooting matrix: {doc}`08_troubleshooting_data_loading`
- If browser-only loading is too slow or unreliable for your dataset → {doc}`04_server_tutorial`
- If you are working in notebooks and want programmatic interaction → {doc}`05_jupyter_tutorial`
- Writing a Cellucid-ready file to open here → {doc}`../../python_package/c_data_preparation_api/01_prepare_overview` (Python) or {doc}`../../r_package/c_data_preparation_api/01_cellucid_prepare_overview` (R)
