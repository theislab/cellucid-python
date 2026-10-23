# Share Without Running a Server (Publish a Dataset Catalog)

This tutorial shows how to publish AnnData files so anyone can open them in
Cellucid from one link, **without** running a server.

You’ll learn two closely related workflows:

1) **Public catalog (recommended): a GitHub repository or any static host**
   - You publish a folder holding your AnnData files and a small `datasets.json`.
   - Anyone opens Cellucid with `?data=<that folder>` and picks a dataset.

2) **Your own sample list (optional): run or deploy the web app yourself**
   - You point the web app's **Sample dataset** list at your catalog.
   - Useful for offline demos, an intranet, or your own deployment.

To publish a single file, you need no catalog at all: a link to the `.h5ad`,
`.zarr.zip` or Zarr store folder opens it ({doc}`14_open_from_link`). If you
just want to open your own local files right now, use
{doc}`03_browser_file_picker_tutorial`. For a complete, inspectable public
repository with three tiny datasets, use {doc}`11_custom_dataset_repository`.

## At A Glance

**Audience**
- Wet lab / beginner: follow the checklist; you do not need to understand file formats.
- Computational users: pay attention to dataset ids, file layout and host limits.
- Developers: this is the exact `datasets.json` contract the web app reads.

**Time**
- One dataset in a public repository: ~15–30 minutes (mostly uploading)

**Prerequisites**
- Python with `cellucid` installed (or R with `cellucid`), to write
  Cellucid-ready AnnData files; any AnnData file works, too
- For a public catalog: a GitHub account or a static host that allows
  cross-origin reads (CORS)

## Important: Privacy and Sharing

Publishing files on a public host publishes your dataset.

Before you do this:
- Verify there is **no private metadata** in `obs` (e.g. patient ids, internal notes).
- Remove or rename columns that should not be public; `prepare(obs_keys=...)`
  writes only the columns you name.

If you need **private** sharing:
- Use **Server Mode** ({doc}`04_server_tutorial`) on a private machine + VPN/SSH tunnel.
- Or send the file and let each collaborator open it with **Files**.

## The Catalog Folder Layout

A **dataset catalog** is a folder holding `datasets.json` and the AnnData files
it lists:

```text
catalog/                         # the folder you link to
├── datasets.json                # required
├── pbmc_demo.zarr/              # a Zarr store folder (recommended on static hosts)
├── lung_atlas.h5ad              # an .h5ad file
├── organoid.zarr.zip            # a zipped Zarr store
└── pbmc_demo.cellucid-session   # optional: a dataset's default view
```

**Key rules**
- `datasets.json` sits in the folder you link to.
- Every path in it is **relative** to that folder, and stays below it.
- Files can sit in subfolders (`data/pbmc_demo.zarr/`).

## Step 1 — Write Cellucid-ready AnnData Files

Any AnnData file opens, but a file written by `cellucid.prepare()` (Python) or
`cellucid_prepare()` (R) reads fastest over a network and carries its own name,
description, citation and launch defaults in `uns['cellucid']`.

```python
import cellucid

cellucid.prepare(
    "catalog/pbmc_demo.zarr",
    format="zarr",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    dataset_id="pbmc_demo",
    dataset_name="PBMC demo",
    dataset_description="PBMC scRNA-seq",
    latent_key="X_pca",
)
```

The arguments you should understand, even as a beginner:

- `dataset_id`: a **stable** identifier; sessions and community annotations
  are keyed by it. 1–180 ASCII letters, digits, `.`, `_` or `-`, starting
  with a letter or digit and not ending with `.`. Avoid changing it after you
  publish ({doc}`06_dataset_identity_why_it_matters`).
- `dataset_name` / `dataset_description`: what the dataset list and the
  dataset panel show.
- `format` (required; the path ends in its suffix): `"zarr"` is a folder of
  many small files, for GitHub (no file over 100 MB, no LFS); `"h5ad"` is one
  file, for single-file hosts (Zenodo, Figshare, Hugging Face, S3, a lab
  server), local files, `cellucid serve` and Jupyter; `"zarr.zip"` is the
  Zarr store in one file, wherever an `.h5ad` goes
  ({ref}`which-format`).
- `vector_fields` (optional): per-cell vectors drawn by the velocity overlay
  ({doc}`../i_vector_field_velocity/index`).

Every argument: {doc}`../../python_package/c_data_preparation_api/index`.

## Step 2 — Write `datasets.json`

`datasets.json` names the datasets, where each file is, and which one opens
first:

```json
{
  "schema": 2,
  "default": "pbmc_demo",
  "datasets": [
    { "id": "pbmc_demo", "path": "pbmc_demo.zarr/" },
    { "id": "lung_atlas", "path": "lung_atlas.h5ad" }
  ]
}
```

- `schema` is the integer `2`.
- `default` is the id of the dataset the catalog opens on.
- Each `id` is the dataset's stable id. When a file states its own id in
  `uns['cellucid']['id']` (every file `prepare()` writes does), the catalog
  must list it under that same id, or the file is refused when it is opened.
- Each `path` ends in `/` for a Zarr store folder, `.h5ad` for an HDF5 file,
  or `.zarr.zip` for a zipped store. Write each part in its canonical URL
  encoding (`cell%20atlas.h5ad`), without `.`, `..`, a query or a fragment.

The dataset list shows each dataset's name and cell count, read from its
file: names, descriptions and citations are never written in the catalog.
The full contract, including the optional default session, is in
{doc}`07_file_format_expectations`.

## Step 3 — Validate Locally (Before You Upload)

- **One file:** open it with the **Zarr folder**, **H5AD** or **Zarr ZIP**
  button under **Files**. If it opens there, it will open from a link.
- **The whole catalog:** serve the folder with any static server that sends
  CORS headers, then open `?data=http://127.0.0.1:8000/` in Cellucid:

  ```bash
  cd catalog
  python -m http.server 8000   # add CORS headers if Cellucid runs on another origin
  ```

This avoids debugging “hosting problems” that are actually file problems.

```{figure} ../../../_static/screenshots/web_app/app-overview-cell-type.png
:alt: Cellucid web app with the sidebar open and a single-cell embedding colored by cell type.
:width: 1440px

A loaded dataset in Cellucid: the sidebar controls the active view while the categorical legend maps directly to the colored points.
```

## Step 4 — Publish to a Public GitHub Repository

### 4.1 Create a repository

A data-only repository is simplest, e.g. `my-lab/cellucid-data`, with the
catalog folder at its root or in a folder such as `catalog/`.

### 4.2 GitHub constraints (do not skip)

- **Public repository only**: Cellucid reads files through
  `raw.githubusercontent.com`.
- **Branches**: shorthand without a branch reads the repository's default
  branch (whatever it is called); Cellucid never probes branch names. Name
  another branch explicitly:
  `owner/repo@my-branch/catalog`, or
  `https://github.com/owner/repo/tree/my-branch/catalog`.
- **File size**: GitHub refuses files over 100 MiB. A Zarr store folder is
  many small chunk files, so it fits where a large `.h5ad` does not; see
  {doc}`15_host_on_github` for the steps, the limits and the Git LFS
  trade-off.

### 4.3 What to enter in Cellucid

Expand **Session**, type the path into the **Open from link** field, and press
**Connect**:

- `owner/repo` (if `datasets.json` is at the repository root)
- `owner/repo/catalog` (if your catalog is in a folder)
- `owner/repo@my-branch/catalog` (on another branch)

```{figure} ../../../_static/screenshots/data_loading/connect-github-catalog.png
:alt: Close-up of the Open from link control with an owner/repo path typed into the text field and a mouse pointer resting on the Connect button beside it.
:width: 496px

The **Open from link** field takes the catalog folder (or one file in it).
Press **Connect** beside it.
```

Cellucid then reads `datasets.json`, describes every listed file from its
metadata (name, cell and gene counts), shows the dataset list, and opens the
catalog's default dataset; each other dataset is opened when you choose it.

**Shareable URL**

```text
https://www.cellucid.com/?data=gh:owner/repo/catalog&dataset=pbmc_demo
```

## Optional: Your Own Sample Dataset List

The **Sample dataset** list reads the catalog named by the page's
`<meta name="cellucid-sample-catalog" content="...">` element in
`cellucid/index.html` (www.cellucid.com names
`https://theislab.github.io/cellucid-datasets/anndata/`). One `sampleCatalog`
query field replaces it for a single launch:

```text
http://localhost:8000/?sampleCatalog=http://localhost:9000/catalog/
```

The value is an absolute `http(s)` folder address ending in `/`. This is a
development and deployment hook; to share a catalog with others, use a data
link, `?data=https://host/catalog/` ({doc}`14_open_from_link`).

```{note}
A catalog on another origin must allow the Cellucid origin with CORS headers,
including range requests (`Range` request header; `Content-Range` and
`Content-Length` exposed), so the browser can read the files in place.
```

## Edge Cases

- **Changed dataset ids after publishing**: links break; saved sessions and
  community annotations no longer match.
- **A file whose `uns['cellucid']['id']` differs from its catalog id**: the
  file is refused when chosen; list it under the id it states.
- **Wrong repository path or branch**: Cellucid reports that the folder holds
  no `datasets.json`, naming the address it read.
- **A Zarr store without consolidated metadata** on a static host: write it
  with `prepare()` or call `zarr.consolidate_metadata()` before publishing.
- **GitHub Pages and hidden files**: Jekyll skips dot-files such as `.zmetadata`
  and `.zgroup`; publish a `.nojekyll` file at the site root.

## Troubleshooting (Large, Practical)

### Symptom: “Open from link says the folder holds no `datasets.json`”

**Likely causes (ordered)**
1) The path points at the wrong folder.
2) The catalog is on another branch than you think.
3) You published the files but no `datasets.json`.

**How to confirm**
- Read the notification. It names the folder Cellucid read:

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/ holds no datasets.json, so it is not a dataset catalog. …
  ```

- Open that folder's `datasets.json` in a browser tab. If it 404s there too,
  the path or branch is wrong; if it loads there but not in Cellucid, suspect
  the network or CORS.

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: Two stacked Cellucid notifications, the upper naming the raw.githubusercontent.com URL that returned not found and the lower explaining that the GitHub repository could not be reached and suggesting a check for a typo.
:width: 776px

What a wrong repository path looks like. The upper notification hands you the
URL to test by hand.
```

---

### Symptom: “Invalid datasets.json: the catalog is missing 'schema' …”

The catalog has no `"schema"` key. A catalog holds exactly `"schema": 2`,
`"default"` and `"datasets"`, and lists AnnData files (each written as a
Cellucid-ready AnnData with `prepare()`).

---

### Symptom: “The dataset list shows up, but choosing a dataset fails”

**Likely causes**
- A `path` names a file or folder that does not exist (404).
- The file states another id than its catalog entry.
- The file is not an AnnData file, or a Zarr store lacks its metadata.

**Fix**
- Compare each `path` with the published names; the notification names the
  dataset and the reason.

---

### Symptom: “It loads, but gene search is slow”

**Likely causes**
- The expression matrix is stored cell-major (CSR): the first gene waits for
  its index to be built.

**Fix**
- Write the file with `prepare()`, which stores the matrix gene-major (CSC).

---

### Symptom: “The vector field overlay toggle is disabled / no fields appear”

**Likely causes**
1) The file has no `obsm['<field>_<basis>']` array beside its embedding `obsm['X_<basis>']` (`velocity_umap` beside `X_umap`), with the same number of columns.
2) The field is bound to another embedding than the one the view shows.

**Fix**
- Write the vectors with `prepare(..., vector_fields=...)`, or name them
  `<field>_<basis>` in `obsm` after their embedding `X_<basis>`.
- Show an embedding the field is bound to; the overlay info line names them.

See {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`.

---

### Symptom: “GitHub blocks my push / files are too big”

**Fix**
- Write the dataset as a `.zarr` folder: each chunk is a small file.
- Or host the files elsewhere (Zenodo, Hugging Face, an object store) and link
  to them, or use server mode.

## Next Steps

- Want to open a file from your computer (no publishing)? → {doc}`03_browser_file_picker_tutorial`
- Want a complete repository to copy and adapt? → {doc}`11_custom_dataset_repository`
- Data on another machine? → {doc}`04_server_tutorial`
- File and catalog expectations: {doc}`07_file_format_expectations`
- If loading fails: {doc}`08_troubleshooting_data_loading`
