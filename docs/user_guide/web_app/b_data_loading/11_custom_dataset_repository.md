# Publish a Custom Dataset Repository

A public Cellucid dataset repository is useful when you want a colleague,
reviewer, class, or workshop participant to open your datasets without
installing Python or waiting for you to keep a server running. The repository
is ordinary static content: a **dataset catalog**, one `datasets.json` beside
the files it lists.

The small
[cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets)
repository is the reference implementation. It contains a catalog of five
deterministic synthetic datasets: three AnnData Zarr folders written by
`cellucid.prepare(..., format="zarr")`, a CITE-seq MuData Zarr folder written by
`cellucid.prepare(data=mdata, format="zarr")`, and the same CITE-seq cells as a
Seurat object saved as `saveRDS()` saves it; the 2-D AnnData again as an
`.h5ad` and a `.zarr.zip`; the `.nojekyll` file GitHub Pages needs, a wet-lab-friendly README, one root-level
`generate_datasets.py` file that reproducibly renews them, and the policy files
any public repository needs (`LICENSE`, `CONTRIBUTING.md`, `SECURITY.md`,
`SUPPORT.md`, `CODE_OF_CONDUCT.md`). Deliberately, that is all: there are no
tests, validation harness, package scaffolding, workflows, `.github/`
directory, or site assets to copy. The README walks from one ordinary Python
`prepare(...)` call through the GitHub upload and optional GitHub Pages setup,
and names the exact catalog path and dataset links to try.

:::{important}
Publishing data in a public repository publishes the coordinates, metadata,
gene values, graphs, vector fields, and provenance that you include. Inspect
the files for private identifiers before you push them. Use
{doc}`04_server_tutorial` for data that must remain private.
:::

## Try the reference repository first

This takes about a minute and proves that the browser can reach GitHub raw
content from your network.

1. Open [Cellucid](https://www.cellucid.com/).
2. Expand **Session** in the left sidebar.
3. In the **Link or server** bar under **Open data**, paste the repository's
   catalog path, `theislab/cellucid-demo-custom-datasets/catalog`.
4. Click **Open** (or press Enter).
5. Wait for the connection message, then choose one of the five entries in
   the dataset dropdown above the bar, whose label now reads
   "Dataset from <host>".

```{figure} ../../../_static/screenshots/data_loading/connect-github-catalog.png
:alt: The Link or server bar with the GitHub path "theislab/cellucid-demo-custom-datasets/catalog" typed into it and the pointer on its Open button.
:width: 496px

The **Link or server** bar with the reference path entered. Type the catalog
folder (or one file in it) and press **Open**.
```

The shorthand reads the repository's default branch (`main`, `master` or any
other), as GitHub resolves `HEAD`, without a request to find it; Cellucid never
probes branch names.

On success the dataset dropdown lists the catalog's datasets in place of the
samples, labelled "Dataset from <host>", and the catalog's default one opens.
One line appears under the bar reading **Link**, the address and
**Disconnect** — that line is how you tell a live connection from a path you
merely typed. **Disconnect** brings the samples back.

### Five tiny datasets, five useful checks

| Exact dataset ID | Contents | A useful first question |
|---|---|---|
| `synthetic-cell-types-2d` | 72 cells, 8 genes, four `obs` fields, two 2D embeddings, and a weighted neighbour graph | Do the three `cell_type` categories agree with their marker-gene patterns? |
| `synthetic-development-3d` | 96 cells, 10 genes, five `obs` fields, a 3D embedding `X_branches`, a weighted neighbour graph, and the velocity field `velocity_branches` drawn over it | Does the vector field follow both developmental branches when you orbit the view? |
| `synthetic-trajectory-1d` | 48 cells, 6 genes, four `obs` fields, a 1D embedding `X_trajectory`, a weighted neighbour graph, and the velocity field `velocity_trajectory` drawn over it | How do `pseudotime` and gene activity change along the line? |
| `synthetic-citeseq` | 60 cells measured twice: a MuData of the modalities `rna` (8 genes) and `prot` (4 proteins), `cell_type` and `sample`, `X_umap` and `X_pca`, a weighted neighbour graph | Under **Data mapping › Genes › Modality**, does `CD19-ADT` in `prot` mark the cells `SYN-B-MARKER` marks in `rna`? ({doc}`19_multimodal_mudata`) |
| `synthetic-citeseq-seurat` | the same 60 cells as a Seurat object (`.rds`): assays `RNA` and `ADT`, reductions `pca` and `umap`, the graph `RNA_snn` | Does the object open as it was saved, with its assays as modalities? ({doc}`17_open_r_objects`) |

Each dataset uses the navigation mode expected for its geometry: Planar in 1D
and 2D, Orbit in 3D. A catalog lists any file Cellucid opens, so the Seurat
`.rds` is listed like the Zarr folders; it is read once when it is listed or
opened, where the Zarr folders are read in place.

A dataset-specific link removes the final selection step:

```text
https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog&dataset=synthetic-development-3d
```

Without `dataset=`, a link opens on the dataset the catalog names as its
`default`.

## The repository contract

The repository can contain source code and explanatory files, but the path you
give Cellucid must resolve to a **catalog folder** with this shape:

```text
my-cellucid-datasets/
├── README.md
├── generate_datasets.py
└── catalog/
    ├── datasets.json
    ├── study-a.zarr/
    ├── study-b.zarr/
    └── study-a.cellucid-session    # optional default view of study-a
```

The generator is the reference repository's sole reproducibility script. Your
own data repository may use a notebook or another preparation workflow
instead; Cellucid itself only reads the catalog folder.

`datasets.json` is the exact catalog. Its top-level object contains only
`schema`, `default`, and `datasets`:

```json
{
  "schema": 2,
  "default": "study-a",
  "datasets": [
    {
      "id": "study-a",
      "path": "study-a.zarr/",
      "session": {
        "path": "study-a.cellucid-session",
        "bytes": 3256,
        "sha256": "9465cfc1cd2dca59c0e9cc9c3ff3f4a26aac8dc57004eb51a641a2e437693c80"
      }
    },
    { "id": "study-b", "path": "study-b.zarr/" }
  ]
}
```

The rules are deliberately strict:

- `schema` is exactly `2`.
- `default` is an exact dataset id present in `datasets`.
- `datasets` is a non-empty array of entries holding `id`, `path` and,
  optionally, `session`; nothing else.
- Every `id` is unique and follows the portable id rule (1–180 ASCII letters,
  digits, `.`, `_` or `-`, starting with a letter or digit, not ending with
  `.`). A file that states its own id (`uns['cellucid']['id']`) is listed under
  exactly that id.
- Every `path` is unique, relative to the catalog folder, cannot escape it, and
  names a Zarr store folder (ending in `/`), an `.h5ad`, or a `.zarr.zip` — or
  a file in any other format Cellucid reads, by its own ending (`.h5mu`,
  `.loom`, `.rds`, …).
- `session`, when present, names a `.cellucid-session` file with its exact
  `bytes` and lowercase `sha256`.

Names, descriptions, citations and counts are not in the catalog: the dataset
list reads them from each file's metadata and `uns['cellucid']` when the
catalog is opened, and a dataset's file is opened when it is chosen. The full
contract, including `uns['cellucid']`, is in
{doc}`07_file_format_expectations`.

## Write the AnnData files

### From Python

Use `cellucid.prepare(...)` once per dataset. Keep every array aligned to the
same cell order and give each dataset a stable portable id.

```python
import cellucid

cellucid.prepare(
    "catalog/study-a.zarr",
    format="zarr",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    dataset_id="study-a",
    dataset_name="Study A",
    dataset_description="A public, documented study.",
    citation_text="Author et al. (2026)",
    citation_url="https://doi.org/10.0000/example",
)
```

The complete guide explains every argument:
{doc}`../../python_package/c_data_preparation_api/index`.

Then write `datasets.json` by hand (it is a few lines), listing each file under
the id it states.

### From R

`cellucid-r`'s `cellucid_prepare()` writes the same Cellucid-ready AnnData.
Begin with {doc}`../../r_package/index`.

### From any AnnData file

Any `.h5ad` or AnnData Zarr store can be listed as it is. A file that states no
id is listed under the catalog's id; to carry its name, description and
citation, add `uns['cellucid']` (or write it again with `prepare()`, which also
stores the expression matrix gene-major for fast gene reads).

## Validate before publishing

Test the exact files, not just the script that was meant to create them.

### 1. Test one file in the browser

Open Cellucid, expand **Session**, click **Open…** under **File or folder**,
then **Folder…** (or **File…** for an `.h5ad` or `.zarr.zip`), and choose
`catalog/study-a.zarr`. Confirm:

- the reported cell and gene counts;
- every intended dimension and its default navigation mode;
- representative categorical and continuous `obs` fields;
- at least one gene;
- the neighbour graph, if the file has one;
- each vector field in each bound dimension.

Repeat for every dataset. Opening a file does not read `datasets.json`, so this
step isolates problems inside one file.

### 2. Test the complete catalog locally

Serve the catalog folder with any static server that sends CORS headers, and
open `?data=http://127.0.0.1:8000/` in Cellucid. Select each catalog entry in
turn and repeat the checks above.

### 3. Inspect the public raw catalog

After publishing, open this pattern in a browser:

```text
https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
```

It must return the JSON file itself, not an HTML sign-in page, a Git LFS
pointer, or a 404 response.

### 4. Exercise the real share link

Use a private/incognito window that is not relying on your local files. Open:

```text
https://www.cellucid.com/?data=gh:<owner>/<repo>/<path>&dataset=<exact-dataset-id>
```

Check the browser console and Network panel as well as the canvas. A convincing
test reaches a stable render, loads one metadata field and one gene, and has no
failed requests or uncaught console errors.

## Publish and address the right branch

1. Create a **public** GitHub repository.
2. Commit the catalog folder: `datasets.json` and every file it lists.
3. Push to the branch you intend to serve.
4. Open the raw catalog URL and then the Cellucid share URL.

For the repository's default branch, use the compact form:

```text
owner/repo/path/to/catalog
```

For a different branch, tag or commit, put `@ref` immediately after the
repository name:

```text
owner/repo@release-2026/path/to/catalog
```

Do not write `owner/repo/release-2026/path/to/catalog`; that means the
`release-2026/path/to/catalog` folder on the default branch. Cellucid also
accepts an exact GitHub tree URL or exact `raw.githubusercontent.com` URL, both
of which already identify the branch.

If changing files in place would make an old result ambiguous, publish a new
dataset id or a new branch and share a new exact URL. Saved sessions also check
dataset identity, so stable ids should describe stable cell ordering and
meaning.

GitHub Pages (with its required `.nojekyll` file), what each host costs per
request, the size limits, updating a published store, and Git LFS are in
{doc}`15_host_on_github`. A Pages-hosted catalog is linked by its folder
address, ending in `/`:

```text
https://www.cellucid.com/?data=https%3A%2F%2F<owner>.github.io%2F<repo>%2Fcatalog%2F&dataset=<exact-dataset-id>
```

## Keep the download small and predictable

The browser reads only what each view needs, but the choices made when writing
the files still matter:

- write only the genes and metadata fields the public visualization needs
  (`prepare(gene_identifiers=..., obs_keys=...)`);
- keep the expression matrix gene-major (CSC), as `prepare()` writes it, so one
  gene is one contiguous read;
- prefer coarse categorical fields to high-cardinality ones, whose colours are not distinguishable;
- preserve the original full analysis object outside the public files.

GitHub blocks files larger than 100 MiB, so publish each dataset as a Zarr
store folder of small chunk files ({doc}`15_host_on_github` explains why, with
the limits and the Git LFS trade-off). For datasets beyond GitHub's repository
and Pages sizes, use a static object host or {doc}`04_server_tutorial`.

For a static host or CDN (a bucket, a lab web server), publish the same catalog
folder with correct CORS headers (including range requests), then link to the
encoded folder and dataset:

```text
https://www.cellucid.com/?data=https%3A%2F%2Fdata.example.org%2Fcatalog%2F&dataset=study-a
```

The GitHub shorthand is specific to public GitHub raw content; it is not a
generic private-repository authentication mechanism. Every other way to open
data from a link (CELLxGENE, GEO, a single `.h5ad` file on any host) is in
{doc}`14_open_from_link`.

## Fast diagnosis

The most common mistake is pointing Cellucid at the repository instead of the
catalog folder inside it. That produces one notification, *The linked data
could not be opened: https://raw.githubusercontent.com/…/ holds no
datasets.json, so it is not a dataset catalog. …*, naming the exact raw folder
it read, which you can paste straight into a browser tab.

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: A Cellucid notification saying the linked data could not be opened because the GitHub folder holds no datasets.json, after a repository path rather than its catalog folder was opened.
:width: 760px

Connecting to a repository with the catalog folder left off. Opening the URL
the notification names yourself confirms whether the problem is your path,
your branch, or your network.
```

| Symptom | Check first |
|---|---|
| … holds no datasets.json, so it is not a dataset catalog | The value must point at the catalog folder (or link one file or Zarr store folder instead); open `datasets.json` under the folder the notification names. |
| `Invalid datasets.json: …` | The message names the key or value: an unknown key, a missing `default`, a path that names neither a store folder nor a file Cellucid reads, a repeated id or path. |
| A dataset cannot be listed or opened | The message names the dataset and the reason: a missing file (404), an id that differs from the file's `uns['cellucid']['id']`, or a file in no format Cellucid reads. |
| A request returns tiny text | Check whether Git LFS stored a pointer instead of the file. |
| A Zarr store opens locally but not from Pages | Add `.nojekyll`, so Pages publishes the store's dot-files. |
| It works locally but not on the lab network | Test `raw.githubusercontent.com`; a firewall or content blocker may require another host or server mode. |
| A saved session refuses the dataset | Confirm the exact dataset id, cell count, and cell ordering are the same as when it was saved. |

Continue with {doc}`08_troubleshooting_data_loading` for browser, CORS, memory,
embedding, and server diagnostics. The reference repository's
[README](https://github.com/theislab/cellucid-demo-custom-datasets#readme)
cross-references this guide and gives the shortest Python-to-GitHub path,
including its single-file regeneration workflow, browser checks, privacy,
updates, multiple datasets, and GitHub Pages.
