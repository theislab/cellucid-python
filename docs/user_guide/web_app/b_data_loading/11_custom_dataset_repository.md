# Publish a Custom Dataset Repository

A public Cellucid dataset repository is useful when you want a colleague,
reviewer, class, or workshop participant to open a prepared dataset without
installing Python or waiting for you to keep a server running. The repository
is ordinary static content: one exact catalog plus one or more complete
prepared dataset directories.

The small
[cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets)
repository is the reference implementation. It contains three deterministic
synthetic datasets and a wet-lab-friendly README. Deliberately, that is all:
there are no repository scripts, tests, workflows, or site assets to copy.
The README walks from one ordinary Python `prepare(...)` call through the
GitHub upload and optional GitHub Pages setup.

:::{important}
Publishing prepared data in a public repository publishes the coordinates,
metadata, gene values, graphs, vector fields, and provenance that you include.
Inspect the generated files for private identifiers before you push them.
Use {doc}`04_server_tutorial` for data that must remain private.
:::

## Try the reference repository first

This takes about a minute and proves that the browser can reach GitHub raw
content from your network.

1. Open [Cellucid](https://www.cellucid.com/).
2. Expand **Session** in the left sidebar.
3. In **GitHub data:**, paste this exact value:

   ```text
   theislab/cellucid-demo-custom-datasets/exports
   ```

4. Click **Load**.
5. Wait for the connection message, then choose one of the three entries under
   **Sample datasets:**.

The shorthand above resolves deterministically to the `main` branch. Cellucid
does not probe `master`, `gh-pages`, or any other branch.

### Three tiny datasets, three useful checks

| Exact dataset ID | Contents | A useful first question |
|---|---|---|
| `synthetic-cell-types-2d` | 72 cells, 8 genes, four `obs` fields, a 2D embedding, and 144 weighted edges | Do the three `cell_type` categories agree with their marker-gene patterns? |
| `synthetic-development-3d` | 96 cells, 10 genes, five `obs` fields, a 3D embedding, 187 weighted edges, and 3D `velocity_umap` vectors | Does the vector field follow both developmental branches when you orbit the view? |
| `synthetic-trajectory-1d` | 48 cells, 6 genes, four `obs` fields, a 1D embedding, 93 weighted edges, and 1D `velocity_umap` vectors | How do `pseudotime` and gene activity change along the line? |

Each dataset uses the navigation mode expected for its geometry: Planar in 1D
and 2D, Orbit in 3D.

Dataset-specific links remove the final selection step:

- [Open the 2D cell-type example](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-cell-types-2d)
- [Open the 3D development example](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-development-3d)
- [Open the 1D trajectory example](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-trajectory-1d)

A GitHub catalog with exactly one dataset loads that sole entry automatically.
A catalog with multiple datasets requires an exact `dataset=` ID in a startup
URL. An interactive connection without `dataset=` stays connected and asks the
user to choose under **Sample datasets:**; it never guesses the first entry.

## The repository contract

The repository can contain source code and explanatory files, but the path you
give Cellucid must resolve to an **exports root** with this shape:

```text
my-cellucid-datasets/
├── README.md
└── exports/
    ├── datasets.json
    ├── study-a/
    │   ├── dataset_identity.json
    │   ├── points_2d.bin.gz
    │   └── ...
    └── study-b/
        ├── dataset_identity.json
        ├── points_3d.bin.gz
        └── ...
```

`datasets.json` is the exact catalog. Its top-level object contains only
`version`, `default`, and `datasets`:

```json
{
  "version": 1,
  "default": "study-a",
  "datasets": [
    {
      "id": "study-a",
      "path": "study-a/",
      "name": "Study A",
      "description": "A public, documented study.",
      "n_cells": 40000,
      "n_genes": 2500
    },
    {
      "id": "study-b",
      "path": "study-b/"
    }
  ]
}
```

The rules are deliberately strict:

- `version` is exactly `1`.
- `default` is a non-empty exact dataset ID present in `datasets`.
- `datasets` is a non-empty array.
- Every entry has exactly one unique `id` and one unique, safe relative
  directory `path`.
- Every `path` ends in `/`; it is relative to the exports root and cannot
  escape it.
- The only optional catalog fields are `name`, `description`, `n_cells`, and
  `n_genes`.
- If an optional field is present, it exactly matches the canonical value in
  that dataset's `dataset_identity.json`.
- Every listed directory is a complete current prepared generation. A stray or
  partial listed directory rejects the connection instead of producing a
  partly working catalog.

At connection time, Cellucid fetches and validates `datasets.json` and every
listed `dataset_identity.json`. It then fetches the selected dataset's
manifests, coordinates, metadata fields, genes, graph, and vector payloads as
the interaction requires them. This is why a catalog can remain quick to open
without downloading every gene from every dataset.

For the prepared-directory contract itself, including exact embedding and
vector-field names, see
{doc}`07_folder_file_format_expectations_high_level_link_to_spec`.

## Build the prepared directories

### From Python

Use `cellucid.prepare(...)` once per dataset. Keep every array aligned to the
same cell order, give the dataset a stable portable ID, and write each dataset
to its own directory beneath `exports/`.

```python
from cellucid import prepare

prepare(
    latent_space=adata.obsm["X_umap_2d"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    out_dir="exports/study-a",
    dataset_name="Study A",
    dataset_id="study-a",
    obs_categorical_dtype="uint16",
    X_umap_2d=adata.obsm["X_umap_2d"],
)
```

This short example shows the repository boundary; the complete preparation
guide explains categorical dtypes, continuous quantization, sparse expression,
connectivity invariants, vector fields, compression, and deterministic
provenance: {doc}`../../python_package/c_data_preparation_api/01_prepare_export_overview`.

After all dataset directories exist, generate the catalog:

```python
from cellucid.prepare_data import generate_datasets_manifest

generate_datasets_manifest(
    "./exports",
    default_dataset="study-a",
)
```

The helper validates every candidate directory before atomically replacing
`datasets.json`. `default_dataset` is required; this prevents an accidental
directory order from becoming public behavior.

### From R

`cellucid-r` writes the same prepared format. Begin with
{doc}`../../r_package/c_data_preparation_api/01_cellucid_prepare_overview`, then
follow {doc}`../../r_package/d_viewing_loading/02_host_exports_for_sharing` to
assemble and validate an exports root. You can use the Python manifest helper
above, or write the exact v1 JSON from R after checking every identity.

### From AnnData without a static export

Do not put an H5AD file or a Zarr store directly in this catalog. The GitHub
loader consumes prepared Cellucid directories. For an H5AD file, a Zarr
directory, or an in-memory `AnnData`, use the Python
{doc}`server <04_server_tutorial>` or {doc}`Jupyter <05_jupyter_tutorial>`
path. The browser file picker separately accepts a small H5AD or one complete
Zarr v2 ZIP; see {doc}`03_browser_file_picker_tutorial`.

## Validate before publishing

Test the exact artifacts, not just the script that was meant to create them.

### 1. Test one prepared folder in the browser

Open Cellucid, expand **Session**, click **Prepared** under **Local data:**, and
select `exports/study-a/`. Confirm:

- the reported cell and gene counts;
- every intended dimension and its default navigation mode;
- representative categorical and continuous `obs` fields;
- at least one gene;
- connectivity, if declared;
- each vector field in each declared dimension.

Repeat for every dataset. Local folder loading does not use `datasets.json`, so
this step isolates problems inside a prepared dataset.

### 2. Test the complete catalog through the local server

```bash
cellucid serve ./exports --no-browser
```

Open the exact **Viewer URL** printed by Cellucid. It ends in
`/?source=remote`. Select each catalog entry in turn and repeat the checks
above. The server rejects an incomplete exports root before presenting it as a
working collection.

### 3. Inspect the public raw catalog

After publishing, open this pattern in a browser:

```text
https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
```

For the reference repository, the exact URL is:

```text
https://raw.githubusercontent.com/theislab/cellucid-demo-custom-datasets/main/exports/datasets.json
```

It must return the JSON file itself, not an HTML sign-in page, a Git LFS
pointer, or a 404 response.

### 4. Exercise the real share link

Use a private/incognito window that is not relying on your local files. Open:

```text
https://www.cellucid.com/?github=<owner>/<repo>/<path>&dataset=<exact-dataset-id>
```

Check the browser console and Network panel as well as the canvas. A convincing
test reaches a stable render, loads one metadata field and one gene, and has no
failed dataset requests or uncaught console errors.

## Publish and address the right branch

1. Create a **public** GitHub repository.
2. Commit `exports/datasets.json` and every file beneath each listed dataset
   path.
3. Push the complete generation to the branch you intend to serve.
4. Open the raw catalog URL and then the Cellucid share URL.

For a repository on `main`, use the compact form:

```text
owner/repo/path/to/exports
```

For a different branch, put `@branch` immediately after the repository name:

```text
owner/repo@release-2026/path/to/exports
```

Do not write `owner/repo/release-2026/path/to/exports`; that means the
`release-2026/path/to/exports` directory on `main`. Cellucid also accepts an
exact GitHub tree URL or exact `raw.githubusercontent.com` URL, both of which
already identify the branch.

If changing files in place would make an old result ambiguous, publish a new
dataset ID or a new branch and share a new exact URL. Saved sessions also check
dataset identity, so stable IDs should describe stable cell ordering and
meaning.

### Optional: enable GitHub Pages

The GitHub loader above does not require Pages. Enable it when you also want a
stable static exports URL:

1. Open the repository's **Settings → Pages**.
2. Under **Build and deployment**, choose **Deploy from a branch**.
3. Choose `main`, `/(root)`, and **Save**.
4. Wait for the deployment, then open:

   ```text
   https://<owner>.github.io/<repo>/exports/datasets.json
   ```

The response must be the catalog JSON. A Pages-hosted Cellucid link uses that
exports root explicitly:

```text
https://www.cellucid.com/?exportsBaseUrl=https%3A%2F%2F<owner>.github.io%2F<repo>%2Fexports%2F&dataset=<exact-dataset-id>
```

## Keep the download small and predictable

Prepared datasets are already arranged for on-demand browser access, but the
choices made during preparation still matter:

- export only genes and metadata fields that the public visualization needs;
- use sparse gene expression when the source is sparse;
- use gzip and choose quantization deliberately;
- avoid high-cardinality categorical fields that produce unusable legends;
- keep `datasets.json` and every identity concise;
- preserve the original full-precision analysis object outside the public
  visualization export.

GitHub blocks ordinary Git files larger than 100 MiB and recommends keeping
repositories small; see
[About large files on GitHub](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).
Git LFS stores a pointer in ordinary repository content, so it is not a
transparent replacement for files the GitHub raw loader fetches. If one
prepared payload approaches that boundary, use {doc}`04_server_tutorial` or a
static object host designed for large public downloads.

For a static host or CDN, publish the same exports root with correct CORS and
content-encoding headers, then provide an exact encoded root and dataset:

```text
https://www.cellucid.com/?exportsBaseUrl=https%3A%2F%2Fdata.example.org%2Fexports%2F&dataset=study-a
```

The GitHub shorthand is specific to public GitHub raw content; it is not a
generic private-repository authentication mechanism.

## Fast diagnosis

| Symptom | Check first |
|---|---|
| `datasets.json not found` | The value must point at the exports root, not one dataset directory; open the raw catalog URL. |
| `Invalid datasets.json` | Remove unsupported keys, require `default`, end each `path` in `/`, and make IDs and paths unique. |
| A dataset is reported invalid | Open its `dataset_identity.json`; confirm the catalog's optional name, description, and counts match it exactly. |
| The catalog connects but no dataset opens from a link | Add `dataset=<exact-dataset-id>` for a catalog with multiple entries. |
| A binary request returns tiny text | Check whether Git LFS stored a pointer instead of the prepared bytes. |
| It works locally but not on the lab network | Test `raw.githubusercontent.com`; a firewall or content blocker may require another host or server mode. |
| A saved session refuses the dataset | Confirm the exact dataset ID, cell count, and cell ordering are the same generation. |

Continue with {doc}`08_troubleshooting_data_loading` for browser, CORS, memory,
embedding, and server diagnostics. The reference repository's
[README](https://github.com/theislab/cellucid-demo-custom-datasets#readme)
gives the shortest Python-to-GitHub path, including browser checks, privacy,
updates, multiple datasets, and GitHub Pages.
