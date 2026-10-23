# Host a Dataset on GitHub

A public GitHub repository is a free, permanent home for a dataset that anyone
can open in Cellucid from a link: no server, no sign-in, no upload service.
Publish it as a **Zarr store folder**. This page explains why, walks through
the steps, and gives the numbers measured against GitHub's own hosts.

:::{important}
A public repository publishes everything in it: coordinates, metadata, gene
values, graphs and provenance. Check the files for private identifiers before
you push. For data that must stay private, use {doc}`04_server_tutorial`.
:::

## Why a Zarr folder

GitHub refuses any file over 100 MiB and warns above 50 MiB
([About large files on GitHub](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)),
so the `.h5ad` (or `.zarr.zip`) of a real atlas cannot live in an ordinary
repository. A Zarr store folder holds the same AnnData as many small files:

- **Every file is small.** `prepare()` and `cellucid_prepare()` write chunks
  of at most a few MB. The 561,947-cell Suo atlas is 1,961 files, 198 MB in
  all, the largest 1.6 MB.
- **It is read in place.** Both GitHub hosts that serve a repository's files
  (GitHub Pages and `raw.githubusercontent.com`) answer over HTTP/2, allow
  other sites to read them (CORS) and serve byte ranges, so Cellucid fetches
  only the chunks each view needs, up to 100 at once on one connection. A
  small file is also quick to fetch the first time: GitHub's edge fetches a
  large file whole before it answers the first request for any part of it
  (0.8–3.2 s measured for 16–66 MB files on `raw.githubusercontent.com`).
- **It costs nothing to serve.** Ordinary repository files use no Git LFS
  bandwidth (see [Git LFS](#git-lfs-an-alternative-with-a-quota)).

The written format is the standard AnnData Zarr layout that `anndata.read_zarr`
reads; its consolidated metadata lets a browser learn every array with one
request, since a web host cannot list a folder.

## Step by step

### 1. Write the store

From Python:

```python
import cellucid

cellucid.prepare(
    "atlas.zarr",
    format="zarr",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="my-atlas",
    dataset_name="My atlas",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
)
```

From R:

```r
cellucid_prepare(
  "atlas.zarr",
  format = "zarr",
  embeddings = list(X_umap = umap),
  dataset_id = "my-atlas",
  dataset_name = "My atlas",
  obs = obs,
  gene_expression = expr
)
```

The default layout is the one to publish: Zarr v2, the expression matrix
stored gene-major (CSC) in chunks of 65,536 values, every other array in
chunks of 262,144, byte-shuffled and deflated (the expression values
unshuffled), each chunk its own file, with consolidated metadata.
[How fast it reads](#how-fast-it-reads) explains why no other layout is faster
on GitHub. Every argument is described in
{doc}`../../python_package/c_data_preparation_api/index` and
{doc}`../../r_package/index`.

### 2. Add `.nojekyll` (for GitHub Pages)

A Zarr store keeps its metadata in dot-files (`.zmetadata`, `.zgroup`,
`.zattrs`, `.zarray`) and AnnData names each table's index `_index`. GitHub
Pages runs Jekyll by default, which publishes neither files whose names start
with `.` nor folders whose names start with `_`, so a store served from Pages
without the marker below is missing its metadata and its cell names. An empty
`.nojekyll` file at the root of the published branch turns Jekyll off:

```bash
touch .nojekyll
```

`raw.githubusercontent.com` serves every file regardless; the marker only
matters for Pages.

### 3. Commit and push

Use Git from the command line: the github.com upload page takes at most
25 MiB and 100 files at a time.

```bash
git add .nojekyll atlas.zarr
git commit -m "Add atlas"
git push
```

### 4. Open it in Cellucid

Paste either address into **Session › Open from link** and press **Connect**:

| Address | Served by | Example |
|---|---|---|
| `owner/repo/path/atlas.zarr` | `raw.githubusercontent.com`, the repository's default branch; works for any public repository | `my-lab/atlas-data/atlas.zarr` |
| `owner/repo@branch/path/atlas.zarr` | the same host, another branch, tag or commit | `my-lab/atlas-data@v1/atlas.zarr` |
| `https://owner.github.io/repo/path/atlas.zarr/` | GitHub Pages, once enabled (step 5) | `https://my-lab.github.io/atlas-data/atlas.zarr/` |

A `github.com/owner/repo/tree/<branch>/…` link to the folder works too. To
share it, put the address in a Cellucid link (`gh:` before the short form; a
Pages address URL-encoded):

```text
https://www.cellucid.com/?data=gh:my-lab/atlas-data/atlas.zarr
https://www.cellucid.com/?data=https%3A%2F%2Fmy-lab.github.io%2Fatlas-data%2Fatlas.zarr%2F
```

To try one first: the reference repository publishes a 72-cell dataset as a
Zarr folder,
[`gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr`](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr),
and the same dataset as one file,
[`gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad`](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad).

For several datasets in one repository, publish a dataset catalog
(`datasets.json` beside the stores) and link its folder instead:
{doc}`11_custom_dataset_repository`.

### 5. Optional: enable GitHub Pages

Pages answers a first request faster than `raw.githubusercontent.com` (below).
Open the repository's **Settings › Pages**, choose **Deploy from a branch**,
pick the branch and `/ (root)` (or `/docs`, and then leave `docs/` out of the
Pages address), and **Save**. After the deployment finishes, the store is at
`https://<owner>.github.io/<repo>/<path>/atlas.zarr/`; check that
`…/atlas.zarr/.zmetadata` opens in a browser tab.

## How fast it reads

Measured on 22 October 2026 from Munich (a line delivering 4.3 MB/s; GitHub's
nearest edge 19 ms away), with real requests to GitHub's hosts:

| | GitHub Pages | `raw.githubusercontent.com` |
|---|---|---|
| Protocol | HTTP/2 (no HTTP/3) | HTTP/2 (no HTTP/3) |
| Requests in flight on one connection | 100 | 100 |
| First request for a file (edge cache miss), median | 150 ms | 240 ms |
| A file the edge has cached (for 10 min on Pages, 5 min on raw) | 35–45 ms | 35–45 ms |
| First request for a large single file (16–66 MB), which the edge fetches whole | 0.7 s (25 MB) | 0.8–3.2 s |
| Uncached small files per second, 100 in flight | 416 | 317 |
| Bytes per second | as fast as the line | as fast as the line |
| Bursts of thousands of requests | one 503 in 3,000 (Cellucid asks again) | no refusal in 10,000 |

Each view therefore costs one or two round trips plus its bytes over your
line. For Suo (561,947 cells, 8,192 genes) written with the default layout
(198 MB), Cellucid's reader takes these times over links simulated from those
measurements (HTTP/2, 100 requests per connection, every CORS preflight
refused; measured 23 October 2026):

| What | Bytes | 150 ms round trip, 5 MB/s | 80 ms, 50 MB/s |
|---|---|---|---|
| Open (the store's structure) | 0.1 MB | 0.3 s | 0.2 s |
| Every `obs` field | 1.3 MB | 0.4 s | 0.2 s |
| 2-D UMAP | 3.7 MB | 0.9 s | 0.2 s |
| First gene | 0.1 MB | 0.2 s | 0.09 s |
| Next gene | 0.1 MB | 0.2 s | 0.08 s |
| 99 genes at once | 3.6 MB | 0.9 s | 0.2 s |
| Neighbour graph | 66.5 MB | 13.8 s | 1.8 s |
| Every gene in two cell types (a differential expression test) | 62.5 MB | 13.0 s | 1.7 s |

Over GitHub's HTTP/2 the folder reads as fast as one `.h5ad` or `.zarr.zip`
of the same data (46.1 s against 46.9 and 46.7 s for every step above and a
few more at 5 MB/s; 8.0 s against 8.3 and 8.2 s at 50 MB/s). Over HTTP/1.1,
which most other hosts speak, a browser fetches six files at a time and the
same folder took 46.8 s against 9.0 s at 50 MB/s: there, publish one file.

Everything but the first round trip is your line, not GitHub. A sharded Zarr
v3 store (many chunks per file) was measured too: it reads fewer files when a
view needs the whole matrix, but every gene from a shard not yet touched waits
for the shard's index first (0.35 s instead of 0.09 s at 80 ms), so the default
unsharded layout is the faster one for what Cellucid reads.

## Limits

| Limit | Value | What it means here |
|---|---|---|
| One file | blocked over 100 MiB, warned over 50 MiB | a Zarr store written by `prepare()` never comes close |
| Repository | ideally under 1 GB, strongly recommended under 5 GB | about the size of four Suo-sized atlases per repository |
| Published Pages site | at most 1 GB | everything in the published branch or folder counts |
| Pages bandwidth | 100 GB per month (soft) | about 400 complete Suo downloads; most visitors read far less |
| Pages builds | 10 per hour (soft) | batch your pushes |

Sources:
[large files](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github),
[GitHub Pages limits](https://docs.github.com/en/pages/getting-started-with-github-pages/github-pages-limits).
Beyond these, publish on a static object host or use {doc}`04_server_tutorial`.

## Updating a published store

The edges cache each file for up to 10 minutes (Pages) or 5 minutes (raw), so
for a few minutes after a push a reader may see old and new files side by
side. Write an update as a new store (`atlas-v2.zarr`) or on a new branch or
tag, and share that address; a link with `@<tag>` or `@<commit>` never
changes under its readers.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| "… no Zarr store is at this address …" | Nothing is published at that path: the owner, repository, branch or path is wrong, the repository is private, or the push has not finished. Open `https://raw.githubusercontent.com/<owner>/<repo>/HEAD/<path>/atlas.zarr/.zmetadata` in a browser tab to check. From a Pages address, the most common cause is a missing `.nojekyll` (step 2). |
| Opens from `owner/repo/...` but not from the Pages address | No `.nojekyll`: Pages left out the dot-files and `_index` folders. Add it (step 2) and wait for the deployment. |
| "… has no consolidated metadata (.zmetadata) …" | The store was written without consolidated metadata. Write it with `prepare()`/`cellucid_prepare()`, or run `zarr.consolidate_metadata` on it. |
| "… was not found at raw.githubusercontent.com (HTTP 404)" | A single `.h5ad` or `.zarr.zip` link names a file that is not there; check it as above. |
| Files on another branch are not found | Without `@branch` the repository's default branch is read; add `@branch`. |
| A file read returns a few lines of text | The file was committed with Git LFS: its pointer is in the repository, the file on the LFS host. Cellucid follows the pointer for a single `.h5ad` or `.zarr.zip`; a Zarr folder must not be stored in LFS. |
| Works at home, not on the institute network | A proxy or firewall blocks `raw.githubusercontent.com` or `github.io`; try the other address, or {doc}`04_server_tutorial`. |

More in {doc}`08_troubleshooting_data_loading`.

## Git LFS: an alternative with a quota

A single `.h5ad` or `.zarr.zip` over 100 MiB can be committed with
[Git LFS](https://git-lfs.com/). Cellucid opens it from the same link
(`owner/repo/path/atlas.h5ad`): it reads the pointer, then reads the file in
place from `media.githubusercontent.com`, which serves LFS content over HTTP/2
with ranges and CORS (a first read of a file there took about 0.5 s).

The cost is bandwidth. Every byte any visitor's browser reads counts against
the **repository owner's** LFS bandwidth: 10 GiB a month on GitHub Free, Pro
and Free for organizations, 250 GiB on Team and Enterprise Cloud. When it is
used up, GitHub either bills the overage or, without a payment method or with
a zero budget, disables Git LFS for the owner's account or organization until
the next month, for every repository it owns
([Git LFS billing](https://docs.github.com/en/billing/concepts/product-billing/git-lfs)).
One visitor reading Suo's neighbour graph and a few dozen genes is already
about 70 MB. Prefer a Zarr folder: it is read as fast and costs no LFS
bandwidth at all.
