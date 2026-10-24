# File format expectations

**Audience:** everyone (beginners get “what should I have?”, experts get exact file/key expectations)  
**Time:** 10–15 minutes  
**What you’ll learn:** what Cellucid expects of an `.h5ad`, an AnnData Zarr
store, what a file says about itself in `uns['cellucid']`, and the exact
`datasets.json` contract of a dataset catalog
**Prerequisites:** none (helpful: you know which loading option you plan to use)

---

## The four things Cellucid can load (practical view)

Cellucid reads **AnnData**. Every path (the file buttons, a dropped file, a
link, the sample list, a server, a notebook) points the viewer at one of these
**inputs**:

1) **AnnData `.h5ad`**
2) **AnnData Zarr store folder** (`.zarr`, Zarr v2 or v3)
3) **Zipped AnnData Zarr store** (`.zarr.zip`, or `.zip` holding one store)
4) **A dataset catalog**: a folder holding `datasets.json` and the AnnData
   files (1–3) it lists

A **Cellucid-ready AnnData** is not a fifth format: it is an ordinary
`.h5ad` or Zarr store that `cellucid.prepare()` (Python) or `cellucid_prepare()`
(R) wrote for fast reading (the expression matrix gene-major, arrays chunked
for range requests) with `uns['cellucid']` filled in. The exact contract is
{doc}`../../python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification`;
how to write one is in {doc}`../../python_package/c_data_preparation_api/index`
(Python) and {doc}`../../r_package/c_data_preparation_api/index` (R). anndata reads it like
any other file, and Cellucid reads any other file like it.

---

(uns-cellucid)=
## What a file says about itself: `uns['cellucid']`

An AnnData file may carry a mapping `uns['cellucid']` that tells Cellucid its
identity and the choices it opens with. Every key is optional, and a file
without the mapping opens exactly as any AnnData file does.

| Key | Value | Used for |
|---|---|---|
| `schema` | the integer `1` | the version of this mapping; a later version is refused, asking for a newer Cellucid |
| `id` | the dataset's stable id: 1–180 ASCII letters, digits, `.`, `_` or `-`, starting with a letter or digit, not ending with `.` | the dataset id sessions and community annotations key on ({doc}`06_dataset_identity_why_it_matters`) |
| `name` | text | the name shown in the dataset list and panel |
| `description` | text (may be empty) | the dataset panel's description |
| `citation` | a mapping `{text, url}` (`url` optional) | how to cite the dataset |
| `source` | a mapping `{name, url}` (`url` optional) | where the data comes from (the panel's **Source** and **URL**) |
| `defaults` | a mapping of launch defaults (below) | the choices the dataset opens with |
| `outlierQuantiles` | a mapping `{latentKey, minPoints, fields}`: the latent space and category size they were computed with, and one float32 array per categorical obs column (NaN where a cell has none), written by `prepare(latent_key=...)` / `cellucid_prepare(latent_key = ...)` | the outlier filter, which reads one column's array on its slider's first move instead of computing it from the latent space ({doc}`../e_filtering/02_outlier_filtering_per_active_field`) |
| `session` | written by `apply_cellucid_session_to_anndata()` in Python | never read by the web app |

Rules the reader applies:

- Text is stored as text; a one-element list of text is read as that text.
  It must not start or end with white space or hold a character
  with no glyph (a control or zero-width character).
- An address (`citation.url`, `source.url`) is an absolute `http(s)` URL.
- A key outside this table, or a wrong value, refuses the file with a message
  naming it, so a misspelt default never passes unnoticed.

`defaults` holds the launch defaults, each optional: `embedding1d`,
`embedding2d`, `embedding3d` (an `obsm` key, or `{key, columns}` with 0-based
columns), `knn` (an `obsp` key), `expression` (`X`, `raw/X` or
`layers/<name>`), `geneIdColumn` (a var column), `latentKey` (an `obsm` key),
`vectorFieldDefault` (a vector field id), `obsKeys` (a list of obs columns),
`normalizeEmbeddings` (true or false), `centroidMinPoints` (a positive
integer: the fewest cells a category needs to get outlier quantiles). These
are the launch choices `cellucid serve`, `serve_anndata()` and `show_anndata()`
take ({doc}`04_server_tutorial`).

How the file's own values combine with the route that opens it:

- **Id**: the file's id is the dataset id. A catalog lists the file under an
  id, which must equal the file's (else the file is refused); a server that
  names an id serves the file under it.
- **Name**: for a file you open yourself (from disk, a link or a catalog), the
  file's `name` is shown, else its `uns['title']` (as CELLxGENE writes it),
  else the file name. A server that names the dataset shows its name.
- **Description, citation, source**: the file's, whatever route opens it.
- **Launch defaults**: the file's, with every choice a server names on its
  command line taking precedence.

Beside `uns['cellucid']` and `uns['title']`, the one other part of `uns`
Cellucid reads is the colours Scanpy keeps for a categorical obs column,
`uns['<column>_colors']`: the column opens in them
({doc}`../d_fields_coloring_legends/03_color_by_behavior`). They are read with
the column's categories when the column is first shown: opening a file reads
no category list, only how many categories each column has, so a column of
one category per cell (barcodes stored as a categorical) costs nothing until
it is shown.

---

(dataset-catalog-contract)=
## A dataset catalog: `datasets.json`

A dataset catalog is a folder holding `datasets.json` and the AnnData files it
lists. The **Sample dataset** list reads one (www.cellucid.com reads
`https://theislab.github.io/cellucid-datasets/anndata/`), and so does a link to
a folder ({doc}`14_open_from_link`).

```json
{
  "schema": 2,
  "default": "pancreas",
  "datasets": [
    {
      "id": "pancreas",
      "path": "pancreas.zarr/",
      "session": {
        "path": "pancreas.cellucid-session",
        "bytes": 3422,
        "sha256": "14eeb50c90eacd56a78704b2dd03756fb8147cdc2dfc7754d4059d6f6d7753c1"
      }
    },
    { "id": "suo", "path": "suo.h5ad" }
  ]
}
```

- `schema` is the integer `2`. A catalog without it is refused, naming the
  three keys a catalog must hold.
- `default` is the id of the dataset the catalog opens on.
- `datasets` lists one dataset or more. Ids are distinct, and so are paths.
  No other key is read; an unknown key is refused.
- `id`: the dataset's stable id, with the same rule as `uns['cellucid']['id']`.
  When the file states an id, it must be this one.
- `path`: where the dataset is, relative to the folder: a **Zarr store folder**
  (ending in `/`), an **`.h5ad`** file, or a **zipped Zarr store** (`.zarr.zip`
  or `.zip`). Each part is written in canonical URL encoding (as
  `encodeURIComponent` writes it, e.g. `cell%20atlas.h5ad`), without `.`, `..`,
  a query or a fragment, and the path stays below the folder.
- `session` (optional): the dataset's **default view**, a published
  `.cellucid-session` file — a **Save** of the view passed through the publish
  step, which keeps the seven chunks a default carries and binds it to the
  catalog's dataset id (see {doc}`../l_sessions_sharing/04_official_sample_states`) —
  with its exact byte length
  (`bytes`) and lowercase SHA-256 (`sha256`). It is applied when the dataset
  is opened from the catalog ({doc}`../l_sessions_sharing/04_official_sample_states`).
  The length and digest pin it to the catalog: a session file that differs is
  refused rather than applied.

The dataset list shows each dataset's name and cell count, read from its
file's metadata and `uns` alone (no matrix, column or index is read); the
dataset's own file is opened when it is chosen. A Zarr store folder on a web
host needs consolidated metadata (`.zmetadata` for Zarr v2, the root
`zarr.json` of v3), which `prepare()` writes.

To compute `bytes` and `sha256` of a session file:

```python
import hashlib
from pathlib import Path

data = Path("pancreas.cellucid-session").read_bytes()
print(len(data), hashlib.sha256(data).hexdigest())
```

The worked public example is {doc}`11_custom_dataset_repository`; the shorter
publishing workflow is {doc}`02_local_demo_tutorial`.

---

## AnnData H5AD and Zarr (what Cellucid reads)

The browser reads one H5AD file, one ZIP containing a complete Zarr store
(Zarr v2 or v3, sharded or not), a Zarr store folder from your disk (the
**Open… › Zarr folder…**, or a dropped folder), or a Zarr store folder on a web
host (through its consolidated metadata). Server mode and Jupyter accept H5AD paths, Zarr
store folders (v2 or v3), zipped Zarr stores, or an in-memory AnnData object.
The server serves a path as it is on disk, and the browser reads it with the
same reader; an in-memory object is written once to a temporary H5AD.

### Compression the browser reads

- **H5AD:** uncompressed, gzip (with or without shuffle), LZF (h5py's own
  filter) and Fletcher-32 checksums. A dataset compressed with a filter outside
  these (Blosc, Zstandard, LZ4, bitshuffle, bzip2, SZIP, N-bit, scale-offset) is
  refused when it is read, naming the filter; the rest of the file still opens.
- **Zarr:** uncompressed, gzip, zlib, Blosc (blosclz, LZ4, LZ4HC, zlib,
  Zstandard; byte and bit shuffle), LZ4, Zstandard (Zarr v3's default), CRC-32C,
  variable-length strings, and Zarr v3 sharding. An array using another codec
  (bz2, LZMA, Snappy, ZFP, delta, fixed-scale-offset, …) is refused when it is
  read, naming the codec.
- **Served from Python:** `cellucid serve`, `serve_anndata()` and
  `show_anndata()` check every array's filters and codecs before they start.
  A store holding an array the browser cannot decode is read once with
  anndata (and its codec plugins) and served as a decodable temporary H5AD;
  the terminal names the arrays.

### Minimum requirements (every loading route)

You must have at least one `obsm` array Cellucid can draw: dense, numeric, one
row per cell, at least one column (a one-dimensional array is one column), with
finite values that fit in 32-bit floats and not the same coordinates for every
cell.

Which array each dimension shows by default:
- `obsm["X_umap_1d"]`, `obsm["X_umap_2d"]` and `obsm["X_umap_3d"]` fill their
  own dimensions; a suffixed key and its column count must agree exactly.
- When a file declares none of the three and carries the plain `X_umap` that
  `sc.tl.umap()` writes, that array fills the dimension of its width (1, 2 or 3
  columns), so an ordinary Scanpy object loads as written.
- When the file has no UMAP at all, a deterministic **best guess** fills the
  dimensions instead (t-SNE, force layouts, spatial, …, and last the first two
  columns of a wide array such as `X_pca`). See {doc}`13_data_mapping`.

This holds for every loading route: the browser file picker, `cellucid serve`,
and the notebook viewers. In the browser any other array can then be turned
on under **Session › Data mapping** and shown in any view; `cellucid serve` takes
`--embedding-1d/--embedding-2d/--embedding-3d KEY[:c0,c1]`.

### Files from older AnnData releases

Files written before anndata 0.8 (2022) open as written, from the file pickers,
a dropped file, a link and `cellucid serve` alike:

- **anndata 0.7**, and every file SeuratDisk's `Convert()` writes: dataframes
  whose categorical columns keep their categories under `__categories`, arrays
  without encoding attributes, and (Zarr) a store without a root `.zattrs`;
- **anndata 0.6 and earlier**: `obs` and `var` as one record per cell
  (compound datasets) with categories in `uns`, a compound `obsm`, `raw` stored
  as `raw.X` / `raw.var`, and sparse matrices named by `h5sparse_format`.

Every column a file stores is shown, including those a SeuratDisk file leaves
out of its `column-order` (anndata itself drops those). Ordered categories and
missing values keep their meaning. Nothing needs converting first.
`prepare()` and `show_anndata(adata)` take an AnnData you have already read, so
they hold the columns anndata read; pass the file's path to `cellucid serve`,
`serve_anndata` or `show_anndata` to keep the unlisted ones.

### Optional (but commonly expected)

- **Obs fields** (metadata for coloring/filtering): `adata.obs[...]`
- **Gene metadata** (for gene search / naming): `adata.var[...]` and/or `adata.var.index`
- **Gene expression**: `adata.X` (dense or sparse)
- **Connectivity**: `adata.obsp["connectivities"]` (KNN graph)

### Vector fields (optional overlays)

Vector fields are optional per-cell displacement vectors (e.g. RNA velocity, drift) that Cellucid can render as an animated overlay.

If you expect vector fields, the key requirements are:
- **Correct width**: exactly as many columns as the embedding the vectors are drawn over.
- **Correct row order**: vectors must match the same cell order as points/embeddings.
- **Discoverable naming / metadata**: so the viewer can find the fields.

See the overlay UI docs:
- {doc}`../i_vector_field_velocity/index`

#### AnnData naming convention (`obsm`)

Cellucid discovers vector fields from `adata.obsm` keys using scVelo's
convention:

- `<field>_<basis>`: the field `<field>` drawn over the embedding `X_<basis>`

Examples:
- `velocity_umap` over `X_umap`, `velocity_umap_3d` over `X_umap_3d`
- `T_fwd_umap`, `T_bwd_umap` over `X_umap`
- `velocity_pca` over `X_pca`, `velocity_branches` over `X_branches`

Notes:
- A key that starts with `X_` (`X_umap`, `X_umap_3d`) is an embedding, never a
  vector field.
- The key is split at one `_` into a field and a basis whose `X_<basis>` the
  file holds; when several splits match, the longest basis wins. A key with no
  match is an ordinary `obsm` array, so `velocity_umap_2d` beside a bare
  `X_umap` (no `X_umap_2d`) is not a vector field: name it `velocity_umap`.
- A field is bound to every enabled view of `X_<basis>` when both arrays have
  the same number of columns. The keys of one field on every embedding are one
  field: `velocity_umap` and `velocity_umap_3d` are both `velocity`.
- The overlay dropdown shows the fields bound to the embedding the view
  **shows**; a view showing another embedding (for example a t-SNE) has none.

#### How to confirm vector fields exist (quick checks)

- Print `adata.obsm.keys()` and look for keys like `velocity_umap` beside `X_umap`.

Naming and schema details:
- {doc}`../i_vector_field_velocity/01_what_vector_fields_are_user_facing`

---

## Common “format mismatch” problems (and what they mean)

### “Cellucid found nothing to plot”

Meaning:
- No `obsm` array can be drawn: every entry is a table, a sparse array, has
  another row count, or holds a value that cannot be drawn. The message lists
  each entry and why.

Fix:
- Add an embedding (for example `sc.tl.umap(adata)`) or fix the array the
  message names. A file without a UMAP but with any drawable array opens on a
  best guess ({doc}`13_data_mapping`).

### “Fields list is empty”

Meaning:
- `obs` is empty, or `uns['cellucid']['defaults']['obsKeys']` names only some columns.

Fix:
- Ensure `adata.obs` actually contains columns and is not all missing.

### “Gene search returns nothing”

Meaning:
- There is no gene expression (`X`), or the genes were named from a `var` column whose vocabulary nobody is searching in.

Fix:
- Provide `adata.X` and name the genes from the `var` column your readers will type.
- In Jupyter/server mode, pass `gene_id_column=...` if needed.

---

## Edge cases and limits (read if you’re debugging weirdness)

- **NaN/Inf in embeddings**: refused, naming the cell, before the dataset is
  shown.
- **Duplicate gene names**: made unique as AnnData's `var_names_make_unique()`
  does (`MATR3`, `MATR3-1`), and reported.
- **Huge categorical fields**: any number of categories loads and the legend scrolls, but palettes cannot keep 50k–100k colours distinct; consider aggregating.
- **Browser memory**: there is no fixed size ceiling for an `.h5ad` or Zarr
  file. Each array read is held to the size its metadata declares, and one the
  browser cannot allocate is refused as too large, naming the bytes. A server
  hands the file to the browser, which reads only what it asks for.
- **Byte order**: never a compatibility problem. H5AD is byte-order-agnostic,
  because HDF5 converts to host order while reading; Zarr records its own byte
  order per array and the reader swaps when it has to.

---

## Next steps

- Pick a workflow: {doc}`01_loading_options_overview`
- Load locally without Python: {doc}`03_browser_file_picker_tutorial`
- Serve an `.h5ad` or Zarr store from Python: {doc}`04_server_tutorial`
- Publish a public dataset catalog: {doc}`11_custom_dataset_repository`
- Debug loading failures: {doc}`08_troubleshooting_data_loading`
