# How the viewer reads AnnData (`show_anndata()` / `serve_anndata()`)

This page explains how Cellucid visualizes **AnnData directly**: any `.h5ad`,
Zarr store, or in-memory object, prepared with `prepare()` or not.

Typical uses:
- you’re iterating inside an analysis notebook,
- you want a quick “does this dataset look sane?” view,
- or you have a large `.h5ad` or Zarr store and want to view it in place.

For the fastest reads and the easiest sharing, write a Cellucid-ready AnnData
once with `prepare()` and open that file here:
{doc}`../c_data_preparation_api/index`.

## At a glance

**Audience**
- Computational users working with Scanpy/AnnData
- Wet lab users who received a `.h5ad` or Zarr store from a collaborator

**Prerequisites**
- `pip install cellucid`
- A `.h5ad` file, a Zarr store (a folder, or zipped as `.zip`), or an
  in-memory `AnnData`

## How it works

The viewer reads the AnnData store itself, in place, with HTTP range requests:
only the bytes a view needs are transferred, and nothing is converted first.
Python computes nothing per request. It names the dataset, serves
the store's bytes under `/_cellucid/anndata/<dataset_id>.h5ad` (or
`.zarr.zip`, or `.zarr/<key>` for a store folder), states the launch choices
you passed, and hosts the viewer channel — over the notebook's own connection
for `show_anndata()` (no port), over HTTP and a WebSocket for
`serve_anndata()`, which also serves the viewer's web build. The
route is named by the dataset id, so the local file name is never published.

- An `.h5ad` file, a Zarr store folder (Zarr v2 with `.zgroup`, or v3 with
  `zarr.json`), and a zipped Zarr store written by `prepare()` (cells in
  priority order) are served as they are on disk, never copied.
- Any other AnnData file, and an in-memory `AnnData`, are served as a Zarr
  store with their cells in **priority order** (the layout `prepare()`
  writes), each chunk made from the source when the viewer asks for it: a
  subsample the viewer opens is then the first rows of every array and the
  first entries of every gene, so only its cells are sent, whatever the
  source's layout. An in-memory object is copied once, before the call
  returns, into memory-mapped files removed when the server stops or Python
  exits (its sparse expression matrices laid out gene-major, so a gene is one
  range: see below); a file is read where it is, its matrices keeping their
  form (an array stored in chunks -- compressed, or a sparse matrix anndata
  wrote -- is decoded once into memory-mapped files beside the store as its
  rows are first read, since rows in priority order lie all over the file;
  removed when the server stops). Original positions are the ones `prepare()` gives: an object read
  from a file in priority order, its `uns['cellucid']['priority']` kept and
  still describing its cells, is put back in its original order first; any
  other object or file is taken in its own row order, and a
  `cellucid_original_position` column with no record that describes the
  cells (a subset, a reordering, a column alone) names nothing and a warning
  says so -- so selections and labels pulled back land on the same cells.
- An AnnData store holding an array stored with a codec the browser cannot
  decode (an HDF5 filter other than gzip, shuffle, Fletcher-32, or LZF, such
  as Blosc, Zstandard, szip, or scaleoffset; a Zarr codec such as bz2, lzma,
  or zfpy) is decoded on the kernel through that store, chunk by chunk as the
  viewer asks. The terminal names those arrays (`Browser cannot decode: …`,
  then `Serving: a Zarr store in priority order, made as the browser reads
  it`). A MuData, Loom or 10x file is decoded on the kernel the same way.
- A MuData file (`.h5mu`, MuData Zarr) not written by `prepare()`, a Loom file
  and a 10x Cell Ranger `outs/` folder are served in priority order too
  (`Serving: a Zarr store in priority order, …`; a MuData Zarr store for a
  MuData and for a 10x matrix of several feature types), holding what the
  viewer presents from the file itself, so a subsample transfers only its
  cells. A Loom file with no embedding, and a 10x `.h5` alone, are refused
  saying what to do.

## Quick start

### Notebook (recommended)

```python
from cellucid import show_anndata

viewer = show_anndata(
    "data.h5ad",
    height=600,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# A Zarr store (folder or .zip) or an in-memory AnnData works the same way.
```

### Terminal

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1
```

### Python (script)

```python
from cellucid import serve_anndata

serve_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
)  # blocks
```

`dataset_name` and `dataset_id` default to the store's own
`uns['cellucid']['name']` and `uns['cellucid']['id']`. An object without them
needs both, and serving it without them stops with:

```text
Error: dataset_name and dataset_id are required: this AnnData declares no Cellucid identity of its own (uns['cellucid']['id'] and ['name']). Name the dataset, for example dataset_name="My study", dataset_id="my-study-v1" (--dataset-name/--dataset-id on the command line).
```

Sessions and community annotations are keyed by the dataset id.

## Minimum AnnData requirements (the “must haves”)

### 1) An embedding in `adata.obsm`

Cellucid shows one embedding per dimension slot (1D, 2D, 3D) and fills the
slots from `adata.obsm`, in this order:

1. **The slots you name.** `embedding_keys={"2d": "X_tsne"}` (or
   `--embedding-2d X_tsne`) shows `obsm['X_tsne']` in 2D. A wider array is read at
   the 0-based columns you name: `embedding_keys={"2d": ("X_pca", [0, 1])}`, or
   `--embedding-2d X_pca:0,1`.
2. **The default rule, for every slot you did not name.** `X_umap_1d`,
   `X_umap_2d`, and `X_umap_3d` each fill their own slot. When none of them is
   present, one UMAP under its conventional name fills the slot of its own
   column count: the first 1- to 3-column one of `X_umap` (the key
   `sc.tl.umap()` writes), `X_UMAP` (a SingleCellExperiment's `UMAP`),
   `X_wnn.umap` and `X_ref.umap` (a Seurat object's). So the ordinary Scanpy
   object is shown exactly as Scanpy left it, and an object converted from R
   as R named it.
3. **A deterministic best guess, only when neither filled any slot.** Cellucid
   then ranks the other `obsm` arrays by name (UMAP, then t-SNE, force-directed
   layouts, PHATE, PaCMAP/TriMAP, spatial, diffusion maps, other names, and
   latent spaces such as PCA or scVI last) and fills each slot with the
   best-ranked array of exactly that width. If no array has 1 to 3 columns, the
   best-ranked wider one is shown in 2D at columns 0 and 1. The viewer marks a
   guessed slot as guessed under **Session › Data mapping**.

The viewer applies these rules when it opens the store, the same rules it
applies to a file opened in the browser. If nothing in `obsm` can be plotted,
the viewer stops with a message that lists what `obsm` holds and why each entry
was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (3,696) and at least one numeric column. obsm holds: X_df (table), spatial (3,696 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again.
```

Those are the embeddings each dimension shows first, its *defaults*. The
viewer's data mapping also *enables* every other usable `obsm` array 1 to 3
columns wide beside them; see
{doc}`../../web_app/b_data_loading/13_data_mapping`.

The complete rules are in the section on which `obsm` keys are read, below.

### 2) (Optional but common) Gene expression in `adata.X`

- If `adata.X` is missing, Cellucid can still show embeddings + obs fields.
- Gene search / expression coloring will not work without expression values.

### 3) (Optional) Cell metadata in `adata.obs`

Cellucid derives obs fields automatically:
- numeric → continuous
- categorical/string/bool → categorical

Anything else — a `datetime64` collection date, a timedelta, a period — cannot
be shown, and the viewer refuses it naming that column
(`Observation field '<key>' has unsupported dtype '…'`). Either convert it
(`adata.obs["collection_date"].astype(str)`) or name the columns you want:

```python
show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    obs_keys=["cell_type", "leiden", "total_counts"],
)
```

`obs_keys` is accepted by {func}`~cellucid.show_anndata`,
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataViewer`,
{class}`~cellucid.AnnDataServer`, and `cellucid serve … --obs-key`. It also
fixes the order fields appear in, and a column left out is not shown. Every
name must be an `obs` column of the store; one that is not stops the server
before it binds, listing the columns `obs` holds.

## Which `obsm` keys are read (detailed, for debugging)

### Slots you name: `embedding_keys` and `--embedding-<n>d`

`embedding_keys` maps a slot — `"1d"`, `"2d"`, or `"3d"` — to an `obsm` key, or
to a `(key, columns)` pair whose columns are 0-based and one per axis:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    embedding_keys={"2d": "X_tsne", "3d": ("X_pca", [0, 1, 2])},
)
```

```bash
cellucid serve /path/to/data.h5ad --dataset-name "My study" --dataset-id my-study-v1 \
  --embedding-2d X_tsne --embedding-3d X_pca:0,1,2
```

A key without columns is read at columns `0 … d-1`, so a 2D slot reads the first
two columns of any array at least two columns wide. A named array that does not
exist, is not a 2-D array, is too narrow for its slot, or repeats a column stops
the server before it binds, with that reason — the check reads the store's
metadata only. The viewer checks the values when it reads them, and refuses a
selected column holding a non-finite value. A slot you name replaces that
slot's default; slots you do not name still follow the default rule.

### The default rule

| `obsm` key | Slot it fills |
| --- | --- |
| `X_umap_1d` | 1D, columns `[0]` |
| `X_umap_2d` | 2D, columns `[0, 1]` |
| `X_umap_3d` | 3D, columns `[0, 1, 2]` |
| `X_umap`, `X_UMAP`, `X_wnn.umap`, `X_ref.umap` (the first that is 1, 2, or 3 columns wide), only when none of the three is present | the slot of its own width |

Any combination of the three dimensional keys may be present, and each fills its
own slot. If any of them exists, `X_umap` never joins them: an object holding
both `X_umap_3d` and `X_umap` serves 3D only. Of the single UMAPs only the
first fills a slot: Scanpy's `X_umap` comes before a SingleCellExperiment's
`X_UMAP`, and another spelling (`X_Umap`) is no default, only the best guess's
first choice. `sc.tl.umap(adata)` writes two
columns and `sc.tl.umap(adata, n_components=3)` writes three, so both are
readable as written. A default key that cannot be plotted — a non-finite
coordinate, or no spread at all — is refused rather than skipped, because the
key says the array was meant to be drawn.

### The best guess

When neither `embedding_keys` nor the default rule fills any slot, every usable
`obsm` array — dense, numeric, one row per cell, at least one column — is a
candidate, except vector fields (`<field>_<basis>` beside an `X_<basis>`, such
as `velocity_umap` beside `X_umap`) and the `latent_key` array.
Candidates are ranked by the first of these name patterns their lowercased key
contains, ties broken by key in code-point order:

| Rank | Key contains |
| --- | --- |
| 0 | `umap` |
| 1 | `tsne`, `t-sne`, `t_sne` |
| 2 | `draw_graph`, `force` |
| 3 | `phate` |
| 4 | `pacmap`, `trimap` |
| 5 | `spatial` |
| 6 | `diffmap` |
| 7 | none of these |
| 8 | `pca`, `scvi`, `latent`, `harmony`, `x_emb` |

For each of 3D, 2D, and 1D, the slot takes the best-ranked candidate whose width
is exactly that dimension. Only when no candidate of width 1 to 3 fills any slot
is the 2D slot filled from the best-ranked wider array at columns `[0, 1]` — so a
10-column `X_umap` is shown as `X_umap[0,1]`. A candidate whose selected columns
turn out to be non-finite or flat is skipped for the next one. The ranking is
deterministic, the same object always resolves the same way, and the same cases
are pinned for every implementation in
`tests/fixtures/data_mapping_conformance.json`.

### Nothing in your object is renamed or written

Resolution is a read. Cellucid does not assign anything into `adata.obsm` and
does not write to the `.h5ad` or Zarr store you served: the server only reads
its bytes. **Session › Data mapping** names the array each slot shows, under
its embedding id (`X_umap`, or `X_pca[0,1]` when only some columns are read).

### When nothing can be plotted

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (3,696) and at least one numeric column. obsm holds: X_df (table), spatial (3,696 × 2, NaN at cell 17). Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again.
```

Each entry says what the array is and, for a dense array, why it was not used:
the first non-finite cell, a row count that is not the cell count, a dtype that
is not numeric, or its role (`a vector field`, `the latent space`). `obsm holds:
nothing.` means `obsm` is empty. Compute an embedding, or name one with
`embedding_keys`.

### Where the same rule applies

The viewer applies this resolution to every AnnData store it reads: one handed
over by {func}`~cellucid.show_anndata`, {func}`~cellucid.serve_anndata`,
{class}`~cellucid.AnnDataViewer`, {class}`~cellucid.AnnDataServer`, or
`cellucid serve` on a path it detects as direct AnnData input, and one opened in
the browser. {func}`~cellucid.add_transition_drift_to_obsm` uses the default
rule only, because it computes drift in a UMAP.

```{important}
A file written by `prepare()` is read by this same resolution. Its
`embedding_defaults` are stored as `uns['cellucid']['defaults']`, which the
viewer applies before its default rule, and an `embedding_keys` you pass here
overrides them for this server. See
{doc}`../c_data_preparation_api/03_embeddings_and_coordinates`.
```

## Loading and memory behavior (critical for large datasets)

### A path: `.h5ad`, Zarr store folder, or zipped Zarr store

The server reads the store's metadata once, before it binds (cell and gene
counts and the names of the stored elements — never an array), and then only
serves bytes. It holds no matrix in memory. The browser reads what each view
needs: an embedding, an `obs` column, a gene, the neighbour graph when a view
draws it. A row-major (CSR) expression matrix is indexed by gene once, in the
browser, the first time a gene is read; a column-major (CSC) matrix needs no
index.

### In-memory `AnnData`

An in-memory object is served as a store in priority order, made from it as
the viewer reads it; your object is not modified. A cell-major (CSR) sparse
`X` is laid out gene-major once, so one gene is one contiguous read.

The viewer shows the object **as it was when you called** `show_anndata()` /
`serve_anndata()`: before the call returns, the object's arrays are copied
once into memory-mapped files in the system temporary folder — every sparse
expression matrix (`X`, layers, `raw.X`) gene-major in priority order, every
dense one in priority order, the `obsm`/`obsp` matrices — so what the viewer
reads is one consistent snapshot. Replacing an element afterwards
(`adata.X = …`) or changing an array in place
(`np.log1p(adata.X.data, out=adata.X.data)`) does not show, and never mixes
into what is served; view the object again to show a change. The copy costs
a pass over each array (about two seconds per 100 million stored values) and
the disk space of the arrays, never a second copy in memory; the files are
removed when the server stops or Python exits, and by the next server when
the kernel was killed first. The matrices are stored as `prepare()` stores
them, by the same functions: float32 unless integer counts need float64,
duplicate coordinates summed in stored order, stored zeros dropped; a value
that is not finite (NaN, infinity) is refused, saying which matrix. The
measurements are in {doc}`14_performance_scaling_and_lazy_loading`.

## Gene IDs, expression source, and repeated names

Cellucid needs a mapping from “what you type in the gene search box” → a column index in the matrix.

By default, expression is read from `adata.X` and genes are named by
`var.index`. An object without `X` reads `adata.raw.X`, else its only layer,
else its best-named layer (`norm` or `log` first, `counts` and `raw_counts`
last), and reports that choice as guessed — the rule a file opened in the
browser follows. `expression_source` (`--expression`) chooses the matrix
explicitly: `"raw"` reads `adata.raw.X` and names genes from `adata.raw.var`,
and `"layer:<name>"` reads `adata.layers[<name>]`, which shares `var`. A source
the store does not have stops the server before it binds. The viewer can switch to any
other matrix or gene-name column afterwards, under **Session › Data mapping ›
Genes** (see {doc}`../../web_app/b_data_loading/13_data_mapping`).

`gene_id_column` (`--gene-id-column`) names genes from a column of that gene
axis instead of its index. Any non-blank string names that exact column, so
`gene_id_column="index"` selects a column literally named `"index"`:

```python
viewer = show_anndata(
    "data.h5ad",
    gene_id_column="gene_symbols",
    expression_source="raw",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

Repeated and empty names are accepted. An empty value (`""`, `None`, or NaN)
takes that row's index name, and a repeated name is made unique the way
anndata's `var_names_make_unique()` does it — the second `MATR3` is shown as
`MATR3-1`. The viewer reports what changed, for example
`Gene names: 1 gene has no 'gene_symbols' value and uses its var index name. 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.`
The rule is the one `prepare()` and `cellucid_prepare()` report when they
write a file; see
{ref}`python_package-make-unique-gene-names`.

`HLA-DRB1/2` and `Gene A` are ordinary names, and `Gene` and `gene` are two
distinct genes.

## Neighbour graph and vector fields: read when drawn

The neighbour graph and the vector-field overlay cost nothing at startup. The
viewer reads the graph only when a view draws it, and a vector field only when
the overlay shows it. The launch options choose what they open on:

| Capability | Reads | Python | Terminal |
|---|---|---|---|
| Neighbor graph | the `adata.obsp` matrix named, else `obsp['connectivities']` | `knn_key="connectivities"` | `--knn connectivities` |
| Vector fields | `adata.obsm['<field>_<basis>']` over `adata.obsm['X_<basis>']` | `vector_field_default="velocity"` | `--vector-field-default velocity` |

Neither choice limits the viewer: **Session › Data mapping** can choose any
`obsp` graph and bind any `obsm` array that fits an embedding (see
{doc}`/user_guide/web_app/b_data_loading/13_data_mapping`).

`prepare()` takes its own `connectivities=` and `vector_fields=` arguments. See
{doc}`../c_data_preparation_api/07_connectivities_knn_graph` and
{doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Vector fields (velocity/drift overlays)

Cellucid can render per-cell displacement vectors as an animated overlay. The
viewer binds the fields below by their keys, reads one when the overlay shows
it, and any other `obsm` array can be bound under
**Session › Data mapping › Vector fields**.

### Naming convention (scVelo's `<field>_<basis>`)

A vector field is an `adata.obsm` key `<field>_<basis>` beside the embedding
`adata.obsm['X_<basis>']` it is drawn over, for any basis — the names scVelo
writes:

- `velocity_umap` (over `X_umap`, as many columns as it)
- `velocity_umap_3d` (over `X_umap_3d`)
- `T_fwd_umap`, `T_bwd_umap` (CellRank-style drift over `X_umap`)
- `velocity_pca` (over `X_pca`), `velocity_branches` (over `X_branches`)

Rules:
- a key that starts with `X_` is an embedding, never a vector field; any other
  key is split at one `_` into a field and a basis whose `X_<basis>` is in the
  same `obsm`, and the longest basis wins (exact, case-sensitive match; see
  {doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`).
  A key with no such split is an ordinary array: `velocity_umap_2d` beside a
  bare `X_umap` (no `X_umap_2d`) is not a vector field — name it
  `velocity_umap`
- a field is bound to every enabled view of `X_<basis>` when both arrays have
  the same number of columns, and read at the view's columns (`velocity_pca`
  over `X_pca` shown at columns `[0, 1]`); a vector drawn over another layout
  points the wrong way, so no other binding is made by default
- vectors must match `n_cells` and have as many columns as the embedding they
  are bound to
- vectors are scaled by that embedding's own normalization factor, into the same
  coordinate space as its points
- `vector_field_default` names the field the overlay opens on by its field
  id, the `<field>` part (`velocity` for `velocity_umap` and
  `velocity_umap_3d`, which are one field), when the object declares more
  than one; pass it to `show_anndata(...)`,
  `serve_anndata(...)`, `AnnDataViewer(...)`, or `AnnDataServer(...)`. A name
  `obsm` does not declare stops the server before it binds, listing the fields
  it does declare; so does a field none of whose embeddings the viewer shows
  when the store opens (the embeddings chosen or found by the default rule,
  and every other 1-3 column array that is not the latent space), since the
  overlay could not draw it
- an object that declares no field by this rule simply has no overlay

For example, an object containing both `velocity_umap` and
`T_fwd_umap` needs:

```python
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    vector_field_default="velocity",
)
```

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The Pancreas sample in 2D with the velocity overlay running over the cells and the Vector field controls open in the sidebar.
:width: 1440px

The 3,696-cell Pancreas sample in 2D (**Planar**) with its velocity field
drawn as particle trails; the **Vector field** controls (field, density,
speed, trail, size, opacity, colour scheme, Sync with LOD) are open beside it. A `velocity_umap` key is the field
`velocity`.
```

## Connectivity (KNN graph)

On the direct-AnnData path the viewer opens on `obsp['connectivities']` when
the store has it, and reads it only when a view draws it. `knn_key` names
another `obsp` matrix:

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
```

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1 \
  --knn connectivities
```

`knn_key` is accepted by {func}`~cellucid.show_anndata`,
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataViewer`, and
{class}`~cellucid.AnnDataServer`. From the terminal it is `--knn KEY`. A graph written by
`sc.pp.neighbors(adata, key_added="pca15")` is `knn_key="pca15_connectivities"`.

### A matrix that is not there, or not a graph, is an error

Naming a matrix the store does not have stops the server before it binds,
rather than serving a dataset silently missing the thing you asked for:

```text
Error: knn_key names obsp['connectivities'], which this store does not have. obsp holds: []. Compute the neighbour graph with sc.pp.neighbors(adata), or leave knn_key unset.
```

The viewer checks the matrix when it reads it, and never transforms it.
Scanpy's `obsp['distances']` is a directed k-nearest-neighbour matrix, so the
viewer refuses it with that reason:

```text
obsp['distances'] is not symmetric — choose a connectivities matrix.
```

### What the matrix must satisfy

- the matrix must already be square, finite, non-negative, exactly symmetric
  in topology and weight, and zero on the diagonal;
- sparse inputs must not store explicit zeros or duplicate coordinates;
- Cellucid preserves exact Float64 weights and rejects asymmetric, directed,
  negative, non-finite, or self-edge-containing graphs rather than
  transforming them;
- a valid zero-edge matrix is still an explicitly present graph,
- large graphs can be expensive to render in the browser.

```{note}
`prepare()` stores one graph, `obsp['connectivities']`, from its
`connectivities=` argument; the viewer's default rule finds it, so a prepared
file needs no `knn_key`.
```

## What startup prints

The server runs four numbered steps, from `cellucid serve` and
from {func}`~cellucid.serve_anndata` alike. None of them reads an array, so
startup takes about the same time for a small object and an atlas:

```text
$ cellucid serve /path/to/pbmc_demo.h5ad \
      --dataset-name "PBMC demo" --dataset-id pbmc-demo --no-browser \
      --knn connectivities --embedding-2d X_umap:1,0
[1/4] Detecting format...
      Path: /path/to/pbmc_demo.h5ad
      Format: AnnData (.h5ad)
      ✓ Format detected

[2/4] Reading metadata...
      Cells: 4,096
      Genes: 512
      obsm: X_umap
      obs fields: 3
      obsp: connectivities
      layers: none
      ✓ Metadata read

[3/4] Checking launch choices...
      embedding2d: {'key': 'X_umap', 'columns': [1, 0]}
      knn: connectivities
      ✓ Launch choices checked

[4/4] Starting server...
```

Reading the transcript:

- Step 1 names the format and layout: `AnnData (.h5ad)`,
  `AnnData Zarr, zipped`, or `AnnData Zarr, Zarr v2 folder` /
  `AnnData Zarr, Zarr v3 folder` (other formats by their own label, such as
  `MatrixMarket folder, a folder of 3 files`). A store holding an array the
  browser cannot decode adds `Browser cannot decode` and `Serving` lines here.
- Step 2 lists what the store holds, from its metadata.
- Step 3 lists each launch choice you set, as the viewer receives it; a choice
  you left unset is not listed, and the viewer takes it from the store's
  `uns['cellucid']['defaults']`, else its default rule. A Zarr store folder
  without consolidated metadata adds `Consolidated metadata: built in memory`.
- Step 4 establishes the viewer's web build and binds the socket. The full
  transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

### The URLs in the banner

`Local URL` and `Viewer URL` are the origin a browser on the serving machine
opens. When the server is bound to every interface — `--host 0.0.0.0`, `::`, or
an empty host — those two lines still show the loopback origin, because a
wildcard is a statement about which interfaces accept connections and not an
address anything can open. The addresses another machine can use are printed
separately:

A wildcard bind prints the loopback origin plus the machine's own name; see {doc}`12_remote_servers_ssh_tunneling_and_cloud`.

The second block appears only when this machine has a name that resolves to
something other than loopback. A non-loopback bind has no authentication, so
reach it through a tunnel rather than exposing the port:
{doc}`12_remote_servers_ssh_tunneling_and_cloud`, and, for the compute-node
case where the tunnel has to terminate on a login node,
{doc}`17_hpc_slurm_and_compute_node_serving`.

## Troubleshooting

### “dataset_name and dataset_id are required”

**Likely cause**
- the store declares no `uns['cellucid']['id']` and `['name']`

**Fix**
- pass `dataset_name=` and `dataset_id=` (`--dataset-name` / `--dataset-id`)

### “Cellucid found nothing to plot”

**Likely causes**
- you didn’t compute an embedding yet
- every `obsm` array is a table, sparse, non-numeric, or has a non-finite
  coordinate in the columns that would be plotted
- `n_cells` mismatch between `adata.obsm[...]` and `adata.n_obs`

**How to confirm**
- read the `obsm holds:` list: each entry names the array and why it was not used

**Fix**
- compute an embedding (`sc.tl.umap(adata)`), repair the named array, or name a
  usable one with `embedding_keys={"2d": "<key>"}` / `--embedding-2d <key>`

### The wrong embedding is shown

**Likely cause**
- the object has no `X_umap*` key, so a slot was guessed; **Session › Data
  mapping** marks it as guessed

**Fix**
- name the array you want with `embedding_keys` or `--embedding-<n>d`

### “knn_key names obsp[…], which this store does not have”

**Likely causes**
- the neighbor graph was never computed, was computed on a different object, or
  was stored under another key (`sc.pp.neighbors(..., key_added=...)`)

**Fix**
- run `sc.pp.neighbors(adata)` before serving (and re-save the file if you serve
  a path), name the key the object holds, or leave `knn_key` unset

### “obsp[…] is not symmetric — choose a connectivities matrix”

**Likely cause**
- the named matrix is a directed matrix such as `distances`

**Fix**
- name the symmetric `connectivities` matrix instead

### The viewer offers no edges

**Likely cause**
- the store has no `obsp['connectivities']`, and no `knn_key` named another
  graph

**Fix**
- pass `knn_key="<key>"` (`--knn <key>`) naming a symmetric `obsp` graph, or
  choose one under **Session › Data mapping**

### “Observation field ‘…’ has unsupported dtype”

**Likely causes**
- `adata.obs` carries a `datetime64`, timedelta, or period column

**Fix**
- convert the column to a string, or pass `obs_keys=[...]` naming the columns
  to show (`--obs-key` from the CLI, repeated once per column)

### “Gene ‘X’ not found in var”

**Likely causes**
- you’re searching by gene symbol but `var.index` contains Ensembl IDs (or vice versa)
- wrong `gene_id_column`, or a gene axis other than the one you expected
  (`expression_source="raw"` names genes from `raw.var`)
- the name was repeated in `var` and is shown as `NAME-1`

**Fix**
- set `gene_id_column=...` (or `--gene-id-column`) to match your gene IDs

### `Gene cannot be shown` — the gene was found and refused

**Likely cause**
- the column holds an infinity. `NaN` is drawn as the neutral grey meaning “no
  value here”, but `+Inf` or `-Inf` has no position on a colour scale, so the
  viewer refuses the gene or continuous `obs` column with a counted diagnosis
  naming the first offending cells. It is refused one column at a time, as you
  select it.

**Fix**
- repair the values in the object and reload, or run `prepare()` once — it
  scans the whole matrix and counts every non-finite value in a single failure.

The repair one-liners are in {doc}`15_troubleshooting_viewing`; `prepare()`'s
refusals are in {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.

For all other issues: {doc}`15_troubleshooting_viewing`.

## Next steps

- Preparing a Cellucid-ready AnnData (fastest reads): {doc}`../c_data_preparation_api/index`
- Server details + endpoints: {doc}`09_server_mode_advanced`
