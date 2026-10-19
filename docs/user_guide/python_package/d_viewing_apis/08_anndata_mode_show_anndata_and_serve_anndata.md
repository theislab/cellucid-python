# AnnData mode (`show_anndata()` / `serve_anndata()`)

This page explains how Cellucid can visualize **AnnData directly**, without running `prepare()` first.

Use AnnData mode when:
- you’re iterating inside an analysis notebook,
- you want a quick “does this dataset look sane?” view,
- or you have a large `.h5ad` and want read-only-backed access without exporting.

If you want maximum performance + easiest sharing, prefer export mode:
{doc}`07_exported_directory_mode_show_and_serve`.

## At a glance

**Audience**
- Computational users working with Scanpy/AnnData
- Wet lab users who received a `.h5ad`/`.zarr` from a collaborator

**Prerequisites**
- `pip install cellucid`
- A `.h5ad` / `.zarr` / in-memory `AnnData`

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
# A .zarr path or in-memory AnnData uses the same required identity arguments.
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
   present, a bare `X_umap` — the key `sc.tl.umap()` writes — fills the slot of
   its own column count, so the ordinary Scanpy object is served exactly as
   Scanpy left it.
3. **A deterministic best guess, only when neither filled any slot.** Cellucid
   then ranks the other `obsm` arrays by name (UMAP, then t-SNE, force-directed
   layouts, PHATE, PaCMAP/TriMAP, spatial, diffusion maps, other names, and
   latent spaces such as PCA or scVI last) and fills each slot with the
   best-ranked array of exactly that width. If no array has 1 to 3 columns, the
   best-ranked wider one is shown in 2D at columns 0 and 1. The startup report
   marks a guessed slot, for example `Embeddings: 2D from obsm['X_tsne'] (guessed)`.

If nothing in `obsm` can be plotted, construction fails with a message that
lists what `obsm` holds and why each entry was not used:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (3,696) and at least one numeric column. obsm holds: X_df (table), spatial (3,696 × 2, NaN at cell 17).
```

Those are the embeddings each dimension shows first, its *defaults*. The data
mapping the server publishes (`adapter.get_data_mapping()`) also *enables*
every other usable `obsm` array 1 to 3 columns wide beside them; see
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
be served, and construction fails naming that column. Either convert it
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
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataAdapter`, and
`cellucid serve … --obs-key`. It also fixes the order fields appear in, and a
column left out is not served at all.

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
exist, is a table or a sparse matrix, does not have one row per cell, is too
narrow for its slot, repeats a column, or holds a non-finite value in a selected
column fails construction with that reason. A slot you name replaces that slot's
default; slots you do not name still follow the default rule.

### The default rule

| `obsm` key | Slot it fills |
| --- | --- |
| `X_umap_1d` | 1D, columns `[0]` |
| `X_umap_2d` | 2D, columns `[0, 1]` |
| `X_umap_3d` | 3D, columns `[0, 1, 2]` |
| `X_umap`, only when none of the three is present | the slot of its own width, when that is 1, 2, or 3 |

Any combination of the three dimensional keys may be present, and each fills its
own slot. If any of them exists, `X_umap` never joins them: an object holding
both `X_umap_3d` and `X_umap` serves 3D only. `sc.tl.umap(adata)` writes two
columns and `sc.tl.umap(adata, n_components=3)` writes three, so both are
readable as written. A default key that cannot be plotted — a non-finite
coordinate, or no spread at all — fails construction rather than being skipped,
because the key says the array was meant to be drawn.

When the default rule fills any slot, the result is exactly what earlier releases
served: the same arrays, the same coordinates, and no guessing.

### The best guess

When neither `embedding_keys` nor the default rule fills any slot, every usable
`obsm` array — dense, numeric, one row per cell, at least one column — is a
candidate, except `<field>_umap_<n>d` vector fields and the `latent_key` array.
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
does not write to the `.h5ad` or `.zarr` you served — an `.h5ad` is opened
read-only-backed, so it could not. The array that fills each slot is reported
during startup, for example `Embeddings: 2D from obsm['X_umap']`, and again in
step 3 under its embedding id (`2D X_umap`, or `2D X_pca[0,1]` when only some
columns are read), so the transcript always names what is drawn.

### When nothing can be plotted

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell (3,696) and at least one numeric column. obsm holds: X_df (table), spatial (3,696 × 2, NaN at cell 17).
```

Each entry says what the array is and, for a dense array, why it was not used:
the first non-finite cell, a row count that is not the cell count, a dtype that
is not numeric, or its role (`a vector field`, `the latent space`). `obsm holds:
nothing.` means `obsm` is empty. Compute an embedding, or name one with
`embedding_keys`.

### Where the same rule applies

This resolution is used by {func}`~cellucid.show_anndata`,
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataViewer`,
{class}`~cellucid.AnnDataAdapter` (including `AnnDataAdapter.from_file`), and
`cellucid serve` on a path it detects as direct AnnData input.
{func}`~cellucid.add_transition_drift_to_obsm` uses the default rule only,
because it computes drift in a UMAP. The web app's own `.h5ad` and `.zarr`
readers apply the default rule.

```{important}
`prepare()` is a different path and shares none of this. It never looks at
`obsm` keys: you pass every embedding you want exported, by id, in
`embeddings={...}`, with {class}`~cellucid.Embedding` naming the columns of a
wide array. See {doc}`../c_data_preparation_api/03_embeddings_and_coordinates`.
```

## Loading and memory behavior (critical for large datasets)

### `.h5ad` (always read-only-backed)

When you pass a `.h5ad` path:
- Cellucid opens it with `anndata.read_h5ad(..., backed="r")`.
- This is usually the best option for large datasets if you don’t want to export.

### `.zarr` (eager)

Cellucid calls `anndata.read_zarr`, which materializes the AnnData arrays in
memory. Use a `.h5ad` path when read-only-backed access is required.

### In-memory `AnnData`

In-memory datasets are convenient but can be risky for large matrices:
- they can trigger expensive conversions (sparse format changes),
- and “accidentally load everything” becomes more likely.

Rule of thumb: if your dataset is big enough that you worry about RAM, view from a file.

## Gene IDs, expression source, and repeated names

Cellucid needs a mapping from “what you type in the gene search box” → a column index in the matrix.

By default, expression is read from `adata.X` and genes are named by
`var.index`. `expression_source` (`--expression`) chooses another matrix:
`"raw"` reads `adata.raw.X` and names genes from `adata.raw.var`, and
`"layer:<name>"` reads `adata.layers[<name>]`, which shares `var`. A source the
object does not have fails construction.

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
anndata's `var_names_make_unique()` does it — the second `MATR3` is served as
`MATR3-1`. Startup reports what changed, for example
`Gene names: 1 gene has no 'gene_symbols' value and uses its var index name. 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'.`
The rule is the one both export writers apply; see
{ref}`python_package-make-unique-gene-names`.

```{warning}
Served gene names must still be drawable exactly as stored: a name carrying a
control or zero-width character, or edge whitespace, fails construction. Payload
routes are integer indices, so `HLA-DRB1/2` and `Gene A` are served fine, and
`Gene` and `gene` are two distinct genes.
```

## Two optional capabilities, both asked for

The neighbor graph and the vector-field overlay are read only when you ask for
them. Everything else — points, obs fields, gene expression, centroids — is
always served.

| Capability | Reads | Python | Terminal |
|---|---|---|---|
| Neighbor graph | the `adata.obsp` matrix you name | `knn_key="connectivities"` | `--knn connectivities` |
| Vector fields | `adata.obsm['<field>_umap_<n>d']` | `serve_vector_fields=True` | `--vector-fields` |

Both are off for the same reason: each is read and validated **in full** before
the server binds its socket, because the manifests the viewer fetches first
declare what they contain. On a large object that is the difference between a
server that opens in a second and one that takes minutes, and most sessions never
turn either overlay on. Neither choice affects anything else the viewer can do.

The startup report names which state you are in, for each of them separately —
served, present in the object but not asked for, or absent. See
{doc}`/user_guide/web_app/b_data_loading/04_server_tutorial` for a full
transcript.

`prepare()` is unchanged by any of this: it takes `connectivities=` and
`vector_fields=` arguments, so passing the data has always been the ask. See
{doc}`../c_data_preparation_api/07_connectivities_knn_graph` and
{doc}`../c_data_preparation_api/08_vector_fields_velocity_displacement`.

## Vector fields (velocity/drift overlays)

Cellucid can render per-cell displacement vectors as an animated overlay, when
you ask for it with `serve_vector_fields=True` or `--vector-fields`.

### Naming convention (UMAP basis)

Required dimension-suffixed keys in `adata.obsm`:

- `velocity_umap_2d` (shape `(n_cells, 2)`)
- `velocity_umap_3d` (shape `(n_cells, 3)`)
- `T_fwd_umap_2d`, `T_bwd_umap_2d` (CellRank-style drift)

An unsuffixed key such as `velocity_umap` is not discovered as a vector field.

Rules:
- a `<field>_umap_<n>d` key is bound to the `<n>`D slot only when that slot shows
  `X_umap_<n>d` or `X_umap`, because a UMAP vector drawn over another layout
  points the wrong way. A key whose slot shows something else is not served
- vectors must match `n_cells` and the dimension of the embedding they are bound to
- vectors are scaled by that embedding's own normalization factor, into the same
  coordinate space as its points
- when the AnnData object declares more than one field ID, pass the exact
  `vector_field_default` to `show_anndata(...)`, `serve_anndata(...)`, or
  `AnnDataViewer(...)`; construction raises `ValueError` if it is omitted or
  does not name an available field
- `vector_field_default` selects among *served* fields, so naming one without
  asking for vector fields is refused rather than ignored
- an object that declares no field in this grammar simply serves none: asking for
  vector fields asks for whatever `obsm` declares, and declaring nothing is an
  ordinary result rather than an error. A *malformed* declaration still fails

For example, an object containing both `velocity_umap_2d` and
`T_fwd_umap_2d` needs:

```python
viewer = show_anndata(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    vector_field_default="velocity_umap",
)
```

```{figure} ../../../_static/screenshots/vector_field_velocity/overlay-controls.png
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with its exported `velocity` field; the field, density, speed,
trail, size, opacity, palette, and LOD controls are visible. The sample is a
prepared export, so its field id is the one `prepare()` was given; served
straight from AnnData, a `velocity_umap_2d` key is the field `velocity_umap`.
```

## Connectivity (KNN graph)

On the direct-AnnData path the neighbor graph is **opt-in, and off by default**.
A `.h5ad` straight out of `sc.pp.neighbors()` serves without a graph unless you
name the `obsp` matrix to draw:

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
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataViewer`,
{class}`~cellucid.AnnDataServer`, {class}`~cellucid.AnnDataAdapter`, and
`AnnDataAdapter.from_file`. From the terminal it is `--knn KEY`, and like every
other AnnData-only flag it is rejected when the path is a prepared export, which
already declares whatever artifacts it holds. A graph written by
`sc.pp.neighbors(adata, key_added="pca15")` is `knn_key="pca15_connectivities"`.

### What “off” means

- No `obsp` matrix is read. Whether `connectivities` exists is checked, because
  the startup report says so, but that is a key lookup and never touches the
  matrix.
- `dataset_identity.json` reports `stats.has_connectivity` as `false` and
  `stats.n_edges` as `null`.
- `connectivity_manifest.json` answers `404 No connectivity data`.
- `connectivity/edges.src.bin`, `connectivity/edges.dst.bin`, and
  `connectivity/edges.weights.f64.bin` answer `404`.
- The viewer reads a dataset that has no graph, exactly as it reads one whose
  object never had a graph, and offers no edge drawing. Points, obs fields,
  gene expression, centroids, and vector fields are unaffected.

### Why off is the default

Building the edge list is the single longest part of starting a server on a
large object. A 50-neighbor graph over millions of cells is hundreds of
millions of stored neighbors, and validating symmetry, deduplicating the stored
neighbors into edges, and ordering them costs minutes and several times the
graph's own memory. Most sessions never draw the graph.

### What “on” costs, and when you pay it

Naming a `knn_key` reads and validates the **whole** matrix before the server
binds its socket. This is not deferrable: the first document the viewer fetches
is the manifest, and the manifest declares the exact edge count and the neighbor
maximum, both of which are properties of the finished edge list. So the cost
lands during startup — between `Loading AnnData` and the printed URL — rather
than the first time someone switches edges on.

For a sparse matrix of 50,000,000 or more stored neighbors, Cellucid logs a
warning naming that count and the cell count before the validation starts, so a
long wait is announced rather than merely observed.

### A matrix that is not there, or not a graph, is an error

Naming a matrix the object does not have fails at construction rather than
serving a dataset silently missing the thing you asked for:

```text
knn_key names obsp['connectivities'], which this object does not have. obsp holds: []. Compute the neighbour graph with sc.pp.neighbors(adata), or serve without a graph.
```

The matrix is checked when it is read, and never transformed. Scanpy's
`obsp['distances']` is a directed k-nearest-neighbour matrix, so it is refused
with that reason:

```text
obsp['distances'] is not symmetric — choose a connectivities matrix. The neighbour graph must be undirected: Connectivity matrix must be exactly symmetric in topology and weight.
```

Any other contract violation fails the same way, before the socket exists, with
the violation quoted after `obsp['<key>'] cannot be served as the neighbour graph:`.

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
Only the direct-AnnData path chooses a matrix at serve time. `prepare()` takes
one graph through its `connectivities=` argument, and the prepared-export server
serves whatever artifacts the export contains.
```

## What startup prints

The direct-AnnData server runs five numbered steps, from `cellucid serve` and
from {func}`~cellucid.serve_anndata` alike. On a large object nearly all of the
time is spent inside step 2, which reports the adapter build as it happens
instead of going quiet:

Startup prints five numbered steps; the transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

Reading the transcript:

- `Embeddings: 2D from obsm['X_umap']` in step 2 names the exact array each
  slot resolved to. A declared object prints, for example,
  `Embeddings: 2D from obsm['X_umap_2d'], 3D from obsm['X_umap_3d']`; an array
  read at some of its columns adds them (`2D from obsm['X_pca'] columns [3, 1]`),
  and a guessed slot ends in `(guessed)`.
- Building the manifests and the categorical centroids is step 4. It used to
  run silently after step 3 reported success, which made a long pause look like
  a stall; it is now its own reported step.
- The prepared-export server still runs three steps. Five steps means direct
  AnnData input, and so does the `/?anndata=true` viewer query.

### The three connectivity states

The `Connectivity` line in step 3 distinguishes a graph that is served from one
this object holds but was not asked to serve:

| Printed | Meaning |
| --- | --- |
| `yes (762,984 edges)` | the graph was asked for, validated, and is served |
| `not served (obsp['connectivities'] is present; pass knn_key='connectivities', or --knn connectivities, to draw it)` | the object has a graph and this server was not asked for it |
| `no` | `adata.obsp` has no `connectivities` at all |

With `--knn connectivities`, step 2 gains the two lines that account for the wait,
and step 3 reports the edge count:

```text
[2/5] Loading AnnData...
      Mode: read-only backed h5ad
      Obs columns: classifying 16
      Embeddings: resolving obsm keys
      Embeddings: 2D from obsm['X_umap']
      Vector fields: scanning obsm
      Vector fields: 0 declared
      Connectivity: reading obsp['connectivities']
      Connectivity: 762,984 edges, 30 neighbors at most
      ✓ File opened
```

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

## Troubleshooting (AnnData mode)

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
- the object has no `X_umap*` key, so a slot was guessed; the startup report
  marks it `(guessed)`

**Fix**
- name the array you want with `embedding_keys` or `--embedding-<n>d`

### “knn_key names obsp[…], which this object does not have”

**Likely causes**
- the neighbor graph was never computed, was computed on a different object, or
  was stored under another key (`sc.pp.neighbors(..., key_added=...)`)

**Fix**
- run `sc.pp.neighbors(adata)` before serving (and re-save the file if you serve
  a path), name the key the object holds, or serve without a graph

### “obsp[…] is not symmetric — choose a connectivities matrix”

**Likely cause**
- the named matrix is a directed matrix such as `distances`

**Fix**
- name the symmetric `connectivities` matrix instead

### The viewer offers no edges, and the object has a graph

**Likely cause**
- connectivity is off by default on the direct-AnnData path. Step 3 of startup
  says so: `Connectivity: not served (obsp['connectivities'] is present; …)`

**Fix**
- pass `knn_key="connectivities"`, or `--knn connectivities` from the terminal,
  and expect the whole graph to be validated before the URL is printed

### “obs field ‘…’ has unsupported dtype”

**Likely causes**
- `adata.obs` carries a `datetime64`, timedelta, or period column

**Fix**
- convert the column to a string, or pass `obs_keys=[...]` naming the columns
  to serve (`--obs-key` from the CLI, repeated once per column)

### “Gene ‘X’ not found in var”

**Likely causes**
- you’re searching by gene symbol but `var.index` contains Ensembl IDs (or vice versa)
- wrong `gene_id_column`, or a gene axis other than the one you expected
  (`expression_source="raw"` names genes from `raw.var`)
- the name was repeated in `var` and is served as `NAME-1`

**Fix**
- set `gene_id_column=...` (or `--gene-id-column`) to match your gene IDs

### `Gene cannot be shown` — the gene was found and refused (HTTP 422)

**Likely cause**
- the column is not entirely finite float32. This server publishes finite
  float32 only, because a colour scale has no position for an infinity, so a
  gene or continuous `obs` column holding `NaN`, `±Inf`, or a magnitude float32
  cannot hold is answered **422** with a counted JSON diagnosis rather than
  served. It is refused one column at a time, as you select it.

**Fix**
- repair the values in the object and reload, or run `prepare()` once — the
  export scans every gene and names all the affected ones in a single failure.

The status, the exact body, and the repair one-liners are in
{doc}`15_troubleshooting_viewing`; the export-side counterpart is in
{doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`.

For all other issues: {doc}`15_troubleshooting_viewing`.

## Next steps

- Export mode (best performance): {doc}`07_exported_directory_mode_show_and_serve`
- Server details + endpoints: {doc}`09_server_mode_advanced`
