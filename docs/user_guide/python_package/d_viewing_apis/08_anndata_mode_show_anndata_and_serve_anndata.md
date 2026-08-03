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

### 1) UMAP embedding in `adata.obsm`

Each of these keys names its own dimension, and always decides it:

- `X_umap_1d` (shape `(n_cells, 1)`)
- `X_umap_2d` (shape `(n_cells, 2)`)
- `X_umap_3d` (shape `(n_cells, 3)`)

**An object that declares none of them is still served when it carries a bare
`X_umap`** — the key `sc.tl.umap()` writes. That array is read at the dimension
its own column count states: two columns are a 2D embedding, three columns are
a 3D embedding, one column is a 1D embedding. The ordinary Scanpy object is
therefore served exactly as Scanpy left it, with no rename and no edit to your
object.

Declare the dimensional keys when you want more than one dimension available at
once, or when the coordinates you want drawn are not the ones `X_umap` holds:

```python
adata.obsm["X_umap_2d"] = adata.obsm["X_umap"]
```

Use `X_umap_3d` if you ran `sc.tl.umap(adata, n_components=3)` and want those
coordinates under a key of their own. When you serve a `.h5ad` or `.zarr` path
rather than an in-memory object, re-save the file after the assignment.

If nothing resolves, construction fails with
“No UMAP embedding Cellucid can read was found in adata.obsm.”. The message
lists the `obsm` keys the object does have; names the shape of `X_umap` when
its width is the reason nothing resolved; and prints an exact statement to run
for every other array that could be an embedding. The complete rules are in the
section on which `obsm` keys are read, below.

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

### The keys that name their own dimension

| `obsm` key | Required shape | Dimension served |
| --- | --- | --- |
| `X_umap_1d` | `(n_cells, 1)` | 1D |
| `X_umap_2d` | `(n_cells, 2)` | 2D |
| `X_umap_3d` | `(n_cells, 3)` | 3D |

Any combination of the three may be present, and each is served at the
dimension its key names. The highest dimension present is the one the viewer
opens in. A key whose array does not carry exactly the column count its name
declares is an error, not a reinterpretation: `X_umap_3d` holding
`(n_cells, 2)` fails construction instead of being served as 2D.

### The bare `X_umap`, resolved by its own width

When none of `X_umap_1d`, `X_umap_2d`, or `X_umap_3d` is present and `obsm`
carries `X_umap`, that one array is read at the dimension its column count
states:

| `adata.obsm["X_umap"].shape` | Read as | Result |
| --- | --- | --- |
| `(n_cells, 1)` | 1D | served as the 1D embedding |
| `(n_cells, 2)` | 2D | served as the 2D embedding |
| `(n_cells, 3)` | 3D | served as the 3D embedding |
| any other width | — | construction fails, naming the shape |

`sc.tl.umap(adata)` writes two columns and `sc.tl.umap(adata, n_components=3)`
writes three, so both are readable as written. Only one dimension comes out of
this path, because there is only one array: to offer 2D *and* 3D in the same
session, declare `X_umap_2d` and `X_umap_3d`.

### An explicit key wins, and the bare key then stays out

If any of `X_umap_1d`, `X_umap_2d`, or `X_umap_3d` exists anywhere in `obsm`,
the set of embeddings is exactly the dimensional keys present, and `X_umap`
never joins it. An object holding both `X_umap_3d` and `X_umap` serves 3D only.
That is deliberate: a writer who named a dimension named it for the whole
object, and a second array of unknown provenance is not promoted alongside it.

### Nothing in your object is renamed or written

Resolution is a read. Cellucid does not assign `X_umap_1d`, `X_umap_2d`, or
`X_umap_3d` into `adata.obsm`, does not rename `X_umap`, and does not write to
the `.h5ad` or `.zarr` you served — an `.h5ad` is opened read-only-backed, so it
could not. The key that resolved is reported during startup, for example
`Embeddings: 2D from obsm['X_umap']`, so the transcript always names the array
being drawn.

### When no width is readable

An `X_umap` of some other width is a latent space someone named after a plot,
and no rename makes it drawable, so the message names its shape:

```text
No UMAP embedding Cellucid can read was found in adata.obsm. Cellucid reads the exact keys 'X_umap_1d', 'X_umap_2d', and 'X_umap_3d', and reads a bare 'X_umap' as the dimension its own column count states. Available obsm keys: ['X_pca', 'X_umap'].
  adata.obsm['X_umap'] has shape (3696, 10), and Cellucid draws 1, 2, or 3 dimensions, so that array names a dimension no viewer renders. Assign the columns you mean to draw under the key for their own count.
```

Assign the columns you actually want drawn under the key for their own count:

```python
adata.obsm["X_umap_2d"] = umap_coordinates_2d
```

When some *other* `obsm` array has `n_cells` rows and 1, 2, or 3 columns, the
message adds a `Fix:` block containing one copy-pasteable statement per
candidate, each naming that array's shape.

### Where the same rule applies

This resolution is the direct-AnnData rule. It is used by
{func}`~cellucid.show_anndata`, {func}`~cellucid.serve_anndata`,
{class}`~cellucid.AnnDataViewer`, {class}`~cellucid.AnnDataAdapter` (including
`AnnDataAdapter.from_file`), `cellucid serve` on a path it detects as direct
AnnData input,
and {func}`~cellucid.add_transition_drift_to_obsm`, which picks the embedding
its drift is computed in. The web app's own `.h5ad` and `.zarr` readers apply
the same rule, so an object that serves also opens in the browser.

```{important}
`prepare()` is a different path and shares none of this. It never looks at
`obsm` keys: the coordinates arrive as arrays, in the `X_umap_1d=`,
`X_umap_2d=`, and `X_umap_3d=` arguments, and each argument names the dimension
of the array you hand it. Passing `adata.obsm["X_umap"]` as `X_umap_2d=` is the
whole of “resolution” on the export path. See
{doc}`../c_data_preparation_api/03_embeddings_and_coordinates`.
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

## Gene IDs (`gene_id_column`) and duplicates

Cellucid needs a mapping from “what you type in the gene search box” → a column index in the matrix.

By default, `gene_id_column=None` uses `var.index`. Any non-blank string names
that exact `var` column, so `gene_id_column="index"` selects a column literally
named `"index"`.

If your gene IDs are in a column (e.g. `"gene_symbols"`), pass:

```python
viewer = show_anndata(
    "data.h5ad",
    gene_id_column="gene_symbols",
    dataset_name="My study",
    dataset_id="my-study-v1",
)
```

```{warning}
AnnData construction rejects duplicate gene names, and names carrying
characters with no glyph (control, zero-width, or edge whitespace). Payload
routes are integer indices, so `HLA-DRB1/2` and `Gene A` are served fine, and
`Gene` and `gene` are two distinct genes. Choose a distinct name in `var.index`
or the selected `gene_id_column`.
```

## Two optional capabilities, both asked for

The neighbor graph and the vector-field overlay are read only when you ask for
them. Everything else — points, obs fields, gene expression, centroids — is
always served.

| Capability | Reads | Python | Terminal |
|---|---|---|---|
| Neighbor graph | `adata.obsp['connectivities']` | `serve_connectivity=True` | `--connectivity` |
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
- vectors must match `n_cells` and the current dimension (2D vs 3D)
- vectors are scaled into the same normalized coordinate space as the embedding points
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
:alt: The real Pancreatic endocrinogenesis sample in its 2D Planar view with the velocity_umap particle overlay and controls visible.
:width: 1440px

The real 3,696-cell Pancreas sample in Build 2026-07-27.1, rendered in 2D
Planar mode with `velocity_umap`; the field, density, speed, trail, size,
opacity, palette, and LOD controls are visible.
```

## Connectivity (KNN graph)

On the direct-AnnData path the neighbor graph is **opt-in, and off by default**.
A `.h5ad` straight out of `sc.pp.neighbors()` serves without a graph unless you
ask for one.

```python
viewer = show_anndata(
    "data.h5ad",
    dataset_name="My study",
    dataset_id="my-study-v1",
    serve_connectivity=True,
)
```

```bash
cellucid serve /path/to/data.h5ad \
  --dataset-name "My study" \
  --dataset-id my-study-v1 \
  --connectivity
```

`serve_connectivity` is accepted by {func}`~cellucid.show_anndata`,
{func}`~cellucid.serve_anndata`, {class}`~cellucid.AnnDataViewer`,
{class}`~cellucid.AnnDataServer`, {class}`~cellucid.AnnDataAdapter`, and
`AnnDataAdapter.from_file`. From the terminal it is `--connectivity`, and like
every other AnnData-only flag it is rejected when the path is a prepared
export, which already declares whatever artifacts it holds.

### What “off” means

- `adata.obsp["connectivities"]` is never read. Whether the key exists is
  checked, because the startup report says so, but that is a key lookup and
  never touches the matrix.
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
graph's own memory. Most sessions never draw the graph, and until this was
opt-in every one of them paid for it.

### What “on” costs, and when you pay it

Passing `serve_connectivity=True` reads and validates the **whole** graph
before the server binds its socket. This is not deferrable: the first document
the viewer fetches is the manifest, and the manifest declares the exact edge
count and the neighbor maximum, both of which are properties of the finished
edge list. So the cost lands during startup — between `Loading AnnData` and the
printed URL — rather than the first time someone switches edges on.

For a sparse matrix of 50,000,000 or more stored neighbors, Cellucid logs a
warning naming that count and the cell count before the validation starts, so a
long wait is announced rather than merely observed.

### Asking for a graph that is not there is an error

`serve_connectivity=True` (or `--connectivity`) on an object whose `obsp` has no
`connectivities` fails at construction rather than serving a dataset silently
missing the thing you asked for:

```text
Connectivity was asked for, and adata.obsp has no 'connectivities' matrix to serve. Compute the neighbor graph with sc.pp.neighbors(adata) before serving, or serve without asking for connectivity.
```

An invalid graph fails the same way, before the socket exists, with the exact
contract violation quoted.

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
Only the direct-AnnData path changed. `prepare()` was already opt-in through its
`connectivities=` argument, and the prepared-export server is unchanged: an
export either contains the connectivity artifacts or it does not, and there is
nothing to ask for at serve time.
```

## What startup prints

The direct-AnnData server runs five numbered steps, from `cellucid serve` and
from {func}`~cellucid.serve_anndata` alike. On a large object nearly all of the
time is spent inside step 2, which reports the adapter build as it happens
instead of going quiet:

Startup prints five numbered steps; the transcript is in {doc}`../../web_app/b_data_loading/04_server_tutorial`.

Reading the transcript:

- `Embeddings: 2D from obsm['X_umap']` in step 2 names the exact array each
  dimension resolved to, so the key that decided the embedding is never a
  guess. A declared object prints, for example,
  `Embeddings: 2D from obsm['X_umap_2d'], 3D from obsm['X_umap_3d']`.
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
| `not served (obsp['connectivities'] is present; pass serve_connectivity=True, or --connectivity, to draw it)` | the object has a graph and this server was not asked for it |
| `no` | `adata.obsp` has no `connectivities` at all |

With `--connectivity`, step 2 gains the two lines that account for the wait,
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

### “No UMAP embedding Cellucid can read was found in adata.obsm”

**Likely causes**
- you didn’t compute UMAP yet
- `obsm` carries `X_umap`, but with a width no viewer draws (a 10-column latent
  space named after a plot, not 1, 2, or 3 coordinates)
- keys are present but have wrong shape — `X_umap_3d` holding two columns
- `n_cells` mismatch between `adata.obsm[...]` and `adata.n_obs`

**How to confirm**
- print `adata.obsm.keys()` and the shapes of candidate arrays

**Fix**
- run the statement the error prints, e.g.
  `adata.obsm["X_umap_2d"] = adata.obsm["X_umap"]`
- when the message names a shape, assign the columns you actually want drawn
  under the key for their own count
- or compute UMAP (Scanpy) and store explicit keys (`X_umap_2d` / `X_umap_3d`)

### “Connectivity was asked for, and adata.obsp has no ‘connectivities’ matrix to serve”

**Likely causes**
- `serve_connectivity=True` or `--connectivity` was passed for an object whose
  neighbor graph was never computed, or was computed on a different object

**Fix**
- run `sc.pp.neighbors(adata)` before serving (and re-save the file if you serve
  a path), or serve without asking for connectivity

### The viewer offers no edges, and the object has a graph

**Likely cause**
- connectivity is off by default on the direct-AnnData path. Step 3 of startup
  says so: `Connectivity: not served (obsp['connectivities'] is present; …)`

**Fix**
- pass `serve_connectivity=True`, or `--connectivity` from the terminal, and
  expect the whole graph to be validated before the URL is printed

### “obs field ‘…’ has unsupported dtype”

**Likely causes**
- `adata.obs` carries a `datetime64`, timedelta, or period column

**Fix**
- convert the column to a string, or pass `obs_keys=[...]` naming the columns
  to serve (`--obs-key` from the CLI, repeated once per column)

### “Gene ‘X’ not found in var”

**Likely causes**
- you’re searching by gene symbol but `var.index` contains Ensembl IDs (or vice versa)
- wrong `gene_id_column`

**Fix**
- set `gene_id_column=...` to match your gene IDs

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
