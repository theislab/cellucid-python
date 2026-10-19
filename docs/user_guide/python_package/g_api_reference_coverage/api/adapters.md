# Adapters (AnnData → Cellucid data model)

```{eval-rst}
.. currentmodule:: cellucid
```

Adapters are the “server-side glue” that let Cellucid serve AnnData directly (without exporting first).

Most users never need to instantiate an adapter manually:
- use {func}`~cellucid.show_anndata` (notebook) or {func}`~cellucid.serve_anndata` / `cellucid serve …` (browser tab)

If you’re debugging, extending, or integrating Cellucid into custom servers, {class}`~cellucid.AnnDataAdapter` is the primary public adapter.

---

## Fast path (for developers)

```python
from cellucid import AnnDataAdapter

adapter = AnnDataAdapter(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
# AnnDataAdapter.from_file uses the same required identity arguments.

identity = adapter.get_dataset_identity()
obs_manifest = adapter.get_obs_manifest()
print(identity.get("name"), len(obs_manifest.get("fields", [])))

adapter.close()
```

---

## Practical path (what an adapter does)

### It emulates the exported on-disk format

The web viewer expects files like:
- `dataset_identity.json`
- `obs_manifest.json`
- `embeddings/<index>_<dim>d.bin`, one per served embedding
- `var/<index>.values.f32` (payload files are named by index, never by gene name)

In AnnData mode, the adapter serves these as **virtual endpoints** computed from AnnData on demand.

### Loading behavior (important for large datasets)

- `.h5ad` paths are always opened read-only-backed so gene expression columns
  are fetched on demand.
- `.zarr` paths are loaded eagerly with `anndata.read_zarr`.
- In-memory AnnData uses whatever you already loaded into RAM.

### Which `obsm` keys become embeddings

The adapter serves one embedding per dimension slot. `embedding_keys` names the
array of a slot — `{"2d": "X_tsne"}`, or `{"2d": ("X_pca", [0, 1])}` with
0-based columns. Every other slot follows the default rule: `X_umap_1d`,
`X_umap_2d`, and `X_umap_3d` each fill their own slot, and an object that
declares none of them but carries the bare `X_umap` that `sc.tl.umap()` writes
has it fill the slot of its own column count (1, 2, or 3). When neither fills
any slot, the adapter makes a deterministic best guess from the other `obsm`
arrays and records which slots were guessed. The adapter renames nothing and
mutates nothing on your object.

- `adapter.embeddings` lists the served embeddings, each with its canonical id
  (`X_umap`, or `X_pca[0,1]` when only some columns are read), `obsm` key,
  columns, dimension, and payload index.
- When nothing can be plotted, construction raises
  `Cellucid found nothing to plot…`, listing every `obsm` entry and why it was
  not used.

### The data mapping it serves

`adapter.get_data_mapping()` returns the data mapping the adapter serves. Its
`embeddings` are `{"enabled": [...], "defaults": {"1d", "2d", "3d"}}`: the
served slots are the defaults, and every other usable `obsm` array 1 to 3
columns wide is enabled beside them at all of its columns (not vector fields,
the latent key or arrays the best guess refused; see
`cellucid.data_mapping.default_enabled_embeddings`), each as an `obsm` key and
0-based columns. Then come `knn`,
`geneNames` (`"index"`, or `"column:<name>"` for `gene_id_column`),
`expression` (`"X"`, `"raw"`, `"layer:<name>"`, or `None` without expression),
`latent`, and the served vector bindings. It is the shape a session fingerprint
records as `dataMapping`. `adapter.get_mapping_candidates()` returns the catalog
the server publishes at `/_cellucid/mapping/candidates.json`: that mapping as
`defaults`, plus every `obsm`, `obsp`, layer and `var` column with its shape,
dtype and usability, read from metadata only.

## Connectivity is asked for, not assumed

`knn_key` is keyword-only and `None` by default, on
{class}`~cellucid.AnnDataAdapter` and on `AnnDataAdapter.from_file` alike.

| `knn_key` | What the adapter does |
|---|---|
| `None` (default) | No `obsp` matrix is read. `has_connectivity()` returns `False`, the dataset identity reports `has_connectivity: false` and `n_edges: null`, and there is no connectivity manifest to serve. |
| an `obsp` key, such as `"connectivities"` | That complete matrix is read and validated during construction, before the adapter can be used, and its manifest declares the edge count and the neighbor maximum. A matrix that is not symmetric is refused. |

Off is the default because building the edge list is the single longest part of
startup on a large object — a 50-neighbor graph over millions of cells is
hundreds of millions of stored neighbors — and most sessions never draw the
graph. Validation is not deferred when you do ask for it: the manifest the
viewer fetches first has to state the edge count, so the whole graph is read up
front rather than on the first edge request.

Naming an `obsp` key the object does not hold raises during construction. It is
an explicit request, so an unmeetable one is an error rather than a quietly
graph-less dataset.

```python
from cellucid import AnnDataAdapter

adapter = AnnDataAdapter(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    knn_key="connectivities",
)
print(adapter.has_connectivity())
adapter.close()
```

### Build progress

The adapter is silent by default (`quiet=True`), so using it as a library
prints nothing. A server passes its own `quiet` through, because opening a
large object spends minutes inside the constructor: with reporting on, the
build names `Obs columns`, `Embeddings` (first `resolving obsm keys`, then the
exact key resolved per dimension), `Vector fields`, and `Connectivity` as each
happens. See {doc}`server` for the full startup report.

---

## API reference

```{eval-rst}
.. autoclass:: AnnDataAdapter
   :members:
   :show-inheritance:
```

---

## Edge cases (do not skip)

- If no `obsm` array can be plotted, construction fails before any route exists.
- Repeated gene names are made unique (`MATR3`, `MATR3-1`) and empty values in
  `gene_id_column` take the row's index name. Names carrying characters with no
  glyph are rejected during adapter construction. Names are not routes, so
  nothing else about them is constrained.
- If `adata.X` is CSR, the adapter may materialize a CSC copy for efficient column access (memory trade-off).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Gene expression lookup is very slow”
Fix:
- Prefer a `.h5ad` path for read-only-backed access.
- For repeated access, export with {func}`~cellucid.prepare` instead.

### Symptom: “Cellucid found nothing to plot”
Fix:
- The error lists every `obsm` entry and why it was not used. Repair the array
  it names, compute an embedding, or pass `embedding_keys={"2d": "<key>"}`
  naming a usable array.

### Symptom: “knn_key names obsp['connectivities'], which this object does not have.”
Fix:
- Compute the neighbor graph first (`sc.pp.neighbors(adata)`), or
- build the adapter without `knn_key="connectivities"`, which is the default and
  never touches `obsp` at all.

---

## See also

- {doc}`server` for AnnData servers
- {doc}`export` for creating exported datasets
