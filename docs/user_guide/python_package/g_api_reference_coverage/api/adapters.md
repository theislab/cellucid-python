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
- `points_2d.bin`, `points_3d.bin`
- `var/<index>.values.f32` (payload files are named by index, never by gene name)

In AnnData mode, the adapter serves these as **virtual endpoints** computed from AnnData on demand.

### Loading behavior (important for large datasets)

- `.h5ad` paths are always opened read-only-backed so gene expression columns
  are fetched on demand.
- `.zarr` paths are loaded eagerly with `anndata.read_zarr`.
- In-memory AnnData uses whatever you already loaded into RAM.

### Which `obsm` keys become embeddings

`X_umap_1d`, `X_umap_2d`, and `X_umap_3d` each name their own dimension and
always decide it.

An object that declares none of those three but carries the bare `X_umap` that
`sc.tl.umap()` writes is read at the dimension its own column count states: a
one-, two-, or three-column array is resolved by its own width and served as
1D, 2D, or 3D. The adapter renames nothing and mutates nothing on your object.

- A bare `X_umap` of any other width — a ten-column latent space named after a
  plot, say — is refused, and the message names its shape.
- One explicit dimensional key anywhere in `obsm` means the writer declared
  dimensions deliberately, so the bare key never joins the resolved set.

The same resolution is used by {func}`~cellucid.serve_anndata`,
{func}`~cellucid.show_anndata`, {class}`~cellucid.AnnDataViewer`, and
`cellucid serve`, because all of them build this adapter.

### Connectivity is asked for, not assumed

`serve_connectivity` is keyword-only and `False` by default, on
{class}`~cellucid.AnnDataAdapter` and on `AnnDataAdapter.from_file` alike.

| `serve_connectivity` | What the adapter does |
|---|---|
| `False` (default) | `adata.obsp['connectivities']` is never read. `has_connectivity()` returns `False`, the dataset identity reports `has_connectivity: false` and `n_edges: null`, and there is no connectivity manifest to serve. |
| `True` | The complete graph is read and validated during construction, before the adapter can be used, and its manifest declares the edge count and the neighbor maximum. |

Off is the default because building the edge list is the single longest part of
startup on a large object — a 50-neighbor graph over millions of cells is
hundreds of millions of stored neighbors — and most sessions never draw the
graph. Validation is not deferred when you do ask for it: the manifest the
viewer fetches first has to state the edge count, so the whole graph is read up
front rather than on the first edge request.

Asking for connectivity on an object whose `obsp` holds no `connectivities`
matrix raises during construction. It is an explicit request, so an unmeetable
one is an error rather than a quietly graph-less dataset.

```python
from cellucid import AnnDataAdapter

adapter = AnnDataAdapter(
    adata,
    dataset_name="My study",
    dataset_id="my-study-v1",
    serve_connectivity=True,
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

- If your embedding keys are missing or have unexpected shapes, the adapter cannot serve `points_*d.bin`.
- Duplicate gene names, and names carrying characters with no glyph, are
  rejected during adapter construction. Names are not routes, so nothing else
  about them is constrained.
- If `adata.X` is CSR, the adapter may materialize a CSC copy for efficient column access (memory trade-off).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Gene expression lookup is very slow”
Fix:
- Prefer a `.h5ad` path for read-only-backed access.
- For repeated access, export with {func}`~cellucid.prepare` instead.

### Symptom: “No embeddings detected”
Fix:
- Ensure `adata.obsm` holds either an explicitly dimensioned embedding
  (`X_umap_1d`, `X_umap_2d`, or `X_umap_3d`) or a bare `X_umap` of 1, 2, or 3
  columns, which is read at the dimension its column count states.
- The error lists the `obsm` keys the object actually has, and when a bare
  `X_umap` is present but too wide it names that array's shape as well.

### Symptom: “Connectivity was asked for, and adata.obsp has no 'connectivities' matrix to serve”
Fix:
- Compute the neighbor graph first (`sc.pp.neighbors(adata)`), or
- build the adapter without `serve_connectivity=True`, which is the default and
  never touches `obsp` at all.

---

## See also

- {doc}`server` for AnnData servers
- {doc}`export` for creating exported datasets
