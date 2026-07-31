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
- Ensure you have an explicitly dimensioned embedding in `adata.obsm`
  (`X_umap_1d`, `X_umap_2d`, or `X_umap_3d`).

---

## See also

- {doc}`server` for AnnData servers
- {doc}`export` for creating exported datasets
