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

adapter = AnnDataAdapter(adata)  # in-memory
# or: adapter = AnnDataAdapter.from_file("data.h5ad")

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
- `var/<gene>.values.f32.bin`

In AnnData mode, the adapter serves these as **virtual endpoints** computed from AnnData on demand.

### Lazy loading behavior (important for large datasets)

- `.h5ad` can be served in *backed* mode so gene expression columns are fetched on demand.
- `.zarr` is inherently chunked/lazy.
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
- Duplicate gene IDs can make gene lookup ambiguous; prefer stable, unique identifiers.
- If `adata.X` is CSR, the adapter may materialize a CSC copy for efficient column access (memory trade-off).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Gene expression lookup is very slow”
Fix:
- Prefer serving a backed `.h5ad` or `.zarr` over in-memory dense matrices.
- For repeated access, export with {func}`~cellucid.prepare` instead.

### Symptom: “No embeddings detected”
Fix:
- Ensure you have an embedding in `adata.obsm` with a supported key (e.g. `X_umap`, `X_umap_2d`, `X_umap_3d`).

---

## See also

- {doc}`server` for AnnData servers
- {doc}`export` for creating exported datasets
