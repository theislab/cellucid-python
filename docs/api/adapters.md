# Data Adapters

```{eval-rst}
.. currentmodule:: cellucid
```

Adapter classes for reading and transforming AnnData into the Cellucid format. These are used internally by the server and viewer classes.

---

## AnnDataAdapter

Adapter for reading AnnData objects and providing data in the format expected by the Cellucid viewer. Handles lazy loading for `.h5ad` and `.zarr` files.

```python
from cellucid import AnnDataAdapter

# Create adapter from various sources
adapter = AnnDataAdapter(adata)           # In-memory
adapter = AnnDataAdapter("data.h5ad")     # HDF5 file
adapter = AnnDataAdapter("data.zarr")     # Zarr store

# Access data
positions = adapter.get_positions()
metadata = adapter.get_metadata()
expression = adapter.get_gene_expression("CD3D")
```

```{eval-rst}
.. autoclass:: AnnDataAdapter
   :members:
   :show-inheritance:
```

---

## See Also

- {func}`~cellucid.prepare` - Export AnnData to static files
- {class}`~cellucid.AnnDataServer` - Uses adapter internally for serving
