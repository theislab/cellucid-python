# Server Functions

```{eval-rst}
.. currentmodule:: cellucid
```

Functions and classes for running a local HTTP visualization server.

---

## Functions

### serve_anndata

Serve an AnnData object via HTTP for viewing in a browser. Supports in-memory objects, `.h5ad` files, and `.zarr` stores.

```python
from cellucid import serve_anndata

# Starts server and opens browser
serve_anndata(adata)
serve_anndata("data.h5ad")
serve_anndata("data.zarr")
```

```{eval-rst}
.. autofunction:: serve_anndata
```

---

### serve

Serve a pre-exported Cellucid dataset via HTTP. Use this when you've already run {func}`~cellucid.prepare`.

```python
from cellucid import prepare, serve

prepare(adata, "./my_export")
serve("./my_export")
```

```{eval-rst}
.. autofunction:: serve
```

---

## Classes

### AnnDataServer

Server class for AnnData with fine-grained control over the server lifecycle.

```python
from cellucid import AnnDataServer

server = AnnDataServer(adata)
server.start()
# ... do other work ...
server.stop()
```

```{eval-rst}
.. autoclass:: AnnDataServer
   :members:
   :show-inheritance:
```

---

### CellucidServer

Server class for pre-exported data with fine-grained control.

```python
from cellucid import CellucidServer

server = CellucidServer("./my_export")
server.start()
# ... do other work ...
server.stop()
```

```{eval-rst}
.. autoclass:: CellucidServer
   :members:
   :show-inheritance:
```

---

## See Also

- {func}`~cellucid.show_anndata` - Display in Jupyter instead of browser
- {func}`~cellucid.prepare` - Export data for static hosting
