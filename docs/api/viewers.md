# Viewer Classes

```{eval-rst}
.. currentmodule:: cellucid
```

High-level classes for embedding the Cellucid viewer in Jupyter notebooks. These provide more control than the simple {func}`~cellucid.show_anndata` and {func}`~cellucid.show` functions.

---

## AnnDataViewer

Class-based interface for visualizing AnnData objects in Jupyter. Provides access to the underlying widget for customization.

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(adata)
viewer.show()

# Access widget properties
viewer.width = 1000
viewer.height = 800
```

```{eval-rst}
.. autoclass:: AnnDataViewer
   :members:
   :show-inheritance:
```

---

## CellucidViewer

Class-based interface for visualizing pre-exported data in Jupyter.

```python
from cellucid import CellucidViewer

viewer = CellucidViewer("./my_export")
viewer.show()
```

```{eval-rst}
.. autoclass:: CellucidViewer
   :members:
   :show-inheritance:
```

---

## See Also

- {func}`~cellucid.show_anndata` - Simple function interface for AnnData
- {func}`~cellucid.show` - Simple function interface for exported data
- {class}`~cellucid.AnnDataServer` - Serve via HTTP instead of embedding
