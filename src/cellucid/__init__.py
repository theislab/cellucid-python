"""
Cellucid - 3D UMAP WebGL visualization toolkit.

Utilities for exporting AnnData/array data for the 3D UMAP web viewer,
serving data locally or remotely, and integrating with Jupyter notebooks.

Main exports:
    - export_data_for_web: Export data to cellucid format
    - serve: Serve a dataset directory via HTTP
    - CellucidServer: Server class for more control
    - CellucidViewer: Jupyter notebook integration
"""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("cellucid")
except _metadata.PackageNotFoundError:  # pragma: no cover - handled in dev installs
    __version__ = "0.0.0"

from .prepare_data import export_data_for_web
from .server import serve, CellucidServer

# Lazy import for jupyter module (has optional dependencies)
def __getattr__(name):
    if name == "CellucidViewer":
        from .jupyter import CellucidViewer
        return CellucidViewer
    if name == "show":
        from .jupyter import show
        return show
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "export_data_for_web",
    "serve",
    "CellucidServer",
    "CellucidViewer",
    "show",
    "__version__",
]
