"""Utilities for exporting AnnData/array data for the 3D UMAP web viewer."""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("cellucid")
except _metadata.PackageNotFoundError:  # pragma: no cover - handled in dev installs
    __version__ = "0.0.0"

from .prepare_data import export_data_for_web

__all__ = ["export_data_for_web", "__version__"]
