"""
Cellucid - 1D, 2D, and 3D single-cell WebGL visualization toolkit.

Writes Cellucid-ready AnnData files, serves AnnData files and objects to the
Cellucid web viewer (which reads them itself), and shows the viewer in
notebooks as a widget.

Main exports:
    - prepare: Write a Cellucid-ready AnnData (.h5ad, .zarr or .zarr.zip)
    - restore_original_order: Put a Cellucid-ready file's cells back in the order given
    - serve_anndata / AnnDataServer: Serve an AnnData file, store or object
    - show_anndata / Widget: Show an AnnData file, store or object in a notebook
    - AnnDataViewer: The same viewer API, in a browser tab of its own
"""

from importlib import metadata as _metadata
from typing import TYPE_CHECKING, Any

__version__ = _metadata.version("cellucid")

if TYPE_CHECKING:
    # Type-checking-only re-exports of every name in ``__all__``.
    #
    # ``TYPE_CHECKING`` is ``False`` at runtime, so this block never executes and
    # the lazy ``__getattr__`` below still keeps ``import cellucid`` free of
    # numpy, pandas, scipy, and anndata. A static checker does read it, and
    # without it the module ``__getattr__`` erases the public API: every export
    # resolves to ``Any``, so ``prepare(nonexistent_argument=123)`` checks clean
    # even though the package ships ``py.typed``.
    from .anndata_server import AnnDataServer, serve_anndata
    from .anndata_session import (
        PullSummary,
        Selection,
        ViewerResult,
        apply_cellucid_session_to_anndata,
        pull,
        read_gmt,
        write_gmt,
    )
    from .jupyter import AnnDataViewer, Widget, show_anndata
    from .prepare_data import prepare, restore_original_order
    from .session_bundle import CellucidSessionBundle
    from .subsample import subsample_cells
    from .vector_fields import add_transition_drift_to_obsm, compute_transition_drift
    from .web_cache import clear_web_cache, get_web_cache_dir


# All imports are lazy to ensure fast CLI startup.
# Heavy dependencies (numpy, scipy, pandas) are only loaded when needed.
def __getattr__(name: str) -> Any:
    # Writing a Cellucid-ready AnnData (imports numpy, pandas, scipy)
    if name == "prepare":
        from .prepare_data import prepare

        return prepare
    if name == "restore_original_order":
        from .prepare_data import restore_original_order

        return restore_original_order

    # AnnData in a notebook
    if name == "show_anndata":
        from .jupyter import show_anndata

        return show_anndata
    if name == "Widget":
        from .jupyter._widget import Widget

        return Widget
    if name == "AnnDataViewer":
        from .jupyter import AnnDataViewer

        return AnnDataViewer

    # AnnData server
    if name == "serve_anndata":
        from .anndata_server import serve_anndata

        return serve_anndata
    if name == "AnnDataServer":
        from .anndata_server import AnnDataServer

        return AnnDataServer

    # Session bundles + AnnData mutation
    if name == "CellucidSessionBundle":
        from .session_bundle import CellucidSessionBundle

        return CellucidSessionBundle
    if name == "apply_cellucid_session_to_anndata":
        from .anndata_session import apply_cellucid_session_to_anndata

        return apply_cellucid_session_to_anndata

    # The round trip: results back from the viewer and its exports
    if name in {"pull", "write_gmt", "read_gmt", "PullSummary", "Selection", "ViewerResult"}:
        from . import anndata_session

        return getattr(anndata_session, name)
    if name == "clear_web_cache":
        from .web_cache import clear_web_cache

        return clear_web_cache
    if name == "get_web_cache_dir":
        from .web_cache import get_web_cache_dir

        return get_web_cache_dir

    # Subsampling: the cells a Cellucid subsample keeps
    if name == "subsample_cells":
        from .subsample import subsample_cells

        return subsample_cells

    # Vector field utilities (CellRank drift, etc.)
    if name == "compute_transition_drift":
        from .vector_fields import compute_transition_drift

        return compute_transition_drift
    if name == "add_transition_drift_to_obsm":
        from .vector_fields import add_transition_drift_to_obsm

        return add_transition_drift_to_obsm

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    # Writing a Cellucid-ready AnnData
    "prepare",
    "restore_original_order",
    # Viewing AnnData
    "show_anndata",
    "Widget",
    "AnnDataViewer",
    "serve_anndata",
    "AnnDataServer",
    # Sessions
    "CellucidSessionBundle",
    "apply_cellucid_session_to_anndata",
    # The round trip
    "pull",
    "write_gmt",
    "read_gmt",
    "PullSummary",
    "Selection",
    "ViewerResult",
    "clear_web_cache",
    "get_web_cache_dir",
    # Subsampling
    "subsample_cells",
    # Vector field utilities
    "compute_transition_drift",
    "add_transition_drift_to_obsm",
    # Version
    "__version__",
]
