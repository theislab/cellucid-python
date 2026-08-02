"""
Cellucid Jupyter Integration

Provides seamless integration with Jupyter notebooks for visualizing
cellucid datasets directly in notebook cells.

Usage:
    from cellucid.jupyter import CellucidViewer, show, show_anndata

    # Quick visualization
    viewer = show_anndata(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )

    # Control the viewer from Python
    viewer.highlight_cells([100, 200, 300], color="#ff0000")
    viewer.set_color_by("cell_type")

    # React to user interactions
    @viewer.on_selection
    def handle_selection(event):
        print(f"Selected {len(event['cells'])} cells")

Communication (bidirectional, works in all environments):
    Python → Frontend:
        - viewer.highlight_cells()
        - viewer.set_color_by()
        - viewer.set_visibility()
        - viewer.clear_highlights()
        - viewer.reset_view()
        - viewer.send_message()

    Frontend → Python:
        - @viewer.on_selection
        - @viewer.on_hover
        - @viewer.on_click
        - @viewer.on_ready
        - @viewer.on_command_error
        - @viewer.on_message

    Both directions work in Jupyter, JupyterLab, Google Colab, and VSCode.
    Frontend → Python uses HTTP POST to the local data server.

The module is a package. Its submodules are:

``_wire``
    Exact value checks for every message crossing the Python/iframe boundary,
    standard library only.
``_hooks``
    :class:`HookRegistry` and :class:`ViewerState`.
``_base``
    :class:`BaseViewer`, notebook-context detection, the HTTP event-routing
    bridge, and interpreter-exit cleanup.
``_exported`` / ``_anndata``
    :class:`CellucidViewer` and :func:`show` for a prepared export directory;
    :class:`AnnDataViewer` and :func:`show_anndata` for an AnnData object.
"""

from ._anndata import AnnDataViewer, show_anndata
from ._base import BaseViewer, cleanup_all
from ._exported import CellucidViewer, show
from ._hooks import HookCallback, HookRegistry, ViewerState

__all__ = [
    "AnnDataViewer",
    "BaseViewer",
    "CellucidViewer",
    "HookCallback",
    "HookRegistry",
    "ViewerState",
    "cleanup_all",
    "show",
    "show_anndata",
]
