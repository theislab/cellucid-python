"""
Cellucid Jupyter Integration

Provides seamless integration with Jupyter notebooks for visualizing
cellucid datasets directly in notebook cells.

Usage:
    from cellucid.jupyter import show_anndata

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

Communication (bidirectional):
    Python → Frontend:
        - viewer.highlight_cells()
        - viewer.set_color_by()
        - viewer.set_visibility()
        - viewer.clear_highlights()
        - viewer.reset_view()
        - viewer.set_data_mapping()
        - viewer.request_session_bundle()
        - viewer.send_message()

    Frontend → Python:
        - @viewer.on_selection
        - @viewer.on_hover
        - @viewer.on_click
        - @viewer.on_ready
        - @viewer.on_command_error
        - @viewer.on_data_mapping
        - @viewer.on_message

    ``show_anndata()`` returns a :class:`Widget` (an anywidget): both
    directions, and the data the viewer reads, travel through the notebook's
    own connection, so it works in JupyterLab, Jupyter Notebook 7, VS Code,
    Google Colab, marimo and Databricks with no port or proxy. The widget's
    messages reach Python between cells. :class:`AnnDataViewer` is the same
    API for a viewer in a browser tab of its own, served by this machine
    (its page talks back over a WebSocket at ``/_cellucid/ws``).

The module is a package. Its submodules are:

``_wire``
    Exact value checks for every message on the viewer channel.
``_hooks``
    :class:`HookRegistry` and :class:`ViewerState`.
``_viewer``
    :class:`AnnDataViewer` (the viewer API, and a browser-tab viewer) and
    :func:`show_anndata`, notebook-context detection, the viewer-channel
    bridge, and interpreter-exit cleanup.
``_widget``
    :class:`Widget`, the notebook widget, and the viewer channel over its
    comm. Imported when first used: it needs ``anywidget``.
"""

from typing import TYPE_CHECKING, Any

from ._hooks import HookCallback, HookRegistry, ViewerState
from ._viewer import AnnDataViewer, SessionBundleRequest, cleanup_all, show_anndata

if TYPE_CHECKING:
    from ._widget import Widget


def __getattr__(name: str) -> Any:
    # The widget imports anywidget (and ipywidgets): loaded on first use, so
    # the rest of the package starts without them.
    if name == "Widget":
        from ._widget import Widget

        return Widget
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "AnnDataViewer",
    "HookCallback",
    "HookRegistry",
    "SessionBundleRequest",
    "ViewerState",
    "Widget",
    "cleanup_all",
    "show_anndata",
]
