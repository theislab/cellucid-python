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
        - viewer.get_session_bundle()
        - viewer.send_message()

    Frontend → Python:
        - @viewer.on_selection
        - @viewer.on_hover
        - @viewer.on_click
        - @viewer.on_ready
        - @viewer.on_command_error
        - @viewer.on_data_mapping
        - @viewer.on_message

    Both directions travel on one WebSocket the viewer opens to the local
    data server (``/_cellucid/ws``), so they work wherever the page in the
    notebook output can reach that server and open a WebSocket to it:
    classic Jupyter and JupyterLab on the same machine, an SSH-forwarded
    port, or a proxy that passes WebSockets (jupyter-server-proxy, with
    ``client_server_url=`` and ``allowed_hosts=``). VS Code's notebook webview
    and Google Colab's proxy have not been tested.

The module is a package. Its submodules are:

``_wire``
    Exact value checks for every message on the viewer channel.
``_hooks``
    :class:`HookRegistry` and :class:`ViewerState`.
``_viewer``
    :class:`AnnDataViewer` and :func:`show_anndata` for an AnnData file, store
    or object, notebook-context detection, the viewer-channel bridge, and
    interpreter-exit cleanup.
"""

from ._hooks import HookCallback, HookRegistry, ViewerState
from ._viewer import AnnDataViewer, cleanup_all, show_anndata

__all__ = [
    "AnnDataViewer",
    "HookCallback",
    "HookRegistry",
    "ViewerState",
    "cleanup_all",
    "show_anndata",
]
