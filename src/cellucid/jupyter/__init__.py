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

    # Control the viewer from Python: every action of the viewer is a method
    viewer.highlight_cells(cells=[100, 200, 300], color="#ff0000")
    viewer.color_by(field={"source": "obs", "key": "cell_type"})
    request = viewer.fields_list()  # a read: request.result() in a later cell

    # React to user interactions
    @viewer.on_selection
    def handle_selection(event):
        print(f"Selected {len(event['cells'])} cells")

Communication (bidirectional):
    Python → Frontend:
        - one method per action of the catalog (viewer.color_by(),
          viewer.highlight_cells(), viewer.history_get(), ...), each returning
          an ActionRequest; viewer.act(type, **params) for any of them
        - viewer.push() (a pandas column as a field)
        - viewer.request_session()
        - viewer.load_session() (a session file's trackers, added beside the person's)
        - with viewer.turn("QC"): (every action inside is one turn of the history)
        - viewer.history (LiveHistory: the history's actions, its latest change)
        - viewer.send_message()

    Frontend → Python:
        - ActionRequest.result(): an action's outcome, a read's result
          (ActionRequest.unrecorded, and an ActionNotRecordedWarning, say why
          a change made no history step)
        - @viewer.on_selection
        - @viewer.on_hover
        - @viewer.on_click
        - @viewer.on_ready
        - @viewer.on_command_error
        - @viewer.on_data_mapping
        - @viewer.on_history_change
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
``_actions``
    The action catalog (``cellucid/action_catalog.json``, the web app's
    document byte for byte) and the one conversion of every action's
    parameters; :data:`OMIT`.
``_actions_generated``
    :class:`ViewerActions`: one method per action, written by
    ``scripts/generate_actions.py`` from the catalog (never by hand).
``_hooks``
    :class:`HookRegistry`, :data:`HOOK_EVENTS` (the events a hook takes) and
    :class:`ViewerState`.
``_history``
    :class:`LiveHistory` (``viewer.history``) and :class:`LoadSessionRequest`.
``_viewer``
    :class:`AnnDataViewer` (the viewer API, and a browser-tab viewer),
    :class:`ActionRequest`, :class:`ActionError`,
    :class:`ActionNotRecordedWarning`, :class:`SessionRequest` and
    :func:`show_anndata`, notebook-context detection, the viewer-channel
    bridge, and interpreter-exit cleanup.
``_widget``
    :class:`Widget`, the notebook widget, and the viewer channel over its
    comm. Imported when first used: it needs ``anywidget``.
"""

from typing import TYPE_CHECKING, Any

from ._actions import OMIT, Omit, action_catalog
from ._actions_generated import ViewerActions
from ._history import LiveHistory, LoadSessionRequest
from ._hooks import HOOK_EVENTS, HookCallback, HookRegistry, ViewerState
from ._viewer import (
    ActionError,
    ActionNotRecordedWarning,
    ActionRequest,
    AnnDataViewer,
    SessionRequest,
    cleanup_all,
    show_anndata,
)

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
    "HOOK_EVENTS",
    "OMIT",
    "ActionError",
    "ActionNotRecordedWarning",
    "ActionRequest",
    "AnnDataViewer",
    "HookCallback",
    "HookRegistry",
    "LiveHistory",
    "LoadSessionRequest",
    "Omit",
    "SessionRequest",
    "ViewerActions",
    "ViewerState",
    "Widget",
    "action_catalog",
    "cleanup_all",
    "show_anndata",
]
