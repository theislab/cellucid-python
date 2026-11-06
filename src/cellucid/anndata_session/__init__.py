"""Apply one step of a Cellucid session to :class:`anndata.AnnData`, and pull results.

The module is a package. Its public surface is the entry points and the
records they return; everything else is split by responsibility across
private submodules:

``_schema``
    Every closed key set of the state-domain values an application reads,
    written once.
``_records``
    The internal plan records (a column to write, a column a step holds).
``_primitives``
    The exact value checks the rest of the package is built on.
``_cells``
    Which cells a step's state knows (the cells the viewer had open).
``_labels`` / ``_highlights`` / ``_fields``
    The planners that turn renames, highlight groups and user-defined fields
    into renamed and prepared ``obs`` columns.
``_mudata``
    A MuData as the viewer showed it.
``_apply``
    The dataset checks (fingerprint, cell order, the fields a step uses),
    the plan (:class:`ApplyPlan`, one row per change), the atomic apply
    (:class:`ApplyReport`) and its record in ``uns['cellucid']['history_applied']``.
``_columns``
    The ``obs`` columns Cellucid writes (one rule per kind) and the one writer
    sessions and pulls share.
``_results``
    The viewer's results (highlights, fields, scores, gene sets), the web
    app's exports (obs-column CSV, GMT, session files), and :func:`pull`, the
    one call that writes them onto an AnnData or MuData.
"""

from ._apply import ApplyPlan, ApplyReport, apply_cellucid_session_to_anndata, history_applied
from ._results import PullSummary, Selection, ViewerResult, pull, read_gmt, write_gmt

__all__ = [
    "ApplyPlan",
    "ApplyReport",
    "history_applied",
    "PullSummary",
    "Selection",
    "ViewerResult",
    "apply_cellucid_session_to_anndata",
    "pull",
    "read_gmt",
    "write_gmt",
]
