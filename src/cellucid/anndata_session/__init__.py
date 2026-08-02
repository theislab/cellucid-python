"""Apply one exact Cellucid session contract to :class:`anndata.AnnData`.

The module is a package. Its public surface is the entry point and the summary
it returns; everything else is split by responsibility across private
submodules:

``_schema``
    Every closed key set, chunk profile, and chunk-id prefix the contract
    declares, written once.
``_records``
    :class:`ApplySummary` and the internal plan records.
``_primitives``
    The exact value checks the rest of the package is built on.
``_chunks``
    Validation of the chunk listing a session declares.
``_highlights`` / ``_fields``
    The planners that turn highlight pages and user-defined fields into
    prepared ``obs`` columns.
``_apply``
    Fingerprint validation, ``uns`` recording, and the atomic application.
"""

from ._apply import apply_cellucid_session_to_anndata
from ._records import ApplySummary

__all__ = ["ApplySummary", "apply_cellucid_session_to_anndata"]
