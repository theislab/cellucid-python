"""
The public entry point and the checks it performs itself.

:func:`apply_cellucid_session_to_anndata` validates the bundle's dataset
fingerprint against the AnnData it was handed, asks the planners for every
column the session declares, and only then writes: each plan is materialized
after all of them are built, so a session that fails halfway leaves the AnnData
untouched. :func:`_build_uns` records what was applied under ``adata.uns``.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

from .._anndata_files import is_mudata
from .._plural import plural
from ..data_mapping import (
    MODALITY_KEY,
    mapped_var_count,
    require_mapping_matches_anndata,
    require_mapping_modality,
    resolve_gene_names,
)
from ..session_bundle import CellucidSessionBundle
from ..session_bundle import _validate_fingerprint as _validate_bundle_fingerprint
from ._chunks import _validate_current_chunk_inventory
from ._fields import _plan_user_defined_fields
from ._highlights import _plan_highlights
from ._mudata import MuDataSessionView
from ._primitives import (
    _require_exact_bool,
    _require_nonempty_string,
    _require_nonnegative_int,
)
from ._records import ApplySummary, _ColumnPlan


def _validate_fingerprint(
    fingerprint: Any,
    adata: Any,
    *,
    expected_dataset_id: str,
) -> tuple[dict[str, Any], Any]:
    """Verify the fingerprint against ``adata`` before anything is applied.

    Returns the fingerprint and what the session is checked and planned
    against: ``adata`` itself, or for a MuData the cells and the modality its
    data mapping names, as the viewer showed them (:class:`MuDataSessionView`).

    The order is the viewer's: identity and cell count first, then the data
    mapping -- every key it names must exist on ``adata`` -- and then the gene
    count of the axis that mapping reads. ``cellOrder`` is validated for shape
    only. Its digest is taken over the normalized coordinates the viewer drew,
    and this bridge applies a session to an AnnData whose row order the caller
    controls, so re-deriving it would compare against coordinates the caller
    never asked to be drawn; enforcing that identity is the viewer's job.
    """
    try:
        fp = _validate_bundle_fingerprint(fingerprint)
    except ValueError as error:
        raise ValueError(f"session datasetFingerprint is invalid: {error}") from error
    dataset_id = _require_nonempty_string(
        fp["datasetId"],
        label="datasetFingerprint.datasetId",
    )
    cell_count = _require_nonnegative_int(
        fp["cellCount"],
        label="datasetFingerprint.cellCount",
    )
    var_count = _require_nonnegative_int(
        fp["varCount"],
        label="datasetFingerprint.varCount",
    )

    mismatches: list[str] = []
    if dataset_id != expected_dataset_id:
        mismatches.append(
            f"datasetId {dataset_id!r} != expected_dataset_id {expected_dataset_id!r}"
        )
    n_obs = getattr(adata, "n_obs", None)
    if cell_count != n_obs:
        mismatches.append(f"cellCount {cell_count} != adata.n_obs {n_obs}")
    if mismatches:
        raise ValueError("Dataset fingerprint mismatch: " + "; ".join(mismatches))

    mapping = fp["dataMapping"]
    if is_mudata(adata):
        require_mapping_modality(mapping, [str(name) for name in adata.mod])
        adata = MuDataSessionView(adata, mapping[MODALITY_KEY])
        mapping = {key: value for key, value in mapping.items() if key != MODALITY_KEY}
    require_mapping_matches_anndata(adata, mapping)
    expected_var_count = mapped_var_count(adata, mapping)
    if var_count != expected_var_count:
        axis = (
            "adata.raw.n_vars"
            if mapping["expression"] == "raw"
            else "genes on a mapping without expression"
            if mapping["expression"] is None
            else "adata.n_vars"
        )
        raise ValueError(
            f"Dataset fingerprint mismatch: varCount {var_count} != {axis} {expected_var_count}"
        )
    return copy.deepcopy(fp), adata


def _require_session_fields(adata: Any, requires: Any, mapping: dict[str, Any]) -> None:
    """Refuse an AnnData that lacks a field the session uses, naming all of them.

    ``requires`` lists the obs columns and the genes -- named under the
    session's data mapping, as the viewer named them -- that a view shows or
    filters on or a highlight group is bound to.
    """
    if not isinstance(requires, dict) or set(requires) != {"obs", "var"}:
        raise ValueError("bundle.manifest.requires must have exactly obs and var")
    obs_keys = set(adata.obs.columns)
    missing_obs = [key for key in requires["obs"] if key not in obs_keys]
    missing_var: list[str] = []
    if requires["var"]:
        gene_names = set(resolve_gene_names(adata, mapping).names)
        missing_var = [key for key in requires["var"] if key not in gene_names]
    if not missing_obs and not missing_var:
        return
    parts: list[str] = []
    if missing_obs:
        parts.append(f"the obs {plural(len(missing_obs), 'column')} " + ", ".join(repr(key) for key in missing_obs))
    if missing_var:
        parts.append(f"the {plural(len(missing_var), 'gene')} " + ", ".join(repr(key) for key in missing_var))
    raise ValueError(
        "This session uses " + " and ".join(parts) + ", which this dataset does not have."
    )


def _build_uns(
    adata: Any,
    *,
    bundle: Any,
    fingerprint: dict[str, Any],
    expected_dataset_id: str,
    highlight_meta: Any | None,
    highlights_uns: dict[str, Any] | None,
    overlays: Any | None,
    plans: list[_ColumnPlan],
) -> dict[str, Any]:
    raw_uns = getattr(adata, "uns", None)
    if not isinstance(raw_uns, dict):
        raise TypeError("adata.uns must be a dictionary")
    next_uns = copy.deepcopy(raw_uns)
    cellucid = next_uns.get("cellucid")
    if cellucid is None:
        cellucid = {}
        next_uns["cellucid"] = cellucid
    if not isinstance(cellucid, dict):
        raise TypeError("adata.uns['cellucid'] must be a dictionary if present")
    session = cellucid.get("session")
    if session is not None and not isinstance(session, dict):
        raise TypeError("adata.uns['cellucid']['session'] must be a dictionary if present")
    next_session: dict[str, Any] = {
        "manifest": copy.deepcopy(bundle.manifest),
        "dataset_fingerprint": fingerprint,
        "applied": {
            "expected_dataset_id": expected_dataset_id,
            "contract": "exact",
        },
        "materialized_fields": {plan.name: copy.deepcopy(plan.metadata) for plan in plans},
        "chunks": {},
    }
    if highlight_meta is not None:
        next_session["chunks"]["highlights/meta"] = highlight_meta
    if overlays is not None:
        next_session["chunks"]["core/field-overlays"] = overlays
    if highlights_uns is not None:
        next_session["highlights"] = highlights_uns
    cellucid["session"] = next_session
    return next_uns


def apply_cellucid_session_to_anndata(
    bundle: CellucidSessionBundle | str | Path,
    adata: Any,
    *,
    expected_dataset_id: str,
    inplace: bool = False,
    add_highlights: bool = True,
    highlights_prefix: str = "cellucid_highlight__",
    add_user_defined_fields: bool = True,
    user_defined_prefix: str = "",
    include_deleted_user_defined_fields: bool = False,
    store_uns: bool = True,
    return_summary: bool = False,
) -> Any | tuple[Any, ApplySummary]:
    """Apply a validated session atomically to the exact matching AnnData dataset.

    Every check -- dataset identity, shape, the cell-order record's form (its
    digest is the viewer's to compare), the obs columns and
    genes the session uses, chunk inventory, bounds, and target-column
    conflicts -- runs before ``adata`` (or its copy) is changed, so a failed
    application leaves nothing half-written.

    Parameters
    ----------
    bundle : CellucidSessionBundle, str, or Path
        The session, or the path of a ``.cellucid-session`` file.
    adata : anndata.AnnData or mudata.MuData
        The dataset the session was saved on, in the same cell order. For a
        MuData the session is checked against what the viewer showed -- the
        cells with every modality placed among them, and the modality its
        data mapping names -- and the added columns go to the global
        ``mdata.obs``; ``mdata.uns`` records the session.
    expected_dataset_id : str
        Required. The dataset id the session must have been saved on.
    inplace : bool
        Change ``adata`` itself instead of a copy. A backed AnnData needs
        ``inplace=True``.
    add_highlights : bool
        Add one boolean ``obs`` column per highlight group.
    highlights_prefix : str
        Prefix of each highlight column, followed by the group id.
    add_user_defined_fields : bool
        Add the session's user-defined cell fields as ``obs`` columns.
    user_defined_prefix : str
        Prefix of each user-defined field column, followed by its key.
    include_deleted_user_defined_fields : bool
        Also add user-defined fields that were deleted in the session.
    store_uns : bool
        Record the manifest, dataset fingerprint, and added columns in
        ``adata.uns['cellucid']['session']``.
    return_summary : bool
        Also return an ``ApplySummary`` listing the added ``obs`` columns.

    Returns
    -------
    anndata.AnnData or tuple of (anndata.AnnData, ApplySummary)
        The changed AnnData (``adata`` itself when ``inplace=True``), with the
        summary when ``return_summary=True``.

    Raises
    ------
    TypeError
        An argument of the wrong type, or a malformed bundle.
    ValueError
        The session does not match ``adata`` or ``expected_dataset_id``, uses
        a field ``adata`` lacks, holds an unknown or invalid chunk, or would
        overwrite an existing ``obs`` column.
    """
    expected_dataset_id = _require_nonempty_string(
        expected_dataset_id,
        label="expected_dataset_id",
    )
    _require_exact_bool(inplace, label="inplace")
    _require_exact_bool(add_highlights, label="add_highlights")
    _require_exact_bool(add_user_defined_fields, label="add_user_defined_fields")
    _require_exact_bool(
        include_deleted_user_defined_fields,
        label="include_deleted_user_defined_fields",
    )
    _require_exact_bool(store_uns, label="store_uns")
    _require_exact_bool(return_summary, label="return_summary")
    if not isinstance(highlights_prefix, str):
        raise TypeError("highlights_prefix must be a string")
    if not isinstance(user_defined_prefix, str):
        raise TypeError("user_defined_prefix must be a string")

    if isinstance(bundle, str | Path):
        bundle = CellucidSessionBundle(bundle)
    for attribute in (
        "dataset_fingerprint",
        "manifest",
        "list_chunk_ids",
        "decode_chunk",
    ):
        if not hasattr(bundle, attribute):
            raise TypeError(f"bundle must expose {attribute!r}")
    if not hasattr(adata, "obs") or not hasattr(adata, "var") or not hasattr(adata, "uns"):
        raise TypeError("adata must be an AnnData-compatible object")
    is_backed = getattr(adata, "isbacked", None)
    if type(is_backed) is not bool:
        raise TypeError("adata.isbacked must be exactly True or False")
    if is_backed and not inplace:
        raise ValueError(
            "A backed AnnData target cannot be applied with inplace=False; "
            "pass an explicitly materialized AnnData object or choose inplace=True"
        )

    fingerprint, shown = _validate_fingerprint(
        bundle.dataset_fingerprint,
        adata,
        expected_dataset_id=expected_dataset_id,
    )
    _require_session_fields(shown, bundle.manifest["requires"], fingerprint["dataMapping"])
    raw_chunk_ids = bundle.list_chunk_ids()
    if not isinstance(raw_chunk_ids, list) or any(
        not isinstance(item, str) for item in raw_chunk_ids
    ):
        raise TypeError("bundle.list_chunk_ids() must return a list of strings")
    if len(raw_chunk_ids) != len(set(raw_chunk_ids)):
        raise ValueError("Session bundle contains duplicate chunk ids")
    inventory = _validate_current_chunk_inventory(bundle, raw_chunk_ids)
    chunk_ids = inventory.ids

    existing_columns = set(adata.obs.columns)
    highlight_plans, highlights_uns, highlight_meta = _plan_highlights(
        bundle,
        chunk_ids,
        n_obs=int(shown.n_obs),
        obs=shown.obs,
        prefix=highlights_prefix,
        existing_columns=existing_columns,
        materialize=add_highlights,
        metadata_by_id=inventory.metadata_by_id,
    )
    field_plans, overlays = _plan_user_defined_fields(
        bundle,
        chunk_ids,
        shown,
        mapping=fingerprint["dataMapping"],
        prefix=user_defined_prefix,
        existing_columns=existing_columns,
        materialize=add_user_defined_fields,
        include_deleted=include_deleted_user_defined_fields,
        metadata_by_id=inventory.metadata_by_id,
    )
    plans = [*highlight_plans, *field_plans]

    next_obs = adata.obs.copy(deep=True)
    for plan in plans:
        next_obs[plan.name] = plan.values
    next_uns = (
        _build_uns(
            adata,
            bundle=bundle,
            fingerprint=fingerprint,
            expected_dataset_id=expected_dataset_id,
            highlight_meta=highlight_meta,
            highlights_uns=highlights_uns,
            overlays=overlays,
            plans=plans,
        )
        if store_uns
        else None
    )

    target = adata if inplace else adata.copy()
    target.obs = next_obs
    if next_uns is not None:
        target.uns = next_uns

    summary = ApplySummary(added_obs_columns=[plan.name for plan in plans])
    if return_summary:
        return target, summary
    return target
