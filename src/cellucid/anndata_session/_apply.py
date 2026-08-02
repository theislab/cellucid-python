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

from ..session_bundle import CellucidSessionBundle
from ._chunks import _validate_current_chunk_inventory
from ._fields import _plan_user_defined_fields
from ._highlights import _plan_highlights
from ._primitives import (
    _require_exact_bool,
    _require_exact_keys,
    _require_nonempty_string,
    _require_nonnegative_int,
)
from ._records import ApplySummary, _ColumnPlan


def _validate_fingerprint(
    fingerprint: Any,
    adata: Any,
    *,
    expected_dataset_id: str,
) -> dict[str, Any]:
    # cellOrder ties the session to the cell ordering it was saved against. It
    # is accepted but not required, and never re-derived here: this applies a
    # session to an AnnData whose row order the caller controls, and the
    # viewer's digest is over exported coordinates an AnnData need not contain.
    # Enforcing that identity is the viewer's job; accepting both shapes is what
    # lets the viewer add the field without stranding existing session files.
    expected_keys = {"sourceType", "datasetId", "cellCount", "varCount"}
    fp = _require_exact_keys(
        fingerprint,
        expected_keys,
        label="session datasetFingerprint",
        optional_keys={"cellOrder"},
    )
    _require_nonempty_string(fp["sourceType"], label="datasetFingerprint.sourceType")
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
    n_vars = getattr(adata, "n_vars", None)
    if cell_count != n_obs:
        mismatches.append(f"cellCount {cell_count} != adata.n_obs {n_obs}")
    if var_count != n_vars:
        mismatches.append(f"varCount {var_count} != adata.n_vars {n_vars}")
    if mismatches:
        raise ValueError("Dataset fingerprint mismatch: " + "; ".join(mismatches))
    return copy.deepcopy(fp)


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
    """Apply a validated session atomically to the exact matching AnnData dataset."""
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

    fingerprint = _validate_fingerprint(
        bundle.dataset_fingerprint,
        adata,
        expected_dataset_id=expected_dataset_id,
    )
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
        n_obs=int(adata.n_obs),
        obs_index=adata.obs_names,
        prefix=highlights_prefix,
        existing_columns=existing_columns,
        materialize=add_highlights,
        metadata_by_id=inventory.metadata_by_id,
    )
    field_plans, overlays = _plan_user_defined_fields(
        bundle,
        chunk_ids,
        adata,
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
