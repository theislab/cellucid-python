"""
Apply a Cellucid `.cellucid-session` bundle to an `anndata.AnnData`.

Session bundles are treated as untrusted input:
- bounds checks on indices/codes
- dataset mismatch policies
- opt-in destructive mutations
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

from .session_bundle import CellucidSessionBundle
from .session_codecs import decode_delta_uvarint, decode_user_defined_codes

logger = logging.getLogger("cellucid.anndata_session")

DatasetMismatchPolicy = Literal["error", "warn_skip", "skip"]
ColumnConflictPolicy = Literal["error", "overwrite", "suffix"]


def _ensure_cellucid_uns(uns: dict[str, Any]) -> dict[str, Any]:
    root = uns.setdefault("cellucid", {})
    if not isinstance(root, dict):
        raise TypeError("adata.uns['cellucid'] must be a dict if present")
    return root


def _safe_key(s: str) -> str:
    key = re.sub(r"[^0-9a-zA-Z_]+", "_", (s or "").strip())
    key = re.sub(r"_+", "_", key).strip("_")
    if not key:
        return "cellucid"
    if key[0].isdigit():
        return f"_{key}"
    return key


def _resolve_column_name(
    existing: set[str],
    name: str,
    policy: ColumnConflictPolicy,
) -> str:
    if name not in existing:
        return name

    if policy == "error":
        raise ValueError(f"Column already exists: {name}")
    if policy == "overwrite":
        return name
    if policy != "suffix":
        raise ValueError(f"Unknown column conflict policy: {policy}")

    i = 2
    while True:
        candidate = f"{name}__{i}"
        if candidate not in existing:
            return candidate
        i += 1


@dataclass(frozen=True)
class ApplySummary:
    added_obs_columns: list[str]
    added_var_columns: list[str]
    skipped_due_to_mismatch: bool
    mismatch_reasons: list[str]


def apply_cellucid_session_to_anndata(
    bundle: CellucidSessionBundle | str,
    adata: "Any",
    *,
    inplace: bool = False,
    dataset_mismatch: DatasetMismatchPolicy = "warn_skip",
    expected_dataset_id: str | None = None,
    add_highlights: bool = True,
    highlights_prefix: str = "cellucid_highlight__",
    add_user_defined_fields: bool = True,
    user_defined_prefix: str = "",
    include_deleted_user_defined_fields: bool = False,
    store_uns: bool = True,
    column_conflict: ColumnConflictPolicy = "suffix",
    return_summary: bool = False,
) -> "Any" | tuple["Any", ApplySummary]:
    """
    Apply a `.cellucid-session` bundle onto an AnnData object.

    By default returns the (possibly copied) `adata`.
    If `return_summary=True`, returns `(adata, summary)`.
    """
    try:
        import pandas as pd  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError("apply_cellucid_session_to_anndata requires pandas") from e

    if isinstance(bundle, str):
        bundle = CellucidSessionBundle(bundle)

    if not inplace:
        adata = adata.copy()

    mismatch_reasons: list[str] = []
    fp = bundle.dataset_fingerprint or {}

    fp_cells = fp.get("cellCount")
    fp_vars = fp.get("varCount")
    fp_dataset_id = fp.get("datasetId")

    if isinstance(fp_cells, int) and fp_cells != getattr(adata, "n_obs", None):
        mismatch_reasons.append(f"cellCount {fp_cells} != adata.n_obs {getattr(adata, 'n_obs', None)}")
    if isinstance(fp_vars, int) and fp_vars != getattr(adata, "n_vars", None):
        mismatch_reasons.append(f"varCount {fp_vars} != adata.n_vars {getattr(adata, 'n_vars', None)}")
    if expected_dataset_id is not None and fp_dataset_id is not None and fp_dataset_id != expected_dataset_id:
        mismatch_reasons.append(f"datasetId {fp_dataset_id!r} != expected_dataset_id {expected_dataset_id!r}")

    has_mismatch = len(mismatch_reasons) > 0
    skip_dataset_dependent = False

    if has_mismatch:
        if dataset_mismatch == "error":
            raise ValueError("Dataset fingerprint mismatch: " + "; ".join(mismatch_reasons))
        if dataset_mismatch == "warn_skip":
            logger.warning("Dataset fingerprint mismatch; skipping dataset-dependent chunks: %s", mismatch_reasons)
            skip_dataset_dependent = True
        elif dataset_mismatch == "skip":
            skip_dataset_dependent = True
        else:
            raise ValueError(f"Unknown dataset_mismatch policy: {dataset_mismatch}")

    added_obs: list[str] = []
    added_var: list[str] = []

    if store_uns:
        cellucid_uns = _ensure_cellucid_uns(adata.uns)
        session_uns = cellucid_uns.setdefault("session", {})
        if not isinstance(session_uns, dict):
            raise TypeError("adata.uns['cellucid']['session'] must be a dict if present")
        session_uns["manifest"] = bundle.manifest
        session_uns["dataset_fingerprint"] = fp
        session_uns["applied"] = {
            "dataset_mismatch_policy": dataset_mismatch,
            "expected_dataset_id": expected_dataset_id,
            "skip_dataset_dependent": skip_dataset_dependent,
        }

    chunk_ids = set(bundle.list_chunk_ids())

    # ---------------------------------------------------------------------
    # Highlights → adata.obs boolean columns
    # ---------------------------------------------------------------------
    if add_highlights and not skip_dataset_dependent and "highlights/meta" in chunk_ids:
        meta = bundle.decode_chunk("highlights/meta")
        if store_uns:
            _ensure_cellucid_uns(adata.uns).setdefault("session", {}).setdefault("chunks", {})["highlights/meta"] = meta

        pages = meta.get("pages") if isinstance(meta, dict) else None
        if isinstance(pages, list):
            existing_cols = set(adata.obs.columns)
            highlights_uns = None
            if store_uns:
                session = _ensure_cellucid_uns(adata.uns).setdefault("session", {})
                highlights_uns = session.setdefault("highlights", {})
                if not isinstance(highlights_uns, dict):
                    highlights_uns = None

            for page in pages:
                if not isinstance(page, dict):
                    continue
                page_id = page.get("id")
                page_name = page.get("name")
                groups = page.get("highlightedGroups") or []
                if not isinstance(groups, list):
                    continue

                for group in groups:
                    if not isinstance(group, dict):
                        continue
                    group_id = group.get("id")
                    if not isinstance(group_id, str) or not group_id:
                        continue

                    membership_chunk_id = f"highlights/cells/{group_id}"
                    if membership_chunk_id in chunk_ids:
                        raw = bundle.decode_chunk(membership_chunk_id)
                        indices = decode_delta_uvarint(
                            raw,
                            max_count=int(getattr(adata, "n_obs", 0)),
                            max_index=int(getattr(adata, "n_obs", 0)) - 1,
                        )
                    else:
                        indices = np.empty(0, dtype=np.uint32)

                    base_name = f"{highlights_prefix}{_safe_key(group_id)}"
                    col_name = _resolve_column_name(existing_cols, base_name, column_conflict)
                    existing_cols.add(col_name)

                    mask = np.zeros(int(getattr(adata, "n_obs", 0)), dtype=bool)
                    if indices.size > 0:
                        mask[indices] = True
                    adata.obs[col_name] = pd.Series(mask, index=adata.obs_names)
                    added_obs.append(col_name)

                    if isinstance(highlights_uns, dict):
                        highlights_uns.setdefault("groups", {})[group_id] = {
                            "obs_column": col_name,
                            "page_id": page_id,
                            "page_name": page_name,
                            "group": group,
                        }

    # ---------------------------------------------------------------------
    # User-defined categorical fields → adata.obs/adata.var
    # ---------------------------------------------------------------------
    if add_user_defined_fields and not skip_dataset_dependent and "core/field-overlays" in chunk_ids:
        overlays = bundle.decode_chunk("core/field-overlays")
        if store_uns:
            _ensure_cellucid_uns(adata.uns).setdefault("session", {}).setdefault("chunks", {})[
                "core/field-overlays"
            ] = overlays

        udf = overlays.get("userDefinedFields") if isinstance(overlays, dict) else None
        if isinstance(udf, list):
            existing_obs = set(adata.obs.columns)
            existing_var = set(getattr(adata, "var", pd.DataFrame()).columns)
            for field in udf:
                if not isinstance(field, dict):
                    continue
                if field.get("kind") != "category":
                    continue
                if field.get("isDeleted") is True and not include_deleted_user_defined_fields:
                    continue

                field_id = field.get("id")
                if not isinstance(field_id, str) or not field_id:
                    continue

                source = field.get("source")
                target = "var" if source == "var" else "obs"

                codes_chunk_id = f"user-defined/codes/{field_id}"
                if codes_chunk_id not in chunk_ids:
                    continue

                raw = bundle.decode_chunk(codes_chunk_id)
                codes = decode_user_defined_codes(raw)

                if target == "obs":
                    expected_len = int(getattr(adata, "n_obs", 0))
                    if codes.shape[0] != expected_len:
                        logger.warning(
                            "Skipping user-defined field %s: codes length %s != adata.n_obs %s",
                            field_id,
                            codes.shape[0],
                            expected_len,
                        )
                        continue
                else:
                    expected_len = int(getattr(adata, "n_vars", 0))
                    if codes.shape[0] != expected_len:
                        logger.warning(
                            "Skipping user-defined var field %s: codes length %s != adata.n_vars %s",
                            field_id,
                            codes.shape[0],
                            expected_len,
                        )
                        continue

                categories = field.get("categories") if isinstance(field.get("categories"), list) else []
                categories = [str(c) for c in categories]

                base_key = str(field.get("key") or field_id)
                col_base = f"{user_defined_prefix}{_safe_key(base_key)}"

                if target == "obs":
                    col_name = _resolve_column_name(existing_obs, col_base, column_conflict)
                    existing_obs.add(col_name)
                else:
                    col_name = _resolve_column_name(existing_var, col_base, column_conflict)
                    existing_var.add(col_name)

                # Pandas uses -1 for NA in categorical codes; sanitize out-of-range codes.
                codes_i32 = codes.astype(np.int32, copy=False)
                if categories:
                    invalid = (codes_i32 < 0) | (codes_i32 >= len(categories))
                    if invalid.any():
                        codes_i32 = codes_i32.copy()
                        codes_i32[invalid] = -1
                    cat = pd.Categorical.from_codes(codes_i32, categories=categories, ordered=False)
                else:
                    # No categories list; store raw codes as integers.
                    cat = pd.Series(codes_i32)

                if target == "obs":
                    adata.obs[col_name] = pd.Series(cat, index=adata.obs_names)
                    added_obs.append(col_name)
                else:
                    adata.var[col_name] = pd.Series(cat, index=adata.var_names)
                    added_var.append(col_name)

    summary = ApplySummary(
        added_obs_columns=added_obs,
        added_var_columns=added_var,
        skipped_due_to_mismatch=skip_dataset_dependent,
        mismatch_reasons=mismatch_reasons,
    )
    if return_summary:
        return adata, summary
    return adata
