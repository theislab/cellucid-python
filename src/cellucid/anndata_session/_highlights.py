"""
Planning the ``obs`` columns that carry a session's highlight pages.

One boolean column is planned per enabled highlight group, named from the page
and group labels, with the group's declared metadata recorded alongside it. The
cell indices come from the cell-set chunk the group names
(:func:`cellucid.session_codecs.decode_cell_set`), and every index is checked
against the AnnData's cell count before any column is built. A category group
may be stored relative to its category's cells: those are read from
``adata.obs[fieldKey] == categoryName``, and a column whose category no longer
holds the cells the session recorded is refused.
"""

from __future__ import annotations

import copy
import re
from typing import Any

import numpy as np
import pandas as pd

from ..session_codecs import cell_set_needs_reference, decode_cell_set
from ._primitives import (
    _require_category_primitive,
    _require_exact_bool,
    _require_exact_keys,
    _require_field_key,
    _require_finite_number,
    _require_highlight_identity,
    _require_nonempty_string,
    _require_nonnegative_int,
    _reserve_column,
)
from ._records import _ColumnPlan
from ._schema import (
    _HIGHLIGHT_CATEGORY_GROUP_KEYS,
    _HIGHLIGHT_CELLS_PREFIX,
    _HIGHLIGHT_GROUP_REQUIRED_KEYS,
    _HIGHLIGHT_GROUP_TYPES,
    _HIGHLIGHT_PAGE_KEYS,
    _HIGHLIGHT_RANGE_GROUP_KEYS,
    _HIGHLIGHT_ROOT_KEYS,
)


def _validate_highlight_group_metadata(group: dict[str, Any], *, label: str) -> None:
    group_type = group.get("type")
    if group_type not in _HIGHLIGHT_GROUP_TYPES:
        raise ValueError(f"{label}.type is unsupported")
    if group_type == "category":
        expected_keys = _HIGHLIGHT_CATEGORY_GROUP_KEYS
    elif group_type == "range":
        expected_keys = _HIGHLIGHT_RANGE_GROUP_KEYS
    else:
        expected_keys = _HIGHLIGHT_GROUP_REQUIRED_KEYS
    _require_exact_keys(group, expected_keys, label=label)
    _require_highlight_identity(group["id"], prefix="highlight_", label=f"{label}.id")
    _require_nonempty_string(group["label"], label=f"{label}.label")
    _require_exact_bool(group["enabled"], label=f"{label}.enabled")
    cell_count = _require_nonnegative_int(group["cellCount"], label=f"{label}.cellCount")
    if cell_count < 1:
        raise ValueError(f"{label}.cellCount must be positive")
    if group_type in {"category", "range"}:
        _require_field_key(group["fieldKey"], label=f"{label}.fieldKey")
        if group["fieldSource"] not in {"obs", "var"}:
            raise ValueError(f"{label}.fieldSource must be 'obs' or 'var'")
    if group_type == "category":
        _require_category_primitive(
            group["categoryName"],
            label=f"{label}.categoryName",
        )
    elif group_type == "range":
        range_min = _require_finite_number(group["rangeMin"], label=f"{label}.rangeMin")
        range_max = _require_finite_number(group["rangeMax"], label=f"{label}.rangeMax")
        if range_min > range_max:
            raise ValueError(f"{label}.rangeMin must not exceed rangeMax")


def _category_reference(obs: pd.DataFrame, group: dict[str, Any], *, label: str) -> np.ndarray:
    """The cells of a category group's category, ascending.

    The viewer stores a category group relative to these cells when that is
    smaller; the decoder checks their count and digest against the ones it
    recorded, so a column whose category changed is refused, never misread.
    """
    if group["type"] != "category" or group["fieldSource"] != "obs":
        raise ValueError(
            f"{label} is stored relative to a category, and it is not a category "
            "group on an obs field"
        )
    field_key = group["fieldKey"]
    matching = [index for index, key in enumerate(obs.columns) if key == field_key]
    if len(matching) != 1:
        raise ValueError(
            f"{label} is stored relative to the category {group['categoryName']!r} of "
            f"obs field {field_key!r}, which must exist exactly once in adata.obs"
        )
    values = obs.iloc[:, matching[0]].astype(object).to_numpy()
    category = group["categoryName"]
    members = np.fromiter(
        (type(value) is type(category) and value == category for value in values),
        dtype=bool,
        count=values.size,
    ) if isinstance(category, bool) else values == category
    return np.flatnonzero(np.asarray(members, dtype=bool)).astype(np.uint32)


def _plan_highlights(
    bundle: Any,
    chunk_ids: set[str],
    *,
    n_obs: int,
    obs: pd.DataFrame,
    prefix: str,
    existing_columns: set[str],
    materialize: bool,
    metadata_by_id: dict[str, dict[str, Any]],
) -> tuple[list[_ColumnPlan], dict[str, Any] | None, Any | None]:
    membership_prefix = "highlights/cells/"
    membership_ids = {chunk_id for chunk_id in chunk_ids if chunk_id.startswith(membership_prefix)}
    if "highlights/meta" not in chunk_ids:
        if membership_ids:
            raise ValueError("Session contains highlight membership chunks without highlights/meta")
        return [], None, None

    raw_meta = bundle.decode_chunk("highlights/meta")
    meta = _require_exact_keys(
        raw_meta,
        _HIGHLIGHT_ROOT_KEYS,
        label="highlights/meta",
    )
    pages = meta["pages"]
    if not isinstance(pages, list):
        raise TypeError("highlights/meta.pages must be an array")
    if not pages:
        raise ValueError("highlights/meta.pages must be a non-empty array")
    active_page_id = _require_highlight_identity(
        meta["activePageId"],
        prefix="page_",
        label="highlights/meta.activePageId",
    )

    plans: list[_ColumnPlan] = []
    referenced_memberships: set[str] = set()
    page_ids: set[str] = set()
    group_ids: set[str] = set()
    ordered_group_ids: list[str] = []
    highlights_uns: dict[str, Any] = {"groups": {}}

    for page_index, raw_page in enumerate(pages):
        page = _require_exact_keys(
            raw_page,
            _HIGHLIGHT_PAGE_KEYS,
            label=f"highlights/meta.pages[{page_index}]",
        )
        page_id = _require_highlight_identity(
            page["id"],
            prefix="page_",
            label=f"highlights/meta.pages[{page_index}].id",
        )
        if page_id in page_ids:
            raise ValueError(f"Duplicate highlight page id: {page_id}")
        page_ids.add(page_id)
        page_name = _require_nonempty_string(
            page["name"],
            label=f"highlights/meta.pages[{page_index}].name",
        )
        page_color = _require_nonempty_string(
            page["color"],
            label=f"highlights/meta.pages[{page_index}].color",
        )
        if re.fullmatch(r"#[0-9A-Fa-f]{6}", page_color) is None:
            raise ValueError(
                f"highlights/meta.pages[{page_index}].color must be a six-digit hex color"
            )
        groups = page["highlightedGroups"]
        if not isinstance(groups, list):
            raise TypeError(
                f"highlights/meta.pages[{page_index}].highlightedGroups must be an array"
            )

        for group_index, raw_group in enumerate(groups):
            if not isinstance(raw_group, dict):
                raise TypeError(
                    f"highlights/meta.pages[{page_index}].highlightedGroups"
                    f"[{group_index}] must be an object"
                )
            group_label = f"highlights/meta.pages[{page_index}].highlightedGroups[{group_index}]"
            _validate_highlight_group_metadata(raw_group, label=group_label)
            group_id = raw_group["id"]
            if group_id in group_ids:
                raise ValueError(f"Duplicate highlight group id: {group_id}")
            group_ids.add(group_id)
            ordered_group_ids.append(group_id)
            membership_chunk_id = f"{membership_prefix}{group_id}"
            cell_count = raw_group["cellCount"]

            if membership_chunk_id not in chunk_ids:
                raise ValueError(
                    f"Highlight group {group_id!r} is missing required chunk "
                    f"{membership_chunk_id!r}"
                )
            membership_meta = metadata_by_id[membership_chunk_id]
            if membership_meta["label"] != f"Highlight cells: {raw_group['label']}":
                raise ValueError(
                    f"{membership_chunk_id} label must match highlight group {group_id!r}"
                )
            raw_membership = bundle.decode_chunk(membership_chunk_id)
            if not isinstance(raw_membership, bytes | bytearray | memoryview):
                raise TypeError(f"{membership_chunk_id} must decode to binary bytes")
            if membership_meta["uncompressedBytes"] != len(raw_membership):
                raise ValueError(f"{membership_chunk_id} uncompressedBytes must match its payload")
            reference = (
                _category_reference(obs, raw_group, label=group_label)
                if cell_set_needs_reference(raw_membership)
                else None
            )
            indices = decode_cell_set(
                raw_membership,
                cell_count=n_obs,
                reference=reference,
            )
            if len(indices) != cell_count:
                raise ValueError(
                    f"{membership_chunk_id} contains {len(indices)} indices, "
                    f"but metadata declares {cell_count}"
                )
            referenced_memberships.add(membership_chunk_id)

            column_name = f"{prefix}{group_id}"
            if materialize:
                _reserve_column(existing_columns, column_name)
                mask = np.zeros(n_obs, dtype=bool)
                mask[indices] = True
                plans.append(
                    _ColumnPlan(
                        name=column_name,
                        values=pd.Series(mask, index=obs.index, copy=False),
                        metadata={
                            "kind": "highlight",
                            "group_id": group_id,
                            "page_id": page_id,
                        },
                    )
                )
            highlights_uns["groups"][group_id] = {
                "obs_column": column_name if materialize else None,
                "page_id": page_id,
                "page_name": page_name,
                "group": copy.deepcopy(raw_group),
            }

    if active_page_id not in page_ids:
        raise ValueError(
            f"highlights/meta.activePageId {active_page_id!r} is not in the page inventory"
        )
    if membership_ids != referenced_memberships:
        unexpected = sorted(membership_ids - referenced_memberships)
        raise ValueError(
            "Session contains undeclared highlight membership chunks: " + ", ".join(unexpected)
        )
    manifest_membership_order = [
        chunk_id for chunk_id in metadata_by_id if chunk_id.startswith(_HIGHLIGHT_CELLS_PREFIX)
    ]
    expected_membership_order = [
        f"{_HIGHLIGHT_CELLS_PREFIX}{group_id}" for group_id in ordered_group_ids
    ]
    if manifest_membership_order != expected_membership_order:
        raise ValueError("Highlight membership chunks must match the exact metadata group order")
    return plans, highlights_uns, copy.deepcopy(meta)
