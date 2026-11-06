"""
Planning the ``obs`` columns that carry a session's highlight groups.

The ``highlights`` value at a step (``app/state/domains/highlights.js``) is
``{active, pages: [{id, name, color, groups: [{id, type, label, enabled,
cells, field?, category?, range?}]}]}``: each group's cells are original
positions, ascending. One boolean column is planned per group of every page,
named ``cellucid_highlight__<page>__<group>`` from the page's name and the
group's label (:func:`highlight_column`), with the group's page, label, tool
and (for a category or range group) the field it was made from recorded
beside it. A cell the session does not know (one the viewer did not open) is
``<NA>``.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from ._columns import membership_column
from ._primitives import (
    _require_bool,
    _require_category_primitive,
    _require_cell_list,
    _require_exact_keys,
    _require_field_key,
    _require_finite_number,
    _require_list,
    _require_record,
    _require_source,
    _require_text,
)
from ._records import _Candidate
from ._schema import (
    CATEGORY_GROUP_KEYS,
    FIELD_REF_KEYS,
    GROUP_KEYS,
    GROUP_TYPES,
    HIGHLIGHTS_KEYS,
    PAGE_KEYS,
    RANGE_GROUP_KEYS,
    RANGE_KEYS,
)

_COLOR_RE = re.compile(r"#[0-9A-Fa-f]{6}\Z")
_NOT_A_WORD = re.compile(r"[^A-Za-z0-9]+")
#: The start of every highlight group's column.
HIGHLIGHT_PREFIX = "cellucid_highlight__"


def column_part(text: str, identifier: str) -> str:
    """A name as part of a column name: each run of characters other than
    ASCII letters and digits one underscore, none at either end; ``identifier``
    (the record's id) when nothing is left."""
    part = _NOT_A_WORD.sub("_", text).strip("_")
    return part if part else identifier


def unique_column(name: str, taken: set[str]) -> str:
    """``name``, or ``name_2``, ``name_3``... -- the first not in ``taken`` --
    added to ``taken``: two groups whose names read alike get two columns."""
    candidate = name
    number = 2
    while candidate in taken:
        candidate = f"{name}_{number}"
        number += 1
    taken.add(candidate)
    return candidate


def highlight_column(page: str, page_id: str, label: str, group_id: str) -> str:
    """A highlight group's column: ``cellucid_highlight__<page>__<group>``."""
    return f"{HIGHLIGHT_PREFIX}{column_part(page, page_id)}__{column_part(label, group_id)}"


def _group_metadata(group: dict[str, Any], *, label: str) -> dict[str, Any]:
    """A group checked, and what a column of it records (no cells)."""
    group_type = group.get("type") if isinstance(group, dict) else None
    if group_type not in GROUP_TYPES:
        raise ValueError(f"This session's {label}.type is not a highlight tool this version reads.")
    keys = CATEGORY_GROUP_KEYS if group_type == "category" else RANGE_GROUP_KEYS if group_type == "range" else GROUP_KEYS
    _require_exact_keys(group, keys, label=label)
    metadata: dict[str, Any] = {
        "label": _require_text(group["label"], label=f"{label}.label"),
        "type": group_type,
        "enabled": _require_bool(group["enabled"], label=f"{label}.enabled"),
    }
    if group_type in ("category", "range"):
        field = _require_exact_keys(group["field"], FIELD_REF_KEYS, label=f"{label}.field")
        metadata["field"] = {
            "source": _require_source(field["source"], label=f"{label}.field.source"),
            "key": _require_field_key(field["key"], label=f"{label}.field.key"),
        }
    if group_type == "category":
        metadata["category"] = _require_category_primitive(group["category"], label=f"{label}.category")
    if group_type == "range":
        bounds = _require_exact_keys(group["range"], RANGE_KEYS, label=f"{label}.range")
        low = _require_finite_number(bounds["min"], label=f"{label}.range.min")
        high = _require_finite_number(bounds["max"], label=f"{label}.range.max")
        if low > high:
            raise ValueError(f"This session's {label}.range has its min above its max.")
        metadata["range"] = {"min": bounds["min"], "max": bounds["max"]}
    return metadata


def _plan_highlights(
    value: Any,
    *,
    n_obs: int,
    index: pd.Index,
    known: np.ndarray | None,
) -> list[_Candidate]:
    """One column per highlight group of every page, checked whole before any is built."""
    highlights = _require_exact_keys(value, HIGHLIGHTS_KEYS, label="highlights value")
    pages = _require_list(highlights["pages"], label="highlights value pages")
    page_ids: set[str] = set()
    group_ids: set[str] = set()
    checked: list[tuple[str, np.ndarray, dict[str, Any]]] = []
    for page_index, raw_page in enumerate(pages):
        page_label = f"highlights value pages[{page_index}]"
        page = _require_exact_keys(raw_page, PAGE_KEYS, label=page_label)
        page_id = _require_text(page["id"], label=f"{page_label}.id")
        if page_id in page_ids:
            raise ValueError(f"This session's highlights value holds the page {page_id!r} twice.")
        page_ids.add(page_id)
        page_name = _require_text(page["name"], label=f"{page_label}.name")
        color = _require_text(page["color"], label=f"{page_label}.color")
        if _COLOR_RE.match(color) is None:
            raise ValueError(f"This session's {page_label}.color must be a six-digit hex colour.")
        for group_index, raw_group in enumerate(_require_list(page["groups"], label=f"{page_label}.groups")):
            group_label = f"{page_label}.groups[{group_index}]"
            group = _require_record(raw_group, label=group_label)
            metadata = _group_metadata(group, label=group_label)
            group_id = _require_text(group["id"], label=f"{group_label}.id")
            if group_id in group_ids:
                raise ValueError(f"This session's highlights value holds the group {group_id!r} twice.")
            group_ids.add(group_id)
            cells = _require_cell_list(group["cells"], n_obs, label=f"{group_label}.cells")
            checked.append((group_id, cells, {
                "kind": "highlight", "group_id": group_id, "page_id": page_id, "page_name": page_name, **metadata,
            }))
    active = highlights["active"]
    if active is not None and active not in page_ids:
        raise ValueError(f"This session's highlights value shows the page {active!r}, which it does not hold.")
    taken: set[str] = set()
    candidates: list[_Candidate] = []
    for group_id, cells, metadata in checked:
        name = unique_column(highlight_column(metadata["page_name"], metadata["page_id"], metadata["label"], group_id), taken)
        values = membership_column(n_obs, cells, open_cells=known, index=index)
        candidates.append(_Candidate(
            kind="highlight", domain="highlights", name=name, values=values,
            cells=int(np.count_nonzero(values.fillna(False).to_numpy(dtype=bool))), note=None, metadata=metadata,
        ))
    return candidates
