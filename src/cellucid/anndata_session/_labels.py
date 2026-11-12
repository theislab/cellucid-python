"""
Planning a session's labels: renamed obs columns and renamed categories.

The ``field-overlays`` value at a step names fields and categories as the
dataset does (``registries/rename-registry.js`` ``toEntries``):

- ``renamedFields`` ``[{source, key, name}]``: an obs column (``source``
  ``"obs"``) is renamed to ``name``, keeping its place (a ``column-renamed``
  row); a gene's display name (``"var"``) is recorded, never applied to data;
- ``renamedCategories`` ``[{source, field, category, name}]``: in the obs
  column ``field``, the category ``category`` (a label, as JavaScript tells
  labels apart: ``1`` is not ``"1"``) reads ``name`` (a
  ``categories-renamed`` row, or ``categories-merged`` when ``name`` is a
  category the column has, or another rename of the column takes it too). A
  categorical column's categories are renamed in place; a text column has its
  values renamed;
- ``deletedFields`` ``[{source, key, purged}]``: recorded, never applied to
  data (a ``field-deleted`` row; the column stays).

Each becomes one row of the plan, applied (``write``), already applied
(``unchanged``: the object holds the new name, and an earlier apply recorded
the rename), or not applied with the reason (``not-applied``): a column or
category the object does not have, a new name the object already uses for a
column it did not get from Cellucid.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .._exact_integers import as_exact_text, holds_unsafe_integers
from ._primitives import (
    _category_identity,
    _require_bool,
    _require_category_primitive,
    _require_exact_keys,
    _require_field_key,
    _require_list,
    _require_source,
)
from ._records import _ColumnPlan
from ._schema import DELETED_FIELD_KEYS, RENAMED_CATEGORY_KEYS, RENAMED_FIELD_KEYS


def _check_entries(overlays: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """The renames and deletions of a field-overlays value, checked."""
    renamed_fields: list[dict[str, Any]] = []
    seen_fields: set[tuple[str, str]] = set()
    for position, raw in enumerate(_require_list(overlays["renamedFields"], label="field-overlays value renamedFields")):
        label = f"field-overlays value renamedFields[{position}]"
        entry = _require_exact_keys(raw, RENAMED_FIELD_KEYS, label=label)
        source = _require_source(entry["source"], label=f"{label}.source")
        key = _require_field_key(entry["key"], label=f"{label}.key")
        name = _require_field_key(entry["name"], label=f"{label}.name")
        if (source, key) in seen_fields:
            raise ValueError(f"This session's field-overlays value renames the {source} field {key!r} twice.")
        seen_fields.add((source, key))
        renamed_fields.append({"source": source, "key": key, "name": name})
    renamed_categories: list[dict[str, Any]] = []
    seen_categories: set[tuple[str, str, tuple[str, Any]]] = set()
    for position, raw in enumerate(
        _require_list(overlays["renamedCategories"], label="field-overlays value renamedCategories")
    ):
        label = f"field-overlays value renamedCategories[{position}]"
        entry = _require_exact_keys(raw, RENAMED_CATEGORY_KEYS, label=label)
        source = _require_source(entry["source"], label=f"{label}.source")
        field = _require_field_key(entry["field"], label=f"{label}.field")
        category = _require_category_primitive(entry["category"], label=f"{label}.category")
        new_label = _require_category_primitive(entry["name"], label=f"{label}.name")
        identity = (source, field, _category_identity(category))
        if identity in seen_categories:
            raise ValueError(
                f"This session's field-overlays value renames the category {category!r} of {field!r} twice."
            )
        seen_categories.add(identity)
        renamed_categories.append({"source": source, "field": field, "category": category, "name": new_label})
    deleted_fields: list[dict[str, Any]] = []
    for position, raw in enumerate(_require_list(overlays["deletedFields"], label="field-overlays value deletedFields")):
        label = f"field-overlays value deletedFields[{position}]"
        entry = _require_exact_keys(raw, DELETED_FIELD_KEYS, label=label)
        deleted_fields.append({
            "source": _require_source(entry["source"], label=f"{label}.source"),
            "key": _require_field_key(entry["key"], label=f"{label}.key"),
            "purged": _require_bool(entry["purged"], label=f"{label}.purged"),
        })
    return renamed_fields, renamed_categories, deleted_fields


#: Why a label is not applied, in the words Python and R both use.
GENE_NAME_NOTE = "A gene's display name is recorded, never applied to the object."
GENE_CATEGORY_NOTE = "A gene has no categories."
DELETION_NOTE = "A deletion is recorded, never applied to the object: the column stays."


def absent_column_note(column: str) -> str:
    return f'The object has no obs column "{column}".'


def modality_column_note(column: str) -> str:
    return f'"{column}" is a column of a MuData modality; labels apply to the global obs.'


def taken_column_note(column: str) -> str:
    return f'The object already has an obs column "{column}" that Cellucid did not write; it is never replaced.'


def numbers_note(column: str) -> str:
    return f'The obs column "{column}" holds numbers, not categories.'


def no_category_note(column: str, category: Any) -> str:
    return f'The obs column "{column}" has no category "{label_text(category)}".'


def label_text(value: Any) -> str:
    """A category label as text, as JavaScript writes it (``1``, ``1.5``, ``true``)."""
    if isinstance(value, bool | np.bool_):
        return "true" if value else "false"
    if isinstance(value, float | np.floating) and float(value).is_integer():
        return str(int(value))
    return str(value)


def _absent_column(column: str, obs: pd.DataFrame, modality_columns: set[str]) -> str | None:
    """Why a label of the obs column ``column`` cannot apply, or None when it can."""
    if column in modality_columns and column not in obs.columns:
        return modality_column_note(column)
    if column not in obs.columns:
        return absent_column_note(column)
    return None


def _renamed_categorical(column: pd.Series, renames: list[tuple[Any, Any]]) -> tuple[pd.Series, tuple[int, ...]]:
    """A categorical column with categories renamed (merged where two read
    alike, in the place of the first), and for each category after the place
    of the first category before it was made from (its colors follow it)."""
    categories = list(column.cat.categories)
    labels = list(categories)
    for category, name in renames:
        identity = _category_identity(category)
        for index, value in enumerate(categories):
            if _category_identity(value) == identity:
                labels[index] = name
    merged: list[Any] = []
    first: list[int] = []
    position: dict[tuple[str, Any], int] = {}
    remap = np.empty(len(labels) + 1, dtype=np.int64)
    remap[-1] = -1
    for index, label in enumerate(labels):
        identity = _category_identity(label)
        if identity not in position:
            position[identity] = len(merged)
            merged.append(label)
            first.append(index)
        remap[index] = position[identity]
    codes = remap[column.cat.codes.to_numpy().astype(np.int64)]
    values = pd.Categorical.from_codes(codes, categories=pd.Index(merged, dtype=object), ordered=bool(column.cat.ordered))
    return pd.Series(values, index=column.index), tuple(first)


def _renamed_text(column: pd.Series, renames: list[tuple[Any, Any]]) -> pd.Series:
    """A text column with the values each rename names replaced."""
    out = column.copy()
    values = column.to_numpy(dtype=object)
    for category, name in renames:
        identity = _category_identity(category)
        hits = np.fromiter((_category_identity(value) == identity for value in values), dtype=bool, count=values.size)
        # A column of pandas' text dtype holds text only.
        out[hits] = str(name) if isinstance(column.dtype, pd.StringDtype) else name
    return out


def _shown_column(column: pd.Series) -> pd.Series:
    """A column as the viewer reads its labels: a boolean one's values are
    its categories (``false``, ``true``), and an integer one holding
    identifiers past 2^53 is their exact decimal text (``_exact_integers``);
    any other column as it is. A label renamed in either is written as text."""
    if column.dtype.kind == "b" or isinstance(column.dtype, pd.BooleanDtype):
        return column.astype(object).where(column.notna(), None)
    if holds_unsafe_integers(column):
        return as_exact_text(column)
    return column


def _labels_of(column: pd.Series) -> list[Any] | None:
    """A column's categories (a categorical's, a text or boolean column's values,
    identifiers past 2^53 as text: ``_shown_column``), or None for numbers."""
    column = _shown_column(column)
    if isinstance(column.dtype, pd.CategoricalDtype):
        return list(column.cat.categories)
    if column.dtype.kind in "OSU" or pd.api.types.is_string_dtype(column.dtype):
        return list(pd.unique(column.dropna()))
    return None


def _count(column: pd.Series, category: Any) -> int:
    identity = _category_identity(category)
    values = _shown_column(column).to_numpy(dtype=object)
    return int(sum(1 for value in values if value is not None and not (isinstance(value, float) and np.isnan(value))
                   and _category_identity(value) == identity))


@dataclass(frozen=True)
class _Labels:
    """The labels of a step as they apply to one object's ``obs``: the plan's
    rows (in their order), the columns renamed (old name to new), the
    categorical or text columns rebuilt with renamed categories (by their
    name before any rename), and the renames an apply records."""

    rows: list[dict[str, Any]]
    rename: dict[str, str]
    replaced: list[_ColumnPlan]
    renamed: list[dict[str, Any]]


def _row(kind: str, name: str, status: str, **extra: Any) -> dict[str, Any]:
    return {"kind": kind, "domain": "field-overlays", "name": name, "column": None, "from": [], "cells": None,
            "status": status, "note": None, **extra}


def _plan_labels(
    value: dict[str, Any],
    obs: pd.DataFrame,
    *,
    modality_columns: set[str],
    earlier: list[dict[str, Any]],
) -> _Labels:
    """The renames and deletions of a field-overlays value as rows against
    ``obs``; ``earlier`` lists the renames earlier applies recorded."""
    renamed_fields, renamed_categories, deleted_fields = _check_entries(value)
    present = [str(column) for column in obs.columns]
    earlier_columns = {(entry["from"], entry["to"]) for entry in earlier if entry["kind"] == "column"}
    earlier_categories = {(entry["column"], entry["from"], entry["to"]) for entry in earlier if entry["kind"] == "category"}
    rows: list[dict[str, Any]] = []
    rename: dict[str, str] = {}
    renamed: list[dict[str, Any]] = []

    # Columns: every applicable rename at once (two swapped names apply).
    applicable = [entry for entry in renamed_fields
                  if entry["source"] == "obs" and entry["name"] != entry["key"] and entry["key"] in present]
    sources = {entry["key"] for entry in applicable}
    kept = set(present) - sources
    targets: set[str] = set()
    for entry in renamed_fields:
        key, name = entry["key"], entry["name"]
        if entry["source"] == "var":
            rows.append(_row("column-renamed", name, "not-applied", **{"from": [key], "note": GENE_NAME_NOTE}))
            continue
        if name == key:
            continue
        reason = _absent_column(key, obs, modality_columns)
        if reason is None:
            if name in kept or name in targets:
                rows.append(_row("column-renamed", name, "not-applied", **{"from": [key], "note": taken_column_note(name)}))
                continue
            targets.add(name)
            rename[key] = name
            renamed.append({"kind": "column", "column": None, "from": key, "to": name})
            rows.append(_row("column-renamed", name, "write", **{"from": [key]}))
        elif name in present and (key, name) in earlier_columns:
            rows.append(_row("column-renamed", name, "unchanged", **{"from": [key]}))
        else:
            rows.append(_row("column-renamed", name, "not-applied", **{"from": [key], "note": reason}))

    # Categories, in the column each is in once the columns are renamed.
    earlier_names = {old: new for old, new in earlier_columns}
    by_column: dict[str, list[tuple[Any, Any]]] = {}
    targets_of: dict[str, list[tuple[str, Any]]] = {}
    for entry in renamed_categories:
        if entry["source"] == "obs":
            targets_of.setdefault(entry["field"], []).append(_category_identity(entry["name"]))
    for entry in renamed_categories:
        field, category, name = entry["field"], entry["category"], entry["name"]
        if entry["source"] == "var":
            rows.append(_row("categories-renamed", label_text(name), "not-applied", **{
                "column": field, "from": [label_text(category)], "note": GENE_CATEGORY_NOTE,
            }))
            continue
        before = field if field in present else earlier_names.get(field, field)
        after = rename.get(before, before)
        reason = _absent_column(before, obs, modality_columns)
        if reason is not None:
            rows.append(_row("categories-renamed", label_text(name), "not-applied", **{
                "column": field, "from": [label_text(category)], "note": reason,
            }))
            continue
        column = obs[before]
        labels = _labels_of(column)
        if labels is None:
            rows.append(_row("categories-renamed", label_text(name), "not-applied", **{
                "column": after, "from": [label_text(category)], "note": numbers_note(after),
            }))
            continue
        identities = [_category_identity(item) for item in labels]
        has_from = _category_identity(category) in identities
        has_to = _category_identity(name) in identities
        shared = targets_of[field].count(_category_identity(name)) > 1
        # Into a category the column holds beside it, or the same as another
        # rename of the column: two categories become one. Once applied, the
        # source is gone and only the shared target tells.
        merged = shared or (has_from and has_to and _category_identity(category) != _category_identity(name))
        kind = "categories-merged" if merged else "categories-renamed"
        if has_from:
            by_column.setdefault(before, []).append((category, name))
            renamed.append({"kind": "category", "column": after, "from": label_text(category), "to": label_text(name)})
            rows.append(_row(kind, label_text(name), "write", **{
                "column": after, "from": [label_text(category)], "cells": _count(column, category),
            }))
        elif has_to and (after, label_text(category), label_text(name)) in earlier_categories:
            rows.append(_row(kind, label_text(name), "unchanged", **{"column": after, "from": [label_text(category)]}))
        else:
            rows.append(_row(kind, label_text(name), "not-applied", **{
                "column": after, "from": [label_text(category)], "note": no_category_note(after, category),
            }))
    replaced: list[_ColumnPlan] = []
    for column_name, renames in by_column.items():
        column = obs[column_name]
        if isinstance(column.dtype, pd.CategoricalDtype):
            values, first = _renamed_categorical(column, renames)
            replaced.append(_ColumnPlan(name=column_name, values=values, metadata={}, sources=first))
        else:
            values = _renamed_text(_shown_column(column), renames)
            if column.dtype.kind == "b" or isinstance(column.dtype, pd.BooleanDtype):
                # A boolean column with a label renamed is text: the others as the viewer names them.
                values = values.map(lambda value: label_text(value) if isinstance(value, bool | np.bool_) else value)
            replaced.append(_ColumnPlan(name=column_name, values=values, metadata={}))

    for entry in deleted_fields:
        rows.append(_row("field-deleted", entry["key"], "not-applied", note=DELETION_NOTE))
    return _Labels(rows=rows, rename=rename, replaced=replaced, renamed=renamed)
