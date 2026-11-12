"""
The colors a step shows categorical fields in: the user-defined fields' and
the dataset's columns' recolored in the view in focus.

The ``colors`` value at a step (``app/state/domains/colors.js``) holds, for
each view (live first, then kept views), the fields whose colors differ from
their opening colors; a categorical field's entry lists each category
recolored (``colors``: ``{category, color}`` by label, then labels held for
cells not open yet). The ``views`` value names the view in focus (``focus``,
its position; 0 is the live view).

A user field's category is shown, in the view in focus, in the color the
person or an agent gave it there, else in the color it opens with: its
``palette`` (a field made from highlight pages takes its pages'), else
Cellucid's own category colors by position
(:func:`cellucid._category_colors.default_category_color`). An apply writes
those colors with the field's column (``uns['<column>_colors']``). A column
that came with the dataset and was recolored in the view in focus is written
with the colors it is shown in as well (:func:`plan_column_colors`, one
``column-colors`` row each). Every other color the step holds (a continuous
field's color map or range, a gene's, another view's) is how the data was
shown, not applied.
cellucid-r reads them by the same rule (``R/session-apply.R``
``.apply_focused_recolors``), and so does the web app's "Save AnnData with
this state…" (``app/history/exports/anndata-copy.js``).
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd

from .._category_colors import default_category_color
from .._plural import plural
from .._scanpy_colors import colors_key, text_colors
from ._columns import _next_colors
from ._labels import absent_column_note, label_text
from ._primitives import (
    _category_identity,
    _refuse,
    _require_category_primitive,
    _require_list,
    _require_record,
    _require_source,
)
from ._records import _ColumnPlan

#: A color as a session writes it.
_HEX_COLOR = re.compile(r"#[0-9a-f]{6}")

#: Recolored categories of one user field, by ``(source, field id)``: each
#: category's identity (``_category_identity``) to its ``#rrggbb``.
Recolors = dict[tuple[str, str], dict[tuple[str, Any], str]]

#: The colors row's note when a field's colors went with its column.
COLORS_NOTE = (
    "The colors of the user fields written and of the dataset columns recolored go with their columns "
    "(uns['<column>_colors']); every other color is how the data was shown: not applied to the object."
)


def _focused_entries(colors: Any, views: Any) -> list[tuple[str, str, str | None, dict[tuple[str, Any], str]]]:
    """Each field the view in focus recolors categories of, in the value's
    order: ``(source, key, user id or None, {category identity: #rrggbb})``
    (``key`` None for a user field)."""
    if colors is None:
        return []
    record = _require_record(colors, label="colors value")
    per_view = _require_list(record.get("views"), label="colors value views")
    focus = 0
    if views is not None:
        given = _require_record(views, label="views value").get("focus")
        if isinstance(given, bool) or not isinstance(given, int) or given < 0:
            raise _refuse("views value focus", "must be a non-negative whole number.")
        focus = given
    if focus >= len(per_view):
        raise _refuse("colors value", f"holds {len(per_view)} {plural(len(per_view), 'view')}, and the view in focus is view {focus}.")
    entries = _require_list(per_view[focus], label=f"colors value views[{focus}]")
    out: list[tuple[str, str, str | None, dict[tuple[str, Any], str]]] = []
    for position, raw in enumerate(entries):
        label = f"colors value views[{focus}][{position}]"
        entry = _require_record(raw, label=label)
        field = _require_record(entry.get("field"), label=f"{label}.field")
        source = _require_source(field.get("source"), label=f"{label}.field.source")
        user = field.get("user")
        key = field.get("key")
        if user is not None:
            if not isinstance(user, str) or not user:
                raise _refuse(f"{label}.field.user", "must be a non-empty text.")
        elif not isinstance(key, str) or not key:
            raise _refuse(f"{label}.field.key", "must be a non-empty text.")
        if "colors" not in entry:
            continue
        recolored: dict[tuple[str, Any], str] = {}
        for index, item in enumerate(_require_list(entry["colors"], label=f"{label}.colors")):
            pair = _require_record(item, label=f"{label}.colors[{index}]")
            if set(pair) != {"category", "color"}:
                raise _refuse(f"{label}.colors[{index}]", "must name category and color only.")
            category = _require_category_primitive(pair["category"], label=f"{label}.colors[{index}].category")
            color = pair["color"]
            if not isinstance(color, str) or not _HEX_COLOR.fullmatch(color):
                raise _refuse(f"{label}.colors[{index}].color", "must be one #rrggbb color.")
            recolored[_category_identity(category)] = color
        out.append((source, key if isinstance(key, str) and user is None else "", user, recolored))
    return out


def focused_recolors(colors: Any, views: Any) -> Recolors:
    """The categories of user fields recolored in the view in focus.

    ``colors`` and ``views`` are the step's values (None when the session
    holds none: no category recolored, the live view in focus).
    """
    return {(source, user): recolored for source, _key, user, recolored in _focused_entries(colors, views)
            if user is not None}


def focused_column_recolors(colors: Any, views: Any) -> dict[str, dict[tuple[str, Any], str]]:
    """The categories of the dataset's obs columns recolored in the view in
    focus, by the column's key as the dataset names it (a column renamed in
    the session keeps its key), in the value's order. A gene's colors (a
    color map, a range) are never categories: not here."""
    return {key: recolored for source, key, user, recolored in _focused_entries(colors, views)
            if user is None and source == "obs" and recolored}


def shown_category_colors(categories: list[Any], palette: tuple[str, ...] | None,
                          recolored: dict[tuple[str, Any], str] | None) -> tuple[str, ...]:
    """The colors a categorical user field's categories are shown in, one
    ``#rrggbb`` a category in its order: the view in focus's recolors, else
    its ``palette``, else Cellucid's own colors by position."""
    given = recolored or {}
    return tuple(
        given.get(_category_identity(category))
        or (palette[index] if palette is not None else default_category_color(index))
        for index, category in enumerate(categories)
    )


# ---------------------------------------------------------------------------
# A dataset column the session recolored

#: A color as a file may hold it and the viewer reads it (``utils/hex-color.js``
#: ``bytesOfHex``): ``#rgb``, ``#rrggbb`` or ``#rrggbbaa``, either case.
_FILE_COLOR = re.compile(r"#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6}|[0-9a-fA-F]{8})")

#: The plan row of a dataset column's colors.
COLUMN_COLORS_KIND = "column-colors"


def _drawn_color(text: str) -> str:
    """A file's color as the viewer draws it: ``#rrggbb`` in lowercase
    (``#rgb`` widened, the alpha of ``#rrggbbaa`` dropped)."""
    digits = "".join(digit * 2 for digit in text[1:]) if len(text) == 4 else text[1:7]
    return "#" + digits.lower()


def opening_palette(held: Any, count: int) -> list[str] | None:
    """The colors the ``count`` categories of a dataset column open with in the
    viewer from the colors Scanpy keeps for it (``held``,
    ``uns['<column>_colors']``): the first ``count``, as drawn, when there
    are that many and each is a hex color (``palettes.js``
    ``scanpyCategoryPalette``); None otherwise, and the column opens in
    Cellucid's own colors by position."""
    colors = text_colors(held)
    if colors is None or count == 0 or len(colors) < count:
        return None
    if not all(_FILE_COLOR.fullmatch(color) for color in colors[:count]):
        return None
    return [_drawn_color(color) for color in colors[:count]]


def modality_colors_note(column: str) -> str:
    return f'"{column}" is a column of a MuData modality; colors apply to the global obs.'


def not_categorical_colors_note(column: str) -> str:
    return f'The obs column "{column}" is not categorical; Scanpy keeps colors for a categorical column only.'


def missing_categories_note(column: str, labels: Sequence[str]) -> str:
    quoted = ", ".join(f'"{label}"' for label in labels)
    return (f"Not written: the {plural(len(labels), 'color')} of {quoted}, "
            f'which the obs column "{column}" does not have.')


@dataclass(frozen=True)
class ColumnColors:
    """The dataset columns a step recolored, against one object: the plan's
    rows (``column-colors``, domain ``colors``, in the colors value's order)
    and the colors to write, by the column's name after the apply's renames."""

    rows: list[dict[str, Any]]
    recolor: dict[str, tuple[str, ...]]


def plan_column_colors(
    recolored: Mapping[str, Mapping[tuple[str, Any], str]],
    target: Any,
    *,
    replaced: Sequence[_ColumnPlan],
    rename: Mapping[str, str],
    earlier_names: Mapping[str, str],
    modality_columns: set[str],
) -> ColumnColors:
    """Each dataset column the view in focus recolored categories of
    (:func:`focused_column_recolors`, by the key the dataset names it),
    written with the colors it is shown in: each category, in the column's
    order once the apply's labels are applied, takes the color the person or
    an agent gave it, else the color it opens with -- the colors Scanpy keeps
    for the column (:func:`opening_palette`; a category merged by the apply
    takes its first source's), else Cellucid's own colors by position.

    A row is ``write`` when those differ from the colors the column holds
    after the apply's other changes, ``unchanged`` when they are the same,
    ``not-applied`` (with the reason) for a column the object lacks, a
    MuData modality's column, a column that is not categorical, or one with
    none of the categories recolored."""
    obs = target.obs
    replaced_by_name = {plan.name: plan for plan in replaced}
    would_hold = _next_colors(target, [], [], list(replaced), dict(rename))
    rows: list[dict[str, Any]] = []
    recolor: dict[str, tuple[str, ...]] = {}

    def row(name: str, status: str, note: str | None = None) -> dict[str, Any]:
        return {"kind": COLUMN_COLORS_KIND, "domain": "colors", "name": name, "column": None, "from": [],
                "cells": None, "status": status, "note": note}

    for key, given in recolored.items():
        before = key if key in obs.columns else earlier_names.get(key, key)
        after = rename.get(before, before)
        if before not in obs.columns:
            note = modality_colors_note(key) if key in modality_columns else absent_column_note(key)
            rows.append(row(after, "not-applied", note))
            continue
        column = obs[before]
        if not isinstance(column.dtype, pd.CategoricalDtype):
            rows.append(row(after, "not-applied", not_categorical_colors_note(after)))
            continue
        relabeled = replaced_by_name.get(before)
        count = len(column.cat.categories)
        categories = list(relabeled.values.cat.categories) if relabeled is not None else list(column.cat.categories)
        sources = relabeled.sources if relabeled is not None and relabeled.sources is not None else tuple(range(count))
        palette = opening_palette(target.uns.get(colors_key(before)), count)
        identities = [_category_identity(category) for category in categories]
        shown = tuple(
            given.get(identity) or (palette[source] if palette is not None else default_category_color(source))
            for identity, source in zip(identities, sources, strict=True)
        )
        present = set(identities)
        missing = [label_text(identity[1]) for identity in given if identity not in present]
        if len(missing) == len(given):
            rows.append(row(after, "not-applied", missing_categories_note(after, missing)))
            continue
        held = colors_key(after)
        holds = would_hold[held] if held in would_hold else target.uns.get(held)
        if text_colors(holds) == list(shown):
            rows.append(row(after, "unchanged", missing_categories_note(after, missing) if missing else None))
            continue
        recolor[after] = shown
        rows.append(row(after, "write", missing_categories_note(after, missing) if missing else None))
    return ColumnColors(rows=rows, recolor=recolor)
