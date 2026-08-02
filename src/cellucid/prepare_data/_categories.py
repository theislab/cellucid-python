"""
Category label identity and storage.

A label is drawn verbatim wherever the viewer names it, so this module decides
both what a label is allowed to be -- an exact JSON scalar, distinct after that
representation, and text with a glyph for every character -- and how many of
them a codes payload can carry.
"""

import math
from collections.abc import Iterable, Sequence
from numbers import Integral, Real
from typing import NamedTuple

import numpy as np

from .._contracts import _WHITESPACE_RUN_PATTERN, JsonScalar, _display_text_defect


def _require_displayable_category_labels(
    categories: Sequence[JsonScalar],
    *,
    field_key: str,
) -> None:
    """Require every category label to read on screen as the value it stores.

    A label is drawn verbatim in the legend, the field selector, and every
    exported figure, so a character with no glyph makes the stored value differ
    from the one the reader sees. Two failures follow from that, and both are
    rejected here rather than repaired:

    * a label carrying an invisible character -- ``'Liver '`` reads as
      ``'Liver'`` but never matches it;
    * two labels that a whitespace-collapsing renderer draws identically --
      ``'T cell'`` and ``'T  cell'`` are two legend rows with one appearance.

    Trimming instead of rejecting would rewrite an annotation the caller never
    asked to change, and would merge two distinct categories into one, moving
    cells between them without saying so.
    """
    defects: list[str] = []
    for label in categories:
        if not isinstance(label, str):
            continue
        if not label:
            defects.append("'' is empty")
            continue
        defect = _display_text_defect(label)
        if defect is not None:
            defects.append(f"{label!r} {defect}")
    if defects:
        shown = defects[:10]
        remainder = len(defects) - len(shown)
        listed = "\n".join(f"  - {defect}" for defect in shown)
        if remainder:
            listed += f"\n  - ... and {remainder} more"
        counted = "1 label" if len(defects) == 1 else f"{len(defects)} labels"
        raise ValueError(
            f"Categorical field {field_key!r} has {counted} the viewer "
            f"cannot show as written:\n{listed}\n"
            "Labels are drawn verbatim in the legend, the field selector, and "
            "exported figures, so these characters are invisible on screen while "
            "still making the label a different value from the one it looks like. "
            "Cellucid does not clean them for you: trimming a label rewrites your "
            "annotation and can merge two categories into one, moving cells "
            "between them. Clean the column and export again, for example:\n"
            f"    obs[{field_key!r}] = obs[{field_key!r}].astype('string').str.strip()\n"
            f"    obs[{field_key!r}] = obs[{field_key!r}].astype('category')"
        )

    collapsed_to_label: dict[str, str] = {}
    for label in categories:
        if not isinstance(label, str):
            continue
        collapsed = _WHITESPACE_RUN_PATTERN.sub(" ", label)
        existing = collapsed_to_label.get(collapsed)
        if existing is not None:
            raise ValueError(
                f"Categorical field {field_key!r} labels {existing!r} and {label!r} "
                "are stored as different categories but are drawn identically: a "
                "run of whitespace collapses to a single space in the legend, the "
                "field selector, and exported figures. Rename one of them so the "
                "two categories can be told apart, then export again."
            )
        collapsed_to_label[collapsed] = label


def _json_category_values(
    values: Iterable[object],
    *,
    field_key: str,
) -> list[JsonScalar]:
    """Preserve one exact JSON-scalar identity for each category label."""
    categories: list[str | bool | int | float] = []
    seen: set[tuple[str, object]] = set()
    for raw_value in values:
        value = raw_value.item() if isinstance(raw_value, np.generic) else raw_value
        token: tuple[str, object]
        if isinstance(value, bool):
            token = ("boolean", value)
        elif isinstance(value, str):
            token = ("string", value)
        elif isinstance(value, Integral):
            integer_value = int(value)
            if abs(integer_value) > 9_007_199_254_740_991:
                raise ValueError(
                    f"Categorical field {field_key!r} contains integer label "
                    f"{integer_value!r} outside JavaScript's exact integer range."
                )
            value = integer_value
            token = ("number", value)
        elif isinstance(value, Real) and math.isfinite(value):
            value = float(value)
            token = ("number", value)
        else:
            raise ValueError(f"Categorical field {field_key!r} labels must be finite JSON scalars.")
        if token in seen:
            raise ValueError(
                f"Categorical field {field_key!r} labels must be unique "
                "after exact JSON representation."
            )
        seen.add(token)
        categories.append(value)
    _require_displayable_category_labels(categories, field_key=field_key)
    return categories


class _CategoricalCodesStorage(NamedTuple):
    """How one categorical field's codes are stored on disk."""

    dtype: type[np.uint8] | type[np.uint16]
    dtype_str: str
    missing_value: int
    codes_extension: str


def _declared_categorical_storage(
    obs_categorical_dtype: str,
    n_categories: int,
    *,
    field_key: str,
) -> _CategoricalCodesStorage:
    """Resolve the codes storage a caller declared, and require the field to fit.

    ``obs_categorical_dtype`` is the caller's choice rather than something
    derived from the data, so the only thing left to settle is whether a
    field's categories fit the declared width. uint8 spends code 255 on the
    missing marker and so carries 255 categories (0-254); uint16 spends 65,535
    the same way and carries 65,535.

    The exporter resolves the storage here twice for each field -- once in the
    pass that records what will be written, once in the loop that writes it --
    so a field cannot be summarised as one width and then encoded as another.
    """
    if obs_categorical_dtype == "uint8":
        if n_categories > 255:
            raise ValueError(
                f"Field {field_key!r} has {n_categories} categories, "
                "but uint8 supports at most 255."
            )
        return _CategoricalCodesStorage(np.uint8, "uint8", 255, "u8")
    if n_categories > 65_535:
        raise ValueError(
            f"Field {field_key!r} has {n_categories} categories, "
            "but uint16 supports at most 65,535."
        )
    return _CategoricalCodesStorage(np.uint16, "uint16", 65535, "u16")
