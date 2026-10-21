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

from .._contracts import (
    _WHITESPACE_RUN_PATTERN,
    CATEGORICAL_CODE_WIDTHS,
    JsonScalar,
    _display_text_defect,
)


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

    dtype: type[np.uint8] | type[np.uint16] | type[np.uint32]
    dtype_str: str
    missing_value: int
    codes_extension: str


def _storage(dtype: type[np.uint8] | type[np.uint16] | type[np.uint32]) -> _CategoricalCodesStorage:
    dtype_str = np.dtype(dtype).name
    width = CATEGORICAL_CODE_WIDTHS[dtype_str]
    return _CategoricalCodesStorage(dtype, dtype_str, width.missing, width.extension)


_CATEGORICAL_CODES_STORAGE = {
    storage.dtype_str: storage
    for storage in (_storage(np.uint8), _storage(np.uint16), _storage(np.uint32))
}
CATEGORICAL_CODE_DTYPES = tuple(_CATEGORICAL_CODES_STORAGE)


def _categorical_storage_for_count(n_categories: int) -> _CategoricalCodesStorage:
    """The narrowest codes storage that carries ``n_categories`` categories.

    Writers that choose the width themselves (``cellucid serve``,
    ``show_anndata``) use this, so uint32 appears only for a field with more
    than 65,535 categories and every narrower field keeps its narrower codes.
    """
    for storage in _CATEGORICAL_CODES_STORAGE.values():
        if n_categories <= storage.missing_value:
            return storage
    raise ValueError(
        f"A categorical field has {n_categories} categories; categorical codes "
        "carry at most 4,294,967,295 (uint32)."
    )


def _declared_categorical_storage(
    obs_categorical_dtype: str,
    n_categories: int,
    *,
    field_key: str,
) -> _CategoricalCodesStorage:
    """Resolve the codes storage a caller declared, and require the field to fit.

    ``obs_categorical_dtype`` is the caller's choice rather than something
    derived from the data, so the only thing left to settle is whether a
    field's categories fit the declared width. Each width spends its largest
    value on the missing marker: uint8 carries 255 categories, uint16 65,535
    and uint32 4,294,967,295.

    The exporter resolves the storage here twice for each field -- once in the
    pass that records what will be written, once in the loop that writes it --
    so a field cannot be summarised as one width and then encoded as another.
    """
    storage = _CATEGORICAL_CODES_STORAGE[obs_categorical_dtype]
    if n_categories > storage.missing_value:
        wider = _categorical_storage_for_count(n_categories).dtype_str
        raise ValueError(
            f"Field {field_key!r} has {n_categories:,} categories, but "
            f"{storage.dtype_str} supports at most {storage.missing_value:,}. "
            f"Pass obs_categorical_dtype={wider!r}."
        )
    return storage
