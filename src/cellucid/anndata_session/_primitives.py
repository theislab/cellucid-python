"""
Exact checks on the session values an application reads.

A session file's values are checked against their hashes when the file is
read (:class:`cellucid.CellucidSession`), which proves they are what was
saved, not that they are what this version understands. These checks are
exact rather than permissive: booleans must be booleans, key sets must match
exactly, and a field name must be one the viewer accepts. Each takes a
``label`` so the refusal names the value and field it was reading, and each
refuses with ``ValueError``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Collection
from typing import Any

import numpy as np

from ..session import _Array


def _refuse(label: str, what: str) -> ValueError:
    return ValueError(f"This session's {label} {what}")


def _require_bool(value: Any, *, label: str) -> bool:
    if type(value) is not bool:
        raise _refuse(label, "must be true or false.")
    return value


def _require_text(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise _refuse(label, "must be a non-empty text.")
    return value


def _require_optional_text(value: Any, *, label: str) -> str | None:
    if value is None:
        return None
    return _require_text(value, label=label)


def _require_record(value: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise _refuse(label, "must be a record.")
    return value


def _require_list(value: Any, *, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise _refuse(label, "must be a list.")
    return value


def _require_exact_keys(value: Any, expected: Collection[str], *, label: str) -> dict[str, Any]:
    record = _require_record(value, label=label)
    missing = sorted(set(expected) - set(record))
    unknown = sorted(set(record) - set(expected))
    if missing or unknown:
        details = [*(["missing " + ", ".join(missing)] if missing else []),
                   *(["unknown " + ", ".join(unknown)] if unknown else [])]
        raise _refuse(label, f"has other fields than this version reads ({'; '.join(details)}).")
    return record


def _require_finite_number(value: Any, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
        raise _refuse(label, "must be a finite number.")
    return float(value)


def _require_source(value: Any, *, label: str) -> str:
    if value not in ("obs", "var"):
        raise _refuse(label, 'must be "obs" or "var".')
    return str(value)


#: What the viewer refuses in a field name (``StateValidator.validateFieldKey``
#: in ``app/utils/state-validator.js``): control characters and zero-width
#: characters. A ``:`` is allowed -- ATAC feature names use it.
_FIELD_KEY_CONTROL = re.compile("[\u0000-\u001f\u007f-\u009f]")
_FIELD_KEY_INVISIBLE = re.compile("[​⁠﻿]")


def _require_field_key(value: Any, *, label: str) -> str:
    """A field name as the viewer accepts one when a field is named or renamed."""
    key = _require_text(value, label=label)
    if key.strip() != key:
        raise _refuse(label, "cannot have leading or trailing whitespace.")
    for pattern, what in ((_FIELD_KEY_CONTROL, "control"), (_FIELD_KEY_INVISIBLE, "zero-width")):
        found = pattern.search(key)
        if found is not None:
            raise _refuse(label, f"cannot contain the {what} character U+{ord(found.group(0)):04X}.")
    return key


def _require_category_primitive(value: Any, *, label: str) -> str | bool | int | float:
    """A category label: a string, finite number, or boolean (``validateCategoryLabel``)."""
    if isinstance(value, str | bool):
        return value
    if isinstance(value, int | float) and math.isfinite(value):
        return value
    raise _refuse(label, "must be a string, finite number, or boolean.")


def _category_identity(value: Any) -> tuple[str, Any]:
    """A category label's identity as JavaScript tells labels apart: ``1``
    and ``"1"`` differ, ``1`` and ``1.0`` do not."""
    if isinstance(value, bool | np.bool_):
        return ("boolean", bool(value))
    if isinstance(value, str):
        return ("string", value)
    if isinstance(value, int | float | np.integer | np.floating):
        return ("number", float(value))
    return ("other", repr(value))


def _require_category_list(value: Any, *, label: str) -> list[str | bool | int | float]:
    items = _require_list(value, label=label)
    if not items:
        raise _refuse(label, "must list at least one category.")
    out: list[str | bool | int | float] = []
    seen: set[tuple[str, Any]] = set()
    for index, item in enumerate(items):
        exact = _require_category_primitive(item, label=f"{label}[{index}]")
        identity = _category_identity(exact)
        if identity in seen:
            raise _refuse(label, f"lists the category {exact!r} twice.")
        seen.add(identity)
        out.append(exact)
    return out


def _require_array(value: Any, types: Collection[str], *, label: str) -> _Array:
    """A typed array of one of ``types`` (JavaScript type names)."""
    if not isinstance(value, _Array) or value.type not in types:
        raise _refuse(label, f"must be a {' or '.join(sorted(types))}.")
    return value


def _require_cell_list(value: Any, n_obs: int, *, label: str) -> np.ndarray:
    """Cells by original position: an ascending ``Uint32Array``, each below ``n_obs``."""
    cells = _require_array(value, ("Uint32Array",), label=label).data
    if cells.size > 1 and not bool(np.all(cells[1:] > cells[:-1])):
        raise _refuse(label, "must list cells strictly ascending.")
    if cells.size and int(cells[-1]) >= n_obs:
        raise _refuse(label, f"names the cell {int(cells[-1])}, and the dataset has {n_obs}.")
    return cells
