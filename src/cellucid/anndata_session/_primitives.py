"""
Exact value checks shared by every part of the session reader.

A session bundle is untrusted input, so these are exact rather than permissive:
booleans must be booleans and not truthy values, key sets must match exactly
rather than be contained, and identifiers must match the wire pattern. Each
takes a ``label`` so the failure names the field the caller was reading.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np

from ._schema import _CATEGORY_LABEL_MAX_LENGTH, _FIELD_KEY_MAX_LENGTH, _WIRE_ID_RE


def _require_exact_bool(value: Any, *, label: str) -> bool:
    if type(value) is not bool:
        raise TypeError(f"{label} must be exactly True or False")
    return value


def _require_nonempty_string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a non-empty string")
    return value


def _require_wire_id(value: Any, *, label: str) -> str:
    identifier = _require_nonempty_string(value, label=label)
    if _WIRE_ID_RE.fullmatch(identifier) is None:
        raise ValueError(
            f"{label} must use only ASCII letters, digits, '.', '_', or '-' "
            "and be at most 180 characters"
        )
    return identifier


def _require_nonnegative_int(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    if value < 0:
        raise ValueError(f"{label} must be non-negative")
    return int(value)


def _require_exact_keys(
    value: Any,
    expected: set[str],
    *,
    label: str,
    optional_keys: set[str] | None = None,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    allowed = expected | (optional_keys or set())
    if not expected <= set(value) or not set(value) <= allowed:
        missing = sorted(expected - set(value))
        unknown = sorted(set(value) - allowed)
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unknown:
            details.append("unknown " + ", ".join(unknown))
        raise ValueError(f"{label} has invalid fields ({'; '.join(details)})")
    return value


def _require_field_key(value: Any, *, label: str) -> str:
    key = _require_nonempty_string(value, label=label)
    if key.strip() != key:
        raise ValueError(f"{label} cannot have leading or trailing whitespace")
    if len(key) > _FIELD_KEY_MAX_LENGTH:
        raise ValueError(f"{label} exceeds {_FIELD_KEY_MAX_LENGTH} Unicode code points")
    if ":" in key:
        raise ValueError(f"{label} cannot contain ':'")
    return key


def _reserve_column(existing: set[str], name: str) -> None:
    if name in existing:
        raise ValueError(f"Column already exists: {name}")
    existing.add(name)


def _require_optional_string(value: Any, *, label: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string or null")
    return value


def _require_category_list(
    value: Any,
    *,
    label: str,
    nonempty: bool,
) -> list[str | bool | int | float]:
    if not isinstance(value, list):
        raise TypeError(f"{label} must be an array")
    if nonempty and not value:
        raise ValueError(f"{label} must be a non-empty array")
    output: list[str | bool | int | float] = []
    identities: set[tuple[str, Any]] = set()
    for index, item in enumerate(value):
        exact = _require_category_primitive(item, label=f"{label}[{index}]")
        if isinstance(exact, str) and len(exact) > _CATEGORY_LABEL_MAX_LENGTH:
            raise ValueError(
                f"{label}[{index}] exceeds {_CATEGORY_LABEL_MAX_LENGTH} Unicode code points"
            )
        identity: tuple[str, Any]
        if isinstance(exact, bool):
            identity = ("boolean", exact)
        elif isinstance(exact, str):
            identity = ("string", exact)
        else:
            identity = ("number", float(exact))
        if identity in identities:
            raise ValueError(f"{label} must contain unique exact values")
        identities.add(identity)
        output.append(exact)
    return output


def _require_finite_number(value: Any, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{label} must be a finite number")
    if not np.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


def _require_highlight_identity(value: Any, *, prefix: str, label: str) -> str:
    identity = _require_nonempty_string(value, label=label)
    if re.fullmatch(rf"{re.escape(prefix)}[1-9]\d*", identity) is None:
        raise ValueError(f"{label} must use the current {prefix}<positive integer> identity")
    return identity


def _require_category_primitive(value: Any, *, label: str) -> str | bool | int | float:
    if isinstance(value, str | bool):
        return value
    if isinstance(value, int | float) and not isinstance(value, bool):
        _require_finite_number(value, label=label)
        return value
    raise TypeError(f"{label} must be a string, finite number, or boolean")
