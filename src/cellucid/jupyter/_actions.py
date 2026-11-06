"""
The action catalogue, as the viewer API sends actions.

Everything a person can do in the viewer, and everything an agent or a
notebook reads from it, is an action of one catalogue, defined once in the web
app (``assets/js/app/actions/``). ``cellucid/action_catalog.json`` is that
catalogue's document byte for byte (pinned in ``tests/test_action_catalog.py``);
:class:`cellucid.jupyter.AnnDataViewer` has one generated method per action
(``_actions_generated.py``) and the generic :meth:`~AnnDataViewer.act`, both of
which check and convert parameters here, by the catalogue, the same way for
every action:

- a parameter the action does not have, or a required one missing, is refused
  at once, naming the parameters it has;
- an ``x-cells`` parameter (cells by original position) takes rows, a boolean
  mask or obs names, and travels as ``uint32``;
- an ``x-values`` parameter (one number per cell) travels as its value type;
- anything else becomes JSON (NumPy scalars and arrays, tuples and pandas
  indexes as lists), and the viewer checks it against the action's schema.

Values are checked further by the viewer, which answers with the action's
outcome or the one sentence that says why it did not apply it.
"""

from __future__ import annotations

import difflib
import json
from collections.abc import Callable, Mapping, Sequence
from functools import cache
from importlib import resources
from typing import Any, Final, TypeAlias

import numpy as np
import pandas as pd

from .._plural import plural
from ._wire import _MAX_CELL_INDEX

#: Cells an action names: 0-based rows of the viewer's ``obs`` (ints, a NumPy
#: integer array), a boolean mask of every cell, or obs names (a list of
#: names, a ``pandas.Index``); a boolean ``pandas.Series`` names the cells of
#: its index that are True.
CellsLike: TypeAlias = Sequence[int] | Sequence[bool] | Sequence[str] | np.ndarray | pd.Index | pd.Series

#: One number per cell: a NumPy array (``uint8``/``uint16``/``uint32`` codes or
#: floats) or a list of numbers.
ValuesLike: TypeAlias = Sequence[float] | np.ndarray | pd.Series


class Omit:
    """The default of an optional action parameter: leave it out, and the
    viewer applies the action's own default. Its one value is ``OMIT``."""

    __slots__ = ()
    _instance: Omit | None = None

    def __new__(cls) -> Omit:
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:
        return "OMIT"

    def __reduce__(self) -> tuple[type[Omit], tuple[()]]:
        return (Omit, ())


OMIT: Final = Omit()



@cache
def action_catalog() -> dict[str, dict[str, Any]]:
    """Every action of the catalogue this package was made with, by type."""
    text = resources.files("cellucid").joinpath("action_catalog.json").read_text(encoding="utf-8")
    document = json.loads(text)
    if document.get("schema") != 1:
        raise RuntimeError("cellucid/action_catalog.json is not an action catalogue this package reads")
    return {entry["type"]: entry for entry in document["actions"]}


def tool_name(action: str) -> str:
    """The method an action is called by (``filter.category.visibility`` ->
    ``filter_category_visibility``), as the web app's ``toolNameOf`` makes it."""
    return action.replace(".", "_").replace("-", "_")


def require_action(action: object) -> dict[str, Any]:
    """The catalogue's entry of ``action`` (its dotted type), or a refusal that
    names the nearest ones."""
    if type(action) is not str:
        raise TypeError(f"An action is named by its dotted type, such as 'color.by', not {action!r}")
    catalog = action_catalog()
    entry = catalog.get(action)
    if entry is not None:
        return entry
    near = difflib.get_close_matches(action, list(catalog), n=3)
    hint = f"; did you mean {', '.join(repr(name) for name in near)}?" if near else ""
    raise ValueError(
        f"There is no action {action!r}{hint} The viewer's actions are listed in "
        "cellucid.jupyter.action_catalog() (or as the viewer's methods)."
    )


def _annotated(schema: Mapping[str, Any], key: str) -> bool:
    return schema.get(key) is True or any(_annotated(sub, key) for sub in schema.get("oneOf", []))


def _fewest_items(schema: Mapping[str, Any], key: str) -> int:
    """The ``minItems`` of the array form of an ``x-cells``/``x-values`` schema."""
    if schema.get(key) is True:
        return int(schema.get("minItems", 0))
    return max((_fewest_items(sub, key) for sub in schema.get("oneOf", []) if _annotated(sub, key)), default=0)


def _json_param(value: Any, label: str) -> Any:
    """A parameter value as JSON: NumPy scalars and arrays, tuples and pandas
    indexes become their Python values and lists."""
    if value is None or type(value) in (bool, int, float, str):
        return value
    if isinstance(value, np.generic):
        return _json_param(value.item(), label)
    if isinstance(value, np.ndarray):
        return [_json_param(item, label) for item in value.tolist()]
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError(f"{label} keys must be text, not {key!r}")
            out[key] = _json_param(item, f"{label}.{key}")
        return out
    if isinstance(value, str | bytes):
        raise TypeError(f"{label} must be text, not bytes")
    tolist = getattr(value, "tolist", None)
    if callable(tolist) and not isinstance(value, Sequence):
        # A pandas Index or Series.
        return _json_param(tolist(), label)
    if isinstance(value, Sequence):
        return [_json_param(item, f"{label}[{index}]") for index, item in enumerate(value)]
    raise TypeError(f"{label} must be JSON-like (text, numbers, lists, dicts), not {type(value).__name__}")


def _values_param(value: Any, label: str) -> np.ndarray:
    """One number per cell, as the value type it travels as: unsigned codes
    as they are, floats as float32/float64, booleans as uint8, other integers
    as float64 (exact up to 2**53)."""
    to_numpy = getattr(value, "to_numpy", None)
    array = np.asarray(to_numpy() if callable(to_numpy) else value)
    if array.ndim != 1:
        raise TypeError(f"{label} must be one-dimensional: one number per cell")
    kind, size = array.dtype.kind, array.dtype.itemsize
    if kind == "b":
        return array.astype(np.uint8)
    if kind == "u" and size <= 4:
        return array
    if kind == "f" and size in (4, 8):
        return array
    if kind == "f" and size == 2:
        return array.astype(np.float32)
    if kind in "iuf":
        return array.astype(np.float64)
    try:
        return array.astype(np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{label} must be numbers, not {array.dtype}") from error


def action_params(
    action: str,
    params: Mapping[str, Any],
    *,
    cells: Callable[[Any, str], np.ndarray],
) -> dict[str, Any]:
    """Check ``params`` against ``action``'s parameters and convert them for
    the wire. ``OMIT`` values are left out; ``cells(value, label)`` turns an
    ``x-cells`` value into original positions (``uint32``)."""
    entry = require_action(action)
    schema = entry["params"]
    properties: dict[str, Any] = schema.get("properties", {})
    out: dict[str, Any] = {}
    for name, value in params.items():
        if value is OMIT:
            continue
        sub = properties.get(name)
        if sub is None:
            known = ", ".join(properties) or "none"
            raise TypeError(f"The action {action!r} has no parameter {name!r}; its parameters are: {known}")
        label = f"{tool_name(action)}({name}=...)"
        if value is None:
            out[name] = None
        elif _annotated(sub, "x-cells") or _annotated(sub, "x-values"):
            key = "x-cells" if _annotated(sub, "x-cells") else "x-values"
            array = cells(value, label) if key == "x-cells" else _values_param(value, label)
            fewest = _fewest_items(sub, key)
            if array.size < fewest:
                what = "cell" if key == "x-cells" else "value"
                raise ValueError(f"{label} needs at least {fewest} {plural(fewest, what)}")
            out[name] = array
        else:
            out[name] = _json_param(value, label)
    missing = [name for name in schema.get("required", []) if name not in out]
    if missing:
        raise TypeError(
            f"The action {action!r} needs {', '.join(repr(name) for name in missing)}"
        )
    return out


def cell_positions(
    value: Any, *, label: str, n_obs: int | None, obs_names: Callable[[], pd.Index | None]
) -> np.ndarray:
    """Cells as original positions (``uint32``), in the order given, each once.

    ``value`` is rows (ints, a NumPy integer array), a boolean mask of every
    cell (a list, a NumPy or pandas boolean array), a boolean ``pandas.Series``
    (its index's cells that are True), or obs names (text, a ``pandas.Index``).
    ``obs_names()`` gives the viewer's obs names (None: this process cannot
    read them); it is asked only for cells named by name, since reading the
    names of a large file costs what rows never need.
    """
    if isinstance(value, str | bytes):
        raise TypeError(f"{label} must list the cells (rows, a mask or obs names), not be one string")
    if isinstance(value, pd.Series) and value.dtype == bool:
        value = value.index[value.to_numpy()]
    if isinstance(value, pd.Series):
        value = value.to_numpy()
    if isinstance(value, np.ndarray | pd.Index):
        array = np.asarray(value)
    else:
        items = list(value)
        kinds = {bool if type(item) is bool else int if isinstance(item, int | np.integer) else str if isinstance(item, str)
                 else type(item) for item in items}
        if len(kinds) > 1 or not kinds <= {bool, int, str}:
            raise TypeError(f"{label} must be rows (ints), a mask (bools) or obs names (text), not a mix")
        array = np.asarray(items, dtype=object if kinds == {str} else None)
    if array.ndim != 1:
        raise TypeError(f"{label} must be one-dimensional")
    if array.size == 0:
        positions = np.empty(0, dtype=np.int64)
    elif array.dtype.kind == "b":
        if n_obs is not None and array.size != n_obs:
            raise ValueError(
                f"{label} is a mask of {array.size} {plural(array.size, 'cell')}; the viewer's dataset has {n_obs}"
            )
        positions = np.flatnonzero(array)
    elif array.dtype.kind in "OUS":
        names = obs_names()
        if names is None:
            raise ValueError(f"{label} names cells by obs name, and this viewer's dataset has none this process can read; name them by row")
        if not names.is_unique:
            raise ValueError(
                f"{label}: this viewer's dataset repeats obs names, so a name does not say which cell it is; "
                "name the cells by row"
            )
        positions = names.get_indexer(pd.Index(array))
        if (positions < 0).any():
            missing = array[positions < 0]
            raise ValueError(f"{label}: {missing.size} of the cells named are not in the viewer's dataset (first: {missing[0]!r})")
    elif array.dtype.kind in "iu":
        positions = array.astype(np.int64)
    else:
        raise TypeError(f"{label} must be rows, a boolean mask or obs names, not {array.dtype}")
    if positions.size:
        if int(positions.min()) < 0:
            raise ValueError(f"{label} must be non-negative rows")
        if int(positions.max()) > _MAX_CELL_INDEX or (n_obs is not None and int(positions.max()) >= n_obs):
            limit = "2**32" if n_obs is None else f"the dataset's {n_obs} {plural(n_obs, 'cell')}"
            raise ValueError(f"{label} names row {int(positions.max())}, beyond {limit}")
        if np.unique(positions).size != positions.size:
            raise ValueError(f"{label} names a cell more than once")
    return positions.astype(np.uint32)
