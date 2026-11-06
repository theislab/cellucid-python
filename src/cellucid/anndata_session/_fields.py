"""
Planning the ``obs`` columns that carry a session's user-defined fields.

The ``field-overlays`` value at a step (``app/state/domains/field-overlays.js``)
lists every field a person or a notebook made (``userFields``), each one of:

- categorical: ``categories`` and ``codes``, a per-cell column ``{cells,
  values}`` by original position (``cells`` null: every cell); a code is a
  category's index, and its code type's largest value marks a cell with no
  category. The column is a pandas ``category`` with those categories in
  order, missing where a cell has no category;
- continuous with its own values (a column a notebook pushed): ``values``, a
  per-cell column of ``Float32Array`` values, NaN for no value. The column is
  ``float32``, NaN where a cell has no value;
- a continuous copy of another field (``sourceField``, no values of its own):
  the viewer reads its values from the field it copies, and so does this: an
  obs column, or a gene named the way the session's data mapping names it
  (its gene-name set on its gene axis, made unique with the shared rule),
  copied from the expression matrix the mapping reads.

Each column is named by the field's key (as shown). A field deleted in the
session, and a purged one (which holds no values), is not applied, and says
so. A cell the session does not know (one the viewer did not open) is
missing.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .._plural import plural
from ..data_mapping import mapped_expression_matrix, resolve_gene_names
from ._columns import categorical_column, continuous_column, without_unknown
from ._primitives import (
    _require_array,
    _require_bool,
    _require_category_list,
    _require_cell_list,
    _require_exact_keys,
    _require_field_key,
    _require_list,
    _require_optional_text,
    _require_record,
    _require_source,
    _require_text,
)
from ._records import _Candidate
from ._schema import (
    CODE_TYPES,
    COLUMN_KEYS,
    OVERLAP_STRATEGIES,
    SOURCE_FIELD_KEYS,
    USER_CATEGORY_KEYS,
    USER_FIELD_KEYS,
    USER_VALUES_KEYS,
)


def _source_field(value: Any, *, label: str) -> dict[str, Any]:
    record = _require_record(value, label=label)
    if "sourceKey" not in record or set(record) - SOURCE_FIELD_KEYS:
        raise ValueError(f"This session's {label} must name sourceKey, and may name sourceIndex and kind only.")
    out: dict[str, Any] = {"sourceKey": _require_text(record["sourceKey"], label=f"{label}.sourceKey")}
    if "sourceIndex" in record:
        index = record["sourceIndex"]
        if isinstance(index, bool) or not isinstance(index, int) or index < 0:
            raise ValueError(f"This session's {label}.sourceIndex must be a non-negative integer.")
        out["sourceIndex"] = index
    if "kind" in record:
        out["kind"] = _require_text(record["kind"], label=f"{label}.kind")
    return out


def _operation(value: Any, *, label: str) -> None:
    if value is None:
        return
    record = _require_record(value, label=label)
    _require_text(record.get("type"), label=f"{label}.type")


def _created_at(value: Any, *, label: str) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int | float) or not np.isfinite(value) or value < 0:
        raise ValueError(f"This session's {label} must be a finite non-negative number or null.")


def _column_part(value: Any, types: tuple[str, ...], n_obs: int, *, label: str) -> tuple[np.ndarray | None, np.ndarray, str]:
    """A per-cell column ``{cells, values}``: its cells (None: every cell),
    its values, and their typed-array type."""
    column = _require_exact_keys(value, COLUMN_KEYS, label=label)
    cells = None if column["cells"] is None else _require_cell_list(column["cells"], n_obs, label=f"{label}.cells")
    values = _require_array(column["values"], types, label=f"{label}.values")
    expected = n_obs if cells is None else cells.size
    if values.data.size != expected:
        held = int(values.data.size)
        raise ValueError(
            f"This session's {label} holds {held} {plural(held, 'value')} for "
            + ("every one of the dataset's" if cells is None else "its")
            + f" {expected} {plural(expected, 'cell')}."
        )
    return cells, values.data, values.type


def _numeric_series(values: Any, *, n_obs: int, index: pd.Index, label: str) -> pd.Series:
    array = np.asarray(values.toarray()) if sparse.issparse(values) else np.asarray(values)
    if array.ndim == 2 and array.shape[1] == 1:
        array = array[:, 0]
    if array.ndim != 1 or array.shape[0] != n_obs:
        raise ValueError(f"{label} must contain exactly {n_obs} cell values")
    if (
        not np.issubdtype(array.dtype, np.number)
        or np.issubdtype(array.dtype, np.bool_)
        or np.issubdtype(array.dtype, np.complexfloating)
    ):
        raise ValueError(f"{label} must contain real numeric values")
    # Copied as the viewer copies them (``overlay-public.js``): NaN, the
    # stored missing value -- and every cell a MuData modality did not
    # measure -- stays NaN.
    return pd.Series(array.astype(np.float64 if array.dtype == np.float64 else np.float32), index=index)


def _copied_values(shown: Any, *, mapping: dict[str, Any], source: str, source_field: dict[str, Any], label: str) -> pd.Series:
    """The values of a continuous copy, read from the field it copies."""
    source_key = source_field["sourceKey"]
    n_obs = int(shown.n_obs)
    if source == "obs":
        matching = [position for position, key in enumerate(shown.obs.columns) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} copies the obs column {source_key!r}, which must exist exactly once")
        source_index = matching[0]
        values = shown.obs.iloc[:, source_index].to_numpy(copy=True)
    else:
        # A gene is named the way the viewer named it: by the data mapping's
        # gene-name set on its gene axis, made unique with the shared rule, so
        # the second of two 'MATR3' rows is 'MATR3-1' here exactly as it was
        # in the session. Its values come from the matrix the mapping read.
        gene_names = resolve_gene_names(shown, mapping).names
        matching = [position for position, key in enumerate(gene_names) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} copies the gene {source_key!r}, which must exist exactly once")
        source_index = matching[0]
        try:
            matrix = mapped_expression_matrix(shown, mapping)
        except ValueError as error:
            raise ValueError(f"{label} cannot copy gene values: {error}") from error
        values = matrix[:, source_index]
    declared = source_field.get("sourceIndex")
    if declared is not None and declared != source_index:
        raise ValueError(
            f"{label} sourceIndex {declared} does not match exact source {source_key!r} at index {source_index}"
        )
    return _numeric_series(values, n_obs=n_obs, index=shown.obs_names, label=label)


def _check_category_definition(field: dict[str, Any], *, label: str) -> list[Any]:
    categories = _require_category_list(field["categories"], label=f"{label}.categories")
    _require_list(field["sourcePages"], label=f"{label}.sourcePages")
    if field["overlapStrategy"] not in OVERLAP_STRATEGIES:
        raise ValueError(f"This session's {label}.overlapStrategy is not one this version reads.")
    _require_optional_text(field["overlapLabel"], label=f"{label}.overlapLabel")
    if field["intersectionLabels"] is not None:
        _require_record(field["intersectionLabels"], label=f"{label}.intersectionLabels")
    _require_optional_text(field["uncoveredLabel"], label=f"{label}.uncoveredLabel")
    return categories


#: Why a field the session made is not applied.
DELETED_NOTE = "Deleted in the session."
PURGED_NOTE = "Deleted for good in the session: it holds no values."


def _plan_user_fields(
    fields: Any,
    shown: Any,
    *,
    mapping: dict[str, Any],
    known: np.ndarray | None,
) -> list[_Candidate]:
    """One column per user field, every field checked before any column is built."""
    n_obs = int(shown.n_obs)
    index = shown.obs_names
    ids: set[str] = set()
    candidates: list[_Candidate] = []
    for position, raw in enumerate(_require_list(fields, label="field-overlays value userFields")):
        label = f"field-overlays value userFields[{position}]"
        record = _require_record(raw, label=label)
        kind = record.get("kind")
        if kind not in ("category", "continuous"):
            raise ValueError(f"This session's {label}.kind must be \"category\" or \"continuous\".")
        purged = _require_bool(record.get("purged"), label=f"{label}.purged")
        copy_of = kind == "continuous" and record.get("sourceField") is not None
        keys = (
            USER_FIELD_KEYS if purged or copy_of
            else USER_CATEGORY_KEYS if kind == "category"
            else USER_VALUES_KEYS
        )
        field = _require_exact_keys(record, keys, label=label)
        field_id = _require_text(field["id"], label=f"{label}.id")
        if field_id in ids:
            raise ValueError(f"This session's field-overlays value holds the user field {field_id!r} twice.")
        ids.add(field_id)
        source = _require_source(field["source"], label=f"{label}.source")
        key = _require_field_key(field["key"], label=f"{label}.key")
        deleted = _require_bool(field["deleted"], label=f"{label}.deleted")
        if purged and not deleted:
            raise ValueError(f"This session's {label} is purged and not deleted.")
        _created_at(field["createdAt"], label=f"{label}.createdAt")
        _operation(field["operation"], label=f"{label}.operation")
        source_field = None if field["sourceField"] is None else _source_field(field["sourceField"], label=f"{label}.sourceField")
        if copy_of and source_field is None:
            raise ValueError(f"This session's {label} copies no field.")
        metadata: dict[str, Any] = {
            "kind": kind, "field_id": field_id, "key": key, "source": source, "deleted": deleted,
            "source_field": source_field,
        }
        what = f"The user field {key!r} ({field_id})"

        def unapplied(note: str, metadata: dict[str, Any] = metadata, key: str = key) -> _Candidate:
            return _Candidate(kind="user-field", domain="field-overlays", name=key, values=None, cells=None, note=note, metadata=metadata)

        if purged:
            candidates.append(unapplied(PURGED_NOTE))
            continue
        if kind == "category":
            categories = _check_category_definition(field, label=label)
            cells, codes, code_type = _column_part(field["codes"], tuple(CODE_TYPES), n_obs, label=f"{label}.codes")
            missing = CODE_TYPES[code_type]
            if len(categories) > missing:
                raise ValueError(f"This session's {label} has more categories than its code type ({code_type}) can name.")
            codes64 = codes.astype(np.int64)
            invalid = np.flatnonzero((codes64 >= len(categories)) & (codes64 != missing))
            if invalid.size:
                count = len(categories)
                raise ValueError(
                    f"This session's {label}.codes holds the code {int(codes64[invalid[0]])}, and the field has "
                    f"{count} {plural(count, 'category', 'categories')}."
                )
            if deleted:
                candidates.append(unapplied(DELETED_NOTE))
                continue
            values = categorical_column(codes, categories, cells=cells, n_obs=n_obs, missing_code=missing, index=index)
        elif copy_of:
            assert source_field is not None
            if deleted:
                candidates.append(unapplied(DELETED_NOTE))
                continue
            values = _copied_values(shown, mapping=mapping, source=source, source_field=source_field, label=what)
        else:
            cells, numbers, _ = _column_part(field["values"], ("Float32Array",), n_obs, label=f"{label}.values")
            if bool(np.isinf(numbers).any()):
                raise ValueError(f"This session's {label}.values holds an infinity; a field holds finite values or NaN.")
            if deleted:
                candidates.append(unapplied(DELETED_NOTE))
                continue
            values = continuous_column(numbers, cells=cells, n_obs=n_obs, index=index)
        known_values = without_unknown(values, known)
        candidates.append(_Candidate(
            kind="user-field", domain="field-overlays", name=key, values=known_values,
            cells=int(known_values.notna().sum()), note=None, metadata=metadata,
        ))
    return candidates
