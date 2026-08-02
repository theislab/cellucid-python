"""
Planning the ``obs`` columns that carry a session's user-defined fields.

A user-defined field is either categorical -- codes plus categories, decoded
from its own chunk -- or a continuous alias of an existing column. Both are
validated against the field overlay the session declares before a column is
built, including the operation that produced it, its creation timestamp, and,
for categoricals, the category list and the overlap strategy that resolved
cells belonging to more than one source page.
"""

from __future__ import annotations

import copy
from collections.abc import Sequence
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy import sparse

from ..session_codecs import decode_user_defined_codes
from ._primitives import (
    _require_category_list,
    _require_exact_bool,
    _require_exact_keys,
    _require_field_key,
    _require_nonempty_string,
    _require_nonnegative_int,
    _require_optional_string,
    _require_wire_id,
    _reserve_column,
)
from ._records import _ColumnPlan
from ._schema import (
    _CATEGORY_FIELD_KEYS,
    _CONTINUOUS_FIELD_KEYS,
    _OVERLAP_STRATEGIES,
    _OVERLAY_ROOT_KEYS,
    _USER_DEFINED_CODES_PREFIX,
)


def _validate_source_field(value: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    allowed = {"sourceKey", "sourceIndex", "kind"}
    if "sourceKey" not in value or set(value) - allowed:
        raise ValueError(f"{label} has invalid fields")
    source_key = _require_nonempty_string(
        value["sourceKey"],
        label=f"{label}.sourceKey",
    )
    output: dict[str, Any] = {"sourceKey": source_key}
    if "sourceIndex" in value:
        output["sourceIndex"] = _require_nonnegative_int(
            value["sourceIndex"],
            label=f"{label}.sourceIndex",
        )
    if "kind" in value:
        output["kind"] = _require_nonempty_string(
            value["kind"],
            label=f"{label}.kind",
        )
    return output


def _validate_operation(value: Any, *, label: str) -> None:
    if value is None:
        return
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object or null")
    _require_nonempty_string(value.get("type"), label=f"{label}.type")


def _numeric_series(
    values: Any,
    *,
    n_obs: int,
    obs_index: pd.Index,
    label: str,
) -> pd.Series:
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
        raise TypeError(f"{label} must contain real numeric values")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} must contain only finite values")
    return pd.Series(array.copy(), index=obs_index)


def _materialize_continuous_alias(
    adata: Any,
    *,
    source: str,
    source_field: dict[str, Any],
    label: str,
) -> pd.Series:
    source_key = source_field["sourceKey"]
    n_obs = int(adata.n_obs)
    obs_index = adata.obs_names
    if source == "obs":
        matching = [index for index, key in enumerate(adata.obs.columns) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} source obs field {source_key!r} must exist exactly once")
        source_index = matching[0]
        values = adata.obs.iloc[:, source_index].to_numpy(copy=True)
    else:
        matching = [index for index, key in enumerate(adata.var_names) if key == source_key]
        if len(matching) != 1:
            raise ValueError(f"{label} source gene {source_key!r} must exist exactly once")
        source_index = matching[0]
        if getattr(adata, "X", None) is None:
            raise ValueError(f"{label} cannot copy gene values because adata.X is absent")
        values = adata.X[:, source_index]

    declared_index = source_field.get("sourceIndex")
    if declared_index is not None and declared_index != source_index:
        raise ValueError(
            f"{label} sourceIndex {declared_index} does not match exact source "
            f"{source_key!r} at index {source_index}"
        )
    return _numeric_series(
        values,
        n_obs=n_obs,
        obs_index=obs_index,
        label=label,
    )


def _validate_category_metadata(
    field: dict[str, Any],
    *,
    label: str,
) -> list[str | bool | int | float]:
    categories = _require_category_list(
        field["categories"],
        label=f"{label}.categories",
        nonempty=True,
    )
    codes_length = _require_nonnegative_int(
        field["codesLength"],
        label=f"{label}.codesLength",
    )
    codes_type = field["codesType"]
    if codes_type not in {"Uint8Array", "Uint16Array"}:
        raise ValueError(f"{label}.codesType must be Uint8Array or Uint16Array")
    if not isinstance(field["centroidsByDim"], dict):
        raise TypeError(f"{label}.centroidsByDim must be an object")
    normalized_dims = field["normalizedDims"]
    if (
        not isinstance(normalized_dims, list)
        or any(
            isinstance(item, bool) or not isinstance(item, int) or item not in {1, 2, 3}
            for item in normalized_dims
        )
        or len(normalized_dims) != len(set(normalized_dims))
    ):
        raise ValueError(f"{label}.normalizedDims must be unique dimensions 1, 2, or 3")
    if not isinstance(field["sourcePages"], list):
        raise TypeError(f"{label}.sourcePages must be an array")
    if field["overlapStrategy"] not in _OVERLAP_STRATEGIES:
        raise ValueError(f"{label}.overlapStrategy is invalid")
    overlap_label = _require_optional_string(
        field["overlapLabel"],
        label=f"{label}.overlapLabel",
    )
    if overlap_label == "":
        raise ValueError(f"{label}.overlapLabel must be non-empty or null")
    if field["intersectionLabels"] is not None and not isinstance(
        field["intersectionLabels"], dict
    ):
        raise TypeError(f"{label}.intersectionLabels must be an object or null")
    uncovered_label = _require_optional_string(
        field["uncoveredLabel"],
        label=f"{label}.uncoveredLabel",
    )
    if uncovered_label == "":
        raise ValueError(f"{label}.uncoveredLabel must be non-empty or null")
    if field["sourceField"] is not None:
        _validate_source_field(field["sourceField"], label=f"{label}.sourceField")
    _validate_operation(field["operation"], label=f"{label}.operation")
    if codes_length > 0 and codes_type == "Uint8Array" and len(categories) > 255:
        raise ValueError(f"{label} has too many categories for Uint8Array codes")
    return categories


def _validate_field_created_at(value: Any, *, label: str) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{label} must be a finite non-negative number or null")
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{label} must be a finite non-negative number or null")


def _validate_string_record(value: Any, *, label: str) -> None:
    if not isinstance(value, dict):
        raise TypeError(f"{label} must be an object")
    for key, item in value.items():
        _require_nonempty_string(key, label=f"{label} key")
        _require_nonempty_string(item, label=f"{label} value for {key!r}")


def _validate_unique_string_array(value: Any, *, label: str) -> list[str]:
    if not isinstance(value, list):
        raise TypeError(f"{label} must be an array")
    output: list[str] = []
    seen: set[str] = set()
    for index, item in enumerate(value):
        exact = _require_nonempty_string(item, label=f"{label}[{index}]")
        if exact in seen:
            raise ValueError(f"{label} duplicates {exact!r}")
        seen.add(exact)
        output.append(exact)
    return output


def _plan_user_defined_fields(
    bundle: Any,
    chunk_ids: set[str],
    adata: Any,
    *,
    prefix: str,
    existing_columns: set[str],
    materialize: bool,
    include_deleted: bool,
    metadata_by_id: dict[str, dict[str, Any]],
) -> tuple[list[_ColumnPlan], Any | None]:
    codes_prefix = "user-defined/codes/"
    code_chunk_ids = {chunk_id for chunk_id in chunk_ids if chunk_id.startswith(codes_prefix)}
    if "core/field-overlays" not in chunk_ids:
        if code_chunk_ids:
            raise ValueError(
                "Session contains user-defined code chunks without core/field-overlays"
            )
        return [], None

    raw_overlays = bundle.decode_chunk("core/field-overlays")
    overlays = _require_exact_keys(
        raw_overlays,
        _OVERLAY_ROOT_KEYS,
        label="core/field-overlays",
    )
    renames = _require_exact_keys(
        overlays["renames"],
        {"fields", "categories"},
        label="core/field-overlays.renames",
    )
    _validate_string_record(
        renames["fields"],
        label="core/field-overlays.renames.fields",
    )
    _validate_string_record(
        renames["categories"],
        label="core/field-overlays.renames.categories",
    )
    deleted_fields = _require_exact_keys(
        overlays["deletedFields"],
        {"deleted", "purged"},
        label="core/field-overlays.deletedFields",
    )
    deleted = set(
        _validate_unique_string_array(
            deleted_fields["deleted"],
            label="core/field-overlays.deletedFields.deleted",
        )
    )
    purged = set(
        _validate_unique_string_array(
            deleted_fields["purged"],
            label="core/field-overlays.deletedFields.purged",
        )
    )
    if not purged.issubset(deleted):
        raise ValueError("core/field-overlays.deletedFields.purged must be a subset of deleted")
    fields = overlays["userDefinedFields"]
    if not isinstance(fields, list):
        raise TypeError("core/field-overlays.userDefinedFields must be an array")

    plans: list[_ColumnPlan] = []
    field_ids: set[str] = set()
    referenced_code_chunks: set[str] = set()
    categorical_field_ids: list[str] = []
    n_obs = int(adata.n_obs)
    for index, raw_field in enumerate(fields):
        if not isinstance(raw_field, dict):
            raise TypeError(f"core/field-overlays.userDefinedFields[{index}] must be an object")
        label = f"core/field-overlays.userDefinedFields[{index}]"
        kind = raw_field.get("kind")
        expected_keys = _CATEGORY_FIELD_KEYS if kind == "category" else _CONTINUOUS_FIELD_KEYS
        field = _require_exact_keys(raw_field, expected_keys, label=label)
        if kind not in {"category", "continuous"}:
            raise ValueError(f"{label}.kind must be 'category' or 'continuous'")
        field_id = _require_wire_id(field["id"], label=f"{label}.id")
        if field_id in field_ids:
            raise ValueError(f"Duplicate user-defined field id: {field_id}")
        field_ids.add(field_id)
        source = field["source"]
        if source not in {"obs", "var"}:
            raise ValueError(f"{label}.source must be 'obs' or 'var'")
        field_key = _require_field_key(field["key"], label=f"{label}.key")
        is_deleted = _require_exact_bool(
            field["isDeleted"],
            label=f"{label}.isDeleted",
        )
        is_purged = _require_exact_bool(
            field["isPurged"],
            label=f"{label}.isPurged",
        )
        if is_purged and not is_deleted:
            raise ValueError(f"{label} purged fields must also be deleted")
        _validate_field_created_at(field["createdAt"], label=f"{label}.createdAt")
        should_materialize = materialize and (include_deleted or not is_deleted) and not is_purged

        column_name = f"{prefix}{field_key}"
        if kind == "continuous":
            source_field = _validate_source_field(
                field["sourceField"],
                label=f"{label}.sourceField",
            )
            _validate_operation(field["operation"], label=f"{label}.operation")
            if should_materialize:
                _reserve_column(existing_columns, column_name)
                values = _materialize_continuous_alias(
                    adata,
                    source=source,
                    source_field=source_field,
                    label=f"Continuous field {field_id!r}",
                )
                plans.append(
                    _ColumnPlan(
                        name=column_name,
                        values=values,
                        metadata={
                            "kind": "continuous",
                            "field_id": field_id,
                            "source": source,
                            "source_field": source_field,
                        },
                    )
                )
            continue

        categories = _validate_category_metadata(field, label=label)
        categorical_field_ids.append(field_id)
        codes_length = field["codesLength"]
        codes_type = field["codesType"]
        code_chunk_id = f"{codes_prefix}{field_id}"
        if codes_length != n_obs:
            raise ValueError(
                f"User-defined categorical field {field_id!r} has codesLength "
                f"{codes_length}, expected {n_obs} cell-aligned codes"
            )
        if n_obs == 0:
            if code_chunk_id in chunk_ids:
                raise ValueError(
                    f"Empty user-defined field {field_id!r} must not contain a codes chunk"
                )
            codes = np.empty(0, dtype=np.uint8 if codes_type == "Uint8Array" else np.uint16)
        else:
            if code_chunk_id not in chunk_ids:
                raise ValueError(
                    f"User-defined categorical field {field_id!r} is missing "
                    f"required chunk {code_chunk_id!r}"
                )
            code_meta = metadata_by_id[code_chunk_id]
            if code_meta["label"] != f"User-defined codes: {field_key}":
                raise ValueError(
                    f"{code_chunk_id} label must match user-defined field {field_id!r}"
                )
            raw_codes = bundle.decode_chunk(code_chunk_id)
            if not isinstance(raw_codes, bytes | bytearray | memoryview):
                raise TypeError(f"{code_chunk_id} must decode to binary bytes")
            if code_meta["uncompressedBytes"] != len(raw_codes):
                raise ValueError(f"{code_chunk_id} uncompressedBytes must match its payload")
            codes = decode_user_defined_codes(
                raw_codes,
                expected_length=codes_length,
                expected_codes_type=codes_type,
            )
            referenced_code_chunks.add(code_chunk_id)
            expected_dtype = np.dtype(np.uint8 if codes_type == "Uint8Array" else np.uint16)
            if codes.dtype != expected_dtype:
                raise ValueError(f"{code_chunk_id} dtype {codes.dtype} does not match {codes_type}")
            if len(codes) != codes_length:
                raise ValueError(
                    f"{code_chunk_id} contains {len(codes)} codes, expected {codes_length}"
                )

        category_codes = codes.astype(np.int64, copy=True)
        missing_sentinel = (
            (255 if codes_type == "Uint8Array" else 65_535)
            if field["uncoveredLabel"] is None
            else None
        )
        if missing_sentinel is not None:
            category_codes[category_codes == missing_sentinel] = -1
        invalid = np.flatnonzero((category_codes < -1) | (category_codes >= len(categories)))
        if invalid.size:
            position = int(invalid[0])
            code = int(category_codes[position])
            raise ValueError(
                f"User-defined categorical field {field_id!r} contains code "
                f"{code} at position {position}, but declares {len(categories)} categories"
            )

        if should_materialize:
            _reserve_column(existing_columns, column_name)
            categorical = pd.Categorical.from_codes(
                cast("Sequence[int]", category_codes),
                categories=pd.Index(categories),
                ordered=False,
            )
            plans.append(
                _ColumnPlan(
                    name=column_name,
                    values=pd.Series(categorical, index=adata.obs_names),
                    metadata={
                        "kind": "category",
                        "field_id": field_id,
                        "source": source,
                    },
                )
            )

    if code_chunk_ids != referenced_code_chunks:
        unexpected = sorted(code_chunk_ids - referenced_code_chunks)
        raise ValueError(
            "Session contains undeclared user-defined code chunks: " + ", ".join(unexpected)
        )
    # The AnnData bridge deliberately leaves core/state opaque. It proves
    # completeness and exact field-overlay order inside each priority bucket;
    # the web restore owner proves which active live/snapshot fields are eager.
    for priority in ("eager", "lazy"):
        actual_order = [
            chunk_id.removeprefix(_USER_DEFINED_CODES_PREFIX)
            for chunk_id, meta in metadata_by_id.items()
            if (chunk_id.startswith(_USER_DEFINED_CODES_PREFIX) and meta["priority"] == priority)
        ]
        expected_order = [
            field_id for field_id in categorical_field_ids if field_id in actual_order
        ]
        if actual_order != expected_order:
            raise ValueError(
                f"User-defined {priority} code chunks must match the exact "
                "field-overlay inventory order"
            )
    return plans, copy.deepcopy(overlays)
