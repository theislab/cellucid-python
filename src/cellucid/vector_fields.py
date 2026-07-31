"""Exact-current vector-field validation and transition-drift utilities."""

from __future__ import annotations

import math
import re
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from scipy import sparse

if TYPE_CHECKING:
    import anndata


_VECTOR_KEY_PATTERN = re.compile(r"^(?P<field>.+_umap)_(?P<dimension>[123])d$")

SUPPORTED_EMBEDDING_KEYS: tuple[str, ...] = ("X_umap_1d", "X_umap_2d", "X_umap_3d")


def is_vector_field_declaration_key(value: object) -> bool:
    """Return whether an AnnData ``obsm`` key exactly declares a vector field."""
    return isinstance(value, str) and _VECTOR_KEY_PATTERN.fullmatch(value) is not None


def embedding_declaration_hints(obsm: Any, *, n_cells: int) -> list[str]:
    """Return one exact rename per undeclared obsm array that could be an embedding.

    ``sc.tl.umap()`` writes ``adata.obsm['X_umap']`` without a dimensional
    suffix, so the most common first failure is one declaration away from
    working. Each hint is a complete, copy-pasteable statement.
    """
    hints: list[str] = []
    for key in obsm:
        if not isinstance(key, str):
            continue
        if key in SUPPORTED_EMBEDDING_KEYS or is_vector_field_declaration_key(key):
            continue
        shape = getattr(obsm[key], "shape", None)
        if not isinstance(shape, tuple) or len(shape) != 2:
            continue
        rows, columns = shape
        if int(rows) != n_cells or int(columns) not in (1, 2, 3):
            continue
        hints.append(
            f'adata.obsm["X_umap_{int(columns)}d"] = adata.obsm[{key!r}]'
            f"  # {key!r} has shape {(int(rows), int(columns))}"
        )
    # Scanpy's own key is the one a first-time user almost always means.
    hints.sort(key=lambda hint: (0 if "'X_umap'" in hint else 1, hint))
    return hints


@dataclass(frozen=True)
class ValidatedVectorFields:
    """Validated vector arrays and their exact source keys."""

    fields: dict[str, dict[int, np.ndarray]]
    source_keys: dict[str, dict[int, str]]
    default_field: str | None


def _require_nonempty_string(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a native non-empty string.")
    if not value:
        raise ValueError(f"{label} must be a native non-empty string.")
    return value


def _dense_array(values: object) -> np.ndarray:
    if sparse.issparse(values):
        return np.asarray(cast(sparse.spmatrix, values).toarray())
    return np.asarray(values)


def _finite_float32_matrix(
    values: object,
    *,
    label: str,
    n_rows: int | None,
    declared_dimension: int | None,
) -> tuple[np.ndarray, int]:
    if values is None:
        raise ValueError(f"{label} must not be None.")
    array = _dense_array(values)
    if array.dtype.kind not in {"i", "u", "f"}:
        raise TypeError(f"{label} must contain real numeric values.")
    if array.ndim == 1 and declared_dimension == 1:
        # A declared 1D quantity carries one coordinate per cell, so a plain
        # length-n vector already is the (n, 1) column it declares. AnnData
        # performs exactly this reshape when such an array is put in obsm.
        array = array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"{label} must be a 2D array, got shape {array.shape}.")
    if array.shape[0] <= 0:
        raise ValueError(f"{label} must contain at least one row.")
    if n_rows is not None and array.shape[0] != n_rows:
        raise ValueError(f"{label} must have exactly {n_rows} rows, got {array.shape[0]}.")
    if declared_dimension is None:
        dimension = int(array.shape[1])
        if dimension not in (1, 2, 3):
            raise ValueError(f"{label} must have exactly 1, 2, or 3 columns, got {array.shape[1]}.")
    else:
        dimension = declared_dimension
        if array.shape[1] != dimension:
            raise ValueError(
                f"{label} must have exactly {dimension} columns, got {array.shape[1]}."
            )
    if not np.isfinite(array).all():
        raise ValueError(f"{label} must contain only finite values.")
    with np.errstate(over="ignore", invalid="ignore"):
        float32_array = array.astype(np.float32, copy=False)
    if not np.isfinite(float32_array).all():
        raise ValueError(f"{label} contains values outside the finite float32 range.")
    return float32_array, dimension


def validate_vector_fields(
    values: Mapping[Any, Any] | None,
    *,
    n_cells: int,
    available_dimensions: Collection[int],
    vector_field_default: object = None,
) -> ValidatedVectorFields:
    """Validate the complete mapping used by highest-dimension field metadata."""
    if values is not None and not isinstance(values, Mapping):
        raise TypeError("vector_fields must be a mapping of vector field keys to arrays.")

    if vector_field_default is not None:
        declared_default = _require_nonempty_string(
            vector_field_default,
            label="vector_field_default",
        )
    else:
        declared_default = None

    if not values:
        if declared_default is not None:
            raise ValueError("vector_field_default was provided but no vector fields exist.")
        return ValidatedVectorFields({}, {}, None)

    allowed_dimensions = set(available_dimensions)
    fields: dict[str, dict[int, np.ndarray]] = {}
    source_keys: dict[str, dict[int, str]] = {}

    for raw_key, raw_values in values.items():
        key = _require_nonempty_string(
            raw_key,
            label="vector field key",
        )
        match = _VECTOR_KEY_PATTERN.fullmatch(key)
        if match is None:
            raise ValueError(
                f"vector field key {key!r} must exactly match '<field>_umap_<1|2|3>d'."
            )

        field_id = match.group("field")
        dimension_text = match.group("dimension")
        declared_dimension = int(dimension_text)
        array, dimension = _finite_float32_matrix(
            raw_values,
            label=f"Vector field {key!r}",
            n_rows=n_cells,
            declared_dimension=declared_dimension,
        )

        if dimension not in allowed_dimensions:
            raise ValueError(f"Vector field {key!r} requires a matching {dimension}D embedding.")
        if dimension in fields.setdefault(field_id, {}):
            other_key = source_keys[field_id][dimension]
            raise ValueError(
                f"Vector field keys {other_key!r} and {key!r} both declare "
                f"the same {dimension}D field {field_id!r}."
            )
        fields[field_id][dimension] = array
        source_keys.setdefault(field_id, {})[dimension] = key

    ordered_fields = {
        field_id: {dimension: fields[field_id][dimension] for dimension in sorted(fields[field_id])}
        for field_id in sorted(fields)
    }
    ordered_source_keys = {
        field_id: {
            dimension: source_keys[field_id][dimension]
            for dimension in sorted(source_keys[field_id])
        }
        for field_id in sorted(source_keys)
    }

    if declared_default is None:
        if len(ordered_fields) > 1:
            raise ValueError(
                "vector_field_default is required when more than one vector field is declared."
            )
        default_field = next(iter(ordered_fields))
    else:
        if declared_default not in ordered_fields:
            raise ValueError(
                f"vector_field_default {declared_default!r} does not match an "
                f"available field: {list(ordered_fields)}."
            )
        default_field = declared_default

    return ValidatedVectorFields(
        fields=ordered_fields,
        source_keys=ordered_source_keys,
        default_field=default_field,
    )


def scale_vector_field(
    vectors: np.ndarray,
    *,
    scale_factor: object,
    label: str,
) -> np.ndarray:
    """Scale vectors using one required finite positive scale."""
    if isinstance(scale_factor, bool) or not isinstance(scale_factor, Real):
        raise ValueError(f"{label} requires a finite positive scale factor.")
    numeric_scale = float(scale_factor)
    if not math.isfinite(numeric_scale) or numeric_scale <= 0:
        raise ValueError(f"{label} requires a finite positive scale factor.")
    scaled = vectors.astype(np.float64) * numeric_scale
    if not np.isfinite(scaled).all():
        raise ValueError(f"{label} scaling produced non-finite values.")
    with np.errstate(over="ignore", invalid="ignore"):
        float32_scaled = scaled.astype(np.float32)
    if not np.isfinite(float32_scaled).all():
        raise ValueError(f"{label} scaling produced values outside the finite float32 range.")
    return float32_scaled


def _transition_matrix(
    values: np.ndarray | sparse.spmatrix,
    *,
    n_cells: int,
) -> np.ndarray | sparse.csr_matrix:
    if sparse.issparse(values):
        matrix = cast(sparse.spmatrix, values).tocsr()
        if matrix.dtype.kind not in {"i", "u", "f"}:
            raise TypeError("transition_matrix must contain real numeric values.")
        if matrix.shape != (n_cells, n_cells):
            raise ValueError(
                f"transition_matrix must have shape ({n_cells}, {n_cells}), got {matrix.shape}."
            )
        if not np.isfinite(matrix.data).all():
            raise ValueError("transition_matrix must contain only finite values.")
        if np.any(matrix.data < 0):
            raise ValueError("transition_matrix values must be non-negative.")
        return matrix.astype(np.float64)

    matrix = np.asarray(values)
    if matrix.dtype.kind not in {"i", "u", "f"}:
        raise TypeError("transition_matrix must contain real numeric values.")
    if matrix.ndim != 2 or matrix.shape != (n_cells, n_cells):
        raise ValueError(
            f"transition_matrix must have shape ({n_cells}, {n_cells}), got {matrix.shape}."
        )
    if not np.isfinite(matrix).all():
        raise ValueError("transition_matrix must contain only finite values.")
    if np.any(matrix < 0):
        raise ValueError("transition_matrix values must be non-negative.")
    return matrix.astype(np.float64, copy=False)


def compute_transition_drift(
    transition_matrix: np.ndarray | sparse.spmatrix,
    embedding: np.ndarray,
    *,
    normalize_rows: bool,
) -> np.ndarray:
    """Compute exact per-cell transition drift in embedding space."""
    if type(normalize_rows) is not bool:
        raise TypeError("normalize_rows must be exactly True or False.")
    embedding_array, _dimension = _finite_float32_matrix(
        embedding,
        label="embedding",
        n_rows=None,
        declared_dimension=None,
    )
    n_cells = int(embedding_array.shape[0])
    matrix = _transition_matrix(transition_matrix, n_cells=n_cells)
    working_embedding = embedding_array.astype(np.float64)
    product = np.asarray(matrix @ working_embedding, dtype=np.float64)

    if normalize_rows:
        row_sums = np.asarray(matrix.sum(axis=1), dtype=np.float64).reshape(-1)
        if not np.isfinite(row_sums).all() or np.any(row_sums == 0):
            raise ValueError(
                "Every transition_matrix row must have a nonzero finite row sum "
                "when normalize_rows=True."
            )
        product = product / row_sums[:, None]

    drift = product - working_embedding
    drift_array, _dimension = _finite_float32_matrix(
        drift,
        label="computed transition drift",
        n_rows=n_cells,
        declared_dimension=int(embedding_array.shape[1]),
    )
    return drift_array


def _umap_embeddings(
    adata: anndata.AnnData,
) -> dict[int, tuple[str, np.ndarray]]:
    embeddings: dict[int, tuple[str, np.ndarray]] = {}
    n_cells = int(adata.n_obs)
    for dimension in (1, 2, 3):
        key = f"X_umap_{dimension}d"
        if key not in adata.obsm:
            continue
        array, _resolved_dimension = _finite_float32_matrix(
            adata.obsm[key],
            label=f"Embedding {key!r}",
            n_rows=n_cells,
            declared_dimension=dimension,
        )
        embeddings[dimension] = (key, array)

    if not embeddings:
        message = (
            "AnnData must contain one or more exact UMAP embedding keys: "
            "'X_umap_1d', 'X_umap_2d', or 'X_umap_3d'. "
            f"Available obsm keys: {list(adata.obsm.keys())}."
        )
        hints = embedding_declaration_hints(adata.obsm, n_cells=n_cells)
        if hints:
            joined = "\n".join(f"    {hint}" for hint in hints)
            message += (
                "\n  Fix: declare an existing array under the key for its "
                f"column count:\n{joined}"
            )
        raise ValueError(message)
    return embeddings


def add_transition_drift_to_obsm(
    adata: anndata.AnnData,
    transition_matrix: np.ndarray | sparse.spmatrix,
    *,
    basis: str = "umap",
    field_prefix: str = "T_fwd",
    dim: int | None = None,
    normalize_rows: bool,
    overwrite: bool = False,
) -> str:
    """Compute transition drift and store one sanctioned UMAP vector key."""
    basis = _require_nonempty_string(basis, label="basis")
    if basis != "umap":
        raise ValueError("basis must be exactly 'umap'.")
    field_prefix = _require_nonempty_string(
        field_prefix,
        label="field_prefix",
    )
    if type(overwrite) is not bool:
        raise TypeError("overwrite must be exactly True or False.")

    embeddings = _umap_embeddings(adata)
    if dim is None:
        if len(embeddings) != 1:
            raise ValueError("dim must be explicit when more than one UMAP embedding exists.")
        dimension = next(iter(embeddings))
    else:
        if isinstance(dim, bool) or not isinstance(dim, Integral):
            raise TypeError("dim must be an integer dimension.")
        dimension = int(dim)
        if dimension not in (1, 2, 3):
            raise ValueError("dim must be exactly 1, 2, or 3.")
        if dimension not in embeddings:
            raise ValueError(f"No matching {dimension}D UMAP embedding exists.")

    embedding = embeddings[dimension][1]
    drift = compute_transition_drift(
        transition_matrix,
        embedding,
        normalize_rows=normalize_rows,
    )

    field_id = f"{field_prefix}_umap"
    output_key = f"{field_id}_{dimension}d"
    if output_key in adata.obsm and not overwrite:
        raise KeyError(
            f"adata.obsm already contains key {output_key!r}; set overwrite=True to replace it."
        )

    validate_vector_fields(
        {output_key: drift},
        n_cells=int(adata.n_obs),
        available_dimensions=embeddings,
    )
    adata.obsm[output_key] = drift
    return output_key
