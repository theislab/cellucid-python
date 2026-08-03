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

#: ``sc.tl.umap()`` writes this key, without a dimension in its name.
UNSUFFIXED_EMBEDDING_KEY: str = "X_umap"


def is_vector_field_declaration_key(value: object) -> bool:
    """Return whether an AnnData ``obsm`` key exactly declares a vector field."""
    return isinstance(value, str) and _VECTOR_KEY_PATTERN.fullmatch(value) is not None


def _servable_column_count(array: object) -> int | None:
    """Return an array's column count when Cellucid can serve that width."""
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple):
        return None
    # AnnData stores a plain ``(n_cells,)`` vector as the single column it is.
    if len(shape) == 1:
        return 1
    if len(shape) != 2:
        return None
    columns = int(shape[1])
    return columns if columns in (1, 2, 3) else None


def resolve_embedding_keys(obsm: Any) -> dict[int, str]:
    """Return every servable ``{dimension: obsm key}`` this object declares.

    A key that names its own dimension -- ``X_umap_1d``, ``X_umap_2d``,
    ``X_umap_3d`` -- always decides that dimension. When an object declares
    none of them, ``sc.tl.umap()``'s own ``X_umap`` resolves to the dimension
    its column count states, so the ordinary Scanpy object is servable as
    written. A dimensional key anywhere in ``obsm`` means the writer declared
    dimensions deliberately, and ``X_umap`` is then left to the vector-field
    and hint paths rather than silently joining a declared set.
    """
    declared = {
        dimension: key
        for dimension, key in ((dimension, f"X_umap_{dimension}d") for dimension in (1, 2, 3))
        if key in obsm
    }
    if declared:
        return declared
    if UNSUFFIXED_EMBEDDING_KEY not in obsm:
        return {}
    columns = _servable_column_count(obsm[UNSUFFIXED_EMBEDDING_KEY])
    if columns is None:
        return {}
    return {columns: UNSUFFIXED_EMBEDDING_KEY}


def unservable_unsuffixed_embedding_shape(obsm: Any) -> tuple[int, ...] | None:
    """Return ``X_umap``'s shape when its width is the reason nothing resolved.

    A ten-column ``X_umap`` is a latent space someone named after a plot. It
    produces no rename hint, because no rename makes ten columns servable, so
    the error message has to name the shape itself.
    """
    if UNSUFFIXED_EMBEDDING_KEY not in obsm:
        return None
    if _servable_column_count(obsm[UNSUFFIXED_EMBEDDING_KEY]) is not None:
        return None
    shape = getattr(obsm[UNSUFFIXED_EMBEDDING_KEY], "shape", None)
    if not isinstance(shape, tuple):
        return None
    return tuple(int(extent) for extent in shape)


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


def missing_embedding_message(obsm: Any, *, n_cells: int, headline: str) -> str:
    """Return one message naming the rule, the object, and the exact remedy."""
    message = (
        f"{headline} Cellucid reads the exact keys 'X_umap_1d', 'X_umap_2d', and "
        "'X_umap_3d', and reads a bare 'X_umap' as the dimension its own column "
        "count states. "
        f"Available obsm keys: {list(obsm.keys())}."
    )
    unservable_shape = unservable_unsuffixed_embedding_shape(obsm)
    if unservable_shape is not None:
        message += (
            f"\n  adata.obsm['X_umap'] has shape {unservable_shape}, and Cellucid "
            "draws 1, 2, or 3 dimensions, so that array names a dimension no "
            "viewer renders. Assign the columns you mean to draw under the key "
            "for their own count."
        )
    hints = embedding_declaration_hints(obsm, n_cells=n_cells)
    if hints:
        joined = "\n".join(f"    {hint}" for hint in hints)
        message += (
            "\n  Fix: declare an existing array under the key for its column "
            f"count, then re-run (re-save the file if you serve a path):\n{joined}"
        )
    return message


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
    """Validate one matrix as float32-publishable and return it in float64.

    The payload is float32, but the rounding to float32 belongs to the one
    place that writes the bytes. Rounding here as well would round twice --
    once on the input and once on the scaled result -- and every component
    whose two roundings disagree lands one ULP away from the correctly rounded
    value. cellucid-r scales the caller's doubles and rounds once, so a second
    rounding here is also what made the two writers disagree on vectors/*.
    """
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
        if not np.isfinite(array.astype(np.float32, copy=False)).all():
            raise ValueError(f"{label} contains values outside the finite float32 range.")
    return array.astype(np.float64, copy=False), dimension


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
    """Scale vectors using one required finite positive scale.

    This is the only float32 rounding a vector payload gets: the caller's own
    values are scaled at their own precision and the product is rounded once,
    exactly as cellucid-r does.
    """
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
    """Compute exact per-cell transition drift in embedding space.

    The drift is returned in float64 because it is an input to the export, not
    a payload: rounding it here would round once more than the export does.
    """
    if type(normalize_rows) is not bool:
        raise TypeError("normalize_rows must be exactly True or False.")
    working_embedding, _dimension = _finite_float32_matrix(
        embedding,
        label="embedding",
        n_rows=None,
        declared_dimension=None,
    )
    n_cells = int(working_embedding.shape[0])
    matrix = _transition_matrix(transition_matrix, n_cells=n_cells)
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
        declared_dimension=int(working_embedding.shape[1]),
    )
    return drift_array


def _umap_embeddings(
    adata: anndata.AnnData,
) -> dict[int, tuple[str, np.ndarray]]:
    embeddings: dict[int, tuple[str, np.ndarray]] = {}
    n_cells = int(adata.n_obs)
    for dimension, key in sorted(resolve_embedding_keys(adata.obsm).items()):
        array, _resolved_dimension = _finite_float32_matrix(
            adata.obsm[key],
            label=f"Embedding {key!r}",
            n_rows=n_cells,
            declared_dimension=dimension,
        )
        embeddings[dimension] = (key, array)

    if not embeddings:
        raise ValueError(
            missing_embedding_message(
                adata.obsm,
                n_cells=n_cells,
                headline="AnnData must contain one UMAP embedding Cellucid can read.",
            )
        )
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
