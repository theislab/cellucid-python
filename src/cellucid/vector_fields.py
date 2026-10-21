"""Vector-field validation and transition-drift utilities.

A vector field is a per-cell displacement drawn over one embedding. In an
export it is bound to embedding *ids*: ``vector_fields={"velocity": {"X_umap":
v2, "X_umap_3d": v3}}`` binds the field ``velocity`` to two items, and each
array is ``(n_cells, dimension of that item)``. The AnnData server builds the
same mapping from the ``<field>_umap_<d>d`` keys in ``obsm`` with the default
binding rule in :mod:`cellucid.data_mapping`.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import TYPE_CHECKING, cast

import numpy as np
from scipy import sparse

from .continuous_payload_diagnosis import describe_non_finite
from .data_mapping import (
    VECTOR_KEY_PATTERN,
    default_embedding_slots,
    describe_obsm,
    obsm_entries,
)

if TYPE_CHECKING:
    import anndata


@dataclass(frozen=True)
class ValidatedVectorFields:
    """Validated vector arrays, keyed field id -> embedding id.

    Fields are in code-point order and each field's embeddings follow the
    order of the ``embedding_dimensions`` mapping it was validated against, so
    the payload indices derived from both orders are the ones every writer
    uses. Arrays are float64; the one float32 rounding happens when they are
    scaled.
    """

    fields: dict[str, dict[str, np.ndarray]]
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
        raise ValueError(
            f"{label} must contain only finite values. "
            + describe_non_finite(array, label)
        )
    with np.errstate(over="ignore", invalid="ignore"):
        narrowed = array.astype(np.float32, copy=False)
        if not np.isfinite(narrowed).all():
            raise ValueError(
                f"{label} contains values outside the finite float32 range. "
                + describe_non_finite(array, label)
            )
    # Underflow is the other end of the same range, and this path was missing it
    # while every other float32 payload checked both ends: a nonzero component
    # below the smallest float32 subnormal becomes 0.0, which is finite, so an
    # overflow-only check published a vector whose component had silently become
    # zero. cellucid-r refuses exactly this in `.write_float32_matrix_row_major`,
    # and one export format cannot mean two things.
    if np.any((array != 0) & (narrowed == 0)):
        raise ValueError(
            f"{label} contains values outside the finite float32 range. "
            + describe_non_finite(array, label)
        )
    return array.astype(np.float64, copy=False), dimension


def validate_vector_fields(
    values: object,
    *,
    n_cells: int,
    embedding_dimensions: Mapping[str, int],
    vector_field_default: object = None,
    label_for: Callable[[str, str], str] | None = None,
) -> ValidatedVectorFields:
    """Validate every field's bindings against the embeddings they name.

    ``values`` maps a field id to ``{embedding id: array}``. Each array must be
    ``(n_cells, d)`` for that embedding's dimension ``d`` (a plain vector is
    accepted for a one-dimensional embedding) and hold finite values that are
    representable as float32. ``vector_field_default`` is required when more
    than one field is declared.
    """
    if values is not None and (isinstance(values, str | bytes) or not isinstance(values, Mapping)):
        raise TypeError(
            "vector_fields must map each field id to a mapping of embedding ids to arrays, "
            "for example {'velocity': {'X_umap': velocity_2d}}."
        )

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
        return ValidatedVectorFields({}, None)

    order = {embedding_id: position for position, embedding_id in enumerate(embedding_dimensions)}
    fields: dict[str, dict[str, np.ndarray]] = {}
    for raw_field, raw_bindings in values.items():
        field_id = _require_nonempty_string(raw_field, label="vector field id")
        if isinstance(raw_bindings, str | bytes) or not isinstance(raw_bindings, Mapping):
            raise TypeError(
                f"Vector field {field_id!r} must map embedding ids to arrays, for example "
                f"{{{field_id!r}: {{'X_umap': vectors}}}}."
            )
        if not raw_bindings:
            raise ValueError(f"Vector field {field_id!r} binds no embedding.")
        bound: dict[str, np.ndarray] = {}
        for raw_embedding, raw_values in raw_bindings.items():
            embedding_id = _require_nonempty_string(
                raw_embedding,
                label=f"Vector field {field_id!r} embedding id",
            )
            if embedding_id not in embedding_dimensions:
                raise ValueError(
                    f"Vector field {field_id!r} is bound to embedding {embedding_id!r}, "
                    f"which is not one of the exported embeddings: {list(embedding_dimensions)}."
                )
            label = (
                label_for(field_id, embedding_id)
                if label_for is not None
                else f"Vector field {field_id!r} on embedding {embedding_id!r}"
            )
            array, _dimension = _finite_float32_matrix(
                raw_values,
                label=label,
                n_rows=n_cells,
                declared_dimension=embedding_dimensions[embedding_id],
            )
            bound[embedding_id] = array
        fields[field_id] = {
            embedding_id: bound[embedding_id]
            for embedding_id in sorted(bound, key=lambda embedding_id: order[embedding_id])
        }

    ordered_fields = {field_id: fields[field_id] for field_id in sorted(fields)}
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

    return ValidatedVectorFields(fields=ordered_fields, default_field=default_field)


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
        raise ValueError(
            f"{label} scaling produced non-finite values. "
            + describe_non_finite(scaled, f"{label} scaled by {numeric_scale:g}")
        )
    with np.errstate(over="ignore", invalid="ignore"):
        float32_scaled = scaled.astype(np.float32)
    if not np.isfinite(float32_scaled).all():
        raise ValueError(
            f"{label} scaling produced values outside the finite float32 range. "
            + describe_non_finite(scaled, f"{label} scaled by {numeric_scale:g}")
        )
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
            raise ValueError(
                "transition_matrix must contain only finite values. "
                + describe_non_finite(
                    matrix.data, "transition_matrix", positions=False
                )
            )
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
        raise ValueError(
            "transition_matrix must contain only finite values. "
            + describe_non_finite(matrix, "transition_matrix")
        )
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
    """Return the UMAP embedding of every dimension the default rule resolves."""
    n_cells = int(adata.n_obs)
    slots = default_embedding_slots(obsm_entries(describe_obsm(adata.obsm, n_obs=n_cells)))
    embeddings: dict[int, tuple[str, np.ndarray]] = {}
    for slot, ref in slots.items():
        if ref is None:
            continue
        dimension = int(slot[0])
        source = np.asarray(
            adata.obsm[ref["key"]].toarray()
            if sparse.issparse(adata.obsm[ref["key"]])
            else adata.obsm[ref["key"]]
        )
        if source.ndim == 1:
            source = source.reshape(-1, 1)
        array, _resolved_dimension = _finite_float32_matrix(
            source[:, ref["columns"]] if source.ndim == 2 else source,
            label=f"Embedding {ref['key']!r}",
            n_rows=n_cells,
            declared_dimension=dimension,
        )
        embeddings[dimension] = (ref["key"], array)

    if not embeddings:
        raise ValueError(
            "AnnData must contain a UMAP embedding Cellucid can read: obsm['X_umap_1d'], "
            "obsm['X_umap_2d'], obsm['X_umap_3d'], or a 1- to 3-column obsm['X_umap']. "
            f"obsm holds: {list(adata.obsm.keys())}."
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

    # The key is the one the default binding rule reads back: it binds to the
    # UMAP embedding of the same dimension when the dataset is served.
    if VECTOR_KEY_PATTERN.fullmatch(output_key) is None:
        raise ValueError(f"{output_key!r} is not a '<field>_umap_<1|2|3>d' vector key.")
    adata.obsm[output_key] = drift
    return output_key
