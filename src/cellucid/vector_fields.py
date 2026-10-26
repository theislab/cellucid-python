"""Transition-drift vector fields, computed into ``obsm``.

A vector field is a per-cell displacement drawn over one embedding. In an
AnnData it is the ``obsm`` array ``<field>_<basis>``, which the viewer draws
as the field ``<field>`` over the embedding ``X_<basis>`` -- scVelo's naming
convention (``velocity_umap`` over ``X_umap``, ``T_fwd_branches`` over
``X_branches``). These utilities compute such a field from a cell-cell
transition matrix (CellRank's, for example) and store it under that key.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
from scipy import sparse

from ._non_finite import describe_non_finite
from ._plural import plural
from .data_mapping import EMBEDDING_KEY_PREFIX, parse_vector_field_key

if TYPE_CHECKING:
    import anndata


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
    """Validate one matrix as storable in float32 and return it in float64.

    The rounding to float32 belongs to the one place that writes the values,
    so the array keeps the caller's precision until then.
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
        raise ValueError(f"{label} must have exactly {n_rows} {plural(n_rows, 'row')}, got {array.shape[0]}.")
    if declared_dimension is None:
        dimension = int(array.shape[1])
        if dimension < 1:
            raise ValueError(f"{label} must have at least one column.")
    else:
        dimension = declared_dimension
        if array.shape[1] != dimension:
            raise ValueError(
                f"{label} must have exactly {dimension} {plural(dimension, 'column')}, got {array.shape[1]}."
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
    # Underflow is the other end of the same range: a nonzero component below
    # the smallest float32 subnormal becomes 0.0, which is finite, so an
    # overflow-only check would store a vector whose component had silently
    # become zero.
    if np.any((array != 0) & (narrowed == 0)):
        raise ValueError(
            f"{label} contains values outside the finite float32 range. "
            + describe_non_finite(array, label)
        )
    return array.astype(np.float64, copy=False), dimension


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

    The drift of cell ``i`` is the transition-weighted mean position of its
    successors minus its own position: ``T @ E - E``, with each row of ``T``
    first divided by its sum when ``normalize_rows`` is True.

    Parameters
    ----------
    transition_matrix : numpy.ndarray or scipy.sparse matrix
        Non-negative, finite ``(n_cells, n_cells)`` cell-cell transition
        matrix (CellRank's, for example).
    embedding : numpy.ndarray
        Finite ``(n_cells, n_dims)`` coordinates the drift is computed in; a
        1D array is not accepted. Values must be representable in float32.
    normalize_rows : bool
        Required. True divides every row of ``transition_matrix`` by its sum
        (every row sum must then be nonzero); False uses it as given.

    Returns
    -------
    numpy.ndarray
        ``(n_cells, n_dims)`` float64 drift, at the precision it was computed
        in; it is rounded to float32 once, when it is written.

    Raises
    ------
    TypeError
        A non-bool ``normalize_rows``, or non-numeric input.
    ValueError
        A shape mismatch, a negative or non-finite value, a value outside the
        float32 range, or a zero row sum with ``normalize_rows=True``.
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


def add_transition_drift_to_obsm(
    adata: anndata.AnnData,
    transition_matrix: np.ndarray | sparse.spmatrix,
    *,
    basis: str = "umap",
    field_prefix: str = "T_fwd",
    normalize_rows: bool,
    overwrite: bool = False,
) -> str:
    """Compute transition drift over ``obsm['X_<basis>']`` and store it as ``<field_prefix>_<basis>``.

    The drift is computed in the coordinates of the embedding ``X_<basis>``,
    every column of it, and stored as ``obsm['<field_prefix>_<basis>']``, as
    wide as the embedding: the key the viewer binds to that embedding by
    default (``T_fwd_umap`` over ``X_umap``, ``T_fwd_umap_3d`` over
    ``X_umap_3d``, ``T_fwd_branches`` over ``X_branches``).

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData holding ``obsm['X_<basis>']``; changed in place.
    transition_matrix : numpy.ndarray or scipy.sparse matrix
        Non-negative, finite ``(n_obs, n_obs)`` cell-cell transition matrix.
    basis : str
        The embedding ``obsm['X_<basis>']`` the drift is computed in.
    field_prefix : str
        The vector field name: the drift is stored as
        ``obsm['<field_prefix>_<basis>']``.
    normalize_rows : bool
        Required; see :func:`compute_transition_drift`.
    overwrite : bool
        Replace an existing ``obsm['<field_prefix>_<basis>']``.

    Returns
    -------
    str
        The ``obsm`` key written.

    Raises
    ------
    KeyError
        ``obsm['<field_prefix>_<basis>']`` exists and ``overwrite`` is False.
    ValueError
        ``obsm['X_<basis>']`` is missing, the key would not read back as the
        field ``field_prefix`` on that embedding, or the inputs are invalid
        (see :func:`compute_transition_drift`).
    """
    basis = _require_nonempty_string(basis, label="basis")
    field_prefix = _require_nonempty_string(field_prefix, label="field_prefix")
    if type(overwrite) is not bool:
        raise TypeError("overwrite must be exactly True or False.")

    embedding_key = EMBEDDING_KEY_PREFIX + basis
    if embedding_key not in adata.obsm:
        raise ValueError(
            f"basis={basis!r} names obsm[{embedding_key!r}], which this AnnData does not have. "
            f"obsm holds: {list(adata.obsm.keys())}."
        )
    output_key = f"{field_prefix}_{basis}"
    # The key must read back as this field on this embedding: the viewer
    # binds it by the vector-field key rule, which a field named like an
    # embedding, or another embedding with a longer matching basis, defeats.
    parsed = parse_vector_field_key(output_key, {*adata.obsm.keys(), output_key})
    if parsed is None or parsed.field != field_prefix or parsed.basis != basis:
        raise ValueError(
            f"obsm[{output_key!r}] would not read as the field {field_prefix!r} on "
            f"obsm[{embedding_key!r}]: "
            + (
                "a key that starts with 'X_' names an embedding, never a vector field."
                if parsed is None
                else f"it reads as the field {parsed.field!r} on obsm[{parsed.embedding_key!r}]."
            )
            + " Choose another field_prefix."
        )
    if output_key in adata.obsm and not overwrite:
        raise KeyError(
            f"adata.obsm already contains key {output_key!r}; set overwrite=True to replace it."
        )

    values = adata.obsm[embedding_key]
    source = np.asarray(
        cast(sparse.spmatrix, values).toarray() if sparse.issparse(values) else values
    )
    if source.ndim == 1:
        source = source.reshape(-1, 1)
    embedding, _dimension = _finite_float32_matrix(
        source,
        label=f"Embedding {embedding_key!r}",
        n_rows=int(adata.n_obs),
        declared_dimension=None,
    )
    drift = compute_transition_drift(
        transition_matrix,
        embedding,
        normalize_rows=normalize_rows,
    )
    adata.obsm[output_key] = drift
    return output_key
