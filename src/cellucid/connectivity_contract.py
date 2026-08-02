"""Exact connectivity validation and canonical edge-pair construction."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Any

import numpy as np
from scipy import sparse

from ._byte_order import little_endian_payload_dtype

CONNECTIVITY_MANIFEST_FILENAME = "connectivity_manifest.json"
CONNECTIVITY_BINARY_DIRNAME = "connectivity"
CONNECTIVITY_SOURCES_PATH = "connectivity/edges.src.bin"
CONNECTIVITY_DESTINATIONS_PATH = "connectivity/edges.dst.bin"
CONNECTIVITY_WEIGHTS_PATH = "connectivity/edges.weights.f64.bin"


@dataclass(frozen=True)
class ConnectivityEdgePairs:
    """One validated, canonical connectivity payload."""

    sources: np.ndarray
    destinations: np.ndarray
    weights: np.ndarray
    n_edges: int
    max_neighbors: int
    index_dtype: str
    index_bytes: int


def connectivity_index_storage(n_cells: int) -> tuple[type[np.unsignedinteger], str, int]:
    """Return the exact browser-supported index storage for one cell axis."""
    if not isinstance(n_cells, Integral) or isinstance(n_cells, bool | np.bool_) or n_cells <= 0:
        raise ValueError("n_cells must be a positive integer.")
    if n_cells <= 65_536:
        return np.uint16, "uint16", 2
    if n_cells <= 4_294_967_296:
        return np.uint32, "uint32", 4
    raise ValueError(
        "Connectivity cannot exceed 4,294,967,296 cells because the current "
        "browser contract supports at most the complete uint32 index domain."
    )


def build_connectivity_manifest(
    *,
    n_cells: int,
    n_edges: int,
    max_neighbors: int,
    index_dtype: str,
    index_bytes: int,
    compression: int | None,
) -> dict[str, int | str | None]:
    """Build the sole canonical weighted edge manifest."""
    _, expected_index_dtype, expected_index_bytes = connectivity_index_storage(n_cells)
    if index_dtype != expected_index_dtype or index_bytes != expected_index_bytes:
        raise ValueError("Connectivity index storage does not match the complete cell axis.")
    if type(n_edges) is not int or n_edges < 0:
        raise ValueError("n_edges must be a non-negative integer.")
    if n_edges > n_cells * (n_cells - 1) // 2:
        raise ValueError("n_edges exceeds the complete undirected cell graph.")
    if (
        type(max_neighbors) is not int
        or max_neighbors < 0
        or max_neighbors >= n_cells
        or (n_edges == 0 and max_neighbors != 0)
        or (n_edges > 0 and max_neighbors == 0)
    ):
        raise ValueError("max_neighbors does not match a valid undirected graph.")
    if compression is not None and (type(compression) is not int or not 1 <= compression <= 9):
        raise ValueError("compression must be None or one integer from 1 to 9.")

    suffix = ".gz" if compression is not None else ""
    return {
        "format": "edge_pairs",
        "n_cells": int(n_cells),
        "n_edges": n_edges,
        "max_neighbors": max_neighbors,
        "index_bytes": index_bytes,
        "index_dtype": index_dtype,
        "sourcesPath": f"{CONNECTIVITY_SOURCES_PATH}{suffix}",
        "destinationsPath": f"{CONNECTIVITY_DESTINATIONS_PATH}{suffix}",
        "weightsPath": f"{CONNECTIVITY_WEIGHTS_PATH}{suffix}",
        "weight_dtype": "float64",
        "weight_bytes": 8,
        "compression": compression,
    }


def _validate_numeric_dtype(dtype: np.dtype[Any]) -> None:
    if np.issubdtype(dtype, np.complexfloating):
        raise ValueError("Connectivity values must be real.")
    if not (np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.bool_)):
        raise ValueError("Connectivity values must be numeric or boolean.")


def _validated_float64_weights(
    values: np.ndarray,
    *,
    reject_zero: bool,
) -> np.ndarray:
    """Validate weights and return their exact little-endian Float64 form."""
    _validate_numeric_dtype(values.dtype)
    if not np.all(np.isfinite(values)):
        raise ValueError("Connectivity values must all be finite.")
    if np.any(values < 0):
        raise ValueError("Connectivity values must all be non-negative.")
    if reject_zero and np.any(values == 0):
        raise ValueError("Sparse connectivity must not contain stored zero entries.")

    with np.errstate(over="ignore", invalid="ignore"):
        weights = values.astype(little_endian_payload_dtype(np.float64), copy=False)
    if not np.all(np.isfinite(weights)):
        raise ValueError("Connectivity values must remain finite when represented as Float64.")

    if np.issubdtype(values.dtype, np.integer):
        large_mask = values > 9_007_199_254_740_992
        if np.any(large_mask):
            original_large = values[large_mask].reshape(-1)
            float_large = weights[large_mask].reshape(-1)
            if any(
                int(original) != int(converted)
                for original, converted in zip(
                    original_large,
                    float_large,
                    strict=True,
                )
            ):
                raise ValueError("Connectivity values must be exactly representable as Float64.")
    elif np.issubdtype(values.dtype, np.floating) and values.dtype.itemsize > 8:
        if not np.array_equal(
            weights.astype(values.dtype),
            values,
        ):
            raise ValueError("Connectivity values must be exactly representable as Float64.")

    return weights


def _validate_unique_sparse_coordinates(connectivities: sparse.spmatrix) -> None:
    coordinates = connectivities.tocoo(copy=False)
    rows = np.asarray(coordinates.row)
    columns = np.asarray(coordinates.col)
    if rows.size < 2:
        return

    column_count = np.uint64(connectivities.shape[1])
    coordinate_keys = rows.astype(np.uint64, copy=True)
    coordinate_keys *= column_count
    coordinate_keys += columns.astype(np.uint64, copy=False)
    coordinate_keys.sort()
    duplicate_offsets = np.flatnonzero(coordinate_keys[1:] == coordinate_keys[:-1])
    if duplicate_offsets.size:
        duplicate_offset = int(duplicate_offsets[0])
        duplicate_key = coordinate_keys[duplicate_offset]
        raise ValueError(
            "Connectivity contains duplicate sparse coordinates at "
            f"({int(duplicate_key // column_count)}, "
            f"{int(duplicate_key % column_count)})."
        )


def _sparse_upper_pairs(
    connectivities: sparse.spmatrix,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _validate_unique_sparse_coordinates(connectivities)
    coordinates = connectivities.tocoo(copy=False)
    rows = np.asarray(coordinates.row, dtype=np.int64)
    columns = np.asarray(coordinates.col, dtype=np.int64)
    weights = _validated_float64_weights(
        np.asarray(coordinates.data),
        reject_zero=True,
    )

    if np.any(rows == columns):
        raise ValueError("Connectivity diagonal values must all be exactly zero.")

    matrix = sparse.coo_matrix(
        (weights, (rows, columns)),
        shape=connectivities.shape,
    ).tocsr()
    matrix.sort_indices()
    if (matrix != matrix.transpose()).nnz != 0:
        raise ValueError("Connectivity matrix must be exactly symmetric in topology and weight.")

    upper_mask = rows < columns
    return (
        rows[upper_mask],
        columns[upper_mask],
        weights[upper_mask],
    )


def _dense_upper_pairs(
    connectivities: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = np.asarray(connectivities)
    if matrix.ndim != 2:
        raise ValueError(f"Connectivity matrix must be two-dimensional; got shape {matrix.shape}.")
    weights_matrix = _validated_float64_weights(
        matrix,
        reject_zero=False,
    )

    if np.any(np.diag(weights_matrix) != 0):
        raise ValueError("Connectivity diagonal values must all be exactly zero.")
    if not np.array_equal(weights_matrix, weights_matrix.transpose()):
        raise ValueError("Connectivity matrix must be exactly symmetric in topology and weight.")

    sources, destinations = np.nonzero(np.triu(weights_matrix, k=1))
    return (
        np.asarray(sources, dtype=np.int64),
        np.asarray(destinations, dtype=np.int64),
        weights_matrix[sources, destinations],
    )


def validate_connectivity_edges(
    connectivities: Any,
    *,
    n_cells: int,
) -> ConnectivityEdgePairs:
    """Validate one exact graph and construct its canonical edge sequence."""
    dtype, dtype_name, index_bytes = connectivity_index_storage(n_cells)

    shape = getattr(connectivities, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 2 or shape != (n_cells, n_cells):
        raise ValueError(
            "Connectivity matrix shape must exactly match the cell axis: "
            f"expected {(n_cells, n_cells)}, got {shape}."
        )

    if sparse.issparse(connectivities):
        sources, destinations, weights = _sparse_upper_pairs(connectivities)
    else:
        sources, destinations, weights = _dense_upper_pairs(connectivities)

    if sources.size:
        order = np.lexsort((destinations, sources))
        sources = sources[order]
        destinations = destinations[order]
        weights = weights[order]

    if (
        sources.shape != destinations.shape
        or sources.shape != weights.shape
        or np.any(weights <= 0)
    ):
        raise RuntimeError("Validated connectivity edge payloads are not aligned.")

    degrees: np.ndarray = np.zeros(n_cells, dtype=np.uint64)
    if sources.size:
        np.add.at(degrees, sources, 1)
        np.add.at(degrees, destinations, 1)
    max_neighbors = int(degrees.max()) if degrees.size else 0

    little_endian_index_dtype = little_endian_payload_dtype(dtype)
    stored_sources: np.ndarray = sources.astype(
        little_endian_index_dtype,
        copy=False,
    )
    stored_destinations: np.ndarray = destinations.astype(
        little_endian_index_dtype,
        copy=False,
    )
    stored_weights = weights.astype(little_endian_payload_dtype(np.float64), copy=False)
    return ConnectivityEdgePairs(
        sources=stored_sources,
        destinations=stored_destinations,
        weights=stored_weights,
        n_edges=int(stored_sources.size),
        max_neighbors=max_neighbors,
        index_dtype=dtype_name,
        index_bytes=index_bytes,
    )
