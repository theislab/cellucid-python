"""
Checking what :func:`cellucid.prepare` was given and turning it into elements.

Everything here runs before anything is written -- but a MuData modality's
matrix, which is checked as it is converted, one modality at a time, while
the output is written beside its target under a temporary name (a failure
then removes it, and nothing reaches the target path). The rules are cellucid-r's
(``R/inputs.R``, ``R/anndata-elements.R``), so the two writers accept the same
inputs and produce the same AnnData.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any, cast
from urllib.parse import urlparse

import numpy as np
import pandas as pd
from scipy import sparse

from .._console import console_print
from .._contracts import _require_display_text
from .._non_finite import describe_non_finite
from .._plural import plural
from ..data_mapping import (
    EMBEDDING_KEY_PREFIX,
    column_gene_names,
    make_unique_names,
    parse_vector_field_key,
    unique_names_note,
)
from ._store import ZARR_METADATA_NAMES

#: The largest integer float32 holds exactly. Integer counts above it are
#: stored as float64, so every count reads back as it was given.
FLOAT32_EXACT_INTEGER_LIMIT = 2**24

_FLOAT32_MAX = float(np.finfo(np.float32).max)
#: Values tested at a time for being whole numbers (:func:`_exact_dtype`).
_INTEGRAL_BLOCK = 1 << 20

# Characters that cannot appear in a file name on every platform, plus the
# path separators and control characters.
_FORBIDDEN_KEY_CHARACTERS = re.compile(r'[/\\<>:"|?*\x00-\x1f]')
_WINDOWS_DEVICE_NAMES = frozenset(
    {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)), *(f"LPT{i}" for i in range(1, 10))}
)


def _plural(count: int, singular: str, plural: str) -> str:
    return singular if count == 1 else plural


# =============================================================================
# Element keys
# =============================================================================


def require_element_key(key: object, what: str) -> str:
    """One name that becomes a path in the store: an obs or var column, an obsm key.

    In a Zarr store each one is a folder on disk, or an entry name in a
    ``.zarr.zip`` someone may unpack, so it obeys the rules a portable path
    component obeys on every platform. The same rule holds for an ``.h5ad``,
    so the dataset can be converted to Zarr wherever it is copied.
    """
    if not isinstance(key, str) or not key:
        raise TypeError(f"{what} names must be non-empty strings; got {key!r}.")
    problem = None
    if _FORBIDDEN_KEY_CHARACTERS.search(key):
        problem = (
            "contains a character that cannot appear in a file name on every "
            'platform (/ \\ < > : " | ? * or a control character)'
        )
    elif key in (".", "..") or key in ZARR_METADATA_NAMES:
        problem = "is reserved by the Zarr store layout"
    elif key.endswith((" ", ".")):
        problem = "ends with a space or a dot, which Windows removes from file names"
    elif key.split(".", 1)[0].upper() in _WINDOWS_DEVICE_NAMES:
        problem = "is a reserved device name on Windows"
    if problem is not None:
        raise ValueError(
            f"{what} {key!r} {problem}. It becomes a path in the store, so rename it."
        )
    return key


def require_distinct_element_keys(keys: Sequence[str], what: str) -> None:
    """Names side by side in one group must differ by more than letter case.

    On macOS and Windows ``Sex`` and ``sex`` are the same folder, and one
    would overwrite the other.
    """
    seen: dict[str, str] = {}
    for key in keys:
        folded = key.casefold()
        if folded in seen:
            first = seen[folded]
            if first == key:
                raise ValueError(f"{what} name {key!r} is used twice.")
            raise ValueError(
                f"{what} names {first!r} and {key!r} differ only in letter case. They "
                "become folders of one Zarr group, and on macOS and Windows those are "
                "the same folder, so rename one of them."
            )
        seen[folded] = key


def _require_frame(frame: object, what: str, *, rows: int, rows_of: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{what} must be a pandas DataFrame, got {type(frame).__name__}.")
    if len(frame) != rows:
        raise ValueError(
            f"{what} has {len(frame):,} {_plural(len(frame), 'row', 'rows')}, but {rows_of} "
            f"{_plural(rows, 'has', 'have')} {rows:,}."
        )
    columns = list(frame.columns)
    for column in columns:
        require_element_key(column, f"{what} column")
    if "_index" in columns:
        raise ValueError(
            f"{what} column '_index' is reserved: AnnData stores the row names under "
            "that name. Rename the column."
        )
    require_distinct_element_keys(columns, f"{what} column")
    return frame


# =============================================================================
# Numbers
# =============================================================================


def _require_real(array: np.ndarray, label: str) -> None:
    if array.dtype.kind not in {"b", "i", "u", "f"}:
        raise TypeError(f"{label} must hold real numbers, got dtype {array.dtype}.")


def _require_finite_float32(values: np.ndarray, label: str, *, positions: bool = True) -> None:
    """Every value finite and within the float32 range.

    Decided from the smallest and largest value alone -- a NaN makes both NaN,
    an infinity is one of them -- so the check allocates nothing beside the
    values (a whole-matrix temporary would double the memory a large matrix
    needs).
    """
    if values.size == 0:
        return
    if values.dtype.kind == "b":
        within = True
    else:
        low, high = float(values.min()), float(values.max())
        within = bool(np.isfinite(low) and np.isfinite(high) and max(-low, high) <= _FLOAT32_MAX)
    if not within:
        raise ValueError(
            f"{label} must hold finite numbers within the float32 range. "
            + describe_non_finite(values, label, positions=positions)
        )


def cell_matrix(values: object, label: str) -> np.ndarray:
    """One per-cell array of real numbers as float32, two-dimensional.

    A plain ``(n_cells,)`` vector is one column. Sparse input is densified.
    """
    if values is None:
        raise ValueError(f"{label} has no values.")
    if sparse.issparse(values):
        values = cast(Any, values).toarray()
    array = np.asarray(values)
    _require_real(array, label)
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"{label} must be a 2-D array or a vector, got shape {array.shape}.")
    _require_finite_float32(array, label)
    return np.ascontiguousarray(array, dtype=np.float32)


# =============================================================================
# obsm: the embeddings and the vector fields
# =============================================================================


def resolve_embeddings(embeddings: object) -> tuple[dict[str, np.ndarray], int]:
    """Every embedding becomes ``obsm[key]``, whole: every column, as given.

    The viewer picks the columns it plots and normalizes them itself.
    """
    if isinstance(embeddings, str | bytes) or not isinstance(embeddings, Mapping):
        raise TypeError(
            "embeddings must map each obsm key to an array, for example "
            "embeddings={'X_umap': adata.obsm['X_umap']}."
        )
    if not embeddings:
        raise ValueError("embeddings must contain at least one embedding.")
    keys = list(embeddings)
    for key in keys:
        require_element_key(key, "Embedding")
    require_distinct_element_keys(keys, "Embedding")
    matrices: dict[str, np.ndarray] = {}
    n_cells: int | None = None
    first = ""
    for key in keys:
        matrix = cell_matrix(embeddings[key], f"Embedding {key!r}")
        if matrix.shape[1] < 1:
            raise ValueError(f"Embedding {key!r} has no columns.")
        if n_cells is None:
            n_cells, first = int(matrix.shape[0]), key
        elif matrix.shape[0] != n_cells:
            raise ValueError(
                f"All embeddings must have one row per cell: {first!r} has {n_cells:,} {plural(n_cells, 'row')} "
                f"and {key!r} has {matrix.shape[0]:,}."
            )
        matrices[key] = matrix
    assert n_cells is not None
    if n_cells < 1:
        raise ValueError("Embeddings must contain at least one cell.")
    return matrices, n_cells


def resolve_vector_fields(
    vector_fields: object,
    matrices: Mapping[str, np.ndarray],
    n_cells: int,
) -> tuple[dict[str, np.ndarray], list[str]]:
    """Each vector field becomes ``obsm['<field>_<basis>']``, scVelo's name for it.

    A field bound to the embedding ``X_<basis>`` is stored as
    ``obsm['<field>_<basis>']`` (``velocity`` on ``X_umap`` is
    ``velocity_umap``, on ``X_branches`` ``velocity_branches``), which the
    viewer draws as the field ``<field>`` over that embedding: the vector-field
    key rule (:func:`cellucid.data_mapping.parse_vector_field_key`). Each array
    has one column per column of its embedding. Returns the arrays by obsm key,
    and the field names.
    """
    if vector_fields is None:
        return {}, []
    if isinstance(vector_fields, str | bytes) or not isinstance(vector_fields, Mapping):
        raise TypeError(
            "vector_fields must map each field to a mapping of embeddings to arrays, "
            "for example vector_fields={'velocity': {'X_umap': velocity_umap}}."
        )
    arrays: dict[str, np.ndarray] = {}
    named: dict[str, tuple[str, str]] = {}
    for field, bindings in vector_fields.items():
        if not isinstance(field, str) or not field:
            raise TypeError(f"Vector field names must be non-empty strings; got {field!r}.")
        if isinstance(bindings, str | bytes) or not isinstance(bindings, Mapping) or not bindings:
            raise TypeError(
                f"Vector field {field!r} must map embeddings to arrays, for example "
                f"{{{field!r}: {{'X_umap': vectors}}}}."
            )
        for embedding, values in bindings.items():
            label = f"Vector field {field!r} on embedding {embedding!r}"
            if embedding not in matrices:
                raise ValueError(
                    f"{label}: {embedding!r} is not one of the embeddings: {list(matrices)}."
                )
            if not embedding.startswith(EMBEDDING_KEY_PREFIX) or embedding == EMBEDDING_KEY_PREFIX:
                raise ValueError(
                    f"{label}: the viewer binds a vector field to the embedding 'X_<basis>' "
                    f"it is named for ('velocity_umap' to 'X_umap'), so the embedding needs "
                    f"an 'X_<basis>' key. Rename {embedding!r}, for example to "
                    f"'X_{embedding}'."
                )
            array = cell_matrix(values, label)
            width = int(matrices[embedding].shape[1])
            if array.shape[1] != width:
                raise ValueError(
                    f"{label} has {array.shape[1]} {plural(array.shape[1], 'column')}, and {embedding!r} has {width}: "
                    "a vector field has one component per embedding column."
                )
            if array.shape[0] != n_cells:
                raise ValueError(
                    f"{label} must have one row per cell ({n_cells:,}), got {array.shape[0]:,}."
                )
            basis = embedding[len(EMBEDDING_KEY_PREFIX) :]
            key = f"{field}_{basis}"
            if key in arrays:
                other_field, other_embedding = named[key]
                raise ValueError(
                    f"{label} and vector field {other_field!r} on {other_embedding!r} would "
                    f"both be stored as obsm[{key!r}]. Rename one field."
                )
            require_element_key(key, "Vector field obsm key")
            arrays[key] = array
            named[key] = (field, embedding)
    collisions = {key.casefold() for key in arrays} & {key.casefold() for key in matrices}
    if collisions:
        raise ValueError(
            f"A vector field would be stored as obsm[{sorted(collisions)[0]!r}], which is also "
            "the key of an embedding. Rename the field or the embedding."
        )
    obsm_keys = {*matrices, *arrays}
    for key, (field, embedding) in named.items():
        parsed = parse_vector_field_key(key, obsm_keys)
        if parsed is None or parsed.field != field or parsed.embedding_key != embedding:
            reads_as = (
                "an embedding: a key that starts with 'X_' is never a vector field"
                if parsed is None
                else f"the field {parsed.field!r} on {parsed.embedding_key!r}"
            )
            raise ValueError(
                f"Vector field {field!r} on embedding {embedding!r} would be stored as "
                f"obsm[{key!r}], which the viewer reads as {reads_as}. Rename the field."
            )
    return arrays, list(vector_fields)


# =============================================================================
# obs and var columns
# =============================================================================


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA:
        return True
    return isinstance(value, float) and value != value


def _string_column(values: np.ndarray, label: str) -> Any:
    """Text: a categorical when a value repeats or is missing, else a string array.

    This is anndata's own ``strings_to_categoricals`` rule. Categories are in
    code-point order; a missing value is code -1.
    """
    # Vectorized: pandas finds the candidates, and the few it flags are held
    # to this module's own rule (None, pandas NA, a float NaN); one pass per
    # value in Python cost half a second per million cells and text column.
    missing = np.asarray(pd.isna(values), dtype=bool)
    flagged = np.flatnonzero(missing)
    if flagged.size:
        missing[flagged] = [_is_missing(values[index]) for index in flagged]
    present = values[~missing]
    if present.size and pd.api.types.infer_dtype(present, skipna=False) != "string":
        value = next(item for item in present if not isinstance(item, str))
        raise TypeError(
            f"{label} mixes text with {type(value).__name__} values, which has no AnnData "
            "encoding. Convert it first, for example with .astype(str)."
        )
    distinct = set(present.tolist())
    if missing.any() or len(distinct) < len(present):
        labels = values.copy() if missing.any() else values
        if missing.any():
            labels[missing] = None
        return pd.Categorical(labels, categories=sorted(distinct), ordered=False)
    return np.asarray(values, dtype=object)


def _frame_column(series: pd.Series, label: str) -> Any:
    """How one data-frame column is stored, by anndata's own encodings."""
    dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        categorical = pd.Categorical(series)
        kind = categorical.categories.dtype.kind
        if kind == "O" and not all(isinstance(value, str) for value in categorical.categories):
            raise TypeError(
                f"{label} has categories that are not all text. Convert them first, for "
                "example with .cat.rename_categories(str)."
            )
        if kind not in {"O", "i", "u", "f", "b"}:
            raise TypeError(f"{label} has {categorical.categories.dtype} categories, which anndata cannot store.")
        return categorical
    if isinstance(dtype, pd.BooleanDtype | pd.core.arrays.integer.IntegerDtype):
        array = series.array
        if not series.isna().any():
            return series.to_numpy(dtype=dtype.numpy_dtype)
        return array
    if isinstance(dtype, pd.core.arrays.floating.FloatingDtype):
        series = series.astype(dtype.numpy_dtype)
        dtype = series.dtype
    if isinstance(dtype, pd.StringDtype) or (isinstance(dtype, np.dtype) and dtype.kind == "O"):
        return _string_column(series.to_numpy(dtype=object), label)
    if isinstance(dtype, np.dtype) and dtype.kind == "b":
        return series.to_numpy()
    if isinstance(dtype, np.dtype) and dtype.kind in {"i", "u"}:
        return series.to_numpy()
    if isinstance(dtype, np.dtype) and dtype.kind == "f":
        values = series.to_numpy()
        infinite = np.flatnonzero(np.isinf(values))
        if infinite.size:
            raise ValueError(
                f"{label} holds {infinite.size:,} infinite "
                f"{_plural(int(infinite.size), 'value', 'values')} (first at row {int(infinite[0]):,}). "
                "A colour scale has no position for an infinity; replace it with NaN, which is "
                "stored as a missing value."
            )
        return values
    raise TypeError(
        f"{label} is a {dtype} column, which has no AnnData encoding Cellucid writes. Convert "
        "it first, for example with .astype(str) or .astype(float)."
    )


def _frame_columns(frame: pd.DataFrame, what: str) -> dict[str, Any]:
    return {
        str(column): _frame_column(frame[column], f"{what} column {column!r}")
        for column in frame.columns
    }


def _row_names(frame: pd.DataFrame, what: str) -> list[str] | None:
    """The frame's own row names, or None for pandas' default ``0 .. n-1`` index."""
    index = frame.index
    if isinstance(index, pd.RangeIndex) and index.start == 0 and index.step == 1:
        return None
    names: list[str] = index.tolist()
    for position, name in enumerate(names):
        if not isinstance(name, str) or not name:
            raise TypeError(
                f"{what} row names must be non-empty text; row {position:,} is {name!r}. "
                f"Convert them first, for example with {what}.index = {what}.index.astype(str)."
            )
    return names


def resolve_obs(obs: object, n_cells: int) -> tuple[np.ndarray, dict[str, Any]]:
    """The obs index (the cell names) and the stored obs columns.

    The cells are named by ``obs.index``, or ``"0"`` to ``"n - 1"`` (as
    anndata names them) when obs is absent or has pandas' default index.
    Repeated names are made unique as anndata's ``obs_names_make_unique``
    makes them, and the change is reported.
    """
    if obs is None:
        return np.asarray([str(i) for i in range(n_cells)], dtype=object), {}
    frame = _require_frame(obs, "obs", rows=n_cells, rows_of="the embeddings")
    names = _row_names(frame, "obs")
    if names is None:
        names = [str(i) for i in range(n_cells)]
    unique = make_unique_names(names)
    if unique.renamed:
        note = unique_names_note(filled_from_index=0, renamed=unique.renamed, source=None)
        console_print(f"Cell names: {note}")
    return np.asarray(unique.names, dtype=object), _frame_columns(frame, "obs")


# =============================================================================
# X and var
# =============================================================================


def _exact_dtype(values: np.ndarray) -> type[np.floating[Any]]:
    """float32, unless the values are integer counts float32 cannot hold exactly."""
    if values.size == 0:
        return np.float32
    if values.dtype.kind in {"i", "u"}:
        integral = True
    elif values.dtype.kind == "f" and values.dtype.itemsize > 4:
        # In blocks, so the test's temporaries stay small beside the values.
        flat = values.reshape(-1)
        integral = all(
            bool(np.all(np.floor(block) == block))
            for block in (flat[start : start + _INTEGRAL_BLOCK] for start in range(0, flat.size, _INTEGRAL_BLOCK))
        )
    else:
        return np.float32
    if integral and max(-float(values.min()), float(values.max())) > FLOAT32_EXACT_INTEGER_LIMIT:
        return np.float64
    return np.float32


def expression_dtype(values: np.ndarray) -> type[np.floating[Any]]:
    """The type an expression matrix's values are stored in: float32, or
    float64 for integer counts above 2^24 (what ``prepare()`` writes, and what
    a server presents for an object; :func:`stored_expression`)."""
    return _exact_dtype(values)


def require_expression_values(values: np.ndarray, label: str) -> None:
    """Expression values ``prepare()`` stores: real, finite, within the float32 range."""
    _require_real(values, label)
    _require_finite_float32(values, label, positions=False)


def _canonical_in_place(matrix: Any) -> None:
    """Sorted, unique indices and no stored zero, a duplicate coordinate summed
    in stored order (scipy's own sort orders a coordinate's values by value, and
    a sum of three or more can then differ in its last bit from the one R's
    Matrix gives and the browser takes: a stable sort keeps stored order)."""
    if not matrix.has_sorted_indices and matrix.nnz > 0:
        lengths = np.diff(matrix.indptr)
        lines = np.repeat(np.arange(len(lengths), dtype=np.int64), lengths)
        order = np.lexsort((matrix.indices, lines))
        matrix.indices = matrix.indices[order]
        matrix.data = matrix.data[order]
        matrix.has_sorted_indices = True
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()


def _in_type(matrix: Any, dtype: Any) -> Any:
    """``matrix``'s values in ``dtype``; a value nearer zero than the type holds becomes no value."""
    if matrix.data.dtype == dtype:
        return matrix
    stored = matrix.data.astype(dtype)
    typed = type(matrix)((stored, matrix.indices, matrix.indptr), shape=matrix.shape)
    if not np.all(stored):
        typed = type(matrix)((stored, matrix.indices.copy(), matrix.indptr.copy()), shape=matrix.shape)
        typed.eliminate_zeros()
    return typed


def canonical_block(block: Any, dtype: Any) -> Any:
    """Rows (CSR) or genes (CSC) of an expression matrix, a copy, as ``prepare()``
    stores them (:func:`canonical_csc` then :func:`stored_expression`): duplicate
    coordinates summed, stored zeros dropped, indices sorted, values in ``dtype``."""
    matrix = block.copy()
    _canonical_in_place(matrix)
    return _in_type(matrix, dtype)


def canonical_csc(value: Any) -> sparse.csc_matrix:
    """``value`` as a CSC matrix with sorted, unique indices and no stored zero.

    The caller's matrix is never changed. ``csc_matrix()`` of a CSC input
    shares its ``data``/``indices``/``indptr``, and the canonicalising steps
    work in place on those arrays, so a CSC input that needs any of them is
    copied first; one that is already canonical is used as it is (no copy),
    and any other format is converted into new arrays.
    """
    matrix = sparse.csc_matrix(value)
    shared = sparse.issparse(value) and value.format == "csc"
    if shared and (not matrix.has_canonical_format or not np.all(matrix.data)):
        matrix = matrix.copy()
    _canonical_in_place(matrix)
    return matrix


def _is_backed_sparse(value: object) -> bool:
    """An AnnData opened with ``backed="r"``: its ``X`` is anndata's on-disk sparse dataset."""
    return type(value).__module__.startswith("anndata") and getattr(value, "format", None) in ("csr", "csc")


#: Rows read from a backed CSR matrix at a time while it is turned gene-major
#: (:func:`_backed_csc`): a scheduling choice that keeps each block small
#: beside the gene-major result, not a limit on the matrix.
_BACKED_ROWS_PER_BLOCK = 1 << 16


def _backed_csc(dataset: Any) -> sparse.csc_matrix:
    """A backed AnnData's ``X`` as an in-memory CSC matrix, the file read twice by row blocks.

    A CSR matrix on disk is turned gene-major without ever being held whole
    in memory: one pass counts each gene's values, the next places each block's
    values at their genes' next free positions. Memory holds the CSC result
    and one block, where reading the CSR matrix first would hold both whole.
    """
    if dataset.format == "csc":
        return sparse.csc_matrix(dataset.to_memory())
    n_rows, n_columns = (int(extent) for extent in dataset.shape)
    blocks = [(start, min(n_rows, start + _BACKED_ROWS_PER_BLOCK)) for start in range(0, n_rows, _BACKED_ROWS_PER_BLOCK)]
    counts = np.zeros(n_columns, dtype=np.int64)
    dtype: np.dtype[Any] | None = None
    for start, stop in blocks:
        block = sparse.csr_matrix(dataset[start:stop])
        counts += np.bincount(block.indices, minlength=n_columns)
        dtype = block.data.dtype
    indptr = np.zeros(n_columns + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    nnz = int(indptr[-1])
    index_dtype = np.int64 if max(n_rows, nnz) > np.iinfo(np.int32).max else np.int32
    data = np.empty(nnz, dtype=dtype if dtype is not None else np.float32)
    indices = np.empty(nnz, dtype=index_dtype)
    following = indptr[:-1].copy()
    for start, stop in blocks:
        block = sparse.csr_matrix(dataset[start:stop])
        block.sort_indices()
        rows = np.repeat(np.arange(start, stop, dtype=index_dtype), np.diff(block.indptr))
        order = np.argsort(block.indices, kind="stable")
        columns = block.indices[order]
        per_column = np.bincount(block.indices, minlength=n_columns)
        first = np.cumsum(per_column) - per_column
        positions = following[columns] + (np.arange(columns.size) - first[columns])
        data[positions] = block.data[order]
        indices[positions] = rows[order]
        following += per_column
    return sparse.csc_matrix((data, indices, indptr), shape=(n_rows, n_columns))


def resolve_expression(gene_expression: object, n_cells: int) -> sparse.csc_matrix | np.ndarray:
    """X, cells by genes, gene-major: CSC when sparse, dense otherwise.

    The values keep the type they were given; :func:`stored_expression` gives
    them the type they are stored in once ``gene_identifiers`` has chosen the
    genes. A stored zero is no value and is dropped.
    """
    label = "gene_expression"
    if _is_backed_sparse(gene_expression):
        gene_expression = _backed_csc(gene_expression)
    if sparse.issparse(gene_expression):
        matrix = canonical_csc(gene_expression)
        _require_real(matrix.data, label)
        _require_finite_float32(matrix.data, label, positions=False)
    else:
        if gene_expression is None:
            raise ValueError(f"{label} has no values.")
        array = np.asarray(gene_expression)
        _require_real(array, label)
        if array.ndim != 2:
            raise ValueError(f"{label} must be a 2-D cells x genes matrix, got shape {array.shape}.")
        _require_finite_float32(array, label)
        matrix = array
    if matrix.shape[0] != n_cells:
        raise ValueError(
            f"{label} has {matrix.shape[0]:,} {_plural(matrix.shape[0], 'row', 'rows')}, but the "
            f"embeddings have {n_cells:,} {plural(n_cells, 'cell')}. gene_expression is cells x genes; transpose a "
            "genes x cells matrix first."
        )
    return matrix


def stored_expression(expression: sparse.csc_matrix | np.ndarray) -> sparse.csc_matrix | np.ndarray:
    """X in the type it is stored in: float32, or float64 for integer counts above 2^24.

    Decided on the genes that are stored -- after ``gene_identifiers`` has
    chosen them -- so a large count in a gene left out does not widen the
    genes kept. Every count of the stored genes reads back as it was given.
    """
    values = np.asarray(expression.data if sparse.issparse(expression) else expression)
    dtype = _exact_dtype(values)
    if values.dtype == dtype:
        return expression
    if isinstance(expression, sparse.csc_matrix):
        # Only the values change type: the new matrix shares the index arrays
        # (scipy's astype() would copy them too), which may be the caller's
        # own, unless a value nearer zero than the smallest float32 became
        # zero, which is no value and is dropped, on index arrays of its own.
        result: sparse.csc_matrix = _in_type(expression, dtype)
        return result
    return expression.astype(dtype)


def resolve_genes(
    expression: sparse.csc_matrix | np.ndarray,
    var: object,
    var_gene_id_column: str | None,
    gene_identifiers: Sequence[str] | None,
) -> tuple[sparse.csc_matrix | np.ndarray, np.ndarray, dict[str, Any]]:
    """The genes: X narrowed to ``gene_identifiers``, and the var index and columns."""
    n_genes = int(expression.shape[1])
    if var is None:
        frame = pd.DataFrame(index=pd.RangeIndex(n_genes))
    else:
        frame = _require_frame(var, "var", rows=n_genes, rows_of="gene_expression (columns)")
    if var_gene_id_column is not None and (
        not isinstance(var_gene_id_column, str) or var_gene_id_column not in frame.columns
    ):
        raise ValueError(
            f"var_gene_id_column names var[{var_gene_id_column!r}], which var does not have. "
            f"var holds: {[str(column) for column in frame.columns]}."
        )

    names = _row_names(frame, "var")
    if names is None:
        if var_gene_id_column is None:
            raise ValueError(
                "The genes have no names: give var an index of gene identifiers, or pass "
                "var_gene_id_column naming the var column that holds them."
            )
        names = column_gene_names(frame, var_gene_id_column).names
    unique = make_unique_names(names)
    if unique.renamed:
        note = unique_names_note(filled_from_index=0, renamed=unique.renamed, source=None)
        console_print(f"Gene names: {note}")
    frame = frame.set_axis(pd.Index(unique.names, dtype=object), axis=0)

    if gene_identifiers is not None:
        if isinstance(gene_identifiers, str | bytes) or not isinstance(gene_identifiers, Sequence):
            raise TypeError("gene_identifiers must be a sequence of gene names.")
        wanted = list(gene_identifiers)
        if not all(isinstance(name, str) and name for name in wanted):
            raise TypeError("gene_identifiers must hold non-empty strings.")
        if len(set(wanted)) != len(wanted):
            raise ValueError("gene_identifiers must not contain duplicate identifiers.")
        shown = (
            unique.names
            if var_gene_id_column is None
            else column_gene_names(frame, var_gene_id_column).names
        )
        position = {name: index for index, name in enumerate(shown)}
        absent = [name for name in wanted if name not in position]
        if absent:
            more = f" and {len(absent) - 5:,} more" if len(absent) > 5 else ""
            raise ValueError(
                f"gene_identifiers contains identifiers not found in var: {absent[:5]}{more}. "
                "Identifiers are matched against the gene names the viewer shows "
                "(var_gene_id_column, or the var index), after repeated names were made unique."
            )
        positions = np.asarray([position[name] for name in wanted], dtype=np.int64)
        expression = expression[:, positions]
        if sparse.issparse(expression):
            expression = sparse.csc_matrix(expression)
        frame = frame.iloc[positions]

    if var_gene_id_column is not None:
        # The column is stored as it is; the viewer fills its empty values
        # from the index and makes repeated names unique when it reads it, and
        # this says ahead of time what that will change.
        shown_names = column_gene_names(frame, var_gene_id_column)
        note = unique_names_note(
            filled_from_index=shown_names.filled_from_index,
            renamed=shown_names.renamed,
            source=var_gene_id_column,
        )
        if note is not None:
            console_print(f"Gene names from var column {var_gene_id_column!r}: {note}")
    index = np.asarray(frame.index.tolist(), dtype=object)
    return expression, index, _frame_columns(frame, "var")


# =============================================================================
# obsp: the neighbour graph
# =============================================================================


def resolve_connectivities(connectivities: object, n_cells: int) -> sparse.csr_matrix:
    """``obsp['connectivities']``: a symmetric CSR float32 cells x cells matrix.

    Weights are finite and non-negative, the diagonal is zero, and the matrix
    is exactly symmetric; a stored zero is not an edge and is dropped.
    """
    label = "connectivities"
    if sparse.issparse(connectivities):
        graph = sparse.csr_matrix(cast(Any, connectivities))
    else:
        array = np.asarray(connectivities)
        if array.ndim != 2:
            raise ValueError(f"{label} must be a cells x cells matrix, got shape {array.shape}.")
        graph = sparse.csr_matrix(array)
    _require_real(graph.data, label)
    if graph.shape != (n_cells, n_cells):
        raise ValueError(
            f"{label} must be {n_cells:,} x {n_cells:,} (one row and one column per cell), "
            f"got {graph.shape[0]:,} x {graph.shape[1]:,}."
        )
    graph = graph.astype(np.float64)
    graph.sum_duplicates()
    if not np.isfinite(graph.data).all():
        raise ValueError(
            "Connectivity weights must all be finite. "
            + describe_non_finite(graph.data, label, positions=False)
        )
    if (graph.data < 0).any():
        raise ValueError("Connectivity weights must all be non-negative.")
    graph.eliminate_zeros()
    if graph.diagonal().any():
        raise ValueError("Connectivity diagonal values must all be exactly zero.")
    if (graph != graph.T).nnz:
        raise ValueError("Connectivity matrix must be exactly symmetric.")
    _require_finite_float32(graph.data, label, positions=False)
    graph = graph.astype(np.float32)
    graph.sort_indices()
    return graph


# =============================================================================
# uns['cellucid']
# =============================================================================


def _require_url(value: object, label: str) -> str:
    text = _require_display_text(value, label=label, allow_empty=False)
    parsed = urlparse(text)
    if parsed.scheme not in ("http", "https") or not parsed.netloc:
        raise ValueError(f"{label} must be an absolute http(s) address; got {text!r}.")
    return text


def linked_text(
    text: object,
    url: object,
    *,
    text_label: str,
    url_label: str,
    text_key: str,
) -> dict[str, str] | None:
    """``{text_key: text, "url": url}``, or None when neither is given; a url needs its text."""
    if text is None:
        if url is not None:
            raise ValueError(f"{text_label} is required whenever {url_label} is supplied.")
        return None
    linked = {text_key: _require_display_text(text, label=text_label, allow_empty=False)}
    if url is not None:
        linked["url"] = _require_url(url, url_label)
    return linked
