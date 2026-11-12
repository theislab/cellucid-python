"""
Each cell's total count: ``obs['total_counts']`` and ``uns['cellucid']['totalCounts']``.

A Cellucid-ready AnnData whose ``X`` holds counts records each cell's total
count -- its library size, the ``t`` of library-size normalization
(``x * 10,000 / t``, scanpy's ``normalize_total(target_sum=1e4)``) -- so the
viewer normalizes a gene from one column instead of reading every gene once
to sum them. The record says what the totals cover:

- ``"all genes"``: the caller's own column (``prepare(total_counts=...)``),
  each cell's total over every gene measured, before any gene was selected
  (scanpy's ``calculate_qc_metrics`` run before subsetting, Seurat's
  ``nCount_RNA``). At least the sum of the stored genes, cell by cell.
- ``"source panel"``: the caller's own column, each cell's total over every
  gene of the source's gene panel (``prepare(total_counts=...,
  total_counts_over="source panel")``): a targeted or integration panel
  (the built-in atlases' 8,192 genes), more genes than are stored but not the
  whole transcriptome. Checked as ``"all genes"`` is.
- ``"stored genes"``: the sum of each row of ``X`` as written, which
  :func:`stored_row_totals` computes when no column is named.

``X`` holds counts when every value it stores is a whole number of at least
zero (an integer type, or floats that are all whole). Summed in float64, so
the sum of integer counts below 2^53 is exact and both writers store the same
bytes (an R double and a NumPy float64 are both ``<f8``).

cellucid-r computes and checks the same in ``R/total-counts.R``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import numpy as np
from scipy import sparse

from .._plural import plural
from .._scratch import Scratch
from .._source_arrays import is_lazy_sparse, sparse_source
from ._gene_major import PASS_ENTRIES, GeneMajor
from ._priority import ORIGINAL_POSITION_COLUMN

__all__ = [
    "OVER_ALL_GENES",
    "OVER_SOURCE_PANEL",
    "OVER_STORED_GENES",
    "TOTAL_COUNTS_COLUMN",
    "require_counts",
    "require_totals_column",
    "stored_row_totals",
    "total_counts_column",
    "total_counts_columns",
    "total_counts_over",
    "total_counts_over_by_modality",
    "with_total_counts",
]

#: The obs column the writers compute totals into.
TOTAL_COUNTS_COLUMN = "total_counts"
#: What a record's totals cover.
OVER_ALL_GENES = "all genes"
OVER_SOURCE_PANEL = "source panel"
OVER_STORED_GENES = "stored genes"
#: What a caller's named column may cover (``total_counts_over=``).
NAMED_OVER = (OVER_ALL_GENES, OVER_SOURCE_PANEL)

#: How far below the stored genes' sum a named total may be and still be one
#: (a float32 total of a cell rounds above 2^24).
_RELATIVE_TOLERANCE = 1e-6


def _counts_block(values: np.ndarray) -> bool:
    """Whether every value of a block is a whole number of at least zero."""
    if values.size == 0 or values.dtype.kind == "b":
        return True
    if values.dtype.kind in "iu":
        return bool(values.min() >= 0)
    return bool(values.min() >= 0) and bool(np.all(np.floor(values) == values))


def stored_row_totals(matrix: Any, n_rows: int) -> np.ndarray | None:
    """Each row's sum (float64) when the matrix holds counts, else None.

    ``matrix`` is an expression matrix as :func:`._inputs.resolve_expression`
    gives it (a dense array, any scipy sparse matrix, a backed or lazily read
    one) or as :func:`._inputs.stored_expression` lays it out
    (:class:`._gene_major.GeneMajor`, a dense array). Read one block of
    :data:`._gene_major.PASS_ENTRIES` values at a time. Integer and boolean
    values are whole numbers by type (as the stored type is decided, and as
    cellucid-r reads a logical matrix: 0 and 1).
    """
    if matrix is None:
        return None
    totals = np.zeros(n_rows, dtype=np.float64)
    if isinstance(matrix, np.ndarray):
        if matrix.dtype.kind not in "biuf":
            return None
        rows = max(1, PASS_ENTRIES // max(1, matrix.shape[1]))
        for start in range(0, matrix.shape[0], rows):
            block = np.asarray(matrix[start : start + rows])
            if not _counts_block(block):
                return None
            totals[start : start + rows] = block.sum(axis=1, dtype=np.float64)
        return totals
    if isinstance(matrix, GeneMajor):
        data, indices, nnz = matrix.data, matrix.indices, matrix.nnz
        if np.dtype(data.dtype).kind not in "biuf":
            return None
        for start in range(0, nnz, PASS_ENTRIES):
            values = np.asarray(data[start : start + PASS_ENTRIES])
            if not _counts_block(values):
                return None
            rows = np.asarray(indices[start : start + PASS_ENTRIES])
            totals += np.bincount(rows, weights=values.astype(np.float64), minlength=n_rows)
        return totals
    if not (sparse.issparse(matrix) or is_lazy_sparse(matrix) or hasattr(matrix, "group") or hasattr(matrix, "attrs")):
        raise TypeError(f"{type(matrix).__name__} is not an expression matrix.")
    scratch = Scratch()
    try:
        source = sparse_source(matrix, scratch)
        if np.dtype(source.data.dtype).kind not in "biuf":
            return None
        pointers = np.asarray(source.indptr(), dtype=np.int64) if source.form == "csr" else None
        for start in range(0, source.nnz, PASS_ENTRIES):
            end = min(source.nnz, start + PASS_ENTRIES)
            values = np.asarray(source.data.span(start, end))
            if not _counts_block(values):
                return None
            if pointers is None:
                rows = np.asarray(source.indices.span(start, end))
            else:
                rows = np.searchsorted(pointers, np.arange(start, end, dtype=np.int64), side="right") - 1
            totals += np.bincount(rows, weights=values.astype(np.float64), minlength=n_rows)
        return totals
    finally:
        scratch.close()


def total_counts_column(value: object) -> str | None:
    """``total_counts=`` of an AnnData: the obs column holding each cell's total count over all genes."""
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise TypeError(
            "total_counts names the obs column holding each cell's total count over all genes "
            "(computed before any gene was selected), for example total_counts='total_counts'."
        )
    return value


def total_counts_columns(value: object) -> dict[str, str]:
    """``total_counts=`` of a MuData: each modality's column of its cells' totals."""
    if value is None:
        return {}
    if not isinstance(value, Mapping) or not value:
        raise TypeError(
            "With a MuData, total_counts maps a modality to the column of each cell's total "
            "count over all its genes, for example {'rna': 'total_counts'} (a column of the "
            "MuData's obs, or 'rna:total_counts' for one of mdata.mod['rna'].obs)."
        )
    named: dict[str, str] = {}
    for modality, column in value.items():
        if not isinstance(modality, str) or not isinstance(column, str) or not column:
            raise TypeError("total_counts maps modality names to column names (text).")
        named[modality] = column
    return named


def total_counts_over(value: object, *, named: bool) -> str:
    """``total_counts_over=``: what the column ``total_counts=`` names covers
    (``"all genes"`` when None); refused without a named column, since the
    totals the writer sums itself are over the stored genes."""
    if value is None:
        return OVER_ALL_GENES
    if value not in NAMED_OVER:
        raise ValueError(
            f"total_counts_over must be {NAMED_OVER[0]!r} (every gene measured) or {NAMED_OVER[1]!r} "
            f"(every gene of the source's gene panel, more than are stored), got {value!r}."
        )
    if not named:
        raise ValueError(
            "total_counts_over says what the column total_counts names covers, and no total_counts was "
            "given; name the column of each cell's total count with total_counts=."
        )
    return str(value)


def total_counts_over_by_modality(value: object, named: Mapping[str, str]) -> dict[str, str]:
    """``total_counts_over=`` of a MuData: what each named modality's column
    covers -- one value for every one, or a mapping of modality to value
    (a modality left out: ``"all genes"``)."""
    if value is None or isinstance(value, str):
        over = total_counts_over(value, named=bool(named))
        return {modality: over for modality in named}
    if not isinstance(value, Mapping):
        raise TypeError(
            "With a MuData, total_counts_over is 'all genes' or 'source panel' for every column "
            "total_counts names, or a mapping of modality to one of them, for example {'rna': 'source panel'}."
        )
    unnamed = sorted(str(modality) for modality in value if modality not in named)
    if unnamed:
        raise ValueError(
            f"total_counts_over names {unnamed}, for which total_counts names no column: {sorted(named)}."
        )
    return {modality: total_counts_over(value.get(modality), named=True) for modality in named}


def require_totals_column(
    values: Any, name: str, sums: np.ndarray, cell_names: np.ndarray, *, measured: np.ndarray | None = None
) -> None:
    """A named column holds each cell's total, at least the stored genes' sum."""
    label = f"total_counts names obs column {name!r}, which"
    if not isinstance(values, np.ndarray) or values.dtype.kind not in "iuf":
        raise ValueError(f"{label} is not a column of numbers.")
    numbers = values.astype(np.float64)
    if measured is not None:
        numbers = numbers[measured]
    bad = np.flatnonzero(~np.isfinite(numbers) | (numbers < 0))
    if bad.size:
        raise ValueError(
            f"{label} has {bad.size:,} {plural(int(bad.size), 'cell')} with no total of 0 or more "
            f"(first {cell_names[int(bad[0])]!r}: {float(numbers[int(bad[0])])!r})."
        )
    short = np.flatnonzero(numbers < sums * (1 - _RELATIVE_TOLERANCE) - _RELATIVE_TOLERANCE)
    if short.size:
        cell = int(short[0])
        raise ValueError(
            f"{label} is below the sum of the genes stored for {short.size:,} "
            f"{plural(int(short.size), 'cell')} (first {cell_names[cell]!r}: {numbers[cell]:g} against "
            f"{sums[cell]:g}), so it is not each cell's total count over all genes. Name the column "
            "of total counts (scanpy's calculate_qc_metrics writes total_counts), not a log of it."
        )


def require_counts(sums: np.ndarray | None, *, has_matrix: bool) -> np.ndarray:
    """``sums`` (:func:`stored_row_totals`): a named total needs a matrix of counts."""
    if not has_matrix:
        raise ValueError(
            "total_counts names each cell's total count of gene_expression, and no "
            "gene_expression was supplied."
        )
    if sums is None:
        raise ValueError(
            "total_counts names each cell's total count, and the expression matrix does not "
            "hold counts (a value is not a whole number of at least 0): the totals normalize "
            "counts, so write the counts as X (for example adata.layers['counts'])."
        )
    return sums


def with_total_counts(
    obs: dict[str, Any],
    sums: np.ndarray | None,
    named: str | None,
    *,
    cell_names: np.ndarray,
    field_name: Callable[[str], str] | None = None,
    has_matrix: bool,
    over: str = OVER_ALL_GENES,
) -> tuple[dict[str, Any], dict[str, str] | None]:
    """``obs`` with the cells' totals, and the record of them (None: none).

    ``sums`` are :func:`stored_row_totals` of the matrix as written (None: it
    holds no counts, or there is none), each row a cell of ``obs`` in its
    order; ``named`` the caller's column (``total_counts=``), recorded as
    covering ``over`` (:func:`total_counts_over`); ``field_name``
    the name the viewer shows a column of ``obs`` under (a MuData modality's:
    ``<modality>:<column>``), the column itself when None.
    """
    def show(column: str) -> str:
        return column if field_name is None else field_name(column)

    if named is not None:
        counts = require_counts(sums, has_matrix=has_matrix)
        if named not in obs:
            raise ValueError(f"total_counts names obs column {named!r}, which obs does not hold: {sorted(obs)}.")
        require_totals_column(obs[named], named, counts, cell_names)
        return obs, {"field": show(named), "over": over}
    if sums is None:
        return obs, None
    taken = [column for column in obs if column.lower() == TOTAL_COUNTS_COLUMN]
    if taken:
        # The data's own column: recorded only when it is exactly the stored genes' sums.
        own = obs[taken[0]]
        if (
            taken[0] == TOTAL_COUNTS_COLUMN
            and isinstance(own, np.ndarray)
            and own.dtype.kind in "iuf"
            and np.array_equal(own.astype(np.float64), sums)
        ):
            return obs, {"field": show(TOTAL_COUNTS_COLUMN), "over": OVER_STORED_GENES}
        return obs, None
    # Before the original positions, which a file in priority order keeps last.
    written = {key: values for key, values in obs.items() if key != ORIGINAL_POSITION_COLUMN}
    written[TOTAL_COUNTS_COLUMN] = sums
    if ORIGINAL_POSITION_COLUMN in obs:
        written[ORIGINAL_POSITION_COLUMN] = obs[ORIGINAL_POSITION_COLUMN]
    return written, {"field": show(TOTAL_COUNTS_COLUMN), "over": OVER_STORED_GENES}
