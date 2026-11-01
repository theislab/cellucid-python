"""
Arrays of a source, read by rows or ranges in any order.

What the served store (:mod:`cellucid._served_store`) and the gene-major
layout (:mod:`cellucid.prepare_data._gene_major`) read a matrix or an array
of one row per cell through, whatever holds it: numpy in memory, an h5py
dataset, a Zarr array, anndata's backed sparse dataset. Rows in priority order
are spread over the whole source, so a read asks for thousands of rows (or
short ranges: a cell-major matrix's rows) from every part of it; they are
taken with one fancy index from a view of the whole array that answers any
row at memory speed: the array itself in memory; an uncompressed HDF5 dataset
stored in one piece mapped where it lies in its file (nothing copied); any
other (chunked, compressed) decoded into a memory-mapped file in a scratch
folder (:mod:`cellucid._scratch`), each source chunk once, when a read first
needs it.
"""

from __future__ import annotations

import itertools
import threading
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse

from ._scratch import Scratch

__all__ = [
    "Rows",
    "SparseSource",
    "Vector",
    "encoding_of",
    "is_group",
    "is_lazy_sparse",
    "is_sparse_group",
    "range_index",
    "sparse_source",
]

_ENCODING = "encoding-type"


def encoding_of(node: Any) -> str | None:
    """A stored element's ``encoding-type``, or None."""
    value = node.attrs.get(_ENCODING) if hasattr(node, "attrs") else None
    return value.decode() if isinstance(value, bytes) else value


def is_group(node: Any) -> bool:
    """An h5py or Zarr group (a stored element), never an in-memory table."""
    return hasattr(node, "keys") and hasattr(node, "attrs") and not hasattr(node, "dtype") and not hasattr(node, "iloc")


def is_sparse_group(node: Any) -> bool:
    """A stored CSR or CSC matrix (a group of ``data``, ``indices``, ``indptr``)."""
    return is_group(node) and encoding_of(node) in ("csr_matrix", "csc_matrix")


def range_index(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """The positions of ranges ``[starts[i], ends[i])``, one after another."""
    lengths = ends - starts
    total = int(lengths.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64)
    offsets = np.repeat(starts - np.concatenate([[0], np.cumsum(lengths)[:-1]]), lengths)
    return offsets + np.arange(total, dtype=np.int64)


#: Names of the decoded copies in a scratch folder.
_SOURCE_NAMES = itertools.count()


def _mapped_in_place(array: Any) -> np.ndarray | None:
    """An uncompressed HDF5 dataset stored in one piece, mapped where it lies in its
    file (read-only; nothing copied), or None for any other array."""
    try:
        import h5py
    except ImportError:  # pragma: no cover - h5py is a dependency
        return None
    if not isinstance(array, h5py.Dataset) or array.chunks is not None or array.compression is not None:
        return None
    # A file read through the operating system as it lies on disk, named by
    # its driver's identifier: ``File.driver`` names the default driver
    # "windows" on Windows (HDF5's Windows driver is the sec2 driver there).
    driver = array.file.id.get_access_plist().get_driver()
    if driver not in (h5py.h5fd.SEC2, h5py.h5fd.STDIO) or array.dtype.kind not in "biuf" or array.size == 0:
        return None
    offset = array.id.get_offset()
    if offset is None:
        return None
    return np.memmap(array.file.filename, dtype=array.dtype, mode="r", offset=offset, shape=array.shape)


class Rows:
    """One 1-D or 2-D array of a source, read by rows in any order.

    The view of the whole array is the array itself in memory, an
    uncompressed HDF5 dataset mapped in place, or -- any other -- a
    memory-mapped file in ``scratch`` each source chunk is decoded into when a
    read first needs it (the chunks a read needs found with numpy, under one
    lock), never a step per row or range. The decoded copy costs the disk
    space of the array as decoded, and is removed with the scratch folder. A
    mapped view (in place, or the decoded copy) is let go by :meth:`release`
    (the scratch folder closing): Windows neither removes nor lets another
    program replace a file that is mapped.
    """

    def __init__(self, array: Any, scratch: Scratch | None) -> None:
        self.array = array
        self.shape = tuple(int(extent) for extent in array.shape)
        self.length = self.shape[0]
        self.dtype = np.dtype(array.dtype)
        chunks = getattr(array, "chunks", None)
        self._scratch = scratch
        self._lock = threading.Lock()
        #: Every row at memory speed, once there (None until a chunked source is first read).
        self._view: np.ndarray | None = array if isinstance(array, np.ndarray) else _mapped_in_place(array)
        #: The rows one source chunk holds, and which chunks are decoded into the view.
        self.chunk = int(chunks[0]) if self._view is None and chunks else self.length
        self._ready: np.ndarray | None = None
        #: Every source chunk decoded: rows are taken with no check.
        self._complete = False
        if scratch is not None and isinstance(self._view, np.memmap):
            scratch.keep(self)

    def _whole(self, rows: np.ndarray) -> np.ndarray:
        """The view of the whole array, every source chunk ``rows`` lie in decoded into it."""
        if self._view is not None and (self._ready is None or self._complete):
            return self._view
        needed = np.unique(rows // self.chunk) if len(rows) else np.zeros(0, dtype=np.int64)
        return self._decoded(needed)

    def _decoded(self, needed: np.ndarray) -> np.ndarray:
        with self._lock:
            if self._view is None:
                if self.array is None:
                    raise RuntimeError("The served store is closed.")
                if self._scratch is None:
                    raise RuntimeError("A chunked source array is read without a scratch folder.")
                self._view = np.lib.format.open_memmap(
                    self._scratch.path() / f"source-{next(_SOURCE_NAMES)}.npy", mode="w+", dtype=self.dtype, shape=self.shape
                )
                self._ready = np.zeros(-(-self.length // self.chunk), dtype=bool)
                self._scratch.keep(self)
            assert self._ready is not None
            for index in needed[~self._ready[needed]].tolist():
                start, end = index * self.chunk, min(self.length, (index + 1) * self.chunk)
                self._view[start:end] = self.array[start:end]
                self._ready[index] = True
            self._complete = bool(self._ready.all())
            return self._view

    def release(self) -> None:
        """Let the mapped view go, the file it maps closed (its scratch folder
        is closing). A row asked after is refused."""
        with self._lock:
            self.array = None
            self._view = None
            self._ready = None
            self._complete = False

    def take(self, rows: np.ndarray, columns: slice = slice(None)) -> np.ndarray:
        """Rows ``rows`` (any order, repeats allowed), and ``columns`` of a 2-D array."""
        rows = np.asarray(rows, dtype=np.int64)
        if len(rows) == 0:
            empty = np.zeros((0, *self.shape[1:]), dtype=self.dtype)
            return empty if len(self.shape) == 1 else empty[:, columns]
        whole = self._whole(rows)
        return whole[rows] if len(self.shape) == 1 else whole[rows, columns]

    def span(self, start: int, end: int) -> np.ndarray:
        """Rows ``[start, end)``, one slice of the view (a pass in stored order)."""
        if end <= start:
            return np.zeros((0, *self.shape[1:]), dtype=self.dtype)
        if self._view is None or not (self._ready is None or self._complete):
            first, last = start // self.chunk, (end - 1) // self.chunk
            self._decoded(np.arange(first, last + 1, dtype=np.int64))
        assert self._view is not None
        return self._view[start:end]


class Vector(Rows):
    """One 1-D array of a source (a sparse matrix's indices or data), read by ranges."""

    def ranges(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """The values of ranges ``[starts[i], ends[i])``, one after another."""
        starts = np.asarray(starts, dtype=np.int64)
        ends = np.asarray(ends, dtype=np.int64)
        if len(starts) == 1:
            return np.asarray(self.span(int(starts[0]), int(ends[0])))
        return self.take(range_index(starts, ends))


@dataclass(frozen=True)
class SparseSource:
    """A compressed sparse matrix of a source: its form, shape and three arrays.

    ``matrix`` is the scipy matrix itself when it is held in memory: its lines
    are then gathered by scipy (one pass in compiled code, no index array).
    """

    form: str
    shape: tuple[int, int]
    indptr: Callable[[], np.ndarray]
    indices: Vector
    data: Vector
    matrix: Any = None

    @property
    def nnz(self) -> int:
        return self.indices.length

    @property
    def n_lines(self) -> int:
        """Lines of the major axis: rows of a CSR matrix, columns of a CSC one."""
        return self.shape[0] if self.form == "csr" else self.shape[1]

    @property
    def extent(self) -> int:
        """Positions of the minor axis: columns of a CSR matrix, rows of a CSC one."""
        return self.shape[1] if self.form == "csr" else self.shape[0]

    def lines(self, lines: np.ndarray, pointers: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Major lines ``lines`` (any order): their lengths, and their entries'
        minor positions and values, line after line, each in stored order.
        ``pointers`` is the source's ``indptr`` (read once by the caller)."""
        lines = np.asarray(lines, dtype=np.int64)
        starts, ends = pointers[lines], pointers[lines + 1]
        lengths = (ends - starts).astype(np.int64)
        if len(lines) and int(lines[-1]) - int(lines[0]) == len(lines) - 1 and bool(np.all(np.diff(lines) == 1)):
            # Lines one after another: one slice of each array (a view, in memory).
            first, last = int(starts[0]), int(ends[-1])
            return lengths, self.indices.span(first, last), self.data.span(first, last)
        if self.matrix is not None:
            picked = self.matrix[lines] if self.form == "csr" else self.matrix[:, lines]
            return lengths, np.asarray(picked.indices), np.asarray(picked.data)
        return lengths, self.indices.ranges(starts, ends), self.data.ranges(starts, ends)


def is_lazy_sparse(value: object) -> bool:
    """A dask array of sparse blocks: anndata's lazily read sparse matrix
    (``anndata.experimental.read_lazy``, ``read_elem_lazy``)."""
    return hasattr(value, "dask") and hasattr(value, "compute") and sparse.issparse(getattr(value, "_meta", None))


def _decoded_lazy(value: Any, scratch: Scratch | None) -> SparseSource:
    """A dask array of sparse blocks decoded once, block of rows by block of
    rows (its own row chunks), into a CSR matrix in memory-mapped files in
    ``scratch``: never the whole matrix in memory (computing it would hold
    every value, and scipy's concatenation of its blocks a second copy). Each
    block's entries are kept as stored; the layout makes every line canonical."""
    if scratch is None:
        raise RuntimeError("A lazily read sparse matrix is read without a scratch folder.")
    n_rows, n_columns = (int(extent) for extent in value.shape)
    index_dtype = np.dtype(np.int64 if n_columns > np.iinfo(np.int32).max else np.int32)
    data_dtype = np.dtype(value.dtype)
    indptr = np.zeros(n_rows + 1, dtype=np.int64)
    name = next(_SOURCE_NAMES)
    paths = (scratch.path() / f"lazy-{name}-indices.bin", scratch.path() / f"lazy-{name}-data.bin")
    row = 0
    with open(paths[0], "wb") as indices, open(paths[1], "wb") as data:
        for rows in value.chunks[0]:
            block = value[row : row + int(rows)].compute()
            block = block.tocsr() if sparse.issparse(block) else sparse.csr_matrix(block)
            indptr[row + 1 : row + int(rows) + 1] = indptr[row] + np.asarray(block.indptr[1:], dtype=np.int64)
            np.asarray(block.indices).astype(index_dtype, copy=False).tofile(indices)
            np.asarray(block.data).astype(data_dtype, copy=False).tofile(data)
            row += int(rows)
    nnz = int(indptr[-1])

    def mapped(path: Any, dtype: np.dtype) -> np.ndarray:
        return np.memmap(path, dtype=dtype, mode="r", shape=(nnz,)) if nnz else np.zeros(0, dtype=dtype)

    return SparseSource(
        form="csr",
        shape=(n_rows, n_columns),
        indptr=lambda: indptr,
        indices=Vector(mapped(paths[0], index_dtype), scratch),
        data=Vector(mapped(paths[1], data_dtype), scratch),
    )


def sparse_source(value: Any, scratch: Scratch | None) -> SparseSource:
    """A sparse matrix as a :class:`SparseSource`: scipy CSR or CSC in memory
    (another format converted to CSR), an h5py or Zarr sparse group,
    anndata's backed sparse dataset (read through its stored group), or a dask
    array of sparse blocks (decoded once into ``scratch``, :func:`_decoded_lazy`)."""
    if is_lazy_sparse(value):
        return _decoded_lazy(value, scratch)
    group = getattr(value, "group", None)
    if group is not None and hasattr(value, "format"):
        value = group
    if sparse.issparse(value):
        matrix = value if value.format in ("csr", "csc") else sparse.csr_matrix(value)
        return SparseSource(
            form=matrix.format,
            shape=(int(matrix.shape[0]), int(matrix.shape[1])),
            indptr=lambda: np.asarray(matrix.indptr),
            indices=Vector(np.asarray(matrix.indices), None),
            data=Vector(np.asarray(matrix.data), None),
            matrix=matrix,
        )
    if not is_sparse_group(value):
        raise TypeError(f"{type(value).__name__} is not a sparse matrix.")
    form = "csr" if encoding_of(value) == "csr_matrix" else "csc"
    shape = [int(extent) for extent in value.attrs["shape"]]
    pointers = value["indptr"]
    return SparseSource(
        form=form,
        shape=(shape[0], shape[1]),
        indptr=lambda: np.asarray(pointers[...]),
        indices=Vector(value["indices"], scratch),
        data=Vector(value["data"], scratch),
    )
