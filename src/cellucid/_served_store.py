"""
The store a server presents for data the kernel holds: priority order, made
as it is read.

When Python serves an AnnData object, or an AnnData file (``.h5ad``, Zarr)
that is not in priority order -- or a MuData file (a MuData Zarr store, each
modality's rows in the order of its cells), a Loom file or a 10x Cell Ranger
matrix (as the browser presents them: :mod:`cellucid._served_formats`) --
the viewer reads a Zarr v2 store that holds the same AnnData with its cells
in **priority order** -- the layout
``prepare()`` writes (:mod:`cellucid.prepare_data._priority`; the
Cellucid-ready AnnData specification): row ``r`` is the cell of the ``r``-th
lowest ``uniform`` rank under :data:`SERVED_PRIORITY_SEED`,
``obs['cellucid_original_position']`` names each row's original position --
the one ``prepare()`` gives the cell, by the one rule for a stored order
(:func:`~cellucid.prepare_data._priority.stored_order`): an object in a
verified stored order is put back in its original order first, and any
other source's rows are each their own original position (a position
column with no record that describes the cells names nothing, and is
replaced, with the warning ``prepare()`` gives) -- ``uns['cellucid']['priority']``
records the seed, and a sparse gene-major ``X`` comes with
``varm['cellucid_prefix_counts']``. A subsample the viewer opens with that
seed is then the first ``k`` rows of every per-cell array, and the first
entries of every gene: the server is asked for, and sends, bytes in
proportion to ``k`` whatever the source's own layout.

Nothing of a file is copied that a read can take from it in place (an
uncompressed array stored in one piece); an array stored in chunks is decoded
once into the store's scratch folder as its rows are first read; an object
held in memory is copied once, as the store is made (its snapshot):

- ``obs``, ``var``, ``uns``, ``varm``, ``varp`` (and ``raw``'s ``var`` and
  ``varm``) are written once, uncompressed, by anndata into a store in memory
  (``obs`` in priority order, the stored outlier quantiles with it);
- every array with one row per cell -- ``X``, ``layers``, ``raw/X``,
  ``obsm``, ``obsp`` -- is a plan, and each of its chunks is made when it is
  asked for, from the source (h5py or Zarr arrays on disk: an uncompressed
  HDF5 array stored in one piece mapped where it lies, any other decoded once,
  source chunk by source chunk as served chunks first need them, into a
  memory-mapped file beside the store, so rows in any order are one fancy
  index at memory speed; an object's arrays from their snapshot):

  ======================  ==================================================
  source                  presented as
  ======================  ==================================================
  dense rows              the same array, rows in priority order
  CSR (cell-major)        CSR, rows in priority order (entries move, never
                          re-sorted: a row keeps its own entries)
  CSC (gene-major)        CSC, each gene's row positions renumbered and
                          sorted, so its first entries are its prefix
  ``obsp`` graphs         both axes in priority order
  ======================  ==================================================

- an object held in memory is served from a snapshot taken as the store is
  made, into memory-mapped files beside the store (removed with it, and by a
  later store when the process that made them was killed): its sparse
  expression matrices (``X``, layers, ``raw/X``; CSR, the usual scanpy
  object, or CSC) laid out gene-major in priority order -- the browser reads
  a gene of a cell-major matrix only after indexing every stored value -- with
  their values stored as ``prepare()`` stores them, by the same functions
  (duplicates summed, stored zeros dropped, float32 unless integer counts need
  float64, a non-finite value refused: :class:`_GeneMajor`); its dense arrays
  copied in priority order (an expression matrix in ``prepare()``'s type);
  its sparse ``obsm``/``obsp`` matrices copied with duplicates summed and
  stored zeros dropped. The snapshot costs a pass over each array (about two
  seconds per 100 million stored values) and the disk space of the arrays,
  never a second copy in memory.

Chunks are uncompressed -- the viewer reads them over loopback or a notebook's
comm, where compressing on the kernel would cost more than the bytes saved --
and cut as ``prepare()`` cuts them, by the constants of
:mod:`cellucid.prepare_data._store` (per-cell arrays -- ``obs``, ``obsm``
rows, the stored outlier quantiles, MuData's cell maps -- by
``CELL_CHUNK_ELEMENTS``, an expression matrix's entries by
``SLICED_CHUNK_ELEMENTS``, anything read whole by ``WHOLE_CHUNK_ELEMENTS``;
what anndata writes in memory is re-cut to them): the store a server presents
holds the arrays ``prepare()`` writes, chunk for chunk, as cellucid-r's
served store holds ``cellucid_prepare()``'s.

The viewer shows the object as it was when the store was made: replacing an
element afterwards (``adata.X = ...``) or changing an array in place
(``np.log1p(adata.X.data, out=adata.X.data)``) does not show, and never mixes
into what is served. View the object again to show a change.
"""

from __future__ import annotations

import contextlib
import itertools
import json
import os
import secrets
import shutil
import tempfile
import threading
import warnings
import weakref
from collections import OrderedDict
from collections.abc import Callable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from scipy import sparse

from ._plural import plural
from .prepare_data._priority import (
    ORIGINAL_POSITION_COLUMN,
    PREFIX_COUNTS_KEY,
    CellPriority,
    in_original_order,
    prefix_levels,
    stored_order,
    without_positions,
)
from .prepare_data._store import (
    SLICED_CHUNK_ELEMENTS,
    WHOLE_CHUNK_ELEMENTS,
    _cell_chunks,
    _row_chunks,
    _strip_chunks,
    _vector_chunks,
)

if TYPE_CHECKING:
    import anndata

    from ._anndata_files import AnnDataContents, AnnDataStore
    from ._served_formats import Presented, PresentedModalities

__all__ = ["SERVED_PRIORITY_SEED", "PriorityStore", "file_order"]

#: The seed of a served store's priority order: the subsampling
#: specification's default seed, so the viewer's default subsample of a
#: served dataset is a prefix.
SERVED_PRIORITY_SEED = 0

#: The genes (or rows) whose renumbered, sorted entries one served store
#: keeps (bytes): a long gene spans several chunks of the served store, and
#: is sorted once for all of them while it is kept.
LINE_CACHE_BYTES = 128 * 2**20

#: Stored values one step of a pass over a matrix holds at a time.
_PASS_ENTRIES = 2**23

_ENCODING = "encoding-type"


def _encoding(kind: str, version: str) -> dict[str, str]:
    return {_ENCODING: kind, "encoding-version": version}


# =============================================================================
# Reading a source
# =============================================================================


class _Cache:
    """A least-recently-used cache bounded by the bytes of its values."""

    def __init__(self, budget: int) -> None:
        self._budget = budget
        self._entries: OrderedDict[Any, tuple[Any, int]] = OrderedDict()
        self._bytes = 0
        self._lock = threading.Lock()

    def get(self, key: Any, make: Callable[[], Any], size: Callable[[Any], int]) -> Any:
        with self._lock:
            found = self._entries.get(key)
            if found is not None:
                self._entries.move_to_end(key)
                return found[0]
        value = make()
        nbytes = size(value)
        if nbytes > self._budget:
            return value
        with self._lock:
            if key not in self._entries:
                self._entries[key] = (value, nbytes)
                self._bytes += nbytes
                while self._bytes > self._budget:
                    _, (_, dropped) = self._entries.popitem(last=False)
                    self._bytes -= dropped
        return value


def _range_index(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """The positions of ranges ``[starts[i], ends[i])``, one after another."""
    lengths = ends - starts
    total = int(lengths.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64)
    offsets = np.repeat(starts - np.concatenate([[0], np.cumsum(lengths)[:-1]]), lengths)
    return offsets + np.arange(total, dtype=np.int64)


class _Rows:
    """One 1-D or 2-D array of a source, read by rows in any order.

    Rows in priority order are spread over the whole source, so a served
    chunk asks for thousands of rows (or short ranges, a cell-major matrix's
    rows) from every part of it. They are taken with one fancy index
    (``np.take``) from a view of the whole array that answers any row at
    memory speed: the array itself in memory; an uncompressed HDF5 dataset
    stored in one piece mapped where it is in its file (nothing copied); any
    other (chunked, compressed: h5py, zarr) decoded into a memory-mapped file
    in the store's scratch folder, each source chunk once, when a served
    chunk first needs it (the chunks a served chunk needs found with numpy,
    under one lock), never a step per row or range. The decoded copy costs
    the disk space of the array as decoded, and is removed with the store.
    """

    def __init__(self, array: Any, scratch: _Scratch | None) -> None:
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

    def _whole(self, rows: np.ndarray) -> np.ndarray:
        """The view of the whole array, every source chunk ``rows`` lie in decoded into it."""
        if self._view is not None and (self._ready is None or self._complete):
            return self._view
        needed = np.unique(rows // self.chunk)
        with self._lock:
            if self._view is None:
                if self._scratch is None:
                    raise RuntimeError("A chunked source array is read without the store's scratch folder.")
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
        """Let the decoded copy go (the store is closing; its folder is removed next)."""
        with self._lock:
            if self._ready is not None:
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


#: Names of the decoded copies in a store's scratch folder.
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
    if array.file.driver not in ("sec2", "stdio") or array.dtype.kind not in "biuf" or array.size == 0:
        return None
    offset = array.id.get_offset()
    if offset is None:
        return None
    return np.memmap(array.file.filename, dtype=array.dtype, mode="r", offset=offset, shape=array.shape)


class _Vector(_Rows):
    """One 1-D array of a source (a sparse matrix's indices or data), read by ranges."""

    def ranges(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """The values of ranges ``[starts[i], ends[i])``, one after another."""
        starts = np.asarray(starts, dtype=np.int64)
        ends = np.asarray(ends, dtype=np.int64)
        if len(starts) == 1 and isinstance(self.array, np.ndarray):
            return np.asarray(self.array[int(starts[0]) : int(ends[0])])
        return self.take(_range_index(starts, ends))


@dataclass(frozen=True)
class _Sparse:
    """A compressed sparse matrix of a source: its form, shape and three arrays."""

    form: str
    shape: tuple[int, int]
    indptr: Callable[[], np.ndarray]
    indices: _Vector
    data: _Vector

    @property
    def nnz(self) -> int:
        return self.indices.length


# =============================================================================
# The arrays a served store makes on request
# =============================================================================


class _Planned:
    """A Zarr v2 array whose chunks are made when they are asked for."""

    shape: tuple[int, ...]
    chunks: tuple[int, ...]
    dtype: np.dtype
    attrs: dict[str, Any] | None = None

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        """The values of chunk ``index`` (smaller at the array's edge)."""
        raise NotImplementedError

    def metadata(self) -> dict[str, Any]:
        return {
            "zarr_format": 2,
            "shape": list(self.shape),
            "chunks": list(self.chunks),
            "dtype": self.dtype.str,
            "compressor": None,
            "fill_value": False if self.dtype.kind == "b" else 0,
            "filters": None,
            "order": "C",
            "dimension_separator": ".",
        }

    def chunk_bytes(self, index: tuple[int, ...]) -> bytes | None:
        """Chunk ``index`` as Zarr v2 stores it (padded to the whole chunk), or None outside the grid."""
        grid = [-(-extent // size) for extent, size in zip(self.shape, self.chunks, strict=True)]
        if len(index) != len(grid) or any(not 0 <= i < n for i, n in zip(index, grid, strict=True)):
            return None
        values = np.ascontiguousarray(self.block(index), dtype=self.dtype)
        if values.shape != tuple(self.chunks):
            padded = np.zeros(self.chunks, dtype=self.dtype)
            padded[tuple(slice(0, extent) for extent in values.shape)] = values
            values = padded
        return values.tobytes()


class _Values(_Planned):
    """An array computed once, when it is first read (a matrix's pointers, its prefix counts)."""

    def __init__(self, shape: tuple[int, ...], dtype: Any, chunks: tuple[int, ...], make: Callable[[], np.ndarray]) -> None:
        self.shape = shape
        self.dtype = np.dtype(dtype)
        self.chunks = chunks
        self._make = make
        self._values: np.ndarray | None = None
        self._lock = threading.Lock()

    def values(self) -> np.ndarray:
        with self._lock:
            if self._values is None:
                self._values = np.asarray(self._make(), dtype=self.dtype)
            return self._values

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        values = self.values()
        return values[tuple(slice(i * c, (i + 1) * c) for i, c in zip(index, self.chunks, strict=True))]


class _DenseRows(_Planned):
    """A dense per-cell array with its rows in priority order."""

    def __init__(self, array: Any, order: np.ndarray | None, chunks: tuple[int, ...], scratch: _Scratch, *, encoded: bool = True) -> None:
        self.array = array
        self.source = _Rows(array, scratch)
        self.order = order
        self.shape = tuple(int(extent) for extent in array.shape)
        self.dtype = np.dtype(array.dtype)
        self.chunks = chunks
        self.attrs = _encoding("array", "0.2.0") if encoded else None

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        start = index[0] * self.chunks[0]
        end = min(self.shape[0], start + self.chunks[0])
        columns = (
            slice(index[1] * self.chunks[1], min(self.shape[1], (index[1] + 1) * self.chunks[1]))
            if len(self.shape) == 2
            else slice(None)
        )
        if self.order is None:
            return np.asarray(self.array[start:end] if len(self.shape) == 1 else self.array[start:end, columns])
        return self.source.take(self.order[start:end], columns)


class _DenseGeneMajor(_Planned):
    """A dense expression matrix of a file stored in chunks of whole rows (or
    of more genes than a served strip), laid out once, on its first read,
    gene-major in priority order into a memory-mapped file.

    Read through the rows, every served strip of a gene would decode every
    chunk of the matrix (300,000 x 200 in 1,000 x 200 gzip chunks: 1.2 s a
    gene, each time); laid out once -- one pass over the file by blocks of
    whole rows, memory one block -- a strip is one slice of the file it wrote.
    """

    def __init__(self, array: Any, priority: CellPriority | None, chunks: tuple[int, int], folder: Callable[[], Path], name: str) -> None:
        self.array = array
        self.priority = priority
        self.shape = (int(array.shape[0]), int(array.shape[1]))
        self.dtype = np.dtype(array.dtype)
        self.chunks = chunks
        self.attrs = _encoding("array", "0.2.0")
        self._folder = folder
        self._name = name
        self._lock = threading.Lock()
        self._laid: np.ndarray | None = None

    def _lay_out(self) -> np.ndarray:
        with self._lock:
            if self._laid is None:
                n_obs, n_vars = self.shape
                out = np.lib.format.open_memmap(
                    self._folder() / f"{self._name}-dense.npy", mode="w+", dtype=self.dtype, shape=(n_vars, n_obs)
                )
                step = max(1, _PASS_ENTRIES // max(1, n_vars))
                for start in range(0, n_obs, step):
                    block = np.asarray(self.array[start : start + step])
                    rows = np.arange(start, start + len(block)) if self.priority is None else self.priority.new_row[start : start + len(block)]
                    out[:, rows] = block.T
                out.flush()
                self._laid = out
            return self._laid

    def release(self) -> None:
        with self._lock:
            self._laid = None

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        laid = self._lay_out()
        rows = slice(index[0] * self.chunks[0], (index[0] + 1) * self.chunks[0])
        genes = slice(index[1] * self.chunks[1], (index[1] + 1) * self.chunks[1])
        return np.ascontiguousarray(laid[genes, rows].T)


def _reads_whole_rows(array: Any, strip: tuple[int, int]) -> bool:
    """Whether a stored dense matrix's chunks span more genes than a served strip (rows read whole for any gene)."""
    chunks = getattr(array, "chunks", None)
    if isinstance(array, np.ndarray) or not chunks or len(chunks) != 2:
        return False
    return int(chunks[1]) > strip[1]


class _Entries(_Planned):
    """The ``indices`` or ``data`` of a served sparse matrix."""

    def __init__(self, matrix: _ServedSparse, part: str, chunk: int) -> None:
        self.matrix = matrix
        self.part = part
        self.shape = (matrix.nnz,)
        self.dtype = matrix.indices_dtype if part == "indices" else matrix.source.data.dtype
        self.chunks = _vector_chunks(matrix.nnz, chunk)

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        start = index[0] * self.chunks[0]
        indices, data = self.matrix.entries(start, min(self.shape[0], start + self.chunks[0]))
        return indices if self.part == "indices" else data


class _ServedSparse:
    """A source's compressed sparse matrix as served: its major axis in
    priority order (``major``: the source line of each served line) and/or its
    minor axis renumbered (``minor``: the served position of each source
    position, entries sorted by it within each line), or cut to a window of
    it (``window``: minor positions ``[lo, hi)``, each line's entries there
    one run, as in a matrix whose lines are sorted: a Cell Ranger matrix's
    features of one type)."""

    def __init__(
        self,
        source: _Sparse,
        *,
        major: np.ndarray | None,
        minor: np.ndarray | None,
        cache: _Cache,
        window: tuple[int, int] | None = None,
    ) -> None:
        self.source = source
        self.form = source.form
        self.major = major
        self.minor = minor
        self.window = window
        major_extent = source.shape[0] if self.form == "csr" else source.shape[1]
        #: The source's pointers, read once, when first needed.
        self._source = _Values((major_extent + 1,), np.int64, (major_extent + 1,), source.indptr)
        if window is None:
            self.shape = source.shape
            self.nnz = source.nnz
            self._runs: tuple[np.ndarray, np.ndarray] | None = None
        else:
            lo, hi = window
            self.shape = (source.shape[0], hi - lo) if self.form == "csr" else (hi - lo, source.shape[1])
            # Where each line's entries in the window are: one pass over the indices, now (the store's shape needs them).
            self._runs = _window_runs(source, self._source.values(), lo, hi)
            self.nnz = int((self._runs[1] - self._runs[0]).sum())
        minor_extent = self.shape[1] if self.form == "csr" else self.shape[0]
        self.indices_dtype = np.dtype(np.int64 if minor_extent > np.iinfo(np.int32).max else np.int32)
        self.indptr_dtype = np.dtype(np.int64 if self.nnz >= 2**31 else np.int32)
        self._cache = cache
        self.pointers = _Values(
            (major_extent + 1,), self.indptr_dtype, _vector_chunks(major_extent + 1, WHOLE_CHUNK_ELEMENTS), self._served_indptr
        )

    def _source_indptr(self) -> np.ndarray:
        return self._source.values()

    def _bounds(self, lines: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The source entry runs ``[starts, ends)`` of source lines ``lines``."""
        if self._runs is not None:
            return self._runs[0][lines], self._runs[1][lines]
        pointers = self._source_indptr()
        return pointers[lines], pointers[lines + 1]

    def _served_indptr(self) -> np.ndarray:
        if self.major is None and self._runs is None:
            return self._source_indptr()
        lines = np.arange(len(self._source_indptr()) - 1, dtype=np.int64)
        starts, ends = self._bounds(lines if self.major is None else self.major)
        return np.concatenate([[0], np.cumsum(ends - starts, dtype=np.int64)])

    def _line(self, line: int) -> tuple[np.ndarray, np.ndarray]:
        """One served line's entries (sorted by their renumbered positions), kept while it is long."""

        def make() -> tuple[np.ndarray, np.ndarray]:
            source_line = np.asarray([line if self.major is None else int(self.major[line])], dtype=np.int64)
            starts, ends = self._bounds(source_line)
            return self._renumbered(self.source.indices.ranges(starts, ends), self.source.data.ranges(starts, ends), ends - starts)

        indices, data = self._cache.get(("line", id(self), line), make, lambda v: v[0].nbytes + v[1].nbytes)
        return np.asarray(indices), np.asarray(data)

    def _renumbered(self, indices: np.ndarray, data: np.ndarray, lengths: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.window is not None:
            return np.asarray(indices, dtype=np.int64) - self.window[0], data
        if self.minor is None:
            return indices, data
        positions = self.minor[indices]
        # Within each line by position: one sort of line * extent + position
        # (stable, as lines are kept in order and a repeated position keeps
        # its stored order), a third of what sorting by two keys costs.
        extent = len(self.minor)
        lines = np.repeat(np.arange(len(lengths), dtype=np.int64) * extent, lengths)
        order = np.argsort(lines + positions, kind="stable")
        return positions[order], data[order]

    def entries(self, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        """Served entries ``[start, end)``: their minor positions and values. A
        served chunk's ``indices`` and ``data`` are asked for apart, and made
        together once (kept for the other)."""
        indices, data = self._cache.get(
            ("entries", id(self), start, end), lambda: self._entries(start, end), lambda v: v[0].nbytes + v[1].nbytes
        )
        return np.asarray(indices), np.asarray(data)

    def _entries(self, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        pointers = self.pointers.values()
        first = int(np.searchsorted(pointers, start, side="right")) - 1
        last = int(np.searchsorted(pointers, end, side="left"))
        lines = np.arange(first, last, dtype=np.int64)
        source_lines = lines if self.major is None else self.major[first:last]
        starts, ends = self._bounds(source_lines)
        offset = start - int(pointers[first])
        if self.minor is not None and int((ends - starts).max(initial=0)) >= SLICED_CHUNK_ELEMENTS:
            # A long line (a gene of a gene-major matrix) is sorted once and
            # kept for the chunks it spans.
            parts = [self._line(int(line)) for line in lines]
            indices = np.concatenate([part[0] for part in parts])
            data = np.concatenate([part[1] for part in parts])
        else:
            indices, data = self._renumbered(
                self.source.indices.ranges(starts, ends), self.source.data.ranges(starts, ends), ends - starts
            )
        return indices[offset : offset + end - start], data[offset : offset + end - start]

    def prefix_counts(self, levels: np.ndarray) -> np.ndarray:
        """``varm['cellucid_prefix_counts']`` of a served gene-major matrix: one pass over its genes."""
        assert self.form == "csc"
        pointers = self._source_indptr()
        n_vars = self.shape[1]
        counts = np.zeros((n_vars, len(levels)), dtype=np.uint32)
        gene = 0
        while gene < n_vars:
            # Genes holding about _PASS_ENTRIES values at a time.
            end = int(np.searchsorted(pointers, pointers[gene] + _PASS_ENTRIES, side="right")) - 1
            end = min(n_vars, max(gene + 1, end))
            starts, ends = pointers[gene:end], pointers[gene + 1 : end + 1]
            rows = self.source.indices.ranges(starts, ends)
            if self.minor is not None:
                rows = self.minor[rows]
            counts[gene:end] = _counts_below(rows, ends - starts, levels)
            gene = end
        return counts


def _window_runs(source: _Sparse, pointers: np.ndarray, lo: int, hi: int) -> tuple[np.ndarray, np.ndarray]:
    """Each line's run of entries at minor positions ``[lo, hi)``: ``(starts, ends)``, one pass over the indices.

    Lines must hold their entries sorted (as Cell Ranger writes them), so a
    window's entries are one run of each line.
    """
    n_lines = len(pointers) - 1
    starts = np.empty(n_lines, dtype=np.int64)
    ends = np.empty(n_lines, dtype=np.int64)
    line = 0
    while line < n_lines:
        stop = int(np.searchsorted(pointers, pointers[line] + _PASS_ENTRIES, side="right")) - 1
        stop = min(n_lines, max(line + 1, stop))
        first, last = int(pointers[line]), int(pointers[stop])
        indices = np.asarray(source.indices.ranges(np.asarray([first]), np.asarray([last])), dtype=np.int64)
        local = pointers[line : stop + 1] - first
        if len(indices) > 1:
            rising = np.diff(indices) > 0
            # A step down is allowed only where a line ends.
            breaks = np.zeros(len(indices) - 1, dtype=bool)
            inner = local[1:-1]
            breaks[inner[(inner > 0) & (inner < len(indices))] - 1] = True
            if not np.all(rising | breaks):
                raise ValueError(
                    "This matrix holds its entries unsorted within a cell, so one feature type's columns cannot be "
                    "read as one run: read it with scanpy (sc.read_10x_h5(path, gex_only=False)), write it as .h5ad "
                    "and serve that."
                )
        lengths = np.diff(local)
        owner = np.repeat(np.arange(stop - line, dtype=np.int64), lengths)
        below_lo = np.bincount(owner, weights=indices < lo, minlength=stop - line).astype(np.int64)
        below_hi = np.bincount(owner, weights=indices < hi, minlength=stop - line).astype(np.int64)
        starts[line:stop] = pointers[line:stop] + below_lo
        ends[line:stop] = pointers[line:stop] + below_hi
        line = stop
    return starts, ends


def _counts_below(rows: np.ndarray, lengths: np.ndarray, levels: np.ndarray) -> np.ndarray:
    """For lines of ``lengths`` entries at ``rows`` (any order within a line):
    each line's entries below each level."""
    n_lines, n_levels = len(lengths), len(levels)
    if n_levels == 0:
        return np.zeros((n_lines, 0), dtype=np.uint32)
    # An entry at row r counts toward level j when r < levels[j]: its first level is
    # the number of levels <= r.
    first = np.searchsorted(levels, rows, side="right")
    lines = np.repeat(np.arange(n_lines, dtype=np.int64), lengths)
    histogram = np.bincount(lines * (n_levels + 1) + first, minlength=n_lines * (n_levels + 1))
    return np.cumsum(histogram.reshape(n_lines, n_levels + 1), axis=1)[:, :n_levels].astype(np.uint32)


class _GeneMajor:
    """An expression matrix held in memory (sparse, any format), laid out once,
    gene-major in priority order, into memory-mapped files in ``folder``.

    The layout is the store's snapshot of the matrix: made when the store is
    made, before the caller has the object back, so the bytes served are the
    matrix as it was then whatever is done to the object later (an in-place
    ``np.log1p(adata.X.data, out=adata.X.data)`` is not seen; view the object
    again to show it). Its values are stored as ``prepare()`` stores them, by
    the same functions (:mod:`cellucid.prepare_data._inputs`): duplicate
    coordinates summed, stored zeros dropped, float32 unless integer counts
    need float64, a non-finite value refused. Each gene's entries land in
    ascending served row; the prefix counts are its entries below each level.

    A cell-major (CSR) matrix is read twice by blocks of rows in priority
    order (counting, then placing); a gene-major (CSC) one twice by blocks of
    genes. Memory holds one block beside the files; the files hold the matrix
    once (its stored values, at 4 or 8 bytes, and a 4- or 8-byte row each).
    """

    def __init__(
        self,
        matrix: Any,
        priority: CellPriority,
        folder: Callable[[], Path],
        name: str,
        stop: threading.Event,
        levels: np.ndarray | None,
        label: str,
    ) -> None:
        from .prepare_data._inputs import expression_dtype, require_expression_values

        self._stop = stop
        self._lock = threading.Lock()
        matrix = matrix if matrix.format in ("csr", "csc") else sparse.csr_matrix(matrix)
        n_obs, n_vars = matrix.shape
        self.shape = (int(n_obs), int(n_vars))
        values = np.asarray(matrix.data)
        require_expression_values(values, label)
        self.dtype = np.dtype(expression_dtype(values))
        self.indices_dtype = np.dtype(np.int64 if n_obs > np.iinfo(np.int32).max else np.int32)
        self._levels = np.zeros(0, dtype=np.int64) if levels is None else np.asarray(levels, dtype=np.int64)
        self._laid: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
        lay_out = self._lay_out_rows if matrix.format == "csr" else self._lay_out_genes
        self._laid = lay_out(matrix, priority, folder(), name)
        self.nnz = int(self._laid[0][-1])

    def laid_out(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        with self._lock:
            if self._laid is None:
                raise RuntimeError("The served store is closed.")
            return self._laid

    def release(self) -> None:
        """Let go of the memory-mapped arrays, so their files can be removed (Windows keeps a mapped file)."""
        with self._lock:
            self._laid = None

    def counts(self) -> np.ndarray:
        return self._counts

    def _files(self, folder: Path, name: str, nnz: int) -> tuple[np.ndarray, np.ndarray]:
        def mapped(part: str, dtype: Any) -> np.ndarray:
            if nnz == 0:
                return np.zeros(0, dtype=dtype)
            return np.lib.format.open_memmap(folder / f"{name}-{part}.npy", mode="w+", dtype=dtype, shape=(nnz,))

        return mapped("indices", self.indices_dtype), mapped("data", self.dtype)

    def _checked(self) -> None:
        if self._stop.is_set():
            raise RuntimeError("The served store was closed while its matrix was laid out.")

    def _lay_out_rows(
        self, matrix: sparse.csr_matrix, priority: CellPriority, folder: Path, name: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        from .prepare_data._inputs import canonical_block

        n_obs, n_vars = self.shape
        order = priority.order
        lengths = np.diff(np.asarray(matrix.indptr, dtype=np.int64))[order]
        served_end = np.cumsum(lengths, dtype=np.int64)
        levels = self._levels
        # Steps of rows ending at every level, so the cursors there are the prefix counts.
        boundaries = np.unique(np.concatenate([levels[levels < n_obs], [n_obs]]))

        def steps() -> Iterator[tuple[int, int]]:
            row = 0
            for boundary in boundaries.tolist():
                while row < boundary:
                    self._checked()
                    base = 0 if row == 0 else int(served_end[row - 1])
                    stop = int(np.searchsorted(served_end, base + _PASS_ENTRIES, side="right"))
                    stop = min(boundary, max(row + 1, stop))
                    yield row, stop
                    row = stop

        def block(row: int, stop: int) -> sparse.csc_matrix:
            # The step's rows in priority order, by gene (each gene's rows ascending).
            result: sparse.csc_matrix = canonical_block(sparse.csc_matrix(matrix[order[row:stop]]), self.dtype)
            return result

        per_gene = np.zeros(n_vars, dtype=np.int64)
        for row, stop in steps():
            per_gene += np.diff(block(row, stop).indptr)
        indptr = np.concatenate([[0], np.cumsum(per_gene, dtype=np.int64)])
        out_rows, out_data = self._files(folder, name, int(indptr[-1]))
        cursor = indptr[:-1].copy()
        counts = np.zeros((n_vars, len(levels)), dtype=np.uint32)
        level_at = {int(level): j for j, level in enumerate(levels.tolist())}
        for row, stop in steps():
            part = block(row, stop)
            step = np.diff(part.indptr)
            within = np.arange(part.nnz, dtype=np.int64) - np.repeat(part.indptr[:-1].astype(np.int64), step)
            places = np.repeat(cursor, step) + within
            out_rows[places] = part.indices.astype(np.int64) + row
            out_data[places] = part.data
            cursor += step
            if stop in level_at:
                counts[:, level_at[stop]] = cursor - indptr[:-1]
        for j, level in enumerate(levels.tolist()):
            if level >= n_obs:
                counts[:, j] = per_gene
        return self._finished(indptr, out_rows, out_data, counts)

    def _lay_out_genes(
        self, matrix: sparse.csc_matrix, priority: CellPriority, folder: Path, name: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        from .prepare_data._inputs import canonical_block

        n_obs, n_vars = self.shape
        pointers = np.asarray(matrix.indptr, dtype=np.int64)

        def steps() -> Iterator[tuple[int, int]]:
            gene = 0
            while gene < n_vars:
                self._checked()
                end = int(np.searchsorted(pointers, pointers[gene] + _PASS_ENTRIES, side="right")) - 1
                end = min(n_vars, max(gene + 1, end))
                yield gene, end
                gene = end

        def block(gene: int, end: int) -> sparse.csc_matrix:
            part: sparse.csc_matrix = canonical_block(matrix[:, gene:end], self.dtype)
            return part

        per_gene = np.zeros(n_vars, dtype=np.int64)
        for gene, end in steps():
            per_gene[gene:end] = np.diff(block(gene, end).indptr)
        indptr = np.concatenate([[0], np.cumsum(per_gene, dtype=np.int64)])
        out_rows, out_data = self._files(folder, name, int(indptr[-1]))
        counts = np.zeros((n_vars, len(self._levels)), dtype=np.uint32)
        for gene, end in steps():
            part = block(gene, end)
            lengths = np.diff(part.indptr)
            rows = priority.new_row[part.indices]
            owner = np.repeat(np.arange(end - gene, dtype=np.int64), lengths)
            # Each gene's entries by served row.
            order = np.lexsort((rows, owner))
            first, last = int(indptr[gene]), int(indptr[end])
            out_rows[first:last] = rows[order]
            out_data[first:last] = part.data[order]
            counts[gene:end] = _counts_below(rows, lengths, self._levels)
        return self._finished(indptr, out_rows, out_data, counts)

    def _finished(
        self, indptr: np.ndarray, out_rows: np.ndarray, out_data: np.ndarray, counts: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        for array in (out_rows, out_data):
            if isinstance(array, np.memmap):
                array.flush()
        self._counts = counts
        return indptr, out_rows, out_data


class _GeneMajorEntries(_Planned):
    def __init__(self, matrix: _GeneMajor, part: str) -> None:
        self.matrix = matrix
        self.part = part
        self.shape = (matrix.nnz,)
        self.dtype = np.dtype(matrix.indices_dtype if part == "indices" else matrix.dtype)
        self.chunks = _vector_chunks(matrix.nnz, SLICED_CHUNK_ELEMENTS)

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        _, rows, data = self.matrix.laid_out()
        values = rows if self.part == "indices" else data
        return values[index[0] * self.chunks[0] : (index[0] + 1) * self.chunks[0]]


# =============================================================================
# Elements of a source
# =============================================================================


def _encoding_of(node: Any) -> str | None:
    value = node.attrs.get(_ENCODING) if hasattr(node, "attrs") else None
    return value.decode() if isinstance(value, bytes) else value


def _is_group(node: Any) -> bool:
    """An h5py or Zarr group (a stored element), never an in-memory table."""
    return hasattr(node, "keys") and hasattr(node, "attrs") and not hasattr(node, "dtype") and not hasattr(node, "iloc")


def _sparse_of_memory(matrix: Any) -> _Sparse:
    form = "csr" if matrix.format == "csr" else "csc"
    return _Sparse(
        form=form,
        shape=(int(matrix.shape[0]), int(matrix.shape[1])),
        indptr=lambda: np.asarray(matrix.indptr),
        indices=_Vector(np.asarray(matrix.indices), None),
        data=_Vector(np.asarray(matrix.data), None),
    )


def _sparse_of_group(group: Any, scratch: _Scratch) -> _Sparse:
    form = "csr" if _encoding_of(group) == "csr_matrix" else "csc"
    shape = [int(extent) for extent in group.attrs["shape"]]
    pointers = group["indptr"]
    return _Sparse(
        form=form,
        shape=(shape[0], shape[1]),
        indptr=lambda: np.asarray(pointers[...]),
        indices=_Vector(group["indices"], scratch),
        data=_Vector(group["data"], scratch),
    )


def _is_sparse_group(node: Any) -> bool:
    return _is_group(node) and _encoding_of(node) in ("csr_matrix", "csc_matrix")


def _as_element(value: Any) -> Any:
    """An element as the store reads it: numpy, scipy CSR/CSC, or an h5py/Zarr node.

    anndata's backed sparse datasets are read through their stored group; a
    dask array is computed; a ``numpy.matrix`` is a plain array.
    """
    group = getattr(value, "group", None)
    if group is not None and hasattr(value, "format"):
        return group
    if sparse.issparse(value):
        return value if value.format in ("csr", "csc") else sparse.csr_matrix(value)
    if isinstance(value, np.matrix):
        return np.asarray(value)
    if hasattr(value, "compute"):
        return np.asarray(value.compute())
    return value


def _is_dense_array(value: Any) -> bool:
    """A numeric numpy, h5py or Zarr array (read by rows)."""
    return hasattr(value, "dtype") and hasattr(value, "shape") and not _is_group(value) and np.dtype(value.dtype).kind in "biuf"


def _take_rows(value: Any, order: np.ndarray) -> Any:
    """A per-cell element with its rows in priority order (a table, an array, anything indexable)."""
    if hasattr(value, "iloc"):
        return value.iloc[order]
    if hasattr(value, "dtype") and not isinstance(value, np.ndarray):
        value = np.asarray(value[...])
    return value[order]


def _permute_quantiles(uns: Any, order: np.ndarray, n_obs: int) -> Any:
    """``uns`` with the stored outlier quantiles (one value per cell) in priority order."""
    if not isinstance(uns, Mapping):
        return uns
    cellucid = uns.get("cellucid")
    quantiles = cellucid.get("outlierQuantiles") if isinstance(cellucid, Mapping) else None
    fields = quantiles.get("fields") if isinstance(quantiles, Mapping) else None
    if not isinstance(fields, Mapping):
        return uns
    assert isinstance(cellucid, Mapping) and isinstance(quantiles, Mapping)
    placed = {
        key: (np.asarray(values)[order] if np.ndim(values) >= 1 and len(values) == n_obs else values)
        for key, values in fields.items()
    }
    return {**uns, "cellucid": {**cellucid, "outlierQuantiles": {**quantiles, "fields": placed}}}


def _without(mapping: Mapping[str, Any], key: str) -> dict[str, Any]:
    return {name: value for name, value in mapping.items() if name != key}


#: The prefix of a store's scratch folder in the temporary folder; a file
#: ``owner`` in it names the process that made it.
_SCRATCH_PREFIX = "cellucid-served-"


def _process_alive(pid: int) -> bool:
    """Whether process ``pid`` runs (on this machine), never signalling it."""
    if pid <= 0:
        return False
    if os.name == "nt":  # pragma: no cover - Windows only
        import ctypes

        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        try:
            code = ctypes.c_ulong()
            return bool(kernel32.GetExitCodeProcess(handle, ctypes.byref(code))) and code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def remove_stale_scratch(temporary: Path | None = None) -> list[Path]:
    """Remove the scratch folders of stores whose process is gone (a kernel
    killed before it closed them), and return them. A folder whose owner runs,
    or names none it can read, is left alone."""
    removed = []
    for folder in Path(temporary or tempfile.gettempdir()).glob(f"{_SCRATCH_PREFIX}*"):
        try:
            pid = int((folder / "owner").read_text().strip())
        except (OSError, ValueError):
            continue
        if pid != os.getpid() and not _process_alive(pid):
            shutil.rmtree(folder, ignore_errors=True)
            removed.append(folder)
    return removed


class _Scratch:
    """The folder a store's layouts and snapshots are written to, made when
    first needed and removed with the store (after its memory-mapped arrays
    are let go). It names its process (``owner``), so a later store removes it
    when that process was killed before closing it (:func:`remove_stale_scratch`)."""

    def __init__(
        self,
        closers: list[Callable[[], None]],
        stop: threading.Event,
        layouts: list[_GeneMajor | _DenseGeneMajor | _Rows],
        snapshots: list[np.ndarray],
    ) -> None:
        self._closers = closers
        self._stop = stop
        self._layouts = layouts
        self._snapshots = snapshots
        self._lock = threading.Lock()
        self.folder: Path | None = None

    def path(self) -> Path:
        with self._lock:
            if self._stop.is_set():
                raise RuntimeError("The served store is closed.")
            if self.folder is None:
                remove_stale_scratch()
                folder = Path(tempfile.mkdtemp(prefix=_SCRATCH_PREFIX))
                (folder / "owner").write_text(str(os.getpid()))
                self.folder = folder
                layouts = self._layouts
                snapshots = self._snapshots

                def remove() -> None:
                    for layout in layouts:
                        layout.release()
                    snapshots.clear()
                    shutil.rmtree(folder, ignore_errors=True)

                self._closers.append(remove)
            return self.folder

    def keep(self, layout: _Rows) -> None:
        """Let ``layout``'s memory-mapped file go before the folder is removed."""
        with self._lock:
            self._layouts.append(layout)


# =============================================================================
# The store
# =============================================================================


class PriorityStore:
    """A served AnnData store in priority order (module docstring).

    ``get(key)`` answers one store key: a metadata document, an element
    written in memory, or a chunk made from the source; None for a key the
    store does not have. ``etag`` names this store's version (one per
    store, so a browser never mixes two). Close it to release an open file
    and the files the layout of a cell-major matrix wrote.
    """

    def __init__(self) -> None:
        import zarr

        self._documents: dict[str, Any] = {}
        self._memory = zarr.storage.MemoryStore(store_dict=self._documents)
        self.root = zarr.open_group(self._memory, mode="w", zarr_format=2)
        self.root.attrs.update(_encoding("anndata", "0.1.0"))
        self._planned: dict[str, _Planned] = {}
        self._lines = _Cache(LINE_CACHE_BYTES)
        self._closers: list[Callable[[], None]] = []
        #: Set when the store closes: a layout still running stops at its next step.
        self._stop = threading.Event()
        #: Matrices laid out into the scratch folder (let go before it is removed).
        self._layouts: list[_GeneMajor | _DenseGeneMajor | _Rows] = []
        #: Dense arrays copied from an object held in memory (memory-mapped), let go on close.
        self._snapshots: list[np.ndarray] = []
        #: The folder the layouts write (nothing here refers back to the store,
        #: so dropping the store releases it).
        self._folder = _Scratch(self._closers, self._stop, self._layouts, self._snapshots)
        self._finalizer = weakref.finalize(self, PriorityStore._release, self._closers, self._stop)
        self.etag = f'"served-{secrets.token_hex(12)}"'
        self._consolidated: bytes | None = None
        #: The seed of the store's priority order.
        self.seed = SERVED_PRIORITY_SEED
        self.n_obs = 0
        #: What the store holds: ``"anndata"``, or ``"mudata"`` (modalities under ``mod``).
        self.layout = "anndata"
        #: The rows of each AnnData group (``""``, ``"mod/<name>/"``): its arrays of one row per cell are chunked as such.
        self._groups: dict[str, int] = {}

    # -- building ------------------------------------------------------------

    @classmethod
    def of_anndata(cls, adata: anndata.AnnData) -> PriorityStore:
        """The store of an AnnData object, its cells in priority order.

        Its cells' original positions are the ones ``prepare()`` gives them,
        by the one rule (:func:`~cellucid.prepare_data._priority.stored_order`):
        an object in a verified stored order (read from a file in priority
        order, its record kept) is put back in its original order first; any
        other is taken in its own order, and a position column with no record
        that describes the cells names nothing (it is not served; a warning
        says so, as ``prepare()`` says it).
        """
        if adata.n_obs <= 0:
            raise ValueError("AnnData datasets must contain at least one observation")
        adata = in_original_order(adata)
        store = cls()
        try:
            store._plan_object(adata)
            store._settle()
        except BaseException:
            store.close()
            raise
        return store

    @classmethod
    def of_file(cls, source: AnnDataStore) -> PriorityStore:
        """The store of an AnnData file (``.h5ad``, a Zarr folder or zip).

        A file in priority order already (``uns['cellucid']['priority']``,
        served this way because the browser cannot decode one of its codecs)
        keeps its own order; any other is put in priority order.
        """
        store = cls()
        try:
            store._plan_file(store._open(source))
            store._settle()
        except BaseException:
            store.close()
            raise
        return store

    @classmethod
    def of_mudata_file(cls, source: AnnDataStore) -> PriorityStore:
        """The store of a MuData file (``.h5mu``, a MuData Zarr folder or zip), its cells in priority order.

        The shared cells are in priority order, and each modality's rows in
        the order of its cells (``prepare()``'s MuData layout): a prefix of
        the cells is a prefix of every modality.
        """
        store = cls()
        try:
            store._plan_mudata_file(store._open(source))
            store._settle()
        except BaseException:
            store.close()
            raise
        return store

    @classmethod
    def of_loom(cls, source: AnnDataStore) -> PriorityStore:
        """The store of a Loom file, presented as the browser presents it (:mod:`cellucid._served_formats`)."""
        from ._served_formats import loom_presented

        store = cls()
        try:
            store._plan_presented(loom_presented(store._open(source), source.path.name))
            store._settle()
        except BaseException:
            store.close()
            raise
        return store

    @classmethod
    def of_cell_ranger(cls, source: AnnDataStore) -> PriorityStore:
        """The store of a 10x Cell Ranger ``.h5``, or of an ``outs/`` folder with its analysis
        (:mod:`cellucid._served_formats`): one AnnData, or a MuData of its feature types."""
        from ._served_formats import cell_ranger_presented, outs_analysis, outs_matrix_file

        store = cls()
        try:
            if source.layout == "folder":
                found = outs_matrix_file(source.files or ())
                if found is None:
                    raise ValueError(f"{source.path.name} holds no Cell Ranger count matrix file.")
                path, base = found
                root = store._open_hdf5(source.path / path)
                files = source.files or ()
                presented = cell_ranger_presented(
                    root, source.path.name, analysis=lambda barcodes: outs_analysis(source.path, files, base, barcodes)
                )
            else:
                presented = cell_ranger_presented(store._open(source), source.path.name, analysis=None)
            store._plan_presented(presented)
            store._settle()
        except BaseException:
            store.close()
            raise
        return store

    def _open_hdf5(self, path: Path) -> Any:
        import h5py

        # A larger chunk cache per dataset than h5py's 1 MiB: rows in
        # priority order revisit chunks in no order.
        file = h5py.File(path, "r", rdcc_nbytes=64 * 2**20, rdcc_nslots=100003)
        self._closers.append(file.close)
        return file

    def _open(self, source: AnnDataStore) -> Any:
        from ._anndata_files import _container_of

        if _container_of(source) == "hdf5":
            return self._open_hdf5(source.path)
        import zarr

        if source.layout == "zip":
            zip_store = zarr.storage.ZipStore(source.path, mode="r")
            self._closers.append(zip_store.close)
            return zarr.open_group(zip_store, path=source.zip_root, mode="r")
        return zarr.open_group(str(source.path), mode="r")

    def _write(self, path: str, value: Any) -> None:
        """One element written in memory by anndata (uncompressed: :meth:`_settle`).

        Every chunk is stored, even one that holds only the fill value: zarr
        leaves such a chunk out by default, and a scalar's (a seed of 0)
        has no fill value a reader could take instead.
        """
        import zarr
        from anndata.io import write_elem

        parent, _, key = path.rpartition("/")
        with zarr.config.set({"array.write_empty_chunks": True}):
            write_elem(self.root[parent] if parent else self.root, key, value, dataset_kwargs={"compressor": None})

    def _group(self, path: str, attrs: Mapping[str, Any] | None = None) -> None:
        self.root.create_group(path).attrs.update(dict(attrs) if attrs is not None else _encoding("dict", "0.1.0"))

    def _plan_object(self, adata: anndata.AnnData) -> None:
        n_obs = int(adata.n_obs)
        # In its original order (`of_anndata`): each row is its own original
        # position, and a position column left names nothing.
        obs = adata.obs
        if ORIGINAL_POSITION_COLUMN in obs.columns:
            without_positions({ORIGINAL_POSITION_COLUMN: None}, "obs")
            obs = obs.drop(columns=[ORIGINAL_POSITION_COLUMN])
        priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED)
        assert priority is not None
        levels = prefix_levels(n_obs)
        self.n_obs = n_obs
        self._plan_anndata(
            "",
            obs=obs,
            var=adata.var,
            uns=dict(adata.uns),
            x=None if adata.X is None else _as_element(adata.X),
            layers={key: _as_element(value) for key, value in adata.layers.items()},
            obsm=dict(adata.obsm),
            obsp=dict(adata.obsp),
            varm=dict(adata.varm),
            varp=dict(adata.varp),
            priority=priority,
            levels=levels,
            record=levels,
        )
        raw = adata.raw
        if raw is not None:
            self._group("raw", _encoding("raw", "0.1.0"))
            self._write("raw/var", raw.var)
            self._write("raw/varm", dict(raw.varm))
            self._matrix("raw/X", _as_element(raw.X), priority, levels=None, expression=True)

    def _plan_anndata(
        self,
        prefix: str,
        *,
        obs: Any,
        var: Any,
        uns: Mapping[str, Any],
        x: Any,
        layers: Mapping[str, Any],
        obsm: Mapping[str, Any],
        obsp: Mapping[str, Any],
        varm: Mapping[str, Any],
        varp: Mapping[str, Any],
        priority: CellPriority | None,
        levels: np.ndarray | None,
        record: np.ndarray | None,
    ) -> None:
        """One AnnData at ``prefix`` (``""``, or ``"mod/<name>/"``), its cells in ``priority``'s order.

        ``obs`` and the arrays are the source's, in its own order; ``levels``
        are the prefix levels of its rows (a modality's: its rows below each
        global level); ``record``, when given, the levels ``uns['cellucid']['priority']``
        records (the dataset's own, at the root). ``priority`` None keeps the
        rows as stored (a file in priority order already).
        """
        n_obs = len(obs)
        self._groups[prefix] = n_obs
        if priority is not None:
            obs = obs.iloc[priority.order].copy()
            obs[ORIGINAL_POSITION_COLUMN] = priority.original.astype(np.uint32)
            uns = _permute_quantiles(uns, priority.order, n_obs)
            if record is not None:
                uns = self._with_priority(uns, priority, record)
        self._write(f"{prefix}obs", obs)
        self._write(f"{prefix}var", var)
        self._write(f"{prefix}uns", uns)
        counts = self._matrix(f"{prefix}X", x, priority, levels=levels, expression=True)
        self._group(f"{prefix}layers")
        for key, value in layers.items():
            self._matrix(f"{prefix}layers/{key}", value, priority, levels=None, expression=True)
        self._per_cell(f"{prefix}obsm", obsm, priority)
        self._graphs(f"{prefix}obsp", obsp, priority)
        # Prefix counts the source had are for its own order.
        self._write(f"{prefix}varm", _without(varm, PREFIX_COUNTS_KEY) if priority is not None else dict(varm))
        if counts is not None:
            assert levels is not None
            self._plan_counts(f"{prefix}varm/{PREFIX_COUNTS_KEY}", counts, len(var), levels)
        self._write(f"{prefix}varp", dict(varp))

    def _plan_file(self, root: Any) -> None:
        if _encoding_of(root) != "anndata":
            raise ValueError(
                "This AnnData store predates anndata's current encoding (anndata 0.8): "
                "write it again with a current anndata to serve it."
            )
        n_obs = _rows_of(root)
        if n_obs <= 0:
            raise ValueError("AnnData datasets must contain at least one observation")
        self.n_obs = n_obs
        own = _verified_priority(root, warn=True)
        priority: CellPriority | None = None
        levels: np.ndarray | None = None
        if own is not None:
            # In a verified priority order already: kept as stored.
            self.seed = own
        else:
            # Each row its own original position (a position column the file
            # holds is replaced: it names nothing).
            priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED)
            assert priority is not None
            levels = prefix_levels(n_obs)
        self._plan_file_group(root, "", priority, levels, record=levels)

    def _plan_file_group(
        self, root: Any, prefix: str, priority: CellPriority | None, levels: np.ndarray | None, *, record: np.ndarray | None
    ) -> None:
        """One AnnData group of a file (the root, or a MuData modality) at ``prefix``."""
        from anndata.io import read_elem

        def elements(key: str) -> dict[str, Any]:
            group = root.get(key, {})
            return {name: group[name] for name in group}

        self._plan_anndata(
            prefix,
            obs=read_elem(root["obs"]),
            var=read_elem(root["var"]),
            uns=read_elem(root["uns"]) if "uns" in root else {},
            x=root.get("X"),
            layers=elements("layers"),
            obsm=elements("obsm"),
            obsp=elements("obsp"),
            varm=read_elem(root["varm"]) if "varm" in root else {},
            varp=read_elem(root["varp"]) if "varp" in root else {},
            priority=priority,
            levels=levels,
            record=record,
        )
        if "raw" in root and _is_group(root["raw"]):
            raw = root["raw"]
            self._group(f"{prefix}raw", dict(raw.attrs))
            for key in raw:
                if key == "X":
                    self._matrix(f"{prefix}raw/X", raw["X"], priority, levels=None, expression=True)
                else:
                    self._write(f"{prefix}raw/{key}", read_elem(raw[key]))

    def _plan_presented(self, presented: Presented | PresentedModalities) -> None:
        """A source the browser presents through a layout view (Loom, Cell Ranger), as it presents it."""
        from ._served_formats import Presented

        n_obs = len(presented.obs)
        if n_obs <= 0:
            raise ValueError("Datasets must contain at least one observation")
        priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED)
        assert priority is not None
        levels = prefix_levels(n_obs)
        self.n_obs = n_obs
        if isinstance(presented, Presented):
            self._plan_anndata(
                "",
                obs=presented.obs,
                var=presented.var,
                uns=presented.uns,
                x=presented.x,
                layers=presented.layers,
                obsm=presented.obsm,
                obsp=presented.obsp,
                varm={},
                varp={},
                priority=priority,
                levels=levels,
                record=levels,
            )
            return
        # Several feature types over the same cells: a MuData, every modality measuring every cell.
        names = list(presented.modalities)
        obs = presented.obs.iloc[priority.order].copy()
        obs[ORIGINAL_POSITION_COLUMN] = priority.original.astype(np.uint32)
        self._mudata_root(names)
        self._groups[""] = n_obs
        self._write("obs", obs)
        self._write("uns", self._with_priority({}, priority, levels))
        self._per_cell("obsm", presented.obsm, priority)
        self._group("obsp")
        self._write("var", pd.concat([modality.var for modality in presented.modalities.values()]))
        self._write("varm", {})
        self._write("varp", {})
        self._write("obsmap", {name: np.arange(1, n_obs + 1, dtype=np.uint32) for name in names})
        offsets = np.cumsum([0, *[len(modality.var) for modality in presented.modalities.values()]])
        self._write(
            "varmap",
            {
                name: np.where(
                    (np.arange(offsets[-1]) >= offsets[i]) & (np.arange(offsets[-1]) < offsets[i + 1]),
                    np.arange(offsets[-1]) - offsets[i] + 1,
                    0,
                ).astype(np.uint32)
                for i, name in enumerate(names)
            },
        )
        for name, modality in presented.modalities.items():
            self._group(f"mod/{name}", _encoding("anndata", "0.1.0"))
            self._plan_anndata(
                f"mod/{name}/",
                obs=modality.obs.iloc[:, :0],
                var=modality.var,
                uns={},
                x=modality.x,
                layers=modality.layers,
                obsm={},
                obsp={},
                varm={},
                varp={},
                priority=priority,
                levels=levels,
                record=None,
            )

    def _mudata_root(self, names: list[str]) -> None:
        """The root of a MuData store: shared cells (``axis`` 0), and ``mod`` in order."""
        self.layout = "mudata"
        self.root.attrs.update({"encoding-type": "MuData", "encoding-version": "0.1.0", "axis": 0})
        self._group("mod", {"mod-order": names})

    def _plan_mudata_file(self, root: Any) -> None:
        from anndata.io import read_elem

        from ._anndata_files import _modality_names

        if (_encoding_of(root) or "").lower() != "mudata" or "mod" not in root:
            raise ValueError("This MuData store has no MuData encoding or no mod group of modalities.")
        if int(np.asarray(root.attrs.get("axis", 0))) != 0:
            raise ValueError(
                "This MuData shares features, not cells (axis != 0); Cellucid shows cells. Write it with shared cells (axis=0)."
            )
        obs = read_elem(root["obs"])
        n_obs = len(obs)
        if n_obs <= 0:
            raise ValueError("MuData datasets must contain at least one observation")
        names = _modality_names(root["mod"])
        if not names:
            raise ValueError("The MuData store holds no modality.")
        self.n_obs = n_obs
        own = _verified_priority(root, warn=True)
        priority: CellPriority | None = None
        levels: np.ndarray | None = None
        if own is not None:
            # Served through the store because the browser cannot decode a codec: kept in its own order.
            self.seed = own
        else:
            # Each row its own original position (a position column replaced).
            priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED)
            assert priority is not None
            levels = prefix_levels(n_obs)
        self._mudata_root(names)
        root_attrs = {key: root.attrs[key] for key in root.attrs if key not in ("encoding-type", "encoding-version", "axis")}
        if root_attrs:
            self.root.attrs.update({key: _attribute_value(value) for key, value in root_attrs.items()})
        uns = read_elem(root["uns"]) if "uns" in root else {}
        # A modality's cell map; one with none measured every cell, in order
        # (as the browser reads it: ``views/modalities.js``).
        stored_maps = root.get("obsmap")
        obsmaps: dict[str, np.ndarray] = {}
        for name in names:
            if stored_maps is not None and name in stored_maps:
                obsmaps[name] = np.asarray(stored_maps[name][...])
                continue
            rows = _rows_of(root["mod"][name])
            if rows != n_obs:
                raise ValueError(
                    f"Modality {name!r} has {rows:,} {plural(rows, 'cell')} and the dataset {n_obs:,}, and no cell map "
                    "(obsmap) says which they are."
                )
            obsmaps[name] = np.arange(1, n_obs + 1, dtype=np.uint32)
        if priority is not None:
            obs = obs.iloc[priority.order].copy()
            obs[ORIGINAL_POSITION_COLUMN] = priority.original.astype(np.uint32)
            uns = self._with_priority(_permute_quantiles(uns, priority.order, n_obs), priority, levels)  # type: ignore[arg-type]
        self._groups[""] = n_obs
        self._write("obs", obs)
        self._write("uns", uns)
        self._write("var", read_elem(root["var"]))
        obsm = root.get("obsm", {})
        self._per_cell("obsm", {key: obsm[key] for key in obsm}, priority)
        obsp = root.get("obsp", {})
        self._graphs("obsp", {key: obsp[key] for key in obsp}, priority)
        self._write("varm", read_elem(root["varm"]) if "varm" in root else {})
        self._write("varp", read_elem(root["varp"]) if "varp" in root else {})
        self._write("varmap", read_elem(root["varmap"]) if "varmap" in root else {})
        served_maps: dict[str, np.ndarray] = {}
        for name in names:
            group = root["mod"][name]
            if (_encoding_of(group) or "") != "anndata":
                raise ValueError(
                    f"Modality {name!r} predates anndata's current encoding (anndata 0.8): write the MuData again "
                    "with a current mudata to serve it."
                )
            self._group(f"mod/{name}", {key: _attribute_value(group.attrs[key]) for key in group.attrs})
            if priority is None:
                served_maps[name] = obsmaps[name]
                self._plan_file_group(group, f"mod/{name}/", None, None, record=None)
                continue
            modality, served_map, below = _modality_order(priority, obsmaps[name], _rows_of(group))
            served_maps[name] = served_map
            assert levels is not None
            self._plan_file_group(group, f"mod/{name}/", modality, below[levels], record=None)
        self._write("obsmap", served_maps)

    def _with_priority(self, uns: Mapping[str, Any], priority: CellPriority, levels: np.ndarray) -> dict[str, Any]:
        """``uns`` with ``uns['cellucid']['priority']`` of this store (replacing one the source had)."""
        cellucid = uns.get("cellucid")
        own = dict(cellucid) if isinstance(cellucid, Mapping) else {}
        own["priority"] = priority.record(levels)
        return {**uns, "cellucid": own}

    def _put(self, key: str, document: Any) -> None:
        from zarr.core.buffer import default_buffer_prototype

        self._documents[key] = default_buffer_prototype().buffer.from_bytes(json.dumps(document).encode())

    def _plan(self, path: str, array: _Planned) -> None:
        self._put(f"{path}/.zarray", array.metadata())
        if array.attrs:
            self._put(f"{path}/.zattrs", array.attrs)
        self._planned[path] = array

    def _sparse_group(self, path: str, form: str, shape: tuple[int, int]) -> None:
        self._put(f"{path}/.zgroup", {"zarr_format": 2})
        self._put(f"{path}/.zattrs", {**_encoding(f"{form}_matrix", "0.1.0"), "shape": list(shape)})

    def _plan_sparse(self, path: str, served: _ServedSparse, chunk: int) -> None:
        self._sparse_group(path, served.form, served.shape)
        self._plan(f"{path}/indptr", served.pointers)
        self._plan(f"{path}/indices", _Entries(served, "indices", chunk))
        self._plan(f"{path}/data", _Entries(served, "data", chunk))

    def _matrix(
        self, path: str, value: Any, priority: CellPriority | None, *, levels: np.ndarray | None, expression: bool
    ) -> Callable[[], np.ndarray] | None:
        """Plan one per-cell matrix at ``path``.

        ``expression``: a cells x genes matrix (``X``, a layer, ``raw/X``),
        presented gene-major when it is a cell-major matrix in memory.
        Returns how its prefix counts are made when ``levels`` are given and
        it is served gene-major and sparse.
        """
        if value is None:
            return None
        order = None if priority is None else priority.order
        from ._served_formats import CellRangerX

        if isinstance(value, CellRangerX):
            # Cell-major, as stored: rows in priority order, a feature type's columns as a window.
            matrix = value.matrix
            group = matrix.group
            whole = value.start == 0 and value.count == matrix.n_features
            served = _ServedSparse(
                _Sparse(
                    form="csr",
                    shape=(matrix.n_cells, matrix.n_features),
                    indptr=lambda: np.asarray(group["indptr"][...]),
                    indices=_Vector(group["indices"], self._folder),
                    data=_Vector(group["data"], self._folder),
                ),
                major=order,
                minor=None,
                cache=self._lines,
                window=None if whole else (value.start, value.start + value.count),
            )
            self._plan_sparse(path, served, SLICED_CHUNK_ELEMENTS if expression else WHOLE_CHUNK_ELEMENTS)
            return None
        if sparse.issparse(value) and expression and priority is not None:
            # In memory: a snapshot, stored as prepare() stores it.
            return self._gene_major(path, value, priority, levels)
        if sparse.issparse(value) and priority is not None:
            # A graph or an obsm matrix in memory: a snapshot, canonical.
            value = sparse.csr_matrix(value, copy=True) if value.format == "csr" else sparse.csc_matrix(value, copy=True)
            value.sum_duplicates()
            value.eliminate_zeros()
        if sparse.issparse(value) or _is_sparse_group(value):
            source = _sparse_of_memory(value) if sparse.issparse(value) else _sparse_of_group(value, self._folder)
            renumber = priority is not None and source.form == "csc"
            served = _ServedSparse(
                source,
                major=order if source.form == "csr" else None,
                minor=priority.new_row if renumber and priority is not None else None,
                cache=self._lines,
            )
            self._plan_sparse(path, served, SLICED_CHUNK_ELEMENTS if expression else WHOLE_CHUNK_ELEMENTS)
            if levels is not None and source.form == "csc":
                return lambda: served.prefix_counts(levels)
            return None
        if not _is_dense_array(value):
            raise ValueError(f"{path} is not a matrix Cellucid serves ({type(value).__name__}).")
        if isinstance(value, np.ndarray) and priority is not None:
            # In memory: a snapshot in priority order (rows as served).
            value = self._dense_snapshot(path, value, priority, expression=expression)
            order = None
        shape = tuple(int(extent) for extent in value.shape)
        chunks: tuple[int, ...]
        if len(shape) == 2 and expression and _reads_whole_rows(value, _strip_chunks((shape[0], shape[1]))):
            laid = _DenseGeneMajor(value, priority, _strip_chunks((shape[0], shape[1])), self._folder.path, path.replace("/", "-"))
            self._layouts.append(laid)
            self._plan(path, laid)
            return None
        if len(shape) == 2:
            chunks = _strip_chunks((shape[0], shape[1])) if expression else _row_chunks((shape[0], shape[1]))
        else:
            # One value per cell.
            chunks = _cell_chunks(shape[0])
        self._plan(path, _DenseRows(value, order, chunks, self._folder))
        return None

    def _gene_major(
        self, path: str, matrix: Any, priority: CellPriority, levels: np.ndarray | None
    ) -> Callable[[], np.ndarray] | None:
        """An expression matrix in memory, laid out now (its snapshot): ``_GeneMajor``."""
        laid = _GeneMajor(matrix, priority, self._folder.path, path.replace("/", "-"), self._stop, levels, path)
        self._layouts.append(laid)
        n_obs, n_vars = laid.shape
        self._sparse_group(path, "csc", (n_obs, n_vars))
        self._plan(
            f"{path}/indptr",
            _Values(
                (n_vars + 1,),
                np.int64 if laid.nnz >= 2**31 else np.int32,
                _vector_chunks(n_vars + 1, WHOLE_CHUNK_ELEMENTS),
                lambda: laid.laid_out()[0],
            ),
        )
        self._plan(f"{path}/indices", _GeneMajorEntries(laid, "indices"))
        self._plan(f"{path}/data", _GeneMajorEntries(laid, "data"))
        return None if levels is None else laid.counts

    def _dense_snapshot(self, path: str, value: np.ndarray, priority: CellPriority | None, *, expression: bool) -> np.ndarray:
        """A dense array held in memory, copied now (its snapshot) in priority
        order into a memory-mapped file: an expression matrix in the type
        ``prepare()`` stores it (a non-finite value refused), any other as it is."""
        from .prepare_data._inputs import expression_dtype, require_expression_values

        array = np.asarray(value)
        dtype = array.dtype
        if expression:
            require_expression_values(array, path)
            dtype = np.dtype(expression_dtype(array))
        snapshot = np.lib.format.open_memmap(
            self._folder.path() / f"{path.replace('/', '-')}-dense.npy", mode="w+", dtype=dtype, shape=array.shape
        )
        order = None if priority is None else priority.order
        step = max(1, _PASS_ENTRIES // max(1, int(np.prod(array.shape[1:], dtype=np.int64))))
        for start in range(0, array.shape[0], step):
            rows = slice(start, start + step) if order is None else order[start : start + step]
            snapshot[start : start + step] = array[rows]
        snapshot.flush()
        self._snapshots.append(snapshot)
        return snapshot

    def _plan_counts(self, path: str, make: Callable[[], np.ndarray], n_vars: int, levels: np.ndarray) -> None:
        """``varm['cellucid_prefix_counts']``: one column per level, each read whole (made on first read)."""
        counts = _Values((n_vars, len(levels)), np.uint32, (max(1, min(n_vars, WHOLE_CHUNK_ELEMENTS)), 1), make)
        counts.attrs = _encoding("array", "0.2.0")
        self._plan(path, counts)

    def _per_cell(self, key: str, elements: Mapping[str, Any], priority: CellPriority | None) -> None:
        """``obsm``: dense arrays and sparse matrices planned, anything else (a table) written in priority order."""
        from anndata.io import read_elem

        self._group(key)
        order = None if priority is None else priority.order
        for name, value in elements.items():
            value = _as_element(value)
            path = f"{key}/{name}"
            if sparse.issparse(value) or _is_sparse_group(value) or _is_dense_array(value) and len(value.shape) == 2:
                self._matrix(path, value, priority, levels=None, expression=False)
                continue
            if _is_group(value) or hasattr(value, "attrs") and hasattr(value, "dtype") and not isinstance(value, np.ndarray):
                value = read_elem(value)
            self._write(path, value if order is None else _take_rows(value, order))

    def _graphs(self, key: str, elements: Mapping[str, Any], priority: CellPriority | None) -> None:
        """``obsp``: each graph with both axes in priority order."""
        from anndata.io import read_elem

        self._group(key)
        for name, value in elements.items():
            value = _as_element(value)
            path = f"{key}/{name}"
            if sparse.issparse(value) or _is_sparse_group(value):
                source = _sparse_of_memory(value) if sparse.issparse(value) else _sparse_of_group(value, self._folder)
                served = _ServedSparse(
                    source,
                    major=None if priority is None else priority.order,
                    minor=None if priority is None else priority.new_row,
                    cache=self._lines,
                )
                self._plan_sparse(path, served, WHOLE_CHUNK_ELEMENTS)
                continue
            dense = np.asarray(value if isinstance(value, np.ndarray) else read_elem(value))
            self._write(path, dense if priority is None else dense[np.ix_(priority.order, priority.order)])

    def _settle(self) -> None:
        """Every array anndata wrote uncompressed (its scalars carry zarr's default
        codec), and a numeric one with a fill value."""
        import numcodecs

        for key in [key for key in self._documents if key.endswith("/.zarray") or key == ".zarray"]:
            document = json.loads(self._documents[key].to_bytes())
            if document.get("fill_value") is None and document.get("dtype") != "|O":
                # A numeric array anndata wrote with no fill value (a scalar):
                # zero, as prepare() writes one (every chunk is stored).
                document["fill_value"] = False if document["dtype"].endswith("b1") else 0
                self._put(key, document)
            if document.get("compressor") is None:
                continue
            codec = numcodecs.get_codec(document["compressor"])
            prefix = key[: -len(".zarray")]
            for chunk in [name for name in self._documents if name.startswith(prefix) and "/" not in name[len(prefix) :] and not name[len(prefix) :].startswith(".")]:
                from zarr.core.buffer import default_buffer_prototype

                raw = codec.decode(self._documents[chunk].to_bytes())
                self._documents[chunk] = default_buffer_prototype().buffer.from_bytes(bytes(raw))
            document["compressor"] = None
            self._put(key, document)
        self._chunk_as_prepare()

    def _per_cell_rows(self, path: str) -> int | None:
        """The rows of the AnnData group whose per-cell element ``path`` is in (``obs``, ``obsm``,
        the stored outlier quantiles, MuData's cell maps), else None."""
        if path.startswith("obsmap/"):
            return self._groups.get("")
        for prefix, rows in self._groups.items():
            for area in ("obs/", "obsm/", "uns/cellucid/outlierQuantiles/"):
                if path.startswith(f"{prefix}{area}"):
                    return rows
        return None

    def _chunk_as_prepare(self) -> None:
        """Every array anndata wrote, chunked as ``prepare()`` chunks it (one
        place: :mod:`cellucid.prepare_data._store`): an array of one value or row
        per cell by :data:`~cellucid.prepare_data._store.CELL_CHUNK_ELEMENTS`
        (a subsample reads the first rows), any other vector by
        :data:`~cellucid.prepare_data._store.WHOLE_CHUNK_ELEMENTS` (read whole).
        Arrays the store plans are made in those chunks already."""
        for key in [key for key in self._documents if key.endswith("/.zarray")]:
            path = key[: -len("/.zarray")]
            if path in self._planned:
                continue
            document = json.loads(self._documents[key].to_bytes())
            shape = document["shape"]
            if not shape:
                continue
            rows = self._per_cell_rows(path)
            if rows is not None and shape[0] == rows:
                chunks: tuple[int, ...] = _cell_chunks(rows) if len(shape) == 1 else _row_chunks((shape[0], shape[1]))
            elif len(shape) == 1:
                chunks = _vector_chunks(shape[0], WHOLE_CHUNK_ELEMENTS)
            else:
                continue
            if list(chunks) != document["chunks"]:
                self._rechunk(path, document, chunks)

    def _rechunk(self, path: str, document: dict[str, Any], chunks: tuple[int, ...]) -> None:
        """Lay one uncompressed array anndata wrote out in ``chunks`` (along its first axis)."""
        import numcodecs
        import zarr
        from zarr.core.buffer import default_buffer_prototype

        values = np.asarray(zarr.open_array(self._memory, path=path, mode="r", zarr_format=2)[...])
        for name in [name for name in self._documents if name.startswith(f"{path}/") and "/" not in name[len(path) + 1 :] and not name[len(path) + 1 :].startswith(".")]:
            del self._documents[name]
        document = {**document, "chunks": list(chunks)}
        self._put(f"{path}/.zarray", document)
        # A text array (dtype "|O", vlen-utf8), whatever numpy string type zarr reads it as.
        text = any(item.get("id") == "vlen-utf8" for item in document.get("filters") or ())
        if text:
            values = values.astype(object)
        separator = document.get("dimension_separator", ".")
        size = chunks[0]
        for index in range(-(-values.shape[0] // size)):
            block = values[index * size : (index + 1) * size]
            if len(block) < size:
                # A chunk at the array's edge is stored whole.
                padded = np.empty((size, *values.shape[1:]), dtype=values.dtype)
                padded[...] = "" if text else np.zeros((), dtype=values.dtype)
                padded[: len(block)] = block
                block = padded
            raw = numcodecs.VLenUTF8().encode(np.ascontiguousarray(block, dtype=object)) if text else np.ascontiguousarray(block).tobytes()
            name = separator.join([str(index), *(["0"] * (values.ndim - 1))])
            self._documents[f"{path}/{name}"] = default_buffer_prototype().buffer.from_bytes(bytes(raw))

    # -- serving -------------------------------------------------------------

    def _consolidated_metadata(self) -> bytes:
        if self._consolidated is None:
            metadata = {}
            for key in sorted(self._documents):
                if key.rsplit("/", 1)[-1] in (".zgroup", ".zattrs", ".zarray"):
                    metadata[key] = json.loads(self._documents[key].to_bytes())
            self._consolidated = json.dumps({"zarr_consolidated_format": 1, "metadata": metadata}).encode()
        return self._consolidated

    def get(self, key: str) -> bytes | None:
        """One key of the store, or None when it has none."""
        if key == ".zmetadata":
            return self._consolidated_metadata()
        found = self._documents.get(key)
        if found is not None:
            return bytes(found.to_bytes())
        path, _, name = key.rpartition("/")
        array = self._planned.get(path)
        if array is None or not name:
            return None
        try:
            index = tuple(int(part) for part in name.split("."))
        except ValueError:
            return None
        return array.chunk_bytes(index)

    def contents(self) -> AnnDataContents:
        """What the store holds, from its metadata (as for a store on disk)."""
        import zarr

        from ._anndata_files import _contents_of_root, _mudata_contents

        root = zarr.open_group(self._memory, mode="r", zarr_format=2)
        return _mudata_contents(root) if self.layout == "mudata" else _contents_of_root(root)

    def stored_keys(self) -> Iterator[str]:
        """Every key that has a value: documents, elements written, and planned chunks."""
        yield ".zmetadata"
        yield from list(self._documents)
        for path, array in self._planned.items():
            grid = [range(-(-extent // size)) for extent, size in zip(array.shape, array.chunks, strict=True)]
            for index in _product(grid):
                yield f"{path}/{'.'.join(str(i) for i in index)}"

    @staticmethod
    def _release(closers: list[Callable[[], None]], stop: threading.Event) -> None:
        stop.set()
        while closers:
            closer = closers.pop()
            with contextlib.suppress(Exception):
                closer()

    def close(self) -> None:
        """Close the source file and remove the files this store wrote."""
        self._finalizer()


def _product(grid: list[range]) -> Iterator[tuple[int, ...]]:
    from itertools import product

    if any(len(axis) == 0 for axis in grid):
        return iter(())
    return product(*grid)


def _attribute_value(value: Any) -> Any:
    """An HDF5 or Zarr attribute as JSON holds it."""
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return [_attribute_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _modality_order(priority: CellPriority, obsmap: np.ndarray, n_rows: int) -> tuple[CellPriority, np.ndarray, np.ndarray]:
    """A MuData modality's rows in the order of its cells (``prepare()``'s
    ``CellPriority.modality``): its priority, its served cell map, and its
    rows below each served cell (``below[k]``: the modality rows of the first
    ``k`` cells). Each of its rows is its own original position (a dataset
    in a verified stored order is served as stored)."""
    placed = np.asarray(obsmap)[priority.order]
    measured = placed > 0
    mod_order = placed[measured].astype(np.int64) - 1
    if len(mod_order) != n_rows or len(np.unique(mod_order)) != n_rows:
        raise ValueError("A MuData modality's cell map must place each of its rows at exactly one cell.")
    new_row = np.empty(n_rows, dtype=np.int64)
    new_row[mod_order] = np.arange(n_rows, dtype=np.int64)
    served = np.zeros(len(placed), dtype=np.uint32)
    served[measured] = np.arange(1, n_rows + 1, dtype=np.uint32)
    below = np.concatenate([[0], np.cumsum(measured)]).astype(np.int64)
    modality = CellPriority(seed=priority.seed, order=mod_order, original=mod_order, new_row=new_row)
    return modality, served, below


def _rows_of(root: Any) -> int:
    """The rows of an AnnData group's ``obs``, from its index's length."""
    obs = root["obs"]
    index = obs.attrs.get("_index", "_index")
    index = index.decode() if isinstance(index, bytes) else str(index)
    return int(obs[index].shape[0])


def _verified_priority(root: Any, *, warn: bool) -> int | None:
    """The seed of a store in a verified priority order, None for any other.

    The one rule (:func:`~cellucid.prepare_data._priority.stored_order`,
    the web app's ``readPriorityLayout``): ``uns['cellucid']['priority']``
    with a valid seed and ``specVersion``, ``nObs`` the store's cell count,
    ``levels`` those of that count, and ``obs['cellucid_original_position']``,
    when there, the seed's order. A store whose record does not verify, or a
    position column with no record, is served as cells Cellucid did not
    write (each row its own original position); with ``warn``, the reason is
    given as ``prepare()`` gives it.
    """
    from anndata.io import read_elem

    uns = root.get("uns")
    cellucid = uns.get("cellucid") if uns is not None and _is_group(uns) else None
    group = cellucid.get("priority") if cellucid is not None and _is_group(cellucid) else None
    record = read_elem(group) if group is not None else None
    obs = root["obs"]
    column = np.asarray(read_elem(obs[ORIGINAL_POSITION_COLUMN])) if ORIGINAL_POSITION_COLUMN in obs else None
    order = stored_order(record, _rows_of(root), column)
    if warn and order.notice is not None:
        warnings.warn(order.notice, UserWarning, stacklevel=3)
    if order.positions is None:
        if warn and column is not None:
            without_positions({ORIGINAL_POSITION_COLUMN: None}, "obs")
        return None
    assert record is not None  # a verified order has its record
    return int(np.asarray(record["seed"]).reshape(-1)[0])


def file_order(source: AnnDataStore) -> str:
    """How an AnnData or MuData file stores its cells: ``"priority"`` (a
    verified ``uns['cellucid']['priority']``: :func:`_verified_priority`),
    ``"as-given"`` (anndata's or mudata's current encoding, any other order,
    a record that no longer describes the cells included), or ``"older"`` (an
    AnnData layout before anndata 0.8, which only anndata and the browser
    read)."""
    from ._anndata_files import _open_root

    with _open_root(source) as root:
        encoding = (_encoding_of(root) or "").lower()
        if encoding not in ("anndata", "mudata"):
            return "older"
        return "priority" if _verified_priority(root, warn=False) is not None else "as-given"
