"""
The gene-major layout of an expression matrix: the one pass every writer and
server shares.

``prepare()`` writes a sparse ``X`` gene-major (CSC), each gene's entries
ascending by row, its cells in priority order (:mod:`._priority`), with
``varm['cellucid_prefix_counts']``; the served store (:mod:`cellucid._served_store`)
presents an object's or a file's expression matrices the same way, so the
browser reads one gene as one range. :class:`GeneMajorLayout` makes that
layout from any sparse source (:mod:`cellucid._source_arrays`: scipy CSR or
CSC in memory, a backed AnnData's matrix, an h5py or Zarr group) in any row
order, storing the values as ``prepare()`` stores them:

- every line canonical: entries ascending by position, a position stored
  twice summed in stored order (scipy's own sort orders a coordinate's values
  by value, and a sum of three or more can then differ in its last bit from
  the one R's Matrix gives and the browser takes: a stable sort keeps stored
  order), no stored zero;
- values real and finite within the float32 range (refused otherwise, counted
  by :func:`cellucid._non_finite.describe_non_finite`), stored as float32, or
  as float64 when they are integer counts float32 cannot hold exactly (above
  2^24) -- decided on the genes stored; a value nearer zero than float32
  holds becomes no value.

Two passes, memory one block of :data:`PASS_ENTRIES` stored values beside the
result (and nothing else: no sort of the whole matrix, no 8-byte index per
value):

1. **Counting**, in the source's own stored order (sequential reads): each
   gene's stored values and what decides their type.
2. **Placing**. A cell-major (CSR) source is read by blocks of rows in the
   target order, each block's entries placed at their genes' next free
   positions, so every gene's entries land in ascending row with no sort, and
   the cursors at each prefix level are the prefix counts. A gene-major (CSC)
   source is read by blocks of genes, each gene's rows renumbered and ordered
   by one stable sort of one key (:func:`line_order`).

A gene-major source already canonical, in the target order and of the stored
type is used as it is (``prepare()`` of a Cellucid-ready file read back with
anndata copies nothing), or copied in one sequential pass when a copy is asked
for (a served snapshot). cellucid-r lays a matrix out the same way, in C
(``src/chunks.c``).
"""

from __future__ import annotations

import contextlib
import threading
import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from scipy import sparse

from .._non_finite import describe_non_finite_blocks
from .._scratch import machine_memory_bytes
from .._source_arrays import SparseSource

__all__ = ["PASS_ENTRIES", "GeneMajor", "GeneMajorLayout", "counts_below", "line_order", "position_dtype", "stored_dtype"]

#: Bytes of temporaries a stored value takes while a cell-major block is
#: placed (its row, its new place, the block's CSR and CSC copies).
_BLOCK_BYTES_PER_VALUE = 40

#: Stored values one step of a pass holds at a time, so its temporaries take
#: a 64th of the machine's memory (the share a served store keeps of what it
#: makes, :func:`cellucid._served_store.default_cache_bytes`): 6.7M on a 16 GB
#: machine. A scheduling choice from the machine's real memory, not a limit on
#: the matrix; the layout's time does not depend on it above about a million
#: values (3M cells x 10: 0.66-0.69 s for 2^20 to 2^24), and the bytes never
#: do. Every pass over a matrix steps by it: the layout here, the served
#: store's dense snapshots and layouts and its prefix counts.
PASS_ENTRIES = max(1, machine_memory_bytes() // 64 // _BLOCK_BYTES_PER_VALUE)

#: The largest integer float32 holds exactly. Integer counts above it are
#: stored as float64, so every count reads back as it was given.
FLOAT32_EXACT_INTEGER_LIMIT = 2**24

#: The largest finite float32: a larger magnitude would be stored as an infinity.
FLOAT32_MAX = float(np.finfo(np.float32).max)

#: Values tested at a time for being whole numbers.
_INTEGRAL_BLOCK = 1 << 20


def line_order(lines: np.ndarray, positions: np.ndarray, extent: int) -> np.ndarray:
    """The stable order of entries by line, then position: one stable sort of the
    single int64 key ``line * extent + position``, the order a two-key
    ``lexsort`` gives (entries of one line and position keep their stored
    order) at a thirtieth of its time (40M values: 1.2 s against 44 s)."""
    key = np.asarray(lines, dtype=np.int64) * int(extent) + np.asarray(positions, dtype=np.int64)
    return np.argsort(key, kind="stable")


def position_dtype(extent: int) -> np.dtype:
    """The type of a gene-major matrix's row positions: int32 unless it has 2^31 rows or more."""
    return np.dtype(np.int64 if extent > np.iinfo(np.int32).max else np.int32)


def counts_below(rows: np.ndarray, lengths: np.ndarray, levels: np.ndarray) -> np.ndarray:
    """For lines of ``lengths`` entries at ``rows`` (any order within a line):
    each line's entries below each level (``varm['cellucid_prefix_counts']``)."""
    n_lines, n_levels = len(lengths), len(levels)
    if n_levels == 0:
        return np.zeros((n_lines, 0), dtype=np.uint32)
    # An entry at row r counts toward level j when r < levels[j]: its first
    # level is the number of levels <= r.
    first = np.searchsorted(levels, rows, side="right")
    lines = np.repeat(np.arange(n_lines, dtype=np.int64), lengths)
    histogram = np.bincount(lines * (n_levels + 1) + first, minlength=n_lines * (n_levels + 1))
    return np.cumsum(histogram.reshape(n_lines, n_levels + 1), axis=1)[:, :n_levels].astype(np.uint32)


@dataclass(frozen=True)
class GeneMajor:
    """A cells x genes matrix gene-major: what ``prepare()`` writes as a sparse ``X``.

    ``indptr`` (int64) bounds each gene's entries, ``indices`` their rows
    (:func:`position_dtype`), ascending within each gene, ``data`` their values
    in the stored type; ``prefix_counts`` (``n_vars x len(levels)`` uint32)
    when prefix levels were given. Read as a CSC matrix (``format``, ``nnz``).
    """

    format: ClassVar[str] = "csc"

    shape: tuple[int, int]
    indptr: np.ndarray
    indices: np.ndarray
    data: np.ndarray
    prefix_counts: np.ndarray | None

    @property
    def nnz(self) -> int:
        return int(self.indptr[-1])


def _canonical(lengths: np.ndarray, minor: np.ndarray, values: np.ndarray, extent: int) -> tuple[sparse.csr_matrix, bool]:
    """Lines of entries (``lengths``, then each entry's minor position and
    value, in stored order) as a canonical CSR matrix of those lines: each
    line's entries ascending by position, a position stored twice summed in
    stored order, no stored zero; and whether they were canonical as given
    (the matrix then shares the given arrays, which are never changed)."""
    pointers = np.zeros(len(lengths) + 1, dtype=np.int64)
    np.cumsum(lengths, out=pointers[1:])
    index_dtype = np.int64 if max(extent, int(pointers[-1])) > np.iinfo(np.int32).max else np.int32
    block = sparse.csr_matrix(
        (np.asarray(values), np.asarray(minor, dtype=index_dtype), pointers.astype(index_dtype)),
        shape=(len(lengths), extent),
    )
    if block.nnz == 0:
        return block, True
    owned = False
    if not block.has_sorted_indices:
        # Each line's entries by position; a repeated position keeps its stored order.
        order = line_order(np.repeat(np.arange(len(lengths), dtype=np.int64), lengths), block.indices, extent)
        block = sparse.csr_matrix((block.data[order], block.indices[order], block.indptr), shape=block.shape)
        block.has_sorted_indices = True
        owned = True
    unchanged = not owned
    if not block.has_canonical_format:
        # Summed in stored order (scipy adds a run of one position left to right).
        block = block if owned else block.copy()
        owned = True
        block.sum_duplicates()
        unchanged = False
    if not np.all(block.data):
        block = block if owned else block.copy()
        block.eliminate_zeros()
        unchanged = False
    return block, unchanged


def _typed(block: sparse.csr_matrix, dtype: np.dtype) -> sparse.csr_matrix:
    """``block``'s values in ``dtype``; a value nearer zero than the type holds becomes no value."""
    if block.data.dtype == dtype:
        return block
    typed = sparse.csr_matrix((block.data.astype(dtype), block.indices, block.indptr), shape=block.shape)
    if not np.all(typed.data):
        typed = sparse.csr_matrix((typed.data, typed.indices.copy(), typed.indptr.copy()), shape=block.shape)
        typed.eliminate_zeros()
    return typed


def _stored_in_float32(values: np.ndarray) -> np.ndarray:
    """Which values are not zero once stored (a value nearer zero than float32
    holds is none; an integer count above 2^24, stored as float64, never is)."""
    if values.dtype.kind == "f" and values.dtype.itemsize > 4:
        with np.errstate(over="ignore", under="ignore"):
            return np.asarray(values.astype(np.float32) != 0)
    return np.asarray(values != 0)


class _Stats:
    """What decides the stored type, gathered block by block: the values'
    kind, their range, and (float64) whether every one is a whole number."""

    def __init__(self, dtype: np.dtype) -> None:
        self.kind = dtype.kind
        self.itemsize = dtype.itemsize
        self.low = np.inf
        self.high = -np.inf
        #: Every float64 value so far a whole number.
        self.whole = True

    def add(self, values: np.ndarray, low: float | None = None, high: float | None = None) -> None:
        """Values (``low``/``high``: their range, when known already)."""
        if values.size == 0 or self.kind == "b":
            return
        self.low = min(self.low, float(values.min()) if low is None else low)
        self.high = max(self.high, float(values.max()) if high is None else high)
        if self.kind == "f" and self.itemsize > 4 and self.whole:
            self.whole = all(
                bool(np.all(np.floor(part) == part))
                for part in (values[start : start + _INTEGRAL_BLOCK] for start in range(0, values.size, _INTEGRAL_BLOCK))
            )

    def dtype(self) -> np.dtype:
        """float32, unless the values are integer counts float32 cannot hold exactly."""
        if self.low == np.inf or self.kind not in "iuf":
            return np.dtype(np.float32)
        integral = self.kind in "iu" or (self.kind == "f" and self.itemsize > 4 and self.whole)
        if integral and max(-self.low, self.high) > FLOAT32_EXACT_INTEGER_LIMIT:
            return np.dtype(np.float64)
        return np.dtype(np.float32)


def stored_dtype(values: np.ndarray) -> np.dtype:
    """The type expression values are stored in: float32, unless they are
    integer counts float32 cannot hold exactly (above 2^24), float64."""
    stats = _Stats(np.dtype(values.dtype))
    stats.add(np.asarray(values).reshape(-1))
    return stats.dtype()


def _range(values: np.ndarray) -> tuple[float, float] | None:
    """The smallest and largest value (None for none, or booleans)."""
    if values.size == 0 or values.dtype.kind == "b":
        return None
    return float(values.min()), float(values.max())


def _within_float32(extent: tuple[float, float] | None) -> bool:
    """Every value finite and within the float32 range, decided from the
    smallest and largest alone (a NaN makes both NaN)."""
    if extent is None:
        return True
    low, high = extent
    return bool(np.isfinite(low) and np.isfinite(high) and max(-low, high) <= FLOAT32_MAX)


#: Allocates one output array: ``(part, dtype, length)`` -> an array to fill
#: (``part`` names it: ``"indices"`` or ``"data"``).
Allocate = Callable[[str, np.dtype, int], np.ndarray]


def _in_memory(part: str, dtype: np.dtype, length: int) -> np.ndarray:
    return np.empty(length, dtype=dtype)


class GeneMajorLayout:
    """The gene-major layout of one sparse source in a row order (module docstring).

    ``order[r]`` is the source row stored row ``r`` holds and ``new_row[s]`` the
    stored row of source row ``s`` (both None: rows as stored); ``genes`` the
    source columns stored, in their stored order (None: every one); ``levels``
    the prefix levels counted (None: no prefix counts). ``allocate`` makes each
    output array (in memory by default; a served store maps files, kept in
    :attr:`outputs` while they fill); ``copy`` asks for new arrays even where
    the source could be used as it is (a snapshot). ``as_stored`` keeps the
    values as the source stores them -- no line made canonical, no type
    changed: only moved, a gene's entries ascending by row and a row stored
    twice in stored order -- so the result's size and type are the source's,
    known before any value is read (a file a server presents).

    :meth:`run` makes it -- :meth:`count`, then :meth:`place` -- in the
    caller's thread or another one; while it runs, :attr:`counted` (the
    pointers known), :meth:`ready` (entries placed in their final positions)
    and :meth:`progress` say how far it is, and :meth:`entries_now` makes any
    range from what the placing pass has placed so far. ``stop``, when set,
    ends a run at its next block.
    """

    def __init__(
        self,
        source: SparseSource,
        *,
        order: np.ndarray | None,
        new_row: np.ndarray | None,
        label: str,
        genes: np.ndarray | None = None,
        levels: np.ndarray | None = None,
        allocate: Allocate = _in_memory,
        copy: bool = False,
        as_stored: bool = False,
        stop: threading.Event | None = None,
    ) -> None:
        dtype = np.dtype(source.data.dtype)
        if dtype.kind not in "biuf":
            raise TypeError(f"{label} must hold real numbers, got dtype {dtype}.")
        self.source = source
        self.order = order
        self.new_row = new_row
        self.label = label
        self.genes = genes
        self.levels = np.zeros(0, dtype=np.int64) if levels is None else np.asarray(levels, dtype=np.int64)
        self._count_levels = levels is not None
        self._allocate = allocate
        self._copy = copy
        self.as_stored = as_stored
        self._stop = stop
        n_obs = source.shape[0]
        n_vars = source.shape[1] if genes is None else len(genes)
        self.shape = (int(n_obs), int(n_vars))
        #: The pointers, once counted.
        self.indptr: np.ndarray | None = None
        #: The stored type: the source's (``as_stored``), else known once counted.
        self.dtype: np.dtype | None = dtype if as_stored else None
        #: The output arrays made so far, by part (``"indices"``, ``"data"``).
        self.outputs: dict[str, np.ndarray] = {}
        #: Set once the pointers are known (the counting pass is done).
        self.counted = threading.Event()
        #: Set once :meth:`run` ends, laid out or failed (:attr:`error`).
        self.finished = threading.Event()
        self.error: BaseException | None = None
        self.result: GeneMajor | None = None
        self._canonical = False
        self._per_gene: np.ndarray | None = None
        self._placed_until = 0
        #: The placing pass's progress over a cell-major source, moved on
        #: under :attr:`_turn` after each block is placed: the rows placed (in
        #: the target order) and each gene's next position.
        self._turn = threading.Lock()
        self._rows_placed = 0
        self._cursor: np.ndarray | None = None
        #: The source's pointers that pass read (a file's are read once).
        self._source_pointers: np.ndarray | None = None
        #: When that pass began; set once :meth:`place` begins, and once its
        #: progress is published.
        self._placing_began: float | None = None
        self._placing = threading.Event()
        self._progress_shared = threading.Event()
        #: What :meth:`entries_now` made before the pass ended: the rows it
        #: passed over itself, the entries it copied from the pass's output,
        #: and the seconds it took (all, and the last).
        self.request_rows = 0
        self.request_copied = 0
        self._request_spent = 0.0
        self._request_last = 0.0
        # Progress: stored values gone through by each pass, and when the work began.
        self._work = (1 if as_stored and source.form == "csc" else 2) * max(1, source.nnz)
        self._done = 0
        self._began: float | None = None
        # Source column -> stored gene (-1: not stored).
        if genes is not None:
            self._gene_of: np.ndarray | None = np.full(source.shape[1], -1, dtype=np.int64)
            self._gene_of[np.asarray(genes, dtype=np.int64)] = np.arange(len(genes), dtype=np.int64)
        else:
            self._gene_of = None

    # -- what is known, and how far it is ------------------------------------

    @property
    def nnz(self) -> int | None:
        """The stored values: the source's (``as_stored``, every gene), else known once counted."""
        if self.indptr is not None:
            return int(self.indptr[-1])
        if self.as_stored and self.genes is None:
            return self.source.nnz
        return None

    @property
    def work(self) -> int:
        """Stored values the layout goes through, over its passes (its size among other work)."""
        return self._work

    def ready(self, end: int) -> bool:
        """Whether the stored entries before ``end`` are in their final positions."""
        return self.result is not None or end <= self._placed_until

    def sooner(self) -> bool:
        """Whether the placing pass over a cell-major source will end sooner
        than :meth:`entries_now` makes its rest again: at the pace it has
        placed so far, its seconds left are no more than the last
        range made from its progress took, or than every one took together.
        The first keeps a range asked near the end from passing over the rows
        the pass is about to place; the second bounds a burst of ranges (a
        hundred genes asked at once) to twice what waiting for the pass costs,
        where each would otherwise pass over the rows left again."""
        began, rows = self._placing_began, self._rows_placed
        if self.result is not None or began is None or rows <= 0:
            return False
        working = time.monotonic() - began
        if working <= 0:
            return False
        left = working / rows * (self.shape[0] - rows)
        return left <= self._request_last or left <= self._request_spent

    def progress(self) -> tuple[float, float | None]:
        """The share of the work done, and the seconds left at the pace so far (None before it has a pace)."""
        if self.result is not None:
            return 1.0, 0.0
        share = min(1.0, self._done / self._work)
        if self._began is None or self._done == 0:
            return share, None
        elapsed = time.monotonic() - self._began
        return share, elapsed * (self._work - self._done) / self._done

    def _step(self, values: int) -> None:
        if self._stop is not None and self._stop.is_set():
            raise RuntimeError("The served store was closed while its matrix was laid out.")
        self._done += values

    # -- the run -------------------------------------------------------------

    def run(self) -> GeneMajor:
        """Lay the matrix out (counting it first, unless :meth:`count` has); the
        result (also :attr:`result`), or the error raised (also :attr:`error`)."""
        try:
            self.count()
            return self.place()
        except BaseException as error:
            self.error = error
            raise
        finally:
            self.counted.set()
            self.finished.set()

    def count(self) -> None:
        """The counting pass: each stored gene's values, and their stored type
        (which refuses a value that is not finite within the float32 range)."""
        if self.indptr is not None:
            return
        if self._began is None:
            self._began = time.monotonic()
        pointers = self._pointers()
        per_gene = self._count_as_stored(pointers) if self.as_stored else self._count_canonical(pointers)
        indptr = np.zeros(self.shape[1] + 1, dtype=np.int64)
        np.cumsum(per_gene, out=indptr[1:])
        self._per_gene = per_gene
        self.indptr = indptr
        self.counted.set()

    def entries_now(self, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        """Stored entries ``[start, end)`` -- their rows and values -- exactly
        as the placing pass will place them, made now by a layout of the genes
        they reach into alone: for a reader that asks before the placing pass
        has reached them. A gene-major source: their own values. A cell-major
        one: the entries the pass has placed are copied, and only those of the
        rows it has not reached are made, found by one scan of the source's
        column indices for the genes' columns (:meth:`_scanned`; a source
        whose lines the layout makes canonical: its rows not reached, read
        line by line); before the pass begins, every row's the same way. The
        pass goes on meanwhile: it writes only past the cursors read. After
        :meth:`count`."""
        assert self.indptr is not None and self.dtype is not None
        first = int(np.searchsorted(self.indptr, start, side="right")) - 1
        last = int(np.searchsorted(self.indptr, end, side="left"))
        stored = np.arange(first, last, dtype=np.int64)
        genes = stored if self.genes is None else np.asarray(self.genes, dtype=np.int64)[stored]
        part = GeneMajorLayout(
            self.source, order=self.order, new_row=self.new_row, label=self.label, genes=genes, as_stored=self.as_stored, stop=self._stop
        )
        # Counted already: these genes' share of the pointers, and the type
        # every gene of the matrix is stored in (what this pass places).
        base = int(self.indptr[first])
        part.indptr = self.indptr[first : last + 1] - base
        part._per_gene = np.diff(part.indptr)
        part.dtype = self.dtype
        part.counted.set()
        offset = start - base
        if self.source.form == "csc":
            placed = part.place()
            return np.asarray(placed.indices[offset : offset + end - start]), np.asarray(placed.data[offset : offset + end - start])
        began = time.monotonic()
        # A pass about to publish its progress (its outputs being made) is waited for.
        while self._placing.is_set() and not self._progress_shared.wait(0.01) and not self.finished.is_set():
            pass
        with self._turn:
            rows, cursor = self._rows_placed, None if self._cursor is None else self._cursor[first:last] - base
            result = self.result
        if result is not None:
            return np.asarray(result.indices[start:end]), np.asarray(result.data[start:end])
        if cursor is None:
            # The pass has not begun: every row, from each gene's start.
            rows, cursor = 0, part.indptr[:-1].copy()
        with self._turn:
            if self._source_pointers is None:
                self._source_pointers = self._pointers()
            pointers = self._source_pointers
        indices, data = part._outputs(int(part.indptr[-1]))
        # What the pass has placed, within the range: copied.
        copied = 0
        for low, high in zip(part.indptr[:-1].tolist(), cursor.tolist(), strict=True):
            low, high = max(low, offset), min(high, offset + end - start)
            if high > low:
                indices[low:high] = self.outputs["indices"][base + low : base + high]
                data[low:high] = self.outputs["data"][base + low : base + high]
                copied += high - low
        # The rows it has not reached: placed here, from its cursors.
        if self.as_stored or self._canonical:
            self._scanned(pointers, first, last, rows, cursor, indices, data)
        else:
            part._place_rows(pointers, part.indptr, first_row=rows, cursor=cursor, outputs=(indices, data))
        took = time.monotonic() - began
        with self._turn:
            self.request_rows += self.shape[0] - rows
            self.request_copied += copied
            self._request_spent += took
            self._request_last = took
        return np.asarray(indices[offset : offset + end - start]), np.asarray(data[offset : offset + end - start])

    def _scanned(
        self, pointers: np.ndarray, first: int, last: int, rows: int, cursor: np.ndarray, indices: np.ndarray, data: np.ndarray
    ) -> None:
        """The entries of stored genes ``first .. last - 1`` in stored rows
        ``rows`` on, of a cell-major source whose entries are stored as they
        are (``as_stored``, or every line canonical in the stored type), placed
        at each gene's ``cursor`` (relative to ``first``'s pointer) in
        ``indices``/``data``. One scan of the source's column indices, block
        by block (a comparison a stored value), finds the genes' entries; only
        those are gathered and ordered -- by gene, then stored row, a row's
        entries of one gene in stored order (one stable sort) -- where the
        placing pass reads, splits and re-sorts every value."""
        source = self.source
        n_obs = self.shape[0]
        found: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        for start, end in self._blocks(pointers, None, source.n_lines):
            low, high = int(pointers[start]), int(pointers[end])
            columns = np.asarray(source.indices.span(low, high))
            if self._gene_of is None:
                hits = np.flatnonzero((columns >= first) & (columns < last))
                genes = columns[hits].astype(np.int64) - first
            else:
                stored = self._gene_of[columns]
                hits = np.flatnonzero((stored >= first) & (stored < last))
                genes = stored[hits] - first
            if hits.size == 0:
                continue
            # Each hit's line, and its row as stored.
            row = np.searchsorted(pointers[start : end + 1], low + hits, side="right") - 1 + start
            if self.new_row is not None:
                row = self.new_row[row]
            kept = row >= rows
            if not np.all(kept):
                hits, genes, row = hits[kept], genes[kept], row[kept]
            found.append((genes, row.astype(np.int64), np.asarray(source.data.span(low, high))[hits]))
        if not found:
            return
        genes = np.concatenate([part[0] for part in found])
        row = np.concatenate([part[1] for part in found])
        values = np.concatenate([part[2] for part in found])
        order = np.argsort(genes * n_obs + row, kind="stable")
        genes, row, values = genes[order], row[order], values[order]
        per_gene = np.bincount(genes, minlength=last - first)
        starts = np.zeros(last - first, dtype=np.int64)
        np.cumsum(per_gene[:-1], out=starts[1:])
        places = np.repeat(cursor - starts, per_gene) + np.arange(genes.size, dtype=np.int64)
        indices[places] = row.astype(indices.dtype, copy=False)
        data[places] = values.astype(data.dtype, copy=False)

    def count_prefixes(self) -> np.ndarray:
        """Each stored gene's values below each level of rows in the target
        order (``varm['cellucid_prefix_counts']``), exactly as the placing pass
        records them, by a counting pass over the source alone (nothing
        placed): for a reader that asks before the placing pass has ended.
        After :meth:`count`."""
        assert self.indptr is not None and self.dtype is not None
        pointers = self._pointers()
        n_obs, n_vars = self.shape
        levels = self.levels
        counts = np.zeros((n_vars, len(levels)), dtype=np.uint32)
        if self.source.form == "csc":
            source_lines = None if self.genes is None else np.asarray(self.genes, dtype=np.int64)
            for start, end in self._blocks(pointers, source_lines, n_vars):
                lines = np.arange(start, end, dtype=np.int64) if source_lines is None else source_lines[start:end]
                block = _typed(self._block(pointers, lines)[0], self.dtype)
                rows = block.indices.astype(np.int64)
                if self.new_row is not None:
                    rows = self.new_row[rows]
                counts[start:end] = counts_below(rows, np.diff(block.indptr).astype(np.int64), levels)
            return counts
        # A cell-major source: its lines in the target order, level by level.
        level_at = {int(level): j for j, level in enumerate(levels.tolist())}
        boundaries = np.unique(np.concatenate([levels[levels < n_obs], [n_obs]])).astype(np.int64)
        cursor = np.zeros(n_vars, dtype=np.int64)
        row = 0
        for boundary in boundaries.tolist():
            span = np.arange(row, boundary, dtype=np.int64) if self.order is None else self.order[row:boundary]
            for start, end in self._blocks(pointers, span, boundary - row):
                block = _typed(self._block(pointers, span[start:end])[0], self.dtype)
                genes = block.indices if self._gene_of is None else self._gene_of[block.indices]
                cursor += np.bincount(genes[genes >= 0], minlength=n_vars)
            row = boundary
            if boundary in level_at:
                counts[:, level_at[boundary]] = cursor
        return counts

    def place(self) -> GeneMajor:
        """The placing pass (after :meth:`count`)."""
        assert self.indptr is not None and self.dtype is not None and self._per_gene is not None
        if self._began is None:
            self._began = time.monotonic()
        self._placing.set()
        pointers = self._pointers()
        indptr = self.indptr
        result: GeneMajor | None = None
        if self.source.form == "csc":
            if self._canonical and self.new_row is None and self.genes is None:
                result = self._as_stored(pointers, indptr)
            if result is None:
                result = self._place_genes(pointers, indptr)
        else:
            result = self._place_rows(pointers, indptr, shared=True)
        for array in (result.indices, result.data):
            if isinstance(array, np.memmap):
                array.flush()
        self.result = result
        return result

    def _pointers(self) -> np.ndarray:
        return np.asarray(self.source.indptr(), dtype=np.int64)

    def _blocks(self, pointers: np.ndarray, lines: np.ndarray | None, count: int) -> Iterator[tuple[int, int]]:
        """Steps ``[start, end)`` over ``count`` lines (``lines[i]`` the source line of
        step position ``i``; None: the source's own), about :data:`PASS_ENTRIES` stored values each."""
        lengths = np.diff(pointers) if lines is None else (pointers[lines + 1] - pointers[lines])
        ends = np.cumsum(lengths, dtype=np.int64)
        start = 0
        while start < count:
            base = 0 if start == 0 else int(ends[start - 1])
            end = min(count, max(start + 1, int(np.searchsorted(ends, base + PASS_ENTRIES, side="right"))))
            yield start, end
            start = end

    def _block(self, pointers: np.ndarray, lines: np.ndarray) -> tuple[sparse.csr_matrix, bool]:
        """Source lines ``lines`` as a CSR matrix of those lines: canonical, and
        whether they were canonical as given (``as_stored``: as given)."""
        lengths, minor, values = self.source.lines(lines, pointers)
        if self.as_stored:
            stored_pointers = np.zeros(len(lengths) + 1, dtype=np.int64)
            np.cumsum(lengths, out=stored_pointers[1:])
            return sparse.csr_matrix((values, minor, stored_pointers), shape=(len(lengths), self.source.extent)), True
        return _canonical(lengths, minor, values, self.source.extent)

    def _count_as_stored(self, pointers: np.ndarray) -> np.ndarray:
        """Each stored gene's values, as stored: a gene-major source's own
        pointers, a cell-major one's indices counted (one pass over them)."""
        source = self.source
        if source.form == "csc":
            lengths = np.diff(pointers)
            return lengths.copy() if self.genes is None else lengths[np.asarray(self.genes, dtype=np.int64)]
        per_gene = np.zeros(self.shape[1], dtype=np.int64)
        for start, end in self._blocks(pointers, None, source.n_lines):
            columns = np.asarray(source.indices.span(int(pointers[start]), int(pointers[end])))
            if self._gene_of is not None:
                genes = self._gene_of[columns]
                per_gene += np.bincount(genes[genes >= 0], minlength=self.shape[1])
            else:
                per_gene += np.bincount(columns, minlength=self.shape[1])
            self._step(int(pointers[end] - pointers[start]))
        return per_gene

    def _count_canonical(self, pointers: np.ndarray) -> np.ndarray:
        """Each stored gene's values once canonical, and their stored type;
        notes whether every line was canonical with nothing dropped (the
        source's own arrays then are the result's)."""
        source = self.source
        n_lines = source.n_lines
        csr = source.form == "csr"
        per_gene = np.zeros(self.shape[1], dtype=np.int64)
        stats = _Stats(np.dtype(source.data.dtype))
        float64 = stats.kind == "f" and stats.itemsize > 4
        canonical = True
        for start, end in self._blocks(pointers, None, n_lines):
            block, unchanged = self._block(pointers, np.arange(start, end, dtype=np.int64))
            extent = _range(block.data)
            if not _within_float32(extent):
                blocks = (
                    self._block(pointers, np.arange(first, last, dtype=np.int64))[0].data
                    for first, last in self._blocks(pointers, None, n_lines)
                )
                raise ValueError(
                    f"{self.label} must hold finite numbers within the float32 range. "
                    + describe_non_finite_blocks(blocks, self.label)
                )
            # Zeros are gone; a float64 value nearer zero than float32 holds is none once stored.
            stored = _stored_in_float32(block.data) if float64 else None
            canonical = canonical and unchanged and (stored is None or bool(np.all(stored)))
            lengths = np.diff(block.indptr)
            if csr:
                columns = block.indices if stored is None else block.indices[stored]
                if self._gene_of is not None:
                    genes = self._gene_of[columns]
                    per_gene += np.bincount(genes[genes >= 0], minlength=self.shape[1])
                    kept = self._gene_of[block.indices] >= 0
                    stats.add(block.data[kept])
                else:
                    per_gene += np.bincount(columns, minlength=self.shape[1])
                    stats.add(block.data, *(extent or (None, None)))
            else:
                if stored is None:
                    line_counts = lengths.astype(np.int64)
                else:
                    owner = np.repeat(np.arange(end - start, dtype=np.int64), lengths)
                    line_counts = np.bincount(owner[stored], minlength=end - start)
                if self._gene_of is not None:
                    targets = self._gene_of[start:end]
                    kept = targets >= 0
                    per_gene[targets[kept]] = line_counts[kept]
                    stats.add(block.data[np.repeat(kept, lengths)])
                else:
                    per_gene[start:end] = line_counts
                    stats.add(block.data, *(extent or (None, None)))
            self._step(int(pointers[end] - pointers[start]))
        self.dtype = stats.dtype()
        self._canonical = canonical and np.dtype(source.data.dtype) == self.dtype
        return per_gene

    def _outputs(self, nnz: int) -> tuple[np.ndarray, np.ndarray]:
        assert self.dtype is not None
        self.outputs["indices"] = self._allocate("indices", position_dtype(self.shape[0]), nnz)
        self.outputs["data"] = self._allocate("data", self.dtype, nnz)
        return self.outputs["indices"], self.outputs["data"]

    def _as_stored(self, pointers: np.ndarray, indptr: np.ndarray) -> GeneMajor | None:
        """A canonical gene-major source in the target order, of the stored type:
        its own arrays (or, with ``copy``, one sequential copy of them)."""
        source = self.source
        in_memory = source.matrix is not None and source.indices.dtype == position_dtype(self.shape[0])
        if in_memory and not self._copy:
            indices, data = np.asarray(source.matrix.indices), np.asarray(source.matrix.data)
            self.outputs.update(indices=indices, data=data)
        else:
            indices, data = self._outputs(int(indptr[-1]))
            for start, end in self._blocks(pointers, None, source.n_lines):
                first, last = int(pointers[start]), int(pointers[end])
                indices[first:last] = source.indices.span(first, last)
                data[first:last] = source.data.span(first, last)
                self._placed_until = last
                self._step(last - first)
        counts = None
        if self._count_levels:
            counts = np.zeros((self.shape[1], len(self.levels)), dtype=np.uint32)
            for start, end in self._blocks(pointers, None, source.n_lines):
                first, last = int(indptr[start]), int(indptr[end])
                counts[start:end] = counts_below(indices[first:last], np.diff(indptr[start : end + 1]), self.levels)
        return GeneMajor(self.shape, indptr, indices, data, counts)

    def _place_genes(self, pointers: np.ndarray, indptr: np.ndarray) -> GeneMajor:
        """A gene-major source, by blocks of stored genes: rows renumbered, ordered by one stable sort."""
        assert self.dtype is not None
        n_vars = self.shape[1]
        indices, data = self._outputs(int(indptr[-1]))
        counts = np.zeros((n_vars, len(self.levels)), dtype=np.uint32) if self._count_levels else None
        source_lines = None if self.genes is None else np.asarray(self.genes, dtype=np.int64)
        for start, end in self._blocks(pointers, source_lines, n_vars):
            lines = np.arange(start, end, dtype=np.int64) if source_lines is None else source_lines[start:end]
            block = _typed(self._block(pointers, lines)[0], self.dtype)
            lengths = np.diff(block.indptr).astype(np.int64)
            rows = block.indices.astype(np.int64)
            values = block.data
            if self.new_row is not None:
                rows = self.new_row[rows]
                order = line_order(np.repeat(np.arange(end - start, dtype=np.int64), lengths), rows, self.shape[0])
                rows, values = rows[order], values[order]
            first, last = int(indptr[start]), int(indptr[end])
            indices[first:last] = rows
            data[first:last] = values
            if counts is not None:
                counts[start:end] = counts_below(rows, lengths, self.levels)
            self._placed_until = last
            self._step(int((pointers[lines + 1] - pointers[lines]).sum()))
        return GeneMajor(self.shape, indptr, indices, data, counts)

    def _place_rows(
        self,
        pointers: np.ndarray,
        indptr: np.ndarray,
        *,
        shared: bool = False,
        first_row: int = 0,
        cursor: np.ndarray | None = None,
        outputs: tuple[np.ndarray, np.ndarray] | None = None,
    ) -> GeneMajor:
        """A cell-major source, by blocks of rows in the target order: each
        block's entries at their genes' next free positions (rows ascending
        within every gene, no sort); the cursors at each level are the prefix
        counts. ``shared``: the pass :meth:`entries_now` reads the progress
        of (published under :attr:`_turn` after each block). From
        ``first_row`` on, with the
        genes' ``cursor`` there, into ``outputs``: the rows a reader makes
        that the pass has not reached (no prefix counts)."""
        assert self.dtype is not None
        n_obs, n_vars = self.shape
        indices, data = self._outputs(int(indptr[-1])) if outputs is None else outputs
        cursor = indptr[:-1].copy() if cursor is None else cursor
        levels = self.levels
        counts = np.zeros((n_vars, len(levels)), dtype=np.uint32) if self._count_levels else None
        level_at = {int(level): j for j, level in enumerate(levels.tolist())}
        # Steps of rows ending at every level, so the cursors there are the prefix counts.
        boundaries = np.unique(np.concatenate([levels[levels < n_obs], [n_obs]])).astype(np.int64)
        row = first_row
        if shared:
            with self._turn:
                self._source_pointers = pointers
                self._cursor = cursor
                self._placing_began = time.monotonic()
            self._progress_shared.set()
        for boundary in boundaries[boundaries > first_row].tolist():
            span = np.arange(row, boundary, dtype=np.int64) if self.order is None else self.order[row:boundary]
            for start, end in self._blocks(pointers, span, boundary - row):
                step = self._place_block(pointers, span, row, start, end, cursor, indices, data)
                with self._turn if shared else contextlib.nullcontext():
                    cursor += step
                    self._rows_placed = row + end
            row = boundary
            if counts is not None and boundary in level_at:
                counts[:, level_at[boundary]] = cursor - indptr[:-1]
        return GeneMajor(self.shape, indptr, indices, data, counts)

    def _place_block(
        self, pointers: np.ndarray, span: np.ndarray, row: int, start: int, end: int, cursor: np.ndarray, indices: np.ndarray, data: np.ndarray
    ) -> np.ndarray:
        """Rows ``span[start:end]`` (stored rows ``row + start`` on) placed at
        their genes' cursors; how far each cursor moves past them."""
        assert self.dtype is not None
        n_vars = self.shape[1]
        lines = span[start:end]
        block = _typed(self._block(pointers, lines)[0], self.dtype)
        if self._gene_of is not None:
            genes = self._gene_of[block.indices]
            kept = genes >= 0
            lengths = np.bincount(np.repeat(np.arange(end - start, dtype=np.int64), np.diff(block.indptr))[kept], minlength=end - start)
            block = sparse.csr_matrix((block.data[kept], genes[kept], np.concatenate([[0], np.cumsum(lengths)])), shape=(end - start, n_vars))
        # The step's rows by gene, each gene's rows ascending (a counting sort, in compiled code).
        part = block.tocsc()
        step: np.ndarray = np.diff(part.indptr).astype(np.int64)
        # Each entry's place: its gene's cursor, plus its rank within the gene in the step.
        places = np.repeat(cursor - part.indptr[:-1], step)
        places += np.arange(part.nnz, dtype=np.int64)
        indices[places] = part.indices.astype(indices.dtype, copy=False) + indices.dtype.type(row + start)
        data[places] = part.data
        del places
        self._step(int((pointers[lines + 1] - pointers[lines]).sum()))
        return step
