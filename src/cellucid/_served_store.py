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
verified stored order is taken from its stored rows by the original position
each one records (one permutation, never a copy back first), and any other
source's rows are each their own original position (a position column with
no record that describes the cells names nothing, and is replaced, with the
warning ``prepare()`` gives) -- ``uns['cellucid']['priority']`` records the
seed, and a sparse gene-major ``X`` comes with
``varm['cellucid_prefix_counts']``. A subsample the viewer opens with that
seed is then the first ``k`` rows of every per-cell array, and the first
entries of every gene: the server is asked for, and sends, bytes in
proportion to ``k`` whatever the source's own layout. :attr:`PriorityStore.obs_names`
are the cells' names in their original order, taken with the rest, which
the viewer's Python side names cells by.

Nothing of a file is copied that a read can take from it in place (an
uncompressed array stored in one piece); an array stored in chunks is decoded
once into the store's scratch folder (:mod:`cellucid._scratch`) as its rows
are first read (:mod:`cellucid._source_arrays`):

- ``obs``, ``var``, ``uns``, ``varm``, ``varp`` (and ``raw``'s ``var`` and
  ``varm``) are written once, uncompressed, by anndata into a store in memory
  (``obs`` in priority order, the stored outlier quantiles with it);
- every array with one row per cell -- ``obsm``, ``obsp``, a dense matrix of
  a file -- is a plan, and each of its chunks is made when it is asked for,
  from the source:

  ======================  ==================================================
  source                  presented as
  ======================  ==================================================
  dense rows              the same array, rows in priority order
  CSC (gene-major)        CSC, each gene's row positions renumbered and
                          sorted, so its first entries are its prefix
  ``obsp`` graphs         both axes in priority order
  ======================  ==================================================

- every expression matrix (``X``, layers, ``raw/X``) that is cell-major (a
  CSR file or backed object, the usual scanpy object) or held in memory is
  **laid out gene-major** in priority order, so the browser reads a gene as
  one range (it would read a cell-major matrix's every stored value before
  its first gene), by the one layout pass ``prepare()`` writes with
  (:mod:`cellucid.prepare_data._gene_major`), into memory-mapped files beside
  the store. The layout runs in the store's **background worker**, one matrix
  after another (``X`` first): the store, its metadata and every other
  element are made at once, so the viewer opens and draws without waiting,
  and a chunk of a matrix not laid out yet is made at once from the source --
  the object's snapshot, or the file -- exactly as it will be laid out (a
  pass over a cell-major source's lines, the values of the genes it covers
  of a gene-major one): the viewer never waits for a layout. Only a file's
  matrix the worker has not counted yet (where each gene's values go is
  known once counted) is answered :class:`NotReady` -- ``503`` with
  ``Retry-After``, the time left at the pace so far. ``X`` of an object or a
  file, and every matrix of an object, are laid out from the start; a file's
  layers and ``raw/X`` when first asked for. The layout costs a pass over
  each matrix (about two seconds per 100 million stored values) and its disk
  space, never a second copy in memory.

An object held in memory is served from a **snapshot**: its expression
matrices laid out (values stored as ``prepare()`` stores them, by the same
functions: duplicates summed, stored zeros dropped, float32 unless integer
counts need float64, a non-finite value refused at once), its dense arrays
(``obsm``) copied in priority order by the same worker, its sparse
``obsm``/``obsp`` matrices copied with duplicates summed and stored zeros
dropped. While the worker lays out or copies an array of the object, that
array is **read-only** (``writeable`` off, put back when it ends): an in-place
change (``np.log1p(adata.X.data, out=adata.X.data)``) is refused by numpy
instead of mixing into what is served, and a chunk asked for meanwhile is
taken from it at once. A file's matrix is
laid out as it is stored (its values and type kept). The viewer shows the
object as it was when the store was made: replacing an element afterwards
(``adata.X = ...``) does not show. View the object again to show a change.

Chunks are uncompressed -- the viewer reads them over loopback or a notebook's
comm, where compressing on the kernel would cost more than the bytes saved --
and cut as ``prepare()`` cuts them, by the constants of
:mod:`cellucid.prepare_data._store` (per-cell arrays -- ``obs``, ``obsm``
rows, the stored outlier quantiles, MuData's cell maps -- by
``CELL_CHUNK_ELEMENTS``, anything read whole by ``WHOLE_CHUNK_ELEMENTS``;
``obs`` is cut in them straight from its columns in served order, and what
else anndata writes in memory is re-cut to them), but for an expression
matrix's entries, cut by :data:`SERVED_ENTRY_CHUNK_ELEMENTS` (a request on
this machine costs more than ``prepare()``'s chunk of them moves): the store a
server presents holds the arrays ``prepare()`` writes, value for value, as
cellucid-r's served store holds ``cellucid_prepare()``'s, chunk for chunk.
A server answers a ``Range`` of a chunk with only those bytes made
(:meth:`PriorityStore.body`).
"""

from __future__ import annotations

import contextlib
import json
import logging
import math
import numbers
import secrets
import threading
import traceback
import warnings
import weakref
from collections import OrderedDict
from collections.abc import Callable, Iterator, Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from scipy import sparse

from ._plural import plural
from ._scratch import Scratch, machine_memory_bytes
from ._source_arrays import (
    Rows,
    SparseSource,
    Vector,
    encoding_of,
    is_group,
    is_lazy_sparse,
    is_sparse_group,
    sparse_source,
)
from .prepare_data._gene_major import (
    PASS_ENTRIES,
    GeneMajorLayout,
    counts_below,
    line_order,
    position_dtype,
)
from .prepare_data._priority import (
    ORIGINAL_POSITION_COLUMN,
    PREFIX_COUNTS_KEY,
    CellPriority,
    prefix_levels,
    stored_order,
    stored_order_of,
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

__all__ = ["SERVED_ENTRY_CHUNK_ELEMENTS", "SERVED_PRIORITY_SEED", "NotReady", "PriorityStore", "StoreBody", "file_order"]

logger = logging.getLogger(__name__)

#: The seed of a served store's priority order: the subsampling
#: specification's default seed, so the viewer's default subsample of a
#: served dataset is a prefix.
SERVED_PRIORITY_SEED = 0

#: Values one chunk of a served expression matrix's ``indices`` and ``data``
#: holds (``X``, layers, ``raw/X``; 1 MiB of float32). A request on this
#: machine costs a fixed time whatever its length, at both ends -- measured on
#: an Apple M1 Pro over loopback, six connections: Python's ``http.server``
#: answers at most about 4,700 requests a second (0.21 ms of the interpreter
#: each), cellucid-r's handlers on R's thread about 1,600 (0.6 ms), and a
#: browser's ``fetch`` from a worker moves 46-79 MB/s in requests of 64 KiB
#: (Chromium, Firefox, WebKit), 300-510 MB/s in requests of 1 MiB, and no more
#: in larger ones (4 MiB: 169-590 MB/s) -- so at ``prepare()``'s 2^14 values a
#: chunk a pass over a matrix (an analysis, many genes) was bound by requests.
#: 2^18 values is the least at which every engine's loopback rate has reached
#: its plateau. A viewer reads a served chunk in
#: part (its chunks are stored uncompressed: a ``Range`` of one is only those
#: bytes, made alone, :meth:`_Planned.chunk_bytes`), so a gene, or a
#: subsample's first entries of a gene, reads no more than with smaller
#: chunks. A scheduling choice, measured; cellucid-r's served stores use the
#: same (``.chunk_elements_served``). A file ``prepare()`` writes keeps
#: :data:`~cellucid.prepare_data._store.SLICED_CHUNK_ELEMENTS`: its chunks are
#: compressed, and a chunk is the least a reader of it decodes.
SERVED_ENTRY_CHUNK_ELEMENTS = 2**18

def default_cache_bytes() -> int:
    """What one served store keeps of what it makes (bytes): a 64th of the
    machine's memory: at most the share the web app gives its decoded chunks
    (a 32nd of what its working set leaves of half the device;
    ``data/memory-budget.js``). cellucid-r's served
    stores follow the same rule (``cache_bytes``)."""
    return machine_memory_bytes() // 64


def require_cache_bytes(value: object) -> int:
    """``cache_bytes`` as a served store takes it: None for
    :func:`default_cache_bytes`, else a non-negative whole number of bytes
    (0 keeps nothing). cellucid-r's ``cache_bytes`` is checked the same way."""
    if value is None:
        return default_cache_bytes()
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < 0:
        raise ValueError(f"cache_bytes must be None or one non-negative whole number of bytes, got {value!r}.")
    return int(value)

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

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()
            self._bytes = 0




# =============================================================================
# The background worker: matrices laid out one after another
# =============================================================================


class NotReady(Exception):
    """A key of a file's matrix the worker has not counted yet: the server answers ``503``
    with ``Retry-After`` :attr:`seconds` (None: no pace measured yet, and the
    server asks again in one second, the header's resolution) and the
    message, which says how far it is."""

    def __init__(self, message: str, seconds: float | None) -> None:
        super().__init__(message)
        self.seconds = seconds


class StoreBody:
    """One key of a served store as a server answers it: its length, and its
    bytes ``[start, end)`` made when asked (``read``), so a ``Range`` of a
    chunk makes only those bytes. A key held as bytes reads slices of them."""

    __slots__ = ("read", "size")

    def __init__(self, size: int, read: Callable[[int, int], bytes]) -> None:
        self.size = size
        self.read = read

    @classmethod
    def of(cls, data: bytes) -> StoreBody:
        """Bytes held whole."""
        view = memoryview(data)
        return cls(len(data), lambda start, end: data if start == 0 and end == len(data) else bytes(view[start:end]))


class _ReadOnly:
    """The arrays of an object the worker reads, read-only while it reads them.

    Counted per array: one two matrices share (a layer that is ``X``) is
    writeable again once both are laid out; one that was read-only already is
    left alone. An in-place change meanwhile is refused by numpy ("assignment
    destination is read-only") rather than mixed into the snapshot."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._held: dict[int, tuple[np.ndarray, int]] = {}

    def hold(self, arrays: list[Any]) -> list[np.ndarray]:
        held = []
        with self._lock:
            for array in arrays:
                if not isinstance(array, np.ndarray) or isinstance(array, np.memmap):
                    continue
                entry = self._held.get(id(array))
                if entry is not None:
                    self._held[id(array)] = (entry[0], entry[1] + 1)
                elif array.flags.writeable:
                    array.flags.writeable = False
                    self._held[id(array)] = (array, 1)
                else:
                    continue
                held.append(array)
        return held

    def release(self, arrays: list[np.ndarray]) -> None:
        with self._lock:
            for array in arrays:
                entry = self._held.get(id(array))
                if entry is None:
                    continue
                if entry[1] > 1:
                    self._held[id(array)] = (entry[0], entry[1] - 1)
                    continue
                del self._held[id(array)]
                with contextlib.suppress(ValueError):
                    array.flags.writeable = True

    def release_all(self) -> None:
        with self._lock:
            held = [array for array, _ in self._held.values()]
            self._held.clear()
        for array in held:
            with contextlib.suppress(ValueError):
                array.flags.writeable = True


class _Job:
    """One matrix laid out (or copied) by the store's worker.

    ``work`` does it; ``progress`` says the share done and the seconds left
    at its pace (None before it has one); ``values`` sizes it among the jobs
    ahead of a key that waits; ``held`` are the object's arrays it reads,
    writeable again (``protection``) once it ends, done or failed."""

    def __init__(
        self,
        label: str,
        values: int,
        work: Callable[[], object],
        progress: Callable[[], tuple[float, float | None]],
        *,
        held: list[np.ndarray] | None = None,
        protection: _ReadOnly | None = None,
    ) -> None:
        self.label = label
        self.values = max(1, int(values))
        self._work = work
        self.progress = progress
        self._held = held or []
        self._protection = protection
        self.queued = False
        self.started = threading.Event()
        self.finished = threading.Event()
        self.error: BaseException | None = None
        self.seconds: float | None = None

    def run(self) -> None:
        import time

        self.started.set()
        began = time.monotonic()
        try:
            self._work()
        except BaseException as error:
            # Its frames may hold an array mapping a file of the store, which
            # Windows would keep from being removed while the error is alive.
            traceback.clear_frames(error.__traceback__)
            self.error = error
        finally:
            self.seconds = time.monotonic() - began
            if self._protection is not None:
                self._protection.release(self._held)
            self._held = []
            self.finished.set()

    def require(self) -> None:
        """Raise the error it ended with, as the reason a key cannot be answered."""
        if self.error is not None:
            raise RuntimeError(f"{self.label} could not be laid out for the viewer: {self.error}") from self.error


class _Worker:
    """The store's background worker: one daemon thread running its jobs in
    the order queued (``X`` first). Nothing here refers to the store, so
    dropping the store closes it (:meth:`PriorityStore._release`)."""

    def __init__(self, stop: threading.Event) -> None:
        self._stop = stop
        self._ready = threading.Condition()
        self._queue: list[_Job] = []
        self._thread: threading.Thread | None = None
        #: Stored values a second, from the jobs done so far.
        self._pace: float | None = None

    def queue(self, job: _Job) -> None:
        with self._ready:
            if job.queued or self._stop.is_set():
                return
            job.queued = True
            self._queue.append(job)
            if self._thread is None:
                self._thread = threading.Thread(target=self._loop, name="cellucid-layout", daemon=True)
                self._thread.start()
            self._ready.notify()

    def _loop(self) -> None:
        while True:
            with self._ready:
                while not self._queue and not self._stop.is_set():
                    self._ready.wait()
                if self._stop.is_set():
                    return
                job = self._queue[0]
            job.run()
            with self._ready:
                self._queue.remove(job)
                if job.error is None and job.seconds:
                    self._pace = job.values / job.seconds
            if job.error is not None and not self._stop.is_set():
                logger.warning("%s could not be laid out for the viewer: %s", job.label, job.error)

    def not_ready(self, job: _Job) -> NotReady:
        """Why ``job``'s key waits: how far it is, and the seconds left for it
        and the jobs ahead of it (queued now, if it was not)."""
        self.queue(job)
        with self._ready:
            ahead = list(self._queue[: self._queue.index(job)]) if job in self._queue else []
        share, left = job.progress()
        running = next((other for other in ahead if other.started.is_set()), None)
        pace = self._pace
        if running is not None:
            done, running_left = running.progress()
            if running_left is not None and done < 1:
                pace = running.values * (1 - done) / max(running_left, 1e-9)
        seconds: float | None = 0.0
        for other in ahead:
            other_share, other_left = other.progress()
            if other.started.is_set() and other_left is not None:
                seconds = None if seconds is None else seconds + other_left
            else:
                seconds = None if seconds is None or pace is None else seconds + other.values / pace
        if job.started.is_set() and left is not None:
            seconds = None if seconds is None else seconds + left
        else:
            seconds = None if seconds is None or pace is None else seconds + job.values * (1 - share) / pace
        if job.started.is_set():
            state = f"{share:.0%} done"
        elif ahead:
            state = f"after {len(ahead)} other {plural(len(ahead), 'matrix', 'matrices')}"
        else:
            state = "starting"
        when = "" if seconds is None else f", about {max(1, math.ceil(seconds)):,} s left"
        return NotReady(f"{job.label} is being laid out gene-major for the viewer ({state}{when}).", seconds)

    def wait(self, jobs: list[_Job], timeout: float | None) -> None:
        """Queue every job and wait until each has ended."""
        for job in jobs:
            self.queue(job)
        for job in jobs:
            if not job.finished.wait(timeout):
                raise TimeoutError(f"{job.label} was not laid out within {timeout} s.")
            job.require()

    def close(self) -> None:
        """Stop after the running job's current block (the store is closing)."""
        with self._ready:
            self._ready.notify_all()
            thread = self._thread
        if thread is not None and thread is not threading.current_thread():
            thread.join()
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

    def in_grid(self, index: tuple[int, ...]) -> bool:
        """Whether ``index`` names a chunk of the array."""
        grid = [-(-extent // size) for extent, size in zip(self.shape, self.chunks, strict=True)]
        return len(index) == len(grid) and all(0 <= i < n for i, n in zip(index, grid, strict=True))

    def chunk_length(self) -> int:
        """The bytes of every chunk as Zarr v2 stores it (padded to the whole chunk)."""
        return int(np.prod(self.chunks, dtype=np.int64)) * self.dtype.itemsize

    def block_part(self, index: tuple[int, ...], first: int, last: int) -> np.ndarray:
        """Values ``[first, last)`` of chunk ``index`` of a vector (offsets in
        the chunk, within the array). An array that can make part of a chunk
        alone makes only those."""
        return np.asarray(self.block(index))[first:last]

    def chunk_bytes(self, index: tuple[int, ...], start: int = 0, end: int | None = None) -> bytes:
        """Bytes ``[start, end)`` of chunk ``index`` (``in_grid``) as Zarr v2
        stores it (padded to the whole chunk; ``end`` None: to its end). A
        vector's are made alone, from the values they cover (a ``Range`` of a
        large chunk is answered without making the chunk); any other array's
        chunk is made whole."""
        size = self.chunk_length()
        end = size if end is None else end
        if not 0 <= start <= end <= size:
            raise ValueError(f"No bytes [{start}, {end}) in a chunk of {size}.")
        if len(self.chunks) == 1:
            item = self.dtype.itemsize
            first, last = start // item, -(-end // item)
            values = np.zeros(last - first, dtype=self.dtype)
            held = max(0, min(last, self.shape[0] - index[0] * self.chunks[0]) - first)
            if held > 0:
                values[:held] = self.block_part(index, first, first + held)
            return values.view(np.uint8)[start - first * item : end - first * item].tobytes()
        values = np.ascontiguousarray(self.block(index), dtype=self.dtype)
        if values.shape != tuple(self.chunks):
            padded = np.zeros(self.chunks, dtype=self.dtype)
            padded[tuple(slice(0, extent) for extent in values.shape)] = values
            values = padded
        return values.reshape(-1).view(np.uint8)[start:end].tobytes()


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

    def __init__(self, array: Any, order: np.ndarray | None, chunks: tuple[int, ...], scratch: Scratch, *, encoded: bool = True) -> None:
        self.array = array
        self.source = Rows(array, scratch)
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



class _DenseSnapshot(_Planned):
    """A dense per-cell array of an object held in memory (an ``obsm`` array)
    as served: its rows in priority order, copied once -- its snapshot -- by
    the store's worker into a memory-mapped file, the object's array
    read-only while it is copied (an in-place change is refused by numpy
    rather than mixed in). Until then a chunk is taken from that array at
    once, the rows it will hold: the store is made, and the viewer draws,
    without waiting for the copy (10M cells of a 30-column latent space: 1.2
    GB)."""

    def __init__(
        self,
        array: np.ndarray,
        order: np.ndarray | None,
        chunks: tuple[int, ...],
        scratch: Scratch,
        name: str,
        stop: threading.Event,
    ) -> None:
        self.array: np.ndarray | None = array
        self.order = order
        self.shape = tuple(int(extent) for extent in array.shape)
        self.dtype = np.dtype(array.dtype)
        self.chunks = chunks
        self.attrs = _encoding("array", "0.2.0")
        self._scratch = scratch
        self._name = name
        self._stop = stop
        self._lock = threading.Lock()
        self._copied: np.ndarray | None = None
        self._done = 0
        self._began: float | None = None
        self.job: _Job | None = None
        scratch.keep(self)

    def copy(self) -> None:
        """The snapshot, by blocks of rows in priority order (the worker runs it)."""
        import time

        self._began = time.monotonic()
        array = self.array
        assert array is not None
        out = np.lib.format.open_memmap(
            self._scratch.path() / f"{self._name}-dense.npy", mode="w+", dtype=self.dtype, shape=self.shape
        )
        step = max(1, PASS_ENTRIES // max(1, int(np.prod(self.shape[1:], dtype=np.int64))))
        for start in range(0, self.shape[0], step):
            if self._stop.is_set():
                raise RuntimeError("The served store was closed while an array of it was copied.")
            end = min(self.shape[0], start + step)
            out[start:end] = array[start:end] if self.order is None else array[self.order[start:end]]
            self._done = end
        out.flush()
        with self._lock:
            self._copied = out

    def progress(self) -> tuple[float, float | None]:
        import time

        share = self._done / max(1, self.shape[0])
        if self._copied is not None:
            return 1.0, 0.0
        if self._began is None or self._done == 0:
            return share, None
        elapsed = time.monotonic() - self._began
        return share, elapsed * (self.shape[0] - self._done) / self._done

    def release(self) -> None:
        with self._lock:
            self._copied = None
            self.array = None

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        start = index[0] * self.chunks[0]
        end = min(self.shape[0], start + self.chunks[0])
        columns = (
            slice(index[1] * self.chunks[1], min(self.shape[1], (index[1] + 1) * self.chunks[1]))
            if len(self.shape) == 2
            else slice(None)
        )
        with self._lock:
            copied, array = self._copied, self.array
        if copied is not None:
            return np.asarray(copied[start:end, columns] if len(self.shape) == 2 else copied[start:end])
        if array is None:
            raise RuntimeError("The served store is closed.")
        rows = slice(start, end) if self.order is None else self.order[start:end]
        return np.asarray(array[rows, columns] if len(self.shape) == 2 else array[rows])


class _DenseGeneMajor(_Planned):
    """A dense expression matrix laid out once, gene-major in priority order,
    into a memory-mapped file, by the store's worker: an object's (its
    snapshot, in the type ``prepare()`` stores it), or a file's stored in
    chunks of whole rows (or of more genes than a served strip).

    Read through the rows, every served strip of a gene would decode every
    chunk of the matrix (300,000 x 200 in 1,000 x 200 gzip chunks: 1.2 s a
    gene, each time); laid out once -- one pass over the source by blocks of
    whole rows, memory one block -- a strip is one slice of the file it wrote.
    A strip asked for before then is taken from the source's rows at once
    (the snapshot, held read-only while it is laid out; a file's rows,
    decoded once into the scratch folder): the viewer never waits for it.
    """

    def __init__(
        self,
        array: Any,
        priority: CellPriority | None,
        dtype: np.dtype,
        chunks: tuple[int, int],
        scratch: Scratch,
        name: str,
        stop: threading.Event,
        worker: _Worker,
    ) -> None:
        self.array = array
        self.worker = worker
        self.new_row = None if priority is None or priority.identity else priority.new_row
        #: Each served row's source row (None: the source's own order).
        self.order = None if priority is None or priority.identity else priority.order
        self.source = Rows(array, scratch)
        self.shape = (int(array.shape[0]), int(array.shape[1]))
        self.dtype = np.dtype(dtype)
        self.chunks = chunks
        self.attrs = _encoding("array", "0.2.0")
        self._scratch = scratch
        self._name = name
        self._stop = stop
        self._lock = threading.Lock()
        self._laid: np.ndarray | None = None
        self._done = 0
        self._began: float | None = None
        self.job: _Job | None = None
        scratch.keep(self)

    def lay_out(self) -> None:
        """The pass, by blocks of whole rows (the worker runs it)."""
        import time

        self._began = time.monotonic()
        n_obs, n_vars = self.shape
        out = np.lib.format.open_memmap(
            self._scratch.path() / f"{self._name}-dense.npy", mode="w+", dtype=self.dtype, shape=(n_vars, n_obs)
        )
        step = max(1, PASS_ENTRIES // max(1, n_vars))
        for start in range(0, n_obs, step):
            if self._stop.is_set():
                raise RuntimeError("The served store was closed while its matrix was laid out.")
            block = np.asarray(self.array[start : start + step])
            rows = np.arange(start, start + len(block)) if self.new_row is None else self.new_row[start : start + len(block)]
            out[:, rows] = block.T
            self._done = start + len(block)
        out.flush()
        with self._lock:
            self._laid = out

    def progress(self) -> tuple[float, float | None]:
        import time

        share = self._done / max(1, self.shape[0])
        if self._laid is not None:
            return 1.0, 0.0
        if self._began is None or self._done == 0:
            return share, None
        elapsed = time.monotonic() - self._began
        return share, elapsed * (self.shape[0] - self._done) / self._done

    def release(self) -> None:
        with self._lock:
            self._laid = None
            self.array = None

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        assert self.job is not None
        start = index[0] * self.chunks[0]
        end = min(self.shape[0], start + self.chunks[0])
        genes = slice(index[1] * self.chunks[1], (index[1] + 1) * self.chunks[1])
        if not self.job.finished.is_set():
            # Not laid out yet: the strip's rows taken from the source as they will be laid out.
            with self._lock:
                if self.array is None:
                    raise RuntimeError("The served store is closed.")
            rows = np.arange(start, end, dtype=np.int64) if self.order is None else self.order[start:end]
            return np.ascontiguousarray(self.source.take(rows, genes), dtype=self.dtype)
        self.job.require()
        with self._lock:
            laid = self._laid
        if laid is None:
            raise RuntimeError("The served store is closed.")
        return np.ascontiguousarray(laid[genes, start:end].T)


def _reads_whole_rows(array: Any, strip: tuple[int, int]) -> bool:
    """Whether a stored dense matrix's chunks span more genes than a served strip (rows read whole for any gene)."""
    chunks = getattr(array, "chunks", None)
    if isinstance(array, np.ndarray) or not chunks or len(chunks) != 2:
        return False
    return int(chunks[1]) > strip[1]


class _LaidOut:
    """An expression matrix laid out gene-major by the store's worker
    (:class:`~cellucid.prepare_data._gene_major.GeneMajorLayout`, its output
    files in the scratch folder): its job, and how its arrays are served.

    What is asked for before the worker has placed it is made at once,
    exactly as it will be placed, from what the worker has placed so far and
    the source for the rest -- the snapshot of an object, the file
    (:meth:`GeneMajorLayout.entries_now`, :meth:`GeneMajorLayout.count_prefixes`);
    the viewer waits for a layout only when it ends sooner than its rest
    would be made again (:meth:`GeneMajorLayout.sooner`). Only a file's
    matrix not counted yet is "not ready" (where each gene's values go is
    known once counted): its keys are answered ``503`` with ``Retry-After``."""

    def __init__(self, layout: GeneMajorLayout, label: str, worker: _Worker, held: list[np.ndarray], protection: _ReadOnly) -> None:
        self.layout = layout
        self.worker = worker
        self.job = _Job(label, layout.work, layout.run, layout.progress, held=held, protection=protection)
        self._lock = threading.Lock()
        self._released = False
        #: The entries made last before they were placed (``(start, end)``,
        #: rows, values): a range's indices and data are asked for apart.
        self._made: tuple[tuple[int, int], np.ndarray, np.ndarray] | None = None
        #: The prefix counts made before the worker ended.
        self._counts: np.ndarray | None = None

    def release(self) -> None:
        """Let the mapped output files go (the store is closing)."""
        with self._lock:
            self._released = True
            self.layout.outputs.clear()
            self.layout.result = None

    def pointers(self) -> np.ndarray:
        if not self.layout.counted.is_set():
            raise self.worker.not_ready(self.job)
        self.job.require()
        assert self.layout.indptr is not None
        return self.layout.indptr

    def entries(self, part: str, start: int, end: int, chunk: tuple[int, int]) -> np.ndarray:
        """Entries ``[start, end)`` of ``part``, which lie in the served chunk
        ``chunk`` (``[first, last)``): placed, or -- not placed yet -- made with
        the rest of that chunk (kept for the other part and later ranges of
        it: a pass over a cell-major source makes every gene a chunk reaches)."""
        if not self.layout.ready(end) and self.layout.sooner():
            # The layout ends sooner than its rest would be made again.
            self.job.finished.wait()
        if not self.layout.ready(end):
            if self.job.finished.is_set():
                self.job.require()
            if not self.layout.counted.is_set():
                raise self.worker.not_ready(self.job)
            with self._lock:
                made = self._made
            if made is None or made[0] != chunk:
                made = (chunk, *self.layout.entries_now(*chunk))
                with self._lock:
                    self._made = made
            values = made[1] if part == "indices" else made[2]
            return values[start - chunk[0] : end - chunk[0]]
        with self._lock:
            if self._released:
                raise RuntimeError("The served store is closed.")
            return np.asarray(self.layout.outputs[part][start:end])

    def prefix_counts(self) -> np.ndarray:
        if not self.job.finished.is_set():
            if not self.layout.counted.is_set():
                raise self.worker.not_ready(self.job)
            with self._lock:
                counts = self._counts
            if counts is None:
                counts = self.layout.count_prefixes()
                with self._lock:
                    self._counts = counts
            return counts
        self.job.require()
        result = self.layout.result
        if result is None or result.prefix_counts is None:
            raise RuntimeError("The served store is closed.")
        return result.prefix_counts


class _LaidPointers(_Planned):
    def __init__(self, matrix: _LaidOut, nnz: int) -> None:
        self.matrix = matrix
        n_vars = matrix.layout.shape[1]
        self.shape = (n_vars + 1,)
        self.dtype = np.dtype(np.int64 if nnz >= 2**31 else np.int32)
        self.chunks = _vector_chunks(n_vars + 1, WHOLE_CHUNK_ELEMENTS)

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        return self.matrix.pointers()[index[0] * self.chunks[0] : (index[0] + 1) * self.chunks[0]]


class _LaidEntries(_Planned):
    def __init__(self, matrix: _LaidOut, part: str, nnz: int) -> None:
        self.matrix = matrix
        self.part = part
        layout = matrix.layout
        assert layout.dtype is not None
        self.shape = (nnz,)
        self.dtype = np.dtype(position_dtype(layout.shape[0]) if part == "indices" else layout.dtype)
        self.chunks = _vector_chunks(nnz, SERVED_ENTRY_CHUNK_ELEMENTS)

    def block(self, index: tuple[int, ...]) -> np.ndarray:
        start = index[0] * self.chunks[0]
        return self.block_part(index, 0, min(self.shape[0], start + self.chunks[0]) - start)

    def block_part(self, index: tuple[int, ...], first: int, last: int) -> np.ndarray:
        start = index[0] * self.chunks[0]
        chunk = (start, min(self.shape[0], start + self.chunks[0]))
        return self.matrix.entries(self.part, start + first, start + last, chunk)


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
        return self.block_part(index, 0, min(self.shape[0], start + self.chunks[0]) - start)

    def block_part(self, index: tuple[int, ...], first: int, last: int) -> np.ndarray:
        start = index[0] * self.chunks[0]
        indices, data = self.matrix.entries(start + first, start + last)
        return indices if self.part == "indices" else data


class _ServedSparse:
    """A source's compressed sparse matrix as served: its major axis in
    priority order (``major``: the source line of each served line) and/or its
    minor axis renumbered (``minor``: the served position of each source
    position, entries sorted by it within each line). Each chunk is made when
    it is asked for: the lazy form of a gene-major file's matrix, a graph, or
    an ``obsm`` matrix."""

    def __init__(
        self,
        source: SparseSource,
        *,
        major: np.ndarray | None,
        minor: np.ndarray | None,
        cache: _Cache,
    ) -> None:
        self.source = source
        self.form = source.form
        self.major = major
        self.minor = minor
        major_extent = source.n_lines
        #: The source's pointers, read once, when first needed.
        self._source = _Values((major_extent + 1,), np.int64, (major_extent + 1,), source.indptr)
        self.shape = source.shape
        self.nnz = source.nnz
        minor_extent = source.extent
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
        pointers = self._source_indptr()
        return pointers[lines], pointers[lines + 1]

    def _served_indptr(self) -> np.ndarray:
        if self.major is None:
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
        if self.minor is None:
            return indices, data
        positions = self.minor[indices]
        # Within each line by position, a repeated one in stored order.
        order = line_order(np.repeat(np.arange(len(lengths), dtype=np.int64), lengths), positions, len(self.minor))
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
            # Genes holding about PASS_ENTRIES values at a time.
            end = int(np.searchsorted(pointers, pointers[gene] + PASS_ENTRIES, side="right")) - 1
            end = min(n_vars, max(gene + 1, end))
            starts, ends = pointers[gene:end], pointers[gene + 1 : end + 1]
            rows = self.source.indices.ranges(starts, ends)
            if self.minor is not None:
                rows = self.minor[rows]
            counts[gene:end] = counts_below(rows, ends - starts, levels)
            gene = end
        return counts

# =============================================================================
# Elements of a source
# =============================================================================


def _served_obs(obs: pd.DataFrame, label: str) -> pd.DataFrame:
    """An in-memory object's ``obs`` with its text columns stored as ``prepare()``
    stores them (``_string_column``: a categorical when a value repeats or is
    missing, in code-point order, else a string array), so a served object
    and the same object prepared open alike: an object-dtype column (common
    before Scanpy first plots it) is not left to anndata's string array.
    """
    from .prepare_data._inputs import _string_column

    converted: dict[str, Any] = {}
    for column in obs.columns:
        series = obs[column]
        text = isinstance(series.dtype, pd.StringDtype) or (
            series.dtype == object and pd.api.types.infer_dtype(series, skipna=True) == "string"
        )
        if text:
            converted[column] = _string_column(series.to_numpy(dtype=object), f"{label} column {column!r}")
    if not converted:
        return obs
    served = obs.copy(deep=False)
    for column, value in converted.items():
        served[column] = value
    return served


def _as_element(value: Any, *, expression: bool = False) -> Any:
    """An element as the store reads it: numpy, scipy CSR/CSC, or an h5py/Zarr node.

    anndata's backed sparse datasets are read through their stored group; a
    dask array (anndata's lazily read element) is computed, a sparse one
    staying sparse -- except an expression matrix's (``expression``) dask
    array of sparse blocks, kept lazy and decoded block by block into the
    store's scratch folder as it is laid out
    (:func:`cellucid._source_arrays.sparse_source`), never whole in memory;
    a ``numpy.matrix`` is a plain array.
    """
    group = getattr(value, "group", None)
    if group is not None and hasattr(value, "format"):
        return group
    if expression and is_lazy_sparse(value):
        return value
    if hasattr(value, "compute") and hasattr(value, "dask") and not sparse.issparse(value):
        value = value.compute()
    if sparse.issparse(value):
        return value if value.format in ("csr", "csc") else sparse.csr_matrix(value)
    if isinstance(value, np.matrix):
        return np.asarray(value)
    return value


def _is_dense_array(value: Any) -> bool:
    """A numeric numpy, h5py or Zarr array (read by rows)."""
    return hasattr(value, "dtype") and hasattr(value, "shape") and not is_group(value) and np.dtype(value.dtype).kind in "biuf"


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


def _stored_parts(values: Any) -> list[tuple[str, np.ndarray]]:
    """The per-cell arrays anndata stores for a frame's index or column
    (``write_dataframe`` writes each one's ``._values``), each by the suffix
    of its path in the column: a categorical's codes (its categories are not
    per cell), a nullable column's values and mask (``write_nullable``: a
    missing value stored as 0 or ``""``), else the array itself. A type
    anndata cannot write has been refused by anndata already (the frame's
    head is written first)."""
    if isinstance(values, pd.Categorical):
        return [("/codes", np.asarray(values.codes))]
    if isinstance(values, pd.arrays.IntegerArray | pd.arrays.BooleanArray):
        return [("/values", values.to_numpy(na_value=0, dtype=values.dtype.numpy_dtype)), ("/mask", np.asarray(values.isna()))]
    if isinstance(values, pd.api.extensions.ExtensionArray) and isinstance(values.dtype, pd.StringDtype):
        return [("/values", values.to_numpy(na_value="")), ("/mask", np.asarray(values.isna()))]
    return [("", np.asarray(values))]


def _without(mapping: Mapping[str, Any], key: str) -> dict[str, Any]:
    return {name: value for name, value in mapping.items() if name != key}

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

    def __init__(self, cache_bytes: int | None = None) -> None:
        import zarr

        self._documents: dict[str, Any] = {}
        self._memory = zarr.storage.MemoryStore(store_dict=self._documents)
        self.root = zarr.open_group(self._memory, mode="w", zarr_format=2)
        self.root.attrs.update(_encoding("anndata", "0.1.0"))
        self._planned: dict[str, _Planned] = {}
        #: The genes (or rows) whose renumbered, sorted entries the store
        #: keeps, least recently used first out, within ``cache_bytes``
        #: (:func:`require_cache_bytes`): a long gene spans several chunks of
        #: the store, and is sorted once for all of them while it is kept.
        self._lines = _Cache(require_cache_bytes(cache_bytes))
        self._closers: list[Callable[[], None]] = []
        #: Set when the store closes: a layout still running stops at its next step.
        self._stop = threading.Event()
        #: The folder the layouts and decoded source arrays write; what maps a
        #: file of it is let go on close (nothing here refers back to the
        #: store, so dropping the store releases it).
        self._folder = Scratch(self._stop)
        #: The background worker laying the expression matrices out, the jobs
        #: it has (each matrix's, queued from the start or when first asked
        #: for), and the object's arrays it holds read-only while it reads them.
        self._worker = _Worker(self._stop)
        self._jobs: list[_Job] = []
        self._protection = _ReadOnly()
        self._finalizer = weakref.finalize(
            self,
            PriorityStore._release,
            self._closers,
            self._stop,
            self._worker,
            self._protection,
            self._folder,
            self._planned,
            self._lines,
        )
        self.etag = f'"served-{secrets.token_hex(12)}"'
        self._consolidated: bytes | None = None
        #: The seed of the store's priority order.
        self.seed = SERVED_PRIORITY_SEED
        self.n_obs = 0
        #: What the store holds: ``"anndata"``, or ``"mudata"`` (modalities under ``mod``).
        self.layout = "anndata"
        #: The rows of each AnnData group (``""``, ``"mod/<name>/"``): its arrays of one row per cell are chunked as such.
        self._groups: dict[str, int] = {}
        #: The cells' names in their original order (the order the viewer
        #: names cells in), taken with the rest of what is served; None when
        #: the source has no obs index.
        self.obs_names: pd.Index | None = None

    # -- building ------------------------------------------------------------

    @classmethod
    def of_anndata(cls, adata: anndata.AnnData, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of an AnnData object, its cells in priority order.

        Its cells' original positions are the ones ``prepare()`` gives them,
        by the one rule (:func:`~cellucid.prepare_data._priority.stored_order`):
        an object in a verified stored order (read from a file in priority
        order, its record kept) is taken from its stored rows by the original
        position each records -- the way back to the original order and on to
        the served one is one permutation, so an object in the served order
        already is served as it is; any other is taken in its own order, and
        a position column with no record that describes the cells names
        nothing (it is not served; a warning says so, as ``prepare()`` says it).
        """
        if adata.n_obs <= 0:
            raise ValueError("AnnData datasets must contain at least one observation")
        stored = stored_order_of(adata)
        if stored.notice is not None:
            warnings.warn(stored.notice, UserWarning, stacklevel=2)
        return cls._built(lambda store: store._plan_object(adata, stored.positions), cache_bytes)

    @classmethod
    def of_file(cls, source: AnnDataStore, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of an AnnData file (``.h5ad``, a Zarr folder or zip).

        A file in priority order already (``uns['cellucid']['priority']``,
        served this way because the browser cannot decode one of its codecs)
        keeps its own order; any other is put in priority order.
        """
        return cls._built(lambda store: store._plan_file(store._open(source)), cache_bytes)

    @classmethod
    def of_mudata_file(cls, source: AnnDataStore, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of a MuData file (``.h5mu``, a MuData Zarr folder or zip), its cells in priority order.

        The shared cells are in priority order, and each modality's rows in
        the order of its cells (``prepare()``'s MuData layout): a prefix of
        the cells is a prefix of every modality.
        """
        return cls._built(lambda store: store._plan_mudata_file(store._open(source)), cache_bytes)

    @classmethod
    def of_mudata(cls, mdata: Any, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of a MuData object, as :meth:`of_mudata_file` makes one of
        its file, from a snapshot as :meth:`of_anndata` makes one: no copy
        written, its cells by the original positions ``prepare()`` gives them."""
        stored = stored_order_of(mdata)
        if stored.notice is not None:
            warnings.warn(stored.notice, UserWarning, stacklevel=2)
        return cls._built(lambda store: store._plan_mudata_object(mdata, stored.positions), cache_bytes)

    @classmethod
    def of_loom(cls, source: AnnDataStore, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of a Loom file, presented as the browser presents it (:mod:`cellucid._served_formats`)."""
        from ._served_formats import loom_presented

        return cls._built(lambda store: store._plan_presented(loom_presented(store._open(source), source.path.name)), cache_bytes)

    @classmethod
    def of_cell_ranger(cls, source: AnnDataStore, *, cache_bytes: int | None = None) -> PriorityStore:
        """The store of a 10x Cell Ranger ``.h5``, or of an ``outs/`` folder with its analysis
        (:mod:`cellucid._served_formats`): one AnnData, or a MuData of its feature types."""
        from ._served_formats import cell_ranger_presented, outs_analysis, outs_matrix_file

        def plan(store: PriorityStore) -> None:
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

        return cls._built(plan, cache_bytes)

    @classmethod
    def _built(cls, plan: Callable[[PriorityStore], None], cache_bytes: int | None) -> PriorityStore:
        """A new store planned by ``plan``, or -- when planning fails -- the
        store closed and the error raised. The frames the error passed through
        are cleared first: their variables may hold an array mapping a file
        of the store (a snapshot being copied, a layout being written), which
        Windows would keep from being removed while the error is alive."""
        store = cls(cache_bytes)
        try:
            plan(store)
            store._settle()
        except BaseException as error:
            traceback.clear_frames(error.__traceback__)
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

    def _write_obs(self, path: str, obs: pd.DataFrame, priority: CellPriority | None) -> None:
        """A frame of one row per cell (``obs``), its rows in ``priority``'s
        order and, with an order, each row's original position
        (:data:`ORIGINAL_POSITION_COLUMN`) its last column (None: as stored).

        anndata writes the frame's empty head -- its encodings, column order,
        categories, each array's type -- and every per-cell array it would
        store (the index, a column, a categorical's codes, a nullable column's
        values and mask: :func:`_stored_parts`) is cut straight from the
        column, in served order, in the chunks ``prepare()`` cuts it in: no
        copy of the frame in served order, and no second pass to re-cut what
        anndata chunked otherwise (a 10M-cell ``obs``: seconds before the
        server's first answer)."""
        rows = len(obs)
        order = None if priority is None else priority.order
        head = obs.iloc[:0]
        if priority is not None:
            head = head.assign(**{ORIGINAL_POSITION_COLUMN: np.zeros(0, dtype=np.uint32)})
        self._write(path, head)
        index = json.loads(self._documents[f"{path}/.zattrs"].to_bytes())["_index"]
        parts: list[tuple[str, np.ndarray, np.ndarray | None]] = [
            (f"{path}/{index}{suffix}", values, order) for suffix, values in _stored_parts(obs.index._values)
        ]
        for column, series in obs.items():
            parts.extend((f"{path}/{column}{suffix}", values, order) for suffix, values in _stored_parts(series._values))
        if priority is not None:
            parts.append((f"{path}/{ORIGINAL_POSITION_COLUMN}", priority.original.astype(np.uint32), None))
        chunks = _cell_chunks(rows)
        for array_path, values, rows_order in parts:
            document = json.loads(self._documents[f"{array_path}/.zarray"].to_bytes())
            document = {**document, "shape": [rows], "chunks": list(chunks)}
            self._put(f"{array_path}/.zarray", document)
            size = chunks[0]
            for index_of_chunk in range(-(-rows // size)):
                taken = slice(index_of_chunk * size, (index_of_chunk + 1) * size)
                block = values[taken] if rows_order is None else values[rows_order[taken]]
                self._put_chunk(array_path, document, index_of_chunk, block)

    def _put_chunk(self, path: str, document: Mapping[str, Any], index: int, block: np.ndarray) -> None:
        """Chunk ``index`` (along the first axis) of an uncompressed array as
        Zarr v2 stores it: a text array's (``vlen-utf8``) encoded, any other's
        bytes; one at the array's edge stored whole, padded with zeros (``""``)."""
        import numcodecs
        from zarr.core.buffer import default_buffer_prototype

        chunks = document["chunks"]
        size = chunks[0]
        text = any(item.get("id") == "vlen-utf8" for item in document.get("filters") or ())
        dtype = np.dtype(object) if text else np.dtype(document["dtype"])
        block = np.asarray(block, dtype=dtype)
        if len(block) < size:
            padded = np.empty((size, *block.shape[1:]), dtype=dtype)
            padded[...] = "" if text else np.zeros((), dtype=dtype)
            padded[: len(block)] = block
            block = padded
        raw = numcodecs.VLenUTF8().encode(np.ascontiguousarray(block, dtype=object)) if text else np.ascontiguousarray(block).tobytes()
        name = document.get("dimension_separator", ".").join([str(index), *(["0"] * (len(chunks) - 1))])
        self._documents[f"{path}/{name}"] = default_buffer_prototype().buffer.from_bytes(bytes(raw))

    def _group(self, path: str, attrs: Mapping[str, Any] | None = None) -> None:
        self.root.create_group(path).attrs.update(dict(attrs) if attrs is not None else _encoding("dict", "0.1.0"))

    def _plan_object(self, adata: anndata.AnnData, given: np.ndarray | None) -> None:
        n_obs = int(adata.n_obs)
        # ``given``: each row's original position, for an object in a
        # verified stored order (its position column is that order); None:
        # each row its own, and a position column names nothing.
        obs = adata.obs
        if ORIGINAL_POSITION_COLUMN in obs.columns:
            without_positions({ORIGINAL_POSITION_COLUMN: None}, "obs", verified=given is not None)
            obs = obs.drop(columns=[ORIGINAL_POSITION_COLUMN])
        obs = _served_obs(obs, "obs")
        priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED, given)
        assert priority is not None
        levels = prefix_levels(n_obs)
        self.n_obs = n_obs
        self._plan_anndata(
            "",
            obs=obs,
            var=adata.var,
            uns=dict(adata.uns),
            x=None if adata.X is None else _as_element(adata.X, expression=True),
            layers={key: _as_element(value, expression=True) for key, value in adata.layers.items()},
            obsm=dict(adata.obsm),
            obsp=dict(adata.obsp),
            varm=dict(adata.varm),
            varp=dict(adata.varp),
            priority=priority,
            levels=levels,
            record=levels,
            snapshot=True,
        )
        raw = adata.raw
        if raw is not None:
            self._group("raw", _encoding("raw", "0.1.0"))
            self._write("raw/var", raw.var)
            self._write("raw/varm", dict(raw.varm))
            self._matrix("raw/X", _as_element(raw.X, expression=True), priority, levels=None, expression=True, snapshot=True)

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
        snapshot: bool = False,
    ) -> None:
        """One AnnData at ``prefix`` (``""``, or ``"mod/<name>/"``), its cells in ``priority``'s order.

        ``obs`` and the arrays are the source's, in its own order; ``levels``
        are the prefix levels of its rows (a modality's: its rows below each
        global level); ``record``, when given, the levels ``uns['cellucid']['priority']``
        records (the dataset's own, at the root). ``priority`` None keeps the
        rows as stored (a file in priority order already). ``snapshot``: the
        source is an object in memory (its matrices copied, every one laid out
        from the start).
        """
        n_obs = len(obs)
        self._groups[prefix] = n_obs
        if priority is not None:
            if prefix == "":
                self._name_cells(obs.index, priority)
            uns = _permute_quantiles(uns, priority.order, n_obs)
            if record is not None:
                uns = self._with_priority(uns, priority, record)
        self._write_obs(f"{prefix}obs", obs, priority)
        self._write(f"{prefix}var", var)
        self._write(f"{prefix}uns", uns)
        counts = self._matrix(f"{prefix}X", x, priority, levels=levels, expression=True, snapshot=snapshot, first=True)
        self._group(f"{prefix}layers")
        for key, value in layers.items():
            self._matrix(f"{prefix}layers/{key}", value, priority, levels=None, expression=True, snapshot=snapshot)
        self._per_cell(f"{prefix}obsm", obsm, priority)
        self._graphs(f"{prefix}obsp", obsp, priority)
        # Prefix counts the source had are for its own order.
        self._write(f"{prefix}varm", _without(varm, PREFIX_COUNTS_KEY) if priority is not None else dict(varm))
        if counts is not None:
            assert levels is not None
            self._plan_counts(f"{prefix}varm/{PREFIX_COUNTS_KEY}", counts, len(var), levels)
        self._write(f"{prefix}varp", dict(varp))

    def _plan_file(self, root: Any) -> None:
        if encoding_of(root) != "anndata":
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
            self._name_stored_cells(root, own)
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
        if "raw" in root and is_group(root["raw"]):
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
        self._name_cells(presented.obs.index, priority)
        self._mudata_root(names)
        self._groups[""] = n_obs
        self._write_obs("obs", presented.obs, priority)
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

        if (encoding_of(root) or "").lower() != "mudata" or "mod" not in root:
            raise ValueError("This MuData store has no MuData encoding or no mod group of modalities.")
        if int(np.asarray(root.attrs.get("axis", 0))) != 0:
            raise ValueError(
                "This MuData shares features, not cells (axis != 0); Cellucid shows cells. Write it with shared cells (axis=0)."
            )
        obs = read_elem(root["obs"])
        names = _modality_names(root["mod"])
        if not names:
            raise ValueError("The MuData store holds no modality.")
        if len(obs) <= 0:
            raise ValueError("MuData datasets must contain at least one observation")
        own = _verified_priority(root, warn=True)
        for name in names:
            if (encoding_of(root["mod"][name]) or "") != "anndata":
                raise ValueError(
                    f"Modality {name!r} predates anndata's current encoding (anndata 0.8): write the MuData again "
                    "with a current mudata to serve it."
                )
        stored_maps = root.get("obsmap")
        elements = {key: root.get(key, {}) for key in ("obsm", "obsp")}
        self._plan_mudata(
            obs=obs,
            uns=read_elem(root["uns"]) if "uns" in root else {},
            var=read_elem(root["var"]),
            obsm={key: elements["obsm"][key] for key in elements["obsm"]},
            obsp={key: elements["obsp"][key] for key in elements["obsp"]},
            varm=read_elem(root["varm"]) if "varm" in root else {},
            varp=read_elem(root["varp"]) if "varp" in root else {},
            varmap=read_elem(root["varmap"]) if "varmap" in root else {},
            obsmaps={
                name: np.asarray(stored_maps[name][...]) for name in names if stored_maps is not None and name in stored_maps
            },
            modalities={name: _rows_of(root["mod"][name]) for name in names},
            attrs={key: root.attrs[key] for key in root.attrs if key not in ("encoding-type", "encoding-version", "axis")},
            own=own,
            given=None,
            plan_modality=lambda name, prefix, priority, levels: self._plan_file_group(
                root["mod"][name], prefix, priority, levels, record=None
            ),
            modality_attrs=lambda name: {key: _attribute_value(root["mod"][name].attrs[key]) for key in root["mod"][name].attrs},
        )

    def _plan_mudata_object(self, mdata: Any, given: np.ndarray | None) -> None:
        """A MuData held in memory, planned as its file would be (served from a
        snapshot, as an AnnData object is: :meth:`_plan_object`)."""
        if int(getattr(mdata, "axis", 0)) != 0:
            raise ValueError(
                "This MuData shares features, not cells (axis != 0); Cellucid shows cells. Make one with shared cells (axis=0)."
            )
        names = [str(name) for name in mdata.mod]
        if not names:
            raise ValueError("The MuData holds no modality.")
        if int(mdata.n_obs) <= 0:
            raise ValueError("MuData datasets must contain at least one observation")
        obs = mdata.obs
        if ORIGINAL_POSITION_COLUMN in obs.columns:
            without_positions({ORIGINAL_POSITION_COLUMN: None}, "obs", verified=given is not None)
            obs = obs.drop(columns=[ORIGINAL_POSITION_COLUMN])
        obs = _served_obs(obs, "obs")
        # A modality of a MuData in a verified stored order holds its rows in
        # that order too, each one's original row in its own position column.
        given_rows: dict[str, np.ndarray | None] = {}
        for name in names:
            modality = mdata.mod[name]
            column = ORIGINAL_POSITION_COLUMN in modality.obs.columns
            if given is not None and column:
                from .prepare_data._priority import _require_positions

                given_rows[name] = _require_positions(modality.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), int(modality.n_obs))
            else:
                if column:
                    without_positions({ORIGINAL_POSITION_COLUMN: None}, f"mdata.mod[{name!r}].obs")
                given_rows[name] = None

        def plan_modality(name: str, prefix: str, priority: CellPriority | None, levels: np.ndarray | None) -> None:
            modality = mdata.mod[name]
            self._plan_anndata(
                prefix,
                obs=_served_obs(modality.obs.drop(columns=[ORIGINAL_POSITION_COLUMN], errors="ignore"), f"mdata.mod[{name!r}].obs"),
                var=modality.var,
                uns=dict(modality.uns),
                x=None if modality.X is None else _as_element(modality.X, expression=True),
                layers={key: _as_element(value, expression=True) for key, value in modality.layers.items()},
                obsm=dict(modality.obsm),
                obsp=dict(modality.obsp),
                varm=dict(modality.varm),
                varp=dict(modality.varp),
                priority=priority,
                levels=levels,
                record=None,
                snapshot=True,
            )
            raw = modality.raw
            if raw is not None:
                self._group(f"{prefix}raw", _encoding("raw", "0.1.0"))
                self._write(f"{prefix}raw/var", raw.var)
                self._write(f"{prefix}raw/varm", dict(raw.varm))
                self._matrix(f"{prefix}raw/X", _as_element(raw.X, expression=True), priority, levels=None, expression=True, snapshot=True)

        self._plan_mudata(
            obs=obs,
            uns=dict(mdata.uns),
            var=mdata.var,
            obsm=dict(mdata.obsm),
            obsp=dict(mdata.obsp),
            varm=dict(mdata.varm),
            varp=dict(mdata.varp),
            varmap={name: np.asarray(mdata.varmap[name]) for name in names},
            obsmaps={name: np.asarray(mdata.obsmap[name]) for name in names},
            modalities={name: int(mdata.mod[name].n_obs) for name in names},
            attrs={},
            own=None,
            given=given,
            given_rows=given_rows,
            plan_modality=plan_modality,
            modality_attrs=lambda name: _encoding("anndata", "0.1.0"),
            snapshot=True,
        )

    def _plan_mudata(
        self,
        *,
        obs: pd.DataFrame,
        uns: Mapping[str, Any],
        var: Any,
        obsm: Mapping[str, Any],
        obsp: Mapping[str, Any],
        varm: Any,
        varp: Any,
        varmap: Any,
        obsmaps: Mapping[str, np.ndarray],
        modalities: Mapping[str, int],
        attrs: Mapping[str, Any],
        own: int | None,
        given: np.ndarray | None,
        plan_modality: Callable[[str, str, CellPriority | None, np.ndarray | None], None],
        modality_attrs: Callable[[str], Mapping[str, Any]],
        given_rows: Mapping[str, np.ndarray | None] | None = None,
        snapshot: bool = False,
    ) -> None:
        """One MuData, a file's or an object's: the shared cells in priority
        order (or ``own``: a file's verified order, kept), each modality's rows
        in the order of its cells (``prepare()``'s MuData layout), so a prefix
        of the cells is a prefix of every modality. ``modalities`` names each
        modality's rows; ``given``/``given_rows``, for an object in a verified
        stored order, each cell's and each modality row's original position."""
        n_obs = len(obs)
        names = list(modalities)
        self.n_obs = n_obs
        priority: CellPriority | None = None
        levels: np.ndarray | None = None
        if own is not None:
            # Served through the store because the browser cannot decode a codec: kept in its own order.
            self.seed = own
            self.obs_names = _in_original_order(obs.index, own)
        else:
            # Each row its own original position (a position column replaced),
            # or the one an object in a verified stored order records.
            priority = CellPriority.of(n_obs, SERVED_PRIORITY_SEED, given)
            assert priority is not None
            levels = prefix_levels(n_obs)
        self._mudata_root(names)
        if attrs:
            self.root.attrs.update({key: _attribute_value(value) for key, value in attrs.items()})
        # A modality's cell map; one with none measured every cell, in order
        # (as the browser reads it: ``views/modalities.js``).
        maps: dict[str, np.ndarray] = {}
        for name in names:
            if name in obsmaps:
                maps[name] = np.asarray(obsmaps[name])
                continue
            rows = modalities[name]
            if rows != n_obs:
                raise ValueError(
                    f"Modality {name!r} has {rows:,} {plural(rows, 'cell')} and the dataset {n_obs:,}, and no cell map "
                    "(obsmap) says which they are."
                )
            maps[name] = np.arange(1, n_obs + 1, dtype=np.uint32)
        if priority is not None:
            self._name_cells(obs.index, priority)
            uns = self._with_priority(_permute_quantiles(uns, priority.order, n_obs), priority, levels)  # type: ignore[arg-type]
        self._groups[""] = n_obs
        self._write_obs("obs", obs, priority)
        self._write("uns", uns)
        self._write("var", var)
        self._per_cell("obsm", obsm, priority)
        self._graphs("obsp", obsp, priority)
        self._write("varm", varm)
        self._write("varp", varp)
        self._write("varmap", varmap)
        served_maps: dict[str, np.ndarray] = {}
        for name in names:
            self._group(f"mod/{name}", modality_attrs(name))
            if priority is None:
                served_maps[name] = maps[name]
                plan_modality(name, f"mod/{name}/", None, None)
                continue
            given_of_rows = None if given_rows is None else given_rows.get(name)
            modality, served_map, below = _modality_order(priority, maps[name], modalities[name], given_of_rows)
            served_maps[name] = served_map
            assert levels is not None
            plan_modality(name, f"mod/{name}/", modality, below[levels])
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

    def _name_cells(self, index: Any, priority: CellPriority) -> None:
        """The cells' names in their original order (:attr:`obs_names`): the
        source's ``index`` (its own rows), each at the original position
        ``priority`` gives its row."""
        names = np.asarray(pd.Index(index), dtype=object)
        named = np.empty(len(names), dtype=object)
        named[priority.original] = names[priority.order]
        self.obs_names = pd.Index(named, dtype=object)

    def _name_stored_cells(self, root: Any, seed: int) -> None:
        """:attr:`obs_names` of a file served in its own verified priority order (``seed``)."""
        from anndata.io import read_elem

        obs = root["obs"]
        key = obs.attrs.get("_index", "_index")
        index = read_elem(obs[key.decode() if isinstance(key, bytes) else str(key)])
        self.obs_names = _in_original_order(pd.Index(index), seed)

    def _laid_out_job(self, matrix: _LaidOut | _DenseGeneMajor | _DenseSnapshot, *, first: bool) -> None:
        """Hand a matrix's layout to the worker: now (``first``: an object's,
        or ``X``), or when one of its keys is first asked for."""
        job = matrix.job
        assert job is not None
        self._jobs.append(job)
        if first:
            self._worker.queue(job)

    def _matrix(
        self,
        path: str,
        value: Any,
        priority: CellPriority | None,
        *,
        levels: np.ndarray | None,
        expression: bool,
        snapshot: bool = False,
        first: bool = False,
    ) -> Callable[[], np.ndarray] | None:
        """Plan one per-cell matrix at ``path``.

        ``expression``: a cells x genes matrix (``X``, a layer, ``raw/X``):
        laid out gene-major by the worker when it is cell-major or held in
        memory (``snapshot``: an object's, laid out from the start, its
        arrays read-only meanwhile; ``first``: ``X``, laid out from the
        start; any other when first asked for). Returns how its prefix counts
        are made when ``levels`` are given and it is served gene-major and sparse.
        """
        if value is None:
            return None
        order = None if priority is None or priority.identity else priority.order
        from ._served_formats import CellRangerX

        if isinstance(value, CellRangerX):
            # Cell-major, as stored: a feature type's columns laid out gene-major.
            matrix = value.matrix
            group = matrix.group
            whole = value.start == 0 and value.count == matrix.n_features
            source = SparseSource(
                form="csr",
                shape=(matrix.n_cells, matrix.n_features),
                indptr=lambda: np.asarray(group["indptr"][...]),
                indices=Vector(group["indices"], self._folder),
                data=Vector(group["data"], self._folder),
            )
            genes = None if whole else np.arange(value.start, value.start + value.count, dtype=np.int64)
            return self._laid_out(path, source, priority, levels, genes=genes, snapshot=False, first=first)
        if expression and (sparse.issparse(value) or is_sparse_group(value) or is_lazy_sparse(value)):
            source = sparse_source(value, self._folder)
            if source.form == "csr" or sparse.issparse(value):
                return self._laid_out(path, source, priority, levels, genes=None, snapshot=snapshot, first=first or snapshot)
        if sparse.issparse(value) and priority is not None:
            # A graph or an obsm matrix in memory: a snapshot, canonical.
            value = sparse.csr_matrix(value, copy=True) if value.format == "csr" else sparse.csc_matrix(value, copy=True)
            value.sum_duplicates()
            value.eliminate_zeros()
        if sparse.issparse(value) or is_sparse_group(value):
            source = sparse_source(value, self._folder)
            renumber = priority is not None and not priority.identity and source.form == "csc"
            served = _ServedSparse(
                source,
                major=order if source.form == "csr" else None,
                minor=priority.new_row if renumber and priority is not None else None,
                cache=self._lines,
            )
            self._plan_sparse(path, served, SERVED_ENTRY_CHUNK_ELEMENTS if expression else WHOLE_CHUNK_ELEMENTS)
            if levels is not None and source.form == "csc":
                # Counted by the worker, the first time a subsample asks for them.
                counted: dict[str, np.ndarray] = {}

                def count() -> None:
                    counted["counts"] = served.prefix_counts(levels)

                nnz = served.nnz
                job = _Job(f"{path} (its prefix counts)", nnz, count, lambda: (0.0, None))
                self._jobs.append(job)

                def make() -> np.ndarray:
                    if not job.finished.is_set():
                        raise self._worker.not_ready(job)
                    job.require()
                    return counted["counts"]

                return make
            return None
        if not _is_dense_array(value):
            raise ValueError(f"{path} is not a matrix Cellucid serves ({type(value).__name__}).")
        shape = tuple(int(extent) for extent in value.shape)
        if len(shape) == 2 and expression and (isinstance(value, np.ndarray) or _reads_whole_rows(value, _strip_chunks((shape[0], shape[1])))):
            # An object's dense matrix (its snapshot, in the type prepare()
            # stores it), or a file's stored in whole rows: laid out gene-major.
            dtype = np.dtype(value.dtype)
            held: list[np.ndarray] = []
            if isinstance(value, np.ndarray):
                from .prepare_data._inputs import expression_dtype, require_expression_values

                require_expression_values(value, path)
                dtype = np.dtype(expression_dtype(value))
                held = self._protection.hold([value])
            laid = _DenseGeneMajor(
                value, priority, dtype, _strip_chunks((shape[0], shape[1])), self._folder, path.replace("/", "-"), self._stop, self._worker
            )
            laid.job = _Job(path, shape[0] * shape[1], laid.lay_out, laid.progress, held=held, protection=self._protection)
            self._laid_out_job(laid, first=first or isinstance(value, np.ndarray))
            self._plan(path, laid)
            return None
        chunks: tuple[int, ...]
        if len(shape) == 2:
            chunks = _strip_chunks((shape[0], shape[1])) if expression else _row_chunks((shape[0], shape[1]))
        else:
            # One value per cell.
            chunks = _cell_chunks(shape[0])
        if isinstance(value, np.ndarray) and priority is not None:
            # In memory: a snapshot in priority order (rows as served), copied
            # by the worker, the array read-only meanwhile.
            copied = _DenseSnapshot(value, order, chunks, self._folder, path.replace("/", "-"), self._stop)
            values = int(np.prod(shape, dtype=np.int64))
            copied.job = _Job(f"{path} (its snapshot)", values, copied.copy, copied.progress, held=self._protection.hold([value]), protection=self._protection)
            self._laid_out_job(copied, first=True)
            self._plan(path, copied)
            return None
        self._plan(path, _DenseRows(value, order, chunks, self._folder))
        return None

    def _laid_out(
        self,
        path: str,
        source: SparseSource,
        priority: CellPriority | None,
        levels: np.ndarray | None,
        *,
        genes: np.ndarray | None,
        snapshot: bool,
        first: bool,
    ) -> Callable[[], np.ndarray] | None:
        """An expression matrix laid out gene-major by the worker
        (:class:`~cellucid.prepare_data._gene_major.GeneMajorLayout`) into
        memory-mapped files. An object's (``snapshot``) is stored as
        ``prepare()`` stores it: counted now -- its size and type are in the
        store's metadata, and a value ``prepare()`` refuses is refused at once
        -- and placed by the worker, its arrays read-only meanwhile. A file's
        keeps its values as stored, its size and type known before any value
        is read; counted and placed by the worker."""
        name = path.replace("/", "-")
        folder = self._folder

        def allocate(part: str, dtype: np.dtype, length: int) -> np.ndarray:
            if length == 0:
                return np.zeros(0, dtype=dtype)
            return np.lib.format.open_memmap(folder.path() / f"{name}-{part}.npy", mode="w+", dtype=dtype, shape=(length,))

        order, new_row = (None, None) if priority is None else priority.rows_moved
        layout = GeneMajorLayout(
            source,
            order=order,
            new_row=new_row,
            label=path,
            genes=genes,
            levels=levels,
            allocate=allocate,
            copy=True,
            as_stored=not snapshot,
            stop=self._stop,
        )
        held = self._protection.hold([source.matrix.indptr, source.matrix.indices, source.matrix.data]) if snapshot and source.matrix is not None else []
        try:
            if layout.nnz is None:
                layout.count()
        except BaseException:
            self._protection.release(held)
            raise
        nnz = layout.nnz
        assert nnz is not None
        laid = _LaidOut(layout, path, self._worker, held, self._protection)
        self._folder.keep(laid)
        self._laid_out_job(laid, first=first)
        n_obs, n_vars = layout.shape
        self._sparse_group(path, "csc", (n_obs, n_vars))
        self._plan(f"{path}/indptr", _LaidPointers(laid, nnz))
        self._plan(f"{path}/indices", _LaidEntries(laid, "indices", nnz))
        self._plan(f"{path}/data", _LaidEntries(laid, "data", nnz))
        return None if levels is None else laid.prefix_counts

    def _plan_counts(self, path: str, make: Callable[[], np.ndarray], n_vars: int, levels: np.ndarray) -> None:
        """``varm['cellucid_prefix_counts']``: one column per level, each read whole (made on first read)."""
        counts = _Values((n_vars, len(levels)), np.uint32, (max(1, min(n_vars, WHOLE_CHUNK_ELEMENTS)), 1), make)
        counts.attrs = _encoding("array", "0.2.0")
        self._plan(path, counts)

    def _per_cell(self, key: str, elements: Mapping[str, Any], priority: CellPriority | None) -> None:
        """``obsm``: dense arrays and sparse matrices planned, anything else (a table) written in priority order."""
        from anndata.io import read_elem

        self._group(key)
        order = None if priority is None or priority.identity else priority.order
        for name, value in elements.items():
            value = _as_element(value)
            path = f"{key}/{name}"
            if sparse.issparse(value) or is_sparse_group(value) or _is_dense_array(value) and len(value.shape) == 2:
                self._matrix(path, value, priority, levels=None, expression=False)
                continue
            if is_group(value) or hasattr(value, "attrs") and hasattr(value, "dtype") and not isinstance(value, np.ndarray):
                value = read_elem(value)
            self._write(path, value if order is None else _take_rows(value, order))

    def _graphs(self, key: str, elements: Mapping[str, Any], priority: CellPriority | None) -> None:
        """``obsp``: each graph with both axes in priority order."""
        from anndata.io import read_elem

        self._group(key)
        for name, value in elements.items():
            value = _as_element(value)
            path = f"{key}/{name}"
            if sparse.issparse(value) or is_sparse_group(value):
                moved = priority is not None and not priority.identity
                served = _ServedSparse(
                    sparse_source(value, self._folder),
                    major=priority.order if moved and priority is not None else None,
                    minor=priority.new_row if moved and priority is not None else None,
                    cache=self._lines,
                )
                self._plan_sparse(path, served, WHOLE_CHUNK_ELEMENTS)
                continue
            dense = np.asarray(value if isinstance(value, np.ndarray) else read_elem(value))
            self._write(path, dense if priority is None or priority.identity else dense[np.ix_(priority.order, priority.order)])

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
        import zarr

        values = np.asarray(zarr.open_array(self._memory, path=path, mode="r", zarr_format=2)[...])
        for name in [name for name in self._documents if name.startswith(f"{path}/") and "/" not in name[len(path) + 1 :] and not name[len(path) + 1 :].startswith(".")]:
            del self._documents[name]
        document = {**document, "chunks": list(chunks)}
        self._put(f"{path}/.zarray", document)
        size = chunks[0]
        for index in range(-(-values.shape[0] // size)):
            self._put_chunk(path, document, index, values[index * size : (index + 1) * size])

    # -- serving -------------------------------------------------------------

    def _consolidated_metadata(self) -> bytes:
        if self._consolidated is None:
            metadata = {}
            for key in sorted(self._documents):
                if key.rsplit("/", 1)[-1] in (".zgroup", ".zattrs", ".zarray"):
                    metadata[key] = json.loads(self._documents[key].to_bytes())
            self._consolidated = json.dumps({"zarr_consolidated_format": 1, "metadata": metadata}).encode()
        return self._consolidated

    def body(self, key: str) -> StoreBody | None:
        """One key of the store as its length and its bytes made on request,
        or None when it has none (refused once it is closed): a server answers
        a ``Range`` of a chunk with only those bytes made. A chunk of a matrix
        not laid out yet is made at once from the source; only a chunk of a
        file's matrix the worker has not counted yet raises :class:`NotReady`,
        which says how long it may take."""
        if self._stop.is_set():
            raise RuntimeError("The served store is closed.")
        if key == ".zmetadata":
            return StoreBody.of(self._consolidated_metadata())
        found = self._documents.get(key)
        if found is not None:
            return StoreBody.of(bytes(found.to_bytes()))
        path, _, name = key.rpartition("/")
        array = self._planned.get(path)
        if array is None or not name:
            return None
        try:
            index = tuple(int(part) for part in name.split("."))
        except ValueError:
            return None
        if not array.in_grid(index):
            return None
        return StoreBody(array.chunk_length(), lambda start, end: array.chunk_bytes(index, start, end))

    def get(self, key: str) -> bytes | None:
        """One key of the store whole (:meth:`body`), or None when it has none."""
        found = self.body(key)
        return None if found is None else found.read(0, found.size)

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

    def laid_out(self, timeout: float | None = None) -> None:
        """Lay out every matrix the worker has (those laid out when first asked
        for, too) and wait until each has ended; raise the first error."""
        self._worker.wait(list(self._jobs), timeout)

    @staticmethod
    def _release(
        closers: list[Callable[[], None]],
        stop: threading.Event,
        worker: _Worker,
        protection: _ReadOnly,
        scratch: Scratch,
        planned: dict[str, _Planned],
        lines: _Cache,
    ) -> None:
        """Close the store: the worker stopped (a layout running ends at its
        next block), the object's arrays writeable again, every array that
        maps a file let go -- the layouts, and the planned arrays and cached
        lines that hold the rest (a snapshot, a source array) -- its scratch
        folder removed, then its files closed. Run on :meth:`close`, or when
        the store is dropped (its attributes are still alive then, so they are
        emptied here, not left to be collected after)."""
        stop.set()
        worker.close()
        protection.release_all()
        planned.clear()
        lines.clear()
        scratch.close()
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


def _modality_order(
    priority: CellPriority, obsmap: np.ndarray, n_rows: int, given: np.ndarray | None = None
) -> tuple[CellPriority, np.ndarray, np.ndarray]:
    """A MuData modality's rows in the order of its cells (``prepare()``'s
    ``CellPriority.modality``): its priority, its served cell map, and its
    rows below each served cell (``below[k]``: the modality rows of the first
    ``k`` cells). Each of its rows is its own original position, or the one
    ``given`` names (a modality of an object in a verified stored order)."""
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
    modality = CellPriority(
        seed=priority.seed,
        order=mod_order,
        original=mod_order if given is None else np.asarray(given, dtype=np.int64)[mod_order],
        new_row=new_row,
        identity=bool(np.array_equal(mod_order, np.arange(n_rows))),
    )
    return modality, served, below


def _in_original_order(stored: pd.Index, seed: int) -> pd.Index:
    """Names of rows in a verified priority order (``seed``), put in their original order."""
    priority = CellPriority.of(len(stored), seed)
    assert priority is not None
    names = np.empty(len(stored), dtype=object)
    names[priority.original] = np.asarray(stored, dtype=object)
    return pd.Index(names, dtype=object)


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
    cellucid = uns.get("cellucid") if uns is not None and is_group(uns) else None
    group = cellucid.get("priority") if cellucid is not None and is_group(cellucid) else None
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
        encoding = (encoding_of(root) or "").lower()
        if encoding not in ("anndata", "mudata"):
            return "older"
        return "priority" if _verified_priority(root, warn=False) is not None else "as-given"
