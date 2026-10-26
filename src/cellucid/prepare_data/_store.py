"""
Writing a Cellucid-ready AnnData: the encodings, the chunks and the codecs.

Every element is written with AnnData's own encodings (the anndata 0.8+
element specification: ``encoding-type`` and ``encoding-version`` on every
group and array), so ``anndata.read_h5ad`` and ``anndata.read_zarr`` read the
result without a warning. What this module decides is only *how* each array is
laid out, because the web app reads the file itself with HTTP range requests
and every chunk is one request and one decode:

- arrays the viewer reads by slices -- the ``data`` and ``indices`` of the
  expression matrix, read one gene's range at a time -- are cut into chunks of
  :data:`SLICED_CHUNK_ELEMENTS`, so one gene costs little beyond its own
  values;
- arrays it reads whole -- obs and var columns, the indexes, embeddings,
  vector fields and the neighbour graph -- into chunks of
  :data:`WHOLE_CHUNK_ELEMENTS`, so a whole read is a handful of requests;
- a two-dimensional array read by rows (an embedding) holds every column of a
  block of rows in one chunk; a dense expression matrix holds a tall strip of
  few genes, since it is read one gene at a time.

Every chunk is compressed with DEFLATE (:data:`DEFLATE_LEVEL`), which browsers
inflate natively (``DecompressionStream``): HDF5's gzip filter in an
``.h5ad``, zlib (Zarr v2) or gzip (Zarr v3) in a Zarr store, unsharded. Before
it, every array of values wider than one byte is byte-shuffled, which puts the
similar high bytes of neighbouring values together -- except the expression
values themselves (a sparse ``X``'s ``data``, a dense ``X``), which repeat
whole and compress best as they are. The shuffle is HDF5's shuffle filter in
an ``.h5ad`` and numcodecs' ``Shuffle`` filter in Zarr v2 (a v2 filter is any
numcodecs codec). Zarr v3 has no shuffle codec of its own -- numcodecs'
(``numcodecs.shuffle``) is outside the v3 specification -- so a shuffled v3
array is compressed with the specification's ``blosc`` codec instead: byte
shuffle and zlib at the same level (:func:`_blosc_codec`); every other v3
array is gzip.

Measured on the 561,947-cell, 8,192-gene Suo atlas with the web app's own
reader (``scripts/benchmark-anndata.mjs`` in the web repository; the numbers
are in the specification page): shuffling cut the Zarr store from 243.4 MB to
198 MB -- the gene indices by a third, embeddings and graph weights by 7-9% --
while shuffled expression values of real data grew 14-45% (raw counts,
``normalize_total`` + ``log1p``, proportions); Zstandard (with or without
Blosc) and Blosc-LZ4 were all larger than shuffled deflate and are decoded in
JavaScript rather than natively; chunks of 2^15 sliced elements doubled the
requests of a whole-matrix read, 2^17 and 2^18 fetched 1.7 and 3 times the
bytes of a batch of genes; sharding Zarr v3 cut the requests of a whole-matrix read
but cost every first gene of a shard one more round trip for its index.

The whole output is written beside its target under a temporary name and moved
into place only once it is complete, so an interrupted or failed call never
leaves a half-written file at the path a viewer would open.
"""

from __future__ import annotations

import os
import shutil
import warnings
import zipfile
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Protocol

import numpy as np
import pandas as pd
from scipy import sparse

#: Elements per chunk of an array the viewer reads by slices: the data and
#: indices of the expression matrix, read one gene's range at a time.
SLICED_CHUNK_ELEMENTS = 2**16

#: Elements per chunk of an array the viewer reads whole.
WHOLE_CHUNK_ELEMENTS = 2**18

#: The deflate level of every compressed chunk, zlib's own default. On the
#: dataset above, level 9 made the expression values 6% and everything else
#: under 1% smaller and wrote them twenty times slower; level 4 made the
#: expression values 6% larger. The browser inflates every level equally fast.
DEFLATE_LEVEL = 6

#: The block HDF5 gathers an ``.h5ad``'s metadata in (:func:`_write_h5ad`):
#: the web app's metadata block on a web host (256 KiB,
#: ``REMOTE_POLICY.blockBytes``). The dataset above has 44 KB of metadata;
#: the rest of the block holds the chunk indexes written with the values. More
#: metadata goes on in further blocks; the unused part of the block stays in
#: the file (the dataset above grew by 58 KB).
H5AD_META_BLOCK_BYTES = 256 * 1024

#: The Zarr format ``prepare()`` writes unless asked for the other one, as
#: cellucid-r's ``cellucid_prepare()`` does: every anndata release since 0.8
#: reads it, where Zarr v3 needs anndata 0.11; the web app reads both equally
#: fast.
DEFAULT_ZARR_FORMAT: Literal[2, 3] = 2

#: The timestamp of every ``.zarr.zip`` entry, the earliest a ZIP can carry,
#: so an archive's bytes never depend on when it was written.
_ZIP_EPOCH = (1980, 1, 1, 0, 0, 0)

#: The buffer one ``.zarr.zip`` entry is copied through.
_COPY_BYTES = 1 << 20

#: The metadata documents of a Zarr store, which no element key may shadow.
ZARR_METADATA_NAMES = frozenset({".zarray", ".zattrs", ".zgroup", ".zmetadata", "zarr.json"})


@dataclass(frozen=True)
class OutputTarget:
    """Where the output goes and what it is."""

    path: Path
    #: ``h5ad`` and ``h5mu`` are HDF5 files (an ``.h5mu`` with MuData's
    #: 512-byte user block before the HDF5 signature), ``zarr`` a Zarr store.
    kind: Literal["h5ad", "h5mu", "zarr"]
    #: A Zarr store packed into one stored (uncompressed) ``.zarr.zip``.
    zipped: bool
    zarr_format: Literal[2, 3] | None


@dataclass(frozen=True)
class Elements:
    """Everything one Cellucid-ready AnnData holds, checked and converted.

    ``obs`` and ``var`` map each column name to what is stored for it: a
    :class:`pandas.Categorical`, an object array of strings, a numeric or
    boolean array, or a pandas nullable integer or boolean array.
    """

    x: sparse.csc_matrix | np.ndarray | None
    obs_index: np.ndarray
    obs: dict[str, Any]
    var_index: np.ndarray
    var: dict[str, Any]
    obsm: dict[str, np.ndarray]
    obsp: dict[str, sparse.csr_matrix]
    uns: dict[str, Any]


class StoreContent(Protocol):
    """What a store holds when it is not one AnnData: a MuData (``_mudata``).

    :meth:`write` lays it out under ``root`` with :class:`_Writer` -- an HDF5
    file's root group (``kind`` ``"h5ad"``) or a Zarr group.
    """

    def write(self, root: Any, kind: Literal["h5ad", "zarr"], zarr_format: int | None) -> None: ...


# =============================================================================
# One array, one group
# =============================================================================


def _encoding(kind: str, version: str) -> dict[str, str]:
    return {"encoding-type": kind, "encoding-version": version}


def _reproducible_gzip(level: int) -> Any:
    """Zarr format 3's ``gzip`` codec, writing the same bytes on every run.

    numcodecs' ``GZip`` stamps the time of writing into each chunk's gzip
    header, so two runs over the same input would differ in every chunk.
    This one writes a zero time (as ``gzip -n`` does); the metadata still
    names the standard ``gzip`` codec, and readers see no difference.
    """
    import gzip
    from functools import cached_property

    from numcodecs.compat import ensure_contiguous_ndarray
    from numcodecs.gzip import GZip
    from zarr.codecs import GzipCodec

    class _GZip(GZip):  # type: ignore[misc]
        def encode(self, buf: Any) -> bytes:
            return gzip.compress(ensure_contiguous_ndarray(buf), compresslevel=self.level, mtime=0)

    class _GzipCodec(GzipCodec):
        @cached_property
        def _gzip_codec(self) -> GZip:
            return _GZip(self.level)

    return _GzipCodec(level=level)


def _blosc_codec(itemsize: int) -> Any:
    """Zarr format 3's ``blosc`` codec for a shuffled array of ``itemsize``-byte values.

    Blosc byte-shuffles each block of a chunk and deflates each byte plane of
    it as a zlib stream (``cname`` zlib, ``clevel`` :data:`DEFLATE_LEVEL`),
    which browsers inflate natively. The block size is Blosc's own choice
    (``blocksize`` 0): the whole chunk, or blocks of 1 MiB in a larger one
    (a float64 chunk of 2^18 values is two). On the He atlas
    (71,650 cells) a Zarr v3 store written so came to 29.00 MB, against
    29.11 MB with ``numcodecs.shuffle`` and gzip, and 37.62 MB with gzip and
    no shuffle (the gene indices of ``X`` alone 64% larger).

    zarr-python compresses with numcodecs' Blosc, which it sets to compress in
    the calling thread (``numcodecs.blosc.use_threads = False``), so a chunk
    of several blocks holds them in order and every run writes the same bytes.
    cellucid-r's ``cellucid_prepare()`` writes the same frames: the same
    header, blocks and streams, each stream deflated by the zlib R links.
    """
    from zarr.codecs import BloscCodec

    return BloscCodec(typesize=itemsize, cname="zlib", clevel=DEFLATE_LEVEL, shuffle="shuffle", blocksize=0)


def _vector_chunks(length: int, target: int) -> tuple[int]:
    return (max(1, min(length, target)),)


def _row_chunks(shape: tuple[int, int]) -> tuple[int, int]:
    """Every column of a block of rows: an array read whole, by rows."""
    rows, columns = shape
    columns = max(1, columns)
    return (max(1, min(rows, WHOLE_CHUNK_ELEMENTS // columns)), columns)


def _strip_chunks(shape: tuple[int, int]) -> tuple[int, int]:
    """A tall strip of few columns: a dense matrix read one column at a time."""
    rows, columns = shape
    chunk_rows = max(1, min(rows, WHOLE_CHUNK_ELEMENTS))
    return (chunk_rows, max(1, min(columns, WHOLE_CHUNK_ELEMENTS // chunk_rows)))


class _Writer:
    """Creates the groups and arrays of one store, HDF5 or Zarr.

    In an ``.h5ad`` every group and array is declared first, with its
    attributes, and the values are written afterwards (:meth:`write`): HDF5
    keeps a string attribute's value (every ``encoding-type``) in a global
    heap collection, which it places among the values written so far, or in
    a collection that strings written before it left room in. Declared before
    any value, the attributes' strings sit together right after the metadata,
    where opening reads them with the tree, instead of among ``obs``.
    """

    def __init__(self, root: Any, kind: Literal["h5ad", "zarr"], zarr_format: int | None):
        self.root = root
        self.kind = kind
        self.zarr_format = zarr_format
        #: The values of an ``.h5ad`` still to write: what opening reads, then
        #: the rest in the order declared. Each is a dataset and a function
        #: giving its values, so a converted copy exists only while written.
        self._opening_fills: list[tuple[Any, Callable[[], Any]]] = []
        self._fills: list[tuple[Any, Callable[[], Any]]] = []
        if kind == "zarr":
            if zarr_format == 2:
                from numcodecs import Zlib

                self.codec: Any = Zlib(level=DEFLATE_LEVEL)
            else:
                self.codec = _reproducible_gzip(DEFLATE_LEVEL)

    def _zarr_codecs(self, itemsize: int, *, shuffle: bool = True) -> dict[str, Any]:
        """The filters and compressors of a Zarr array of ``itemsize``-byte values.

        A multi-byte value is byte-shuffled before it is deflated (unless
        ``shuffle`` is false: the expression values, see the module
        docstring), as HDF5's shuffle filter does in an ``.h5ad``: numcodecs'
        ``Shuffle`` filter before zlib in Zarr format 2, the ``blosc`` codec
        (byte shuffle and zlib) in format 3.
        """
        if itemsize <= 1 or not shuffle:
            return {"filters": None, "compressors": self.codec}
        if self.zarr_format == 2:
            from numcodecs import Shuffle

            return {"filters": [Shuffle(elementsize=itemsize)], "compressors": self.codec}
        return {"filters": None, "compressors": _blosc_codec(itemsize)}

    def _fill(self, dataset: Any, values: Callable[[], Any], *, opening: bool) -> None:
        """Write a declared array's values: now in a Zarr store, at the end in an ``.h5ad``."""
        if self.kind == "h5ad":
            (self._opening_fills if opening else self._fills).append((dataset, values))
        else:
            dataset[...] = values()

    def group(self, parent: Any, key: str, attrs: Mapping[str, Any]) -> Any:
        group = parent.create_group(key)
        group.attrs.update(dict(attrs))
        return group

    def dict_group(self, parent: Any, key: str) -> Any:
        return self.group(parent, key, _encoding("dict", "0.1.0"))

    def array(
        self,
        parent: Any,
        key: str,
        values: np.ndarray,
        chunks: tuple[int, ...],
        *,
        dtype: Any = None,
        encoded: bool = True,
        shuffle: bool = True,
        opening: bool = False,
    ) -> None:
        """One numeric or boolean array, compressed in ``chunks``.

        ``dtype`` converts the values as they are written. ``encoded`` adds
        the ``array`` encoding attributes, which every array element carries;
        the three arrays inside a sparse matrix carry none. ``shuffle`` is
        false for the expression values alone (module docstring). ``opening``
        marks an array opening reads (:meth:`write`).
        """
        dtype = np.dtype(values.dtype if dtype is None else dtype)

        def converted() -> np.ndarray:
            return np.ascontiguousarray(values, dtype=dtype)

        if self.kind == "h5ad":
            if values.size:
                # The values are converted to the dataset's type chunk by
                # chunk as they are written (:func:`_write_chunks_in_parallel`),
                # so no whole converted copy is made.
                dataset = parent.create_dataset(
                    key,
                    shape=values.shape,
                    dtype=dtype,
                    chunks=chunks,
                    compression="gzip",
                    compression_opts=DEFLATE_LEVEL,
                    shuffle=shuffle and dtype.itemsize > 1,
                )
                self._fill(dataset, lambda: np.asarray(values), opening=opening)
            else:
                # HDF5 refuses a chunk of an empty dataset.
                dataset = parent.create_dataset(key, shape=values.shape, dtype=dtype)
        else:
            dataset = parent.create_array(
                key,
                shape=values.shape,
                dtype=dtype,
                chunks=chunks,
                fill_value=np.zeros((), dtype=dtype).item(),
                **self._zarr_codecs(dtype.itemsize, shuffle=shuffle),
            )
            if values.size:
                self._fill(dataset, converted, opening=opening)
        if encoded:
            dataset.attrs.update(_encoding("array", "0.2.0"))

    def element(
        self, parent: Any, key: str, value: Any, chunks: tuple[int, ...], *, opening: bool = False
    ) -> None:
        """A string array (as anndata encodes it), or a nullable integer or boolean array.

        A nullable array is a group of two arrays, ``values`` and ``mask``,
        each laid out as every other array of the store. In an ``.h5ad`` a
        string array is declared here and written with the other values,
        exactly as anndata writes one (variable-length UTF-8 strings,
        ``string-array`` 0.2.0). A Zarr format 3 store's is written here as
        anndata writes one, through the array just created: anndata writes
        the values through the array reopened from its metadata, which
        swaps the store's reproducible gzip codec for numcodecs' own. A
        Zarr format 2 store's is written by anndata.
        """
        if isinstance(value, pd.arrays.IntegerArray | pd.arrays.BooleanArray):
            kind = "nullable-integer" if isinstance(value, pd.arrays.IntegerArray) else "nullable-boolean"
            group = self.group(parent, key, _encoding(kind, "0.1.0"))
            self.array(group, "values", value.to_numpy(dtype=value.dtype.numpy_dtype, na_value=0), chunks)
            self.array(group, "mask", np.asarray(value.isna()), chunks)
            return
        if self.kind == "h5ad":
            import h5py

            strings = np.asarray(value, dtype=object)
            # HDF5 refuses a chunk of an empty dataset, and has nothing to
            # compress in one.
            options: dict[str, Any] = (
                {"chunks": chunks, "compression": "gzip", "compression_opts": DEFLATE_LEVEL}
                if len(strings)
                else {}
            )
            dataset = parent.create_dataset(key, shape=strings.shape, dtype=h5py.string_dtype(), **options)
            dataset.attrs.update(_encoding("string-array", "0.2.0"))
            if len(strings):
                self._fill(dataset, lambda: strings, opening=opening)
            return
        if self.zarr_format == 3:
            from zarr.core.dtype import VariableLengthUTF8

            strings = np.asarray(value, dtype=object)
            dataset = parent.create_array(
                key, shape=strings.shape, dtype=VariableLengthUTF8(), chunks=chunks, compressors=self.codec
            )
            dataset.attrs.update(_encoding("string-array", "0.2.0"))
            if len(strings):
                dataset[...] = strings
            return
        from anndata.io import write_elem

        # An empty Zarr array declares the store's codec too: a reader
        # checks the codec before it sees there is no chunk to decode.
        write_elem(parent, key, value, dataset_kwargs={"chunks": chunks, "compressor": self.codec})

    def scalar(self, parent: Any, key: str, value: str | bool | int | float) -> None:
        """One text or number of ``uns``, a 0-dimensional array, which opening reads.

        In an ``.h5ad`` it is stored as anndata stores one (uncompressed: HDF5
        compresses no scalar), declared here and written with the other
        values. anndata's own writer stores a Zarr scalar with zarr-python's
        default codec (Blosc-LZ4 in Zarr v2, Zstandard in v3), not the
        store's, so a Zarr scalar is written here with the store's codec, as
        every other array.
        """
        if self.kind == "h5ad":
            import h5py

            if isinstance(value, str):
                data = np.array(value, dtype=h5py.string_dtype(encoding="utf-8"))
                encoding = _encoding("string", "0.2.0")
            else:
                data = np.array(value)
                encoding = _encoding("numeric-scalar", "0.2.0")
            dataset = parent.create_dataset(key, shape=(), dtype=data.dtype)
            dataset.attrs.update(encoding)
            self._fill(dataset, lambda: data, opening=True)
            return
        codecs: dict[str, Any]
        if isinstance(value, str):
            from zarr.core.dtype import VariableLengthUTF8

            codecs = {"filters": None, "compressors": self.codec}
            if self.zarr_format == 2:
                from numcodecs import VLenUTF8

                codecs["filters"] = [VLenUTF8()]
            dtype: Any = VariableLengthUTF8()
            fill_value: Any = ""
            encoding = _encoding("string", "0.2.0")
        else:
            data = np.asarray(value)
            dtype = data.dtype
            fill_value = np.zeros((), dtype=dtype).item()
            encoding = _encoding("numeric-scalar", "0.2.0")
            codecs = self._zarr_codecs(data.dtype.itemsize)
        dataset = parent.create_array(key, shape=(), dtype=dtype, fill_value=fill_value, **codecs)
        dataset[...] = np.asarray(value)
        dataset.attrs.update(encoding)

    # -- composite elements ---------------------------------------------------

    def uns(self, parent: Any, key: str, value: Any, *, opening: bool = True) -> None:
        """``uns`` and every mapping, list and value in it, laid out as the rest.

        A list is an array read whole, chunked and compressed as an obs
        column is; a text or number is a scalar (:meth:`scalar`). Opening
        reads all of ``uns`` except the stored outlier quantiles, an array
        per categorical obs column read only when outliers are filtered.
        """
        if isinstance(value, Mapping):
            group = self.dict_group(parent, key)
            for name, item in value.items():
                self.uns(group, name, item, opening=opening and not (key == "cellucid" and name == "outlierQuantiles"))
        elif isinstance(value, np.ndarray):
            chunks = _vector_chunks(len(value), WHOLE_CHUNK_ELEMENTS)
            if value.dtype.kind == "O":
                self.element(parent, key, value, chunks, opening=opening)
            else:
                self.array(parent, key, value, chunks, opening=opening)
        else:
            self.scalar(parent, key, value)

    def column(self, parent: Any, key: str, value: Any, *, opening_categories: bool = False) -> None:
        """One dataframe column; ``opening_categories``: a categorical's labels are read as the file opens."""
        length = len(value)
        chunks = _vector_chunks(length, WHOLE_CHUNK_ELEMENTS)
        if isinstance(value, pd.Categorical):
            group = self.group(
                parent, key, {**_encoding("categorical", "0.2.0"), "ordered": bool(value.ordered)}
            )
            self.array(group, "codes", np.asarray(value.codes), chunks)
            categories = value.categories.to_numpy()
            category_chunks = _vector_chunks(len(categories), WHOLE_CHUNK_ELEMENTS)
            if categories.dtype.kind == "O":
                self.element(group, "categories", categories.astype(object), category_chunks, opening=opening_categories)
            else:
                self.array(group, "categories", categories, category_chunks, opening=opening_categories)
        elif isinstance(value, np.ndarray) and value.dtype.kind != "O":
            self.array(parent, key, value, chunks)
        else:
            self.element(parent, key, value, chunks)

    def dataframe(
        self,
        parent: Any,
        key: str,
        index: np.ndarray,
        columns: Mapping[str, Any],
        *,
        opening: bool = False,
        opening_categories: bool = False,
    ) -> None:
        """A dataframe; ``opening``: its index is read when the file opens (the gene
        names); ``opening_categories``: so are its categorical columns' labels (the
        web app lists every obs field's categories as it opens)."""
        group = self.group(
            parent,
            key,
            {
                **_encoding("dataframe", "0.2.0"),
                "_index": "_index",
                "column-order": list(columns),
            },
        )
        self.element(group, "_index", index, _vector_chunks(len(index), WHOLE_CHUNK_ELEMENTS), opening=opening)
        for name, value in columns.items():
            self.column(group, name, value, opening_categories=opening_categories)

    def sparse_pointers(
        self,
        parent: Any,
        key: str,
        matrix: sparse.csc_matrix | sparse.csr_matrix,
        *,
        opening: bool = False,
    ) -> Any:
        """A sparse matrix's group and its ``indptr``; :meth:`sparse_entries` writes the rest."""
        kind = "csc_matrix" if matrix.format == "csc" else "csr_matrix"
        group = self.group(
            parent, key, {**_encoding(kind, "0.1.0"), "shape": [int(extent) for extent in matrix.shape]}
        )
        self.array(
            group,
            "indptr",
            matrix.indptr,
            _vector_chunks(len(matrix.indptr), WHOLE_CHUNK_ELEMENTS),
            dtype=_index_dtype(matrix),
            encoded=False,
            opening=opening,
        )
        return group

    def sparse_entries(
        self, group: Any, matrix: sparse.csc_matrix | sparse.csr_matrix, *, sliced: bool
    ) -> None:
        """A sparse matrix's ``data`` and ``indices``.

        Each array's chunks lie together, in order: a whole-matrix read is a
        few long ranges. (Alternating the two arrays' chunks put a gene's two
        ranges side by side, but made reading one array alone fetch the other
        too: Suo's 1,000 cells over HTTP/1.1 read 113 MB instead of 74 and took
        6.3 s instead of 4.9.) The values of the expression matrix (``sliced``)
        are not shuffled (module docstring).
        """
        chunks = _vector_chunks(matrix.nnz, SLICED_CHUNK_ELEMENTS if sliced else WHOLE_CHUNK_ELEMENTS)
        self.array(group, "data", matrix.data, chunks, encoded=False, shuffle=not sliced)
        self.array(group, "indices", matrix.indices, chunks, dtype=_position_dtype(matrix), encoded=False)

    def write(self, elements: Elements) -> None:
        """Every element, in the order an ``.h5ad`` then holds their values.

        In an ``.h5ad`` everything is declared first (class docstring), so
        the metadata and every attribute come first; then what opening
        reads -- ``uns`` (but its stored outlier quantiles), the gene names
        (the ``var`` index), the pointers of a sparse ``X`` and the labels of
        every categorical ``obs`` column -- so it sits together right after
        them; then the other values in the order declared, the bulk of ``X``
        last. (A Zarr store keeps each chunk in
        its own file or entry, in any order, and is written as declared.)
        """
        root = self.root
        root.attrs.update(_encoding("anndata", "0.1.0"))
        self.uns(root, "uns", elements.uns)
        self.dataframe(root, "var", elements.var_index, elements.var, opening=True)
        x = elements.x
        x_group = (
            None if x is None or isinstance(x, np.ndarray) else self.sparse_pointers(root, "X", x, opening=True)
        )
        self.dataframe(root, "obs", elements.obs_index, elements.obs, opening_categories=True)
        obsm = self.dict_group(root, "obsm")
        for key, values in elements.obsm.items():
            self.array(obsm, key, values, _row_chunks(values.shape))
        obsp = self.dict_group(root, "obsp")
        for key, matrix in elements.obsp.items():
            self.sparse_entries(self.sparse_pointers(obsp, key, matrix), matrix, sliced=False)
        self.dict_group(root, "varm")
        self.dict_group(root, "varp")
        self.dict_group(root, "layers")
        if isinstance(x, np.ndarray):
            self.array(root, "X", x, _strip_chunks(x.shape), shuffle=False)
        elif x is not None:
            self.sparse_entries(x_group, x, sliced=True)
        self.flush()

    def flush(self) -> None:
        """Write the values of an ``.h5ad`` declared so far: what opening reads, then the rest.

        Each array is flushed as soon as it is written: HDF5 keeps chunks
        smaller than its chunk cache until their dataset is flushed or
        closed, and places them only then. (A Zarr store has nothing left
        to write.)
        """
        for dataset, produce in [*self._opening_fills, *self._fills]:
            values = produce()
            if not _write_chunks_in_parallel(dataset, values):
                dataset[...] = values
            dataset.flush()
        self._opening_fills.clear()
        self._fills.clear()


def _chunk_bytes(values: np.ndarray, start: tuple[int, ...], chunks: tuple[int, ...], dtype: np.dtype, shuffle: bool) -> bytes:
    """One HDF5 chunk as its gzip (and shuffle) filters store it.

    The chunk's values in C order, at the dataset's type, padded with zeros
    (the fill value) to the whole chunk at the array's edge, as HDF5 stores an
    edge chunk; byte-shuffled as HDF5's shuffle filter does (every value's
    first byte, then every second byte, ...); deflated as HDF5's deflate
    filter does (a zlib stream at :data:`DEFLATE_LEVEL`).
    """
    import zlib

    selection = tuple(slice(first, first + size) for first, size in zip(start, chunks, strict=True))
    block = np.asarray(values[selection], dtype=dtype)
    if block.shape != chunks:
        padded = np.zeros(chunks, dtype=dtype)
        padded[tuple(slice(0, size) for size in block.shape)] = block
        block = padded
    raw = np.ascontiguousarray(block).view(np.uint8)
    if shuffle and dtype.itemsize > 1:
        raw = np.ascontiguousarray(raw.reshape(-1, dtype.itemsize).T)
    return zlib.compress(raw.tobytes(), DEFLATE_LEVEL)


def _write_chunks_in_parallel(dataset: Any, values: Any) -> bool:
    """Write a numeric ``.h5ad`` array's chunks, compressed on every core.

    HDF5 compresses the chunks of a write one after the other on the calling
    thread, which made deflating the values most of the time ``prepare()``
    takes (93% of it on 200,000 cells). Each chunk is compressed here, as
    HDF5's own filters would (:func:`_chunk_bytes`), by a pool of threads --
    zlib releases the interpreter lock while it deflates -- and stored with
    HDF5's direct chunk write, in grid order, so the file's layout is what
    HDF5 writes. Returns False for an array this does not apply to (text, a
    scalar, an unchunked or otherwise filtered array), which the caller
    writes as before.
    """
    from concurrent.futures import ThreadPoolExecutor
    from itertools import product

    if not isinstance(values, np.ndarray) or values.dtype.kind not in "biuf":
        return False
    chunks = dataset.chunks
    if chunks is None or dataset.compression != "gzip" or dataset.compression_opts != DEFLATE_LEVEL:
        return False
    if dataset.scaleoffset is not None or dataset.fletcher32 or values.shape != dataset.shape:
        return False
    dtype = np.dtype(dataset.dtype)
    shuffle = bool(dataset.shuffle)
    grid = [range(0, length, size) for length, size in zip(dataset.shape, chunks, strict=True)]
    starts = product(*grid)
    workers = max(1, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        # A window of chunks in flight, written in order as each completes:
        # memory holds a few chunks per thread, never the compressed array.
        pending: list[tuple[tuple[int, ...], Any]] = []
        for start in starts:
            pending.append((start, pool.submit(_chunk_bytes, values, start, chunks, dtype, shuffle)))
            if len(pending) >= 4 * workers:
                first, future = pending.pop(0)
                dataset.id.write_direct_chunk(first, future.result(), 0)
        for first, future in pending:
            dataset.id.write_direct_chunk(first, future.result(), 0)
    return True


def _index_dtype(matrix: sparse.csc_matrix | sparse.csr_matrix) -> type[np.signedinteger[Any]]:
    """``indptr``: int32 unless the matrix holds 2^31 stored values or more."""
    return np.int64 if matrix.nnz >= 2**31 else np.int32


def _position_dtype(matrix: sparse.csc_matrix | sparse.csr_matrix) -> type[np.signedinteger[Any]]:
    """``indices``: int32 unless the axis they index (a CSC matrix's rows) has 2^31 entries.

    Decided by the axis, not by the stored-value count: a matrix of more than
    2^31 values over a few million cells keeps 4-byte row positions, half the
    bytes of the second-largest array on disk and over the network, where
    scipy would widen them with ``indptr``. Every reader takes each array's
    own type (anndata, the web app's reader).
    """
    axis = matrix.shape[0] if matrix.format == "csc" else matrix.shape[1]
    return np.int64 if axis > np.iinfo(np.int32).max else np.int32


# =============================================================================
# The store, and its placement
# =============================================================================


def _write_root(
    root: Any, kind: Literal["h5ad", "zarr"], zarr_format: int | None, content: Elements | StoreContent
) -> None:
    if isinstance(content, Elements):
        _Writer(root, kind, zarr_format).write(content)
    else:
        content.write(root, kind, zarr_format)


#: The user block an ``.h5mu`` starts with: mudata reads its first bytes
#: (``MuData``) to recognise the file, and the HDF5 signature follows it.
H5MU_USERBLOCK_BYTES = 512


def _h5mu_header() -> bytes:
    """The text of an ``.h5mu``'s user block, as mudata writes it, naming this writer."""
    from importlib.metadata import version

    return f"MuData (format-version=0.1.0;creator=cellucid;creator-version={version('cellucid')})".encode()


def _write_h5ad(path: Path, content: Elements | StoreContent, *, mudata: bool = False) -> None:
    """An ``.h5ad`` that opens with one read of its start.

    HDF5 gathers metadata -- object headers, attributes, chunk indexes -- in
    blocks of ``meta_block_size`` (:data:`H5AD_META_BLOCK_BYTES`), the first
    one at the start of the file; with HDF5's default 2 KiB blocks it is
    scattered between the arrays, and opening a file waits one round trip per
    level of its tree. Every group and array is declared before any value is
    written (:class:`_Writer`), so the attributes' strings follow the metadata
    block, and then come the values opening reads: ``uns`` (but its outlier
    quantiles), the gene names, the pointers of ``X`` and the labels of every
    categorical ``obs`` column -- 431 KB in all on the dataset above, which
    the web app reads with its first request of the file (``HEAD_BYTES``,
    512 KiB).
    File format 1.10 (``libver``) stores a single chunk's address in its
    dataset's header and the chunks of a fixed-size array in one fixed array
    index, so no chunk needs an index read of its own; HDF5 1.10 and later
    (h5py, anndata, rhdf5, hdf5r, h5wasm) read it. (HDF5's paged file space
    gathered the metadata too, but scattered each array's chunks over pages
    shared with the others: a whole-matrix read over HTTP/1.1 took 224
    requests instead of 33.)

    An ``.h5mu`` (``mudata``) is the same HDF5 file after a user block of
    :data:`H5MU_USERBLOCK_BYTES` holding :func:`_h5mu_header`.
    """
    import h5py

    with h5py.File(
        path,
        "w",
        libver=("v110", "v110"),
        meta_block_size=H5AD_META_BLOCK_BYTES,
        userblock_size=H5MU_USERBLOCK_BYTES if mudata else 0,
    ) as file:
        _write_root(file, "h5ad", None, content)
    if mudata:
        header = _h5mu_header()
        with path.open("r+b") as file:
            file.write(header + b"\0" * (H5MU_USERBLOCK_BYTES - len(header)))


def _write_zarr_folder(path: Path, content: Elements | StoreContent, zarr_format: Literal[2, 3]) -> None:
    import anndata
    import zarr

    root = zarr.open_group(str(path), mode="w", zarr_format=zarr_format)
    # Every array unsharded: one chunk is one file, read with one request.
    with anndata.settings.override(auto_shard_zarr_v3=False):
        _write_root(root, "zarr", zarr_format, content)
    # A browser cannot list a folder: the consolidated metadata is how it
    # learns every array and group with one request. zarr-python warns that
    # consolidated metadata is not yet in the v3 specification; the web app
    # and anndata both read it.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Consolidated metadata is currently not part")
        zarr.consolidate_metadata(str(path), zarr_format=zarr_format)


def _zip_entries(folder: Path) -> list[Path]:
    """The files of a store folder in the order a ``.zarr.zip`` holds them.

    The chunks first, then the metadata documents, the consolidated metadata
    last: a reader opening the archive reads its end -- the central directory
    -- and with it the documents that describe the store, so it opens with
    one read. The chunks of an array follow one another in their grid order.
    """
    files = [path for path in folder.rglob("*") if path.is_file()]

    def order(path: Path) -> tuple[Any, ...]:
        relative = path.relative_to(folder)
        if path.name in ZARR_METADATA_NAMES:
            consolidated = len(relative.parts) == 1
            return (1, consolidated, relative.parts)
        # A chunk: its array's path, then its grid position in numbers
        # ("X/data/10" after "X/data/9"). Format 3 keys a chunk "c/<i>/<j>"
        # under its array ("c" alone for a scalar), format 2 "<i>.<j>".
        parts = relative.parts
        for at in range(len(parts) - 1, -1, -1):
            if parts[at] == "c" and all(part.isdigit() for part in parts[at + 1 :]):
                array, position = parts[:at], tuple(int(part) for part in parts[at + 1 :])
                break
        else:
            array, position = parts[:-1], tuple(int(part) for part in parts[-1].split("."))
        return (0, array, position)

    return sorted(files, key=order)


def _zip_folder(folder: Path, archive: Path) -> None:
    """Pack a store folder into one stored ZIP whose top level is the store.

    Entries are stored, not deflated: every chunk is compressed already, and a
    stored entry is read with one range request at its offset. They are in
    :func:`_zip_entries` order. The bytes depend on the entries alone, as
    cellucid-r's writer makes them: every timestamp is 1980-01-01 00:00 and
    every entry a 0644 regular file, never the time or permissions the folder
    was written with.
    """
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as file:
        for source in _zip_entries(folder):
            entry = zipfile.ZipInfo(source.relative_to(folder).as_posix(), _ZIP_EPOCH)
            entry.create_system = 3
            entry.external_attr = 0o100644 << 16
            entry.file_size = source.stat().st_size
            with source.open("rb") as reader, file.open(entry, "w") as writer:
                shutil.copyfileobj(reader, writer, _COPY_BYTES)


def _remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


@contextmanager
def _staged(target: OutputTarget) -> Iterator[Path]:
    """A temporary sibling of the target, moved into place once complete."""
    parent = target.path.parent
    parent.mkdir(parents=True, exist_ok=True)
    tag = f".{target.path.name}.{os.getpid()}"
    staging = parent / f"{tag}.partial"
    replaced = parent / f"{tag}.replaced"
    for leftover in (staging, replaced):
        _remove(leftover)
    try:
        yield staging
        had_target = target.path.exists() or target.path.is_symlink()
        if had_target:
            os.replace(target.path, replaced)
        try:
            os.replace(staging, target.path)
        except BaseException:
            if had_target:
                os.replace(replaced, target.path)
            raise
        if had_target:
            _remove(replaced)
    finally:
        _remove(staging)


def write_store(target: OutputTarget, content: Elements | StoreContent) -> None:
    """Write ``content`` to ``target``, replacing what was there only once complete."""
    with _staged(target) as staging:
        if target.kind in ("h5ad", "h5mu"):
            _write_h5ad(staging, content, mudata=target.kind == "h5mu")
            return
        assert target.zarr_format is not None
        if not target.zipped:
            _write_zarr_folder(staging, content, target.zarr_format)
            return
        folder = staging.with_name(f"{staging.name}.zarr")
        _remove(folder)
        try:
            _write_zarr_folder(folder, content, target.zarr_format)
            _zip_folder(folder, staging)
        finally:
            _remove(folder)
