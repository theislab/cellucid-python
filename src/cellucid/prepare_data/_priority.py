"""
Priority order: the row order of a Cellucid-ready file.

Every cell has a priority, its ``uniform`` rank in the subsampling
specification (:mod:`cellucid.subsample`): ``rank(seed, p)``, a bijection of
the uint32 positions for one seed, so no two cells share one. A file in
priority order stores its cells by ascending rank: file row ``r`` holds the
cell of the ``r``-th lowest rank. The ``k`` cells a ``uniform`` subsample of
size ``k`` keeps with the file's seed are then its first ``k`` rows, in every
array with one row per cell -- ``obs``, ``obsm``, a dense ``X``, both axes of
``obsp`` -- and, in a gene-major ``X`` whose entries ascend by row within
each gene, the first entries of every gene. A subsample reads one range of
each array and never a byte of a cell it left out.

What the file records (the Cellucid-ready specification):

- ``uns['cellucid']['priority'] = {seed, specVersion, nObs, levels}``: the
  seed of the ranks, the subsampling specification's version they follow,
  the cell count the file was written with, and the prefix sizes ``levels``
  (ascending, :func:`prefix_levels` of ``nObs``) that
  ``varm['cellucid_prefix_counts']`` counts at;
- ``obs['cellucid_original_position']``: each row's position in the order the
  cells were given (uint32), a permutation of ``0 .. n_obs - 1``;
- ``varm['cellucid_prefix_counts']`` beside a sparse gene-major ``X``:
  ``n_vars x len(levels)`` uint32, column ``j`` holding how many of each
  gene's entries lie in rows below ``levels[j]`` (:func:`prefix_counts`).

:func:`cellucid.restore_original_order` puts a file read with anndata back in
the order its cells were given.

The record is the claim, verified before anything relies on it
(:func:`stored_order`, the one rule every reader shares with the web app):
its ``nObs`` is the object's cell count, its ``levels`` those of that count,
and the position column, when there, is the seed's order. anndata keeps
``uns`` through a subset, a concatenation or a reorder, so an object changed
so after it was read carries a record that no longer describes it: it is
then taken as cells Cellucid did not write -- each row its own original
position -- with a warning, never refused and never misread. A position
column without a record names nothing.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy import sparse

from .._plural import plural
from ..subsample import SUBSAMPLE_SPEC_VERSION, _rank_argsort, _ranks, _require_uint32
from ._gene_major import counts_below

if TYPE_CHECKING:
    from ._mudata import Modality

#: The obs column holding each row's original position.
ORIGINAL_POSITION_COLUMN = "cellucid_original_position"

#: The varm key of a sparse gene-major ``X``'s prefix counts.
PREFIX_COUNTS_KEY = "cellucid_prefix_counts"

#: The prefix sizes counted step down from ``n_obs`` by this factor (as an
#: exact fraction: each level is ``floor(previous * 3 / 4)``). Measured on
#: the Suo atlas (561,947 cells, 5,103 genes) and a 1,000,000-cell synthetic
#: file with gene-major ``X`` in 2^16-value chunks: a gene's prefix at 1, 5,
#: 10 and 30% of the cells, read up to the next level, fetched at most 1.01x
#: the bytes of its exact prefix, against 1.04-1.09x with halving levels and
#: 1.00x with 8/7; the table costs 0.09-0.27% of ``X`` (43-45 columns,
#: 1.6-4.1 KB per column stored), against 0.04-0.11% and 0.19-0.58%.
LEVEL_RATIO = (3, 4)


def prefix_levels(n_obs: int) -> np.ndarray:
    """The prefix sizes counted, ascending: ``n_obs`` stepped down by :data:`LEVEL_RATIO` to 1.

    ``n_obs`` itself is not a level (a gene's whole range is its entry count,
    ``X/indptr``); a level repeated by the rounding is kept once.
    """
    numerator, denominator = LEVEL_RATIO
    levels: list[int] = []
    level = int(n_obs)
    while True:
        level = level * numerator // denominator
        if level < 1 or (levels and level == levels[-1]):
            break
        levels.append(level)
    return np.asarray(levels[::-1], dtype=np.int64)


def prefix_counts(matrix: Any, levels: np.ndarray) -> np.ndarray:
    """``varm['cellucid_prefix_counts']``: each gene's entries in rows below each level.

    ``matrix`` is gene-major (``indptr``, ``indices``; a CSC matrix or a
    :class:`._gene_major.GeneMajor`), its rows ascending within each gene or
    not: each entry is counted toward the levels above its row.
    """
    pointers = np.asarray(matrix.indptr, dtype=np.int64)
    return counts_below(np.asarray(matrix.indices), np.diff(pointers), np.asarray(levels, dtype=np.int64))


def _take(values: Any, rows: np.ndarray) -> Any:
    """``values`` (one row per cell) at ``rows``: arrays, categoricals and pandas nullable arrays alike."""
    if isinstance(values, np.ndarray):
        return values[rows]
    if hasattr(values, "take"):
        return values.take(rows)
    return np.asarray(values)[rows]


def _require_positions(values: Any, n_cells: int) -> np.ndarray:
    """A given ``obs['cellucid_original_position']``: a permutation of ``0 .. n_cells - 1``."""
    array = np.asarray(values)
    if array.dtype.kind not in "iu" or array.shape != (n_cells,):
        raise ValueError(
            f"obs[{ORIGINAL_POSITION_COLUMN!r}] holds each cell's original position, a whole "
            f"number per cell; this one is {array.dtype} of shape {array.shape}. Rename the "
            "column, or drop it to store the cells' order as given."
        )
    inside = (array >= 0) & (array < n_cells)
    if not inside.all() or len(np.unique(array)) != n_cells:
        raise ValueError(
            f"obs[{ORIGINAL_POSITION_COLUMN!r}] must hold every position 0 to {n_cells - 1} once: "
            "it records the order the cells were first given in, and a subset or a "
            "concatenation of files no longer has one. Drop the column to store the cells' "
            "order as given."
        )
    return array.astype(np.int64)


def _graph(graph: sparse.csr_matrix, order: np.ndarray, new_row: np.ndarray) -> sparse.csr_matrix:
    """A cells x cells graph with both axes in file order, canonical CSR."""
    rows = graph[order]
    result = sparse.csr_matrix(
        (rows.data, new_row[rows.indices].astype(rows.indices.dtype), rows.indptr), shape=graph.shape
    )
    result.has_sorted_indices = False
    result.sort_indices()
    return result


def _inverse(order: np.ndarray) -> np.ndarray:
    inverse = np.empty(len(order), dtype=np.int64)
    inverse[order] = np.arange(len(order), dtype=np.int64)
    return inverse


@dataclass(frozen=True)
class CellPriority:
    """The row order of one dataset's cells as written, and how each element is put in it.

    ``order[r]`` is the given row file row ``r`` holds, ``new_row`` the
    inverse (the file row of each given row), ``original[r]`` its original
    position. The given rows are the cells in their original order, or an
    object in a verified stored order (:func:`stored_order`): the two
    permutations, back to the original order and on to the file's, are one
    (:meth:`of`), and an object already in the file's order is written as it
    is (``identity``: every element taken without a copy). ``seed`` is the
    priority order's seed; None writes the cells in their original order, with
    no priority record and no position column.
    """

    seed: int | None
    order: np.ndarray
    original: np.ndarray
    new_row: np.ndarray
    #: Every file row is the given row of its number: nothing to reorder.
    identity: bool = False

    @classmethod
    def of(cls, n_cells: int, seed: int | None, given: np.ndarray | None = None) -> CellPriority | None:
        """The order of ``n_cells`` cells for ``seed`` (None: their original order).

        ``given`` is each given row's original position, for an object in a
        verified stored order (:attr:`StoredOrder.positions`); None when each
        row is its own. None when there is nothing to do: no seed, and the
        rows in their original order.
        """
        if seed is None and given is None:
            return None
        n_cells = _require_uint32(n_cells, "the cell count")
        positions = np.arange(n_cells, dtype=np.int64)
        # Ranks never tie for one seed: any sort gives the one order.
        original = positions if seed is None else _rank_argsort(_ranks(positions, seed)).astype(np.int64)
        order = original if given is None else _inverse(np.asarray(given, dtype=np.int64))[original]
        return cls(
            seed=seed, order=order, original=original, new_row=_inverse(order), identity=bool(np.array_equal(order, positions))
        )

    @property
    def recorded(self) -> bool:
        """Whether the file records this order (a priority order: a seed)."""
        return self.seed is not None

    @property
    def rows_moved(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        """``(order, new_row)``, or ``(None, None)`` when no row moves (``identity``)."""
        return (None, None) if self.identity else (self.order, self.new_row)

    def rows(self, values: Any) -> Any:
        return values if self.identity else _take(values, self.order)

    def rows_of_mapping(self, arrays: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        return {key: self.rows(values) for key, values in arrays.items()}

    def obs_columns(self, columns: dict[str, Any]) -> dict[str, Any]:
        """Every column in file order, and (a priority order) the original positions last."""
        placed = {key: self.rows(values) for key, values in columns.items() if key != ORIGINAL_POSITION_COLUMN}
        if self.recorded:
            placed[ORIGINAL_POSITION_COLUMN] = self.original.astype(np.uint32)
        return placed

    def graph(self, graph: sparse.csr_matrix) -> sparse.csr_matrix:
        return graph if self.identity else _graph(graph, self.order, self.new_row)

    def record(self, levels: np.ndarray | None) -> dict[str, Any]:
        """``uns['cellucid']['priority']``."""
        assert self.seed is not None
        return {
            "seed": int(self.seed),
            "specVersion": SUBSAMPLE_SPEC_VERSION,
            "nObs": len(self.order),
            "levels": np.asarray(levels, dtype=np.int64),
        }

    def modality(self, modality: Modality, levels: np.ndarray | None) -> Modality:
        """A MuData modality with its rows ordered by its cells' file rows.

        Its cell map follows the cells (``obsmap`` in file order, renumbered
        to the modality's new rows); its rows' original rows are its own
        given ones (``Modality.given``: a modality of an object in a verified
        stored order) or the rows themselves. In a priority order, its prefix
        levels are its rows of cells below each global level, so a global
        prefix of the cells is a prefix of the modality, and its own
        ``obs['cellucid_original_position']`` holds each row's original row in
        the modality, so its order is restored too.
        """
        placed = np.asarray(modality.obsmap)[self.order]
        measured = placed > 0
        mod_order = placed[measured].astype(np.int64) - 1
        mod_new = _inverse(mod_order)
        obsmap = np.zeros(len(placed), dtype=np.uint32)
        obsmap[measured] = (mod_new[placed[measured].astype(np.int64) - 1] + 1).astype(np.uint32)
        below = np.concatenate([[0], np.cumsum(measured)])
        original = mod_order if modality.given is None else np.asarray(modality.given, dtype=np.int64)[mod_order]
        moved = not np.array_equal(mod_order, np.arange(len(mod_order)))
        obs = {key: _take(values, mod_order) if moved else values for key, values in modality.obs.items() if key != ORIGINAL_POSITION_COLUMN}
        if self.recorded:
            obs[ORIGINAL_POSITION_COLUMN] = original.astype(np.uint32)
        return replace(
            modality,
            obs_index=_take(modality.obs_index, mod_order) if moved else modality.obs_index,
            obs=obs,
            obsm={key: values[mod_order] if moved else values for key, values in modality.obsm.items()},
            obsp={key: _graph(graph, mod_order, mod_new) if moved else graph for key, graph in modality.obsp.items()},
            obsmap=obsmap,
            given=None,
            order=mod_order if moved else None,
            new_row=mod_new if moved else None,
            prefix_levels=None if levels is None else below[np.asarray(levels, dtype=np.int64)].astype(np.int64),
        )


@dataclass(frozen=True)
class StoredOrder:
    """How the rows of an object or a file name its cells (:func:`stored_order`).

    ``positions[r]`` is the original position of row ``r`` when the rows are
    in a verified priority order; None when each row is its own original
    position. ``notice`` says why a priority record the data carries is not
    used (None when there is none, or it verifies).
    """

    positions: np.ndarray | None
    notice: str | None


def _record_scalar(value: Any) -> int | None:
    """One whole number of the record (a NumPy scalar, a 0-d or 1-element array), or None."""
    array = np.asarray(value)
    if array.size != 1 or array.dtype.kind not in "iu":
        return None
    return int(array.reshape(-1)[0])


def stored_order(record: Any, n_obs: int, column: Any | None) -> StoredOrder:
    """The one rule for a priority record: is ``n_obs`` rows' order the one it claims?

    ``record`` is ``uns['cellucid']['priority']`` (None when absent) and
    ``column`` ``obs['cellucid_original_position']`` (None when absent). The
    record verifies when its seed is a uint32, its ``specVersion`` the one
    this package follows, its ``nObs`` is ``n_obs``, its ``levels`` (when
    present) are :func:`prefix_levels` of ``n_obs``, and the column (when
    present) is exactly the seed's order. Then ``positions`` is that order,
    computed from the seed (a record whose column was dropped still
    verifies). Anything else is data Cellucid did not write: ``positions``
    None (each row its own original position), and ``notice`` says why (a
    column without a record names nothing either). The web app reads a file
    by the same record checks (``priority-layout.js``) and checks the column
    where opening reads it, its first chunk of rows, computing the rest from
    the seed (reading every row would cost each opening 4-8 bytes a cell): a
    file subset, combined or reordered after it was written fails both. Only
    a permutation of the rows after the first chunk that leaves the first
    chunk in place is refused here and read in the seed's order by the web;
    no reorder anndata or Scanpy performs does that.
    """
    if record is None:
        if column is None:
            return StoredOrder(None, None)
        return StoredOrder(
            None,
            f"obs[{ORIGINAL_POSITION_COLUMN!r}] has no uns['cellucid']['priority'] beside it, so it names "
            "no order: the cells are taken in the order they are in (each row its own original "
            "position), as the web app reads such a file.",
        )

    def unused(reason: str) -> StoredOrder:
        return StoredOrder(
            None,
            f"uns['cellucid']['priority'] does not describe these cells: {reason}. They are taken "
            "in the order they are in, as cells Cellucid did not write (each row its own original "
            "position); the web app reads such a file the same way.",
        )

    if not isinstance(record, Mapping):
        return unused(f"it is a {type(record).__name__}, not a mapping")
    seed = _record_scalar(record.get("seed"))
    if seed is None or not 0 <= seed <= 0xFFFFFFFF:
        return unused("its seed is not a whole number from 0 to 4294967295")
    version = _record_scalar(record.get("specVersion"))
    if version != SUBSAMPLE_SPEC_VERSION:
        return unused(
            f"it follows version {version} of the subsampling specification, and this package "
            f"follows {SUBSAMPLE_SPEC_VERSION}"
        )
    recorded = _record_scalar(record.get("nObs"))
    if recorded is None:
        return unused("it records no cell count (nObs)")
    if recorded != n_obs:
        return unused(
            f"it was written for {recorded} {plural(recorded, 'cell')}, and there "
            f"{plural(n_obs, 'is', 'are')} {n_obs} (cells were added or removed since)"
        )
    if "levels" in record:
        levels = np.asarray(record["levels"])
        if levels.dtype.kind not in "iu" or not np.array_equal(levels.reshape(-1), prefix_levels(n_obs)):
            return unused(f"its levels are not those of {n_obs} {plural(n_obs, 'cell')}")
    positions = _rank_argsort(_ranks(np.arange(n_obs, dtype=np.int64), seed)).astype(np.int64)
    if column is not None:
        stored = np.asarray(column)
        if stored.dtype.kind not in "iu" or stored.shape != (n_obs,) or not np.array_equal(stored, positions):
            return unused(f"obs[{ORIGINAL_POSITION_COLUMN!r}] is not its seed's order (the cells were reordered since)")
    return StoredOrder(positions, None)


def stored_order_of(data: Any) -> StoredOrder:
    """:func:`stored_order` of an AnnData's or MuData's rows (its global cells)."""
    cellucid_uns = data.uns.get("cellucid") if isinstance(data.uns, Mapping) else None
    record = cellucid_uns.get("priority") if isinstance(cellucid_uns, Mapping) else None
    column = data.obs[ORIGINAL_POSITION_COLUMN].to_numpy() if ORIGINAL_POSITION_COLUMN in data.obs.columns else None
    return stored_order(record, data.n_obs, column)


def without_positions(columns: dict[str, Any], what: str, *, verified: bool = False) -> dict[str, Any]:
    """``columns`` (``obs``) without a position column ``prepare()`` was given.

    The column of an object in a verified stored order (``verified``) is the
    order its rows are in, which ``prepare()`` composes with the order it
    writes (:meth:`CellPriority.of`): it is dropped silently. Any other names
    no order (:func:`stored_order`): written again it would name cells the
    file does not hold where it says. It is dropped, with a warning.
    """
    if ORIGINAL_POSITION_COLUMN not in columns:
        return columns
    if not verified:
        warnings.warn(
            f"{what}[{ORIGINAL_POSITION_COLUMN!r}] is not written: without a priority record that describes "
            "these cells (uns['cellucid']['priority']) it names no order, and prepare() records the "
            "cells' order itself.",
            UserWarning,
            stacklevel=4,
        )
    return {key: values for key, values in columns.items() if key != ORIGINAL_POSITION_COLUMN}


def _rows_of(data: Any, rows: np.ndarray) -> Any:
    """A new in-memory object of ``data``'s ``rows``, in that order (a backed AnnData read into memory)."""
    view = data[rows]
    return view.to_memory() if getattr(data, "isbacked", False) else view.copy()


def restore_original_order(data: Any) -> Any:
    """The cells of a file written in priority order, back in the order they were given.

    A Cellucid-ready file stores its cells by priority (``prepare()``,
    ``priority_seed``), recorded in ``uns['cellucid']['priority']``, each
    row's original position in ``obs['cellucid_original_position']`` (and,
    in a MuData, each modality's rows' in its own ``obs``).
    ``anndata.read_h5ad`` and ``anndata.read_zarr`` (``mudata.read_h5mu``,
    ``mudata.read_zarr``) give the cells in that stored order; this returns a
    new object with every per-cell element -- ``obs``, ``obsm``, ``obsp``
    (both axes), ``X``, ``layers``, and in a MuData every modality's cells --
    in the original order, without the position column,
    ``uns['cellucid']['priority']`` and ``varm['cellucid_prefix_counts']``,
    which describe the stored order only; the stored outlier quantiles
    (``uns['cellucid']['outlierQuantiles']``, one value per cell) follow the
    cells. A backed AnnData (``anndata.read_h5ad(path, backed="r")``) is read
    into memory in that order.

    ``data`` is returned as it is when its rows are each its own original
    position: it has no priority record, or one that no longer describes its
    cells -- it was subset, concatenated or reordered after it was read
    (anndata keeps ``uns`` through all three) -- which a warning reports.
    That is the rule the web app reads such a file by, so a session saved on
    it applies to ``data`` as it is.

    Parameters
    ----------
    data
        An ``anndata.AnnData`` or a ``mudata.MuData``.

    Returns
    -------
    AnnData or MuData
        A new object in the original order, or ``data`` itself.

    Examples
    --------
    >>> adata = anndata.read_zarr("pbmc.zarr")
    >>> adata = cellucid.restore_original_order(adata)
    """
    from .._anndata_files import is_mudata
    from ._data import is_anndata

    if not (is_anndata(data) or is_mudata(data)):
        raise TypeError(f"data must be an anndata.AnnData or a mudata.MuData, got {type(data).__name__}.")
    order = stored_order_of(data)
    if order.notice is not None:
        warnings.warn(order.notice, UserWarning, stacklevel=2)
    if order.positions is None:
        return data
    return _restored(data, order.positions)


def _restored(data: Any, positions: np.ndarray) -> Any:
    """``data`` (rows in verified priority order, ``positions`` each row's original position) in original order."""
    from .._anndata_files import is_mudata

    order = np.argsort(positions, kind="stable")
    restored = _rows_of(data, order)
    restored.obs = restored.obs.drop(columns=[ORIGINAL_POSITION_COLUMN], errors="ignore")
    cellucid_uns = restored.uns.get("cellucid")
    if isinstance(cellucid_uns, dict):
        cellucid_uns.pop("priority", None)
        # The stored outlier quantiles: one value per cell, in the cells' order.
        fields = cellucid_uns.get("outlierQuantiles", {}).get("fields", {})
        for key, values in list(fields.items()):
            if isinstance(values, np.ndarray) and values.shape[:1] == (len(order),):
                fields[key] = values[order]
    for item in [restored, *(restored.mod.values() if is_mudata(restored) else [])]:
        if PREFIX_COUNTS_KEY in item.varm:
            del item.varm[PREFIX_COUNTS_KEY]
    if is_mudata(restored):
        # Each modality's rows in their own original order.
        for name, modality in restored.mod.items():
            if ORIGINAL_POSITION_COLUMN not in modality.obs.columns:
                continue
            rows = _require_positions(modality.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), modality.n_obs)
            modality = _rows_of(modality, np.argsort(rows, kind="stable"))
            modality.obs = modality.obs.drop(columns=ORIGINAL_POSITION_COLUMN)
            restored.mod[name] = modality
        restored.update()
        copied = f"{ORIGINAL_POSITION_COLUMN}"
        restored.obs = restored.obs.drop(columns=[column for column in restored.obs.columns if column.endswith(f":{copied}")])
    return restored
