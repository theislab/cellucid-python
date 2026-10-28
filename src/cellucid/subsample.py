"""Subsampling: which cells of a dataset a subsampled opening keeps.

One algorithm family, specified once (docs: *Subsampling specification*) and
implemented the same in the web app (``assets/js/data/subsample.js``, run while
a dataset opens), here, and in R (``cellucid_subsample_cells()``); the shared
conformance fixture ``subsample_conformance.json`` holds the three to the same
cells. :func:`subsample_cells` gives the cells Cellucid opens for a request, as
original positions (rows of ``obs``), so a notebook can take the same cells the
viewer shows, or choose cells the same way before opening.

Every hash is computed on uint32 values with numpy's wrapping uint32
arithmetic; every float operation of ``spatial`` is one IEEE 754 binary64
operation in the specified order.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from ._plural import plural

if TYPE_CHECKING:
    from numpy.typing import NDArray

__all__ = [
    "SUBSAMPLE_METHODS",
    "SUBSAMPLE_SPEC_VERSION",
    "cell_set_digest",
    "count_from_fraction",
    "decode_cell_set",
    "encode_cell_set",
    "subsample_cells",
]

#: The subsampling methods.
SUBSAMPLE_METHODS = ("uniform", "blocks", "stratified", "spatial")

#: The version of the subsampling specification this module follows: a
#: Cellucid-ready file in priority order records it
#: (``uns['cellucid']['priority']['specVersion']``), since its row order is
#: this specification's ``uniform`` rank.
SUBSAMPLE_SPEC_VERSION = 1

_M32 = 0xFFFFFFFF
_U32 = np.uint32


def _fmix32(h: NDArray[np.uint32]) -> NDArray[np.uint32]:
    h = h ^ (h >> _U32(16))
    h = h * _U32(0x85EBCA6B)
    h = h ^ (h >> _U32(13))
    h = h * _U32(0xC2B2AE35)
    return h ^ (h >> _U32(16))


def _fmix32_int(h: int) -> int:
    h &= _M32
    h ^= h >> 16
    h = (h * 0x85EBCA6B) & _M32
    h ^= h >> 13
    h = (h * 0xC2B2AE35) & _M32
    return h ^ (h >> 16)


def _seed_key(seed: int) -> int:
    return _fmix32_int(seed ^ 0x9E3779B9)


def _ranks(positions: NDArray[Any], seed: int) -> NDArray[np.uint32]:
    with np.errstate(over="ignore"):
        return _fmix32(positions.astype(np.uint32) ^ _U32(_seed_key(seed)))


def _require_uint32(value: Any, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int | np.integer):
        raise TypeError(f"{label} must be an integer; got {value!r}")
    number = int(value)
    if not minimum <= number <= _M32:
        raise ValueError(f"{label} must be an integer in [{minimum}, 2**32 - 1]; got {number}")
    return number


def count_from_fraction(fraction: float, n_obs: int) -> int:
    """The subsample size a fraction of ``n_obs`` cells stands for, exactly.

    The fraction, in ``(0, 1]``, is taken as the rational its shortest
    round-trip decimal names (``repr``, as the web app's ``String(x)`` and
    R's shortest ``%.Ng``), and the size is ``round-half-up(n_obs * fraction)``
    in integer arithmetic, at least 1 and at most ``n_obs``.
    """
    from fractions import Fraction

    n_obs = _require_uint32(n_obs, "n_obs")
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError(f"fraction must be in (0, 1]; got {fraction!r}")
    share = Fraction(repr(float(fraction)))
    k = (2 * n_obs * share.numerator + share.denominator) // (2 * share.denominator)
    return min(n_obs, max(1, k))


def _lowest(positions: NDArray[Any], k: int, seed: int) -> NDArray[np.int64]:
    """The ``k`` positions of lowest rank, ascending (ranks never tie)."""
    if k >= positions.size:
        return np.sort(positions.astype(np.int64))
    ranks = _ranks(positions, seed)
    chosen = np.argpartition(ranks, k - 1)[:k]
    return np.sort(positions[chosen].astype(np.int64))


def _uniform(n_obs: int, k: int, seed: int) -> NDArray[np.int64]:
    return _lowest(np.arange(n_obs, dtype=np.uint32), k, seed)


def _blocks(n_obs: int, k: int, seed: int, block_size: int) -> NDArray[np.int64]:
    count = -(-n_obs // block_size)
    take = min(count, -(-k // block_size) + 1)
    chosen = _lowest(np.arange(count, dtype=np.uint32), take, seed)
    order = chosen[np.argsort(_ranks(chosen, seed), kind="stable")]
    pieces = []
    left = k
    for block in order.tolist():
        if left == 0:
            break
        start = block * block_size
        used = min(block_size, n_obs - start, left)
        pieces.append((start, used))
        left -= used
    pieces.sort()
    return np.concatenate([np.arange(start, start + used, dtype=np.int64) for start, used in pieces])


def _allocate(capacity: list[int], total: int, allocation: str) -> list[int]:
    if total == 0 or not capacity:
        return [0] * len(capacity)
    whole = sum(capacity)
    if total > whole:
        raise ValueError(f"cannot allocate {total} {plural(total, 'cell')} over {whole}")
    if allocation == "proportional":
        out = [total * c // whole for c in capacity]
        remainder = [total * c % whole for c in capacity]
        for s in sorted(range(len(capacity)), key=lambda s: (-remainder[s], s))[: total - sum(out)]:
            out[s] += 1
        return out
    level, used, remaining = 0, 0, len(capacity)
    for c in sorted(capacity):
        step = c - level
        if used + step * remaining > total:
            break
        used += step * remaining
        level = c
        remaining -= 1
    if remaining > 0:
        level += (total - used) // remaining
    out = [min(c, level) for c in capacity]
    left = total - sum(out)
    for s, c in enumerate(capacity):
        if left == 0:
            break
        if c > level:
            out[s] += 1
            left -= 1
    return out


def _split_over_strata(strata: NDArray[np.int64], n_strata: int, k: int, seed: int, allocation: str, minimum: int) -> NDArray[np.int64]:
    sizes = np.bincount(strata, minlength=n_strata)
    occupied = [s for s in range(n_strata) if sizes[s] > 0]
    base = [min(minimum, int(sizes[s])) for s in occupied]
    capacity = [int(sizes[s]) - b for s, b in zip(occupied, base, strict=False)]
    if sum(base) > k:
        raise ValueError(
            f"a subsample of {k} {plural(k, 'cell')} cannot hold {minimum} {plural(minimum, 'cell')} of each of "
            f"{len(occupied)} {plural(len(occupied), 'category', 'categories')} ({sum(base)} {plural(sum(base), 'cell')})"
        )
    extra = iter(_allocate([c for c in capacity if c > 0], k - sum(base), allocation))
    order = np.argsort(strata, kind="stable")
    bounds = np.concatenate([[0], np.cumsum(sizes)])
    kept = []
    for s, b, c in zip(occupied, base, capacity, strict=False):
        quota = b + (next(extra) if c > 0 else 0)
        if quota:
            kept.append(_lowest(order[bounds[s]:bounds[s + 1]].astype(np.uint32), quota, seed))
    return np.sort(np.concatenate(kept)) if kept else np.zeros(0, dtype=np.int64)


def _stratified(codes: NDArray[Any], n_categories: int, k: int, seed: int, allocation: str, minimum: int) -> NDArray[np.int64]:
    values = np.asarray(codes, dtype=np.float64)
    missing = np.isnan(values) | (values < 0)
    if np.any(values[~missing] >= n_categories) or np.any(values[~missing] != np.floor(values[~missing])):
        raise ValueError(f"category codes must be integers in [0, {n_categories}) or missing")
    strata = np.where(missing, n_categories, np.nan_to_num(values, nan=n_categories)).astype(np.int64)
    return _split_over_strata(strata, n_categories + 1, k, seed, allocation, minimum)


def grid_boxes(coordinates: NDArray[Any], grid_size: int) -> NDArray[np.int64]:
    """The grid box of each cell of an embedding (``spatial``), per the specification.

    Coordinates are taken as float32, then each axis' box is
    ``floor((x - min) * grid_size / (max - min))`` in binary64 (``grid_size - 1``
    at the maximum, 0 on an axis of one value); a cell with a non-finite
    coordinate is in box ``grid_size ** dimensions``.
    """
    values = np.asarray(coordinates, dtype=np.float32)
    if values.ndim == 1:
        values = values[:, None]
    dimensions = values.shape[1]
    if dimensions not in (1, 2, 3):
        raise ValueError("a spatial subsample needs a 1-, 2- or 3-dimensional embedding")
    wide = values.astype(np.float64)
    finite = np.all(np.isfinite(wide), axis=1)
    boxes = np.zeros(values.shape[0], dtype=np.int64)
    scale = 1
    for d in range(dimensions):
        column = wide[finite, d]
        low = column.min() if column.size else 0.0
        high = column.max() if column.size else 0.0
        width = high - low
        b: NDArray[np.int64] = np.zeros(values.shape[0], dtype=np.int64)
        if width > 0:
            with np.errstate(invalid="ignore"):
                scaled = (wide[:, d] - low) * float(grid_size)
                b = np.where(finite, np.minimum(np.floor(scaled / width), grid_size - 1), 0).astype(np.int64)
        boxes += b * scale
        scale *= grid_size
    boxes[~finite] = grid_size**dimensions
    return boxes


def _spatial(coordinates: NDArray[Any], k: int, seed: int, grid_size: int) -> NDArray[np.int64]:
    boxes = grid_boxes(coordinates, grid_size)
    occupied, strata = np.unique(boxes, return_inverse=True)
    return _split_over_strata(strata.astype(np.int64), occupied.size, k, seed, "balanced", 0)


def _n_obs_of(data: Any) -> int:
    if isinstance(data, int | np.integer) and not isinstance(data, bool):
        return _require_uint32(data, "n_obs")
    n_obs = getattr(data, "n_obs", None)
    if n_obs is None:
        raise TypeError("data must be an AnnData, a MuData or a cell count")
    return _require_uint32(n_obs, "n_obs")


def _codes_of(data: Any, groupby: Any) -> tuple[NDArray[Any], int]:
    if isinstance(groupby, str):
        if not hasattr(data, "obs"):
            raise TypeError("groupby names an obs column, so data must be an AnnData or a MuData")
        if groupby not in data.obs.columns:
            raise KeyError(f"obs has no column {groupby!r}")
        column = data.obs[groupby]
        if str(column.dtype) == "category":
            return np.asarray(column.cat.codes, dtype=np.int64), len(column.cat.categories)
        # Plain values (text): categories in code-point order, as below.
        groupby = column.to_numpy()
    import pandas as pd

    categorical = pd.Categorical(groupby)
    return np.asarray(categorical.codes, dtype=np.int64), len(categorical.categories)


def _coordinates_of(data: Any, basis: Any, columns: Sequence[int]) -> NDArray[Any]:
    if isinstance(basis, str):
        if not hasattr(data, "obsm"):
            raise TypeError("basis names an obsm key, so data must be an AnnData or a MuData")
        if basis not in data.obsm:
            raise KeyError(f"obsm has no key {basis!r}")
        values = np.asarray(data.obsm[basis])
    else:
        values = np.asarray(basis)
    if values.ndim == 1:
        values = values[:, None]
    return values[:, list(columns)]


def subsample_cells(
    data: Any,
    n_cells: int | None = None,
    *,
    fraction: float | None = None,
    method: Literal["uniform", "blocks", "stratified", "spatial"] = "uniform",
    seed: int = 0,
    groupby: Any = None,
    allocation: Literal["proportional", "balanced"] = "proportional",
    min_per_group: int = 0,
    block_size: int | None = None,
    basis: Any = None,
    columns: Sequence[int] = (0, 1),
    grid_size: int = 32,
) -> NDArray[np.int64]:
    """The cells a Cellucid subsample keeps, as original positions.

    The same request opens the same cells in the web app (Open data, *Cells to
    open*) and in R (``cellucid_subsample_cells()``): every implementation
    follows one specification and is checked against one conformance fixture.

    Parameters
    ----------
    data
        An AnnData, a MuData, or the cell count ``n_obs`` itself (enough for
        ``uniform`` and ``blocks``).
    n_cells
        How many cells to keep; or give ``fraction``. A count of every cell or
        more keeps every cell.
    fraction
        The share of cells to keep, in ``(0, 1]``; the count is
        ``n_obs × fraction`` rounded half up, at least 1, computed exactly
        from the fraction's shortest decimal text (``repr``), never as a
        binary64 product, as the web app and R compute it. Any count from 1 is a
        selection here; the web app opens a subsample of at least 2 cells,
        since one cell has no coordinate variation to draw.
    method
        ``"uniform"`` (every cell alike), ``"blocks"`` (runs of ``block_size``
        neighbouring cells in file order, kept together), ``"stratified"``
        (split over the categories of ``groupby``), or ``"spatial"`` (split
        over the occupied boxes of a ``grid_size`` grid on ``basis``, balanced,
        so sparse regions keep their cells).
    seed
        A uint32. The same seed keeps the same cells; with ``uniform`` and
        ``blocks`` a smaller subsample is inside a larger one.
    groupby
        ``stratified``: an obs column name, or one value per cell. A
        categorical keeps its categories' order; plain values (text) become
        categories in code-point order, as in the web app and R. Cells
        without a category are one more stratum.
    allocation
        ``stratified``: ``"proportional"`` (largest remainder) or
        ``"balanced"`` (equal shares, water-filled).
    min_per_group
        ``stratified``: each category first gets this many cells (or all it has).
    block_size
        ``blocks``: cells per block. The web app derives it from the file's
        chunk layout and shows it; give that value to reproduce its cells.
    basis
        ``spatial``: an obsm key, or an array of coordinates (one row per cell).
    columns
        ``spatial``: the embedding columns used, one to three.
    grid_size
        ``spatial``: boxes per axis.

    Returns
    -------
    numpy.ndarray
        The kept cells' positions (rows of ``obs``), ascending, int64.

    Examples
    --------
    >>> kept = cellucid.subsample_cells(adata, fraction=0.1, seed=0)
    >>> small = adata[kept].copy()
    """
    n_obs = _n_obs_of(data)
    seed = _require_uint32(seed, "seed")
    if method not in SUBSAMPLE_METHODS:
        raise ValueError(f"method must be one of {', '.join(SUBSAMPLE_METHODS)}; got {method!r}")
    if (n_cells is None) == (fraction is None):
        raise TypeError("give exactly one of n_cells and fraction")
    k = count_from_fraction(fraction, n_obs) if fraction is not None else _require_uint32(n_cells, "n_cells", minimum=1)
    if k >= n_obs:
        return np.arange(n_obs, dtype=np.int64)
    if method == "uniform":
        return _uniform(n_obs, k, seed)
    if method == "blocks":
        if block_size is None:
            raise TypeError("blocks needs block_size (the web app shows the one it used)")
        return _blocks(n_obs, k, seed, _require_uint32(block_size, "block_size", minimum=1))
    if method == "stratified":
        if groupby is None:
            raise TypeError("stratified needs groupby")
        if allocation not in ("proportional", "balanced"):
            raise ValueError(f"allocation must be 'proportional' or 'balanced'; got {allocation!r}")
        codes, n_categories = _codes_of(data, groupby)
        if codes.shape[0] != n_obs:
            raise ValueError(
                f"groupby has {codes.shape[0]} {plural(codes.shape[0], 'value')} for {n_obs} {plural(n_obs, 'cell')}"
            )
        return _stratified(codes, n_categories, k, seed, allocation, _require_uint32(min_per_group, "min_per_group"))
    if basis is None:
        raise TypeError("spatial needs basis")
    if not 1 <= len(columns) <= 3 or any(isinstance(c, bool) or not isinstance(c, (int, np.integer)) or c < 0 for c in columns):
        raise ValueError("spatial uses one to three embedding columns, each a non-negative column index")
    grid_size = _require_uint32(grid_size, "grid_size", minimum=1)
    if grid_size ** len(columns) >= 2**53:
        raise ValueError(f"grid_size**{len(columns)} boxes are more than are numbered exactly (2**53)")
    coordinates = _coordinates_of(data, basis, columns)
    if coordinates.shape[0] != n_obs:
        raise ValueError(
            f"basis has {coordinates.shape[0]} {plural(coordinates.shape[0], 'row')} for {n_obs} {plural(n_obs, 'cell')}"
        )
    return _spatial(coordinates, k, seed, _require_uint32(grid_size, "grid_size", minimum=1))


def cell_set_digest(n_obs: int, kept: Sequence[int] | NDArray[Any] | None) -> str:
    """The cell-set digest of the specification: 16 lowercase hex digits.

    ``kept`` is the ascending original positions, or ``None`` for every cell. A
    session saved on a subsample records it (manifest ``cells.digest``).
    """
    n_obs = _require_uint32(n_obs, "n_obs")
    positions = np.arange(n_obs, dtype=np.uint64) if kept is None else np.asarray(kept, dtype=np.uint64)
    k = int(positions.size)
    a = (0x811C9DC5 ^ n_obs) & _M32
    b = (0x9E3779B9 ^ k) & _M32
    for p in positions.tolist():
        a = ((a ^ p) * 0x01000193) & _M32
        x = (b ^ p) & _M32
        b = ((((x << 13) | (x >> 19)) & _M32) * 5 + 0xE6546B64) & _M32
    return f"{_fmix32_int(a ^ k):08x}{_fmix32_int(b ^ n_obs):08x}"


def encode_cell_set(kept: Sequence[int] | NDArray[Any]) -> bytes:
    """A cell set as a session stores it: each gap between ascending positions minus one, LEB128."""
    out = bytearray()
    previous = -1
    for p in np.asarray(kept, dtype=np.int64).tolist():
        if p <= previous:
            raise ValueError("a cell set lists positions strictly ascending")
        gap = p - previous - 1
        previous = p
        while gap >= 0x80:
            out.append((gap & 0x7F) | 0x80)
            gap >>= 7
        out.append(gap)
    return bytes(out)


def decode_cell_set(data: bytes, n_obs: int) -> NDArray[np.int64] | None:
    """The positions of an encoded cell set (``None``: no bytes, every cell)."""
    n_obs = _require_uint32(n_obs, "n_obs")
    if len(data) == 0:
        return None
    positions = []
    previous = -1
    gap = shift = 0
    for byte in data:
        gap |= (byte & 0x7F) << shift
        if byte & 0x80:
            shift += 7
            if shift > 35:
                raise ValueError("a stored cell set is malformed")
            continue
        previous += gap + 1
        if previous >= n_obs:
            raise ValueError(f"a stored cell set names a cell beyond the {n_obs} {plural(n_obs, 'cell')} of the dataset")
        positions.append(previous)
        gap = shift = 0
    if shift:
        raise ValueError("a stored cell set is truncated")
    return np.asarray(positions, dtype=np.int64)
