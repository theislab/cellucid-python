"""Why an array of numbers cannot be stored, counted rather than named.

``x must contain only finite values`` is true and stops there. One NaN in
eighteen million cells and every value NaN produce that same sentence, and
they want opposite responses: the first is a cell to look at, the second is
the wrong matrix. One pass over the array therefore reports how many values
are NaN, how many are each infinity, how many overflow or underflow float32,
and the first positions affected. cellucid-r's ``.describe_non_finite()``
reports the same counts, so an input refused by either writer is refused with
the same information.

The counting is not on the happy path: callers keep their fast all-finite
check and build the description only when it fails.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy import sparse

from ._plural import plural

#: How many offending positions to name. Enough to check, few enough to read.
NON_FINITE_EXAMPLE_LIMIT = 5


@dataclass(frozen=True)
class _NonFiniteCounts:
    """What one pass over an array found, before anyone words it."""

    total: int
    nan: int
    positive_infinity: int
    negative_infinity: int
    float32_overflow: int
    float32_underflow: int
    examples: tuple[int, ...]

    @property
    def offending(self) -> int:
        return (
            self.nan
            + self.positive_infinity
            + self.negative_infinity
            + self.float32_overflow
            + self.float32_underflow
        )

    def counted(self) -> str:
        parts = [
            (self.nan, "NaN"),
            (self.positive_infinity, "+Inf"),
            (self.negative_infinity, "-Inf"),
            (self.float32_overflow, "beyond the float32 range"),
            (self.float32_underflow, "below the smallest float32"),
        ]
        return ", ".join(f"{count:,} {label}" for count, label in parts if count)


def _count_non_finite(
    values: object,
) -> tuple[_NonFiniteCounts, tuple[int, ...]] | None:
    """One pass over any array, or ``None`` when nothing is wrong with it."""
    dense = (
        cast("sparse.spmatrix", values).toarray()
        if sparse.issparse(values)
        else values
    )
    source = np.asarray(dense)
    if source.dtype.kind not in {"i", "u", "f"}:
        # Not a numbers problem. The caller's own type error says it better.
        return None

    with np.errstate(over="ignore", invalid="ignore"):
        narrowed = source.astype(np.float32, copy=False)

    finite = np.isfinite(source)
    overflow = finite & ~np.isfinite(narrowed)
    underflow = (source != 0) & (narrowed == 0)
    offending = ~finite | overflow | underflow
    if not offending.any():
        return None

    flat = np.flatnonzero(offending.reshape(-1))
    counts = _NonFiniteCounts(
        total=int(source.size),
        nan=int(np.isnan(source).sum()),
        positive_infinity=int((source == np.inf).sum()),
        negative_infinity=int((source == -np.inf).sum()),
        float32_overflow=int(overflow.sum()),
        float32_underflow=int(underflow.sum()),
        examples=tuple(int(index) for index in flat[:NON_FINITE_EXAMPLE_LIMIT]),
    )
    return counts, source.shape


def _format_positions(examples: tuple[int, ...], shape: tuple[int, ...]) -> str:
    """The offending positions as the caller indexes them.

    A flat index for a one-dimensional column, ``[row, column]`` for anything
    with axes -- the position that can be pasted back into the failing call.
    cellucid-r prints the same list one-based and column-major, because that is
    how its caller indexes; the report is shared, the convention inside it
    belongs to the reader's language.
    """
    if len(shape) <= 1:
        return ", ".join(f"{index:,}" for index in examples)
    unraveled = np.unravel_index(np.asarray(examples, dtype=np.int64), shape)
    return ", ".join(
        "[" + ", ".join(f"{int(axis[position]):,}" for axis in unraveled) + "]"
        for position in range(len(examples))
    )


def describe_non_finite(
    values: object,
    subject: str,
    *,
    positions: bool = True,
) -> str:
    """Why one array cannot be stored as finite float32, counted rather than named.

    ``positions`` is ``False`` for a sparse matrix, whose stored non-zero
    entries have no index the caller would recognise; the counts are still
    exactly what was found. cellucid-r's ``.describe_non_finite()`` takes the
    same argument for the same reason.

    Returns ``""`` when there is nothing to describe, so a caller can append it
    unconditionally to a message whose own check has already failed.
    """
    counted = _count_non_finite(values)
    if counted is None:
        return ""
    counts, shape = counted
    scope = "values" if positions else "stored non-zero values"
    where = ""
    if positions:
        more = ", ..." if counts.offending > len(counts.examples) else ""
        where = (
            f" First affected {plural(counts.offending, 'position')}: "
            f"{_format_positions(counts.examples, shape)}{more}."
        )
    return (
        f"{subject}: of {counts.total:,} {scope}, {counts.counted()}.{where} "
        "Cellucid stores finite float32, because a colour scale has no "
        "position for an infinity and no value for a magnitude float32 "
        "cannot hold."
    )
