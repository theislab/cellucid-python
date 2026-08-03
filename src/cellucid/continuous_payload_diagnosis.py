"""Why one continuous payload cannot be published, counted rather than named.

Cellucid publishes continuous values — gene expression and continuous ``obs``
columns — as finite float32. A colour scale needs a finite minimum and maximum,
and an infinity has no position on one: it makes the whole field's range
infinite, so every other cell collapses onto a single colour. Both the exporter
and the server therefore refuse a column that is not entirely finite.

Refusing is the easy half. The hard half is saying *why* in a way the person
looking at the screen can act on, and that is what this module exists for. The
refusal used to read

    Gene 'A2ML1' expression must contain only finite values.

which is true, and which the server then discarded entirely: ``_handle_var``
caught it alongside every other exception and answered ``500 Internal server
error``, so the browser could not name the gene, the reason, or anything to do
about it. Measured against a real 18,142,044-cell object, every one of 5,121
genes answered 500 after fifty seconds with a twenty-one byte body.

What replaces it counts. One pass over the column reports how many values are
NaN, how many are each infinity, how many overflow or underflow float32, and
which cells they are — because "one cell out of eighteen million" and "every
cell" are the same message today and want opposite responses.

The counting is not on the happy path. Callers keep their existing fast
all-finite check and build a diagnosis only when it fails, so a healthy column
pays nothing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy import sparse

__all__ = [
    "NON_FINITE_EXAMPLE_LIMIT",
    "NonFinitePayloadDiagnosis",
    "NonFinitePayloadError",
    "describe_non_finite",
    "diagnose_continuous_payload",
]

#: How many offending cell indices to name. Enough to check, few enough to read.
NON_FINITE_EXAMPLE_LIMIT = 5

_GENE_REMEDY = (
    "Cellucid publishes finite float32 only, because a colour scale has no "
    "position for an infinity. Repair the matrix before serving or exporting "
    "it -- `adata.X.data[~np.isfinite(adata.X.data)] = 0` for a sparse matrix, "
    "`np.nan_to_num(adata.X, copy=False)` for a dense one."
)


def _field_remedy(name: str) -> str:
    return (
        "Cellucid publishes finite float32 only, because a colour scale has no "
        "position for an infinity. Repair the column before serving or "
        f"exporting it -- `adata.obs[{name!r}] = "
        f"adata.obs[{name!r}].replace([np.inf, -np.inf], np.nan).fillna(0)` -- "
        "then reload."
    )


@dataclass(frozen=True)
class NonFinitePayloadDiagnosis:
    """One counted description of a continuous payload that cannot be published.

    Attributes
    ----------
    kind
        ``"gene"`` or ``"field"``. The only thing that changes the remedy: a
        gene comes out of ``adata.X`` and a field out of ``adata.obs``.
    name
        The gene symbol or column key, as the reader sees it.
    total
        How many values were examined.
    nan, positive_infinity, negative_infinity
        Counts by kind, in the source values.
    float32_overflow
        Finite in the source, infinite once narrowed to float32.
    float32_underflow
        Nonzero in the source, exactly zero once narrowed to float32. The quiet
        one: it is finite either way, so an overflow-only check would publish a
        column whose every value had silently become zero.
    examples
        The first offending cell indices, ``NaN`` included -- every count above
        is a value that cannot be published.
    """

    kind: str
    name: str
    total: int
    nan: int
    positive_infinity: int
    negative_infinity: int
    float32_overflow: int
    float32_underflow: int
    examples: tuple[int, ...]

    @property
    def offending(self) -> int:
        """How many values cannot be published."""
        return (
            self.nan
            + self.positive_infinity
            + self.negative_infinity
            + self.float32_overflow
            + self.float32_underflow
        )

    def as_payload(self) -> dict:
        """The JSON body a server answers a refused payload request with."""
        return {
            "error": "non_finite_continuous_payload",
            "kind": self.kind,
            "name": self.name,
            "counts": {
                "total": self.total,
                "offending": self.offending,
                "nan": self.nan,
                "positive_infinity": self.positive_infinity,
                "negative_infinity": self.negative_infinity,
                "float32_overflow": self.float32_overflow,
                "float32_underflow": self.float32_underflow,
            },
            "examples": list(self.examples),
            "message": str(self),
        }

    def _counted(self) -> str:
        parts = [
            (self.nan, "NaN"),
            (self.positive_infinity, "+Inf"),
            (self.negative_infinity, "-Inf"),
            (self.float32_overflow, "beyond the float32 range"),
            (self.float32_underflow, "below the smallest float32"),
        ]
        return ", ".join(f"{count:,} {label}" for count, label in parts if count)

    def __str__(self) -> str:
        noun = "Gene" if self.kind == "gene" else "Continuous obs field"
        examples = ", ".join(f"{index:,}" for index in self.examples)
        more = ", ..." if self.offending > len(self.examples) else ""
        remedy = _GENE_REMEDY if self.kind == "gene" else _field_remedy(self.name)
        return (
            f"{noun} {self.name!r} cannot be published: of {self.total:,} "
            f"cells, {self._counted()}. "
            f"First affected cell{'' if len(self.examples) == 1 else 's'}: "
            f"{examples}{more}. {remedy}"
        )


class NonFinitePayloadError(ValueError):
    """A refusal that carries its own counted diagnosis.

    Subclasses ``ValueError`` so every existing ``except ValueError`` around the
    validators keeps working; a caller that wants the counts asks for
    ``.diagnosis``.
    """

    def __init__(self, diagnosis: NonFinitePayloadDiagnosis) -> None:
        super().__init__(str(diagnosis))
        self.diagnosis = diagnosis


def diagnose_continuous_payload(
    values: object,
    *,
    kind: str,
    name: str,
) -> NonFinitePayloadDiagnosis | None:
    """Count every value that cannot be published, or return ``None``.

    ``None`` means the column is publishable, which is how a caller
    distinguishes "the values are wrong" from every other reason a validator
    might have refused -- a non-numeric dtype, or the wrong shape. Those keep
    their own messages rather than being described as a finiteness problem.

    Parameters
    ----------
    values
        Any array-like, dense or sparse.
    kind
        ``"gene"`` or ``"field"``.
    name
        The gene symbol or column key.
    """
    if kind not in {"gene", "field"}:
        raise ValueError('Diagnosis kind must be exactly "gene" or "field".')
    if not isinstance(name, str) or not name:
        raise ValueError("Diagnosis name must be one non-empty string.")

    counted = _count_non_finite(values)
    if counted is None:
        return None
    counts, _shape = counted
    return NonFinitePayloadDiagnosis(
        kind=kind,
        name=name,
        total=counts.total,
        nan=counts.nan,
        positive_infinity=counts.positive_infinity,
        negative_infinity=counts.negative_infinity,
        float32_overflow=counts.float32_overflow,
        float32_underflow=counts.float32_underflow,
        examples=counts.examples,
    )


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
    """One pass over any array, or ``None`` when nothing is wrong with it.

    The server's diagnosis and the exporter's refusal are the same five
    questions asked of the same values, so they are asked once here. What
    differs is who is being told: the server names cells and an AnnData repair,
    the exporter names positions in the array it was handed. That difference is
    wording, and wording is the only thing the two callers do separately.
    """
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
    """Why one exported array cannot be published, counted rather than named.

    ``x must contain only finite values.`` is true and stops there. One NaN in
    eighteen million cells and every value NaN produce that same sentence, and
    they want opposite responses: the first is a cell to look at, the second is
    the wrong matrix. This is the exporter's half of what
    :func:`diagnose_continuous_payload` gives the server, and what
    ``.describe_non_finite()`` gives cellucid-r's caller, so a dataset refused
    by either writer is refused with the same information.

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
        plural = "" if counts.offending == 1 else "s"
        more = ", ..." if counts.offending > len(counts.examples) else ""
        where = (
            f" First affected position{plural}: "
            f"{_format_positions(counts.examples, shape)}{more}."
        )
    return (
        f"{subject} cannot be published: of {counts.total:,} {scope}, "
        f"{counts.counted()}.{where} "
        "Cellucid publishes finite float32 only, because a colour scale has "
        "no position for an infinity and no value for a magnitude float32 "
        "cannot hold."
    )
