#!/usr/bin/env python3
"""Write the shared transition-drift conformance fixture.

    python scripts/generate_transition_drift_conformance.py [OUTPUT]

The fixture (``tests/fixtures/transition_drift_conformance.json``, and a
byte-identical copy in cellucid-r's ``tests/testthat/fixtures/``) pins the
transition drift both packages compute, ``cellucid.compute_transition_drift``
and cellucid-r's ``cellucid_transition_drift()``: ``T @ E - E``, with the
product divided by the row sums of ``T`` first when rows are normalized.

The expected drift is computed here with exact rational arithmetic
(:class:`fractions.Fraction`), independently of numpy, scipy, R and every BLAS,
from inputs drawn with :class:`random.Random`, whose ``random()`` stream is the
same on every Python version. Running this script again writes the same bytes;
cellucid-python's ``tests/test_transition_drift_conformance.py`` checks that.

Two kinds of case:

``exact``
    Every transition value is a multiple of 2**-4 and every coordinate a
    multiple of 2**-8, small enough that every product and every partial sum
    of every row is exactly representable in float64, whatever order a BLAS,
    a sparse kernel or a fused multiply-add sums them in. Each implementation
    then computes the product and the row sums exactly, and the two remaining
    operations -- one division, one subtraction -- are single IEEE 754
    operations, correctly rounded everywhere. The expected drift is those two
    operations applied here in float64, and ``bound`` is 0: an implementation
    must reproduce it bit for bit.

``general``
    Arbitrary float64 transition values (CellRank-like rows) and float32
    coordinates. Summation order is not part of the contract -- numpy's
    pairwise sums, a BLAS's blocked dot products and R's long-double row sums
    all differ in the last bits -- so the expected drift is the exact result
    rounded once to float64, and ``bound`` is, per value, twice the standard
    forward-error bound of the computation in any summation order (Higham,
    *Accuracy and Stability of Numerical Algorithms*, 2nd ed., section 3.1):
    with ``u = 2**-53``, ``n`` the number of cells, ``A = |T| @ |E|`` and
    ``s`` the row sum (1 without normalization), the computed drift is within
    ``(2n + 3) u (A / s + |E|)`` of the exact one to first order, the expected
    value's own rounding included; ``bound`` doubles that to cover the
    second-order terms.
"""

from __future__ import annotations

import json
import random
import struct
import sys
from fractions import Fraction
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPOSITORY_ROOT / "tests" / "fixtures" / "transition_drift_conformance.json"

UNIT_ROUNDOFF = Fraction(1, 2**53)
EXACT_TRANSITION_BITS = 4
EXACT_COORDINATE_BITS = 8

ABOUT = (
    "Transition drift T @ E - E (the product divided by the row sums of T first "
    "when normalizeRows), for cellucid-python's compute_transition_drift and "
    "cellucid-r's cellucid_transition_drift(). Written by cellucid-python's "
    "scripts/generate_transition_drift_conformance.py with exact rational "
    "arithmetic, independently of both. Matrices are row-major; a sparse "
    "transition matrix is 0-based coordinates in row-major order. An exact case "
    "is computed without rounding up to the division and the subtraction, so "
    "bound is 0 and the drift is matched bit for bit; a general case gives the "
    "exact drift rounded to float64 and, per value, a bound twice the forward "
    "error of the computation in any summation order: 2 (2n + 3) u (A / s + |E|) "
    "with u = 2**-53, n cells, A = |T| @ |E| and s the row sum (1 without "
    "normalization). An implementation passes when |drift - expected| <= bound "
    "for every value, with the matrix given dense and given sparse."
)


def _float32(value: float) -> float:
    return float(struct.unpack("<f", struct.pack("<f", value))[0])


class _Case:
    def __init__(
        self,
        name: str,
        *,
        exact: bool,
        normalize_rows: bool,
        sparse: bool,
        rows: list[list[tuple[int, float]]],
        embedding: list[list[float]],
    ) -> None:
        self.name = name
        self.exact = exact
        self.normalize_rows = normalize_rows
        self.sparse = sparse
        self.rows = rows
        self.embedding = embedding

    def document(self) -> dict[str, object]:
        n_cells = len(self.embedding)
        width = len(self.embedding[0])
        drift: list[float] = []
        bound: list[float] = []
        for i, row in enumerate(self.rows):
            row_sum = sum((Fraction(value) for _, value in row), Fraction(0))
            if self.normalize_rows and row_sum == 0:
                raise ValueError(f"{self.name}: row {i} sums to 0")
            for k in range(width):
                coordinate = Fraction(self.embedding[i][k])
                product = sum(
                    (Fraction(value) * Fraction(self.embedding[j][k]) for j, value in row),
                    Fraction(0),
                )
                magnitude = sum(
                    (abs(Fraction(value) * Fraction(self.embedding[j][k])) for j, value in row),
                    Fraction(0),
                )
                if self.exact:
                    self._require_exact(magnitude, row_sum)
                    quotient = float(product)
                    if self.normalize_rows:
                        quotient = quotient / float(row_sum)
                    drift.append(quotient - self.embedding[i][k])
                    bound.append(0.0)
                    continue
                scale = row_sum if self.normalize_rows else Fraction(1)
                drift.append(float(product / scale - coordinate))
                bound.append(
                    float(
                        2
                        * (2 * n_cells + 3)
                        * UNIT_ROUNDOFF
                        * (magnitude / scale + abs(coordinate))
                    )
                )
        if self.sparse:
            entries = [(i, j, value) for i, row in enumerate(self.rows) for j, value in row if value != 0]
            transition: dict[str, object] = {
                "kind": "sparse",
                "rows": [i for i, _, _ in entries],
                "columns": [j for _, j, _ in entries],
                "values": [value for _, _, value in entries],
            }
        else:
            dense = [[0.0] * n_cells for _ in range(n_cells)]
            for i, row in enumerate(self.rows):
                for j, value in row:
                    dense[i][j] = value
            transition = {"kind": "dense", "values": [value for row in dense for value in row]}
        return {
            "name": self.name,
            "exact": self.exact,
            "normalizeRows": self.normalize_rows,
            "nCells": n_cells,
            "width": width,
            "transition": transition,
            "embedding": [value for row in self.embedding for value in row],
            "drift": drift,
            "bound": bound,
        }

    def _require_exact(self, magnitude: Fraction, row_sum: Fraction) -> None:
        # Every partial sum of the row is a multiple of the grid below and no
        # larger than the sum of magnitudes, so it is exact in float64 when
        # that sum spans fewer than 53 bits of the grid.
        grid = 2 ** (EXACT_TRANSITION_BITS + EXACT_COORDINATE_BITS)
        if magnitude * grid >= 2**53 or row_sum * 2**EXACT_TRANSITION_BITS >= 2**53:
            raise ValueError(f"{self.name}: inputs too wide to be summed exactly")


def _exact_transition(generator: random.Random) -> float:
    return int(generator.random() * 256) / 2**EXACT_TRANSITION_BITS


def _exact_coordinate(generator: random.Random) -> float:
    return (int(generator.random() * 2**15) - 2**14) / 2**EXACT_COORDINATE_BITS


def _distinct_columns(generator: random.Random, n_cells: int, count: int) -> list[int]:
    chosen: list[int] = []
    while len(chosen) < count:
        column = int(generator.random() * n_cells)
        if column not in chosen:
            chosen.append(column)
    return sorted(chosen)


def _exact_case(
    name: str,
    seed: int,
    *,
    n_cells: int,
    width: int,
    normalize_rows: bool,
    sparse: bool,
    zero_row: int | None = None,
) -> _Case:
    generator = random.Random(seed)
    rows: list[list[tuple[int, float]]] = []
    for i in range(n_cells):
        if i == zero_row:
            rows.append([])
            continue
        columns = _distinct_columns(generator, n_cells, max(1, n_cells // 2)) if sparse else list(range(n_cells))
        row = [(j, _exact_transition(generator)) for j in columns]
        if all(value == 0 for _, value in row):
            row[0] = (row[0][0], 1.0)
        rows.append([(j, value) for j, value in row if value != 0 or not sparse])
    embedding = [[_exact_coordinate(generator) for _ in range(width)] for _ in range(n_cells)]
    return _Case(name, exact=True, normalize_rows=normalize_rows, sparse=sparse, rows=rows, embedding=embedding)


def _general_case(
    name: str,
    seed: int,
    *,
    n_cells: int,
    width: int,
    normalize_rows: bool,
    sparse: bool,
    neighbours: int,
    scale: float,
) -> _Case:
    generator = random.Random(seed)
    rows: list[list[tuple[int, float]]] = []
    for _ in range(n_cells):
        columns = _distinct_columns(generator, n_cells, neighbours) if sparse else list(range(n_cells))
        values = [generator.random() if sparse or generator.random() < 0.7 else 0.0 for _ in columns]
        if all(value == 0 for value in values):
            values[0] = generator.random() + 0.5
        if not normalize_rows:
            # Already row-stochastic, as CellRank's transition matrix is: the
            # rows sum to 1 up to their own rounding.
            total = sum(values)
            values = [value / total for value in values]
        rows.append([(j, value) for j, value in zip(columns, values, strict=True) if value != 0 or not sparse])
    embedding = [
        [_float32((2 * generator.random() - 1) * scale) for _ in range(width)] for _ in range(n_cells)
    ]
    return _Case(name, exact=False, normalize_rows=normalize_rows, sparse=sparse, rows=rows, embedding=embedding)


def fixture_text() -> str:
    """The fixture, exactly as it is written to disk."""
    cases = [
        _exact_case("exact, dense, rows normalized", 1, n_cells=5, width=2, normalize_rows=True, sparse=False),
        _exact_case("exact, sparse, rows normalized", 2, n_cells=8, width=3, normalize_rows=True, sparse=True),
        _exact_case("exact, dense, one column, as given", 3, n_cells=4, width=1, normalize_rows=False, sparse=False),
        _exact_case(
            "exact, sparse, an empty row, as given",
            4,
            n_cells=6,
            width=2,
            normalize_rows=False,
            sparse=True,
            zero_row=2,
        ),
        _exact_case("exact, sparse, ten columns, rows normalized", 5, n_cells=7, width=10, normalize_rows=True, sparse=True),
        _general_case(
            "general, dense, rows normalized",
            6,
            n_cells=24,
            width=2,
            normalize_rows=True,
            sparse=False,
            neighbours=0,
            scale=10.0,
        ),
        _general_case(
            "general, sparse neighbours, rows normalized",
            7,
            n_cells=60,
            width=3,
            normalize_rows=True,
            sparse=True,
            neighbours=8,
            scale=10.0,
        ),
        _general_case(
            "general, sparse row-stochastic, as given",
            8,
            n_cells=50,
            width=2,
            normalize_rows=False,
            sparse=True,
            neighbours=6,
            scale=10.0,
        ),
        _general_case(
            "general, sparse, coordinates in the ten thousands",
            9,
            n_cells=40,
            width=2,
            normalize_rows=True,
            sparse=True,
            neighbours=5,
            scale=1.0e4,
        ),
    ]
    document = {"_about": ABOUT, "version": 1, "cases": [case.document() for case in cases]}
    return json.dumps(document, indent=1, ensure_ascii=True, allow_nan=False) + "\n"


def main(argv: list[str]) -> int:
    output = Path(argv[1]) if len(argv) > 1 else DEFAULT_OUTPUT
    output.write_text(fixture_text(), encoding="utf-8", newline="\n")
    print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
