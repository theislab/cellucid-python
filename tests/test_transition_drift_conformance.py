"""The shared transition-drift conformance fixture, run against the Python helper.

``tests/fixtures/transition_drift_conformance.json`` is kept byte for byte in
cellucid-python and cellucid-r (``tests/testthat/fixtures/``), and both
repositories pin the same SHA-256. It is written by
``scripts/generate_transition_drift_conformance.py`` with exact rational
arithmetic, independently of numpy, scipy and R, so
:func:`cellucid.compute_transition_drift` and cellucid-r's
``cellucid_transition_drift()`` are each held to the same expected drift:

- an ``exact`` case bit for bit (``bound`` 0): its product and row sums are
  exact in float64 in any summation order, and the division and subtraction
  that remain are correctly rounded everywhere;
- a ``general`` case within ``bound``, twice the forward-error bound of the
  computation in any summation order, since numpy, a BLAS and R sum in
  different orders.

Every case is run with the transition matrix dense and sparse.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy import sparse

from cellucid import compute_transition_drift

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
FIXTURE = REPOSITORY_ROOT / "tests" / "fixtures" / "transition_drift_conformance.json"
GENERATOR = REPOSITORY_ROOT / "scripts" / "generate_transition_drift_conformance.py"

#: The same digest is pinned in cellucid-r
#: (tests/testthat/test-transition-drift.R). Update both together.
FIXTURE_SHA256 = "0c3a15648c21a6441f2e530d25a2837298dacd5113942f3ef280b048a1ea81e8"


def _cases() -> list[dict[str, Any]]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))["cases"]


def test_the_fixture_is_the_pinned_generator_output() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    spec = importlib.util.spec_from_file_location("generate_transition_drift_conformance", GENERATOR)
    assert spec is not None and spec.loader is not None
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    assert generator.fixture_text().encode("utf-8") == FIXTURE.read_bytes()
    document = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert document["version"] == 1
    cases = document["cases"]
    assert sum(case["exact"] for case in cases) == 5
    assert sum(not case["exact"] for case in cases) == 4


@pytest.mark.parametrize("case", _cases(), ids=lambda case: case["name"])
@pytest.mark.parametrize("representation", ["dense", "sparse"])
def test_the_drift_is_the_conformance_drift(case: dict[str, Any], representation: str) -> None:
    n_cells, width = case["nCells"], case["width"]
    transition = case["transition"]
    if transition["kind"] == "dense":
        matrix: Any = np.asarray(transition["values"], dtype=np.float64).reshape(n_cells, n_cells)
    else:
        matrix = sparse.csr_matrix(
            (transition["values"], (transition["rows"], transition["columns"])),
            shape=(n_cells, n_cells),
        ).toarray()
    if representation == "sparse":
        matrix = sparse.csr_matrix(matrix)
    embedding = np.asarray(case["embedding"], dtype=np.float64).reshape(n_cells, width)
    expected = np.asarray(case["drift"], dtype=np.float64).reshape(n_cells, width)
    bound = np.asarray(case["bound"], dtype=np.float64).reshape(n_cells, width)

    drift = compute_transition_drift(matrix, embedding, normalize_rows=case["normalizeRows"])

    assert drift.dtype == np.float64
    if case["exact"]:
        assert drift.tobytes() == expected.tobytes()
    else:
        assert np.all(np.abs(drift - expected) <= bound)
