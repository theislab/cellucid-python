"""The latent-space outlier rule ``prepare()``, ``cellucid serve`` and the web app share.

The web app's ``data/latent-outliers.js`` reproduces these numbers bit for bit
(its ``outlier-quantile-exact-contract`` test compares against bytes this
package writes), so the rule is stated exactly: float32 values widened to
float64, centroid sums accumulated in cell order, squared distances
accumulated in column order, and the share at that distance or closer as
float32.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from cellucid.prepare_data._centroids import _compute_latent_space_quantiles


def _reference(latent: np.ndarray, codes: np.ndarray, n_categories: int, min_points: int) -> np.ndarray:
    """The stated rule, one scalar at a time."""
    values = latent.astype(np.float32)
    out = np.full(len(codes), np.nan, dtype=np.float32)
    for code in range(n_categories):
        members = [cell for cell in range(len(codes)) if codes[cell] == code]
        n = len(members)
        if n < min_points:
            continue
        squared = [0.0] * n
        for column in range(values.shape[1]):
            total = 0.0
            for cell in members:
                total += float(values[cell, column])
            centre = total / n
            for k, cell in enumerate(members):
                difference = float(values[cell, column]) - centre
                squared[k] += difference * difference
        distances = [math.sqrt(value) for value in squared]
        ordered = sorted(distances)
        for k, cell in enumerate(members):
            rank = sum(1 for value in ordered if value <= distances[k])
            out[cell] = np.float32(rank) / np.float32(n)
    return out


def test_ties_share_the_larger_rank_and_small_or_missing_categories_are_nan():
    latent = np.array([[0, 0], [1, 0], [1, 0], [4, 0], [7, 7], [8, 8], [3, 3]], dtype=np.float32)
    codes = np.array([0, 0, 0, 0, 1, 1, -1])
    quantiles = _compute_latent_space_quantiles(latent, codes, ["a", "b"], min_points=3)
    assert quantiles.dtype == np.float32
    assert quantiles[:4].tolist() == [0.75, 0.5, 0.5, 1.0]
    assert np.isnan(quantiles[4:]).all()


def test_the_computation_is_the_stated_rule_bit_for_bit():
    rng = np.random.default_rng(7)
    latent = (rng.normal(0.0, 3.0, size=(400, 7)) * 1e3).astype(np.float32)
    codes = rng.integers(-1, 4, size=400)
    expected = _reference(latent, codes, 4, 10)
    quantiles = _compute_latent_space_quantiles(latent, codes, ["a", "b", "c", "d"], min_points=10)
    assert quantiles.tobytes() == expected.tobytes()


def test_a_float64_latent_is_read_as_the_float32_it_is_stored_as():
    rng = np.random.default_rng(11)
    latent = rng.normal(size=(60, 3))
    codes = np.arange(60) % 2
    from_float64 = _compute_latent_space_quantiles(latent, codes, ["a", "b"], min_points=5)
    from_float32 = _compute_latent_space_quantiles(latent.astype(np.float32), codes, ["a", "b"], min_points=5)
    assert from_float64.tobytes() == from_float32.tobytes()


#: The shared fixture, kept byte for byte in cellucid-r
#: (tests/testthat/fixtures/latent_outlier_conformance.json) and the web app
#: (tests/fixtures/latent-outlier-conformance.json) under the same pin.
#: Regenerate with scripts/generate_latent_outlier_conformance.py and copy it
#: to both.
CONFORMANCE = Path(__file__).parent / "fixtures" / "latent_outlier_conformance.json"
CONFORMANCE_SHA256 = "46ea30291fc2d1fc51fc2021837eb8914062f8080b6ad51995cacaa1875b9d80"


def test_the_shared_conformance_cases_give_the_pinned_bytes():
    assert hashlib.sha256(CONFORMANCE.read_bytes()).hexdigest() == CONFORMANCE_SHA256
    cases = json.loads(CONFORMANCE.read_text(encoding="utf-8"))["cases"]
    assert len(cases) == 4
    for case in cases:
        latent = np.array(case["latent"], dtype=np.float64).reshape(-1, case["width"])
        assert latent.astype(np.float32).astype(np.float64).tobytes() == latent.tobytes(), case["name"]
        quantiles = _compute_latent_space_quantiles(
            latent,
            np.array(case["codes"]),
            [str(index) for index in range(case["nCategories"])],
            min_points=case["minPoints"],
        )
        assert quantiles.astype("<f4").tobytes().hex() == case["quantilesHex"], case["name"]
