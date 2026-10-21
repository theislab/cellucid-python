"""Write the latent-outlier conformance fixture shared by cellucid-python, cellucid-r and the web app.

``tests/fixtures/latent_outlier_conformance.json`` holds latent spaces (float32
values, written as the exact numbers they are), category codes, and the
outlier-quantile bytes ``prepare()`` writes for them. cellucid-r
(``tests/testthat/fixtures/latent_outlier_conformance.json``) and the web app
(``tests/fixtures/latent-outlier-conformance.json``) keep the same bytes and
the same SHA-256 pin, and each checks its own implementation against them.

Run from the repository root: ``python scripts/generate_latent_outlier_conformance.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from cellucid.prepare_data._centroids import _compute_latent_space_quantiles

TARGET = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "latent_outlier_conformance.json"


def _case(name: str, latent: np.ndarray, codes: np.ndarray, n_categories: int, min_points: int) -> dict:
    latent = np.asarray(latent, dtype=np.float32)
    quantiles = _compute_latent_space_quantiles(
        latent, np.asarray(codes), [str(index) for index in range(n_categories)], min_points=min_points
    )
    return {
        "name": name,
        "width": int(latent.shape[1]),
        # Row-major, one float32 value each, as the exact double it widens to.
        "latent": [float(value) for value in latent.reshape(-1)],
        "codes": [int(code) for code in codes],
        "nCategories": n_categories,
        "minPoints": min_points,
        "quantilesHex": quantiles.astype("<f4").tobytes().hex(),
    }


def main() -> None:
    rng = np.random.default_rng(20261019)
    # Pairs mirrored about a far-off centre: their distances tie exactly in
    # real arithmetic, so rounding -- the order of summation -- decides the
    # ranks. A float32 or pairwise summation gives other bytes here.
    centre = np.random.default_rng(0).normal(0.0, 1.0, size=(1, 12)) * 1e4
    offsets = np.random.default_rng(0).integers(-50, 50, size=(10, 12)).astype(np.float64) * 0.37
    mirrored = np.concatenate([centre + offsets, centre - offsets])
    cases = [
        _case(
            "mirrored pairs whose ranks only the stated summation order decides",
            mirrored,
            np.zeros(len(mirrored), dtype=int),
            1,
            2,
        ),
        _case(
            "ties share the larger rank; a small category and a missing code are NaN",
            np.array([[0, 0], [1, 0], [1, 0], [4, 0], [7, 7], [8, 8], [3, 3]]),
            np.array([0, 0, 0, 0, 1, 1, -1]),
            2,
            3,
        ),
        _case(
            "a wide space at a large scale, where the order of summation shows",
            rng.normal(0.0, 3.0, size=(300, 9)) * 1e3 + 1e5,
            rng.integers(-1, 5, size=300),
            5,
            10,
        ),
        _case(
            "one column, many exact ties",
            rng.integers(0, 4, size=(120, 1)).astype(np.float64),
            np.arange(120) % 3,
            3,
            2,
        ),
    ]
    document = {
        "_about": (
            "Latent-space outlier quantiles: float32 values widened to float64, centroid sums "
            "in cell order, squared distances in column order, the share at that distance or "
            "closer as float32. Written by cellucid-python's "
            "scripts/generate_latent_outlier_conformance.py; kept byte for byte in cellucid-r "
            "and the web app."
        ),
        "version": 1,
        "cases": cases,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
