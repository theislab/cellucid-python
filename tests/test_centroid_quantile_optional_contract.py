"""``centroid_outlier_quantile=None`` is a documented input, so it must type-check.

``prepare()`` accepts ``None`` and then publishes every categorical field with
empty centroid lists, which
``docs/user_guide/python_package/c_data_preparation_api/04_obs_cell_metadata.md``
states and the R peer matches with ``NULL``. The package ships ``py.typed``, so
the annotation is a promise enforced inside *someone else's* project: while the
parameter was annotated ``float``, the documented call was an ``[arg-type]``
error for every consumer that type-checks. These tests pin both halves — the
runtime behaviour and the accepted type set — so the annotation cannot silently
narrow back to ``float`` or widen to something that accepts anything at all.
"""

from __future__ import annotations

import inspect
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare
from cellucid.prepare_data._generation import _prepare_generation

CONSUMER_MODULE = "uses_prepare.py"

# ``_categoricalFields`` entries are positional: the sixth element is the
# ``centroidsByDim`` mapping, keyed by embedding dimension.
CENTROIDS_BY_DIM = 5


def _export(out_dir: Path, quantile: float | None) -> dict:
    """Publish one two-cell-per-category export and return its obs manifest."""
    prepare(
        latent_space=np.array(
            [[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0]],
            dtype=np.float32,
        ),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["a", "a", "b", "b"])}),
        X_umap_2d=np.array(
            [[0.0, 0.0], [0.1, 0.0], [5.0, 5.0], [5.1, 5.0]],
            dtype=np.float32,
        ),
        out_dir=out_dir,
        obs_categorical_dtype="uint8",
        dataset_id="centroid-quantile-contract",
        dataset_name="Centroid quantile contract",
        centroid_outlier_quantile=quantile,
        centroid_min_points=1,
        force=True,
    )
    return json.loads((out_dir / "obs_manifest.json").read_text(encoding="utf-8"))


def test_none_publishes_every_categorical_field_with_empty_centroids(tmp_path: Path) -> None:
    manifest = _export(tmp_path / "exports", None)

    assert manifest["centroid_outlier_quantile"] is None

    fields = manifest["_categoricalFields"]
    assert fields, "the export must hold the categorical field this asserts about"
    for field in fields:
        centroids_by_dim = field[CENTROIDS_BY_DIM]
        assert centroids_by_dim, "the dimension keys survive; only the centroids are dropped"
        assert all(centroids == [] for centroids in centroids_by_dim.values()), centroids_by_dim


def test_a_quantile_still_publishes_centroids(tmp_path: Path) -> None:
    """The control run: empty lists are what ``None`` buys, not what every run does."""
    manifest = _export(tmp_path / "exports", 0.95)

    assert manifest["centroid_outlier_quantile"] == 0.95

    fields = manifest["_categoricalFields"]
    assert fields
    for field in fields:
        centroids_by_dim = field[CENTROIDS_BY_DIM]
        assert centroids_by_dim
        assert all(centroids for centroids in centroids_by_dim.values()), centroids_by_dim


def test_both_forwarding_signatures_admit_none() -> None:
    """``prepare`` and the internal generation step must agree on the type."""
    for function in (prepare, _prepare_generation):
        annotation = inspect.signature(function).parameters["centroid_outlier_quantile"].annotation
        assert annotation == float | None, f"{function.__name__}: {annotation}"


def test_a_consumer_type_checks_the_documented_none_call(tmp_path: Path) -> None:
    """Checked from outside the package, exactly as a dependent project would."""
    pytest.importorskip("mypy")
    consumer_dir = tmp_path / "consumer"
    consumer_dir.mkdir()
    (consumer_dir / CONSUMER_MODULE).write_text(
        "import cellucid\n"
        "\n"
        "cellucid.prepare(\n"
        '    obs_categorical_dtype="uint8",\n'
        '    dataset_name="Demo",\n'
        '    dataset_id="demo",\n'
        "    centroid_outlier_quantile=None,\n"
        ")\n"
        "cellucid.prepare(\n"
        '    obs_categorical_dtype="uint8",\n'
        '    dataset_name="Demo",\n'
        '    dataset_id="demo",\n'
        '    centroid_outlier_quantile="0.95",\n'
        ")\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--cache-dir",
            str(tmp_path / "mypy_cache"),
            CONSUMER_MODULE,
        ],
        capture_output=True,
        text=True,
        cwd=consumer_dir,
        check=False,
    )
    output = completed.stdout + completed.stderr

    assert 'has incompatible type "None"' not in output, output
    # The parameter admits None, not anything: a wrong type stays a consumer error.
    assert (
        'Argument "centroid_outlier_quantile" to "prepare" has incompatible type "str"' in output
    ), output
    assert "[arg-type]" in output, output
    assert "Found 1 error" in output, output
