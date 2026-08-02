"""A manifest is validated in memory and then serialized; nothing proved the two
agreed until the bytes were read back.

``json.dumps`` is not a total function into the JSON the browser accepts: a
non-finite float is written as the bare tokens ``NaN``/``Infinity``, which
``json.loads`` happily reads back and ``JSON.parse`` refuses outright, and a
dict key that is not a string is coerced to one. Both publish a file that says
something other than what was validated, and both fail only when the web app
opens the export. The writer therefore re-reads every manifest it writes,
under a parser as strict as the reader's, and requires it to parse back to
exactly the payload that was validated.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data._export import prepare
from cellucid.prepare_data._manifest import (
    _first_json_difference,
    _require_manifest_reads_back,
)


def _prepare(out_dir: Path, **overrides: object) -> None:
    n_cells = 40
    keywords: dict[str, object] = {
        "latent_space": np.arange(n_cells * 3).reshape(n_cells, 3) / n_cells,
        "obs": pd.DataFrame(
            {
                "grp": pd.Categorical(["A"] * 20 + ["B"] * 20),
                "val": np.arange(n_cells) / n_cells,
            }
        ),
        "X_umap_2d": np.column_stack(
            [np.arange(n_cells) / n_cells, np.arange(n_cells)[::-1] / n_cells]
        ),
        "out_dir": out_dir,
        "obs_categorical_dtype": "uint8",
        "dataset_name": "Read back",
        "dataset_id": "readback",
        "centroid_min_points": 2,
        "force": True,
        "created_at": "2020-01-01T00:00:00Z",
    }
    keywords.update(overrides)
    prepare(**keywords)  # type: ignore[arg-type]


@pytest.mark.parametrize("variant", ["no-columns", "empty-obs-keys"])
def test_export_without_observation_fields_declares_object_schemas(
    tmp_path: Path, variant: str
) -> None:
    """``_obsSchemas`` is a JSON object even when the export carries no obs field.

    The web reader calls ``requireRecord`` on it and refuses an array, so the
    empty schema map has to be ``{}``. cellucid-r wrote ``[]`` here and
    published a dataset the viewer then rejected; this pins the Python half of
    the same contract so the two writers cannot drift apart again.
    """
    if variant == "no-columns":
        overrides: dict[str, object] = {"obs": pd.DataFrame(index=range(40))}
    else:
        overrides = {
            "obs": pd.DataFrame({"grp": pd.Categorical(["A"] * 20 + ["B"] * 20)}),
            "obs_keys": [],
        }
    _prepare(tmp_path / "export", **overrides)

    text = (tmp_path / "export" / "obs_manifest.json").read_text(encoding="utf-8")
    assert '"_obsSchemas": {}' in text
    manifest = json.loads(text)
    assert manifest["_obsSchemas"] == {}
    assert isinstance(manifest["_obsSchemas"], dict)
    assert manifest["_continuousFields"] == []
    assert manifest["_categoricalFields"] == []


def test_read_back_refuses_non_standard_json_constants(tmp_path: Path) -> None:
    """``NaN`` and ``Infinity`` are Python JSON, not JSON.

    ``json.loads`` accepts both, so a read-back that used the default parser
    would agree with the writer and still hand the browser a file
    ``JSON.parse`` throws on. The check has to be at least as strict as the
    reader it stands in for.
    """
    path = tmp_path / "manifest.json"
    for value in (float("nan"), float("inf"), float("-inf")):
        payload = {"minValue": value}
        path.write_text(json.dumps(payload), encoding="utf-8")
        # The default Python parser reads the file back without complaint,
        # which is exactly why the writer cannot rely on it.
        assert json.loads(path.read_text(encoding="utf-8")) is not None
        with pytest.raises(RuntimeError, match="does not read back"):
            _require_manifest_reads_back(path, payload)

    path.write_text(json.dumps({"minValue": 0.5}), encoding="utf-8")
    _require_manifest_reads_back(path, {"minValue": 0.5})


_READ_BACK_PAYLOAD: dict[str, object] = {
    "schemas": {},
    "categories": ["only"],
    "n_points": 3,
}


def test_read_back_accepts_the_bytes_that_carry_the_payload(tmp_path: Path) -> None:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(_READ_BACK_PAYLOAD), encoding="utf-8")
    _require_manifest_reads_back(path, _READ_BACK_PAYLOAD)


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("object_became_array", '{"schemas": [], "categories": ["only"], "n_points": 3}'),
        ("array_became_scalar", '{"schemas": {}, "categories": "only", "n_points": 3}'),
        ("key_dropped", '{"schemas": {}, "categories": ["only"]}'),
        ("key_added", '{"schemas": {}, "categories": ["only"], "n_points": 3, "x": 1}'),
        ("key_reordered", '{"categories": ["only"], "schemas": {}, "n_points": 3}'),
        ("value_changed", '{"schemas": {}, "categories": ["only"], "n_points": 4}'),
        ("number_became_string", '{"schemas": {}, "categories": ["only"], "n_points": "3"}'),
        ("not_json", '{"schemas": {},'),
    ],
)
def test_read_back_refuses_a_file_that_says_something_else(
    tmp_path: Path, name: str, text: str
) -> None:
    """Proved against the bytes, not against the encoder.

    Given the payload the writer validated, a file carrying any other shape is
    refused. These are the losses ``json.dumps`` and a truncated write can
    produce between a validated payload and the file a reader opens.
    """
    path = tmp_path / f"{name}.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(RuntimeError, match="does not read back"):
        _require_manifest_reads_back(path, _READ_BACK_PAYLOAD)


def test_first_json_difference_names_the_node_that_disagrees() -> None:
    """The message points at one node, so a failure is actionable."""
    assert _first_json_difference({"a": 1}, {"a": 1}) is None
    assert _first_json_difference([1, 2], [1, 2]) is None
    assert "$.a" in str(_first_json_difference({"a": 1}, {"a": 2}))
    assert "$.a[1]" in str(_first_json_difference({"a": [1, 2]}, {"a": [1, 3]}))
    assert "object" in str(_first_json_difference({"a": {}}, {"a": []}))
    assert "array" in str(_first_json_difference({"a": ["x"]}, {"a": "x"}))
    # A bool is not a number: JSON tells them apart and so does this.
    assert _first_json_difference({"a": True}, {"a": 1}) is not None
    assert _first_json_difference({"a": 1}, {"a": True}) is not None
    # An int and a float of equal value are one JSON number.
    assert _first_json_difference({"a": 1}, {"a": 1.0}) is None


def test_read_back_accepts_every_manifest_the_writer_emits(tmp_path: Path) -> None:
    """The guard must not refuse ordinary exports: one case per optional part."""
    n_cells = 40
    expression = (np.arange(n_cells * 3).reshape(n_cells, 3) / n_cells).astype(
        np.float32
    )
    var = pd.DataFrame(index=["g1", "g2", "g3"])
    connectivities = np.zeros((n_cells, n_cells))
    connectivities[0, 1] = connectivities[1, 0] = 0.5
    connectivities[2, 3] = connectivities[3, 2] = 1.0

    cases: dict[str, dict[str, object]] = {
        "plain": {},
        "quantized": {"obs_continuous_quantization": 8, "var_quantization": 8},
        "genes": {"gene_expression": expression, "var": var},
        "connectivity": {"connectivities": connectivities},
        "compressed": {"compression": 6},
        "no_centroids": {"centroid_outlier_quantile": None},
        "single_category": {
            "obs": pd.DataFrame(
                {
                    "grp": pd.Categorical(["solo"] * n_cells),
                    "val": np.arange(n_cells) / n_cells,
                }
            )
        },
        "one_dimension": {
            "X_umap_2d": None,
            "X_umap_1d": (np.arange(n_cells) / n_cells).reshape(n_cells, 1),
        },
        "vectors": {
            "vector_fields": {
                "velocity_umap_2d": np.column_stack(
                    [np.arange(n_cells) / n_cells, np.arange(n_cells)[::-1] / n_cells]
                )
            },
            "vector_field_default": "velocity_umap",
        },
    }

    for name, overrides in cases.items():
        out_dir = tmp_path / name
        _prepare(out_dir, **overrides)
        for manifest_name in (
            "obs_manifest.json",
            "var_manifest.json",
            "connectivity_manifest.json",
            "dataset_identity.json",
        ):
            manifest_path = out_dir / manifest_name
            if not manifest_path.exists():
                continue
            # Parsed the way the browser parses it: no NaN, no Infinity.
            parsed = json.loads(
                manifest_path.read_text(encoding="utf-8"),
                parse_constant=_reject,
            )
            assert _all_finite(parsed), f"{name}/{manifest_name} carries a non-number"


def _reject(name: str) -> object:
    raise AssertionError(f"non-standard JSON constant {name!r}")


def _all_finite(value: object) -> bool:
    if isinstance(value, dict):
        return all(_all_finite(item) for item in value.values())
    if isinstance(value, list):
        return all(_all_finite(item) for item in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return True
