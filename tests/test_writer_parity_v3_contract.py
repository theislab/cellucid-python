"""Byte parity of the two writers on export format v3.

``tests/fixtures/writer_parity_v3.json`` holds one input -- two cells, five
embeddings (a 1-D vector, three 2D items including a wide array read at
columns ``[3, 1]``, and a 3D item), two vector fields bound to embedding ids,
duplicate and empty gene names in the default column and in an alternative one,
and a one-edge graph -- together with every payload byte and every manifest
the Python writer produces from it. cellucid-r keeps the same file, pins the
same SHA-256, and asserts that ``cellucid_prepare()`` produces the same bytes
and equal manifests, which is what makes the two writers one format.

Every coordinate is dyadic and every axis range a power of two, so the
normalization, the centroids and the vector scaling are exact in float64 and
no summation order can make the writers disagree.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import cellucid
from cellucid import Embedding, prepare

FIXTURE = Path(__file__).parent / "fixtures" / "writer_parity_v3.json"

#: Pinned identically in cellucid-r (tests/testthat/test-writer-parity-v3.R).
FIXTURE_SHA256 = "44779fb1fce0955d6bdbdfc85df9458c8f7edff945cec07a4c9c65fbc7bb9aad"


def _fixture() -> dict[str, Any]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _prepare_arguments(spec: dict[str, Any]) -> dict[str, Any]:
    obs = pd.DataFrame(
        {
            key: (
                pd.Categorical(field["values"])
                if field["kind"] == "category"
                else np.asarray(field["values"], dtype=np.float64)
            )
            for key, field in spec["obs"].items()
        }
    )
    var = pd.DataFrame(spec["var"]["columns"], index=spec["var"]["index"])
    return {
        "latent_space": np.asarray(spec["latent_space"], dtype=np.float64),
        "obs": obs,
        "var": var,
        "gene_expression": np.asarray(spec["gene_expression"], dtype=np.float64),
        "var_gene_id_column": spec["var_gene_id_column"],
        "var_gene_name_columns": spec["var_gene_name_columns"],
        "connectivities": np.asarray(spec["connectivities"], dtype=np.float64),
        "centroid_min_points": spec["centroid_min_points"],
        "obs_categorical_dtype": spec["obs_categorical_dtype"],
        "dataset_name": spec["dataset_name"],
        "dataset_id": spec["dataset_id"],
        "created_at": spec["created_at"],
        "embeddings": {
            key: Embedding(
                np.asarray(item["values"], dtype=np.float64),
                columns=item["columns"],
                label=item["label"],
                source_key=item["source_key"],
            )
            for key, item in spec["embeddings"].items()
        },
        "embedding_defaults": spec["embedding_defaults"],
        "vector_fields": {
            field: {key: np.asarray(values, dtype=np.float64) for key, values in by.items()}
            for field, by in spec["vector_fields"].items()
        },
        "vector_field_default": spec["vector_field_default"],
    }


def test_the_parity_fixture_is_the_pinned_shared_copy() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256


def test_prepare_reproduces_every_parity_payload_byte_and_manifest(tmp_path: Path) -> None:
    fixture = _fixture()
    out = tmp_path / "export"
    prepare(out_dir=out, **_prepare_arguments(fixture["inputs"]))

    written = sorted(
        path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()
    )
    expected_payloads = fixture["expected"]["payloads"]
    expected_manifests = fixture["expected"]["manifests"]
    assert written == sorted([*expected_payloads, *expected_manifests])

    for relative, hex_bytes in expected_payloads.items():
        assert (out / relative).read_bytes() == bytes.fromhex(hex_bytes), relative

    for relative, manifest in expected_manifests.items():
        parsed = json.loads((out / relative).read_text(encoding="utf-8"))
        if relative == "dataset_identity.json":
            assert parsed.pop("cellucid_data_version") == cellucid.__version__
        # Key order is part of what both writers agree on.
        assert json.dumps(parsed) == json.dumps(manifest), relative


def test_the_parity_fixture_exercises_every_v3_surface() -> None:
    """The fixture is only evidence of parity if it reaches every new rule."""
    fixture = _fixture()
    manifests = fixture["expected"]["manifests"]
    identity = manifests["dataset_identity.json"]
    items = identity["embeddings"]["items"]
    assert identity["version"] == 3
    assert [item["dimension"] for item in items.values()] == [1, 2, 2, 2, 3]
    assert any(item["source"]["columns"] != list(range(item["dimension"])) for item in items.values())
    assert identity["embeddings"]["defaults"] == {"1d": "X_line", "2d": "X_umap", "3d": "X_umap_3d"}
    assert {field for field in identity["vector_fields"]["fields"]} == {"drift", "velocity"}

    var_manifest = manifests["var_manifest.json"]
    sets = var_manifest["gene_name_sets"]
    assert len(sets) == 2
    assert all(name_set["renamed"] and name_set["filled_from_index"] for name_set in sets)
    assert sets[1]["names"] is not None and None in sets[1]["names"]

    categorical = manifests["obs_manifest.json"]["_categoricalFields"][0]
    assert list(categorical[5]) == list(items)
    payloads = fixture["expected"]["payloads"]
    assert {path.split("/")[0] for path in payloads} == {
        "connectivity",
        "embeddings",
        "obs",
        "var",
        "vectors",
    }
