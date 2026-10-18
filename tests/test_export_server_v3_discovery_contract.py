"""Exact contract for serving format-v3 exports with ``cellucid serve``.

A directory is an export when it holds ``dataset_identity.json``,
``obs_manifest.json`` and an ``embeddings/`` directory; it is admitted only
with a version-3 identity, so an export written in an earlier format is told to
re-export. Which embedding and vector files it may serve is decided by the
identity: every ``embeddings.items[*].file`` and every vector file, which is
keyed by the embedding id it is drawn over.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare
from cellucid.server._artifacts import _declared_dataset_artifacts
from cellucid.server._datasets import _list_exported_datasets


def _export(out_dir: Path, *, dataset_id: str = "served-v3", compression: int | None = None) -> None:
    plane = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    prepare(
        latent_space=plane,
        obs=pd.DataFrame({"group": pd.Categorical(["a", "b", "a"])}),
        embeddings={"X_umap": plane, "X_tsne": plane * 3.0, "X_line": plane[:, 0]},
        embedding_defaults={"2d": "X_umap"},
        vector_fields={"velocity": {"X_tsne": plane * 0.1, "X_line": plane[:, :1]}},
        out_dir=out_dir,
        dataset_name=f"Served {dataset_id}",
        dataset_id=dataset_id,
        created_at="2026-01-01T00:00:00Z",
        obs_categorical_dtype="uint8",
        centroid_min_points=1,
        compression=compression,
    )


@pytest.mark.parametrize("compression", [None, 6])
def test_the_inventory_is_every_declared_embedding_and_vector_file(
    tmp_path: Path,
    compression: int | None,
) -> None:
    _export(tmp_path, compression=compression)
    suffix = ".gz" if compression is not None else ""

    declared = _declared_dataset_artifacts(tmp_path, include_published_state=False)
    assert {path for path in declared if path.startswith(("embeddings/", "vectors/"))} == {
        f"embeddings/0_1d.bin{suffix}",
        f"embeddings/1_2d.bin{suffix}",
        f"embeddings/2_2d.bin{suffix}",
        f"vectors/0_0.bin{suffix}",
        f"vectors/0_1.bin{suffix}",
    }
    assert declared == {
        path.relative_to(tmp_path).as_posix() for path in tmp_path.rglob("*") if path.is_file()
    }


def test_a_single_export_and_a_root_of_exports_are_discovered(tmp_path: Path) -> None:
    single = tmp_path / "single"
    _export(single)
    assert _list_exported_datasets(single) == [
        {"id": "served-v3", "path": "/", "name": "Served served-v3"}
    ]

    root = tmp_path / "root"
    _export(root / "alpha", dataset_id="alpha")
    _export(root / "beta", dataset_id="beta")
    assert [entry["path"] for entry in _list_exported_datasets(root)] == ["/alpha/", "/beta/"]


def test_a_directory_with_no_export_marker_is_not_an_export(tmp_path: Path) -> None:
    (tmp_path / "notes").mkdir()
    (tmp_path / "points_2d.bin").write_bytes(b"\x00" * 8)
    assert _list_exported_datasets(tmp_path) == []


def test_an_export_in_the_previous_format_is_told_to_re_export(tmp_path: Path) -> None:
    (tmp_path / "dataset_identity.json").write_text(
        json.dumps({"version": 2, "id": "old", "name": "Old"}),
        encoding="utf-8",
    )
    (tmp_path / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (tmp_path / "points_2d.bin").write_bytes(b"\x00" * 8)

    with pytest.raises(
        ValueError,
        match=r"version must be exactly 3\. Re-export the dataset with the current prepare\(\)",
    ):
        _list_exported_datasets(tmp_path)


@pytest.mark.parametrize("state", ["missing", "empty", "file"])
def test_an_export_without_embedding_payloads_is_incomplete(tmp_path: Path, state: str) -> None:
    _export(tmp_path)
    embeddings = tmp_path / "embeddings"
    for path in embeddings.iterdir():
        path.unlink()
    embeddings.rmdir()
    if state == "empty":
        embeddings.mkdir()
    elif state == "file":
        embeddings.write_bytes(b"")

    with pytest.raises(ValueError, match=r"is not a complete exported dataset: it has no embeddings/"):
        _list_exported_datasets(tmp_path)


def test_an_embeddings_directory_that_is_a_symbolic_link_is_refused(tmp_path: Path) -> None:
    export = tmp_path / "export"
    _export(export)
    elsewhere = tmp_path / "elsewhere"
    (export / "embeddings").rename(elsewhere)
    (export / "embeddings").symlink_to(elsewhere, target_is_directory=True)

    with pytest.raises(ValueError, match=r"is not a complete exported dataset"):
        _list_exported_datasets(export)


@pytest.mark.parametrize(
    ("embeddings", "message"),
    [
        ({"files": {"2d": "points_2d.bin"}}, r"embeddings\.items must be a JSON object"),
        ({"items": {}}, r"embeddings\.items must declare an embedding"),
        ({"items": {"X_umap": "embeddings/0_2d.bin"}}, r"embedding 'X_umap' is invalid"),
        ({"items": {"X_umap": {"file": "/abs.bin"}}}, r"one exact relative POSIX path"),
        ({"items": {"X_umap": {"file": "../0_2d.bin"}}}, r"one exact relative POSIX path"),
    ],
)
def test_the_inventory_refuses_a_malformed_embedding_declaration(
    tmp_path: Path,
    embeddings: dict,
    message: str,
) -> None:
    _export(tmp_path)
    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))
    identity["embeddings"] = embeddings
    (tmp_path / "dataset_identity.json").write_text(json.dumps(identity), encoding="utf-8")

    with pytest.raises((TypeError, ValueError), match=message):
        _declared_dataset_artifacts(tmp_path, include_published_state=False)


def test_vector_files_are_read_per_bound_embedding(tmp_path: Path) -> None:
    _export(tmp_path)
    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["vector_fields"]["fields"]["velocity"]["files"] == {
        "X_line": "vectors/0_0.bin",
        "X_tsne": "vectors/0_1.bin",
    }
    identity["vector_fields"]["fields"]["velocity"]["files"]["X_tsne"] = "vectors/../x.bin"
    (tmp_path / "dataset_identity.json").write_text(json.dumps(identity), encoding="utf-8")

    with pytest.raises(ValueError, match=r"vector field 'velocity' file for 'X_tsne'"):
        _declared_dataset_artifacts(tmp_path, include_published_state=False)
