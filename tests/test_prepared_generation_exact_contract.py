from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.cli import _detect_data_format
from cellucid.prepare_data import prepare
from cellucid.server import CORSRequestHandler


def _prepare_kwargs(out_dir: Path, *, dimensions: int) -> dict[str, object]:
    if dimensions == 2:
        embedding = np.array(
            [
                [-3.0, 1.0],
                [0.5, 5.0],
                [8.0, -2.0],
            ],
            dtype=np.float32,
        )
    elif dimensions == 3:
        embedding = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 2.0],
                [0.0, 1.0, 4.0],
            ],
            dtype=np.float32,
        )
    else:
        raise AssertionError("test fixture supports only 2D and 3D")

    return {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        f"X_umap_{dimensions}d": embedding,
        "out_dir": out_dir,
        "dataset_name": "Atomic generation",
        "dataset_id": "atomic-generation",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }


def _snapshot(directory: Path) -> dict[str, bytes]:
    return {
        path.relative_to(directory).as_posix(): path.read_bytes()
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def _write_exported_dataset(directory: Path, dataset_id: str) -> None:
    directory.mkdir()
    (directory / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": dataset_id,
                "name": dataset_id,
            }
        ),
        encoding="utf-8",
    )
    (directory / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (directory / "points_2d.bin").write_bytes(b"\x00" * 8)


def _write_zarr_v2_root(directory: Path) -> None:
    directory.mkdir()
    (directory / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (directory / ".zattrs").write_text("{}", encoding="utf-8")


def test_prepare_exposes_only_canonical_artifact_names() -> None:
    parameters = inspect.signature(prepare).parameters
    assert "obs_manifest_filename" not in parameters
    assert "obs_binary_dirname" not in parameters
    assert "var_manifest_filename" not in parameters
    assert "var_binary_dirname" not in parameters


def test_force_replacement_removes_every_stale_capability(tmp_path: Path) -> None:
    output = tmp_path / "generation"
    initial = _prepare_kwargs(output, dimensions=3)
    initial["var"] = pd.DataFrame(index=["GeneA", "GeneB"])
    initial["gene_expression"] = np.array(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        dtype=np.float32,
    )
    initial["connectivities"] = np.array(
        [[0.0, 0.5, 0.0], [0.5, 0.0, 1.5], [0.0, 1.5, 0.0]],
        dtype=np.float64,
    )
    prepare(**initial)

    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    prepare(**replacement)

    assert (output / "points_2d.bin").is_file()
    assert not (output / "points_3d.bin").exists()
    assert not (output / "var_manifest.json").exists()
    assert not (output / "var").exists()
    assert not (output / "connectivity_manifest.json").exists()
    assert not (output / "connectivity").exists()
    identity = json.loads((output / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["embeddings"]["available_dimensions"] == [2]
    assert identity["stats"]["n_genes"] == 0
    assert identity["stats"]["has_connectivity"] is False


def test_failed_force_generation_preserves_the_prior_generation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import cellucid.prepare_data as prepare_module

    output = tmp_path / "generation"
    prepare(**_prepare_kwargs(output, dimensions=3))
    original = _snapshot(output)
    real_write_binary = prepare_module._write_binary
    calls = 0

    def fail_during_generation(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("synthetic staged write failure")
        return real_write_binary(*args, **kwargs)

    monkeypatch.setattr(prepare_module, "_write_binary", fail_during_generation)
    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    with pytest.raises(OSError, match="synthetic staged write failure"):
        prepare(**replacement)

    assert _snapshot(output) == original
    assert not list(tmp_path.glob(".generation.cellucid-*"))


def test_failed_initial_generation_leaves_no_partial_target(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import cellucid.prepare_data as prepare_module

    output = tmp_path / "generation"

    def fail_write(*_args, **_kwargs):
        raise OSError("synthetic initial write failure")

    monkeypatch.setattr(prepare_module, "_write_binary", fail_write)
    with pytest.raises(OSError, match="synthetic initial write failure"):
        prepare(**_prepare_kwargs(output, dimensions=2))

    assert not output.exists()
    assert not list(tmp_path.glob(".generation.cellucid-*"))


def test_concurrent_generation_for_the_same_target_is_rejected(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    lock_path = tmp_path / ".generation.cellucid.lock"
    lock_path.write_text("active owner\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="generation.*already active"):
        prepare(**_prepare_kwargs(output, dimensions=2))

    assert not output.exists()
    assert lock_path.read_text(encoding="utf-8") == "active owner\n"


def test_publish_failure_restores_the_prior_generation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "generation"
    prepare(**_prepare_kwargs(output, dimensions=3))
    original = _snapshot(output)
    real_rename = Path.rename

    def fail_staging_publish(path: Path, target: Path):
        if (
            path.name.startswith(".generation.cellucid-stage-")
            and Path(target) == output
        ):
            raise OSError("synthetic publish failure")
        return real_rename(path, target)

    monkeypatch.setattr(Path, "rename", fail_staging_publish)
    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    with pytest.raises(OSError, match="synthetic publish failure"):
        prepare(**replacement)

    assert _snapshot(output) == original
    assert not list(tmp_path.glob(".generation.cellucid-*"))


@pytest.mark.parametrize("marker", ["suffix-only", "partial-v2"])
def test_cli_rejects_incomplete_zarr_declarations(
    tmp_path: Path,
    marker: str,
) -> None:
    store = tmp_path / "fixture.zarr"
    store.mkdir()
    if marker == "partial-v2":
        (store / ".zattrs").write_text("{}", encoding="utf-8")
    assert _detect_data_format(store) == "unknown"


def test_cli_accepts_one_complete_zarr_declaration(tmp_path: Path) -> None:
    store = tmp_path / "fixture"
    _write_zarr_v2_root(store)
    assert _detect_data_format(store) == "zarr"


def test_exported_root_rejects_an_orphan_instead_of_omitting_it(
    tmp_path: Path,
) -> None:
    root = tmp_path / "exports"
    root.mkdir()
    _write_exported_dataset(root / "valid", "valid")
    (root / "orphan").mkdir()

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = root
    with pytest.raises(ValueError, match="orphan.*complete exported dataset"):
        handler._list_datasets()
    with pytest.raises(ValueError, match="orphan.*complete exported dataset"):
        _detect_data_format(root)


def test_exported_root_rejects_duplicate_dataset_ids(tmp_path: Path) -> None:
    root = tmp_path / "exports"
    root.mkdir()
    _write_exported_dataset(root / "first", "duplicate")
    _write_exported_dataset(root / "second", "duplicate")

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = root
    with pytest.raises(ValueError, match="duplicate dataset id"):
        handler._list_datasets()
