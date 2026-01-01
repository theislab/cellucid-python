from __future__ import annotations

from pathlib import Path


def test_cli_detect_exported_root_dir(tmp_path: Path) -> None:
    from cellucid.cli import _detect_data_format

    exports_root = tmp_path / "exports"
    exports_root.mkdir()

    ds1 = exports_root / "ds1"
    ds1.mkdir()
    (ds1 / "dataset_identity.json").write_text('{"version": 2, "id": "ds1", "name": "ds1"}', encoding="utf-8")
    (ds1 / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (ds1 / "points_3d.bin").write_bytes(b"\x00")

    assert _detect_data_format(exports_root) == "exported"


def test_server_datasets_use_identity_id(tmp_path: Path) -> None:
    import json

    from cellucid.server import CORSRequestHandler

    dataset_dir = tmp_path / "my_folder_name"
    dataset_dir.mkdir()
    (dataset_dir / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (dataset_dir / "points_3d.bin").write_bytes(b"\x00")
    (dataset_dir / "dataset_identity.json").write_text(
        json.dumps({"version": 2, "id": "stable_dataset_id", "name": "Nice Name"}), encoding="utf-8"
    )

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = dataset_dir

    datasets = handler._list_datasets()
    assert datasets == [{"id": "stable_dataset_id", "path": "/", "name": "Nice Name"}]
