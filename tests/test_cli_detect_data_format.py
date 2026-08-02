from __future__ import annotations

import json
from pathlib import Path
from unittest import mock

import pytest


def test_cli_detect_exported_root_dir(tmp_path: Path) -> None:
    from cellucid.cli import _detect_data_format

    exports_root = tmp_path / "exports"
    exports_root.mkdir()

    ds1 = exports_root / "ds1"
    ds1.mkdir()
    (ds1 / "dataset_identity.json").write_text(
        '{"version": 2, "id": "ds1", "name": "ds1"}', encoding="utf-8"
    )
    (ds1 / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (ds1 / "points_3d.bin").write_bytes(b"\x00")

    assert _detect_data_format(exports_root) == "exported"


def test_unsupported_four_dimensional_points_are_not_detected_as_export(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _detect_data_format
    from cellucid.server import CORSRequestHandler

    dataset_dir = tmp_path / "unsupported"
    dataset_dir.mkdir()
    (dataset_dir / "dataset_identity.json").write_text(
        '{"version": 2, "id": "unsupported", "name": "Unsupported"}',
        encoding="utf-8",
    )
    (dataset_dir / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (dataset_dir / ("points_" + "4d.bin")).write_bytes(b"\x00")

    with pytest.raises(ValueError, match="complete exported dataset"):
        _detect_data_format(dataset_dir)
    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    assert handler._is_dataset_dir(dataset_dir) is False


def test_server_datasets_use_identity_id(tmp_path: Path) -> None:
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


@pytest.mark.parametrize(
    ("identity", "invalid_field"),
    [
        ({"version": 2, "name": "Nice Name"}, "dataset_id"),
        ({"version": 2, "id": "stable-id"}, "dataset_name"),
        ({"version": 2, "id": 7, "name": "Nice Name"}, "dataset_id"),
        ({"version": 2, "id": "stable-id", "name": ""}, "dataset_name"),
    ],
)
def test_server_rejects_malformed_identity_instead_of_using_directory_name(
    tmp_path: Path,
    identity: dict[str, object],
    invalid_field: str,
) -> None:
    from cellucid.server import CORSRequestHandler

    dataset_dir = tmp_path / "must-not-become-identity"
    dataset_dir.mkdir()
    (dataset_dir / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (dataset_dir / "points_3d.bin").write_bytes(b"\x00")
    (dataset_dir / "dataset_identity.json").write_text(
        json.dumps(identity),
        encoding="utf-8",
    )

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = dataset_dir

    with pytest.raises((TypeError, ValueError), match=invalid_field):
        handler._list_datasets()


def test_cli_requires_explicit_identity_for_direct_anndata(tmp_path: Path) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    args = create_parser().parse_args(["serve", str(store), "--quiet"])

    # Both flags are named in one message; reporting only the first would send
    # someone who supplied neither back to the terminal twice for one mistake.
    with pytest.raises(
        ValueError,
        match=r"--dataset-name and --dataset-id are required",
    ):
        _run_serve(args)


def test_cli_passes_exact_identity_without_backed_state_to_zarr(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    args = create_parser().parse_args(
        [
            "serve",
            str(store),
            "--quiet",
            "--dataset-name",
            "CLI fixture",
            "--dataset-id",
            "cli-fixture",
            "--obs-key",
            "cell_type",
            "--obs-key",
            "total_counts",
            "--vector-field-default",
            "drift_umap",
        ]
    )

    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)

    serve.assert_called_once_with(
        data=str(store),
        port=8765,
        host="127.0.0.1",
        open_browser=True,
        quiet=True,
        serve_web_ui=True,
        web_source_url="https://www.cellucid.com",
        web_cache_dir=None,
        allowed_hosts=None,
        latent_key=None,
        dataset_name="CLI fixture",
        dataset_id="cli-fixture",
        obs_keys=["cell_type", "total_counts"],
        vector_field_default="drift_umap",
    )


def _write_cli_export(directory: Path) -> None:
    directory.mkdir()
    (directory / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "cli-export",
                "name": "CLI export",
            }
        ),
        encoding="utf-8",
    )
    (directory / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (directory / "points_2d.bin").write_bytes(b"\x00" * 8)


def test_cli_no_command_and_parse_errors_return_nonzero_without_system_exit(
    capsys,
) -> None:
    from cellucid.cli import main

    assert main([]) == 2
    assert "usage:" in capsys.readouterr().err
    assert main(["--not-a-cellucid-option"]) == 2


def test_cli_runtime_failure_returns_nonzero_without_system_exit(
    tmp_path: Path,
    capsys,
) -> None:
    from cellucid.cli import main

    missing = tmp_path / "does-not-exist"
    assert main(["serve", str(missing), "--quiet"]) == 1
    assert "Path not found" in capsys.readouterr().err


def test_cli_quiet_and_verbose_are_mutually_exclusive() -> None:
    from cellucid.cli import create_parser

    with pytest.raises(SystemExit) as error:
        create_parser().parse_args(
            ["serve", "dataset", "--quiet", "--verbose"]
        )
    assert error.value.code == 2


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--latent-key", "latent"], "--latent-key"),
        (["--dataset-name", "Ignored"], "--dataset-name"),
        (["--dataset-id", "ignored"], "--dataset-id"),
        (["--obs-key", "cell_type"], "--obs-key"),
        (
            ["--vector-field-default", "velocity_umap"],
            "--vector-field-default",
        ),
    ],
)
def test_cli_rejects_every_anndata_only_option_for_prepared_exports(
    tmp_path: Path,
    arguments: list[str],
    message: str,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    export = tmp_path / "export"
    _write_cli_export(export)
    args = create_parser().parse_args(
        ["serve", str(export), "--quiet", *arguments]
    )
    with (
        mock.patch("cellucid.server.serve") as serve,
        pytest.raises(ValueError, match=message),
    ):
        _run_serve(args)
    serve.assert_not_called()


def test_cli_rejects_removed_backed_mode_without_dispatch(tmp_path: Path) -> None:
    from cellucid.cli import main

    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        status = main(
            [
                "serve",
                str(store),
                "--quiet",
                "--dataset-name",
                "CLI fixture",
                "--dataset-id",
                "cli-fixture",
                "--no-backed",
            ]
        )
    assert status == 2
    serve.assert_not_called()


def test_cli_does_not_reclassify_import_error_by_message_text(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    args = create_parser().parse_args(
        [
            "serve",
            str(store),
            "--quiet",
            "--dataset-name",
            "CLI fixture",
            "--dataset-id",
            "cli-fixture",
        ]
    )
    with (
        mock.patch(
            "cellucid.anndata_server.serve_anndata",
            side_effect=ImportError("pip install zarr is unrelated text"),
        ),
        pytest.raises(ImportError, match="unrelated text"),
    ):
        _run_serve(args)


def test_cli_forwards_every_repeated_allowed_host_to_the_prepared_server(
    tmp_path: Path,
) -> None:
    """`--allowed-host` is the operator's explicit reverse-proxy declaration."""
    from cellucid.cli import _run_serve, create_parser

    export = tmp_path / "export"
    _write_cli_export(export)
    args = create_parser().parse_args(
        [
            "serve",
            str(export),
            "--quiet",
            "--allowed-host",
            "hub.example.org",
            "--allowed-host",
            "other.example.org",
        ]
    )

    with mock.patch("cellucid.server.serve") as serve:
        _run_serve(args)

    serve.assert_called_once_with(
        data_dir=str(export),
        port=8765,
        host="127.0.0.1",
        open_browser=True,
        quiet=True,
        serve_web_ui=True,
        web_source_url="https://www.cellucid.com",
        web_cache_dir=None,
        allowed_hosts=["hub.example.org", "other.example.org"],
    )


def test_cli_forwards_every_repeated_allowed_host_to_the_anndata_server(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    args = create_parser().parse_args(
        [
            "serve",
            str(store),
            "--quiet",
            "--dataset-name",
            "CLI fixture",
            "--dataset-id",
            "cli-fixture",
            "--allowed-host",
            "hub.example.org",
        ]
    )

    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)

    assert serve.call_args.kwargs["allowed_hosts"] == ["hub.example.org"]


def test_cli_declares_no_allowed_host_by_default(tmp_path: Path) -> None:
    """No flag must mean exactly the loopback-only rule, not a wider one."""
    from cellucid.cli import _run_serve, create_parser

    export = tmp_path / "export"
    _write_cli_export(export)
    args = create_parser().parse_args(["serve", str(export), "--quiet"])

    with mock.patch("cellucid.server.serve") as serve:
        _run_serve(args)

    assert serve.call_args.kwargs["allowed_hosts"] is None


def test_cli_reports_an_unsupported_allowed_host_without_serving(
    tmp_path: Path,
    capsys,
) -> None:
    from cellucid.cli import main

    export = tmp_path / "export"
    _write_cli_export(export)
    status = main(
        ["serve", str(export), "--quiet", "--allowed-host", "*.example.org"]
    )

    assert status == 1
    assert "wildcard" in capsys.readouterr().err
