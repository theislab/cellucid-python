from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest


def _write_zarr_store(store: Path) -> Path:
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    return store


def test_cli_detects_every_anndata_layout_and_nothing_else(tmp_path: Path) -> None:
    from cellucid.cli import _detect_data_format

    h5ad = tmp_path / "data.h5ad"
    h5ad.write_bytes(b"")
    assert _detect_data_format(h5ad) == "h5ad"
    assert _detect_data_format(_write_zarr_store(tmp_path / "data.zarr")) == "zarr"
    v3 = tmp_path / "v3.zarr"
    v3.mkdir()
    (v3 / "zarr.json").write_text('{"zarr_format":3,"node_type":"group"}', encoding="utf-8")
    assert _detect_data_format(v3) == "zarr"

    # A plain directory of other files is not an AnnData store.
    folder = tmp_path / "folder"
    folder.mkdir()
    (folder / "notes.json").write_text("{}", encoding="utf-8")
    assert _detect_data_format(folder) == "unknown"
    assert _detect_data_format(tmp_path / "missing") == "unknown"
    other = tmp_path / "data.csv"
    other.write_text("a,b\n", encoding="utf-8")
    assert _detect_data_format(other) == "unknown"


def test_cli_reports_a_directory_that_is_not_a_store(tmp_path: Path, capsys) -> None:
    from cellucid.cli import main

    folder = tmp_path / "folder"
    folder.mkdir()
    assert main(["serve", str(folder), "--quiet"]) == 1
    assert "Unable to detect format" in capsys.readouterr().err


def test_cli_requires_an_identity_a_store_does_not_declare(tmp_path: Path) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = _write_zarr_store(tmp_path / "fixture.zarr")
    args = create_parser().parse_args(["serve", str(store), "--quiet"])

    # Both are named in one message; reporting only the first would send
    # someone who supplied neither back to the terminal twice for one mistake.
    with pytest.raises(
        ValueError,
        match=r"dataset_name and dataset_id are required.*--dataset-name/--dataset-id",
    ):
        _run_serve(args)


def test_cli_passes_exact_identity_without_backed_state_to_zarr(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = _write_zarr_store(tmp_path / "fixture.zarr")
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
        normalize_embeddings=None,
        centroid_min_points=None,
        dataset_name="CLI fixture",
        dataset_id="cli-fixture",
        obs_keys=["cell_type", "total_counts"],
        gene_id_column=None,
        embedding_keys=None,
        knn_key=None,
        expression_source=None,
        vector_field_default="drift_umap",
    )


@pytest.mark.parametrize(
    ("flags", "normalize_embeddings", "centroid_min_points"),
    [
        (["--no-normalize-embeddings"], False, None),
        (["--normalize-embeddings", "--centroid-min-points", "3"], True, 3),
    ],
)
def test_cli_forwards_embedding_scaling_and_outlier_category_size(
    tmp_path: Path,
    flags: list[str],
    normalize_embeddings: bool,
    centroid_min_points: int | None,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = _write_zarr_store(tmp_path / "fixture.zarr")
    args = create_parser().parse_args(
        ["serve", str(store), "--quiet", "--dataset-name", "CLI", "--dataset-id", "cli", *flags]
    )

    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)

    assert serve.call_args.kwargs["normalize_embeddings"] is normalize_embeddings
    assert serve.call_args.kwargs["centroid_min_points"] == centroid_min_points


@pytest.mark.parametrize("value", ["0", "-2", "1.5", "ten"])
def test_cli_refuses_a_centroid_min_points_that_is_not_a_positive_integer(
    value: str,
    capsys,
) -> None:
    from cellucid.cli import create_parser

    with pytest.raises(SystemExit) as error:
        create_parser().parse_args(["serve", "dataset", "--centroid-min-points", value])
    assert error.value.code == 2
    assert (
        f"argument --centroid-min-points: expected a positive integer; got {value!r}"
        in capsys.readouterr().err
    )


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


def test_cli_rejects_removed_backed_mode_without_dispatch(tmp_path: Path) -> None:
    from cellucid.cli import main

    store = _write_zarr_store(tmp_path / "fixture.zarr")
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

    store = _write_zarr_store(tmp_path / "fixture.zarr")
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


def test_cli_forwards_every_repeated_allowed_host_to_the_anndata_server(
    tmp_path: Path,
) -> None:
    from cellucid.cli import _run_serve, create_parser

    store = _write_zarr_store(tmp_path / "fixture.zarr")
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

    store = _write_zarr_store(tmp_path / "fixture.zarr")
    args = create_parser().parse_args(["serve", str(store), "--quiet"])

    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)

    assert serve.call_args.kwargs["allowed_hosts"] is None


def test_cli_reports_an_unsupported_allowed_host_without_serving(
    tmp_path: Path,
    capsys,
) -> None:
    from cellucid.cli import main

    store = _write_zarr_store(tmp_path / "fixture.zarr")
    status = main(
        ["serve", str(store), "--quiet", "--allowed-host", "*.example.org"]
    )

    assert status == 1
    assert "wildcard" in capsys.readouterr().err
