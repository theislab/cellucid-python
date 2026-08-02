"""What ``cellucid serve`` prints when it fails, and to whom it is addressed.

The person running ``cellucid serve`` is usually a wet-lab scientist, and the
failures they actually hit -- a forgotten ``--dataset-name``, a port another
program already holds, a path that is not a dataset -- are conditions they can
correct. Printing a Python traceback for those tells them nothing they can act
on and leaks the machine's directory layout into whatever they paste for help.

A defect inside Cellucid is the opposite case: its traceback is the bug report,
so it is printed in full, always, even under ``--quiet``.

These tests hold both halves of that line, and the shape of the output the
documentation tells readers to look for: the last line is a single ``Error:``
line and it carries the whole message.
"""

from __future__ import annotations

import argparse
import errno
import json
import socket
from pathlib import Path
from unittest import mock

import pytest

TRACEBACK_MARKER = "Traceback (most recent call last)"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SERVER_TUTORIAL = (
    REPOSITORY_ROOT / "docs" / "user_guide" / "web_app" / "b_data_loading" / "04_server_tutorial.md"
)


def _prepared_export(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "cli-export",
                "name": "CLI export",
                "description": "",
                "stats": {
                    "n_cells": 1,
                    "n_genes": 0,
                    "n_obs_fields": 0,
                    "n_categorical_fields": 0,
                    "n_continuous_fields": 0,
                    "has_connectivity": False,
                    "n_edges": None,
                },
                "embeddings": {
                    "available_dimensions": [2],
                    "default_dimension": 2,
                    "files": {"2d": "points_2d.bin"},
                },
                "obs_fields": [],
                "export_settings": {
                    "compression": None,
                    "var_quantization": None,
                    "obs_continuous_quantization": None,
                    "obs_categorical_dtype": "uint16",
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "obs_manifest.json").write_text(
        json.dumps(
            {
                "_format": "compact_v1",
                "n_points": 1,
                "centroid_outlier_quantile": None,
                "latent_key": "latent_space",
                "compression": None,
                "_obsSchemas": {},
                "_continuousFields": [],
                "_categoricalFields": [],
            }
        ),
        encoding="utf-8",
    )
    (root / "points_2d.bin").write_bytes(b"\x00" * 8)
    return root


def _zarr_store(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (root / ".zattrs").write_text("{}", encoding="utf-8")
    return root


def _error_lines(stderr: str) -> list[str]:
    return [line for line in stderr.splitlines() if line.startswith("Error:")]


def _assert_operator_report(stderr: str) -> str:
    """One ``Error:`` line, last, and no traceback above it."""
    assert TRACEBACK_MARKER not in stderr, stderr
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, stderr
    errors = _error_lines(stderr)
    assert len(errors) == 1, stderr
    assert lines[-1] == errors[0], stderr
    return errors[0]


def test_a_missing_dataset_identity_names_every_missing_flag_at_once(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    status = main(["serve", str(_zarr_store(tmp_path / "fixture.zarr")), "--quiet"])

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert message.startswith("Error: --dataset-name and --dataset-id are required")
    assert '--dataset-name "My study" --dataset-id my-study-v1' in message


def test_a_single_missing_identity_flag_names_only_that_flag(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    status = main(
        [
            "serve",
            str(_zarr_store(tmp_path / "fixture.zarr")),
            "--quiet",
            "--dataset-name",
            "My study",
        ]
    )

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert message.startswith("Error: --dataset-id is required")
    assert "--dataset-name and" not in message


def test_a_taken_port_is_reported_as_a_port_the_operator_can_change(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    holder = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    holder.bind(("127.0.0.1", 0))
    holder.listen(1)
    taken_port = holder.getsockname()[1]
    try:
        status = main(
            [
                "serve",
                str(_prepared_export(tmp_path / "export")),
                "--quiet",
                "--no-browser",
                "--no-web-ui",
                "--port",
                str(taken_port),
            ]
        )
    finally:
        holder.close()

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert str(taken_port) in message
    assert "already in use" in message
    assert "--port 0" in message


def test_an_unresolvable_host_is_reported_as_a_host_the_operator_can_change(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    status = main(
        [
            "serve",
            str(_prepared_export(tmp_path / "export")),
            "--quiet",
            "--no-browser",
            "--no-web-ui",
            "--host",
            "203.0.113.9",
        ]
    )

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "203.0.113.9" in message
    assert "--host 127.0.0.1" in message


def test_a_missing_path_is_reported_without_a_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    status = main(["serve", str(tmp_path / "does-not-exist"), "--quiet"])

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "Path not found" in message


def test_an_anndata_only_flag_on_a_prepared_export_is_reported_without_a_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    status = main(
        [
            "serve",
            str(_prepared_export(tmp_path / "export")),
            "--quiet",
            "--dataset-name",
            "Ignored",
        ]
    )

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "--dataset-name" in message
    assert "dataset_identity.json" in message


def test_an_unrecognized_input_is_reported_without_a_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    unrecognized = tmp_path / "notes.txt"
    unrecognized.write_text("not a dataset", encoding="utf-8")
    status = main(["serve", str(unrecognized), "--quiet"])

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "Unable to detect format" in message


def test_a_missing_optional_dependency_names_the_package_to_install(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    with mock.patch(
        "cellucid.anndata_server.serve_anndata",
        side_effect=ModuleNotFoundError("No module named 'zarr'", name="zarr"),
    ):
        status = main(
            [
                "serve",
                str(_zarr_store(tmp_path / "fixture.zarr")),
                "--quiet",
                "--dataset-name",
                "My study",
                "--dataset-id",
                "my-study-v1",
            ]
        )

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "pip install zarr" in message


def test_a_headless_machine_without_a_browser_is_told_to_use_no_browser(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    with mock.patch(
        "cellucid.server.serve",
        side_effect=RuntimeError("Could not open the browser for http://127.0.0.1:8765/"),
    ):
        status = main(["serve", str(_prepared_export(tmp_path / "export")), "--quiet"])

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "Could not open the browser" in message


def test_an_internal_defect_keeps_its_traceback_and_asks_for_a_bug_report(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    with mock.patch(
        "cellucid.server.serve",
        side_effect=AttributeError("'NoneType' object has no attribute 'shape'"),
    ):
        status = main(["serve", str(_prepared_export(tmp_path / "export")), "--quiet"])

    assert status == 1
    stderr = capsys.readouterr().err
    assert TRACEBACK_MARKER in stderr
    assert "AttributeError" in stderr
    errors = _error_lines(stderr)
    assert len(errors) == 1, stderr
    assert [line for line in stderr.splitlines() if line.strip()][-1] == errors[0]
    assert "defect in Cellucid" in errors[0]
    assert "https://github.com/theislab/cellucid-python/issues" in errors[0]


def test_an_operator_error_still_records_its_traceback_for_verbose_runs(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Nothing is destroyed; ``--verbose`` still surfaces the full stack."""
    from cellucid.cli import main

    status = main(["serve", str(_zarr_store(tmp_path / "fixture.zarr")), "--verbose"])

    assert status == 1
    stderr = capsys.readouterr().err
    assert TRACEBACK_MARKER in stderr
    errors = _error_lines(stderr)
    assert len(errors) == 1, stderr
    assert [line for line in stderr.splitlines() if line.strip()][-1] == errors[0]


def test_interrupting_the_server_is_not_reported_as_an_error(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from cellucid.cli import main

    with mock.patch("cellucid.server.serve", side_effect=KeyboardInterrupt):
        status = main(["serve", str(_prepared_export(tmp_path / "export")), "--quiet"])

    assert status == 130
    stderr = capsys.readouterr().err
    assert TRACEBACK_MARKER not in stderr
    assert _error_lines(stderr) == []
    assert "Interrupted by user" in stderr


@pytest.mark.parametrize(
    "error",
    [
        AttributeError("internal"),
        IndexError("internal"),
        AssertionError("internal"),
        NameError("internal"),
        ZeroDivisionError("internal"),
    ],
)
def test_every_defect_signature_keeps_its_traceback(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    error: Exception,
) -> None:
    from cellucid.cli import main

    with mock.patch("cellucid.server.serve", side_effect=error):
        status = main(["serve", str(_prepared_export(tmp_path / "export")), "--quiet"])

    assert status == 1
    stderr = capsys.readouterr().err
    assert TRACEBACK_MARKER in stderr
    assert "defect in Cellucid" in _error_lines(stderr)[0]


def test_the_tutorial_quotes_the_line_the_cli_actually_prints(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The transcripts on the server tutorial are generated output, not prose.

    They replaced two screenshots that had gone stale against the code and were
    printing this machine's home directory into a committed image. A quoted
    transcript only stays honest while something compares it to the command.
    """
    from cellucid.cli import _operating_system_message, main

    tutorial = SERVER_TUTORIAL.read_text(encoding="utf-8")

    h5ad = tmp_path / "lung_atlas.h5ad"
    h5ad.write_bytes(b"")
    assert main(["serve", str(h5ad), "--no-browser"]) == 1
    missing_identity = _error_lines(capsys.readouterr().err)[0]
    assert missing_identity in tutorial

    taken_port = _operating_system_message(
        OSError(errno.EADDRINUSE, "Address already in use"),
        argparse.Namespace(host="127.0.0.1", port=8765),
    )
    assert f"Error: {taken_port}" in tutorial

    unbindable_host = _operating_system_message(
        OSError(errno.EADDRNOTAVAIL, "Can't assign requested address"),
        argparse.Namespace(host="203.0.113.9", port=8765),
    )
    assert unbindable_host.split(".")[0] in " ".join(tutorial.split())

    # The page must not teach the behaviour that was removed.
    assert "ignore the traceback above it" not in tutorial
    assert "Errno 48" not in tutorial


@pytest.mark.parametrize(
    "error",
    [
        ValueError("Prepared-data path does not contain one complete current dataset."),
        TypeError("dataset_identity.json must contain a JSON object."),
        KeyError("obs field 'cell_type' not found. Available fields: ['leiden']"),
        OSError("Prepared artifact changed while it was opened."),
    ],
)
def test_every_reported_condition_stays_a_single_actionable_line(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    error: Exception,
) -> None:
    from cellucid.cli import main

    with mock.patch("cellucid.server.serve", side_effect=error):
        status = main(["serve", str(_prepared_export(tmp_path / "export")), "--quiet"])

    assert status == 1
    message = _assert_operator_report(capsys.readouterr().err)
    assert "defect in Cellucid" not in message
    # A ``KeyError`` prints as ``'...'`` unless the CLI unwraps it, and a
    # quoted sentence reads as a bug rather than as an instruction.
    assert not message.startswith("Error: \"'")
