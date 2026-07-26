from __future__ import annotations

import inspect
import json
import re
from pathlib import Path

import pytest

from cellucid._server_base import (
    register_event_callback,
    route_event,
    unregister_event_callback,
)
from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer
from cellucid.jupyter import CellucidViewer, HookRegistry
from cellucid.prepare_data import generate_datasets_manifest
from cellucid.server import CellucidServer

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RETIRED_PATH_LANGUAGE = re.compile(
    r"\bfallback\b|fall(?:ing)?[\s-]+back|best[\s_-]*effort|"
    r"\blegacy\b|\bdeprecated\b|deprecation|backwards?[\s_-]*compat|"
    r"\bshim\b",
    re.IGNORECASE,
)
CURRENT_TEXT_SUFFIXES = {
    ".json",
    ".md",
    ".py",
    ".rst",
    ".sh",
    ".toml",
    ".yaml",
    ".yml",
}


def _identity(dataset_id: str, name: str = "Dataset") -> dict[str, object]:
    return {
        "version": 2,
        "id": dataset_id,
        "name": name,
        "stats": {
            "n_cells": 3,
            "n_genes": 2,
        },
    }


def _write_identity(directory: Path, identity: dict[str, object]) -> None:
    directory.mkdir()
    (directory / "dataset_identity.json").write_text(
        json.dumps(identity),
        encoding="utf-8",
    )


def test_shipped_python_repository_surfaces_contain_only_current_contract_language() -> None:
    roots = [
        REPOSITORY_ROOT / "src",
        REPOSITORY_ROOT / "docs",
        REPOSITORY_ROOT / "scripts",
        REPOSITORY_ROOT / ".github",
    ]
    files = [
        REPOSITORY_ROOT / "README.md",
        REPOSITORY_ROOT / "CHANGELOG.md",
        REPOSITORY_ROOT / "pyproject.toml",
    ]
    for root in roots:
        files.extend(
            path
            for path in root.rglob("*")
            if path.is_file()
            and path.suffix in CURRENT_TEXT_SUFFIXES
            and "_build" not in path.parts
            and ".cellucid-maintenance" not in path.parts
        )

    violations: list[str] = []
    for path in files:
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            start=1,
        ):
            if RETIRED_PATH_LANGUAGE.search(line):
                violations.append(
                    f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: "
                    f"{line.strip()}"
                )

    assert violations == []


def test_mypy_uses_each_ci_interpreter_and_ci_keeps_the_minimum_version() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    mypy_section = pyproject.split("[tool.mypy]", maxsplit=1)[1].split(
        "\n[",
        maxsplit=1,
    )[0]
    assert re.search(r"(?m)^\s*python_version\s*=", mypy_section) is None

    workflow = (REPOSITORY_ROOT / ".github/workflows/test.yml").read_text(
        encoding="utf-8"
    )
    assert 'python: "3.10"' in workflow
    assert "python -m mypy src/cellucid" in workflow


def test_release_version_and_embedding_dimensions_are_coherent() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    version_match = re.search(
        r'(?m)^version = "([^"]+)"\s+# CELLUCID_VERSION$',
        pyproject,
    )
    assert version_match is not None
    version = version_match.group(1)

    changelog = (REPOSITORY_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    citation = (REPOSITORY_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    docs_conf = (REPOSITORY_ROOT / "docs/conf.py").read_text(encoding="utf-8")
    assert f"## [{version}]" in changelog
    assert f"/compare/v{version}...HEAD" in changelog
    assert f"/releases/tag/v{version}" in changelog
    assert f"version: {version}" in citation
    assert f'release = version = "{version}"' in docs_conf
    assert "1D/2D/3D/4D" not in changelog
    assert "Multi-dimensional embedding exports (1D/2D/3D)" in changelog


def test_event_callback_failure_propagates_to_the_http_boundary() -> None:
    viewer_id = "exact-callback"
    viewer_token = "exact-callback-token"

    def fail(_event: dict[str, object]) -> None:
        raise RuntimeError("callback failed")

    register_event_callback(viewer_id, viewer_token, fail)
    try:
        with pytest.raises(RuntimeError, match="callback failed"):
            route_event(viewer_id, viewer_token, {"type": "selection"})
    finally:
        unregister_event_callback(viewer_id, viewer_token)


def test_hook_registry_delivers_each_event_once_and_propagates_failure() -> None:
    hooks = HookRegistry()
    deliveries: list[dict[str, object]] = []
    hooks.register("message", deliveries.append)
    hooks.trigger("message", {"value": 7})
    assert deliveries == [{"event": "message", "value": 7}]

    def fail(_event: dict[str, object]) -> None:
        raise RuntimeError("hook failed")

    hooks.register("selection", fail)
    with pytest.raises(RuntimeError, match="hook failed"):
        hooks.trigger("selection", {"cells": [1]})


def test_send_message_before_display_is_an_explicit_error() -> None:
    viewer = object.__new__(CellucidViewer)
    viewer._displayed = False
    with pytest.raises(RuntimeError, match="display"):
        viewer.send_message({"type": "ping", "requestId": "request-1"})


@pytest.mark.parametrize("server_type", [CellucidServer, AnnDataServer])
def test_starting_a_running_server_is_an_explicit_error(server_type: type) -> None:
    server = object.__new__(server_type)
    server._running = True
    with pytest.raises(RuntimeError, match="already running"):
        server.start()


@pytest.mark.parametrize(
    "path_kind",
    ["wrong-file-extension", "incomplete-zarr-directory"],
)
def test_anndata_file_loader_rejects_undeclared_formats_before_reader_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    path_kind: str,
) -> None:
    import anndata

    calls: list[str] = []

    def forbidden_h5ad(*_args, **_kwargs):
        calls.append("h5ad")
        raise AssertionError("read_h5ad must not run")

    def forbidden_zarr(*_args, **_kwargs):
        calls.append("zarr")
        raise AssertionError("read_zarr must not run")

    monkeypatch.setattr(anndata, "read_h5ad", forbidden_h5ad)
    monkeypatch.setattr(anndata, "read_zarr", forbidden_zarr)

    if path_kind == "wrong-file-extension":
        path = tmp_path / "dataset.bin"
        path.write_bytes(b"not h5ad")
    else:
        path = tmp_path / "dataset.zarr"
        path.mkdir()
        (path / ".zattrs").write_text("{}", encoding="utf-8")

    with pytest.raises(ValueError, match=r"\.h5ad|Zarr"):
        AnnDataAdapter.from_file(
            path,
            dataset_name="Exact loader",
            dataset_id="exact-loader",
        )
    assert calls == []


def test_dataset_manifest_has_one_required_default_and_no_repair_modes() -> None:
    parameters = inspect.signature(generate_datasets_manifest).parameters
    assert "output_filename" not in parameters
    assert parameters["default_dataset"].default is inspect.Parameter.empty


def test_dataset_manifest_rejects_orphan_directories_atomically(
    tmp_path: Path,
) -> None:
    _write_identity(tmp_path / "valid", _identity("valid"))
    (tmp_path / "orphan").mkdir()
    manifest_path = tmp_path / "datasets.json"
    original = b'{"existing":true}\n'
    manifest_path.write_bytes(original)

    with pytest.raises(ValueError, match="orphan.*dataset_identity"):
        generate_datasets_manifest(tmp_path, default_dataset="valid")
    assert manifest_path.read_bytes() == original


def test_dataset_manifest_rejects_missing_or_unknown_default(
    tmp_path: Path,
) -> None:
    _write_identity(tmp_path / "valid", _identity("valid"))

    with pytest.raises(ValueError, match="default_dataset.*missing"):
        generate_datasets_manifest(tmp_path, default_dataset="missing")
    assert not (tmp_path / "datasets.json").exists()


def test_dataset_manifest_rejects_duplicate_dataset_ids(
    tmp_path: Path,
) -> None:
    _write_identity(tmp_path / "first", _identity("duplicate", "First"))
    _write_identity(tmp_path / "second", _identity("duplicate", "Second"))

    with pytest.raises(ValueError, match="duplicate dataset id"):
        generate_datasets_manifest(tmp_path, default_dataset="duplicate")
    assert not (tmp_path / "datasets.json").exists()
