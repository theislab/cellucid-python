from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

import cellucid
from cellucid._viewer_channel import (
    register_event_callback,
    route_event,
    unregister_event_callback,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.jupyter import HOOK_EVENTS, AnnDataViewer, HookRegistry

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
                    f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: {line.strip()}"
                )

    assert violations == []


def test_readme_cross_references_the_current_cellucid_ecosystem() -> None:
    readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
    for url in (
        "https://docs.cellucid.com/en/latest/user_guide/python_package/index.html",
        "https://github.com/theislab/cellucid",
        "https://github.com/theislab/cellucid-r",
        "https://github.com/theislab/cellucid-datasets",
        "https://github.com/theislab/cellucid-demo-custom-datasets",
        "https://github.com/theislab/cellucid-annotation",
    ):
        assert url in readme


def test_readme_renders_completely_outside_the_repository() -> None:
    """PyPI serves the long description with no repository to resolve against."""
    readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
    relative_targets = [
        target
        for target in re.findall(r"\]\(([^)\s]+)\)", readme)
        if not target.startswith(("https://", "http://", "#", "mailto:"))
    ]
    assert relative_targets == []


def test_readme_claims_only_what_the_project_can_evidence() -> None:
    readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
    # No conda-forge feedstock exists for this package.
    assert "conda-forge" not in readme
    assert "anaconda.org" not in readme
    # Unmeasured scale and speed claims must not be advertised.
    for unmeasured in ("millions", "10x faster", "10× faster", "fastest"):
        assert unmeasured.lower() not in readme.lower()


def test_mypy_uses_each_ci_interpreter_and_ci_keeps_the_supported_boundaries() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    mypy_section = pyproject.split("[tool.mypy]", maxsplit=1)[1].split(
        "\n[",
        maxsplit=1,
    )[0]
    assert re.search(r"(?m)^\s*python_version\s*=", mypy_section) is None

    workflow = (REPOSITORY_ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    assert 'python: "3.11"' in workflow
    assert 'python: "3.14"' in workflow
    assert 'python: "3.10"' not in workflow
    assert "python -m mypy src/cellucid" in workflow


def test_ci_checks_out_every_sibling_the_tests_compare_with() -> None:
    """The cross-repository tests (tests/_web_root.py) skip only where a
    sibling is not checked out; CI checks out each one beside this repository
    and requires it, so there they always run (a test reading cellucid-r
    without its checkout failed every job, run 43106597330)."""
    from _web_root import _SIBLINGS

    workflow = (REPOSITORY_ROOT / ".github/workflows/test.yml").read_text(encoding="utf-8")
    assert re.search(r'(?m)^      CELLUCID_SIBLINGS_REQUIRED: "1"$', workflow)
    assert re.search(r"(?m)^          path: cellucid-python$", workflow)
    for repository in _SIBLINGS:
        assert re.search(
            rf"(?m)^          repository: theislab/{re.escape(repository)}\n          path: {re.escape(repository)}$",
            workflow,
        ), repository
    # The web app's sparse checkout holds every tree the tests read of it.
    sparse = re.search(r"repository: theislab/cellucid\n(?:.*\n)*?\s+sparse-checkout: \|\n((?:\s{12}.*\n)+)", workflow)
    assert sparse is not None
    assert sorted(line.strip() for line in sparse.group(1).strip().splitlines()) == [
        "assets/external/",
        "assets/js/",
        "assets/wasm/",
        "assets/workers/",
        "scripts/",
        "tests/fixtures/",
    ]
    assert "actions/setup-node@" in workflow


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
    assert "## [Unreleased]" not in changelog
    assert "[Unreleased]:" not in changelog
    assert f"/releases/tag/v{version}" in changelog
    assert f"version: {version}" in citation
    assert "date-released:" not in citation
    assert f"## [{version}]\n" in changelog
    assert f'release = version = "{version}"' in docs_conf
    assert "1D/2D/3D/4D" not in changelog
    assert "Multi-dimensional embedding exports (1D/2D/3D)" in changelog
    assert cellucid.__doc__ is not None
    assert "1D, 2D, and 3D single-cell WebGL visualization toolkit" in cellucid.__doc__

    validator = REPOSITORY_ROOT / "scripts" / "validate_release.py"
    valid = subprocess.run(
        [sys.executable, str(validator), "--tag", f"v{version}"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert valid.returncode == 0, valid.stderr
    wrong = subprocess.run(
        [sys.executable, str(validator), "--tag", "v999.999.999"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert wrong.returncode != 0
    assert f"must equal v{version}" in wrong.stderr


def test_pypi_release_requires_gated_artifacts_and_trusted_identity() -> None:
    workflow = (REPOSITORY_ROOT / ".github/workflows/pypi-publish.yml").read_text(encoding="utf-8")
    # One trigger: a published GitHub Release. A tag push publishes nothing.
    assert "\non:\n  release:\n    types: [published]\n\npermissions:" in workflow
    assert "push:" not in workflow
    assert "tags:" not in workflow
    assert "workflow_dispatch" not in workflow
    assert "github.ref_name" not in workflow
    assert '--release-event "${GITHUB_EVENT_PATH}"' in workflow
    assert '[[ "${GITHUB_REF}" == "refs/tags/${RELEASE_TAG}" ]]' in workflow
    assert "PYPI_API_TOKEN" not in workflow
    assert "password:" not in workflow
    assert "id-token: write" in workflow
    assert "name: pypi" in workflow
    assert "release-provenance:" in workflow
    assert workflow.count("needs: release-provenance") == 2
    assert workflow.count("fetch-depth: 0") == 1
    assert "refs/heads/main:refs/remotes/origin/main" in workflow
    assert (
        'git merge-base --is-ancestor \\\n            "${GITHUB_SHA}" "refs/remotes/origin/main"'
    ) in workflow
    assert '[[ "$(git rev-parse HEAD)" == "${GITHUB_SHA}" ]]' in workflow
    assert "needs: [quality-gates, artifact-gates]" in workflow
    assert '--tag "${RELEASE_TAG}" \\' in workflow
    assert "--verify-local-tags" in workflow
    assert "python -m pytest" in workflow
    assert "python scripts/normalize_sdist.py dist" in workflow
    assert "python scripts/validate_release.py --sdist dist" in workflow
    assert "python -m twine check --strict dist/*" in workflow
    assert "os: [ubuntu-latest, macos-latest, windows-latest]" in workflow
    assert "artifact: [wheel, sdist]" in workflow


def test_no_workflow_reacts_to_a_tag_push() -> None:
    for path in sorted((REPOSITORY_ROOT / ".github/workflows").glob("*.yml")):
        workflow = path.read_text(encoding="utf-8")
        assert "tags:" not in workflow, path.name


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


def test_hook_events_are_the_protocol_events_and_a_misspelt_one_is_refused() -> None:
    # Every event the viewer sends a notebook takes a hook, except an agent
    # app's own (delivered to its port), plus the catch-all ``message``:
    # cellucid-r's cellucid_on() pins the same names to the same fixture
    # (test-viewer-protocol.R). A misspelt event used to be taken and never fired.
    protocol = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))
    events = {*protocol["events"], *protocol["binaryEvents"]}
    assert set(HOOK_EVENTS) == {event for event in events if not event.startswith("agent_")} | {"message"}

    hooks = HookRegistry()
    for call in (
        lambda: hooks.register("selections", print),
        lambda: hooks.on("selections"),
        lambda: hooks.unregister("selections", print),
        lambda: hooks.clear("selections"),
        lambda: hooks.register("agent_result", print),
    ):
        with pytest.raises(ValueError, match="^event must be one of 'action_result', .*; got '(selections|agent_result)'.$"):
            call()
    viewer = object.__new__(AnnDataViewer)
    for event in ("selections", "message"):
        with pytest.raises(ValueError, match=f"^event must be one of 'action_result', .*; got '{event}'.$"):
            viewer.wait_for_event(event, timeout=0)


def test_send_message_after_stop_is_an_explicit_error() -> None:
    # Before the viewer connects a command waits for it; after the viewer
    # stopped there is nothing left to wait for.
    viewer = object.__new__(AnnDataViewer)
    viewer._message_routing_registered = False
    with pytest.raises(RuntimeError, match="stopped"):
        viewer.send_message({"type": "ping", "requestId": "request-1"})


def test_starting_a_running_server_is_an_explicit_error() -> None:
    server = object.__new__(AnnDataServer)
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

    with pytest.raises(ValueError, match=r"\.h5ad|is not a Zarr store"):
        AnnDataServer(
            path,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Exact loader",
            dataset_id="exact-loader",
        )
    assert calls == []
