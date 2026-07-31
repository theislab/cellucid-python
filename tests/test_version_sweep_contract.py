"""Exact contract for the cross-repository version sweep.

`scripts/sweep_versions.py` is the one command that answers "do all six
repositories declare the same release?". It cannot be a test — the repositories
are separate checkouts with separate CI, and `cellucid-python`'s suite can only
see its own tree — so this module guards the parts that *are* verifiable here:
the declaration table's shape, its coverage of this repository, and the
classification logic that keeps deliberate divergence from being reported as
drift.

The classification matters as much as the sweep. `cellucid/npm/package.json`
sits at `0.0.1` on purpose (it is a name reservation, documented in
`cellucid/npm/publishing.md`), and committed exports stamp the exporter release
that built them, so older exports legitimately record older versions. A sweep
that flagged either would be ignored within a release.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SWEEP_PATH = REPOSITORY_ROOT / "scripts" / "sweep_versions.py"


def _load_sweep() -> ModuleType:
    spec = importlib.util.spec_from_file_location("cellucid_sweep_versions", SWEEP_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


sweep = _load_sweep()


def _package_version() -> str:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'(?m)^version = "([^"]+)"\s+# CELLUCID_VERSION$', pyproject)
    assert match is not None
    return match.group(1)


def test_the_sweep_covers_every_cellucid_repository() -> None:
    assert sweep.REPOSITORIES == (
        "cellucid",
        "cellucid-python",
        "cellucid-r",
        "cellucid-datasets",
        "cellucid-demo-custom-datasets",
        "cellucid-annotation",
    )
    declared = {declaration.repository for declaration in sweep.DECLARATIONS}
    accounted = declared | set(sweep.NO_VERSION_BY_DESIGN)
    assert accounted == set(sweep.REPOSITORIES), (
        "every repository must either declare a version or record why it does "
        f"not; unaccounted: {set(sweep.REPOSITORIES) - accounted}"
    )
    assert set(sweep.NO_VERSION_BY_DESIGN) <= set(sweep.REPOSITORIES)
    assert all(reason for reason in sweep.NO_VERSION_BY_DESIGN.values())


def test_every_declaration_is_classified_and_independent_ones_are_justified() -> None:
    kinds = {sweep.RELEASE, sweep.INDEPENDENT, sweep.ARTIFACT}
    for declaration in sweep.DECLARATIONS:
        assert declaration.kind in kinds, declaration
        assert declaration.repository in sweep.REPOSITORIES, declaration
        if declaration.kind != sweep.RELEASE:
            assert declaration.note, (
                f"{declaration.repository}/{declaration.relative_path} is excluded "
                "from the agreement check, so the reason must be written down"
            )


def test_the_npm_placeholder_is_excluded_from_the_agreement_check() -> None:
    placeholder = [
        declaration
        for declaration in sweep.DECLARATIONS
        if declaration.repository == "cellucid"
        and declaration.relative_path == "npm/package.json"
    ]
    assert len(placeholder) == 1
    assert placeholder[0].kind == sweep.INDEPENDENT
    assert "publishing.md" in placeholder[0].note


def test_this_repository_is_fully_declared_and_internally_consistent() -> None:
    version = _package_version()
    own = [
        declaration
        for declaration in sweep.DECLARATIONS
        if declaration.repository == "cellucid-python"
    ]
    assert own, "the sweep must cover cellucid-python"
    for declaration in own:
        path = REPOSITORY_ROOT / declaration.relative_path
        assert path.is_file(), f"{declaration.relative_path} does not exist"
        text = path.read_text(encoding="utf-8")
        matches = re.findall(declaration.pattern, text)
        assert matches, f"{declaration.relative_path} pattern matched nothing"
        found = matches[:1] if declaration.first_only else matches
        assert set(found) == {version}, (
            f"{declaration.relative_path} declares {sorted(set(found))}, "
            f"not {version}"
        )


def test_changelog_headings_are_read_as_the_current_release_only() -> None:
    changelogs = [
        declaration
        for declaration in sweep.DECLARATIONS
        if declaration.relative_path in {"CHANGELOG.md", "NEWS.md"}
    ]
    assert changelogs
    for declaration in changelogs:
        assert declaration.first_only, (
            f"{declaration.repository}/{declaration.relative_path} lists past "
            "releases; only the topmost heading declares the current version"
        )


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_the_sweep_reports_drift_and_tolerates_deliberate_divergence(tmp_path) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _write(workspace / "cellucid" / "npm" / "package.json", '  "version": "0.0.1",\n')
    _write(
        workspace / "cellucid-r" / "CITATION.cff",
        "version: 1.2.2  # CELLUCID_VERSION\n",
    )
    _write(
        workspace / "cellucid-datasets" / "exports" / "old" / "dataset_identity.json",
        '{"cellucid_data_version": "0.0.9"}\n',
    )
    for name in ("cellucid-demo-custom-datasets", "cellucid-annotation"):
        (workspace / name).mkdir(parents=True, exist_ok=True)

    assert sweep.release_version(workspace) == "1.2.3"
    findings, missing, _ = sweep.collect(workspace)
    assert missing == []

    by_kind: dict[str, set[str]] = {}
    for finding in findings:
        by_kind.setdefault(finding.declaration.kind, set()).add(finding.value)

    # The stale R citation is drift; the npm placeholder and the old export are
    # not, and must never be counted as such.
    assert by_kind[sweep.RELEASE] == {"1.2.3", "1.2.2"}
    assert by_kind[sweep.INDEPENDENT] == {"0.0.1"}
    assert by_kind[sweep.ARTIFACT] == {"0.0.9"}
    assert sweep.report(workspace) == 1


def test_the_sweep_flags_a_partial_workspace_instead_of_passing_quietly(
    tmp_path,
) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _, missing, _ = sweep.collect(workspace)
    assert set(missing) == set(sweep.REPOSITORIES) - {"cellucid-python"}
    # A workspace that cannot see five of six repositories must not look like a
    # clean release gate.
    assert sweep.report(workspace) == 1


def test_the_sweep_fails_when_it_can_not_read_the_release_version(tmp_path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    assert sweep.release_version(workspace) is None
    assert sweep.report(workspace) == 1
