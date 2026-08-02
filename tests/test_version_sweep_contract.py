"""Exact contract for the cross-repository version sweep.

`scripts/sweep_versions.py` is the one command that answers "do all six
repositories declare the same release?". It cannot be a test — the repositories
are separate checkouts with separate CI, and `cellucid-python`'s suite can only
see its own tree — so this module guards the parts that *are* verifiable here:
the declaration table's shape, its coverage of this repository, the marker
discipline that keeps a bump sweep able to find every declaration, and the
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
import subprocess
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


def test_a_marker_carried_off_the_declaring_line_always_states_why() -> None:
    for declaration in sweep.DECLARATIONS:
        if declaration.marker_alternative:
            assert declaration.kind == sweep.RELEASE, (
                f"{declaration.repository}/{declaration.relative_path} is not a "
                "release declaration, so it needs no marker and no exemption"
            )
            assert declaration.marker_reason, (
                f"{declaration.repository}/{declaration.relative_path} does not "
                "carry its marker on the declaring line, so the reason the "
                "ordinary rule cannot apply must be written down"
            )
        else:
            assert not declaration.marker_reason, declaration


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


def _own_declarations() -> list:
    return [
        declaration
        for declaration in sweep.DECLARATIONS
        if declaration.repository == "cellucid-python"
    ]


def test_this_repository_is_fully_declared_and_internally_consistent() -> None:
    version = _package_version()
    own = _own_declarations()
    assert own, "the sweep must cover cellucid-python"
    for declaration in own:
        path = REPOSITORY_ROOT / declaration.relative_path
        assert path.is_file(), f"{declaration.relative_path} does not exist"
        text = path.read_text(encoding="utf-8")
        matches = re.findall(declaration.pattern, text)
        assert matches, f"{declaration.relative_path} pattern matched nothing"
        if declaration.kind != sweep.RELEASE:
            # An artifact stamp records the release that produced it, so it is
            # allowed to lag. It still has to name some release.
            continue
        found = matches[:1] if declaration.first_only else matches
        assert set(found) == {version}, (
            f"{declaration.relative_path} declares {sorted(set(found))}, "
            f"not {version}"
        )


def test_every_release_declaration_here_is_reachable_from_a_bump_grep() -> None:
    """`grep -rn CELLUCID_VERSION` must list every line a bump has to edit."""
    for declaration in _own_declarations():
        if declaration.kind != sweep.RELEASE:
            continue
        path = REPOSITORY_ROOT / declaration.relative_path
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(declaration.pattern, text):
            line = sweep._line_of(text, match.start(1))
            assert sweep._is_marked(text, line, declaration), (
                f"{declaration.relative_path}:{line} declares the release but "
                "carries no CELLUCID_VERSION marker, so a bump sweep inside "
                "this checkout will not list it"
            )
            if declaration.first_only:
                break


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


def _checkout(root: Path) -> None:
    """Make `root` a Git checkout with everything written so far staged.

    The sweep reads the Git index, not the working tree, so a synthetic
    workspace has to be a set of checkouts to be readable at all.
    """
    root.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "-c", "init.defaultBranch=main", "init", "--quiet", str(root)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["git", "-C", str(root), "add", "--all"], check=True, capture_output=True
    )


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
    for name in sweep.REPOSITORIES:
        _checkout(workspace / name)

    assert sweep.release_version(workspace) == "1.2.3"
    result = sweep.collect(workspace)
    assert result.missing_repositories == []
    assert result.unindexed_repositories == []

    by_kind: dict[str, set[str]] = {}
    for finding in result.findings:
        by_kind.setdefault(finding.declaration.kind, set()).add(finding.value)

    # The stale R citation is drift; the npm placeholder and the old export are
    # not, and must never be counted as such.
    assert by_kind[sweep.RELEASE] == {"1.2.3", "1.2.2"}
    assert by_kind[sweep.INDEPENDENT] == {"0.0.1"}
    assert by_kind[sweep.ARTIFACT] == {"0.0.9"}
    assert sweep.report(workspace) == 1


def _marker_workspace(workspace: Path, citation_marker: str) -> None:
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _write(workspace / "cellucid-r" / "CITATION.cff", f"version: 1.2.3{citation_marker}\n")
    for name in sweep.REPOSITORIES:
        _checkout(workspace / name)


def test_a_release_declaration_without_a_marker_fails_the_sweep(tmp_path, capsys) -> None:
    workspace = tmp_path / "workspace"
    _marker_workspace(workspace, citation_marker="")

    citation = [
        finding
        for finding in sweep.collect(workspace).findings
        if finding.path == "cellucid-r/CITATION.cff"
    ]
    assert len(citation) == 1
    # The value agrees with the release, so nothing but the missing marker can
    # be the reason this run fails.
    assert citation[0].value == "1.2.3"
    assert not citation[0].marked

    assert sweep.report(workspace) == 1
    printed = capsys.readouterr().out
    assert "UNMARKED RELEASE DECLARATIONS" in printed
    assert "cellucid-r/CITATION.cff:1" in printed
    assert "1 release declarations carry no CELLUCID_VERSION marker" in printed


def test_the_same_declaration_passes_the_marker_check_once_it_is_marked(
    tmp_path, capsys
) -> None:
    workspace = tmp_path / "workspace"
    _marker_workspace(workspace, citation_marker="  # CELLUCID_VERSION")

    citation = [
        finding
        for finding in sweep.collect(workspace).findings
        if finding.path == "cellucid-r/CITATION.cff"
    ]
    assert len(citation) == 1
    assert citation[0].marked

    sweep.report(workspace)
    assert "carry no CELLUCID_VERSION marker" not in capsys.readouterr().out


def test_a_marker_on_the_following_line_counts_only_where_it_is_declared(
    tmp_path,
) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    # CHANGELOG.md declares `marker_alternative=NEXT_LINE`; CITATION.cff does
    # not, so the identical layout must not satisfy it.
    _write(
        workspace / "cellucid-python" / "CHANGELOG.md",
        "## [1.2.3]\n<!-- CELLUCID_VERSION -->\n",
    )
    _write(
        workspace / "cellucid-r" / "CITATION.cff",
        "version: 1.2.3\n# CELLUCID_VERSION\n",
    )
    for name in sweep.REPOSITORIES:
        _checkout(workspace / name)

    marked = {
        finding.path: finding.marked for finding in sweep.collect(workspace).findings
    }
    assert marked["cellucid-python/CHANGELOG.md"] is True
    assert marked["cellucid-r/CITATION.cff"] is False


def test_the_sweep_reads_the_git_index_and_not_the_working_tree(tmp_path) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _write(workspace / "cellucid-datasets" / ".gitignore", "/exports/_*/\n")
    _write(
        workspace / "cellucid-datasets" / "exports" / "suo" / "dataset_identity.json",
        '{"cellucid_data_version": "1.2.3"}\n',
    )
    for name in sweep.REPOSITORIES:
        _checkout(workspace / name)

    # A local scratch export, ignored and never committed. A working-tree walk
    # would read it and put a value nobody published into a release report.
    _write(
        workspace / "cellucid-datasets" / "exports" / "_test" / "dataset_identity.json",
        '{"cellucid_data_version": "0.0.1a2"}\n',
    )

    artifacts = {
        finding.path
        for finding in sweep.collect(workspace).findings
        if finding.declaration.kind == sweep.ARTIFACT
    }
    assert artifacts == {"cellucid-datasets/exports/suo/dataset_identity.json"}


def test_a_repository_without_a_git_index_can_not_be_swept(tmp_path, capsys) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _write(
        workspace / "cellucid-datasets" / "exports" / "suo" / "dataset_identity.json",
        '{"cellucid_data_version": "1.2.3"}\n',
    )
    for name in sweep.REPOSITORIES:
        if name == "cellucid-datasets":
            (workspace / name).mkdir(parents=True, exist_ok=True)
            continue
        _checkout(workspace / name)

    result = sweep.collect(workspace)
    assert result.unindexed_repositories == ["cellucid-datasets"]
    assert result.missing_repositories == []
    assert not [
        finding for finding in result.findings if finding.path.startswith("cellucid-datasets/")
    ]
    # A directory whose trackedness cannot be established is not a clean sweep.
    assert sweep.report(workspace) == 1
    assert "REPOSITORIES WITHOUT A GIT INDEX" in capsys.readouterr().out


def test_the_sweep_flags_a_partial_workspace_instead_of_passing_quietly(
    tmp_path,
) -> None:
    workspace = tmp_path / "workspace"
    _write(
        workspace / "cellucid-python" / "pyproject.toml",
        'version = "1.2.3"  # CELLUCID_VERSION\n',
    )
    _checkout(workspace / "cellucid-python")
    result = sweep.collect(workspace)
    assert set(result.missing_repositories) == set(sweep.REPOSITORIES) - {
        "cellucid-python"
    }
    # A workspace that cannot see five of six repositories must not look like a
    # clean release gate.
    assert sweep.report(workspace) == 1


def test_the_sweep_fails_when_it_can_not_read_the_release_version(tmp_path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    assert sweep.release_version(workspace) is None
    assert sweep.report(workspace) == 1
