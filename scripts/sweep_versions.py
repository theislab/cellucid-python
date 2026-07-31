#!/usr/bin/env python3
"""Report every version declaration across all six Cellucid repositories.

The `CELLUCID_VERSION` marker comment makes `grep -rn CELLUCID_VERSION` list the
edit set for a bump *inside one checkout*. It cannot answer the question a
release actually asks — do all six repositories agree? — because they are
separate checkouts with separate CI, so no single test suite can see them all.
This script is that answer: run it once from a workspace that holds the
repositories side by side and read the verdict.

    python cellucid-python/scripts/sweep_versions.py

It classifies every declaration it finds, because "all version literals must be
equal" is false here and a sweep that pretends otherwise is useless:

* ``release``     -- must equal the release version. Drift here is a defect.
* ``independent`` -- deliberately unrelated to the release version, with the
                     reason recorded next to it. The npm package is a name
                     reservation pinned at 0.0.1 on purpose (see
                     ``cellucid/npm/publishing.md``); flagging it would train
                     maintainers to ignore the sweep.
* ``artifact``    -- the exporter release that produced a committed export.
                     Older exports legitimately record older versions; these are
                     printed for visibility and never counted as drift.

Repositories that are not checked out are reported as such rather than skipped
silently, so a partial workspace cannot produce a falsely clean verdict.

Exit status is 0 only when every repository was swept and every ``release``
declaration agrees. A partial workspace exits non-zero: an unswept repository
cannot be confirmed, and a release gate that cannot see a repository must not
look like a pass.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

RELEASE = "release"
INDEPENDENT = "independent"
ARTIFACT = "artifact"

REPOSITORIES = (
    "cellucid",
    "cellucid-python",
    "cellucid-r",
    "cellucid-datasets",
    "cellucid-demo-custom-datasets",
    "cellucid-annotation",
)

# Repositories that declare no release version at all, and why. Recorded so the
# sweep can prove the absence is intentional instead of leaving a silent gap.
NO_VERSION_BY_DESIGN = {
    "cellucid": (
        "continuously deployed static site; a deployment is identified by the "
        "build id in the app footer, not a released version (cellucid/CITATION.cff)"
    ),
    "cellucid-annotation": (
        "reference repository layout, not a released artifact; the annotation "
        "contract is versioned by its schema $id suffix (`-v1`), not by a "
        "package version"
    ),
}

# The release version is read from the one place that defines it.
VERSION_SOURCE = ("cellucid-python", "pyproject.toml")
VERSION_SOURCE_PATTERN = r'(?m)^version = "([^"]+)"\s+# CELLUCID_VERSION$'


@dataclass(frozen=True)
class Declaration:
    """One place a version literal is written down."""

    repository: str
    relative_path: str
    pattern: str
    kind: str
    note: str = ""
    glob: bool = False
    # Changelogs list every past release. Only the topmost heading declares the
    # current version; the history below it must stay as written.
    first_only: bool = False


@dataclass(frozen=True)
class Finding:
    declaration: Declaration
    path: str
    line: int
    value: str
    marked: bool


DECLARATIONS: tuple[Declaration, ...] = (
    # -- cellucid (web app) --------------------------------------------------
    Declaration(
        "cellucid",
        "npm/package.json",
        r'(?m)^  "version": "([^"]+)",$',
        INDEPENDENT,
        "npm name-reservation placeholder, deliberately unrelated to the "
        "package releases; documented in cellucid/npm/publishing.md",
    ),
    Declaration(
        "cellucid",
        "npm/publishing.md",
        r"independent of the Python package's version \(`([^`]+)` at time of writing\)",
        RELEASE,
        "cross-reference to the Python release inside the placeholder rationale",
    ),
    # -- cellucid-python -----------------------------------------------------
    Declaration("cellucid-python", "pyproject.toml", VERSION_SOURCE_PATTERN, RELEASE),
    Declaration("cellucid-python", "CITATION.cff", r"(?m)^version: (\S+)", RELEASE),
    Declaration(
        "cellucid-python",
        "CHANGELOG.md",
        r"(?m)^## \[([^\]]+)\]$",
        RELEASE,
        "current release heading",
        first_only=True,
    ),
    Declaration(
        "cellucid-python",
        "docs/conf.py",
        r'(?m)^\s*release = version = "([^"]+)"',
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "scripts/publishing/meta.yaml",
        r'(?m)^\{% set version = "([^"]+)" %\}',
        RELEASE,
        "conda-forge recipe",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/installation.md",
        r'pip install "cellucid==([^"]+)"',
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/h_developer_docs/04_build_install_and_packaging.md",
        r'`version = "([^"]+)"` \(beta\)',
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/h_developer_docs/14_release_process.md",
        r"`cellucid` is currently in beta \(`([^`]+)`\)",
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/c_data_preparation_api/"
        "09_output_format_specification_exports_directory.md",
        r"`cellucid_data_version` is `([^`]+)`",
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        r"Version `([^`]+)` from `packageVersion`",
        RELEASE,
        "R release, declared on the shared documentation site",
    ),
    # -- cellucid-r ----------------------------------------------------------
    Declaration(
        "cellucid-r",
        "DESCRIPTION",
        r"(?m)^Version: (\S+)$",
        RELEASE,
        "DCF cannot carry a comment; the marker is the "
        "Config/cellucid/version-marker field",
    ),
    Declaration(
        "cellucid-r",
        "NEWS.md",
        r"(?m)^# cellucid (\S+)",
        RELEASE,
        "current release heading",
        first_only=True,
    ),
    Declaration("cellucid-r", "CITATION.cff", r"(?m)^version: (\S+)", RELEASE),
    Declaration(
        "cellucid-r",
        "README.md",
        r"Active package version — (\d+\.\d+\.\d+[0-9A-Za-z._-]*)",
        RELEASE,
    ),
    Declaration(
        "cellucid-r",
        "vignettes/installation.Rmd",
        r"source and documentation version is (\d+\.\d+\.\d+[0-9A-Za-z_-]*)",
        RELEASE,
    ),
    # -- cellucid-datasets ---------------------------------------------------
    Declaration(
        "cellucid-datasets",
        "sources/pancreas.json",
        r'(?s)"cellucid_producer": \{.*?"version": "([^"]+)"',
        ARTIFACT,
        "exporter release recorded in the Pancreas build provenance",
    ),
    Declaration(
        "cellucid-datasets",
        "tests/catalog-contract.test.mjs",
        r"(?s)'https://github\.com/theislab/cellucid-python',\s*\n\s*version: '([^']+)'",
        ARTIFACT,
        "test pin of the recorded exporter release",
    ),
    Declaration(
        "cellucid-datasets",
        "exports/*/dataset_identity.json",
        r'"cellucid_data_version": "([^"]+)"',
        ARTIFACT,
        "stamped when the export was built; older exports keep older values",
        glob=True,
    ),
    # -- cellucid-demo-custom-datasets ---------------------------------------
    Declaration(
        "cellucid-demo-custom-datasets",
        "generate_datasets.py",
        r'(?m)^CELLUCID_RELEASE = "([^"]+)"$',
        RELEASE,
    ),
    Declaration(
        "cellucid-demo-custom-datasets",
        "generate_datasets.py",
        r'(?m)^#\s+"cellucid==([^"]+)",$',
        RELEASE,
        "PEP 723 inline script dependency",
    ),
    Declaration(
        "cellucid-demo-custom-datasets",
        "README.md",
        r'python -m pip install "cellucid==([^"]+)"',
        RELEASE,
    ),
    Declaration(
        "cellucid-demo-custom-datasets",
        "README.md",
        r"`(\d+\.\d+\.\d+[0-9A-Za-z._-]*)` is the release these example exports "
        r"were built with",
        RELEASE,
    ),
    Declaration(
        "cellucid-demo-custom-datasets",
        "exports/*/dataset_identity.json",
        r'"cellucid_data_version": "([^"]+)"',
        ARTIFACT,
        "stamped when the export was built",
        glob=True,
    ),
)


def default_workspace() -> Path:
    """The directory that holds every Cellucid checkout side by side."""
    return Path(__file__).resolve().parents[2]


def _resolve(workspace: Path, declaration: Declaration) -> list[Path]:
    root = workspace / declaration.repository
    if declaration.glob:
        return sorted(root.glob(declaration.relative_path))
    candidate = root / declaration.relative_path
    return [candidate] if candidate.is_file() else []


def _line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _is_marked(text: str, line_number: int, declaration: Declaration) -> bool:
    """Whether a bump sweep (`grep -rn CELLUCID_VERSION`) would list this line."""
    lines = text.splitlines()
    if line_number <= len(lines) and "CELLUCID_VERSION" in lines[line_number - 1]:
        return True
    # DESCRIPTION is DCF and cannot carry a comment, so cellucid-r declares the
    # convention with a dedicated field instead.
    return "Config/cellucid/version-marker: CELLUCID_VERSION" in text


def collect(workspace: Path) -> tuple[list[Finding], list[str], list[Declaration]]:
    """Return (findings, missing repositories, declarations with no match)."""
    findings: list[Finding] = []
    missing_repositories = [
        name for name in REPOSITORIES if not (workspace / name).is_dir()
    ]
    unmatched: list[Declaration] = []

    for declaration in DECLARATIONS:
        if declaration.repository in missing_repositories:
            continue
        paths = _resolve(workspace, declaration)
        matched_any = False
        for path in paths:
            text = path.read_text(encoding="utf-8")
            for match in re.finditer(declaration.pattern, text):
                matched_any = True
                line = _line_of(text, match.start(1))
                findings.append(
                    Finding(
                        declaration=declaration,
                        path=path.relative_to(workspace).as_posix(),
                        line=line,
                        value=match.group(1),
                        marked=_is_marked(text, line, declaration),
                    )
                )
                if declaration.first_only:
                    break
        if not matched_any:
            unmatched.append(declaration)

    return findings, missing_repositories, unmatched


def release_version(workspace: Path) -> str | None:
    path = workspace / VERSION_SOURCE[0] / VERSION_SOURCE[1]
    if not path.is_file():
        return None
    match = re.search(VERSION_SOURCE_PATTERN, path.read_text(encoding="utf-8"))
    return match.group(1) if match else None


def _grouped(findings: list[Finding], kind: str) -> dict[str, list[Finding]]:
    grouped: dict[str, list[Finding]] = {}
    for finding in findings:
        if finding.declaration.kind == kind:
            grouped.setdefault(finding.value, []).append(finding)
    return dict(sorted(grouped.items()))


def report(workspace: Path) -> int:
    findings, missing_repositories, unmatched = collect(workspace)
    expected = release_version(workspace)
    out: list[str] = []

    out.append(f"Cellucid version sweep — workspace: {workspace}")
    out.append(f"Repositories swept: {len(REPOSITORIES) - len(missing_repositories)}"
               f" of {len(REPOSITORIES)}")
    source = "/".join(VERSION_SOURCE)
    out.append(f"Release version (source of truth: {source}): {expected or 'UNKNOWN'}")
    out.append("")

    release_groups = _grouped(findings, RELEASE)
    total_release = sum(len(group) for group in release_groups.values())
    out.append(f"RELEASE DECLARATIONS — must all agree ({total_release} found)")
    for value, group in release_groups.items():
        flag = "" if value == expected else "   <-- DRIFT"
        out.append(f"  {value}{flag}")
        for finding in group:
            marker = "marker" if finding.marked else "NO MARKER"
            note = f"  ({finding.declaration.note})" if finding.declaration.note else ""
            out.append(f"    {finding.path}:{finding.line}  [{marker}]{note}")
    out.append("")

    independent_groups = _grouped(findings, INDEPENDENT)
    out.append("INTENTIONALLY INDEPENDENT — excluded from the agreement check")
    for value, group in independent_groups.items():
        for finding in group:
            out.append(f"  {value}  {finding.path}:{finding.line}")
            out.append(f"        {finding.declaration.note}")
    if not independent_groups:
        out.append("  (none)")
    out.append("")

    artifact_groups = _grouped(findings, ARTIFACT)
    out.append(
        "ARTIFACT STAMPS — exporter release recorded in a committed artifact; "
        "older values are expected"
    )
    for value, group in artifact_groups.items():
        out.append(f"  {value}")
        for finding in group:
            out.append(f"    {finding.path}:{finding.line}")
    if not artifact_groups:
        out.append("  (none)")
    out.append("")

    out.append("NO RELEASE VERSION BY DESIGN")
    for name, reason in NO_VERSION_BY_DESIGN.items():
        out.append(f"  {name}")
        out.append(f"        {reason}")
    out.append("")

    unmarked = [
        finding
        for finding in findings
        if finding.declaration.kind == RELEASE and not finding.marked
    ]
    out.append(
        "MISSING CELLUCID_VERSION MARKER — a bump sweep inside the owning "
        f"checkout will not list these ({len(unmarked)} of {total_release})"
    )
    for finding in unmarked:
        out.append(f"  {finding.path}:{finding.line}")
    if not unmarked:
        out.append("  (none)")
    out.append("")

    problems = 0
    if missing_repositories:
        problems += len(missing_repositories)
        out.append("REPOSITORIES NOT CHECKED OUT — not swept, verdict is partial")
        for name in missing_repositories:
            out.append(f"  {name}")
        out.append("")
    if unmatched:
        problems += len(unmatched)
        out.append("DECLARATIONS THAT MATCHED NOTHING — the sweep is out of date")
        for declaration in unmatched:
            out.append(f"  {declaration.repository}/{declaration.relative_path}")
        out.append("")

    drifted = [value for value in release_groups if value != expected]
    if expected is None:
        problems += 1
        out.append("VERDICT: FAIL — could not read the release version")
    elif drifted:
        problems += 1
        out.append(
            f"VERDICT: FAIL — release declarations disagree: "
            f"{', '.join(sorted(set(drifted) | {expected}))}"
        )
    elif missing_repositories:
        out.append(
            f"VERDICT: PARTIAL — {total_release} release declarations agree at "
            f"{expected}, but {len(missing_repositories)} repositories were not swept"
        )
    else:
        out.append(
            f"VERDICT: OK — all {total_release} release declarations across "
            f"{len(REPOSITORIES)} repositories agree at {expected}"
        )

    print("\n".join(out))
    return 1 if problems else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--workspace",
        type=Path,
        default=default_workspace(),
        help="directory containing the Cellucid repositories side by side",
    )
    arguments = parser.parse_args()
    return report(arguments.workspace.resolve())


if __name__ == "__main__":
    sys.exit(main())
