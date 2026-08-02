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
* ``artifact``    -- the exporter release that produced a committed export, or
                     prose that describes one. Older exports legitimately record
                     older values; these are printed for visibility and never
                     counted as drift.

Every ``release`` declaration must also be reachable from inside its own
checkout, so each one carries the marker on its declaring line. A line that
cannot carry a comment — a DCF field, a JSON example, a transcript of a tool's
output — records how it carries the marker instead and why, in
``marker_alternative``/``marker_reason``. An unrecorded gap fails the sweep:
a declaration nobody can grep for is a declaration a bump will miss.

Only files in each repository's Git index are read. A glob that walked the
working tree would pick up ignored scratch directories, so the verdict would
depend on what happened to be lying around and a stale local export could inject
a version into a release report.

Repositories that are not checked out are reported as such rather than skipped
silently, so a partial workspace cannot produce a falsely clean verdict.

Exit status is 0 only when every repository was swept, every ``release``
declaration agrees, and every one of them is marked. A partial workspace exits
non-zero: an unswept repository cannot be confirmed, and a release gate that
cannot see a repository must not look like a pass.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

RELEASE = "release"
INDEPENDENT = "independent"
ARTIFACT = "artifact"

MARKER = "CELLUCID_VERSION"

# `marker_alternative` sentinel: the marker comment sits on the line after the
# declaration instead of on it.
NEXT_LINE = "<next-line>"

# One spelling of "a release version" for every pattern that has no surrounding
# delimiter to bound the match. The suffix must end alphanumerically, so a
# version quoted at the end of a sentence yields `0.9.1` and not `0.9.1.`.
VERSION = r"(\d+\.\d+\.\d+(?:[0-9A-Za-z._-]*[0-9A-Za-z])?)"

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
    # How this declaration carries its `CELLUCID_VERSION` marker when its own
    # line cannot. `NEXT_LINE` puts the comment on the following line; any other
    # value is a literal that must appear somewhere in the file. Empty means the
    # ordinary rule: the marker is on the declaring line.
    marker_alternative: str = ""
    # Why the ordinary rule cannot apply. Required whenever `marker_alternative`
    # is set, so an exemption is always an argument and never a convenience.
    marker_reason: str = ""


@dataclass(frozen=True)
class Finding:
    declaration: Declaration
    path: str
    line: int
    value: str
    marked: bool


@dataclass
class Sweep:
    """Everything one pass over a workspace learned."""

    findings: list[Finding] = field(default_factory=list)
    # Present in `REPOSITORIES` but not checked out here.
    missing_repositories: list[str] = field(default_factory=list)
    # Checked out, but with no readable Git index, so "tracked" is unanswerable
    # and nothing in them can be reported reproducibly.
    unindexed_repositories: list[str] = field(default_factory=list)
    # Declarations whose pattern found nothing in a tracked file.
    unmatched: list[Declaration] = field(default_factory=list)


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
        marker_alternative=NEXT_LINE,
        marker_reason=(
            "the release validator and the docs contract both require the exact "
            "line `## [<version>]`, so the marker goes on the line below it"
        ),
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
        marker_alternative=NEXT_LINE,
        marker_reason=(
            "the recipe validator requires the exact line "
            "`{% set version = \"<version>\" %}`, so the YAML comment goes below it"
        ),
    ),
    Declaration(
        "cellucid-python",
        "scripts/publishing/PUBLISHING.md",
        rf"For release {VERSION}, validate the candidate with:",
        RELEASE,
        "release-runbook heading",
    ),
    Declaration(
        "cellucid-python",
        "scripts/publishing/PUBLISHING.md",
        rf"python scripts/validate_release\.py --tag v{VERSION}",
        RELEASE,
        "the tag a maintainer copies out of the runbook; a bump that misses it "
        "validates the previous release",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/installation.md",
        r'pip install "cellucid==([^"]+)"',
        RELEASE,
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/installation.md",
        rf"lists\n{VERSION}, the normal install and pin commands",
        RELEASE,
        "PyPI availability note",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/installation.md",
        rf"when validating {VERSION} before publication",
        RELEASE,
        "pre-publication source-checkout note",
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
        "docs/user_guide/python_package/h_developer_docs/14_release_process.md",
        rf"with a `v` prefix, such as `v{VERSION}`",
        RELEASE,
        "the tag-naming example",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/h_developer_docs/14_release_process.md",
        rf"Cellucid Python release contract is valid for v{VERSION}\.",
        RELEASE,
        "the exact line the release validator prints",
        marker_alternative=f"<!-- {MARKER} -->",
        marker_reason=(
            "a documented tool transcript has to reproduce the output byte for "
            "byte, and the validator prints no comment; the marker is on the "
            "versioning-policy line of the same page"
        ),
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/python_package/h_developer_docs/14_release_process.md",
        rf"checks out complete Git history, validates `v{VERSION}`",
        RELEASE,
        "the tag the publishing workflow validates",
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
        "docs/user_guide/python_package/c_data_preparation_api/"
        "09_output_format_specification_exports_directory.md",
        r'"cellucid_data_version": "([^"]+)"',
        RELEASE,
        "the worked schema example the prose above it pins",
        marker_alternative=f"<!-- {MARKER} -->",
        marker_reason=(
            "JSON has no comment syntax and the block is validated against the "
            "reader's key set; the prose line that pins this value is marked"
        ),
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        r"Version `([^`]+)` from `packageVersion`",
        RELEASE,
        "R release, declared on the shared documentation site",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        rf"install Cellucid for R, verify version {VERSION}",
        RELEASE,
        "page goal",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        rf"### From CRAN when version {VERSION} is listed",
        RELEASE,
        "CRAN install section heading",
        marker_alternative=NEXT_LINE,
        marker_reason=(
            "an inline comment inside a heading enters the heading text and the "
            "anchor MyST generates from it, so the marker goes on the line below"
        ),
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        rf"When the CRAN package index lists `cellucid` {VERSION}, install it with",
        RELEASE,
        "the CRAN install instruction",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/installation.md",
        rf"If the index does not list version {VERSION} yet",
        RELEASE,
        "the source-install instruction for when CRAN does not list it yet",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/r_package/index.md",
        rf"Install version {VERSION} from the source currently available",
        RELEASE,
        "installation card on the R guide home page",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/web_app/b_data_loading/04_server_tutorial.md",
        rf"(?m)^cellucid v{VERSION}$",
        ARTIFACT,
        "the version banner inside the recorded `cellucid serve` transcripts; it "
        "states what that run printed, so it changes when the transcript is "
        "re-recorded, not at a bump. These were screenshots until a capture was "
        "found to have baked a local path into the image pixels — terminal "
        "output now ships as text, which is also why this pattern is greppable",
    ),
    Declaration(
        "cellucid-python",
        "docs/user_guide/web_app/b_data_loading/10_standard_pancreas_dataset.md",
        rf"the exact Cellucid {VERSION} producer-source digest",
        ARTIFACT,
        "describes the exporter release recorded in "
        "cellucid-datasets/sources/pancreas.json; it changes when the Pancreas "
        "export is rebuilt, not at a bump",
    ),
    # -- cellucid-r ----------------------------------------------------------
    Declaration(
        "cellucid-r",
        "DESCRIPTION",
        r"(?m)^Version: (\S+)$",
        RELEASE,
        "DCF cannot carry a comment; the marker is the "
        "Config/cellucid/version-marker field",
        marker_alternative=f"Config/cellucid/version-marker: {MARKER}",
        marker_reason="DCF has no comment syntax, so the convention is a field",
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
        rf"Active package version — {VERSION}",
        RELEASE,
    ),
    Declaration(
        "cellucid-r",
        "vignettes/installation.Rmd",
        rf"source and documentation version is {VERSION}",
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
        r'(?m)^CELLUCID_RELEASE = "([^"]+)"',
        RELEASE,
    ),
    Declaration(
        "cellucid-demo-custom-datasets",
        "generate_datasets.py",
        r'(?m)^#\s+"cellucid==([^"]+)",',
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
        rf"`{VERSION}` is the release these example exports were built with",
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


def tracked_files(root: Path) -> frozenset[str] | None:
    """Repository-relative paths in ``root``'s Git index, or None if unreadable.

    The sweep reads only these. A working-tree walk would also pick up ignored
    scratch directories — `cellucid-datasets/exports/_*/` is ignored precisely
    so local experiments stay local — and a release verdict that depends on
    uncommitted files is not reproducible from a clean checkout.
    """
    if not (root / ".git").exists():
        return None
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return frozenset(entry for entry in completed.stdout.split("\0") if entry)


def _resolve(
    workspace: Path, declaration: Declaration, tracked: frozenset[str]
) -> list[Path]:
    root = workspace / declaration.repository
    if declaration.glob:
        candidates = sorted(root.glob(declaration.relative_path))
    else:
        candidates = [root / declaration.relative_path]
    return [
        path
        for path in candidates
        if path.is_file() and path.relative_to(root).as_posix() in tracked
    ]


def _line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _is_marked(text: str, line_number: int, declaration: Declaration) -> bool:
    """Whether a bump sweep (`grep -rn CELLUCID_VERSION`) would list this line."""
    lines = text.splitlines()
    if line_number <= len(lines) and MARKER in lines[line_number - 1]:
        return True
    alternative = declaration.marker_alternative
    if not alternative:
        return False
    if alternative == NEXT_LINE:
        return line_number < len(lines) and MARKER in lines[line_number]
    return alternative in text


def collect(workspace: Path) -> Sweep:
    """Read every declaration the workspace can answer for."""
    sweep = Sweep()
    indexes: dict[str, frozenset[str]] = {}
    for name in REPOSITORIES:
        root = workspace / name
        if not root.is_dir():
            sweep.missing_repositories.append(name)
            continue
        tracked = tracked_files(root)
        if tracked is None:
            sweep.unindexed_repositories.append(name)
            continue
        indexes[name] = tracked

    for declaration in DECLARATIONS:
        tracked = indexes.get(declaration.repository)
        if tracked is None:
            continue
        matched_any = False
        for path in _resolve(workspace, declaration, tracked):
            text = path.read_text(encoding="utf-8")
            for match in re.finditer(declaration.pattern, text):
                matched_any = True
                line = _line_of(text, match.start(1))
                sweep.findings.append(
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
            sweep.unmatched.append(declaration)

    return sweep


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
    sweep = collect(workspace)
    findings = sweep.findings
    expected = release_version(workspace)
    out: list[str] = []
    unswept = len(sweep.missing_repositories) + len(sweep.unindexed_repositories)

    out.append(f"Cellucid version sweep — workspace: {workspace}")
    out.append(f"Repositories swept: {len(REPOSITORIES) - unswept}"
               f" of {len(REPOSITORIES)}")
    source = "/".join(VERSION_SOURCE)
    out.append(f"Release version (source of truth: {source}): {expected or 'UNKNOWN'}")
    out.append("Files read: Git-tracked only")
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

    exempted = [
        finding
        for finding in findings
        if finding.declaration.kind == RELEASE and finding.declaration.marker_alternative
    ]
    out.append(
        "MARKERS CARRIED ELSEWHERE — the declaring line cannot hold a comment "
        f"({len(exempted)} of {total_release})"
    )
    for finding in exempted:
        where = (
            "next line"
            if finding.declaration.marker_alternative == NEXT_LINE
            else finding.declaration.marker_alternative
        )
        out.append(f"  {finding.path}:{finding.line}  [{where}]")
        out.append(f"        {finding.declaration.marker_reason}")
    if not exempted:
        out.append("  (none)")
    out.append("")

    unmarked = [
        finding
        for finding in findings
        if finding.declaration.kind == RELEASE and not finding.marked
    ]
    out.append(
        "UNMARKED RELEASE DECLARATIONS — a bump sweep inside the owning "
        f"checkout will not list these ({len(unmarked)} of {total_release})"
    )
    for finding in unmarked:
        out.append(f"  {finding.path}:{finding.line}")
    if not unmarked:
        out.append("  (none)")
    out.append("")

    failures: list[str] = []
    if sweep.unindexed_repositories:
        failures.append(
            f"{len(sweep.unindexed_repositories)} repositories have no readable "
            "Git index"
        )
        out.append("REPOSITORIES WITHOUT A GIT INDEX — trackedness is unanswerable")
        for name in sweep.unindexed_repositories:
            out.append(f"  {name}")
        out.append("")
    if sweep.missing_repositories:
        out.append("REPOSITORIES NOT CHECKED OUT — not swept, verdict is partial")
        for name in sweep.missing_repositories:
            out.append(f"  {name}")
        out.append("")
    if sweep.unmatched:
        failures.append(f"{len(sweep.unmatched)} declarations matched nothing")
        out.append(
            "DECLARATIONS THAT MATCHED NOTHING — the sweep is out of date, or "
            "the file is not tracked"
        )
        for declaration in sweep.unmatched:
            out.append(f"  {declaration.repository}/{declaration.relative_path}")
        out.append("")

    drifted = [value for value in release_groups if value != expected]
    if expected is None:
        failures.append("could not read the release version")
    elif drifted:
        failures.append(
            "release declarations disagree: "
            f"{', '.join(sorted(set(drifted) | {expected}))}"
        )
    if unmarked:
        failures.append(
            f"{len(unmarked)} release declarations carry no {MARKER} marker"
        )

    if failures:
        out.append(f"VERDICT: FAIL — {'; '.join(failures)}")
    elif sweep.missing_repositories:
        out.append(
            f"VERDICT: PARTIAL — {total_release} release declarations agree at "
            f"{expected}, but {len(sweep.missing_repositories)} repositories "
            "were not swept"
        )
    else:
        out.append(
            f"VERDICT: OK — all {total_release} release declarations across "
            f"{len(REPOSITORIES)} repositories agree at {expected} and every "
            f"one of them is reachable from `grep -rn {MARKER}`"
        )

    print("\n".join(out))
    return 1 if failures or sweep.missing_repositories else 0


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
