"""Exact contract for the two ways `cellucid.com` may be written.

Cellucid publishes `www.cellucid.com` (see `cellucid/CNAME` and the app's own
`<link rel="canonical">`), so every URL a reader is meant to open uses the
`www` host. The bare apex `https://cellucid.com` is reserved for identifiers
that must never change once published — JSON Schema `$id` values, the SVG
export namespace, and the provenance strings embedded in exported figures —
plus the CORS origin patterns the Python server matches.

Before this contract existed the split was enforced in one direction only: a
test pinned the documentation header link at the apex and nothing stopped a
frozen identifier from being "normalised" to `www`. This module enforces both
directions:

* every bare-apex occurrence must be a declared identifier or origin pattern,
  so a new navigational link cannot use the apex; and
* every frozen identifier is positively pinned at the apex, so it cannot be
  rewritten to `www`.

The prose statement of the rule lives in `docs/contributing/web_app.md`.
"""

from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"

CANONICAL_APP_ORIGIN = "https://www.cellucid.com"

# `https://www.cellucid.com` does not contain `https://cellucid.com`, so a plain
# substring search finds bare-apex occurrences and nothing else.
BARE_APEX = "https://cellucid.com"

SWEPT_SUFFIXES = frozenset({".md", ".py", ".rst", ".css", ".html", ".txt", ".cff"})
SWEPT_ROOTS = ("docs", "src")
SWEPT_FILES = ("README.md", "CONTRIBUTING.md", "CHANGELOG.md", "CITATION.cff")
EXCLUDED_PARTS = frozenset({"_build", "__pycache__", ".pytest_cache", ".mypy_cache"})

# The page that documents the convention necessarily quotes both forms.
CONVENTION_PAGE = "docs/contributing/web_app.md"

# Frozen identifiers: exact lines that must keep the bare apex. Rewriting any of
# them to `www` invalidates already-published references (a schema `$id` is
# compared as a string) or makes figures exported by two builds disagree about
# their own origin.
FROZEN_APEX_IDENTIFIERS: dict[str, tuple[str, ...]] = {
    "docs/user_guide/web_app/k_figure_export/05_metadata_and_provenance.md": (
        "- `Website`: `https://cellucid.com`",
        '  "generator": "https://cellucid.com",',
        '  "exporter": { "name": "Cellucid", "website": "https://cellucid.com" },',
    ),
    CONVENTION_PAGE: (
        "- `https://cellucid.com/contracts/community-annotation/user-v1.schema.json`",
        "- `https://cellucid.com/contracts/community-annotation/config-v1.schema.json`",
        "- `https://cellucid.com/contracts/community-annotation/merges-v1.schema.json`",
        "- `https://cellucid.com/ns#` (SVG export namespace)",
    ),
}

# `CORSMixin._get_allowed_origin` accepts the apex and any `.cellucid.com`
# subdomain over HTTPS. These lines describe that matcher, so they are neither
# links nor identifiers.
ORIGIN_PATTERN_LINE = "- `https://cellucid.com` and `https://*.cellucid.com`"
ORIGIN_PATTERN_PAGES = (
    "docs/user_guide/python_package/e_jupyter_hooks/"
    "13_security_cors_origins_and_mixed_content.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "06_configuration_env_vars_and_logging.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "17_security_privacy_and_networking.md",
)


def _swept_paths() -> list[Path]:
    paths: list[Path] = []
    for name in SWEPT_FILES:
        candidate = REPOSITORY_ROOT / name
        if candidate.is_file():
            paths.append(candidate)
    for root in SWEPT_ROOTS:
        for candidate in sorted((REPOSITORY_ROOT / root).rglob("*")):
            if not candidate.is_file():
                continue
            if candidate.suffix not in SWEPT_SUFFIXES:
                continue
            if EXCLUDED_PARTS.intersection(candidate.parts):
                continue
            paths.append(candidate)
    return paths


def _bare_apex_occurrences() -> list[tuple[str, int, str]]:
    found: list[tuple[str, int, str]] = []
    for path in _swept_paths():
        relative = path.relative_to(REPOSITORY_ROOT).as_posix()
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if BARE_APEX in line:
                found.append((relative, number, line))
    return found


def test_every_bare_apex_url_is_a_declared_identifier_or_origin_pattern() -> None:
    unclassified: list[str] = []
    for relative, number, line in _bare_apex_occurrences():
        if relative == CONVENTION_PAGE:
            continue
        frozen = {text.strip() for text in FROZEN_APEX_IDENTIFIERS.get(relative, ())}
        if line.strip() in frozen:
            continue
        if relative in ORIGIN_PATTERN_PAGES and line.strip() == ORIGIN_PATTERN_LINE:
            continue
        unclassified.append(f"{relative}:{number}: {line.strip()}")

    assert not unclassified, (
        "bare-apex `https://cellucid.com` is reserved for frozen identifiers and "
        "CORS origin patterns; navigational URLs must use "
        f"{CANONICAL_APP_ORIGIN}. Unclassified occurrences:\n"
        + "\n".join(unclassified)
    )


def test_frozen_identifiers_are_still_written_at_the_bare_apex() -> None:
    for relative, required_lines in FROZEN_APEX_IDENTIFIERS.items():
        contents = (REPOSITORY_ROOT / relative).read_text(encoding="utf-8")
        for required in required_lines:
            assert required in contents, (
                f"{relative} no longer declares the frozen identifier {required!r}. "
                "These are names, not addresses: adding the `www` host to a "
                "published schema $id or to exported figure provenance is a "
                "breaking change, not a normalisation."
            )


def test_navigational_urls_use_the_canonical_www_host() -> None:
    conf = (DOCS_ROOT / "conf.py").read_text(encoding="utf-8")
    assert f'"url": "{CANONICAL_APP_ORIGIN}"' in conf

    index = (DOCS_ROOT / "index.md").read_text(encoding="utf-8")
    assert f"```{{button-link}} {CANONICAL_APP_ORIGIN}" in index

    server_base = (REPOSITORY_ROOT / "src/cellucid/_server_base.py").read_text(
        encoding="utf-8"
    )
    assert f'CELLUCID_WEB_URL = "{CANONICAL_APP_ORIGIN}"' in server_base

    contributing = (REPOSITORY_ROOT / CONVENTION_PAGE).read_text(encoding="utf-8")
    assert f"hosted at `{CANONICAL_APP_ORIGIN}`" in contributing

    readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
    assert CANONICAL_APP_ORIGIN in readme


def test_the_convention_page_states_which_form_is_canonical_for_which_purpose() -> None:
    contributing = (REPOSITORY_ROOT / CONVENTION_PAGE).read_text(encoding="utf-8")
    assert "## Cross-repository conventions" in contributing
    assert "Canonical URL forms" in contributing
    assert "frozen identifier" in contributing
    assert "scripts/sweep_versions.py" in contributing
    assert "CELLUCID_VERSION" in contributing


def test_citation_files_point_at_a_reachable_canonical_landing_page() -> None:
    citation = (REPOSITORY_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert 'url: "https://cellucid.readthedocs.io"' in citation
    assert BARE_APEX not in citation


def test_the_release_version_marker_convention_is_applied_to_this_repository() -> None:
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'(?m)^version = "([^"]+)"\s+# CELLUCID_VERSION$', pyproject)
    assert match is not None
    version = match.group(1)

    citation = (REPOSITORY_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert f"version: {version}  # CELLUCID_VERSION" in citation, (
        "CITATION.cff must carry the CELLUCID_VERSION marker so a bump sweep "
        "lists it, matching cellucid-r/CITATION.cff"
    )
