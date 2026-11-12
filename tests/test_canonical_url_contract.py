"""Exact contract for the two ways `cellucid.com` may be written.

Cellucid publishes `www.cellucid.com` (see `cellucid/CNAME` and the app's own
`<link rel="canonical">`), so every URL a reader is meant to open uses the
`www` host. The bare apex `https://cellucid.com` is reserved for identifiers
that must never change once published — JSON Schema `$id` values, the SVG
export namespace, and the provenance strings embedded in exported figures —
plus the CORS origin patterns the Python server matches.

Before this contract existed the split was enforced in one direction only: a
test pinned the documentation header link at the apex and nothing stopped a
frozen identifier from being "normalized" to `www`. This module enforces both
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
    # The apex is one of the two origins the hosted app is served from, which
    # the server's CORS policy names exactly.
    "src/cellucid/_server_base.py": (
        'HOSTED_APP_ORIGINS = frozenset({CELLUCID_WEB_URL, "https://cellucid.com"})',
    ),
    CONVENTION_PAGE: (
        "- `https://cellucid.com/contracts/community-annotation/user-v3.schema.json`",
        "- `https://cellucid.com/contracts/community-annotation/config-v1.schema.json`",
        "- `https://cellucid.com/contracts/community-annotation/merges-v3.schema.json`",
        "- `https://cellucid.com/ns#` (SVG export namespace)",
    ),
}

# `CORSMixin._get_allowed_origin` accepts the hosted app's two origins exactly
# (`HOSTED_APP_ORIGINS`). This line names them as origins, so it is neither a
# link nor an identifier.
ORIGIN_PATTERN_LINE = "- `https://www.cellucid.com` and `https://cellucid.com` (the hosted app's origins, exactly)"
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
                "breaking change, not a normalization."
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
    assert 'url: "https://docs.cellucid.com"' in citation
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


# -- Documentation host -------------------------------------------------------
# The documentation is served from the Read the Docs canonical domain
# docs.cellucid.com. The former `cellucid.readthedocs.io` host only redirects,
# so no link may name it, and the Sphinx configuration must declare the new host
# as the canonical base and publish a sitemap that robots.txt advertises.

DOCS_ORIGIN = "https://docs.cellucid.com"
DOCS_BASE_URL = f"{DOCS_ORIGIN}/en/latest/"
DOCS_SITEMAP_URL = f"{DOCS_BASE_URL}sitemap.xml"
RETIRED_DOCS_HOST = "cellucid.readthedocs.io"


def test_no_link_names_the_retired_read_the_docs_host() -> None:
    offenders = []
    for path in sorted(REPOSITORY_ROOT.rglob("*")):
        relative = path.relative_to(REPOSITORY_ROOT)
        if not path.is_file() or EXCLUDED_PARTS.intersection(relative.parts):
            continue
        # .claude holds local agent worktrees (other checkouts of this repository, git-excluded): not this tree's files.
        if relative.parts[0] in {".git", ".claude", "build", ".tmp", "data"} or relative == Path(__file__).relative_to(REPOSITORY_ROOT):
            continue
        if path.suffix not in SWEPT_SUFFIXES | {".toml", ".yaml", ".yml", ".ipynb", ".R", ".Rmd"}:
            continue
        if RETIRED_DOCS_HOST in path.read_text(encoding="utf-8", errors="ignore"):
            offenders.append(str(relative))
    assert offenders == []


def test_sphinx_declares_the_canonical_docs_host_and_publishes_a_sitemap() -> None:
    conf = (DOCS_ROOT / "conf.py").read_text(encoding="utf-8")
    assert f'html_baseurl = "{DOCS_BASE_URL}"' in conf
    assert 'sitemap_url_scheme = "{link}"' in conf
    assert '"sphinx_sitemap",' in conf
    assert 'html_extra_path = ["_extra"]' in conf

    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(r'"sphinx-sitemap==\d+\.\d+\.\d+",', pyproject)


def test_docs_pages_carry_a_description_a_link_card_and_a_unique_title() -> None:
    """Search results and link previews of every documentation page.

    sphinxext-opengraph gives each page a meta description and Open Graph tags
    on the canonical host, with the web app's 1200 x 630 card image (the URL
    cellucid's tests/seo-head-contract.test.mjs pins); a title suffix names the
    guide, so the Python, R and web app pages of one name stay distinct.
    """
    conf = (DOCS_ROOT / "conf.py").read_text(encoding="utf-8")
    assert '"sphinxext.opengraph",' in conf
    assert "ogp_site_url = html_baseurl" in conf
    assert 'ogp_image = "https://www.cellucid.com/assets/img/og-preview.jpg"' in conf
    assert '<meta name="twitter:card" content="summary_large_image" />' in conf
    assert 'ogp_social_cards = {"enable": False}' in conf
    for prefix in ("user_guide/web_app/", "user_guide/python_package/", "user_guide/r_package/"):
        assert f'("{prefix}", " — Cellucid ' in conf
    assert 'app.connect("html-page-context", _title_suffix, priority=600)' in conf
    # The reading card (Audience, Time, ...) leaves the doctree before the
    # extension reads it, so 172 descriptions are not one boilerplate.
    assert 'app.connect("html-page-context", _drop_reading_card, priority=400)' in conf
    layout = (DOCS_ROOT / "_templates" / "layout.html").read_text(encoding="utf-8")
    assert "<title>{{ title|striptags|e }}{{ cellucid_title_suffix|e }}</title>" in layout

    index = (DOCS_ROOT / "index.md").read_text(encoding="utf-8")
    front_matter = index.split("---\n")[1]
    assert "myst:\n  html_meta:\n    description: " in front_matter
    assert 'og:title: "Cellucid documentation"' in front_matter

    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(r'"sphinxext-opengraph==\d+\.\d+\.\d+",', pyproject)


def test_no_two_pages_of_one_guide_share_a_title() -> None:
    """A page's <title> is its first heading plus its guide's suffix (conf.py)."""
    import json

    guides = ("user_guide/web_app/", "user_guide/python_package/", "user_guide/r_package/")
    titles: dict[tuple[str, str], list[str]] = {}
    for path in sorted(DOCS_ROOT.rglob("*")):
        relative = path.relative_to(DOCS_ROOT).as_posix()
        if path.suffix not in {".md", ".ipynb"} or relative.startswith(("_", "user_guide/web_app/r_screenshot_checklist")):
            continue
        if path.suffix == ".md":
            text = path.read_text(encoding="utf-8")
        else:
            cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
            text = "\n".join("".join(cell["source"]) for cell in cells if cell["cell_type"] == "markdown")
        heading = re.search(r"(?m)^# (.+)$", text)
        if heading is None:
            continue
        guide = next((prefix for prefix in guides if relative.startswith(prefix)), "")
        titles.setdefault((guide, heading.group(1).strip()), []).append(relative)
    assert {key: pages for key, pages in titles.items() if len(pages) > 1} == {}


def test_docs_publish_llms_txt_with_the_summary_the_web_app_uses() -> None:
    conf = (DOCS_ROOT / "conf.py").read_text(encoding="utf-8")
    assert '"sphinx_llms_txt",' in conf
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(r'"sphinx-llms-txt==\d+\.\d+\.\d+",', pyproject)
    summary = re.search(r"llms_txt_summary = \(\n((?:    \".*\"\n)+)\)", conf)
    assert summary is not None
    text = "".join(re.findall(r'"(.*)"', summary.group(1)))
    # The first two sentences are the web app's llms.txt and JSON-LD wording.
    assert text.startswith(
        "Cellucid is a free, open-source web app for interactive single-cell data "
        "visualization. It renders millions of cells on the GPU in the browser"
    )
    assert f"{CANONICAL_APP_ORIGIN}/" in text


def test_docs_robots_txt_allows_crawling_and_advertises_the_sitemap() -> None:
    robots = (DOCS_ROOT / "_extra" / "robots.txt").read_text(encoding="utf-8")
    lines = [line.strip() for line in robots.splitlines() if line.strip()]
    assert lines == ["User-agent: *", "Disallow:", f"Sitemap: {DOCS_SITEMAP_URL}"]


def test_the_docs_landing_page_links_to_the_app() -> None:
    index = (DOCS_ROOT / "index.md").read_text(encoding="utf-8")
    assert f"{{button-link}} {CANONICAL_APP_ORIGIN}" in index
