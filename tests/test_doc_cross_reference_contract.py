"""Exact contract for how one documentation page references another.

Every internal page reference on the docs site is a `{doc}` role. A bare
backtick is inline literal text: Sphinx does not resolve it, does not check
that the target exists, and does not render a link. The docs build runs with
`-W`, so a `{doc}` role pointing at a missing page fails the build — but a
backticked page name bypasses that check entirely and ships as plain grey text.

That is not hypothetical. A reference to `k_sessions_sharing/index` survived in
the highlighting section long after the directory was renamed to
`l_sessions_sharing/`, and a reference to
`b_data_loading/04_troubleshooting_data_loading` pointed at the server tutorial
for as long as it existed, because nothing was ever asked to resolve either one.

This module fails when a bare inline literal is shaped like a page on this site.
It deliberately does not require the target to exist: a *dead* backticked
reference is the exact defect being prevented, so shape alone is the trigger.

The prose statement of the rule lives in
`docs/user_guide/python_package/h_developer_docs/15_docs_development_and_style_guide.md`.
"""

from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"

EXCLUDED_PARTS = frozenset(
    {"_build", "jupyter_execute", ".ipynb_checkpoints", "__pycache__"}
)

# The page that states the rule has to be able to quote page names as text.
CONVENTION_PAGE = (
    "user_guide/python_package/h_developer_docs/15_docs_development_and_style_guide.md"
)

# The screenshot checklist audits this site page by page, so it names pages as
# filenames on purpose: as abbreviations in inventory tables (`b/09`, `k/index`),
# and as records of pages that were deleted and therefore have no target to link
# to. A link is the wrong shape for both, so the chapter is exempt from this rule
# while every reader-facing chapter stays under it.
AUDIT_SECTION = "user_guide/web_app/r_screenshot_checklist/"


def _is_exempt(relative: str) -> bool:
    return relative == CONVENTION_PAGE or relative.startswith(AUDIT_SECTION)

# A page on this site is `index` or `NN_slug`, optionally behind a path, and
# optionally written with the source extension it has on disk. The extension
# form is the one that slips past a reader: `10_troubleshooting_sessions.md`
# looks like a filename rather than a link, so it was never written as a role
# and Sphinx never resolved it. Both forms are the same defect.
_SEGMENT = r"(?:[a-z]_[a-z0-9_]+|\.\.|[a-z][a-z0-9_]*)"
_LEAF = r"(?:index|\d\d_[a-z0-9_]+)"
_EXTENSION = r"(?:\.(?:md|ipynb))?"
DOC_SHAPED = re.compile(rf"^(?:{_SEGMENT}/)*{_LEAF}{_EXTENSION}$")

# A single-backtick span that is not part of a ``double`` span and is not the
# argument of a role such as {doc}`...` or {ref}`...`.
INLINE_LITERAL = re.compile(r"(?<!`)(?<!\})`([^`\n]+)`(?!`)")
ROLE_PREFIX = re.compile(r"\{[a-z:+_-]+\}$")

FENCE = re.compile(r"^\s*(```+|~~~+)")


def _markdown_pages() -> list[Path]:
    return sorted(
        path
        for path in DOCS_ROOT.rglob("*.md")
        if not EXCLUDED_PARTS & set(path.relative_to(DOCS_ROOT).parts)
    )


def _without_code_blocks(text: str) -> list[str]:
    """Blank out fenced code, keeping line numbers intact."""
    lines: list[str] = []
    fence: str | None = None
    for line in text.splitlines():
        match = FENCE.match(line)
        token = match.group(1)[:3] if match else None
        if fence is None:
            if token is not None:
                fence = token
                lines.append("")
                continue
            lines.append(line)
        else:
            lines.append("")
            if token == fence:
                fence = None
    return lines


def _bare_page_references() -> list[str]:
    offenders: list[str] = []
    for path in _markdown_pages():
        relative = path.relative_to(DOCS_ROOT).as_posix()
        if _is_exempt(str(relative)):
            continue
        for number, line in enumerate(
            _without_code_blocks(path.read_text(encoding="utf-8")),
            1,
        ):
            for match in INLINE_LITERAL.finditer(line):
                if ROLE_PREFIX.search(line[: match.start()]):
                    continue
                literal = match.group(1)
                if DOC_SHAPED.match(literal):
                    offenders.append(f"docs/{relative}:{number}: `{literal}`")
    return offenders


def test_internal_page_references_use_the_doc_role() -> None:
    """No page names another page in bare backticks."""
    offenders = _bare_page_references()
    assert not offenders, (
        "These inline literals name a page on this site. Sphinx never resolves "
        "or validates them, so a wrong target ships silently. Write each as a "
        "{doc} role with a target relative to the linking page:\n  "
        + "\n  ".join(offenders)
    )


def test_the_contract_detects_a_bare_reference() -> None:
    """The shape filter matches page names and not ordinary literals."""
    for page_name in (
        "06_troubleshooting_highlighting",
        "index",
        "../e_filtering/index",
        "a_orientation/04_ui_glossary_terminology",
        "k_sessions_sharing/index",  # dead: still a defect, still caught
    ):
        assert DOC_SHAPED.match(page_name), page_name

    for literal in (
        "vector_fields",  # API parameter that collides with a generated page
        "dataset_identity.json",
        "local-user",
        "prepare()",
        "cellucid/assets/css/README.md",
        "DataState",
    ):
        assert not DOC_SHAPED.match(literal), literal
