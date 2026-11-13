from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"


# A path that exists only on the machine the documentation was written on. The
# placeholders the guides teach with — ``/Users/you/…``, ``/home/you``,
# ``/Users/<username>/`` — are deliberately excluded: they are instructions, not
# captured output. Everything else naming a real home or scratch directory is a
# leak.
_CAPTURE_HOST_PATH = re.compile(
    r"/(?:Users|home)/(?!you\b|<username>|<user>)[A-Za-z0-9._-]+/"
    r"|/private/tmp/"
    r"|\.pyenv/versions/"
)


def test_no_documentation_page_names_a_path_from_the_machine_it_was_written_on() -> None:
    """Committed docs must not carry the author's home, scratch, or pyenv paths.

    Screenshots hide this class of leak in pixels, where no scan reaches, which
    is why terminal transcripts ship as ``text`` blocks. Prose has no such
    excuse, so the same rule is enforced here mechanically. ``docs/**`` is
    walked directly rather than through Sphinx: ``_tooling`` is in
    ``exclude_patterns`` and so is never built, but it is still committed to a
    public repository.
    """
    offenders: list[str] = []
    for page in sorted(DOCS_ROOT.rglob("*.md")):
        if "_build" in page.parts:
            continue
        for number, line in enumerate(page.read_text(encoding="utf-8").splitlines(), 1):
            if _CAPTURE_HOST_PATH.search(line):
                offenders.append(f"{page.relative_to(REPOSITORY_ROOT)}:{number}")
    assert offenders == [], offenders
