"""One spelling: US, in every word a person or an agent reads from this package.

The web app held its strings to US spelling (``color``, ``neighbor``,
``canceled``, ``catalog``) while this package said "Compute the neighbour
graph", its docstrings "catalogue" and "recognise", and its documentation
"neighbour", "greyed", "centred" in prose and in screenshot alt texts
(round 9, R9N-1): the contract lived only in the web repository. The British
word list is shared byte for byte with the web app
(``tests/failure-wording-contract.test.mjs``) and cellucid-r
(``tests/testthat/test-spelling.R``), at the digest pinned below.

Read here: every string literal of ``src/cellucid`` that is prose (words with
a space between, or a one-word sentence such as "Canceled."), its docstrings
and ``#:`` attribute docs, the action catalog's titles and summaries, the
documentation (Markdown, reStructuredText, notebooks, the screenshot
descriptions), the repository's own pages and the unreleased section of the
change log. A word in backticks is a name (``normalise`` of a step recorded
before the rename, ``notifications/cancelled`` of MCP) and is left out, as is
a word that is part of an identifier, a path or a flag.
"""

from __future__ import annotations

import hashlib
import io
import json
import re
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "spelling.json"

#: The same digest is pinned in the web app (tests/failure-wording-contract.test.mjs)
#: and in cellucid-r (tests/testthat/test-spelling.R). Update all three together.
FIXTURE_SHA256 = "f195abcd84ce15a4d39575042ac1ad205a0b4e1b6b8b439d361d4c0dfbbfa77c"

_PATTERNS = [entry["pattern"] for entry in json.loads(FIXTURE.read_text(encoding="utf-8"))["british"]]
#: A British word in prose: never beside a letter, ``_``, ``.``, ``/``, ``\``,
#: ``$`` or ``--`` (an identifier, a path, a flag), nor assigned or called.
BRITISH = re.compile(
    r"(?<![\w.$/\\])(?<!--)(?:" + "|".join(_PATTERNS) + r")(?![\w$/\\]|\.\w|\s*=(?!=)|\()",
    re.IGNORECASE,
)
BACKTICKS = re.compile(r"`+[^`\n]*?`+")
#: A single name in quotes is a name the text quotes (a colormap id `"greys"`, a
#: parameter a message echoes: `has no parameter "colour"`), never its prose.
QUOTED_NAME = re.compile(r"\\?([\"'])[\w-]+\\?\1")
ONE_WORD_SENTENCE = re.compile(r"[A-Za-z]+[.!?…]+")


def british(text: str) -> list[str]:
    """The British spellings in prose ``text``, backtick spans left out."""
    return [match.group(0) for match in BRITISH.finditer(QUOTED_NAME.sub(" ", BACKTICKS.sub(" ", text)))]


def _is_prose(body: str) -> bool:
    return bool(re.search(r"\s", body)) or bool(ONE_WORD_SENTENCE.fullmatch(body))


def _string_body(token: str) -> str:
    match = re.match(r"[rbuRBUfF]*('''|\"\"\"|'|\")", token)
    assert match is not None, token
    quote = match.group(1)
    return token[match.end() : len(token) - len(quote)]


def _python_offenders(path: Path) -> list[str]:
    offenders = []
    source = path.read_text(encoding="utf-8")
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.STRING:
            body = _string_body(token.string)
        elif token.type == getattr(tokenize, "FSTRING_MIDDLE", -1) or (
            token.type == tokenize.COMMENT and token.string.startswith("#:")
        ):
            body = token.string
        else:
            continue
        if token.type != tokenize.COMMENT and not _is_prose(body):
            continue
        for word in british(body):
            offenders.append(f"{path.relative_to(ROOT)}:{token.start[0]}: {word}")
    return offenders


def _text_offenders(path: Path, text: str, first_line: int = 1) -> list[str]:
    return [
        f"{path.relative_to(ROOT)}:{first_line + index}: {word}"
        for index, line in enumerate(text.split("\n"))
        for word in british(line)
    ]


def _notebook_offenders(path: Path) -> list[str]:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    offenders = []
    for number, cell in enumerate(notebook.get("cells", [])):
        source = cell.get("source", "")
        text = "".join(source) if isinstance(source, list) else source
        offenders.extend(f"{entry} (cell {number})" for entry in _text_offenders(path, text))
    return offenders


def _unreleased_changelog() -> str:
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    released = re.search(r"^## \[\d+\.\d+\.\d+\] - \d{4}-\d{2}-\d{2}$", text, re.MULTILINE)
    assert released is not None, "CHANGELOG.md has a released section"
    return text[: released.start()]


def test_the_word_list_is_the_one_shared_with_the_web_app_and_r_at_the_pinned_digest() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    assert json.loads(FIXTURE.read_text(encoding="utf-8"))["version"] == 1


def test_the_check_finds_british_words_in_prose_and_never_in_names() -> None:
    for word, us in [
        ("neighbour", "neighbor"),
        ("Colours", "Colors"),
        ("greyed", "grayed"),
        ("cancelled", "canceled"),
        ("catalogue", "catalog"),
        ("unrecognised", "unrecognized"),
        ("reorganisation", "reorganization"),
        ("learnt", "learned"),
        ("analysed", "analyzed"),
        ("relabelled", "relabeled"),
        ("centred", "centered"),
    ]:
        assert british(f"the {word} one") == [word]
        assert british(f"the {us} one") == []
    assert british("asyncio.CancelledError notifications/cancelled --colour colour= normalise( `cancelled` x.cancelled") == []
    assert british("cancellation, analyses, organism, summaries, precise, raised") == []


def test_the_package_says_every_message_docstring_and_attribute_doc_in_us_spelling() -> None:
    offenders = []
    for path in sorted((ROOT / "src" / "cellucid").rglob("*.py")):
        offenders.extend(_python_offenders(path))
    for path in sorted((ROOT / "src" / "cellucid").rglob("*.json")):
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(r'"(?:title|summary|description)":\s*"((?:[^"\\]|\\.)*)"', text):
            for word in british(match.group(1)):
                offenders.append(f"{path.relative_to(ROOT)}:{text[: match.start()].count(chr(10)) + 1}: {word}")
    for path in sorted((ROOT / "src" / "cellucid").rglob("*.js")):
        for number, line in enumerate(path.read_text(encoding="utf-8").split("\n"), start=1):
            for literal in re.findall(r"(['\"`])((?:(?!\1)[^\\]|\\.)*?)\1", line):
                if _is_prose(literal[1]):
                    offenders.extend(f"{path.relative_to(ROOT)}:{number}: {word}" for word in british(literal[1]))
    assert offenders == []


def test_the_documentation_and_the_repository_pages_are_us_spelling() -> None:
    offenders = []
    docs = ROOT / "docs"
    for path in sorted(docs.rglob("*")):
        if "_build" in path.parts or not path.is_file():
            continue
        if path.suffix in {".md", ".rst", ".txt"}:
            offenders.extend(_text_offenders(path, path.read_text(encoding="utf-8")))
        elif path.suffix == ".ipynb":
            offenders.extend(_notebook_offenders(path))
        elif path.suffix in {".py", ".mjs", ".js", ".json"} and "_tooling" in path.parts:
            # The screenshot tooling: its descriptions are the figures' alt texts.
            text = path.read_text(encoding="utf-8")
            for number, line in enumerate(text.split("\n"), start=1):
                for literal in re.findall(r"(['\"`])((?:(?!\1)[^\\]|\\.)*?)\1", line):
                    if re.search(r"\s", literal[1]):
                        offenders.extend(f"{path.relative_to(ROOT)}:{number}: {word}" for word in british(literal[1]))
    for name in ("README.md", "CONTRIBUTING.md", "SECURITY.md", "SUPPORT.md", "CODE_OF_CONDUCT.md"):
        path = ROOT / name
        if path.is_file():
            offenders.extend(_text_offenders(path, path.read_text(encoding="utf-8")))
    offenders.extend(_text_offenders(ROOT / "CHANGELOG.md", _unreleased_changelog()))
    assert offenders == []


def test_the_tests_pin_the_words_the_package_says() -> None:
    # A test that expects British text forces the source back to it, or fails the day the source
    # is right (round 9: `test_action_catalog.py` still expected "cancelled" after the request
    # said "canceled"). Every prose string of the suite is read: an expectation, and an input the
    # package echoes back. A misspelled name a test feeds in on purpose (an unknown key
    # "colour") is one word, so no prose. This file's own examples are the word list's.
    offenders = []
    for path in sorted((ROOT / "tests").rglob("*.py")):
        if path.name != Path(__file__).name and "fixtures" not in path.parts:
            offenders.extend(_python_offenders(path))
    assert offenders == []
