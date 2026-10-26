"""Every Python file parses on the oldest Python this package supports.

The suite usually runs on a newer Python than ``requires-python`` admits, and
syntax that only a newer Python accepts then passes here and fails on the
oldest one, in CI or for a user. ``ast.parse(..., feature_version=)`` refuses
most of it (``type`` statements, type parameter lists), but not the f-strings
Python 3.12 made legal (PEP 701), which are the likeliest to slip in; those
are found from the tokens of each f-string.
"""

from __future__ import annotations

import ast
import io
import re
import sys
import token
import tokenize
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ("src", "tests", "scripts", "docs")


def _minimum_python() -> tuple[int, int]:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^requires-python\s*=\s*">=\s*(\d+)\.(\d+)', text, re.MULTILINE)
    assert match is not None, "pyproject.toml declares no requires-python lower bound"
    return int(match.group(1)), int(match.group(2))


MINIMUM = _minimum_python()


def _python_files() -> list[Path]:
    files = [
        path
        for directory in SOURCES
        for path in sorted((ROOT / directory).rglob("*.py"))
        if "_build" not in path.relative_to(ROOT).parts
    ]
    assert any(path.name == "_prepare.py" for path in files) and any(
        path == Path(__file__) for path in files
    )
    return files


def _pep701_f_strings(source: str) -> list[str]:
    """The f-strings in ``source`` that only Python 3.12 and later accept.

    Before 3.12 an f-string's replacement fields could not reuse the quote
    that encloses the f-string, hold a backslash or a comment, or, in a
    single-quoted f-string, span lines. The tokenizer of Python 3.12 and later
    splits an f-string into tokens, so each replacement field is seen as the
    code it is.
    """
    found: list[str] = []
    # The quote of each open f-string, innermost last: ', ", ''' or """.
    quotes: list[str] = []
    for item in tokenize.generate_tokens(io.StringIO(source).readline):
        kind, text, line = item.type, item.string, item.start[0]
        if kind == token.FSTRING_END:
            quotes.pop()
            continue
        # A literal part of the outermost f-string is not inside any field.
        if not quotes or (kind == token.FSTRING_MIDDLE and len(quotes) == 1):
            if kind == token.FSTRING_START:
                quotes.append(_quote(text))
            continue
        if kind in (token.STRING, token.FSTRING_START) and _quote(text).startswith(quotes[-1]):
            found.append(
                f"line {line}: a string inside an f-string field reuses the enclosing quote"
            )
        if kind == token.COMMENT:
            found.append(f"line {line}: a comment inside an f-string field")
        if kind == token.NL and len(quotes[-1]) == 1:
            found.append(f"line {line}: a single-quoted f-string field spans lines")
        if "\\" in text:
            found.append(f"line {line}: a backslash inside an f-string field")
        if kind == token.FSTRING_START:
            quotes.append(_quote(text))
    return found


def _quote(text: str) -> str:
    """The opening quote of a string token: ', ", ''' or \"\"\"."""
    body = text.lstrip("rRfFbBuU")
    return body[:3] if body[:3] in ("'''", '"""') else body[:1]


@pytest.mark.parametrize("path", _python_files(), ids=lambda path: str(path.relative_to(ROOT)))
def test_the_file_parses_with_the_oldest_supported_grammar(path: Path) -> None:
    source = path.read_text(encoding="utf-8")
    ast.parse(source, filename=str(path), feature_version=MINIMUM)
    if sys.version_info >= (3, 12) and MINIMUM < (3, 12):
        assert _pep701_f_strings(source) == [], f"{path.relative_to(ROOT)} needs Python 3.12"


@pytest.mark.skipif(
    sys.version_info < (3, 12), reason="the tokenizer splits f-strings from Python 3.12"
)
@pytest.mark.parametrize(
    "source",
    [
        'x = f"{"a"}"\n',
        "x = f'{'a'}'\n",
        'x = f"{f"{1}"}"\n',
        'x = f"{"\\n".join(y)}"\n',
        'x = f"{y # note\n}"\n',
        'x = f"{y +\n 1}"\n',
    ],
)
def test_a_python_312_f_string_is_found(source: str) -> None:
    assert _pep701_f_strings(source)


@pytest.mark.parametrize(
    "source",
    [
        "x = f\"{'a'}\"\n",
        "x = f'{y!r:>{width}}'\n",
        'x = f"""{\n y\n}"""\n',
        "x = f\"{f'{1}'}\"\n",
        'x = f"\\n{y}"\n',
        'x = "{}".format(y)\n',
    ],
)
def test_an_f_string_older_pythons_accept_is_not_found(source: str) -> None:
    ast.parse(source, feature_version=(3, 11))
    if sys.version_info >= (3, 12):
        assert _pep701_f_strings(source) == []
