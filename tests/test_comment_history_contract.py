"""Comments, docstrings, test names and workflow comments hold no history: no finding, audit or plan id,
no round or hunt story, no CI run address and no pointer into the private plans.

The rule is "Comments say what the code does" in CONTRIBUTING.md; the patterns are the
hard classes of the history survey (an id is an exact shape, so a hit is a defect).
"""

from __future__ import annotations

import ast
import io
import re
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TREES = ("src", "tests", "scripts")
WORKFLOWS = ".github/workflows"
JS_TREE = "docs/_tooling"

_NOT_AFTER = r"(?<!kbd )(?<!noqa: )(?<!Shift\+)(?<!Ctrl\+)(?<!Alt\+)(?<!press )(?<!key )(?<!`)"

PATTERNS: dict[str, re.Pattern[str]] = {
    "CI run": re.compile(r"actions/runs/\d+|\brun \d{8,}\b|\b4[0-4]\d{9}\b"),
    "finding id": re.compile(
        _NOT_AFTER
        + r"(?:\bF\d{2,3}\b|\bL\d{1,2}-\d{1,2}\b|\bR\d{1,2}[A-Za-z]{0,2}-\d{1,3}\b"
        + r"|\bRC[A-Z]-\d{1,2}\b|\bCEL-(?:AUDIT-)?\d{2,4}\b|\bCEL-[A-Z]+-\d+\b)"
    ),
    "bug hunt": re.compile(
        r"bug[- ]?hunt|\bround[- ]\d+\b|\blive[- ]pass\b|\blive[- ]\d+\b", re.IGNORECASE
    ),
    "plan pointer": re.compile(r"History v2 ?§|cellucid-private"),
}


def find_history(text: str) -> list[tuple[str, str]]:
    """Return (class, matched text) for every history marker in one comment or docstring."""
    return [(name, m.group(0)) for name, rx in PATTERNS.items() for m in rx.finditer(text)]


def _py_files() -> list[Path]:
    files: list[Path] = []
    for tree in TREES:
        files += [
            p
            for p in (ROOT / tree).rglob("*.py")
            if "__pycache__" not in p.parts and ".egg-info" not in str(p)
        ]
    return sorted(files)


def _workflow_files() -> list[Path]:
    return sorted([*(ROOT / WORKFLOWS).glob("*.yml"), *(ROOT / WORKFLOWS).glob("*.yaml")])


def yaml_comments(source: str) -> list[tuple[int, str]]:
    """(line, text) of each `#` comment of a YAML file: a `#` that opens a line or follows
    white space outside a quoted string."""
    found: list[tuple[int, str]] = []
    for number, text in enumerate(source.splitlines(), start=1):
        quote = ""
        i = 0
        while i < len(text):
            ch = text[i]
            if quote:
                if ch == "\\" and quote == '"':
                    i += 1
                elif ch == quote:
                    quote = ""
            elif ch in "\"'" and (i == 0 or text[i - 1] in " \t:[{,-"):
                quote = ch
            elif ch == "#" and (i == 0 or text[i - 1].isspace()):
                found.append((number, text[i + 1 :]))
                break
            i += 1
    return found


def _js_files() -> list[Path]:
    return sorted(p for p in (ROOT / JS_TREE).rglob("*.mjs") if "node_modules" not in p.parts)


def js_comments(source: str) -> list[tuple[int, str]]:
    """(line, text) of each `//` and `/* */` comment of a JavaScript file: a scan that skips quoted
    strings, template literals (with `${}` nesting) and regular-expression literals."""
    found: list[tuple[int, str]] = []
    i, line, n = 0, 1, len(source)
    braces: list[int] = []  # template-literal nesting: brace depth at each open `${`
    depth = 0
    in_template = False
    last = ""  # last significant character, to tell a regex `/` from a division
    while i < n:
        ch = source[i]
        nxt = source[i + 1] if i + 1 < n else ""
        if in_template:
            if ch == "\\":
                i += 2
                continue
            if ch == "\n":
                line += 1
            if ch == "`":
                in_template = False
                last = "`"
            elif ch == "$" and nxt == "{":
                braces.append(depth)
                depth += 1
                in_template = False
                i += 1
            i += 1
            continue
        if ch == "\n":
            line += 1
        elif ch == "/" and nxt == "/":
            end = source.find("\n", i)
            end = n if end < 0 else end
            found.append((line, source[i + 2 : end]))
            i = end
            continue
        elif ch == "/" and nxt == "*":
            end = source.find("*/", i + 2)
            end = n if end < 0 else end
            body = source[i + 2 : end]
            for offset, text in enumerate(body.split("\n")):
                found.append((line + offset, text))
            line += body.count("\n")
            i = end + 2
            continue
        elif ch in "\"'":
            j = i + 1
            while j < n and source[j] != ch and source[j] != "\n":
                j += 2 if source[j] == "\\" else 1
            i = j + 1
            last = ch
            continue
        elif ch == "`":
            in_template = True
        elif ch == "/" and (last == "" or last in "(,=:[!&|?{};+-*%<>~^"):
            j = i + 1
            in_class = False
            while j < n and source[j] != "\n":
                c = source[j]
                if c == "\\":
                    j += 1
                elif c == "[":
                    in_class = True
                elif c == "]":
                    in_class = False
                elif c == "/" and not in_class:
                    break
                j += 1
            i = j + 1
            last = "/"
            continue
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if braces and depth == braces[-1]:
                braces.pop()
                in_template = True
        if not ch.isspace():
            last = ch
        i += 1
    return found


def scan_js(source: str) -> list[tuple[int, str, str]]:
    """(line, class, text) for each marker in the comments of a JavaScript file."""
    return [(n, c, t) for n, comment in js_comments(source) for c, t in find_history(comment)]


def scan_workflow(source: str) -> list[tuple[int, str, str]]:
    """(line, class, text) for each marker in the comments of a workflow."""
    return [(n, c, t) for n, comment in yaml_comments(source) for c, t in find_history(comment)]


def scan_source(source: str) -> list[tuple[int, str, str]]:
    """(line, class, text) for each marker in the comments, docstrings and test names of a module."""
    hits: list[tuple[int, str, str]] = []
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok.type == tokenize.COMMENT:
            hits += [(tok.start[0], c, t) for c, t in find_history(tok.string)]
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc is not None:
                line = node.body[0].lineno if node.body else 1
                hits += [(line, c, t) for c, t in find_history(doc)]
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith(
            "test"
        ):
            hits += [(node.lineno, c, t) for c, t in find_history(node.name.replace("_", " "))]
            spaced = node.name.replace("_", " ").upper()
            hits += [(node.lineno, c, t) for c, t in find_history(spaced) if c == "finding id"]
    return hits


def test_checker_finds_planted_history() -> None:
    planted = [
        "# F123 was a bug",
        "# see L4-12 for the cause",
        "# R9a-3: wrong order",
        "# CEL-AUDIT-42",
        "# found in round 9",
        "# from the bug hunt",
        "# the live pass showed it",
        "# https://github.com/theislab/cellucid/actions/runs/123456",
        "# run 12345678901 failed",
        "# cellucid-private/plans/open-items.md",
        "# History v2 §3.8",
    ]
    for line in planted:
        assert find_history(line), line
    assert scan_source("def test_f123_case():\n    pass\n")
    assert scan_source("def test_round_9_case():\n    pass\n")
    assert scan_source('def f():\n    """Cause: F123."""\n')
    assert scan_source('x = "F123"  # plain\n') == []
    assert scan_source('x = "F123"\n') == []
    assert find_history("# the key F5 and a 4-byte id") == []


def test_checker_finds_history_in_workflow_comments_only() -> None:
    assert scan_workflow("name: a\n# run 12345678901 failed\nsteps:\n  - run: x # round 9\n") == [
        (2, "CI run", "run 12345678901"),
        (4, "bug hunt", "round 9"),
    ]
    assert scan_workflow("env:\n  A: \"a # round 9\"\n  U: https://x.test/p#round-9\n  B: '#round 9'\n") == []


def test_checker_finds_history_in_javascript_comments_only() -> None:
    planted = (
        "const a = 1; // bug hunt F74\n"
        "/**\n * Crops the panel (round 9).\n */\n"
        "const b = `x ${y /* run 12345678901 */}`;\n"
    )
    assert scan_js(planted) == [
        (1, "finding id", "F74"),
        (1, "bug hunt", "bug hunt"),
        (3, "bug hunt", "round 9"),
        (5, "CI run", "run 12345678901"),
    ]
    code = (
        "const s = 'F74 // round 9';\nconst t = `bug hunt ${1}`;\n"
        "const r = /https?:\\/\\/[/]F74/;\nconst d = a / b; // fine\n"
    )
    assert scan_js(code) == []


def test_javascript_tooling_comments_hold_no_history() -> None:
    problems = []
    for path in _js_files():
        for line, cls, text in scan_js(path.read_text(encoding="utf-8")):
            problems.append(f"{path.relative_to(ROOT)}:{line}: {cls}: {text!r}")
    assert _js_files()
    assert not problems, "history in JavaScript comments:\n" + "\n".join(problems)


def test_workflow_comments_hold_no_history() -> None:
    problems = []
    for path in _workflow_files():
        for line, cls, text in scan_workflow(path.read_text(encoding="utf-8")):
            problems.append(f"{path.relative_to(ROOT)}:{line}: {cls}: {text!r}")
    assert _workflow_files()
    assert not problems, "history in workflow comments:\n" + "\n".join(problems)


def test_comments_docstrings_and_test_names_hold_no_history() -> None:
    problems = []
    for path in _py_files():
        for line, cls, text in scan_source(path.read_text(encoding="utf-8")):
            problems.append(f"{path.relative_to(ROOT)}:{line}: {cls}: {text!r}")
    assert not problems, "history in comments:\n" + "\n".join(problems)
