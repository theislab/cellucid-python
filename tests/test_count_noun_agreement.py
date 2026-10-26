"""A count and its noun agree: "1 cell", "2 cells", never "1 cells".

Every f-string that writes a count before a noun takes the noun from
``cellucid._plural.plural``; ``f"{n} cells"`` reads "1 cells" the day ``n`` is
one, and a hand-rolled ``"" if n == 1 else "s"`` drifts from the shared rule.
Both are refused here for the whole package. What is not a count is not
refused: a label before a plural noun (``f"{label} entries must be …"``) and an
upper-case limit constant, which is always above one.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

from cellucid._plural import plural

SOURCE = Path(__file__).resolve().parents[1] / "src" / "cellucid"

PLURAL_NOUNS = (
    r"cells|points|genes|values|entries|categories|options|columns|rows|pages|items|fields"
    r"|characters|edges|datasets|views|codes|groups|pixels|indices|stores|positions|frames"
    r"|samples|records|files|workers|results|arrays|layers|features|suggestions|requests|steps"
    r"|objects|modalities|embeddings|annotators|voters|votes|chunks|parts|matrices|sets|nodes"
    r"|users|labels|lines|errors|warnings|retries|attempts|bytes|digits|dimensions|axes|keys"
    r"|components|messages|times|reads"
)
NOUN_AFTER = re.compile(r"^\s+(" + PLURAL_NOUNS + r")\b")
LABEL = re.compile(
    r"^(label|name|noun|what|kind|key|path|prefix|owner|context|field|codes_type|position"
    r"|type\(.*\)\.__name__|.*\.dtype)$"
)
CONSTANT = re.compile(r"^[A-Z_][A-Z0-9_]*$")


def test_plural_takes_the_singular_for_one() -> None:
    assert plural(1, "cell") == "cell"
    assert plural(0, "cell") == "cells"
    assert plural(2, "cell") == "cells"
    assert plural(1, "entry", "entries") == "entry"
    assert plural(3, "entry", "entries") == "entries"


def _offenders(tree: ast.AST, where: str) -> list[str]:
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.JoinedStr):
            parts = node.values
            for value, following in zip(parts, parts[1:], strict=False):
                if not (isinstance(value, ast.FormattedValue) and isinstance(following, ast.Constant)):
                    continue
                if not isinstance(following.value, str):
                    continue
                match = NOUN_AFTER.match(following.value)
                if match is None:
                    continue
                expression = ast.unparse(value.value)
                if LABEL.match(expression) or CONSTANT.match(expression):
                    continue
                found.append(f"{where}:{node.lineno}: {{{expression}}} {match.group(1)}")
        if (
            isinstance(node, ast.IfExp)
            and isinstance(node.body, ast.Constant)
            and isinstance(node.orelse, ast.Constant)
            and node.body.value == ""
            and node.orelse.value == "s"
        ):
            found.append(f"{where}:{node.lineno}: {ast.unparse(node)}")
    return found


def test_no_count_is_written_before_a_plural_noun_it_may_not_agree_with() -> None:
    offenders: list[str] = []
    for path in sorted(SOURCE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        offenders += _offenders(tree, str(path.relative_to(SOURCE)))
    assert offenders == [], "use cellucid._plural.plural for the noun"


def test_the_check_catches_a_hard_coded_plural() -> None:
    tree = ast.parse('n = 1\nmessage = f"{n} cells, {n} genes{\'\' if n == 1 else \'s\'}"\n')
    assert len(_offenders(tree, "probe.py")) == 3
