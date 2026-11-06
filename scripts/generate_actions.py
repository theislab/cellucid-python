"""Write ``src/cellucid/jupyter/_actions_generated.py`` from the action catalogue.

The catalogue (``src/cellucid/action_catalog.json``) is the web app's
``catalogDocument`` of every action, byte for byte (written by the web
repository's ``scripts/generate-action-catalog.mjs``; its digest is pinned in
``tests/test_action_catalog.py``). Every action becomes one method of the
viewer (:class:`cellucid.jupyter.AnnDataViewer`, and so of the notebook
widget): its tool name (``color.by`` -> ``color_by``), its parameters as
typed keyword arguments, and a docstring made of its title, summary and the
description of each parameter. Each method calls ``self.act(type, ...)``,
which checks and converts the parameters by the catalogue and sends one
``action`` command: nothing is written by hand per action.

    python scripts/generate_actions.py           # write the module
    python scripts/generate_actions.py --check   # exit 1 when it is stale

``tests/test_action_catalog.py`` fails while the committed module differs
from what this writes.
"""

from __future__ import annotations

import json
import keyword
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = ROOT / "src" / "cellucid" / "action_catalog.json"
MODULE_PATH = ROOT / "src" / "cellucid" / "jupyter" / "_actions_generated.py"

#: What an action of each effect does, in its docstring.
_EFFECT_TEXT = {
    "state": "It changes the data's meaning (fields, labels, filters, highlights, results)",
    "view": "It changes how the cells are shown",
    "layout": "It changes the viewer's layout",
}


def python_name(name: str) -> str:
    """A parameter's keyword in Python: a keyword gets an underscore (``in_``)."""
    return f"{name}_" if keyword.iskeyword(name) or keyword.issoftkeyword(name) else name


def _types_of(schema: dict[str, Any]) -> list[str]:
    declared = schema.get("type")
    if declared is None:
        return []
    return list(declared) if isinstance(declared, list) else [declared]


def _annotated(schema: dict[str, Any], key: str) -> bool:
    return schema.get(key) is True or any(_annotated(sub, key) for sub in schema.get("oneOf", []))


def _literal(values: list[Any]) -> str:
    plain = [value for value in values if value is not None]
    parts = [f"Literal[{', '.join(repr(value) for value in plain)}]"] if plain else []
    if len(plain) != len(values):
        parts.append("None")
    return " | ".join(parts)


def annotation(schema: dict[str, Any]) -> str:
    """The Python type a parameter of this schema takes."""
    if _annotated(schema, "x-cells"):
        nullable = "null" in _types_of(schema) or any("null" in _types_of(sub) for sub in schema.get("oneOf", []))
        return "CellsLike | None" if nullable else "CellsLike"
    if _annotated(schema, "x-values"):
        return "ValuesLike"
    if "const" in schema:
        return _literal([schema["const"]])
    if "enum" in schema:
        return _literal(list(schema["enum"]))
    if "oneOf" in schema:
        parts: list[str] = []
        for sub in schema["oneOf"]:
            for part in annotation(sub).split(" | "):
                if part not in parts:
                    parts.append(part)
        return " | ".join(parts)
    parts = []
    for kind in _types_of(schema):
        if kind == "string":
            parts.append("str")
        elif kind == "integer":
            parts.append("int")
        elif kind == "number":
            parts.append("float")
        elif kind == "boolean":
            parts.append("bool")
        elif kind == "null":
            parts.append("None")
        elif kind == "object":
            parts.append("Mapping[str, Any]")
        elif kind == "array":
            parts.append(f"Sequence[{annotation(schema['items'])}]")
    if not parts:
        return "Any"
    # None last, as Python writes optional types.
    return " | ".join([part for part in parts if part != "None"] + [part for part in parts if part == "None"])


_RST_SPECIAL = re.compile(r"([\\*`|])")
_RST_UNDERSCORE = re.compile(r"(\w)_(?=\W|$)")


def rst(text: str) -> str:
    """Catalogue text as reStructuredText that reads as written."""
    return _RST_UNDERSCORE.sub(r"\1\\_", _RST_SPECIAL.sub(r"\\\1", text))


def _type_words(schema: dict[str, Any]) -> str:
    """A parameter's type, as the docstring says it."""
    if _annotated(schema, "x-cells"):
        return "cells (rows, a boolean mask or obs names)"
    if _annotated(schema, "x-values"):
        return "array of numbers, one per cell"
    if "const" in schema:
        return json.dumps(schema["const"])
    if "enum" in schema:
        return "{" + ", ".join(json.dumps(value) for value in schema["enum"]) + "}"
    if "oneOf" in schema:
        return " or ".join(_type_words(sub) for sub in schema["oneOf"])
    words = {
        "string": "str", "integer": "int", "number": "float", "boolean": "bool",
        "null": "None", "object": "dict", "array": "list",
    }
    return " or ".join(words[kind] for kind in _types_of(schema)) or "any"


def _keys_sentence(schema: dict[str, Any]) -> str:
    """The keys of an object parameter, with which are required."""
    objects = [schema] if "object" in _types_of(schema) else [
        sub for sub in schema.get("oneOf", []) if "object" in _types_of(sub)
    ]
    if len(objects) != 1 or not objects[0].get("properties"):
        return ""
    properties = objects[0]["properties"]
    required = set(objects[0].get("required", []))
    keys = ", ".join(f"``{key}``{'' if key in required else ' (optional)'}" for key in properties)
    return f" Keys: {keys}."


def _wrap(text: str, indent: str, width: int = 96) -> list[str]:
    words = text.split()
    lines: list[str] = []
    line = ""
    for word in words:
        if line and len(indent) + len(line) + 1 + len(word) > width:
            lines.append(indent + line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    if line:
        lines.append(indent + line)
    return lines


_SCALARS = (bool, int, float, str, type(None))


def _default_of(schema: dict[str, Any]) -> str | None:
    """A parameter's Python default: its scalar schema default, else OMIT."""
    if "default" in schema and isinstance(schema["default"], _SCALARS):
        return repr(schema["default"])
    return None


def method_source(action: dict[str, Any]) -> str:
    """One generated method."""
    schema = action["params"]
    properties: dict[str, Any] = schema.get("properties", {})
    required = list(schema.get("required", []))
    optional = [name for name in properties if name not in required]
    signature = ["self", "*"] if properties else ["self"]
    for name in required:
        signature.append(f"{python_name(name)}: {annotation(properties[name])}")
    for name in optional:
        default = _default_of(properties[name])
        if default is None:
            signature.append(f"{python_name(name)}: {annotation(properties[name])} | Omit = OMIT")
        else:
            signature.append(f"{python_name(name)}: {annotation(properties[name])} = {default}")
    tool = action["tool"]
    if not tool.isidentifier() or keyword.iskeyword(tool):
        raise ValueError(f"{action['type']}: {tool!r} is not a method name")

    doc = [f"{rst(action['title'])}.", "", *_wrap(rst(action["summary"]), "")]
    doc.append("")
    if action["effect"] == "read":
        effect = (
            f"Action ``{action['type']}``. It reads only and changes nothing; "
            "``request.result()`` is what it reads."
        )
    elif action["effect"] == "history":
        effect = (
            f"Action ``{action['type']}``. It moves the viewer's history, as the person's Undo, Redo "
            "or a click on a step does, and records no step; ``request.result()`` is its outcome."
        )
    elif action["effect"] == "dataset":
        effect = (
            f"Action ``{action['type']}``. It opens another dataset in place of the one the viewer "
            "shows, as an agent in the viewer does. A notebook's viewer shows the dataset its notebook "
            "serves, so the viewer refuses it from Python: show another dataset with "
            ":func:`cellucid.show_anndata` or :func:`cellucid.serve_anndata`."
        )
    elif action["effect"] == "outside":
        effect = (
            f"Action ``{action['type']}``. It acts outside the viewer's state and history (a file the "
            "person receives in the browser, a publication, or what the browser keeps for later "
            "visits), which Undo does not take back, and records no step; ``request.result()`` is "
            "its outcome."
        )
        if action.get("confirm") == "publishes":
            effect += (
                " It publishes for others to see, so the viewer asks the person first, from Python too, "
                "as when an agent runs it; ``request.result()`` waits for their answer (call it again if "
                "its timeout passes first)."
            )
    else:
        effect = (
            f"Action ``{action['type']}``. {_EFFECT_TEXT[action['effect']]}, as one step of the viewer's "
            "history, which the person can undo; ``request.result()`` is its outcome and "
            "``request.step`` the step."
        )
    doc.extend(_wrap(effect, ""))
    if properties:
        doc.extend(["", "Parameters", "----------"])
        for name in [*required, *optional]:
            sub = properties[name]
            doc.append(f"{python_name(name)} : {rst(_type_words(sub))}")
            text = rst(sub.get("description", "")).strip()
            if text and not text.endswith("."):
                text += "."
            text += _keys_sentence(sub)
            if name in optional:
                default = _default_of(sub)
                text += (
                    " Optional: left out, the viewer decides." if default is None
                    else f" Default ``{default}``."
                )
            doc.extend(_wrap(text.strip() or "See the action's summary.", "    "))
    doc.extend([
        "",
        "Returns",
        "-------",
        "ActionRequest",
        *_wrap(
            "The request, sent at once (or as soon as the viewer connects). In a notebook widget "
            "its answer arrives once the cell has finished: collect it in a later cell.",
            "    ",
        ),
    ])
    body_doc = "\n".join(f"        {line}" if line else "" for line in doc)
    if '"""' in body_doc:
        raise ValueError(f"{action['type']}: its text holds a triple quote")
    args = ", ".join(f"{name}={python_name(name)}" if name.isidentifier() and not keyword.iskeyword(name)
                     else f"**{{{name!r}: {python_name(name)}}}" for name in [*required, *optional])
    call = f'self.act({action["type"]!r}{", " + args if args else ""})'
    lines = [f"    def {tool}("]
    lines.extend(f"        {part}," for part in signature)
    lines.append("    ) -> ActionRequest:")
    lines.append('        r"""' + body_doc.lstrip())
    lines.append('        """')
    lines.append(f"        return {call}")
    return "\n".join(lines)


HEADER = '''"""The viewer's methods, one per action of the catalogue: generated, do not edit.

Written by ``scripts/generate_actions.py`` from ``cellucid/action_catalog.json``
(the web app's action catalogue, digest {digest}). Each method sends one
``action`` command through :meth:`AnnDataViewer.act`, which checks and converts
its parameters by the catalogue, and returns its :class:`ActionRequest`.
"""

from __future__ import annotations

{imports}

if TYPE_CHECKING:
    from ._viewer import ActionRequest


class ViewerActions:
    """Every action of the catalogue as a method (mixed into the viewer)."""

    __slots__ = ()

    if TYPE_CHECKING:

        def act(self, action: str, /, **params: Any) -> ActionRequest: ...
'''


def render(catalog_text: str) -> str:
    """The module's text for a catalogue document."""
    import hashlib

    document = json.loads(catalog_text)
    if document.get("schema") != 1:
        raise ValueError("The action catalogue must be schema 1")
    digest = hashlib.sha256(catalog_text.encode("utf-8")).hexdigest()
    methods = [method_source(action) for action in sorted(document["actions"], key=lambda entry: entry["tool"])]
    body = "\n\n".join(methods)

    def used(*names: str) -> list[str]:
        return [name for name in names if re.search(rf"\b{name}\b", body) is not None]

    groups = [
        ("collections.abc", used("Mapping", "Sequence")),
        ("typing", ["TYPE_CHECKING", "Any", *used("Literal")]),
    ]
    imports = [f"from {module} import {', '.join(names)}" for module, names in groups if names]
    local = used("OMIT", "CellsLike", "Omit", "ValuesLike")
    text = (
        HEADER.replace("{digest}", digest[:16]).replace(
            "{imports}", "\n".join(imports) + ("\n\nfrom ._actions import " + ", ".join(local) if local else "")
        )
        + "\n"
        + body
        + "\n"
    )
    return text


def main(argv: list[str]) -> int:
    check = "--check" in argv
    unknown = [argument for argument in argv if argument != "--check"]
    if unknown:
        print(f"Unknown argument {unknown[0]}; the only one is --check", file=sys.stderr)
        return 2
    text = render(CATALOG_PATH.read_text(encoding="utf-8"))
    current = MODULE_PATH.read_text(encoding="utf-8") if MODULE_PATH.exists() else None
    if current == text:
        print(f"current {MODULE_PATH.relative_to(ROOT)}")
        return 0
    if check:
        print(f"stale   {MODULE_PATH.relative_to(ROOT)}: run python scripts/generate_actions.py")
        return 1
    MODULE_PATH.write_text(text, encoding="utf-8")
    print(f"wrote   {MODULE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
