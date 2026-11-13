"""Write ``src/cellucid/_history_model.py`` from the history model.

The model (``src/cellucid/history_model.schema.json``) is the web app's
``historyModelDocument`` byte for byte (its home is the web repository's
``tests/fixtures/history-model.schema.json``, written by
``scripts/generate-history-model.mjs``; ``tests/test_history_model.py``
compares the two files): every record a
history, a saved state and a session file hold -- trackers, branches, steps,
notes, turns, the journal, the apply plan and its report -- with their fields,
types, enumerations and rules. This writes one frozen dataclass per record
(fields in the model's order, in snake case: ``headId`` is ``head_id``), the
model's constants and enumerations, and the field table the session reader
checks every record against; nothing about a record is written by hand.

Two records are maps rather than dataclasses: ``DomainValues`` (a value hash
per state domain) and ``Reads`` (the same, and ``history``, a branch's tip),
whose keys are domain ids. ``Tracker.steps`` is ``step_records`` here: the
reading model's ``Tracker.steps`` is the steps as a table
(:mod:`cellucid.history`).

    python scripts/generate_history_model.py           # write the module
    python scripts/generate_history_model.py --check   # exit 1 when it is stale
"""

from __future__ import annotations

import json
import keyword
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = ROOT / "src" / "cellucid" / "history_model.schema.json"
MODULE_PATH = ROOT / "src" / "cellucid" / "_history_model.py"

#: Fields whose attribute is not their snake-case name, and why.
_ATTRIBUTE_OVERRIDES = {
    # The reading model's Tracker.steps is the table of the steps.
    ("Tracker", "steps"): "step_records",
}

#: Records read as maps of their keys (domain ids), not as dataclasses.
_MAP_RECORDS = ("DomainValues", "Reads")

#: The model's enumerations, by the field that names them.
_ENUMS = {
    "BY_ORIGINS": ("By", "origin"),
    "HISTORY_HOSTS": ("By", "host"),
    "STEP_STATUSES": ("Step", "status"),
    "STEP_KINDS": ("Step", "kind"),
    "SQUASH_CAUSES": ("Squash", "cause"),
    "TRACKER_START_KINDS": ("Start", "kind"),
    "REPLAY_OUTCOMES": ("ReplayEntry", "outcome"),
    "SHARE_CONTENTS": ("Scope", "content"),
    "SHARE_RESULTS": ("Scope", "results"),
    "SHARE_TIMES": ("Scope", "times"),
    "SHARE_SOURCES": ("Scope", "source"),
    "APPLY_ROW_KINDS": ("ApplyPlanRow", "kind"),
    "APPLY_ROW_STATUSES": ("ApplyPlanRow", "status"),
    "APPLY_PACKAGES": ("Provenance", "package"),
}


def snake(name: str) -> str:
    """``headId`` -> ``head_id``; snake-case keys stay as they are."""
    return re.sub(r"(?<=[a-z0-9])([A-Z])", lambda match: "_" + match.group(1).lower(), name)


def attribute(record: str, key: str) -> str:
    """A field's attribute: its snake-case name, a Python keyword with an
    underscore (``from`` is ``from_``), or its override."""
    if (record, key) in _ATTRIBUTE_OVERRIDES:
        return _ATTRIBUTE_OVERRIDES[(record, key)]
    name = snake(key)
    return f"{name}_" if keyword.iskeyword(name) else name


def spec_of(schema: dict[str, Any], hash_pattern: str) -> dict[str, Any]:
    """A field's spec, as the model writes it (``model.js`` ``specSchema`` read back)."""
    nullable = False
    if "anyOf" in schema:
        branches = schema["anyOf"]
        if len(branches) != 2 or branches[1] != {"type": "null"} or "$ref" not in branches[0]:
            raise ValueError(f"Unread anyOf: {schema}")
        return {"type": "record", "name": branches[0]["$ref"].rsplit("/", 1)[1], "nullable": True}
    if "$ref" in schema:
        return {"type": "record", "name": schema["$ref"].rsplit("/", 1)[1], "nullable": False}
    if "enum" in schema:
        values = list(schema["enum"])
        nullable = None in values
        return {"type": "enum", "values": [value for value in values if value is not None], "nullable": nullable}
    if schema.get("x-time"):
        return {"type": "time", "nullable": True}
    declared = schema.get("type")
    types = list(declared) if isinstance(declared, list) else ([] if declared is None else [declared])
    nullable = "null" in types
    types = [kind for kind in types if kind != "null"]
    if not types:
        if set(schema) - {"description"}:
            raise ValueError(f"Unread spec: {schema}")
        return {"type": "json", "nullable": True}
    kind = types[0]
    if kind == "string":
        if "x-id" in schema:
            return {"type": "id", "kind": schema["x-id"], "nullable": nullable}
        if schema.get("x-name"):
            return {"type": "name", "nullable": nullable}
        if schema.get("pattern") == hash_pattern:
            return {"type": "hash", "nullable": nullable}
        if "pattern" in schema:
            return {"type": "string", "pattern": schema["pattern"], "nullable": nullable}
        if schema.get("minLength") == 1:
            return {"type": "text", "nullable": nullable}
        return {"type": "string", "pattern": None, "nullable": nullable}
    if kind == "boolean":
        return {"type": "boolean", "nullable": nullable}
    if kind == "integer":
        return {"type": "integer", "minimum": schema["minimum"], "maximum": schema.get("maximum"), "nullable": nullable}
    if kind == "number":
        return {"type": "number", "minimum": schema["minimum"], "nullable": nullable}
    if kind == "object":
        return {"type": "object", "nullable": nullable}
    if kind == "array":
        if "prefixItems" in schema:
            return {"type": "pair", "items": [spec_of(item, hash_pattern) for item in schema["prefixItems"]], "nullable": nullable}
        return {"type": "array", "item": spec_of(schema["items"], hash_pattern), "nullable": nullable}
    raise ValueError(f"Unread spec: {schema}")


def annotation(spec: dict[str, Any]) -> str:
    """A field's Python type."""
    kind = spec["type"]
    base = {
        "text": "str", "string": "str", "name": "str", "hash": "str", "id": "str", "enum": "str",
        "boolean": "bool", "integer": "int", "number": "float", "object": "dict[str, Any]", "json": "Any",
    }.get(kind)
    if kind == "time":
        return "float | None"
    if kind == "record":
        base = "dict[str, Any]" if spec["name"] in _MAP_RECORDS else spec["name"]
    elif kind == "array":
        base = f"tuple[{annotation(spec['item'])}, ...]"
    elif kind == "pair":
        base = f"tuple[{', '.join(annotation(item) for item in spec['items'])}]"
    if base is None:
        raise ValueError(f"No annotation for {spec}")
    if base == "Any":
        return base
    return f"{base} | None" if spec.get("nullable") else base


def docstring(text: str, indent: str) -> str:
    escaped = text.replace("\\", "\\\\").replace('"""', '\\"\\"\\"')
    return f'{indent}"""{escaped}"""'


def render(schema_text: str) -> str:
    schema = json.loads(schema_text)
    hash_pattern = schema["$defs"]["AppStamp"]["properties"]["catalog"]["pattern"]
    records: dict[str, list[tuple[str, str, dict[str, Any], bool, str]]] = {}
    for name, definition in schema["$defs"].items():
        required = set(definition["required"])
        fields = []
        # The document's keys are sorted; its required list keeps the model's
        # order, so required fields come in that order, then the optional ones.
        order = [*definition["required"], *sorted(set(definition["properties"]) - required)]
        for key in order:
            field = definition["properties"][key]
            spec = spec_of(field, hash_pattern)
            fields.append((key, attribute(name, key), spec, key not in required, field.get("description", "")))
        records[name] = fields
    constants = schema["x-constants"]
    out: list[str] = []
    out.append('"""')
    out.append("The history model, as Python reads it: one frozen dataclass per record.")
    out.append("")
    out.append("Generated by ``scripts/generate_history_model.py`` from")
    out.append("``history_model.schema.json`` (the web app's ``app/history/model.js``")
    out.append("``historyModelDocument``); do not edit by hand. The reading model")
    out.append("(:mod:`cellucid.history`) and the session reader (:mod:`cellucid.session`)")
    out.append("build every record from these classes and check every record against")
    out.append("``FIELDS``.")
    out.append('"""')
    out.append("")
    out.append("from __future__ import annotations")
    out.append("")
    out.append("import dataclasses")
    out.append("from collections.abc import Mapping as _Mapping")
    out.append("from dataclasses import dataclass")
    out.append("from typing import Any, Final")
    out.append("")
    out.append("#: The model's version (``HISTORY_MODEL_VERSION``).")
    out.append(f"HISTORY_MODEL_VERSION: Final = {schema['version']}")
    out.append("")
    out.append("#: Id forms, as regular expressions (``HISTORY_ID_PATTERNS``).")
    out.append(f"ID_PATTERNS: Final[_Mapping[str, str]] = {schema['x-ids']!r}")
    out.append("#: A value's content hash (SHA-256, lower-case hex).")
    out.append(f"VALUE_HASH_PATTERN: Final = {hash_pattern!r}")
    out.append(f'ROOT_STEP_ID: Final = "{constants["rootStep"]}"')
    out.append(f'FIRST_BRANCH_ID: Final = "{constants["firstBranch"]}"')
    out.append(f'FIRST_BRANCH_NAME: Final = "{constants["firstBranchName"]}"')
    out.append(f"HISTORY_LANE_COUNT: Final = {constants['laneCount']}")
    out.append(f'OUTSIDE_CHANGE_TYPE: Final = "{constants["outsideChangeType"]}"')
    out.append("#: The branch Clear history, keep the state makes.")
    out.append(f'FRESH_START_BRANCH_NAME: Final = "{constants["freshStartBranchName"]}"')
    out.append("#: The root step's label for each start.")
    out.append(f"ROOT_LABELS: Final[_Mapping[str, str]] = {constants['rootLabels']!r}")
    out.append("#: White space in a name, as code points: trimmed from both ends, each inner run one space.")
    codes = ", ".join(f"0x{int(code[2:], 16):04x}" for code in schema["x-names"]["whiteSpace"])
    out.append(f"NAME_WHITE_SPACE: Final = ({codes})")
    out.append("#: The state domains, in the order a full restore applies them.")
    out.append(f"DOMAINS: Final = {tuple(schema['x-domains'])!r}")
    out.append("#: What an action may declare it reads: the domains and ``history``.")
    out.append(f"READ_INPUTS: Final = {tuple(schema['x-readInputs'])!r}")
    out.append("#: The structural edits the journal keeps.")
    out.append(f"JOURNAL_OPS: Final = {tuple(schema['x-journalOps'])!r}")
    out.append("#: The tracker fields a journal entry may put back.")
    out.append(f"TRACKER_PATCH_FIELDS: Final = {tuple(schema['x-trackerPatchFields'])!r}")
    for constant, (record, key) in _ENUMS.items():
        spec = next(entry[2] for entry in records[record] if entry[0] == key)
        out.append(f"{constant}: Final = {tuple(spec['values'])!r}")
    out.append("#: The records read as maps of domain ids (``{domain: value hash}``), not as dataclasses.")
    out.append(f"MAP_RECORDS: Final = {_MAP_RECORDS!r}")
    out.append("")
    out.append("#: Every record's fields in the model's order: ``(key, attribute, spec, optional)``.")
    out.append("FIELDS: Final[_Mapping[str, tuple[tuple[str, str, _Mapping[str, Any], bool], ...]]] = {")
    for name, fields in records.items():
        out.append(f"    {name!r}: (")
        for key, attr, spec, optional, _ in fields:
            out.append(f"        ({key!r}, {attr!r}, {spec!r}, {optional}),")
        out.append("    ),")
    out.append("}")
    out.append("#: What each record is, in the model's words.")
    out.append("DESCRIPTIONS: Final[_Mapping[str, str]] = {")
    for name, definition in schema["$defs"].items():
        out.append(f"    {name!r}: {definition['description']!r},")
    out.append("}")
    out.append("#: The rules each record keeps beyond its fields' types.")
    out.append("INVARIANTS: Final[_Mapping[str, tuple[str, ...]]] = {")
    for name, definition in schema["$defs"].items():
        if "x-invariants" in definition:
            out.append(f"    {name!r}: {tuple(definition['x-invariants'])!r},")
    out.append("}")
    for name, fields in records.items():
        if name in _MAP_RECORDS:
            continue
        out.append("")
        out.append("")
        out.append("@dataclass(frozen=True)")
        out.append(f"class {name}:")
        out.append(docstring(schema["$defs"][name]["description"], "    "))
        out.append("")
        for _key, attr, spec, optional, description in fields:
            kind = annotation(spec)
            if optional:
                kind = f"{kind} | None" if not kind.endswith("| None") and kind != "Any" else kind
            out.append(f"    {attr}: {kind}")
            if description:
                out.append(docstring(description, "    "))
    out.append("")
    out.append("")
    out.append(_BUILDER)
    return "\n".join(out).rstrip() + "\n"


_BUILDER = '''#: Every record class, by its name.
RECORD_CLASSES: Final[_Mapping[str, type]] = {
    name: cls for name, cls in globals().items() if isinstance(cls, type) and name in FIELDS and name not in MAP_RECORDS
}


def build(name: str, value: Any, classes: _Mapping[str, type] | None = None) -> Any:
    """The record ``name`` made from its JSON ``value`` (already checked against
    ``FIELDS``): a dataclass of the record (or of ``classes[name]``, a
    subclass), lists as tuples, a map record as a dict. ``classes`` replaces
    the class of a record wherever it is built from here."""
    chosen = {} if classes is None else classes

    def convert(spec: _Mapping[str, Any], item: Any) -> Any:
        if item is None:
            return None
        kind = spec["type"]
        if kind == "record":
            return build(spec["name"], item, chosen)
        if kind == "array":
            return tuple(convert(spec["item"], entry) for entry in item)
        if kind == "pair":
            return tuple(convert(part, entry) for part, entry in zip(spec["items"], item, strict=True))
        return item

    fields = FIELDS[name]
    if name in MAP_RECORDS:
        specs = {key: spec for key, _, spec, _ in fields}
        return {key: convert(specs[key], item) for key, item in value.items()}
    cls = chosen.get(name, RECORD_CLASSES[name])
    return cls(**{attr: convert(spec, value.get(key)) for key, attr, spec, _ in fields})


def as_json(name: str, record: Any) -> Any:
    """A record as its JSON (the model's keys), the inverse of :func:`build`."""

    def convert(spec: _Mapping[str, Any], item: Any) -> Any:
        if item is None:
            return None
        kind = spec["type"]
        if kind == "record":
            return as_json(spec["name"], item)
        if kind in ("array", "pair"):
            specs = [spec["item"]] * len(item) if kind == "array" else spec["items"]
            return [convert(part, entry) for part, entry in zip(specs, item, strict=True)]
        return item

    fields = FIELDS[name]
    if name in MAP_RECORDS:
        specs = {key: spec for key, _, spec, _ in fields}
        return {key: convert(specs[key], item) for key, item in record.items()}
    out: dict[str, Any] = {}
    for key, attr, spec, optional in fields:
        item = getattr(record, attr)
        if optional and item is None:
            continue
        out[key] = convert(spec, item)
    return out


def replace(record: Any, **changes: Any) -> Any:
    """A copy of a record with ``changes`` (``dataclasses.replace``)."""
    return dataclasses.replace(record, **changes)'''


def main(argv: list[str]) -> int:
    check = "--check" in argv
    unknown = [argument for argument in argv if argument != "--check"]
    if unknown:
        print(f"Unknown argument {unknown[0]}; the only one is --check", file=sys.stderr)
        return 2
    text = render(SCHEMA_PATH.read_text(encoding="utf-8"))
    current = MODULE_PATH.read_text(encoding="utf-8") if MODULE_PATH.exists() else None
    if current == text:
        print(f"current {MODULE_PATH.relative_to(ROOT)}")
        return 0
    if check:
        print(f"stale   {MODULE_PATH.relative_to(ROOT)}: run python scripts/generate_history_model.py")
        return 1
    MODULE_PATH.write_text(text, encoding="utf-8")
    print(f"wrote   {MODULE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
