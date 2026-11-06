"""
The history's rules, as the web app keeps them (``app/history/model.js``).

Every record a session file holds is checked here before anything reads it:
its shape against the model's field table (``_history_model.FIELDS``,
generated from ``history_model.schema.json``), then the rules the types cannot
say (who acted, a step's own rules by its kind, a squashed step's, a group's,
a turn's, a journal entry's, and a tracker's whole tree: one root, branches
that are lines from their fork step to their tip, groups that are runs of
consecutive steps holding whole turns, unique names, the head on the active
branch's path, every note and turn accounted for; and what a session's scope
says it holds). The messages are the web app's, word for word: a reader says
them after "This session file is damaged: ".

Names (branches, groups, trackers, bookmarks) follow the model's rule: a
name's stored form is NFC with the model's white space trimmed and each inner
run one space (:func:`normalize_name`), its key that form lower-cased then NFC
again (:func:`name_key`); a branch, group or tracker cannot be named like an
id of its kind.
"""

from __future__ import annotations

import json
import math
import re
import unicodedata
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

from . import _history_model as model
from .session_codecs import MAX_SAFE_INTEGER


class HistoryModelError(TypeError):
    """A record that breaks the history model; the message names the path and the rule."""


def _fail(message: str) -> HistoryModelError:
    return HistoryModelError(message)


_ID_PATTERNS = {kind: re.compile(pattern) for kind, pattern in model.ID_PATTERNS.items()}
_HASH = re.compile(model.VALUE_HASH_PATTERN)
_WHITE = "".join(chr(code) for code in model.NAME_WHITE_SPACE)
_WHITE_RUN = re.compile(f"[{re.escape(_WHITE)}]+")
_WHITE_EDGES = re.compile(f"^[{re.escape(_WHITE)}]+|[{re.escape(_WHITE)}]+$")
_CONTROL = re.compile("[\u0000-\u001f\u007f-\u009f]")
_ID_EXAMPLES = {
    "tracker": "t-0b7c6f0e-3c1a-4a52-9d1e-5f2b8c4a7e10", "step": "s12", "branch": "b3", "group": "g1", "note": "n4",
    "journal": "j2",
}


def id_number(identifier: str) -> int:
    """The number of an id (``s12`` -> 12)."""
    return int(identifier[1:])


def is_id(kind: str, value: Any) -> bool:
    """Whether ``value`` is an id of ``kind`` (step, branch, group, note, journal, tracker)."""
    return isinstance(value, str) and _ID_PATTERNS[kind].fullmatch(value) is not None


# ---------------------------------------------------------------------------
# Names


def normalize_name(text: str) -> str:
    """A name in its stored form: NFC, white space trimmed, each inner run one space."""
    if not isinstance(text, str):
        raise _fail("A name is text")
    return _WHITE_RUN.sub(" ", _WHITE_EDGES.sub("", unicodedata.normalize("NFC", text)))


def name_key(text: str) -> str:
    """A name's key, what uniqueness compares: its stored form lower-cased, then NFC."""
    return unicodedata.normalize("NFC", normalize_name(text).lower())


_ID_LIKE: dict[str, re.Pattern[str] | None] = {
    "branch": _ID_PATTERNS["branch"], "group": _ID_PATTERNS["group"], "tracker": _ID_PATTERNS["tracker"], "bookmark": None,
}


def name_problem(kind: str, text: str) -> str | None:
    """Why ``text`` cannot name a ``kind`` (branch, group, tracker, bookmark), or None."""
    if kind not in _ID_LIKE:
        raise _fail(f'There is no kind of name "{kind}"')
    name = normalize_name(text)
    if not name:
        return f"A {kind} needs a name."
    if _CONTROL.search(name):
        return f"A {kind} name cannot hold control characters."
    pattern = _ID_LIKE[kind]
    if pattern is not None and pattern.fullmatch(name_key(name)):
        return f'A {kind} cannot be called "{name}": that is the form of a {kind} id. Choose another name.'
    return None


def quoted_list(items: Sequence[str]) -> str:
    """``"a", "b" and "c"`` (``quotedList`` in ``model.js``)."""
    quoted = [f'"{item}"' for item in items]
    return "".join(quoted) if len(quoted) <= 1 else f"{', '.join(quoted[:-1])} and {quoted[-1]}"


def no_branch(ref: str, names: Sequence[str]) -> str:
    return f'There is no branch "{ref}"; the branches are {quoted_list(names)}.'


def no_tracker(ref: str, names: Sequence[str]) -> str:
    return f'There is no tracker "{ref}"; the trackers are {quoted_list(names)}.'


def no_step(step: str, tracker: str) -> str:
    return f'There is no step "{step}" in the tracker "{tracker}"; history_get lists its steps.'


def no_group(ref: str, names: Sequence[str]) -> str:
    if not names:
        return f'There is no group "{ref}"; the tracker has no group.'
    return f'There is no group "{ref}"; the groups are {quoted_list(names)}.'


# ---------------------------------------------------------------------------
# Shapes, from the field table


def _expected(spec: Mapping[str, Any]) -> str:
    kind = spec["type"]
    if kind == "text":
        base = "a non-empty text"
    elif kind == "string":
        base = f"text of the form {spec['pattern']}" if spec.get("pattern") else "text"
    elif kind == "name":
        base = "a name in its stored form (NFC, trimmed, single-spaced)"
    elif kind == "time":
        base = "a time in milliseconds since 1970 or null"
    elif kind == "boolean":
        base = "true or false"
    elif kind == "hash":
        base = "a value hash (64 lower-case hex digits)"
    elif kind == "id":
        base = f'a {spec["kind"]} id like "{_ID_EXAMPLES[spec["kind"]]}"'
    elif kind == "integer":
        base = (
            f"a whole number from {spec['minimum']}" if spec.get("maximum") is None
            else f"a whole number from {spec['minimum']} to {spec['maximum']}"
        )
    elif kind == "number":
        base = f"a number from {spec['minimum']}"
    elif kind == "enum":
        base = f"one of {', '.join(json.dumps(value) for value in spec['values'])}"
    elif kind == "array":
        base = "a list"
    elif kind == "pair":
        base = "a pair"
    elif kind == "record":
        base = f"a {spec['name']} record"
    elif kind == "object":
        base = "a record"
    else:
        base = "a value"
    return f"{base} or null" if spec.get("nullable") and kind != "time" else base


def _is_number(value: Any) -> bool:
    return isinstance(value, int | float) and not isinstance(value, bool)


def _is_safe_integer(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return abs(value) <= MAX_SAFE_INTEGER
    return isinstance(value, float) and math.isfinite(value) and value.is_integer() and abs(value) <= MAX_SAFE_INTEGER


def _item_path(path: str, item: Any, index: int) -> str:
    if isinstance(item, dict) and isinstance(item.get("id"), str):
        return f"{path}[{json.dumps(item['id'], ensure_ascii=False)}]"
    return f"{path}[{index}]"


def check_value(spec: Mapping[str, Any], value: Any, path: str) -> None:
    """One field's value against its spec (``checkValue``)."""
    kind = spec["type"]
    if value is None and (spec.get("nullable") or kind in ("time", "json")):
        return

    def wrong() -> HistoryModelError:
        return _fail(f"{path} must be {_expected(spec)}")

    if kind == "text":
        if not isinstance(value, str) or not value:
            raise wrong()
    elif kind == "string":
        pattern = spec.get("pattern")
        if not isinstance(value, str) or (pattern and re.search(pattern, value) is None):
            raise wrong()
    elif kind == "name":
        if not isinstance(value, str) or not value or normalize_name(value) != value:
            raise wrong()
    elif kind == "time":
        if not _is_number(value) or not math.isfinite(value):
            raise wrong()
    elif kind == "boolean":
        if not isinstance(value, bool):
            raise wrong()
    elif kind == "hash":
        if not isinstance(value, str) or _HASH.fullmatch(value) is None:
            raise wrong()
    elif kind == "id":
        if not is_id(spec["kind"], value):
            raise wrong()
    elif kind == "integer":
        maximum = spec.get("maximum")
        if not _is_safe_integer(value) or value < spec["minimum"] or (maximum is not None and value > maximum):
            raise wrong()
    elif kind == "number":
        if not _is_number(value) or not math.isfinite(value) or value < spec["minimum"]:
            raise wrong()
    elif kind == "enum":
        if isinstance(value, bool) or value not in spec["values"]:
            raise wrong()
    elif kind == "array":
        if not isinstance(value, list):
            raise wrong()
        for index, item in enumerate(value):
            check_value(spec["item"], item, _item_path(path, item, index))
    elif kind == "pair":
        if not isinstance(value, list) or len(value) != 2:
            raise wrong()
        for index, (item_spec, item) in enumerate(zip(spec["items"], value, strict=True)):
            check_value(item_spec, item, f"{path}[{index}]")
    elif kind == "record":
        check_record(spec["name"], value, path)
    elif kind == "object":
        if not isinstance(value, dict):
            raise wrong()
    elif kind == "json":
        return
    else:
        raise _fail(f"{path} has an unknown spec")


def check_record(name: str, value: Any, path: str) -> None:
    """A record's shape: exactly its fields, each its type (``checkRecord``)."""
    fields = model.FIELDS[name]
    if not isinstance(value, dict):
        raise _fail(f"{path} must be a {name} record")
    known = {key for key, _, _, _ in fields}
    for key in value:
        if key not in known:
            raise _fail(f'{path} has an unknown field "{key}"')
    for key, _, spec, optional in fields:
        if key not in value:
            if optional:
                continue
            raise _fail(f'{path} lacks "{key}"')
        check_value(spec, value[key], f"{path}.{key}")


# ---------------------------------------------------------------------------
# Record rules


def check_by(by: Mapping[str, Any], path: str) -> None:
    origin = by["origin"]
    if (origin == "agent") != (by["agent"] is not None):
        raise _fail(f'{path} names an agent exactly when its origin is "agent"')
    if by["model"] is not None and origin != "agent":
        raise _fail(f"{path} names a model only for an agent")
    if (origin == "host") != (by["host"] is not None):
        raise _fail(f'{path} names a host exactly when its origin is "host"')
    if by["turn"] is not None and origin not in ("agent", "host"):
        raise _fail(f"{path} belongs to a turn only for an agent or a host")


def same_by(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    return all(a[key] == b[key] for key, _, _, _ in model.FIELDS["By"])


def _check_note(note: Mapping[str, Any], path: str) -> None:
    check_by(note["by"], f"{path}.by")
    if note["fromReply"] and note["by"]["origin"] != "agent":
        raise _fail(f"{path} is from an agent's reply only when an agent is its author")


def _check_notes(notes: Sequence[Mapping[str, Any]], path: str) -> None:
    for note in notes:
        _check_note(note, f"{path}[{json.dumps(note['id'])}]")


def _same_keys(record: Mapping[str, Any], keys: Iterable[str]) -> bool:
    return sorted(record) == sorted(keys)


def check_step_rules(step: Mapping[str, Any], path: str) -> None:
    """One step's own rules, whatever tree it is in (``checkStepRules``)."""
    check_by(step["by"], f"{path}.by")
    _check_notes(step["notes"], f"{path}.notes")
    number = id_number(step["id"])
    if len(set(step["children"])) != len(step["children"]):
        raise _fail(f"{path} lists a child twice")
    if any(id_number(child) <= number for child in step["children"]):
        raise _fail(f"{path} has a child numbered before it")
    if step["at"] is not None and step["until"] is not None and step["until"] < step["at"]:
        raise _fail(f"{path} ends (until) before it starts (at)")
    if step["status"] == "unapplicable" and step["reason"] is None:
        raise _fail(f"{path} no longer applies and must say why (reason)")
    kind = step["kind"]
    if (step["parent"] is None) != (kind == "root"):
        raise _fail(f'{path} is the root (kind "root") exactly when it has no parent')
    if (kind == "squashed") != (step["squashed"] is not None):
        raise _fail(f"{path} says what it stands for (squashed) exactly when it is squashed")
    if kind != "action" and step["v"] is not None:
        raise _fail(f"{path} records a parameters version (v) only for an action")
    if step["parent"] is None:
        if step["id"] != model.ROOT_STEP_ID:
            raise _fail(f"{path} has no parent, and only the root ({model.ROOT_STEP_ID}) has none")
        if (
            step["type"] is not None or step["params"] is not None or step["code"] is not None or step["touches"]
            or step["before"] or step["reads"] or step["outcome"] is not None or step["coalesceKey"] is not None
            or step["replayOf"] is not None or step["merged"] != 1 or step["ms"] is not None
            or step["status"] != "applied" or step["reason"] is not None or step["group"] is not None
        ):
            raise _fail(f"{path} is the root and records no action")
        return
    if id_number(step["parent"]) >= number:
        raise _fail(f"{path} is numbered before its parent")
    if (kind == "outside") != (step["type"] == model.OUTSIDE_CHANGE_TYPE):
        raise _fail(
            f'{path} is of kind "outside" exactly when it is a change made outside an action ({model.OUTSIDE_CHANGE_TYPE})'
        )
    if (step["type"] == model.OUTSIDE_CHANGE_TYPE) != (step["by"]["origin"] == "outside"):
        raise _fail(
            f'{path} is by "outside" exactly when it is a change made outside an action ({model.OUTSIDE_CHANGE_TYPE})'
        )
    if kind == "squashed":
        if (
            step["type"] is not None or step["params"] is not None or step["code"] is not None
            or step["coalesceKey"] is not None or step["replayOf"] is not None
        ):
            raise _fail(
                f"{path} is squashed and records no action of its own (type, params, code, coalesceKey and replayOf null)"
            )
        check_squash_rules(step["squashed"], f"{path}.squashed")
    else:
        if step["type"] is None:
            raise _fail(f"{path} names no action")
        if step["params"] is None:
            raise _fail(f"{path} records an action and its parameters (params)")
    if kind == "action" and step["v"] is None:
        raise _fail(f"{path} records the version of its action's parameters (v)")
    if len(set(step["touches"])) != len(step["touches"]):
        raise _fail(f"{path} lists a touched domain twice")
    ordered = [domain for domain in model.DOMAINS if domain in step["touches"]]
    if ordered != list(step["touches"]):
        raise _fail(f"{path} lists its touched domains out of domain order")
    if not _same_keys(step["before"], step["touches"]) or not _same_keys(step["after"], step["touches"]):
        raise _fail(f"{path} must hold a value before and after for exactly the domains it touches")


def check_squash_rules(squash: Mapping[str, Any], path: str) -> None:
    """What a squashed step stands for (``checkSquashRules``)."""
    counted = squash["cause"] in ("steps", "clear")
    if counted != (squash["count"] is not None):
        raise _fail(f"{path} counts the steps it replaced (count) exactly when they were recorded and shown (steps, clear)")
    if (squash["cause"] == "share") != (len(squash["by"]) == 0):
        raise _fail(f"{path} names who acted (by) exactly unless it is a shared state")
    _check_distinct_bys(squash["by"], f"{path}.by")
    notes, bookmarks = squash["carried"]["notes"], squash["carried"]["bookmarks"]
    if squash["cause"] != "steps" and (notes or bookmarks):
        raise _fail(f"{path} carries notes or bookmarks only for squashed steps (cause steps)")
    if len(set(notes)) != len(notes):
        raise _fail(f"{path}.carried names a note twice")
    if len({name_key(name) for name in bookmarks}) != len(bookmarks):
        raise _fail(f"{path}.carried names a bookmark twice")


def _check_distinct_bys(bys: Sequence[Mapping[str, Any]], path: str) -> None:
    for index, by in enumerate(bys):
        check_by(by, f"{path}[{index}]")
    for index in range(1, len(bys)):
        if any(same_by(other, bys[index]) for other in bys[:index]):
            raise _fail(f"{path} names one author twice")


def _check_pause_rules(pause: Mapping[str, Any], path: str) -> None:
    """A tracker's paused recording (``checkPauseRules``)."""
    if not pause["by"]:
        raise _fail(f"{path} names who paused it (by)")
    _check_distinct_bys(pause["by"], f"{path}.by")


def _step_count(step: Mapping[str, Any]) -> int:
    """How many recorded steps a step stands for: a squashed step its count (one when none was recorded), any other one."""
    if step["kind"] == "squashed" and step["squashed"]["count"] is not None:
        return int(step["squashed"]["count"])
    return 1


def _check_squash_entry(entry: Mapping[str, Any], path: str) -> None:
    """A squash's journal entry held to the squash rules (``checkSquashEntry``): the
    squashed step takes the id of the range's last step; the range is one line
    of one branch; its count adds up what it replaced; the replaced steps'
    notes and bookmarks are carried; it stays inside one group or leaves none."""
    was = {step["id"]: step for step in entry["before"]["steps"]}
    made = [
        step for step in entry["after"]["steps"]
        if step["kind"] == "squashed" and step["squashed"]["cause"] == "steps" and (
            step["id"] not in was or was[step["id"]]["kind"] != "squashed" or was[step["id"]]["squashed"] != step["squashed"]
        )
    ]
    if len(made) != 1:
        raise _fail(f"{path} makes exactly one squashed step (cause steps)")
    squashed = made[0]
    kept = {step["id"] for step in entry["after"]["steps"]}
    replaced = [step for step in entry["before"]["steps"] if step["id"] not in kept or step["id"] == squashed["id"]]
    last = next((step for step in replaced if step["id"] == squashed["id"]), None)
    if last is None or any(id_number(step["id"]) > id_number(squashed["id"]) for step in replaced):
        raise _fail(f"{path}: the squashed step takes the id of the last step it replaced")
    ordered = sorted(replaced, key=lambda step: id_number(step["id"]))
    if any(step["branch"] != squashed["branch"] or (index > 0 and step["parent"] != ordered[index - 1]["id"])
           for index, step in enumerate(ordered)):
        raise _fail(f"{path} replaces consecutive steps of one branch")
    if squashed["parent"] != ordered[0]["parent"] or list(squashed["children"]) != list(last["children"]):
        raise _fail(f"{path}: the squashed step keeps the first replaced step's parent and the last one's children")
    if squashed["squashed"]["count"] != sum(_step_count(step) for step in replaced):
        raise _fail(f"{path}: its count adds up the steps it replaced (a squashed one counts its own)")
    notes = [
        note for step in replaced
        for note in [*(item["id"] for item in step["notes"]), *(step["squashed"]["carried"]["notes"] if step["kind"] == "squashed" else [])]
    ]
    bookmarks = [
        name for step in replaced
        for name in [*([] if step["bookmark"] is None else [step["bookmark"]]),
                     *(step["squashed"]["carried"]["bookmarks"] if step["kind"] == "squashed" else [])]
    ]
    own = {note["id"] for note in squashed["notes"]}
    if any(note not in own for note in notes) or sorted(set(notes)) != sorted(squashed["squashed"]["carried"]["notes"]):
        raise _fail(f"{path}: the notes of the steps it replaced are carried by the squashed step")
    if sorted(bookmarks) != sorted(squashed["squashed"]["carried"]["bookmarks"]):
        raise _fail(f"{path}: the bookmarks of the steps it replaced are carried by the squashed step")
    groups = {step["group"] for step in replaced}
    if squashed["group"] is not None and groups != {squashed["group"]}:
        raise _fail(f"{path}: a squash inside a group keeps it, and one that covers whole groups leaves none")


def check_group_rules(group: Mapping[str, Any], path: str) -> None:
    """A group's own rules (``checkGroupRules``)."""
    check_by(group["by"], f"{path}.by")
    _check_notes(group["notes"], f"{path}.notes")


def _check_turn_rules(turn: Mapping[str, Any], path: str) -> None:
    check_by(turn["by"], f"{path}.by")
    if turn["by"]["origin"] not in ("agent", "host"):
        raise _fail(f"{path} is an agent's or a notebook's")
    if turn["by"]["turn"] != turn["id"]:
        raise _fail(f"{path}.by.turn must be the turn's own id")
    if turn["startedAt"] is not None and turn["endedAt"] is not None and turn["endedAt"] < turn["startedAt"]:
        raise _fail(f"{path} ends before it starts")
    _check_notes(turn["notes"], f"{path}.notes")


def _check_start_rules(start: Mapping[str, Any], path: str) -> None:
    replay = start["kind"] == "replay"
    linked = replay or start["kind"] == "import"
    if linked != (start["from"] is not None):
        raise _fail(f"{path} names where it comes from (from) exactly for a replay or an import")
    if replay != (start["report"] is not None):
        raise _fail(f"{path} keeps a report exactly for a replay")
    if not replay and start["mapping"] is not None:
        raise _fail(f"{path} keeps a name mapping only for a replay")
    if start["kind"] == "import" and start["from"]["file"] is None:
        raise _fail(f"{path} of an import names its file")


def _check_patch_rules(patch: Mapping[str, Any], path: str) -> None:
    for step in patch["steps"]:
        check_step_rules(step, f"{path}.steps[{json.dumps(step['id'])}]")
    for branch in patch["branches"]:
        check_by(branch["by"], f"{path}.branches[{json.dumps(branch['id'])}].by")
        _check_notes(branch["notes"], f"{path}.branches[{json.dumps(branch['id'])}].notes")
    for group in patch["groups"]:
        check_group_rules(group, f"{path}.groups[{json.dumps(group['id'])}]")
    for index, membership in enumerate(patch["memberships"]):
        if id_number(membership["from"]) > id_number(membership["to"]):
            raise _fail(f"{path}.memberships[{index}] runs from a later step to an earlier one")
    if len({membership["group"] for membership in patch["memberships"]}) != len(patch["memberships"]):
        raise _fail(f"{path}.memberships names one group twice")
    for turn in patch["turns"]:
        _check_turn_rules(turn, f"{path}.turns[{json.dumps(turn['id'])}]")
    for kind in ("steps", "branches", "groups", "turns"):
        ids = [item["id"] for item in patch[kind]]
        if len(set(ids)) != len(ids):
            raise _fail(f"{path}.{kind} holds one record twice")


def _check_journal_rules(entry: Mapping[str, Any], path: str) -> None:
    check_by(entry["by"], f"{path}.by")
    if entry["purged"] != (entry["before"] is None) or entry["purged"] != (entry["after"] is None):
        raise _fail(f"{path} keeps what it replaced and made (before, after) exactly until it is purged")
    if entry["purged"]:
        return
    _check_patch_rules(entry["before"], f"{path}.before")
    _check_patch_rules(entry["after"], f"{path}.after")
    if not _same_keys(entry["before"]["tracker"], entry["after"]["tracker"]):
        raise _fail(f"{path} must replace the same tracker fields it puts back")
    if entry["op"].startswith("group.") and (entry["before"]["steps"] or entry["after"]["steps"]):
        raise _fail(f"{path} records which steps a group held (memberships), never the steps")
    if entry["op"] == "squash":
        _check_squash_entry(entry, path)


def _in_id_order(items: Sequence[Mapping[str, Any]]) -> bool:
    return all(id_number(items[index - 1]["id"]) < id_number(items[index]["id"]) for index in range(1, len(items)))


def _check_unique_names(
    items: Iterable[Mapping[str, Any]],
    kind: str,
    label: Callable[[Mapping[str, Any]], str],
    name_of: Callable[[Mapping[str, Any]], str | None],
    path: str,
) -> None:
    seen: dict[str, Mapping[str, Any]] = {}
    for item in items:
        name = name_of(item)
        if name is None:
            continue
        problem = name_problem(kind, name)
        if problem is not None:
            raise _fail(f"{path}: {problem}")
        key = name_key(name)
        if key in seen:
            raise _fail(f'{path}: {label(seen[key])} and {label(item)} share the name "{name}"')
        seen[key] = item


def validate_tracker(tracker: Mapping[str, Any], path: str = "tracker") -> Mapping[str, Any]:
    """Check one tracker whole: every record's shape and rules, and the tree's (``validateTracker``)."""
    check_record("Tracker", tracker, path)
    return check_tracker_rules(tracker, path)


def check_tracker_rules(tracker: Mapping[str, Any], path: str) -> Mapping[str, Any]:
    """A tracker's rules beyond its records' shapes (``checkTrackerRules``)."""
    where = f'{path} "{tracker["id"]}"'
    check_by(tracker["by"], f"{path}.by")
    _check_start_rules(tracker["start"], f"{path}.start")
    if tracker["removed"] is not None:
        check_by(tracker["removed"]["by"], f"{path}.removed.by")
    _check_notes(tracker["notes"], f"{path}.notes")

    steps: dict[str, Mapping[str, Any]] = {}
    for step in tracker["steps"]:
        if step["id"] in steps:
            raise _fail(f'{where} holds the step "{step["id"]}" twice')
        steps[step["id"]] = step
        check_step_rules(step, f"{path}.steps[{json.dumps(step['id'])}]")
        if id_number(step["id"]) >= tracker["nextStep"]:
            raise _fail(f'{where}: step "{step["id"]}" is numbered past nextStep')
    if not _in_id_order(tracker["steps"]):
        raise _fail(f"{where} lists its steps out of id order")
    if tracker["rootId"] != model.ROOT_STEP_ID:
        raise _fail(f"{where}'s root must be {model.ROOT_STEP_ID}")
    root = steps.get(model.ROOT_STEP_ID)
    if root is None:
        raise _fail(f"{where} has no root step {model.ROOT_STEP_ID}")
    if not same_by(root["by"], tracker["by"]):
        raise _fail(f"{where}'s root is by whoever started it (tracker.by)")
    for step in tracker["steps"]:
        if step["parent"] is not None:
            parent = steps.get(step["parent"])
            if parent is None or step["id"] not in parent["children"]:
                raise _fail(f'{where}: step "{step["id"]}" and its parent disagree')
        for child in step["children"]:
            found = steps.get(child)
            if found is None or found["parent"] != step["id"]:
                raise _fail(f'{where}: step "{step["id"]}" lists a child that is not its own')

    branches: dict[str, Mapping[str, Any]] = {}
    for branch in tracker["branches"]:
        if branch["id"] in branches:
            raise _fail(f'{where} holds the branch "{branch["id"]}" twice')
        branches[branch["id"]] = branch
        check_by(branch["by"], f"{path}.branches[{json.dumps(branch['id'])}].by")
        _check_notes(branch["notes"], f"{path}.branches[{json.dumps(branch['id'])}].notes")
        if id_number(branch["id"]) >= tracker["nextBranch"]:
            raise _fail(f'{where}: branch "{branch["id"]}" is numbered past nextBranch')
    if not _in_id_order(tracker["branches"]):
        raise _fail(f"{where} lists its branches out of id order")
    _check_unique_names(tracker["branches"], "branch", lambda item: f'branch "{item["id"]}"', lambda item: item["name"], where)
    own_steps: dict[str, list[Mapping[str, Any]]] = {branch_id: [] for branch_id in branches}
    for step in tracker["steps"]:
        own = own_steps.get(step["branch"])
        if own is None:
            raise _fail(f'{where}: step "{step["id"]}" is on the branch "{step["branch"]}", which it does not hold')
        own.append(step)
    root_branches = [branch for branch in tracker["branches"] if branch["base"] is None]
    if len(root_branches) != 1:
        raise _fail(f"{where} must have exactly one branch without a fork step (the one holding the root)")
    if root["branch"] != root_branches[0]["id"]:
        raise _fail(f"{where}: the root must be on the branch without a fork step")
    for branch in tracker["branches"]:
        label = f'{where}: branch "{branch["id"]}"'
        own = own_steps[branch["id"]]
        if not own:
            raise _fail(f"{label} has no step of its own")
        if branch["base"] is not None:
            base = steps.get(branch["base"])
            if base is None:
                raise _fail(f"{label} forks from a step the tracker does not hold")
            if base["branch"] == branch["id"]:
                raise _fail(f"{label} forks from one of its own steps")
        firsts = [
            step for step in own if step["parent"] is None or steps[step["parent"]]["branch"] != branch["id"]
        ]
        if len(firsts) != 1 or firsts[0]["parent"] != branch["base"]:
            raise _fail(f"{label}'s steps must form one line starting right after its fork step")
        count = 1
        current = firsts[0]
        while True:
            following = [child for child in current["children"] if steps[child]["branch"] == branch["id"]]
            if len(following) > 1:
                raise _fail(f'{label} splits at step "{current["id"]}"; a branch is one line')
            if not following:
                break
            current = steps[following[0]]
            count += 1
        if count != len(own):
            raise _fail(f"{label}'s steps must form one line starting right after its fork step")
        if current["id"] != branch["tip"]:
            raise _fail(f'{label}\'s tip must be its last step ("{current["id"]}")')

    head = steps.get(tracker["headId"])
    if head is None:
        raise _fail(f"{where}'s head is not one of its steps")
    active = branches.get(tracker["activeBranch"])
    if active is None:
        raise _fail(f"{where}'s active branch is not one of its branches")
    on_path = False
    walker: Mapping[str, Any] | None = steps[active["tip"]]
    while walker is not None:
        if walker["id"] == head["id"]:
            on_path = True
            break
        walker = None if walker["parent"] is None else steps.get(walker["parent"])
    if not on_path:
        raise _fail(f"{where}'s head must be on its active branch's path")
    if active["archived"] and head["branch"] != active["id"]:
        raise _fail(f"{where}'s active branch is archived, so the head must be on it")
    main = branches.get(tracker["mainBranch"])
    if main is None:
        raise _fail(f"{where}'s main branch is not one of its branches")
    if main["archived"]:
        raise _fail(f"{where}'s main branch cannot be archived")
    if tracker["paused"] is not None:
        _check_pause_rules(tracker["paused"], f"{path}.paused")

    groups: dict[str, list[Mapping[str, Any]]] = {}
    for group in tracker["groups"]:
        if group["id"] in groups:
            raise _fail(f'{where} holds the group "{group["id"]}" twice')
        groups[group["id"]] = []
        check_group_rules(group, f"{path}.groups[{json.dumps(group['id'])}]")
        if id_number(group["id"]) >= tracker["nextGroup"]:
            raise _fail(f'{where}: group "{group["id"]}" is numbered past nextGroup')
    if not _in_id_order(tracker["groups"]):
        raise _fail(f"{where} lists its groups out of id order")
    _check_unique_names(tracker["groups"], "group", lambda item: f'group "{item["id"]}"', lambda item: item["name"], where)
    for step in tracker["steps"]:
        if step["group"] is None:
            continue
        members = groups.get(step["group"])
        if members is None:
            raise _fail(f'{where}: step "{step["id"]}" is in the group "{step["group"]}", which it does not hold')
        members.append(step)
    for group_id, members in groups.items():
        if not members:
            raise _fail(f'{where}: group "{group_id}" has no step')
        # Steps are in id order, and along one branch a child's id is above its parent's.
        for index in range(1, len(members)):
            if members[index]["branch"] != members[0]["branch"] or members[index]["parent"] != members[index - 1]["id"]:
                raise _fail(f"{where}: group \"{group_id}\"'s steps must be consecutive steps of one branch")
    # A group never cuts a run of one turn's consecutive steps on its branch: a turn that
    # acts again after someone else's step starts a new run, which the group need not hold.
    for step in tracker["steps"]:
        turn_id = step["by"]["turn"]
        if turn_id is None:
            continue
        parent = None if step["parent"] is None else steps.get(step["parent"])
        if parent is not None and parent["branch"] == step["branch"] and parent["by"]["turn"] == turn_id and parent["group"] != step["group"]:
            raise _fail(
                f'{where}: a group cuts the run of turn "{turn_id}" between steps "{parent["id"]}" and "{step["id"]}"; '
                "a group holds whole runs of a turn"
            )

    bookmarked = [
        {"id": step["id"], "name": name} for step in tracker["steps"]
        for name in [*([] if step["bookmark"] is None else [step["bookmark"]]),
                     *(step["squashed"]["carried"]["bookmarks"] if step["kind"] == "squashed" else [])]
    ]
    _check_unique_names(bookmarked, "bookmark", lambda item: f'step "{item["id"]}"', lambda item: item["name"], where)
    for step in tracker["steps"]:
        if step["kind"] != "squashed":
            continue
        held_notes = {note["id"] for note in step["notes"]}
        if any(note not in held_notes for note in step["squashed"]["carried"]["notes"]):
            raise _fail(f'{where}: step "{step["id"]}" carries a note it does not hold')
    note_ids: set[str] = set()
    every_note = [
        *tracker["notes"],
        *(note for branch in tracker["branches"] for note in branch["notes"]),
        *(note for group in tracker["groups"] for note in group["notes"]),
        *(note for step in tracker["steps"] for note in step["notes"]),
        *(note for turn in tracker["turns"] for note in turn["notes"]),
    ]
    for note in every_note:
        if note["id"] in note_ids:
            raise _fail(f'{where} holds the note "{note["id"]}" twice')
        note_ids.add(note["id"])
        if id_number(note["id"]) >= tracker["nextNote"]:
            raise _fail(f'{where}: note "{note["id"]}" is numbered past nextNote')
    turns: dict[str, Mapping[str, Any]] = {}
    for turn in tracker["turns"]:
        if turn["id"] in turns:
            raise _fail(f'{where} holds the turn "{turn["id"]}" twice')
        turns[turn["id"]] = turn
        _check_turn_rules(turn, f"{path}.turns[{json.dumps(turn['id'])}]")
        if turn["branch"] not in branches:
            raise _fail(f'{where}: turn "{turn["id"]}" acts on a branch the tracker does not hold')
    turn_steps: set[str] = set()
    for step in tracker["steps"]:
        turn_id = step["by"]["turn"]
        if turn_id is None:
            continue
        if turn_id not in turns:
            raise _fail(f'{where}: step "{step["id"]}" belongs to the turn "{turn_id}", which it does not hold')
        turn_steps.add(turn_id)
    for turn in tracker["turns"]:
        # A turn that only acted outside the state (a vote published) keeps that on record.
        if turn["id"] not in turn_steps and len(turn["outside"]) == 0:
            raise _fail(f'{where}: turn "{turn["id"]}" has no step and did nothing outside the state')

    journal_ids: set[str] = set()
    for entry in tracker["journal"]:
        if entry["id"] in journal_ids:
            raise _fail(f'{where} holds the journal entry "{entry["id"]}" twice')
        journal_ids.add(entry["id"])
        if id_number(entry["id"]) >= tracker["nextJournal"]:
            raise _fail(f'{where}: journal entry "{entry["id"]}" is numbered past nextJournal')
        _check_journal_rules(entry, f"{path}.journal[{json.dumps(entry['id'])}]")
    if not _in_id_order(tracker["journal"]):
        raise _fail(f"{where} lists its journal out of id order")
    if tracker["removed"] is not None and not any(
        entry["op"] == "tracker.remove" and not entry["purged"] for entry in tracker["journal"]
    ):
        raise _fail(f"{where} is removed, so its journal must keep the removal (to restore it)")
    return tracker


def _check_trackers(trackers: Sequence[Mapping[str, Any]], active_tracker: str, path: str) -> None:
    ids: set[str] = set()
    for tracker in trackers:
        check_tracker_rules(tracker, "tracker")
        if tracker["id"] in ids:
            raise _fail(f'{path} holds the tracker "{tracker["id"]}" twice')
        ids.add(tracker["id"])
    _check_unique_names(
        [tracker for tracker in trackers if tracker["removed"] is None], "tracker",
        lambda item: f'tracker "{item["id"]}"', lambda item: item["name"], path,
    )
    active = next((tracker for tracker in trackers if tracker["id"] == active_tracker), None)
    if active is None:
        raise _fail(f"{path}'s active tracker is not one of its trackers")
    if active["removed"] is not None:
        raise _fail(f"{path}'s active tracker cannot be removed")


def validate_session(session: Mapping[str, Any]) -> Mapping[str, Any]:
    """Check a session file as read: its header records, every tracker, and
    that its scope names exactly the trackers it holds (``validateSession``)."""
    check_record("Session", session, "session")
    _check_trackers(session["trackers"], session["activeTracker"], "session")
    held = [tracker["id"] for tracker in session["trackers"]]
    if list(session["scope"]["trackers"]) != held:
        raise _fail("session.scope must name exactly the trackers it holds, in order")
    for entry in session["scope"]["branches"]:
        tracker = next((candidate for candidate in session["trackers"] if candidate["id"] == entry["tracker"]), None)
        if tracker is None:
            raise _fail(f'session.scope names branches of a tracker it does not hold ("{entry["tracker"]}")')
        known = {branch["id"] for branch in tracker["branches"]}
        if any(branch not in known for branch in entry["branches"]):
            raise _fail(f'session.scope names a branch the tracker "{entry["tracker"]}" does not hold')
    stamps = [(stamp["type"], stamp["version"]) for stamp in session["actions"]]
    if any(index > 0 and stamps[index - 1] >= stamps[index] for index in range(len(stamps))):
        raise _fail("session.actions must be sorted by type then version, each pair once")
    used: set[str] = set()
    for tracker in session["trackers"]:
        every_step = [
            *tracker["steps"],
            *(step for entry in tracker["journal"] if not entry["purged"] for step in [*entry["before"]["steps"], *entry["after"]["steps"]]),
        ]
        used.update(f'{step["type"]}@{step["v"]}' for step in every_step if step["kind"] == "action")
        if tracker["paused"] is not None and tracker["paused"]["draft"] is not None:
            raise _fail(f'session: tracker "{tracker["id"]}" keeps a paused draft; a file squashes it into a step first')
    stamped = [f"{kind}@{version}" for kind, version in stamps]
    if len(stamped) != len(used) or any(stamp not in used for stamp in stamped):
        raise _fail("session.actions must stamp exactly each action and version its steps use")
    _check_scope_holds(session["scope"], session["trackers"])
    return session


def _check_scope_holds(scope: Mapping[str, Any], trackers: Sequence[Mapping[str, Any]]) -> None:
    """What the share chose is what the file holds (``checkScopeHolds``)."""
    if scope["content"] == "steps" and scope["results"] != "recompute":
        raise _fail("session.scope of steps only holds no results (results is recompute)")
    for tracker in trackers:
        where = f'session: tracker "{tracker["id"]}"'
        if not scope["chats"] and tracker["chats"]:
            raise _fail(f"{where} holds conversations its scope leaves out (chats)")
        if not scope["journal"] and tracker["journal"]:
            raise _fail(f"{where} holds Recently deleted, which its scope leaves out (journal)")
        if not scope["notes"]:
            noted = any(
                item["notes"]
                for item in [tracker, *tracker["branches"], *tracker["groups"], *tracker["steps"], *tracker["turns"]]
            )
            if noted or any(step["bookmark"] is not None for step in tracker["steps"]):
                raise _fail(f"{where} holds notes or bookmarks its scope leaves out (notes)")
    if scope["content"] != "state":
        return
    if len(trackers) != 1 or len(trackers[0]["branches"]) != 1:
        raise _fail("session.scope of a state holds one tracker of one branch")
    tracker = trackers[0]
    steps = tracker["steps"]
    if (
        len(steps) != 2 or steps[1]["kind"] != "squashed" or steps[1]["squashed"]["cause"] != "share"
        or tracker["headId"] != steps[1]["id"]
    ):
        raise _fail("session.scope of a state holds the root and one squashed step (cause share), shown")
    if (
        tracker["turns"] or tracker["groups"] or tracker["notes"] or steps[0]["notes"] or steps[1]["notes"]
        or tracker["journal"] or tracker["chats"]
    ):
        raise _fail("session.scope of a state holds no turn, group, note, journal or conversation")
    # Privacy: a state shared without its history names nobody (no agent, model, notebook, person or turn).
    bys = [tracker["by"], steps[0]["by"], steps[1]["by"], *([] if tracker["paused"] is None else tracker["paused"]["by"])]
    if not all(
        by["origin"] in ("user", "startup") and all(by[key] is None for key in ("agent", "model", "host", "person", "turn"))
        for by in bys
    ):
        raise _fail(
            "session.scope of a state names nobody: every by is the person or the opening, "
            "with no agent, model, notebook, name or turn"
        )
