"""The history model Python reads sessions with, generated from the web app's.

``src/cellucid/history_model.schema.json`` is the web app's
``historyModelDocument`` (``app/history/model.js``) byte for byte, pinned at
the digest the web app and cellucid-r pin. ``src/cellucid/_history_model.py``
is written from it by ``scripts/generate_history_model.py`` and must be
current; the reading model (``cellucid.history``, ``cellucid.CellucidSession``,
the apply plan) exposes every field of every record the schema names; the
names rule and the constants are the schema's.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import importlib.util
import json
import re
import unicodedata
from importlib import resources
from pathlib import Path
from typing import Any

import pytest
from _web_root import sibling_file, web_file

from cellucid import _history_model as model
from cellucid import _history_rules as rules
from cellucid.history import STEP_COLUMNS, Branch, Group, Step, Tracker

ROOT = Path(__file__).resolve().parents[1]
#: The digest of ``history_model.schema.json``; the web app and cellucid-r pin the same one.
HISTORY_MODEL_SHA256 = "bc403ababa0a52bf62a1446b6f3598135b1990cd5a4116818ab26e3696004bc7"
REGENERATE = (
    "Run node scripts/generate-history-model.mjs in the cellucid repository (it writes the three copies and "
    "prints the digest), then python scripts/generate_history_model.py here, and pin the digest in the three suites."
)


def _schema_text() -> str:
    return resources.files("cellucid").joinpath("history_model.schema.json").read_text(encoding="utf-8")


def _generator():
    spec = importlib.util.spec_from_file_location("generate_history_model", ROOT / "scripts" / "generate_history_model.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_schema_is_pinned_and_the_web_apps() -> None:
    data = resources.files("cellucid").joinpath("history_model.schema.json").read_bytes()
    assert hashlib.sha256(data).hexdigest() == HISTORY_MODEL_SHA256, REGENERATE
    assert web_file("tests/fixtures/history-model.schema.json").read_bytes() == data, REGENERATE
    assert sibling_file("cellucid-r", "inst/history_model.schema.json").read_bytes() == data, REGENERATE


def test_the_generated_model_is_current() -> None:
    assert _generator().render(_schema_text()) == (ROOT / "src" / "cellucid" / "_history_model.py").read_text(encoding="utf-8"), (
        "src/cellucid/_history_model.py is stale: run python scripts/generate_history_model.py"
    )
    assert model.HISTORY_MODEL_SHA256 == HISTORY_MODEL_SHA256


def test_every_record_of_the_schema_has_its_fields() -> None:
    schema = json.loads(_schema_text())
    for name, definition in schema["$defs"].items():
        keys = [key for key, _, _, _ in model.FIELDS[name]]
        assert sorted(keys) == sorted(definition["properties"]), name
        assert [key for key, _, _, optional in model.FIELDS[name] if not optional] == definition["required"], name
        if name in model.MAP_RECORDS:
            continue
        cls = model.RECORD_CLASSES[name]
        assert [field.name for field in dataclasses.fields(cls)] == [attr for _, attr, _, _ in model.FIELDS[name]], name


def _snake(key: str) -> str:
    name = re.sub(r"(?<=[a-z0-9])([A-Z])", lambda match: "_" + match.group(1).lower(), key)
    return f"{name}_" if name in ("from",) else name


def test_the_reading_model_exposes_every_field_of_the_schema() -> None:
    """A field the schema has that a reader does not expose fails here."""
    from cellucid import read_session

    session = read_session(Path(__file__).parent / "fixtures" / "history_golden" / "history_golden.cellucid-session")
    tracker = session.tracker()
    objects = {
        "Session": session,
        "Tracker": tracker,
        "Branch": tracker.branches[1],
        "Step": tracker.step("s5"),
        "Group": tracker.groups[0],
        "Squash": tracker.step("s14").squashed,
        "Turn": tracker.turns[1],
        "Note": tracker.notes[0],
        "JournalEntry": tracker.journal[0],
        "Start": session.tracker("Annotation (replay)").start,
        "Scope": session.scope,
        "AppStamp": session.app,
        "ActionStamp": session.actions[0],
        "Requires": session.requires,
        "By": tracker.by,
        "Pause": session.tracker("Annotation (replay)").paused,
        "Carried": tracker.step("s16").squashed.carried,
        "Membership": tracker.journal[2].after.memberships[0],
    }
    for name, obj in objects.items():
        for key, _, _, _ in model.FIELDS[name]:
            assert hasattr(obj, _snake(key)), f"{name}.{key} is not exposed as {_snake(key)}"
    assert isinstance(tracker, Tracker) and isinstance(tracker.branches[0], Branch) and isinstance(tracker.step("s1"), Step)
    assert isinstance(tracker.groups[0], Group)
    assert list(tracker.steps.columns) == list(STEP_COLUMNS)


def test_names_are_compared_as_the_web_app_compares_them() -> None:
    assert rules.normalize_name("  Without  alpha \t") == "Without alpha"
    assert rules.normalize_name("Café") == unicodedata.normalize("NFC", "Café")
    assert rules.name_key("ÉQUIPE  One") == "équipe one"
    assert rules.name_problem("branch", "   ") == "A branch needs a name."
    assert rules.name_problem("branch", "a\u0007b") == "A branch name cannot hold control characters."
    assert rules.name_problem("branch", "B7") == 'A branch cannot be called "B7": that is the form of a branch id. Choose another name.'
    assert rules.name_problem("branch", "B2 cells") is None
    assert rules.name_problem("tracker", "T-cells") is None
    assert rules.name_problem("bookmark", "b7") is None
    assert rules.name_problem("group", "G3") == 'A group cannot be called "G3": that is the form of a group id. Choose another name.'
    assert rules.name_problem("group", "QC") is None
    # The white space set is the schema's, code point for code point.
    schema = json.loads(_schema_text())
    assert [f"U+{code:04X}" for code in model.NAME_WHITE_SPACE] == schema["x-names"]["whiteSpace"]


def test_the_constants_are_the_schemas() -> None:
    schema = json.loads(_schema_text())
    constants = schema["x-constants"]
    assert (
        constants["rootStep"], constants["firstBranch"], constants["firstBranchName"]
    ) == (model.ROOT_STEP_ID, model.FIRST_BRANCH_ID, model.FIRST_BRANCH_NAME)
    assert constants["laneCount"] == model.HISTORY_LANE_COUNT
    assert constants["outsideChangeType"] == model.OUTSIDE_CHANGE_TYPE
    assert constants["freshStartBranchName"] == model.FRESH_START_BRANCH_NAME
    assert model.STEP_KINDS == ("root", "action", "outside", "squashed")
    assert model.SHARE_CONTENTS == ("history", "state", "steps") and model.SHARE_RESULTS == ("include", "recompute")
    assert dict(model.ROOT_LABELS) == constants["rootLabels"]
    assert list(model.DOMAINS) == schema["x-domains"]
    assert list(model.JOURNAL_OPS) == schema["x-journalOps"]
    assert dict(model.ID_PATTERNS) == schema["x-ids"]


def test_a_record_that_breaks_the_model_is_refused_naming_the_path_and_the_rule() -> None:
    by = {"origin": "agent", "agent": None, "model": None, "host": None, "person": None, "turn": None}
    with pytest.raises(rules.HistoryModelError, match='by names an agent exactly when its origin is "agent"'):
        rules.check_record("By", by, "by")
        rules.check_by(by, "by")
    with pytest.raises(rules.HistoryModelError, match='by has an unknown field "colour"'):
        rules.check_record("By", {**by, "colour": 1}, "by")
    with pytest.raises(rules.HistoryModelError, match=r'branch\.lane must be a whole number from 0 to 7'):
        rules.check_record("Branch", {
            "id": "b1", "name": "X", "base": "s1", "tip": "s2", "createdAt": 0, "by": {**by, "origin": "user", "agent": None},
            "lane": 9, "locked": False, "archived": False, "notes": [],
        }, "branch")


# =============================================================================
# The rules of a paused recording, a squash and a group edit (the web app's
# tests/history-model-contract.test.mjs, on the golden history's own records)
# =============================================================================


def _golden_records() -> tuple[Any, list[dict[str, Any]]]:
    from cellucid import read_session

    session = read_session(Path(__file__).parent / "fixtures" / "history_golden" / "history_golden.cellucid-session")
    return session, copy.deepcopy(session._contents.trackers)


def _step(tracker: dict[str, Any], step_id: str) -> dict[str, Any]:
    return next(step for step in tracker["steps"] if step["id"] == step_id)


def _refused(tracker: dict[str, Any], message: str) -> None:
    with pytest.raises(rules.HistoryModelError) as refused:
        rules.check_tracker_rules(tracker, "tracker")
    assert str(refused.value) == message


TRACKER_RULES = [
    pytest.param(1, lambda t: t["paused"].__setitem__("by", []), "tracker.paused names who paused it (by)", id="paused-by-nobody"),
    pytest.param(1, lambda t: t["paused"].__setitem__("by", [t["paused"]["by"][0], t["paused"]["by"][0]]),
                 "tracker.paused.by names one author twice", id="paused-by-twice"),
    pytest.param(0, lambda t: _step(t, "s14")["squashed"]["carried"].__setitem__("bookmarks", ["Kept"]),
                 'tracker.steps["s14"].squashed carries notes or bookmarks only for squashed steps (cause steps)', id="carried-not-steps"),
    pytest.param(0, lambda t: _step(t, "s16")["squashed"]["carried"].__setitem__("notes", ["n6", "n6"]),
                 'tracker.steps["s16"].squashed.carried names a note twice', id="carried-note-twice"),
    pytest.param(0, lambda t: _step(t, "s16")["squashed"]["carried"].__setitem__("bookmarks", ["Before tidying", "BEFORE  tidying"]),
                 'tracker.steps["s16"].squashed.carried names a bookmark twice', id="carried-bookmark-twice"),
    pytest.param(0, lambda t: _step(t, "s16").__setitem__("notes", []),
                 '{where}: step "s16" carries a note it does not hold', id="carried-note-not-held"),
    pytest.param(0, lambda t: _step(t, "s16")["squashed"]["carried"].__setitem__("bookmarks", ["Markers saved"]),
                 '{where}: step "s5" and step "s16" share the name "Markers saved"', id="carried-bookmark-shared"),
    pytest.param(0, lambda t: t["journal"][6]["after"]["steps"][1]["squashed"].__setitem__("count", 3),
                 'tracker.journal["j6"]: its count adds up the steps it replaced (a squashed one counts its own)', id="squash-count"),
    pytest.param(0, lambda t: t["journal"][6]["after"]["steps"][1]["squashed"]["carried"].__setitem__("bookmarks", []),
                 'tracker.journal["j6"]: the bookmarks of the steps it replaced are carried by the squashed step', id="squash-bookmarks"),
    pytest.param(0, lambda t: t["journal"][6]["after"]["steps"][1]["squashed"]["carried"].__setitem__("notes", []),
                 'tracker.journal["j6"]: the notes of the steps it replaced are carried by the squashed step', id="squash-notes"),
    pytest.param(0, lambda t: t["journal"][6]["before"]["steps"][1].__setitem__("branch", "b1"),
                 'tracker.journal["j6"] replaces consecutive steps of one branch', id="squash-branch"),
    pytest.param(0, lambda t: t["journal"][6]["after"]["steps"][1].__setitem__("group", "g0"),
                 'tracker.journal["j6"]: a squash inside a group keeps it, and one that covers whole groups leaves none', id="squash-group"),
    pytest.param(0, lambda t: t["journal"][6]["after"]["steps"][1].__setitem__(
        "squashed", {**t["journal"][6]["after"]["steps"][1]["squashed"], "cause": "clear", "carried": {"notes": [], "bookmarks": []}}),
                 'tracker.journal["j6"] makes exactly one squashed step (cause steps)', id="squash-made"),
    pytest.param(0, lambda t: t["journal"][2]["after"]["memberships"][0].__setitem__("from", "s4"),
                 'tracker.journal["j2"].after.memberships[0] runs from a later step to an earlier one', id="membership-order"),
    pytest.param(0, lambda t: t["journal"][2]["after"]["memberships"].append(dict(t["journal"][2]["after"]["memberships"][0])),
                 'tracker.journal["j2"].after.memberships names one group twice', id="membership-twice"),
    pytest.param(0, lambda t: t["journal"][2]["after"]["steps"].append(copy.deepcopy(_step(t, "s2"))),
                 'tracker.journal["j2"] records which steps a group held (memberships), never the steps', id="group-steps"),
]


@pytest.mark.parametrize(("index", "edit", "message"), TRACKER_RULES)
def test_a_pause_a_squash_and_a_group_edit_keep_their_rules(index: int, edit: Any, message: str) -> None:
    _, trackers = _golden_records()
    tracker = trackers[index]
    rules.check_tracker_rules(tracker, "tracker")
    edit(tracker)
    _refused(tracker, message.replace("{where}", f'tracker "{tracker["id"]}"'))


def _session_record(session: Any, trackers: list[dict[str, Any]]) -> dict[str, Any]:
    manifest = session.manifest
    return {
        "createdAt": manifest["createdAt"], "app": manifest["app"], "actions": manifest["actions"], "scope": manifest["scope"],
        "datasetFingerprint": manifest["datasetFingerprint"], "source": manifest["source"], "requires": manifest["requires"],
        "cells": manifest.get("cells"), "genes": manifest.get("genes"), "activeTracker": manifest["activeTracker"],
        "trackers": trackers,
    }


def test_a_session_stamps_each_action_and_version_and_keeps_no_draft() -> None:
    session, trackers = _golden_records()
    record = _session_record(session, trackers)
    rules.validate_session(record)
    swapped = copy.deepcopy(record)
    swapped["actions"][0], swapped["actions"][1] = swapped["actions"][1], swapped["actions"][0]
    with pytest.raises(rules.HistoryModelError, match="^session.actions must be sorted by type then version, each pair once$"):
        rules.validate_session(swapped)
    extra = copy.deepcopy(record)
    extra["actions"].append({"type": "zzz.never", "version": 1, "params": None})
    with pytest.raises(rules.HistoryModelError, match="^session.actions must stamp exactly each action and version its steps use$"):
        rules.validate_session(extra)
    other_version = copy.deepcopy(record)
    other_version["actions"][0]["version"] = 2
    with pytest.raises(rules.HistoryModelError, match="^session.actions must stamp exactly each action and version its steps use$"):
        rules.validate_session(other_version)
    drafted = copy.deepcopy(record)
    drafted["trackers"][1]["paused"]["draft"] = {"colors": "d" * 64}
    with pytest.raises(rules.HistoryModelError) as refused:
        rules.validate_session(drafted)
    assert str(refused.value) == f'session: tracker "{trackers[1]["id"]}" keeps a paused draft; a file squashes it into a step first'
