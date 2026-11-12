"""
A Cellucid history, read: trackers, branches, groups, steps, notes, turns and the journal.

Every record is the history model's (``history_model.schema.json``, the web
app's ``app/history/model.js``), one frozen dataclass each, generated
(``_history_model.py``): a field of the model is an attribute here in snake
case (``headId`` is ``head_id``; ``from`` is ``from_``). Four records do more
than hold their fields:

- :class:`Tracker`: one tree of steps on the dataset. ``tracker.steps`` is
  every step as a table (one row per step, in step order: the order of their
  ids, parents before children, as the web app and R list them),
  ``tracker.step("s12")`` one step, ``tracker.branch("Main")`` one branch (by
  its name or its id), ``tracker.methods()`` and ``tracker.pseudocode()`` a
  branch as the web app writes it, ``tracker.stale()`` the results whose
  inputs changed since they were made;
- :class:`Branch`: a line of steps from its fork step (``base``) to its tip;
  ``branch.steps`` is its path from the root, as a table;
- :class:`Group`: a named run of consecutive steps of one branch (a folder:
  "QC", "Clustering"); ``group.steps`` is its steps, as a table; each step
  names its group (``step.group``);
- :class:`Step`: one step as recorded, of a ``kind``: ``root``, ``action``
  (its ``type``, ``params`` and their version ``v``), ``outside`` (a change
  made outside an action) or ``squashed`` (several steps merged into one, or
  what changed while recording was paused: no action of its own, and
  ``squashed`` says why and who acted); ``step.state()`` is every state
  domain's value after it, ``step.diff(other)`` what differs between two
  steps, domain by domain.

A step shows its pseudocode as it was recorded (``code``): parameters of an
earlier version (``v`` below the action's version in this package's
catalog) are shown and applied as they were; the web app upgrades them only
to replay them. No picture of a step is ever stored.

The session reader (:class:`cellucid.CellucidSession`) builds them; a live
viewer's history is read the same way (``viewer.history``).
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, fields
from typing import Any, cast

import pandas as pd

from . import _history_model as model
from . import _history_rules as _rules
from ._history_model import (
    ActionStamp,
    ApplyPlanRow,
    AppStamp,
    BranchName,
    BranchTip,
    By,
    Carried,
    CategoryMap,
    CellOrder,
    CellOrderEmbedding,
    DatasetIdentity,
    DatasetName,
    JournalEntry,
    Link,
    Membership,
    NameMap,
    Note,
    Patch,
    Pause,
    Removal,
    Renamed,
    ReplayEntry,
    ReplayOf,
    Requires,
    Scope,
    ScopeBranches,
    Skipped,
    Squash,
    Start,
    TrackerName,
    TrackerPatch,
    Turn,
)
from ._history_model import Mapping as NameMapping
from ._history_rules import HistoryModelError, name_key, normalize_name

__all__ = [
    "STEP_COLUMNS",
    "NOTE_COLUMNS",
    "DOMAIN_NAMES",
    "ActionStamp",
    "AppStamp",
    "ApplyPlanRow",
    "Branch",
    "BranchName",
    "BranchTip",
    "By",
    "Carried",
    "CategoryMap",
    "CellOrder",
    "CellOrderEmbedding",
    "DatasetIdentity",
    "DatasetName",
    "Diff",
    "Group",
    "HistoryModelError",
    "JournalEntry",
    "Link",
    "Membership",
    "NameMap",
    "NameMapping",
    "Note",
    "Patch",
    "Pause",
    "Removal",
    "Renamed",
    "ReplayEntry",
    "ReplayOf",
    "Requires",
    "Scope",
    "ScopeBranches",
    "Skipped",
    "Squash",
    "Start",
    "Step",
    "Tracker",
    "TrackerName",
    "TrackerPatch",
    "Turn",
    "build_tracker",
    "methods_paragraph",
    "name_key",
    "normalize_name",
]

#: The columns of a tracker's step table, in order (the same in R).
STEP_COLUMNS = (
    "id", "branch", "parent", "kind", "type", "v", "label", "code", "origin", "agent", "model", "host", "turn", "at",
    "status", "stale", "bookmark", "group", "n_notes",
)
#: The columns of a tracker's note table, in order (the same in R).
NOTE_COLUMNS = ("id", "on", "target", "text", "origin", "agent", "at", "edited_at", "from_reply")

#: Each state domain's name, as a person reads it (``DOMAINS`` in the web
#: app's ``app/actions/domains.js``).
DOMAIN_NAMES: Mapping[str, str] = {
    "cells": "Open cells",
    "mapping": "Data mapping",
    "expression": "Expression values",
    "field-overlays": "Fields and labels",
    "gene-sets": "Gene sets",
    "views": "Views",
    "active-field": "Active field",
    "filters": "Filters",
    "colors": "Colors",
    "visibility": "Hidden cells",
    "highlights": "Highlights",
    "camera": "Camera",
    "camera-path": "Camera path",
    "settings": "Settings",
    "analysis": "Analysis",
    "layout": "Layout",
}


# ---------------------------------------------------------------------------
# Words, as the web app writes them


def format_count(value: int) -> str:
    """``1200`` -> ``1,200`` (``formatCount``)."""
    return f"{int(value):,}"


def format_count_of(value: int, singular: str, plural: str | None = None) -> str:
    """``"1 cell"``, ``"1,200 cells"`` (``formatCountOf``)."""
    return f"{format_count(value)} {count_noun(value, singular, plural)}"


def count_noun(value: int, singular: str, plural: str | None = None) -> str:
    return singular if value == 1 else (plural if plural is not None else f"{singular}s")


def _js_json(value: Any) -> str:
    """``JSON.stringify`` of plain data, keys in their order (what the web compares)."""
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _same_value(a: Any, b: Any) -> bool:
    """``sameValue`` (``state/domains/same-value.js``): typed arrays by their bytes."""
    array_a = hasattr(a, "digest") and hasattr(a, "data")
    array_b = hasattr(b, "digest") and hasattr(b, "data")
    if array_a or array_b:
        return array_a and array_b and a.digest == b.digest
    if isinstance(a, float) and isinstance(b, float):
        return (math.isnan(a) and math.isnan(b)) or (a == b and math.copysign(1, a) == math.copysign(1, b))
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same_value(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(_same_value(a[key], b[key]) for key in a)
    if isinstance(a, list | dict) or isinstance(b, list | dict):
        return False
    return type(a) is type(b) and a == b if isinstance(a, bool) or isinstance(b, bool) else a == b


def _cells_of(value: Any) -> int:
    data = getattr(value, "data", value)
    return int(len(data))


def _describe_cells(a: Any, b: Any) -> list[str]:
    def opened(open_cells: Any) -> str:
        if open_cells is None:
            return "every cell"
        return f"{format_count_of(open_cells['k'], 'cell')} ({open_cells['method']}, seed {open_cells['seed']})"

    def genes(count: Any) -> str:
        return "the opening's sample" if count is None else format_count_of(count, "cell")

    lines = []
    if _js_json(a["open"]) != _js_json(b["open"]):
        lines.append(f"Cells open: {opened(a['open'])} → {opened(b['open'])}")
    if a["genes"] != b["genes"]:
        lines.append(f"Genes read for: {genes(a['genes'])} → {genes(b['genes'])}")
    return lines


def _describe_user_field(before: Any, after: Any) -> list[str]:
    if before is None:
        return [f'{"Annotation" if after["kind"] == "category" else "Field"} "{after["key"]}" made']
    if after is None:
        return [f'Field "{before["key"]}" removed']
    lines = []
    if before["key"] != after["key"]:
        lines.append(f'Field "{before["key"]}" renamed "{after["key"]}"')
    if before["deleted"] != after["deleted"]:
        lines.append(f'Field "{after["key"]}" {"deleted" if after["deleted"] else "restored"}')
    if before["purged"] != after["purged"] and after["purged"]:
        lines.append(f'Field "{after["key"]}" deleted for good')
    if not _same_value(before.get("categories"), after.get("categories")):
        lines.append(f'Field "{after["key"]}": categories changed')
    elif not _same_value(before.get("codes"), after.get("codes")) or not _same_value(before.get("values"), after.get("values")):
        lines.append(f'Field "{after["key"]}": values changed')
    return lines


def _describe_field_overlays(a: Any, b: Any) -> list[str]:
    lines: list[str] = []
    renamed_a = {f"{entry['source']}:{entry['key']}": entry["name"] for entry in a["renamedFields"]}
    renamed_b = {f"{entry['source']}:{entry['key']}": entry["name"] for entry in b["renamedFields"]}
    for key, name in renamed_b.items():
        if renamed_a.get(key) != name:
            lines.append(f'Field "{key[key.index(":") + 1:]}" shown as "{name}"')
    for key in renamed_a:
        if key not in renamed_b:
            lines.append(f'Field "{key[key.index(":") + 1:]}" shown under its own name')

    def category_id(entry: Any) -> str:
        return _js_json([entry["source"], entry["field"], entry["category"]])

    categories_a = {category_id(entry): entry for entry in a["renamedCategories"]}
    categories_b = {category_id(entry): entry for entry in b["renamedCategories"]}
    for key, entry in categories_b.items():
        was = categories_a.get(key)
        if was is None or was["name"] != entry["name"]:
            lines.append(f'{entry["field"]}: "{_js_text(entry["category"])}" shown as "{_js_text(entry["name"])}"')
    for key, entry in categories_a.items():
        if key not in categories_b:
            lines.append(f'{entry["field"]}: "{_js_text(entry["category"])}" shown under its own name')
    deleted_a = {f"{entry['source']}:{entry['key']}": entry for entry in a["deletedFields"]}
    deleted_b = {f"{entry['source']}:{entry['key']}": entry for entry in b["deletedFields"]}
    for key, entry in deleted_b.items():
        if key not in deleted_a:
            lines.append(f'Field "{entry["key"]}" deleted')
    for key, entry in deleted_a.items():
        if key not in deleted_b:
            lines.append(f'Field "{entry["key"]}" restored')
    user_a = {item["id"]: item for item in a["userFields"]}
    user_b = {item["id"]: item for item in b["userFields"]}
    for key in dict.fromkeys([*user_a, *user_b]):
        before, after = user_a.get(key), user_b.get(key)
        if not _same_value(before, after):
            lines.extend(_describe_user_field(before, after))
    return lines


def _js_text(value: Any) -> str:
    """A label as JavaScript writes it in a template (``1`` stays ``1``, ``True`` is ``true``)."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _describe_gene_sets(a: Any, b: Any) -> list[str]:
    before = {entry["name"]: entry for entry in a["sets"]}
    after = {entry["name"]: entry for entry in b["sets"]}
    lines = []
    for name, entry in after.items():
        was = before.get(name)
        if was is None:
            lines.append(f"Gene set {name} added ({format_count_of(len(entry['genes']), 'gene')})")
        elif _js_json(was) != _js_json(entry):
            lines.append(f"Gene set {name} changed ({format_count_of(len(entry['genes']), 'gene')})")
    for name in before:
        if name not in after:
            lines.append(f"Gene set {name} removed")
    return lines


def _describe_highlights(a: Any, b: Any) -> list[str]:
    lines = []
    pages_a = {page["id"]: page for page in a["pages"]}
    pages_b = {page["id"]: page for page in b["pages"]}
    for page_id, page in pages_b.items():
        was = pages_a.get(page_id)
        if was is None:
            lines.append(f'Highlight page "{page["name"]}" added')
            continue
        if was["name"] != page["name"]:
            lines.append(f'Highlight page "{was["name"]}" renamed "{page["name"]}"')
        if was["color"] != page["color"]:
            lines.append(f'Highlight page "{page["name"]}" recolored')
        groups_a = {group["id"]: group for group in was["groups"]}
        groups_b = {group["id"]: group for group in page["groups"]}
        for group_id, group in groups_b.items():
            before = groups_a.get(group_id)
            if before is None:
                lines.append(f'"{group["label"]}" highlighted ({format_count_of(_cells_of(group["cells"]), "cell")}) on "{page["name"]}"')
            elif before["enabled"] != group["enabled"]:
                lines.append(f'"{group["label"]}" {"shown" if group["enabled"] else "hidden"} on "{page["name"]}"')
            elif before["label"] != group["label"]:
                lines.append(f'"{before["label"]}" renamed "{group["label"]}"')
            elif _cells_of(before["cells"]) != _cells_of(group["cells"]):
                lines.append(f'"{group["label"]}" now {format_count_of(_cells_of(group["cells"]), "cell")}')
        for group_id, group in groups_a.items():
            if group_id not in groups_b:
                lines.append(f'"{group["label"]}" removed from "{page["name"]}"')
    for page_id, page in pages_a.items():
        if page_id not in pages_b:
            lines.append(f'Highlight page "{page["name"]}" removed')
    if a["active"] != b["active"] and b["active"] in pages_b:
        lines.append(f'Highlight page "{pages_b[b["active"]]["name"]}" shown')
    return lines


#: The domains Python and R describe in the web app's words; every other one is said by its name.
_DESCRIBE: Mapping[str, Callable[[Any, Any], list[str]]] = {
    "cells": _describe_cells,
    "field-overlays": _describe_field_overlays,
    "gene-sets": _describe_gene_sets,
    "highlights": _describe_highlights,
}


def describe_domain(domain: str, a: Any, b: Any) -> list[str]:
    """What differs in one domain's value, as lines: the web app's words for
    the open cells, fields and labels, gene sets and highlights; ``"<name>
    changed"`` for any other domain (whose words need the page)."""
    lines = _DESCRIBE[domain](a, b) if domain in _DESCRIBE else []
    return lines if lines else [f"{DOMAIN_NAMES[domain]} changed"]


# ---------------------------------------------------------------------------
# Who took a step, the methods paragraph and the pseudocode


_AGENT_KEYS = {"claude": "claude", "chatgpt": "chatgpt"}
_HOST_NAMES = {"python": "Python", "r": "R"}
#: Domains whose change a methods section reports (``METHOD_DOMAINS``).
METHOD_DOMAINS = frozenset({"mapping", "expression", "field-overlays", "gene-sets", "active-field", "filters", "visibility", "highlights", "analysis"})
_SUBSAMPLE_METHODS = {
    "blocks": "in blocks of neighboring cells",
    "uniform": "uniformly at random",
    "stratified": "stratified by a field",
    "spatial": "spatially",
}


def origin_of(by: By) -> tuple[str, str]:
    """Who acted (a ``By`` record: ``step.by``, ``note.by``): a key for its look
    and a short name (``originOf``)."""
    if by.origin == "outside":
        return "outside", "Not recorded"
    if by.origin == "user":
        return "person", "Person"
    if by.origin == "agent":
        name = by.agent or "Agent"
        return _AGENT_KEYS.get(name.lower(), "agent"), name
    if by.origin == "host":
        host = (by.host or "").lower()
        return (host if host in _HOST_NAMES else "host"), _HOST_NAMES.get(host, "Notebook")
    if by.origin == "session":
        return "session", "Session"
    if by.origin == "replay":
        return "replay", "Replay"
    return "start", "Start"


def _origin_phrase(by: By) -> str | None:
    key, name = origin_of(by)
    if key == "outside":
        return None
    if key == "person":
        return "the user"
    if key in ("claude", "chatgpt", "agent"):
        return f"an AI assistant ({name}{f', {by.model}' if by.model else ''}) at the user's request"
    if key in ("python", "r"):
        return f"a {name} notebook through the Cellucid viewer API"
    if key == "host":
        return "a notebook through the Cellucid viewer API"
    return {"session": "loading a session file", "replay": "replaying a recorded history"}.get(key, "opening the dataset")


def _actors_of(step: model.Step) -> Sequence[By]:
    """Everyone who acted in a step: a squashed step's ``squashed.by`` when it names anyone, else its ``by``."""
    if step.kind == "squashed" and step.squashed is not None and step.squashed.by:
        return step.squashed.by
    return (step.by,)


def squash_phrase(squashed: Squash) -> str:
    """What a squashed step stands for, as a sentence opening (``squashPhrase``)."""
    if squashed.cause == "paused":
        return "Changes made with recording paused"
    if squashed.cause == "steps":
        return f"{format_count_of(squashed.count or 0, 'step')} merged into one"
    if squashed.cause == "clear":
        return "A fresh start from the state reached before the history was cleared"
    if squashed.cause == "share":
        return "A shared state, received without its history"
    raise ValueError(f"A squashed step is made for paused, steps, clear or share, not {squashed.cause}")


def _group_runs(steps: Sequence[model.Step], groups: Sequence[model.Group]) -> list[tuple[model.Group | None, list[model.Step]]]:
    """Consecutive steps that share a group, in order (``groupRuns``)."""
    by_id = {group.id: group for group in groups}
    runs: list[tuple[model.Group | None, list[model.Step]]] = []
    for step in steps:
        if runs and (runs[-1][0].id if runs[-1][0] is not None else None) == step.group:
            runs[-1][1].append(step)
            continue
        if step.group is not None and step.group not in by_id:
            raise ValueError(f"Step {step.id[1:]} names the group \"{step.group}\", which the tracker does not hold")
        runs.append((by_id[step.group] if step.group is not None else None, [step]))
    return runs


def _run_headings(runs: Sequence[tuple[model.Group | None, Sequence[model.Step]]]) -> list[str | None]:
    """Each run's heading (``runHeadings``): its group's name, or "After <group>" /
    "Before <group>" for steps in no group; none on a path with no group."""
    headings: list[str | None] = []
    for index, (group, _steps) in enumerate(runs):
        if group is not None:
            headings.append(group.name)
        elif index > 0:
            headings.append(f"After {cast(model.Group, runs[index - 1][0]).name}")
        else:
            headings.append(f"Before {cast(model.Group, runs[1][0]).name}" if len(runs) > 1 else None)
    return headings


def _trim_end(text: str) -> str:
    return re.sub(r"[.;,\s]+$", "", str(text).strip())


def _in_sentence(text: str) -> str:
    if len(text) < 2:
        return text
    first = re.split(r"\s", text, maxsplit=1)[0]
    if len(first) > 1 and first == first.upper():
        return text
    return text[0].lower() + text[1:]


def number_ranges(numbers: Sequence[int]) -> str:
    """``[1, 2, 3, 5]`` -> ``"1–3 and 5"`` (``numberRanges``)."""
    runs: list[list[int]] = []
    for number in sorted(set(numbers)):
        if runs and number == runs[-1][1] + 1:
            runs[-1][1] = number
        else:
            runs.append([number, number])
    parts = [str(a) if a == b else f"{a}–{b}" for a, b in runs]
    return "".join(parts) if len(parts) <= 1 else f"{', '.join(parts[:-1])} and {parts[-1]}"


def _js_round(value: float) -> int:
    return math.floor(value + 0.5)


def _dataset_sentence(dataset: Mapping[str, Any]) -> str:
    from .session import _js_number

    name = f"“{_trim_end(dataset['name'])}”" if dataset.get("name") else None
    size = " and ".join(
        part for part in (
            format_count_of(dataset["cells"], "cell") if dataset.get("cells") is not None else None,
            format_count_of(dataset["genes"], "gene") if dataset.get("genes") is not None else None,
        ) if part
    )
    citation = _trim_end(dataset["citation"]) if dataset.get("citation") else None
    facts = "; ".join(part for part in (size or None, citation) if part)
    subject = "The dataset" if name is None else f"The dataset {name}"
    open_cells = dataset.get("open_cells")
    how = "was opened in full"
    if open_cells is not None and open_cells["open"] < open_cells["of"]:
        percent = 100 * open_cells["open"] / open_cells["of"]
        # The web app's formatShare: whole from 10 (never 100 short of every
        # cell), two significant figures below.
        share = f"{min(99, _js_round(percent))}%" if percent >= 10 else f"{_js_number(float(f'{percent:.2g}'))}%"
        method = _SUBSAMPLE_METHODS.get(open_cells.get("method") or "")
        seed = f" with seed {open_cells['seed']}" if isinstance(open_cells.get("seed"), int) else ""
        how = (
            f"was opened on {format_count(open_cells['open'])} of its {format_count_of(open_cells['of'], 'cell')} ({share})"
            + (f", chosen {method}{seed}" if method else seed)
        )
    return f"{subject}{f' ({facts})' if facts else ''} {how}."


def _step_item(step: model.Step) -> str:
    if step.kind == "outside":
        return f"{', '.join(DOMAIN_NAMES[d].lower() for d in step.touches if d in METHOD_DOMAINS)} changed outside the recorded steps"
    if step.kind == "squashed":
        domains = ", ".join(DOMAIN_NAMES[d].lower() for d in step.touches)
        return f"{_in_sentence(squash_phrase(cast(Squash, step.squashed)))}{f' ({domains})' if domains else ''}"
    return _in_sentence(_trim_end(step.label))


def _origin_sentence(numbered: Sequence[tuple[model.Step, int]]) -> str | None:
    by_origin: dict[str, list[int]] = {}
    for step, number in numbered:
        for by in _actors_of(step):
            phrase = None if origin_of(by)[0] == "person" else _origin_phrase(by)
            if phrase is None:
                continue
            numbers = by_origin.setdefault(phrase, [])
            if number not in numbers:
                numbers.append(number)
    clauses = []
    for position, (phrase, numbers) in enumerate(by_origin.items()):
        noun = count_noun(len(numbers), "step")
        if position > 0:
            clauses.append(f"{noun} {number_ranges(numbers)} by {phrase}")
        else:
            clauses.append(f"{noun[0].upper()}{noun[1:]} {number_ranges(numbers)} {'was' if len(numbers) == 1 else 'were'} taken by {phrase}")
    if not clauses:
        return None
    if len(clauses) == 1:
        return f"{clauses[0]}."
    return f"{', '.join(clauses[:-1])}, and {clauses[-1]}."


def _list_sentences(numbered: Sequence[tuple[model.Step, int]]) -> list[str]:
    items = "; ".join(f"({number}) {_step_item(step)}" for step, number in numbered)
    n = len(numbered)
    sentences = [f"The following {format_count_of(n, 'step')} {'was' if n == 1 else 'were'} applied{'' if n == 1 else ', in this order'}: {items}."]
    origins = _origin_sentence(numbered)
    if origins is not None:
        sentences.append(origins)
    return sentences


def methods_paragraph(
    steps: Sequence[model.Step], *, dataset: Mapping[str, Any], version: str | None, notes: Sequence[Note] = (),
    groups: Sequence[model.Group] = (),
) -> str:
    """A branch's path of steps as methods text, as the web app writes it
    (``app/history/exports/methods.js``): the tracker's first note (its
    abstract) as a paragraph of its own; the software and its version, the
    dataset, each step that changed the data or what was computed (a change
    made outside an action, a squashed step said as what it stands for), who
    took them when it was not the person, and what is not listed and why.
    When a listed step is in a group, each run of listed steps follows under
    its heading (the group's name, or "Before"/"After" a group)."""
    taken = [step for step in steps if step.kind != "root"]
    flagged = [step for step in taken if step.kind == "action" and step.status == "unapplicable"]
    applied = [step for step in taken if step not in flagged and step.touches]
    listed = [step for step in applied if any(domain in METHOD_DOMAINS for domain in step.touches)]
    presentation = len(applied) - len(listed)
    numbered = [(step, index + 1) for index, step in enumerate(listed)]
    opening = [
        f"Single-cell data were explored interactively in Cellucid{f' (version {version})' if version else ''}.",
        _dataset_sentence(dataset),
    ]
    if not listed:
        opening.append("No recorded step changed the data, its annotation or what was computed.")
    closing = []
    squashed = [number for step, number in numbered if step.kind == "squashed"]
    if len(squashed) == 1:
        closing.append(f"For step {squashed[0]}, only the state it reached was kept, not the operations that led to it.")
    elif squashed:
        closing.append(f"For steps {number_ranges(squashed)}, only the states they reached were kept, not the operations that led to them.")
    if presentation > 0:
        closing.append(
            f"{format_count_of(presentation, 'further step')} only changed the camera, the panel layout, colors or display settings"
            f" and {'is' if presentation == 1 else 'are'} not listed."
        )
    if flagged:
        closing.append(
            f"{format_count_of(len(flagged), 'recorded step')} no longer applied after an earlier step was removed"
            f" and {'is' if len(flagged) == 1 else 'are'} not included."
        )
    closing.append("The complete history, with every parameter, can be reopened from the accompanying Cellucid session file.")

    number_of = {step.id: number for step, number in numbered}
    runs = _group_runs(listed, groups)
    headings = _run_headings(runs)

    def sentences_of(run_steps: Sequence[model.Step]) -> list[str]:
        return _list_sentences([(step, number_of[step.id]) for step in run_steps])

    if all(heading is None for heading in headings):
        paragraphs = [" ".join([*opening, *(s for _group, run in runs for s in sentences_of(run)), *closing])]
    else:
        summary = "The step applied is listed below." if len(listed) == 1 else f"The {format_count_of(len(listed), 'step')} applied are listed below, in order."
        paragraphs = [" ".join([*opening, summary, *closing])]
        paragraphs += [f"{heading}\n{' '.join(sentences_of(run))}" for heading, (_group, run) in zip(headings, runs, strict=True)]
    return "\n\n".join([notes[0].text, *paragraphs] if notes else paragraphs)


def _comment(text: str) -> str:
    return re.sub(r"\s+", " ", str(text)).strip()


def _pseudocode_line(step: model.Step) -> str:
    names = ", ".join(DOMAIN_NAMES[domain] for domain in step.touches)
    if step.kind == "outside":
        return f"# Changed outside the history{f': {names}' if names else ''}"
    if step.kind == "squashed":
        return f"# {squash_phrase(cast(Squash, step.squashed))}{f': {names}' if names else ''} (values only)"
    if step.status == "unapplicable":
        return f"# No longer applies ({_comment(step.reason or 'not applied')}): {step.code}"
    key, name = origin_of(step.by)
    return f"{step.code}" if key == "person" else f"{step.code}  # {name}"


def pseudocode_text(
    steps: Sequence[model.Step], *, dataset_name: str | None, version: str | None, tracker: str, branch: str,
    groups: Sequence[model.Group] = (),
) -> str:
    """A branch's path of steps as pseudocode, from each step's stored line, by
    the web app's rules (``app/history/exports/pseudocode.js``): agents and
    notebooks named in a trailing comment, ``# --- QC ---`` and
    ``# --- end of QC ---`` around a group's steps, ``# Bookmark: <name>``
    after a bookmarked step (and, after a squashed step, each bookmark it
    kept from the steps it replaced, "(on a step merged into this one)"), and a flagged step, a change made outside an
    action and a squashed step (values only) as comments."""
    taken = [step for step in steps if step.kind != "root"]
    app = f"Cellucid {version}" if version else "Cellucid"
    lines = [
        f"# {app}{f' · {_comment(dataset_name)}' if dataset_name else ''}",
        f'# Tracker "{_comment(tracker)}", branch "{_comment(branch)}": {format_count_of(len(taken), "step")}',
    ]
    for group, run in _group_runs(taken, groups):
        name = _comment(group.name) if group is not None else None
        if name is not None:
            lines.append(f"# --- {name} ---")
        for step in run:
            lines.append(_pseudocode_line(step))
            if step.bookmark is not None:
                lines.append(f"# Bookmark: {_comment(step.bookmark)}")
            if step.kind == "squashed" and step.squashed is not None:
                lines.extend(f"# Bookmark: {_comment(name)} (on a step merged into this one)" for name in step.squashed.carried.bookmarks)
        if name is not None:
            lines.append(f"# --- end of {name} ---")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The reading model


def _plain_value(value: Any) -> Any:
    """A value as this package hands it out: typed arrays as read-only numpy arrays."""
    if hasattr(value, "digest") and hasattr(value, "data"):
        return value.data
    if isinstance(value, list):
        return [_plain_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain_value(item) for key, item in value.items()}
    return value


#: The attributes the reading model keeps beside a record's fields (set once, when built).
_OWNER = "_owner"
_READING = "_reading"
_BY_ID = "_steps_by_id"
_STALE = "_stale_cache"


@dataclass(frozen=True)
class _Context:
    """What a tracker reads its values from: the session's values by hash, and the dataset's facts."""

    value: Callable[[str], Any]
    values_included: bool
    dataset: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Step(model.Step):
    """One step as recorded: the root, an action, a change made outside an
    action, or a squashed step (``kind``).

    Every field of the model's ``Step``; ``state()`` is the value of every
    state domain after it, ``diff(other)`` what differs between it and
    another step, ``stale`` why its result no longer matches its inputs (or
    None).
    """

    def _tracker(self) -> Tracker:
        return cast(Tracker, getattr(self, _OWNER))

    def state(self) -> dict[str, Any]:
        """The value of every state domain after this step: ``{domain: value}``.

        A domain's value at a step is what the nearest step on the way from
        the root last left it (the root holds them all). Values are plain
        data, as ``app/state/domains/README.md`` shapes them: typed arrays are
        read-only numpy arrays of the matching type (``Uint32Array`` cells are
        ``uint32``, ``Float32Array`` values ``float32``), a number JSON cannot
        hold is a float (``nan``, ``inf``, ``-0.0``).
        """
        tracker = self._tracker()
        return {domain: _plain_value(value) for domain, value in tracker._raw_state(self.id).items()}

    @property
    def stale(self) -> dict[str, Any] | None:
        """Why this step's result no longer matches its inputs at the tracker's
        head, or None: ``{"cause": step id or None, "domains": [...]}``."""
        return self._tracker()._stale().get(self.id)

    def diff(self, other: Step | str) -> Diff:
        """What differs between the state after this step and after ``other``
        (a step of any tracker, or the id of one of this tracker's)."""
        tracker = self._tracker()
        target = tracker.step(other) if isinstance(other, str) else other
        return Diff.between(self, target)



@dataclass(frozen=True)
class Branch(model.Branch):
    """A line of steps with a unique name, from its fork step (``base``) to its tip.

    Every field of the model's ``Branch``; ``main`` says whether it is its
    tracker's main branch, ``steps`` is its path from the root to its tip as
    a table, ``own_steps`` the steps made on it.
    """

    def _tracker(self) -> Tracker:
        return cast(Tracker, getattr(self, _OWNER))

    @property
    def main(self) -> bool:
        """Whether this is the tracker's main branch (what sharing, replay and load default to)."""
        return self._tracker().main_branch == self.id

    @property
    def path(self) -> tuple[Step, ...]:
        """Its steps from the root to its tip, root first."""
        return self._tracker().path_to(self.tip)

    @property
    def steps(self) -> pd.DataFrame:
        """Its path from the root to its tip, one row per step (``STEP_COLUMNS``)."""
        return self._tracker()._table(self.path)

    @property
    def own_steps(self) -> pd.DataFrame:
        """The steps made on this branch, from its first to its tip."""
        tracker = self._tracker()
        return tracker._table([step for step in tracker.step_records if step.branch == self.id])


@dataclass(frozen=True)
class Group(model.Group):
    """A named run of consecutive steps of one branch: a folder over them,
    nothing merged.

    Every field of the model's ``Group``; ``step_ids`` are its steps (those
    that name it), ``first`` and ``last`` its ends, ``branch`` the branch they
    are on, ``steps`` them as a table.
    """

    def _tracker(self) -> Tracker:
        return cast(Tracker, getattr(self, _OWNER))

    @property
    def step_ids(self) -> tuple[str, ...]:
        """Its steps' ids, in step order."""
        return tuple(step.id for step in self._tracker().step_records if step.group == self.id)

    @property
    def first(self) -> str:
        """Its first step."""
        return self.step_ids[0]

    @property
    def last(self) -> str:
        """Its last step."""
        return self.step_ids[-1]

    @property
    def branch(self) -> str:
        """The branch its steps are on (its id)."""
        return self._tracker().step(self.first).branch

    @property
    def steps(self) -> pd.DataFrame:
        """Its steps, one row each (``STEP_COLUMNS``)."""
        tracker = self._tracker()
        return tracker._table([step for step in tracker.step_records if step.group == self.id])


@dataclass(frozen=True)
class Tracker(model.Tracker):
    """One tree of steps on the dataset, with its branches, turns, notes and journal.

    Every field of the model's ``Tracker`` (``step_records`` holds the
    :class:`Step` records, ``groups`` the :class:`Group` records, ``paused``
    ``None`` while it records, else a :class:`Pause`: since when it was paused
    (``at``), who paused it and everyone who acted while it was (``by``), and
    ``draft``, which a file keeps ``None``); a squashed step's ``squashed``
    has ``carried``, the notes and bookmarks of the steps it replaced
    (:class:`Carried`); ``steps`` is the table of every step, ``step(id)``,
    ``branch(name or id)`` and ``group(name or id)`` one record,
    ``notes_table()`` every note, ``stale()`` the stale results, ``methods()``
    and ``pseudocode()`` a branch in words.
    """

    def _context(self) -> _Context:
        return cast(_Context, getattr(self, _READING))

    # -- records ------------------------------------------------------------

    @property
    def steps(self) -> pd.DataFrame:
        """Every step, one row each, in step order (``STEP_COLUMNS``): its id,
        branch (id), parent, ``kind`` (root, action, outside, squashed), action
        (``type``) and its parameters' version (``v``), label, pseudocode line
        (``code``), who took it (``origin``, ``agent``, ``model``, ``host``,
        ``turn``), when (``at``, a UTC timestamp), whether it still applies
        (``status``), whether its result is stale at the head, its bookmark,
        its group (id) and its number of notes."""
        return self._table(self.step_records)

    def step(self, step: str) -> Step:
        """One step, by its id (``s12``)."""
        found = self._by_id().get(step) if isinstance(step, str) else None
        if found is None:
            raise ValueError(_rules.no_step(str(step), self.name))
        return found

    def branch(self, branch: str | None = None) -> Branch:
        """One branch, by its name (compared by its key) or its id; ``None`` is the active branch."""
        if branch is None:
            branch = self.active_branch
        if not isinstance(branch, str):
            raise TypeError("branch must be a branch's name or id (a string) or None")
        found = next((item for item in self.branches if item.id == branch), None)
        if found is None:
            key = name_key(branch)
            found = next((item for item in self.branches if name_key(item.name) == key), None)
        if found is None:
            raise ValueError(_rules.no_branch(branch, [item.name for item in self.branches]))
        return cast(Branch, found)

    def group(self, group: str) -> Group:
        """One group, by its name (compared by its key) or its id."""
        if not isinstance(group, str):
            raise TypeError("group must be a group's name or id (a string)")
        found = next((item for item in self.groups if item.id == group), None)
        if found is None:
            key = name_key(group)
            found = next((item for item in self.groups if name_key(item.name) == key), None)
        if found is None:
            raise ValueError(_rules.no_group(group, [item.name for item in self.groups]))
        return cast(Group, found)

    def turn(self, turn: str) -> Turn:
        """One turn, by its id."""
        found = next((item for item in self.turns if item.id == turn), None)
        if found is None:
            raise ValueError(f'There is no turn "{turn}" in the tracker "{self.name}".')
        return found

    def path_to(self, step: str) -> tuple[Step, ...]:
        """The steps from the root to ``step``, root first."""
        by_id = self._by_id()
        path: list[Step] = []
        current: Step | None = self.step(step)
        while current is not None:
            path.append(current)
            current = None if current.parent is None else by_id[current.parent]
        return tuple(reversed(path))

    def notes_table(self) -> pd.DataFrame:
        """Every note of the tracker, its branches, groups, steps and turns, in
        note order (``NOTE_COLUMNS``): what it is on (``on``: tracker, branch,
        group, step or turn, and ``target``, its id), its text, its author, when it
        was written and last edited, and whether it was taken from an agent's
        reply."""
        rows = []
        targets = [
            *(("tracker", self.id, note) for note in self.notes),
            *(("branch", branch.id, note) for branch in self.branches for note in branch.notes),
            *(("group", group.id, note) for group in self.groups for note in group.notes),
            *(("step", step.id, note) for step in self.step_records for note in step.notes),
            *(("turn", turn.id, note) for turn in self.turns for note in turn.notes),
        ]
        for on, target, note in sorted(targets, key=lambda entry: _rules.id_number(entry[2].id)):
            rows.append([
                note.id, on, target, note.text, note.by.origin, note.by.agent, _timestamp(note.at),
                _timestamp(note.edited_at), note.from_reply,
            ])
        return pd.DataFrame(rows, columns=list(NOTE_COLUMNS))

    # -- states -------------------------------------------------------------

    def _by_id(self) -> dict[str, Step]:
        cached: dict[str, Step] | None = getattr(self, _BY_ID, None)
        if cached is None:
            cached = {step.id: cast(Step, step) for step in self.step_records}
            object.__setattr__(self, _BY_ID, cached)
        return cached

    def _hashes(self, step: str) -> dict[str, str]:
        hashes: dict[str, str] = {}
        for item in self.path_to(step):
            hashes.update(item.after)
        return hashes

    def _raw_state(self, step: str) -> dict[str, Any]:
        context = self._context()
        if not context.values_included:
            raise ValueError(
                "This session was shared as its steps only, without the values they leave behind: "
                "replay it in Cellucid to see its states."
            )
        return {domain: context.value(value_hash) for domain, value_hash in self._hashes(step).items()}

    def _stale(self) -> dict[str, dict[str, Any]]:
        cached: dict[str, dict[str, Any]] | None = getattr(self, _STALE, None)
        if cached is not None:
            return cached
        path = self.path_to(self.head_id)
        head = self._hashes(self.head_id)
        out: dict[str, dict[str, Any]] = {}
        for index, step in enumerate(path):
            read = {domain: value for domain, value in step.reads.items() if domain != "history"}
            if not read or step.status != "applied":
                continue
            if not all(head.get(domain) == step.after[domain] for domain in step.touches):
                continue
            changed = [domain for domain, value_hash in read.items() if head.get(domain) != value_hash]
            if not changed:
                continue
            cause = next((later for later in path[index + 1 :] if any(domain in changed for domain in later.touches)), None)
            out[step.id] = {
                "cause": None if cause is None else cause.id,
                "domains": [domain for domain in model.DOMAINS if domain in changed],
            }
        object.__setattr__(self, _STALE, out)
        return out

    def stale(self) -> list[dict[str, Any]]:
        """The results stale at the head (``history-v2`` §3.4): each step on
        the head's path whose output is still what the head shows while an
        input it read changed since, with the first later step that changed
        it (``cause``) and the inputs (``domains``), in step order."""
        return [{"step": step, **reason} for step, reason in self._stale().items()]

    # -- words --------------------------------------------------------------

    def methods(self, branch: str | None = None) -> str:
        """The methods text of a branch (main when left out), as the web app writes it:
        the tracker's first note opens it, and group names head its sections."""
        target = self.branch(branch if branch is not None else self.main_branch)
        context = self._context()
        return methods_paragraph(
            target.path, dataset=context.dataset, version=context.dataset.get("version"), notes=self.notes,
            groups=self.groups,
        )

    def pseudocode(self, branch: str | None = None) -> str:
        """A branch's steps (main when left out) as pseudocode, one stored line each,
        groups marked with ``# --- QC ---`` comments."""
        target = self.branch(branch if branch is not None else self.main_branch)
        context = self._context()
        return pseudocode_text(
            target.path, dataset_name=context.dataset.get("name"), version=context.dataset.get("version"),
            tracker=self.name, branch=target.name, groups=self.groups,
        )

    # -- tables -------------------------------------------------------------

    def _table(self, steps: Sequence[model.Step]) -> pd.DataFrame:
        stale = self._stale()
        rows = [
            [
                step.id, step.branch, step.parent, step.kind, step.type, step.v, step.label, step.code, step.by.origin,
                step.by.agent, step.by.model, step.by.host, step.by.turn, _timestamp(step.at), step.status, step.id in stale,
                step.bookmark, step.group, len(step.notes),
            ]
            for step in steps
        ]
        frame = pd.DataFrame(rows, columns=list(STEP_COLUMNS))
        frame["at"] = pd.to_datetime(frame["at"], utc=True)
        frame["v"] = frame["v"].astype("Int64")
        frame["stale"] = frame["stale"].astype(bool)
        frame["n_notes"] = frame["n_notes"].astype("int64")
        return frame


def _timestamp(milliseconds: float | None) -> pd.Timestamp | None:
    return None if milliseconds is None else pd.Timestamp(milliseconds, unit="ms", tz="UTC")


@dataclass(frozen=True)
class Diff:
    """What differs between the states after two steps (of one tracker or two), domain by domain.

    ``rows`` holds ``(domain, name, line)``: the domain's id, its name as a
    person reads it, and one line of what changed, in the web app's words for
    the open cells, fields and labels, gene sets and highlights, and
    ``"<name> changed"`` for the presentation domains. ``table`` is the rows
    as a DataFrame; ``common`` is the two steps' last common step (None
    across trackers).
    """

    a: Step
    b: Step
    common: str | None
    rows: tuple[tuple[str, str, str], ...]

    @property
    def table(self) -> pd.DataFrame:
        return pd.DataFrame(list(self.rows), columns=["domain", "name", "line"])

    @property
    def domains(self) -> tuple[str, ...]:
        """The domains that differ, in domain order."""
        return tuple(dict.fromkeys(domain for domain, _, _ in self.rows))

    @classmethod
    def between(cls, a: Step, b: Step) -> Diff:
        tracker_a, tracker_b = a._tracker(), b._tracker()
        values_a, values_b = tracker_a._hashes(a.id), tracker_b._hashes(b.id)
        rows: list[tuple[str, str, str]] = []
        for domain in model.DOMAINS:
            if values_a.get(domain) == values_b.get(domain):
                continue
            state_a = tracker_a._context().value(values_a[domain]) if tracker_a._context().values_included else None
            state_b = tracker_b._context().value(values_b[domain]) if tracker_b._context().values_included else None
            lines = (
                describe_domain(domain, state_a, state_b)
                if state_a is not None and state_b is not None
                else [f"{DOMAIN_NAMES[domain]} changed"]
            )
            rows.extend((domain, DOMAIN_NAMES[domain], line) for line in lines)
        common = None
        if tracker_a is tracker_b:
            on_a = {step.id for step in tracker_a.path_to(a.id)}
            common = next(step.id for step in reversed(tracker_b.path_to(b.id)) if step.id in on_a)
        return cls(a=a, b=b, common=common, rows=tuple(rows))


def build_tracker(
    record: Mapping[str, Any],
    *,
    value: Callable[[str], Any],
    values_included: bool,
    dataset: Mapping[str, Any],
) -> Tracker:
    """A :class:`Tracker` from its record as the model's JSON (already
    checked): its steps (params and outcome with their typed arrays read in),
    branches, groups, turns, notes, journal and conversations."""
    plain = {
        **record,
        "steps": [{**step, "params": _plain_value(step["params"]), "outcome": _plain_value(step["outcome"])} for step in record["steps"]],
        "journal": [],
    }
    tracker = cast(Tracker, model.build("Tracker", plain, {"Tracker": Tracker, "Branch": Branch, "Group": Group, "Step": Step}))
    journal = tuple(
        model.build("JournalEntry", {
            **entry,
            **{
                part: None if entry[part] is None else {
                    **entry[part],
                    "steps": [{**step, "params": _plain_value(step["params"]), "outcome": _plain_value(step["outcome"])} for step in entry[part]["steps"]],
                }
                for part in ("before", "after")
            },
        })
        for entry in record["journal"]
    )
    object.__setattr__(tracker, "journal", journal)
    object.__setattr__(tracker, _READING, _Context(value=value, values_included=values_included, dataset=dataset))
    for item in (*tracker.step_records, *tracker.branches, *tracker.groups):
        object.__setattr__(item, _OWNER, tracker)
    return tracker


def _field_names(record: str) -> list[str]:
    return [attr for _, attr, _, _ in model.FIELDS[record]]


def record_fields(cls: type) -> list[str]:
    """The attribute names of a record class's fields (the model's, in its order)."""
    return [item.name for item in fields(cls)]
