"""Session files of format 4 (the History v2 shape), made in the tests.

The shared golden fixture (``tests/fixtures/history_golden``) holds every kind
of record the history model has; these helpers write files of the same format
with other content (other datasets, labels, fields, subsamples, branches) and
take files apart to damage them, so the reader's and the applier's rules are
tested on every case. They follow the format the web app writes
(``app/session/README.md``, ``app/session/history-file.js``): steps, parts
and values as JSON cut into pieces, each chunk gzip'd where that saves a
sixteenth of its bytes (else stored as it is), every array once by its
digest, every sub-value of at least 512 characters of skeleton text once as
a part, the values the active head shows and the parts and arrays they name
eager and everything else lazy, every record the model's
(``history_model.schema.json``).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import struct
from pathlib import Path
from typing import Any

import numpy as np

from cellucid._history_model import DOMAINS as DOMAIN_ORDER
from cellucid.session import _ARRAY_TYPES, SESSION_SIGNATURE, _Array, _hash_value, _js_stringify

CREATED_AT = "2026-11-03T12:00:00.000Z"
#: The time of a test tracker's opening (2026-01-01T00:00:00Z).
START = 1_767_225_600_000
#: The skeleton text from which an object or a list is stored as a part (``PART_MIN_CHARS``).
PART_MIN_CHARS = 512
#: The text a JSON chunk is cut after (``CHUNK_TEXT_CHARS``).
CHUNK_TEXT_CHARS = 16 * 1024 * 1024


def typed(type_name: str, values: Any) -> _Array:
    """A typed array as a value holds one (``Uint32Array``, ``Float32Array``...)."""
    data = np.ascontiguousarray(np.asarray(values), dtype=_ARRAY_TYPES[type_name])
    digest = f"{type_name}:{data.size}:{hashlib.sha256(data.tobytes()).hexdigest()}"
    data.flags.writeable = False
    return _Array(type=type_name, digest=digest, data=data)


def cells(values: Any) -> _Array:
    return typed("Uint32Array", values)


def mapping(**changes: Any) -> dict[str, Any]:
    out: dict[str, Any] = {
        "kind": "anndata",
        "embeddings": {
            "enabled": [{"key": "X_umap", "columns": [0, 1]}],
            "defaults": {"1d": None, "2d": {"key": "X_umap", "columns": [0, 1]}, "3d": None},
        },
        "knn": None,
        "geneNames": "index",
        "expression": "X",
        "log1p": {},
        "highlyVariable": None,
        "latent": None,
        "vectorFields": [],
    }
    out.update(changes)
    if "modality" in changes:
        # The canonical key order: modality after knn.
        ordered = {key: out[key] for key in ("kind", "embeddings", "knn", "modality")}
        ordered.update({key: value for key, value in out.items() if key not in ordered})
        out = ordered
    return out


def fingerprint(*, cell_count: int, var_count: int, dataset_id: str | None = "demo", data_mapping: dict[str, Any] | None = None,
                embedding: str = "X_umap", digest: str = "0123456789abcdef") -> dict[str, Any]:
    """A dataset fingerprint; ``digest`` is the cell-order digest of ``embedding``'s
    columns 0 and 1 over the cells (an apply checks it when the object has that
    embedding: ``cellucid.anndata_session._apply.cell_order_digest``)."""
    return {
        "datasetId": dataset_id,
        "name": "Demo",
        "cellCount": cell_count,
        "varCount": var_count,
        "cellOrder": {"embedding": {"key": embedding, "columns": [0, 1]}, "dimension": 2, "digest": digest},
        "dataMapping": mapping() if data_mapping is None else data_mapping,
    }


def overlays(*, renamed_fields: list[Any] | None = None, renamed_categories: list[Any] | None = None,
             deleted_fields: list[Any] | None = None, user_fields: list[Any] | None = None) -> dict[str, Any]:
    return {
        "renamedFields": list(renamed_fields or []),
        "renamedCategories": list(renamed_categories or []),
        "deletedFields": list(deleted_fields or []),
        "userFields": list(user_fields or []),
    }


def user_category(field_id: str, key: str, categories: list[Any], codes: Any, *, cells_of: Any = None,
                  code_type: str = "Uint8Array", deleted: bool = False, source: str = "obs") -> dict[str, Any]:
    return {
        "id": field_id, "source": source, "kind": "category", "key": key, "deleted": deleted, "purged": False,
        "createdAt": 1_767_225_600_000, "sourceField": None, "operation": {"type": "annotation"},
        "categories": categories,
        "codes": {"cells": None if cells_of is None else cells(cells_of), "values": typed(code_type, codes)},
        "sourcePages": [], "overlapStrategy": "first", "overlapLabel": None, "intersectionLabels": None,
        "uncoveredLabel": None, "palette": None,
    }


def user_values(field_id: str, key: str, values: Any, *, cells_of: Any = None, deleted: bool = False) -> dict[str, Any]:
    return {
        "id": field_id, "source": "obs", "kind": "continuous", "key": key, "deleted": deleted, "purged": False,
        "createdAt": 1_767_225_600_000, "sourceField": None, "operation": {"type": "notebook-push"},
        "values": {"cells": None if cells_of is None else cells(cells_of), "values": typed("Float32Array", values)},
    }


def user_copy(field_id: str, key: str, source: str, source_key: str, *, source_index: int | None = None) -> dict[str, Any]:
    source_field: dict[str, Any] = {"kind": "continuous-copy", "sourceKey": source_key}
    if source_index is not None:
        source_field["sourceIndex"] = source_index
    return {
        "id": field_id, "source": source, "kind": "continuous", "key": key, "deleted": False, "purged": False,
        "createdAt": 1_767_225_600_000, "sourceField": source_field, "operation": {"type": "copy"},
    }


def highlights(*pages: dict[str, Any], active: str | None = None) -> dict[str, Any]:
    return {"active": active if active is not None else (pages[0]["id"] if pages else None), "pages": list(pages)}


def page(page_id: str, name: str, *groups: dict[str, Any], color: str = "#2563eb") -> dict[str, Any]:
    return {"id": page_id, "name": name, "color": color, "groups": list(groups)}


def group(group_id: str, label: str, members: Any, *, kind: str = "lasso", enabled: bool = True,
          **extra: Any) -> dict[str, Any]:
    return {"id": group_id, "type": kind, "label": label, "enabled": enabled, "cells": cells(members), **extra}


def state(**changes: Any) -> dict[str, Any]:
    """Every domain an application reads, at its opening value, with ``changes``."""
    out: dict[str, Any] = {
        "cells": {"open": None, "genes": None},
        "mapping": mapping(),
        "field-overlays": overlays(),
        "gene-sets": {"sets": []},
        "highlights": highlights(),
        "camera": {"views": [{"zoom": -0.0, "ratio": math.nan}]},
    }
    for key, value in changes.items():
        out[key.replace("_", "-")] = value
    return out


def _split(value: Any, arrays: dict[str, _Array], parts: dict[str, Any] | None = None, *, root: bool = True) -> Any:
    """A value's skeleton (``splitValue``): typed arrays into ``arrays``, and,
    when ``parts`` is given, every object or list (not the value itself) whose
    skeleton text is at least :data:`PART_MIN_CHARS` into ``parts`` under the
    SHA-256 of that text, bottom-up."""
    if isinstance(value, _Array):
        arrays.setdefault(value.digest, value)
        return {"$array": value.digest}
    if isinstance(value, float) and (not math.isfinite(value) or (value == 0 and math.copysign(1.0, value) < 0)):
        name = "-0" if value == 0 else "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
        return {"$number": name}
    if isinstance(value, list):
        out: Any = [_split(item, arrays, parts, root=False) for item in value]
    elif isinstance(value, dict):
        out = {key: _split(child, arrays, parts, root=False) for key, child in value.items()}
    else:
        return value
    if parts is None or root:
        return out
    text = _js_stringify(out)
    if len(text) < PART_MIN_CHARS:
        return out
    name = hashlib.sha256(text.encode("utf-8")).hexdigest()
    parts.setdefault(name, out)
    return {"$part": name}


def part_names(skeleton: Any) -> list[str]:
    """The part names a skeleton names directly."""
    if isinstance(skeleton, dict):
        if set(skeleton) == {"$part"}:
            return [skeleton["$part"]]
        return [name for child in skeleton.values() for name in part_names(child)]
    if isinstance(skeleton, list):
        return [name for item in skeleton for name in part_names(item)]
    return []


def _reach(skeletons: Any, parts: dict[str, Any]) -> tuple[set[str], set[str]]:
    """Every array and part a set of skeletons reaches, through the parts."""
    arrays: set[str] = set()
    reached: set[str] = set()
    queue = list(skeletons)
    while queue:
        skeleton = queue.pop()
        arrays.update(_digests(skeleton))
        for name in part_names(skeleton):
            if name not in reached:
                reached.add(name)
                queue.append(parts[name])
    return arrays, reached


def pieces(kind: str, items: Any, limit: int = CHUNK_TEXT_CHARS) -> list[bytes]:
    """JSON text cut after ``limit`` characters (``jsonPieces``): a list
    (``kind="list"``, items) or a map (``"map"``, (key, value) pairs); one
    item is never split."""
    out: list[bytes] = []
    texts: list[str] = []
    length = 0
    open_, close = ("[", "]") if kind == "list" else ("{", "}")
    for item in items:
        text = _js_stringify(item) if kind == "list" else f"{_js_stringify(item[0])}:{_js_stringify(item[1])}"
        if texts and length + len(text) > limit:
            out.append(f"{open_}{','.join(texts)}{close}".encode())
            texts, length = [], 0
        texts.append(text)
        length += len(text) + 1
    if texts:
        out.append(f"{open_}{','.join(texts)}{close}".encode())
    return out


def _json(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")


def frame(manifest: dict[str, Any], chunks: list[bytes]) -> bytes:
    """A session file of ``manifest`` and its stored chunks, framed."""
    manifest_bytes = _json(manifest)
    out = bytearray(SESSION_SIGNATURE + struct.pack("<I", len(manifest_bytes)) + manifest_bytes)
    for chunk in chunks:
        out += struct.pack("<I", len(chunk)) + chunk
    return bytes(out)


def unframe(data: bytes) -> tuple[dict[str, Any], list[bytes]]:
    """A session file's manifest and stored chunks."""
    offset = len(SESSION_SIGNATURE)
    length = struct.unpack("<I", data[offset : offset + 4])[0]
    manifest = json.loads(data[offset + 4 : offset + 4 + length])
    offset += 4 + length
    chunks = []
    for _ in manifest["chunks"]:
        size = struct.unpack("<I", data[offset : offset + 4])[0]
        chunks.append(data[offset + 4 : offset + 4 + size])
        offset += 4 + size
    return manifest, chunks


def stored(payload: bytes) -> bytes:
    return gzip.compress(payload, mtime=0)


def store(meta: dict[str, Any], payload: bytes, codec: str = "auto") -> tuple[dict[str, Any], bytes]:
    """A chunk's manifest entry and stored bytes (``storeChunk``): with
    ``codec="auto"`` gzip'd where that saves at least a sixteenth, else as it is."""
    zipped = stored(payload)
    keep = codec == "gzip" or (codec == "auto" and len(zipped) * 16 <= len(payload) * 15)
    body = zipped if keep else payload
    return {**meta, "codec": "gzip" if keep else "none", "storedBytes": len(body), "uncompressedBytes": len(payload)}, body


def payload(entry: dict[str, Any], chunk: bytes) -> bytes:
    """A stored chunk's payload."""
    return gzip.decompress(chunk) if entry["codec"] == "gzip" else chunk


def restore_sizes(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    """Make each chunk's declared sizes its own (after a chunk was changed)."""
    for entry, chunk in zip(manifest["chunks"], chunks, strict=True):
        entry["storedBytes"] = len(chunk)
        entry["uncompressedBytes"] = len(gzip.decompress(chunk)) if entry["codec"] == "gzip" else len(chunk)


def tracker_id(label: str) -> str:
    """A tracker id (``t-<uuid>``) for a test's short name (``"t1"``), the same every time."""
    if re.fullmatch(r"t-[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", label):
        return label
    digest = hashlib.sha256(label.encode("utf-8")).hexdigest()
    return f"t-{digest[:8]}-{digest[8:12]}-4{digest[13:16]}-8{digest[17:20]}-{digest[20:32]}"


def by(origin: str = "user", *, agent: str | None = None, model: str | None = None, host: str | None = None,
       turn: str | None = None) -> dict[str, Any]:
    """A ``By`` record."""
    return {"origin": origin, "agent": agent, "model": model, "host": host, "person": None, "turn": turn}


def note(note_id: str, text: str, author: dict[str, Any] | None = None, *, at: float | None = START) -> dict[str, Any]:
    return {"id": note_id, "text": text, "by": author or by(), "at": at, "editedAt": None, "fromReply": False}



def session_file(
    path: Path,
    *,
    fingerprint_record: dict[str, Any],
    trackers: list[dict[str, Any]],
    active: str | None = None,
    requires: dict[str, list[str]] | None = None,
    cells_record: dict[str, Any] | None = None,
    genes_record: dict[str, Any] | None = None,
    source: dict[str, Any] | None = None,
    app_version: str = "0.9.6",
    content: str = "history",
    times: str = "exact",
    piece_chars: int = CHUNK_TEXT_CHARS,
    codec: str = "auto",
) -> Path:
    """Write a session file of the History v2 shape (every record the model's).

    Each tracker is ``{"id", "name", "root": {domain: value}, "steps": [...],
    "head"?, "branchNames"?: {branch id: name}, "activeBranch"?, "chats"?,
    "notes"?}``; each step ``{"id"?, "parent"?, "type"?, "set": {domain:
    value}, "params"?, "outcome"?, "origin"?, "agent"?, "host"?, "turn"?,
    "branch"?, "status"?, "reason"?, "label"?, "notes"?, "bookmark"?,
    "reads"?, "group"?, "squashed"?}``: its ``after`` the values it sets (only
    those that differ from its parent's are its touches) and its ``before``
    its parent's; its branch ``b0`` unless named; a step with ``squashed``
    (``{cause, count, by}``, ``carried`` empty unless given) is a squashed
    step, with no action. A tracker's ``groups`` are ``[{"id", "name"}]``,
    its ``recording`` ``False`` pauses it (``paused`` since its opening, by
    the person).
    Branches are made from the steps on them, and turns from the steps' turn
    ids. ``content="steps"`` writes the steps only (no values); ``"state"``
    marks a file holding a state without its history (its one tracker must be
    the root and one squashed step of cause ``share``).

    The chunks are laid out as the web app's writer lays them out: JSON cut
    after ``piece_chars`` characters of text, each chunk stored with
    ``codec`` (``"auto"``: gzip where it saves a sixteenth; ``"gzip"``;
    ``"none"``).
    """
    arrays: dict[str, _Array] = {}
    parts: dict[str, Any] = {}
    held: dict[str, Any] = {}
    records: list[tuple[dict[str, Any], list[dict[str, Any]], dict[str, dict[str, str]]]] = []
    types: set[tuple[str, int]] = set()
    for index, tracker in enumerate(trackers):
        identifier = tracker_id(tracker["id"])
        hashes_at: dict[str, dict[str, str]] = {}
        stored_steps: list[dict[str, Any]] = []
        created = START + 1_000_000 * index
        root_after = {}
        for domain in DOMAIN_ORDER:
            if domain in tracker["root"]:
                value_hash = _hash_value(tracker["root"][domain])
                held[value_hash] = tracker["root"][domain]
                root_after[domain] = value_hash
        root = {
            "id": "s0", "parent": None, "children": [], "branch": "b0", "kind": "root", "type": None, "params": None,
            "v": None, "code": None, "label": "Opened", "by": by(), "at": created, "until": created, "merged": 1, "ms": None,
            "touches": [], "before": {}, "after": root_after, "reads": {}, "status": "applied", "reason": None,
            "outcome": None, "coalesceKey": None, "replayOf": None, "squashed": None, "notes": [], "bookmark": None,
            "group": None,
        }
        stored_steps.append(root)
        hashes_at["s0"] = dict(root_after)
        by_id = {"s0": root}
        for number, step in enumerate(tracker.get("steps", []), start=1):
            parent = step.get("parent", stored_steps[-1]["id"])
            after: dict[str, str] = {}
            before: dict[str, str] = {}
            for domain in DOMAIN_ORDER:
                if domain not in step["set"]:
                    continue
                value_hash = _hash_value(step["set"][domain])
                held[value_hash] = step["set"][domain]
                if step.get("status", "applied") == "applied" and hashes_at[parent].get(domain) == value_hash:
                    continue
                after[domain] = value_hash
                before[domain] = hashes_at[parent][domain]
            if step.get("status", "applied") == "unapplicable":
                after = {domain: hashes_at[parent][domain] for domain in after}
            squashed = step.get("squashed")
            if squashed is not None:
                squashed = {"carried": {"notes": [], "bookmarks": []}, **squashed}
            action = None if squashed is not None else step.get("type", "test.step")
            origin = step.get("origin", "agent" if step.get("agent") else "host" if step.get("host") else "user")
            kind = "squashed" if squashed is not None else "outside" if origin == "outside" else "action"
            if kind == "action":
                types.add((action, step.get("v", 1)))
            record = {
                "id": step.get("id", f"s{number}"), "parent": parent, "children": [], "branch": step.get("branch", "b0"),
                "kind": kind, "type": action,
                "params": None if squashed is not None else _split(step.get("params", {}), arrays, parts),
                "v": step.get("v", 1) if kind == "action" else None,
                "code": None if kind != "action" else f"{action.replace('.', '_')}()",
                "label": step.get("label", f"Step {number}"),
                "by": by(origin, agent=step.get("agent"), model=step.get("model"), host=step.get("host"), turn=step.get("turn")),
                "at": created + number * 1000, "until": created + number * 1000, "merged": 1, "ms": 1.5,
                "touches": [domain for domain in DOMAIN_ORDER if domain in after], "before": before, "after": after,
                "reads": step.get("reads", {}), "status": step.get("status", "applied"), "reason": step.get("reason"),
                "outcome": _split(step.get("outcome"), arrays, parts), "coalesceKey": None, "replayOf": None, "squashed": squashed,
                "notes": list(step.get("notes", [])), "bookmark": step.get("bookmark"), "group": step.get("group"),
            }
            by_id[parent]["children"].append(record["id"])
            by_id[record["id"]] = record
            stored_steps.append(record)
            hashes_at[record["id"]] = {**hashes_at[parent], **after}
        head = tracker.get("head", stored_steps[-1]["id"])
        names = tracker.get("branchNames", {})
        branches = []
        for branch_id in dict.fromkeys(step["branch"] for step in sorted(stored_steps, key=lambda item: int(item["id"][1:]))):
            own = [step for step in stored_steps if step["branch"] == branch_id]
            own.sort(key=lambda item: int(item["id"][1:]))
            branches.append({
                "id": branch_id,
                "name": names.get(branch_id, "Main" if branch_id == "b0" else f"Branch {branch_id[1:]}"),
                "base": own[0]["parent"], "tip": own[-1]["id"], "createdAt": own[0]["at"], "by": by(),
                "lane": int(branch_id[1:]) % 8, "locked": False, "archived": False, "notes": [],
            })
        branches.sort(key=lambda item: int(item["id"][1:]))
        turns = []
        for turn_id in dict.fromkeys(step["by"]["turn"] for step in stored_steps if step["by"]["turn"] is not None):
            taken = [step for step in stored_steps if step["by"]["turn"] == turn_id]
            turns.append({
                "id": turn_id, "by": dict(taken[0]["by"]), "label": tracker.get("turnLabels", {}).get(turn_id, "A turn"),
                "conversation": None, "message": None, "startedAt": taken[0]["at"], "endedAt": taken[-1]["at"],
                "branch": taken[0]["branch"], "notes": [],
            })
        notes = list(tracker.get("notes", []))
        every_note = [*notes, *(item for step in stored_steps for item in step["notes"])]
        manifest_record = {
            "id": identifier, "name": tracker["name"],
            "dataset": {"id": "demo", "name": "Demo", "cellOrder": None}, "createdAt": created, "by": by(),
            "start": {"kind": "opened", "from": None, "report": None, "mapping": None}, "removed": None,
            "rootId": "s0", "headId": head,
            "activeBranch": tracker.get("activeBranch", by_id[head]["branch"]), "mainBranch": "b0",
            "paused": None if tracker.get("recording", True) else {"at": created, "by": [by()], "draft": None},
            "branches": branches,
            "groups": [{"id": item["id"], "name": item["name"], "createdAt": created, "by": by(), "notes": []} for item in tracker.get("groups", [])],
            "turns": turns, "notes": notes,
            "nextStep": max(int(step["id"][1:]) for step in stored_steps) + 1,
            "nextBranch": max(int(branch["id"][1:]) for branch in branches) + 1,
            "nextGroup": max([int(item["id"][1:]) + 1 for item in tracker.get("groups", [])], default=0),
            "nextNote": max([int(item["id"][1:]) + 1 for item in every_note], default=0),
            "nextJournal": 0, "steps": [], "journal": None, "chats": [],
        }
        stored_steps.sort(key=lambda item: int(item["id"][1:]))
        records.append((manifest_record, stored_steps, hashes_at))
    active_id = tracker_id(active) if active is not None else records[0][0]["id"]
    active_record = next(record for record in records if record[0]["id"] == active_id)
    referenced = {
        value_hash for _, steps, _ in records for step in steps for value_hash in [*step["before"].values(), *step["after"].values()]
    }
    eager_hashes = list(dict.fromkeys(active_record[2][active_record[0]["headId"]].values()))
    values = content != "steps"
    eager_skeletons = {value_hash: _split(held[value_hash], arrays, parts) for value_hash in eager_hashes} if values else {}
    lazy_skeletons = (
        {value_hash: _split(held[value_hash], arrays, parts) for value_hash in sorted(referenced) if value_hash not in eager_skeletons}
        if values else {}
    )
    eager_arrays, eager_parts = _reach(eager_skeletons.values(), parts)

    chunk_entries: list[tuple[dict[str, Any], bytes]] = []
    for index, (manifest_record, steps, _) in enumerate(records):
        for piece, body in enumerate(pieces("list", steps, piece_chars)):
            chunk_id = f"steps/{index}/{piece}"
            manifest_record["steps"].append(chunk_id)
            chunk_entries.append(({"id": chunk_id, "kind": "steps", "tracker": manifest_record["id"], "priority": "eager"}, body))
    for priority in ("eager", "lazy"):
        own_parts = [(name, skeleton) for name, skeleton in parts.items() if (name in eager_parts) == (priority == "eager")]
        for piece, body in enumerate(pieces("map", own_parts, piece_chars)):
            chunk_entries.append(({"id": f"parts/{priority}/{piece}", "kind": "parts", "priority": priority}, body))
        for digest in [digest for digest in arrays if (digest in eager_arrays) == (priority == "eager")]:
            chunk_entries.append(_array_chunk(arrays[digest], priority))
        own_values = eager_skeletons if priority == "eager" else lazy_skeletons
        for piece, body in enumerate(pieces("map", list(own_values.items()), piece_chars)):
            chunk_entries.append(({"id": f"values/{priority}/{piece}", "kind": "values", "priority": priority}, body))
    for index, tracker in enumerate(trackers):
        for position, (conversation, transcript) in enumerate(tracker.get("chats", [])):
            chunk_id = f"chat/{index}/{position}"
            records[index][0]["chats"].append(chunk_id)
            chunk_entries.append(({"id": chunk_id, "kind": "chat", "tracker": records[index][0]["id"],
                                   "conversation": conversation, "priority": "lazy"}, _json(transcript)))
    metas = []
    chunks = []
    for meta, body in chunk_entries:
        entry, framed = store(meta, body, codec)
        metas.append(entry)
        chunks.append(framed)
    manifest: dict[str, Any] = {
        "createdAt": CREATED_AT,
        "app": {"version": app_version, "catalog": "ab" * 32},
        "actions": [{"type": action, "version": version, "params": "cd" * 32} for action, version in sorted(types)],
        "scope": {
            "trackers": [record[0]["id"] for record in records],
            "branches": [{"tracker": record[0]["id"], "branches": [branch["id"] for branch in record[0]["branches"]]} for record in records],
            "content": content, "results": "include" if values else "recompute", "notes": content != "state",
            "chats": any(tracker.get("chats") for tracker in trackers), "journal": False, "times": times, "source": "full",
        },
        "datasetFingerprint": fingerprint_record,
        "source": source,
        "requires": requires if requires is not None else {"obs": [], "var": []},
        "chrome": {"pressed": {}},
        **({"cells": cells_record} if cells_record is not None else {}),
        **({"genes": genes_record} if genes_record is not None else {}),
        "activeTracker": active_id,
        "trackers": [record[0] for record in records],
        "chunks": metas,
    }
    path.write_bytes(frame(manifest, chunks))
    return path


def _digests(skeleton: Any) -> list[str]:
    if isinstance(skeleton, dict):
        if set(skeleton) == {"$array"}:
            return [skeleton["$array"]]
        return [digest for child in skeleton.values() for digest in _digests(child)]
    if isinstance(skeleton, list):
        return [digest for item in skeleton for digest in _digests(item)]
    return []


def _array_chunk(array: _Array, priority: str) -> tuple[dict[str, Any], bytes]:
    return (
        {"id": f"array/{array.digest}", "kind": "array", "digest": array.digest, "encoding": "raw",
         "priority": priority},
        array.data.tobytes(),
    )


def subsample_record(n_obs: int, k: int, seed: int) -> tuple[dict[str, Any], np.ndarray]:
    """The manifest ``cells`` of a uniform subsample, and its kept cells."""
    import base64

    from cellucid.subsample import cell_set_digest, encode_cell_set, subsample_cells

    kept = subsample_cells(n_obs, k, seed=seed)
    return (
        {
            "spec": {"method": "uniform", "nObs": n_obs, "k": k, "seed": seed},
            "kept": base64.b64encode(encode_cell_set(kept)).decode("ascii"),
            "digest": cell_set_digest(n_obs, kept),
        },
        kept,
    )
