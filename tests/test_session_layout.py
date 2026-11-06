"""Session format 4's chunk layout, read as the web app writes it.

The web app's writer (``app/session/history-file.js``) cuts steps, parts and
values into pieces, stores each chunk gzip'd or as it is, keeps every
sub-value of at least 512 characters of skeleton text once as a part, and
leaves the arrays and parts only step parameters name in the lazy chunks
(``app/session/README.md``, "The file"). These files are written by the test
helper the same way (``_session_files.py``) and then taken apart: parts
shared, cyclic, expanding past what a page holds, unreferenced or renamed,
steps and values in many pieces, chunks stored as they are, and parameters
naming lazy arrays.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import _session_files as files
import pytest

from cellucid import read_session
from cellucid.session import _KIND_ORDER, MAX_VALUE_TEXT

GENES = [f"GENE{index:04d}" for index in range(120)]
OTHER_GENES = [f"OTHER{index:04d}" for index in range(120)]


def _fingerprint() -> dict[str, Any]:
    return files.fingerprint(cell_count=6, var_count=len(GENES), dataset_id="layout")


def _gene_sets(*names: str) -> dict[str, Any]:
    """Gene sets whose gene lists are long enough to be stored as parts."""
    return {"sets": [{"name": name, "genes": OTHER_GENES if name == "Third" else GENES} for name in names]}


def _write(path: Path, **options: Any) -> Path:
    """Three steps whose gene-sets values share the first set, and a step whose
    parameters hold a cell array no shown value holds."""
    root = files.state(gene_sets=_gene_sets("Markers"))
    steps = [
        {"set": {"gene-sets": _gene_sets("Markers", "More")}, "type": "gene-set.add", "params": {"name": "More", "genes": GENES}},
        {"set": {"highlights": files.highlights(files.page("p0", "Page 1", files.group("g0", "Lasso", [0, 1, 2])))},
         "type": "highlight.add", "params": {"cells": files.cells([0, 1, 2, 3])}},
        {"set": {"gene-sets": _gene_sets("Markers", "More", "Third")}, "type": "gene-set.add", "params": {"name": "Third"}},
    ]
    return files.session_file(
        path, fingerprint_record=_fingerprint(),
        trackers=[{"id": "t1", "name": "Tracker 1", "root": root, "steps": steps}], **options,
    )


def _rewrite(path: Path, edit: Any, out: Path) -> Path:
    manifest, chunks = files.unframe(path.read_bytes())
    edit(manifest, chunks)
    out.write_bytes(files.frame(manifest, chunks))
    return out


def _index(manifest: dict[str, Any], chunk_id: str) -> int:
    return next(index for index, entry in enumerate(manifest["chunks"]) if entry["id"] == chunk_id)


def _json_chunk(manifest: dict[str, Any], chunks: list[bytes], chunk_id: str) -> Any:
    index = _index(manifest, chunk_id)
    return json.loads(files.payload(manifest["chunks"][index], chunks[index]))


def _put_json_chunk(manifest: dict[str, Any], chunks: list[bytes], chunk_id: str, payload: Any) -> None:
    index = _index(manifest, chunk_id)
    entry = manifest["chunks"][index]
    manifest["chunks"][index], chunks[index] = files.store(entry, json.dumps(payload, separators=(",", ":")).encode(), entry["codec"])


def _insert_chunk(manifest: dict[str, Any], chunks: list[bytes], meta: dict[str, Any], payload: bytes) -> None:
    """A chunk put where the order puts it (eager before lazy, then by kind)."""
    def key(entry: dict[str, Any]) -> tuple[int, int]:
        return (entry["priority"] == "lazy", _KIND_ORDER.index(entry["kind"]))

    at = next((index for index, entry in enumerate(manifest["chunks"]) if key(entry) > key(meta)), len(manifest["chunks"]))
    entry, body = files.store(meta, payload)
    manifest["chunks"].insert(at, entry)
    chunks.insert(at, body)


def _eager_value_naming_a_part(manifest: dict[str, Any], chunks: list[bytes]) -> tuple[str, str]:
    """An eager value's hash and the first part its skeleton names."""
    values = _json_chunk(manifest, chunks, "values/eager/0")
    for value_hash, skeleton in values.items():
        names = files.part_names(skeleton)
        if names:
            return value_hash, names[0]
    raise AssertionError("no eager value names a part")


def _refused(path: Path) -> str:
    with pytest.raises(ValueError) as refused:
        read_session(path)
    message = str(refused.value)
    assert message.startswith("This session file is damaged: ")
    return message


# =============================================================================
# Parts
# =============================================================================


def test_parts_are_stored_once_and_joined_once(tmp_path: Path) -> None:
    path = _write(tmp_path / "parts.cellucid-session")
    manifest, chunks = files.unframe(path.read_bytes())
    kinds = [entry["kind"] for entry in manifest["chunks"]]
    assert "parts" in kinds
    parts = {name: skeleton for entry, chunk in zip(manifest["chunks"], chunks, strict=True) if entry["kind"] == "parts"
             for name, skeleton in json.loads(files.payload(entry, chunk)).items()}
    # The gene list of "Markers" and "More" is held by the root's, s1's and s3's values and s1's
    # parameters: one part; "Third"'s another.
    assert sorted(skeleton[0] for skeleton in parts.values()) == ["GENE0000", "OTHER0000"]

    session = read_session(path)
    tracker = session.tracker()
    assert tracker.step("s0").state()["gene-sets"] == _gene_sets("Markers")
    assert tracker.step("s3").state()["gene-sets"] == _gene_sets("Markers", "More", "Third")
    # Joined once: every value that names the part holds the same joined sub-value.
    values = session._contents.values
    held = [item["genes"] for value in values.values() if isinstance(value, dict) and value.get("sets")
            for item in value["sets"] if item["name"] != "Third"]
    assert len(held) == 5 and all(item is held[0] for item in held)
    assert tracker.step("s1").params == {"name": "More", "genes": GENES}


def test_a_part_that_names_itself_through_others_is_refused(tmp_path: Path) -> None:
    def cycle(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        _, name = _eager_value_naming_a_part(manifest, chunks)
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        other = "f" * 64
        parts[other] = [{"$part": name}]
        parts[name] = {"loop": {"$part": other}}
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), cycle, tmp_path / "cycle.cellucid-session"))
    assert "its eager values: A stored value names a part inside itself." in message


def test_a_value_that_would_expand_past_a_page_is_refused_before_it_is_expanded(tmp_path: Path) -> None:
    def bomb(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        value_hash, _ = _eager_value_naming_a_part(manifest, chunks)
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        names = [f"{level:064x}" for level in range(24)]
        parts[names[0]] = ["x" * 1000]
        for level in range(1, len(names)):
            parts[names[level]] = [{"$part": names[level - 1]}, {"$part": names[level - 1]}]
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)
        values = _json_chunk(manifest, chunks, "values/eager/0")
        values[value_hash] = {"sets": [{"$part": names[-1]}]}
        _put_json_chunk(manifest, chunks, "values/eager/0", values)

    path = _rewrite(_write(tmp_path / "a.cellucid-session"), bomb, tmp_path / "bomb.cellucid-session")
    started = time.perf_counter()
    message = _refused(path)
    # 2^23 copies of 1000 characters: counted, never made.
    assert time.perf_counter() - started < 10
    assert "A stored value expands to" in message
    assert f"more than a page can hold ({MAX_VALUE_TEXT})" in message
    assert MAX_VALUE_TEXT == 2**29 - 24


def test_a_part_nothing_refers_to_is_refused(tmp_path: Path) -> None:
    def extra(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        parts["0" * 64] = ["unused"]
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), extra, tmp_path / "extra.cellucid-session"))
    assert message.endswith("it holds a part nothing refers to.")


def test_a_parts_name_is_never_trusted(tmp_path: Path) -> None:
    original = _write(tmp_path / "a.cellucid-session")
    manifest, chunks = files.unframe(original.read_bytes())
    _, name = _eager_value_naming_a_part(manifest, chunks)

    def rename(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        # Another name everywhere: the values still hash as written, so the file reads.
        for index, entry in enumerate(manifest["chunks"]):
            if entry["kind"] in ("steps", "parts", "values"):
                text = files.payload(entry, chunks[index]).decode().replace(name, "e" * 64)
                manifest["chunks"][index], chunks[index] = files.store(entry, text.encode(), entry["codec"])

    _same(_reading(_rewrite(original, rename, tmp_path / "renamed.cellucid-session")), _reading(original))

    def tamper(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        parts[name] = ["Tampered", *parts[name][1:]]
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)

    message = _refused(_rewrite(original, tamper, tmp_path / "tampered.cellucid-session"))
    assert message.endswith("its eager values holds a value that does not match its hash.")


def test_a_shown_value_naming_a_lazy_part_is_refused(tmp_path: Path) -> None:
    def lazy(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        _, name = _eager_value_naming_a_part(manifest, chunks)
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        moved = parts.pop(name)
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)
        _insert_chunk(manifest, chunks, {"id": "parts/lazy/9", "kind": "parts", "priority": "lazy"},
                      json.dumps({name: moved}, separators=(",", ":")).encode())

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), lazy, tmp_path / "lazy.cellucid-session"))
    assert "a value of the step it shows names a part that is not among its first chunks" in message


def test_a_part_the_file_lacks_is_refused(tmp_path: Path) -> None:
    def drop(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        _, name = _eager_value_naming_a_part(manifest, chunks)
        parts = _json_chunk(manifest, chunks, "parts/eager/0")
        del parts[name]
        _put_json_chunk(manifest, chunks, "parts/eager/0", parts)

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), drop, tmp_path / "drop.cellucid-session"))
    assert "a value names a part the file does not hold" in message


@pytest.mark.parametrize(
    ("placeholder", "message"),
    [
        pytest.param({"$part": 5}, "A stored value names a part wrongly at value.", id="name"),
        # The placeholder's shape is checked before the part is looked for.
        pytest.param({"$part": "not-a-hash"}, "A stored value names a part wrongly at value.", id="malformed-name"),
        pytest.param({"$part": ""}, "A stored value names a part wrongly at value.", id="empty-name"),
        pytest.param({"$part": "f" * 64}, "a value names a part the file does not hold (ffffffffffff", id="unknown"),
        pytest.param({"$array": ""}, "a value names an array the file does not hold (", id="empty-digest"),
        pytest.param({"$part": "a" * 64, "x": 1}, "A stored value has a malformed placeholder at value.", id="extra-key"),
    ],
)
def test_a_malformed_part_reference_is_refused(tmp_path: Path, placeholder: dict[str, Any], message: str) -> None:
    def malformed(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        value_hash, _ = _eager_value_naming_a_part(manifest, chunks)
        values = _json_chunk(manifest, chunks, "values/eager/0")
        values[value_hash] = {"sets": [placeholder]}
        _put_json_chunk(manifest, chunks, "values/eager/0", values)

    assert message in _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), malformed, tmp_path / "bad.cellucid-session"))


# =============================================================================
# Pieces, codecs and lazy parameters
# =============================================================================


def _reading(path: Path) -> Any:
    session = read_session(path)
    tracker = session.tracker()
    return (
        tracker.steps.drop(columns=[]).to_dict("records"),
        {step: tracker.step(step).state() for step in ("s0", "s1", "s2", "s3")},
        tracker.step("s2").params["cells"].tolist(),
    )


def _same(a: Any, b: Any) -> None:
    assert a[0] == b[0] and a[2] == b[2]
    for step in a[1]:
        assert json.dumps(a[1][step], default=lambda array: array.tolist(), sort_keys=True) == json.dumps(
            b[1][step], default=lambda array: array.tolist(), sort_keys=True)


def test_steps_and_values_in_many_pieces_read_as_one(tmp_path: Path) -> None:
    whole = _write(tmp_path / "whole.cellucid-session")
    cut = _write(tmp_path / "cut.cellucid-session", piece_chars=600)
    manifest, _ = files.unframe(cut.read_bytes())
    ids = [entry["id"] for entry in manifest["chunks"]]
    assert len(manifest["trackers"][0]["steps"]) > 1
    assert manifest["trackers"][0]["steps"] == [chunk_id for chunk_id in ids if chunk_id.startswith("steps/")]
    assert len([chunk_id for chunk_id in ids if chunk_id.startswith("values/")]) > 2
    assert len([chunk_id for chunk_id in ids if chunk_id.startswith("parts/")]) > 1
    _same(_reading(cut), _reading(whole))


def test_chunks_stored_as_they_are_read_the_same(tmp_path: Path) -> None:
    plain = _write(tmp_path / "plain.cellucid-session", codec="none")
    manifest, _ = files.unframe(plain.read_bytes())
    assert {entry["codec"] for entry in manifest["chunks"]} == {"none"}
    assert all(entry["storedBytes"] == entry["uncompressedBytes"] for entry in manifest["chunks"])
    _same(_reading(plain), _reading(_write(tmp_path / "gzip.cellucid-session", codec="gzip")))

    def dishonest(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        manifest["chunks"][0]["uncompressedBytes"] += 1

    message = _refused(_rewrite(plain, dishonest, tmp_path / "dishonest.cellucid-session"))
    assert message.endswith("chunk 0 is stored as it is, at another size.")


def test_parameters_naming_a_lazy_array_are_joined(tmp_path: Path) -> None:
    path = _write(tmp_path / "lazy.cellucid-session")
    manifest, _ = files.unframe(path.read_bytes())
    digest = files.cells([0, 1, 2, 3]).digest
    entry = next(entry for entry in manifest["chunks"] if entry["kind"] == "array" and entry["digest"] == digest)
    assert entry["priority"] == "lazy"
    step = read_session(path).tracker().step("s2")
    assert step.params["cells"].tolist() == [0, 1, 2, 3]
    assert str(step.params["cells"].dtype) == "uint32"


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        pytest.param(lambda manifest, chunks: manifest["trackers"][0].__setitem__("steps", []),
                     "tracker 0 steps must list its steps chunks.", id="no-steps-chunk"),
        pytest.param(lambda manifest, chunks: manifest["trackers"][0].__setitem__("steps", "steps/0/0"),
                     "tracker 0 steps must list its steps chunks.", id="steps-as-text"),
        pytest.param(lambda manifest, chunks: manifest["trackers"][0]["steps"].append("values/eager/0"),
                     "names a steps chunk that is not its own.", id="names-another-chunk"),
        pytest.param(lambda manifest, chunks: manifest["trackers"][0]["steps"].append(manifest["trackers"][0]["steps"][0]),
                     "names a steps chunk that is not its own.", id="names-a-chunk-twice"),
    ],
)
def test_the_layout_is_refused_when_broken(tmp_path: Path, edit: Any, message: str) -> None:
    refused = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), edit, tmp_path / "broken.cellucid-session"))
    assert message in refused


def test_steps_chunks_named_in_another_order_than_the_file_holds_them_are_refused(tmp_path: Path) -> None:
    cut = _write(tmp_path / "cut.cellucid-session", piece_chars=600)

    def swap(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        steps = manifest["trackers"][0]["steps"]
        steps[0], steps[1] = steps[1], steps[0]

    message = _refused(_rewrite(cut, swap, tmp_path / "swapped.cellucid-session"))
    assert message.endswith(f'tracker "{files.tracker_id("t1")}" names its steps chunks in another order than the file holds them.')


def test_a_steps_only_file_holds_the_parts_and_arrays_its_steps_name(tmp_path: Path) -> None:
    """A share of the steps only holds no values, but the arrays and parts its steps' parameters still name."""
    path = _write(tmp_path / "steps.cellucid-session", content="steps")
    manifest, _ = files.unframe(path.read_bytes())
    kinds = {entry["kind"] for entry in manifest["chunks"]}
    assert "values" not in kinds and {"parts", "array"} <= kinds
    tracker = read_session(path).tracker()
    assert tracker.step("s1").params["genes"] == GENES
    assert tracker.step("s2").params["cells"].tolist() == [0, 1, 2, 3]


def test_a_file_without_values_of_the_step_it_shows_is_refused(tmp_path: Path) -> None:
    def drop(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        index = _index(manifest, "values/eager/0")
        manifest["chunks"].pop(index)
        chunks.pop(index)

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), drop, tmp_path / "drop.cellucid-session"))
    assert message.endswith("it holds no values of the step it shows.")


def test_lazy_before_eager_is_out_of_order(tmp_path: Path) -> None:
    def swap(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        index = _index(manifest, "values/eager/0")
        manifest["chunks"].insert(index - 1, manifest["chunks"].pop(index))
        chunks.insert(index - 1, chunks.pop(index))

    message = _refused(_rewrite(_write(tmp_path / "a.cellucid-session"), swap, tmp_path / "swap.cellucid-session"))
    assert 'its chunk "' in message and "is out of order (eager before lazy; steps, parts, array, values, journal, chat)." in message


def test_the_layout_constants_are_the_web_apps() -> None:
    """The kind order and the longest value are ``history-file.js``'s and ``value-parts.js``'s."""
    from test_sessions import _web_source

    writer = _web_source("app/session/history-file.js")
    assert "const KIND_ORDER = Object.freeze(['steps', 'parts', 'array', 'values', 'journal', 'chat']);" in writer
    assert _KIND_ORDER == ("steps", "parts", "array", "values", "journal", "chat")
    assert "export const MAX_VALUE_TEXT = 2 ** 29 - 24;" in _web_source("app/history/value-parts.js")
    assert "export const PART_MIN_CHARS = 512;" in _web_source("app/history/value-parts.js")
    assert files.PART_MIN_CHARS == 512
